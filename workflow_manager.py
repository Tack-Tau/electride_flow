#!/usr/bin/env python3
"""
VASP Workflow Manager - Batch submission with dynamic monitoring
No MongoDB/FireWorks - uses local JSON database for job tracking

Features:
- Sequential VASP workflow: Relax → SC → (PARCHG for semiconductors) → ELF
- PARCHG calculations for semiconductors
- Dynamic job submission with concurrency control

Note: Pre-screening should be done separately with prescreen.py before running this workflow.
"""

import os
import sys
import json
import time
import argparse
import zipfile
import warnings
import subprocess
import shutil
from pathlib import Path
from io import StringIO
from datetime import datetime

from pymatgen.core import Structure
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.ase import AseAtomsAdaptor

try:
    from pyxtal import pyxtal
    PYXTAL_AVAILABLE = True
except ImportError:
    PYXTAL_AVAILABLE = False
    print("WARNING: PyXtal not available. Structures will not be symmetrized.")

warnings.filterwarnings('ignore', category=UserWarning, message='.*POTCAR data with symbol.*')
warnings.filterwarnings('ignore', message='Using UFloat objects with std_dev==0')


def check_electronic_convergence_outcar(outcar_path):
    """
    Check electronic convergence from OUTCAR file.
    
    For timed-out jobs, vasprun.xml is incomplete/corrupted. This function checks
    OUTCAR for "aborting loop because EDIFF is reached" marker, which indicates
    electronic SCF converged in at least one ionic step.
    
    Returns:
        bool: True if electronic convergence was achieved
    """
    if not outcar_path.exists():
        return False
    
    try:
        with open(outcar_path, 'r') as f:
            content = f.read()
        return 'aborting loop because EDIFF is reached' in content
    except Exception:
        return False


def parse_band_gap_from_vasprun(vasprun_path):
    """
    Parse band gap from vasprun.xml using pymatgen's vasprun parser.
    
    Returns:
        tuple: (band_gap, is_semiconductor)
    """
    try:
        vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=True)
        
        if not vr.converged_electronic:
            print(f"    Warning: Electronic SCF not converged")
            return None, False
        
        if vr.efermi is None:
            print(f"    Warning: Fermi energy (e_fermi) not found in vasprun.xml")
            print(f"      This usually means the SC calculation didn't complete properly")
            print(f"      Check: {vasprun_path.parent}/OUTCAR for errors")
            return None, False
        
        band_structure = vr.get_band_structure(line_mode=False, efermi=vr.efermi)
        
        if band_structure is None:
            print(f"    Warning: Band structure not available")
            print(f"      VASP finished but didn't calculate eigenvalues properly")
            return None, False
        
        band_gap_dict = band_structure.get_band_gap()
        if band_gap_dict is None or 'energy' not in band_gap_dict:
            print(f"    Warning: Band gap data not available in band structure")
            return None, False
        
        band_gap = band_gap_dict['energy']
        is_semiconductor = not band_structure.is_metal()
        
        print(f"    e_fermi = {vr.efermi:.4f} eV, band_gap = {band_gap:.4f} eV, is_metal = {band_structure.is_metal()}")
        return band_gap, is_semiconductor
        
    except Exception as e:
        error_msg = str(e).split('\n')[0]
        print(f"    Warning: Could not parse vasprun.xml: {error_msg}")
        print(f"      File: {vasprun_path}")
        return None, False


def get_parchg_incar_settings_band(iband, grid):
    """Get INCAR settings for band-specific PARCHG."""
    return {
        'IBAND': iband,
        'LPARD': True,
        'LAECHG': False,
        'LVHAR': False,
        'LCHARG': False,
        'NGXF': grid[0],
        'NGYF': grid[1],
        'NGZF': grid[2],
        'LWAVE': False,
        'ISTART': 1,
        'ICHARG': 11
    }


def get_parchg_incar_settings_energy(eint, grid):
    """Get INCAR settings for energy-window PARCHG."""
    return {
        'EINT': eint,
        'NBMOD': -3,
        'LPARD': True,
        'LAECHG': False,
        'LVHAR': False,
        'LCHARG': False,
        'NGXF': grid[0],
        'NGYF': grid[1],
        'NGZF': grid[2],
        'LWAVE': False,
        'ISTART': 1,
        'ICHARG': 11
    }


class WorkflowDatabase:
    """Simple JSON-based database for tracking job states."""
    
    def __init__(self, db_path):
        self.db_path = Path(db_path)
        self.data = {'structures': {}, 'config': {}}
        self.load()
    
    def load(self):
        """Load database from JSON file."""
        if self.db_path.exists():
            with open(self.db_path, 'r') as f:
                self.data = json.load(f)
    
    def save(self):
        """Save database to JSON file."""
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = self.db_path.with_suffix('.tmp')
        with open(tmp_path, 'w') as f:
            json.dump(self.data, f, indent=2)
        tmp_path.replace(self.db_path)
    
    def add_structure(self, struct_id, comp_name, struct_idx, base_dir, chemsys=None):
        """Add a new structure to track."""
        self.data['structures'][struct_id] = {
            'composition': comp_name,
            'chemsys': chemsys,
            'structure_idx': struct_idx,
            'state': 'PENDING',
            'relax_job_id': None,
            'sc_job_id': None,
            'elf_job_id': None,
            'parchg_job_ids': [],
            'band_gap': None,
            'is_semiconductor': None,
            'relax_dir': str(base_dir / struct_id / 'Relax'),
            'sc_dir': str(base_dir / struct_id / 'SC'),
            'elf_dir': str(base_dir / struct_id / 'ELF'),
            'parchg_dir': str(base_dir / struct_id / 'PARCHG'),
            'last_updated': datetime.now().isoformat(),
            'error': None
        }
        self.save()
    
    def update_state(self, struct_id, state, **kwargs):
        """Update structure state and additional fields."""
        if struct_id in self.data['structures']:
            self.data['structures'][struct_id]['state'] = state
            self.data['structures'][struct_id]['last_updated'] = datetime.now().isoformat()
            for key, value in kwargs.items():
                self.data['structures'][struct_id][key] = value
            self.save()
    
    def get_structure(self, struct_id):
        """Get structure data."""
        return self.data['structures'].get(struct_id)
    
    def get_by_state(self, state):
        """Get all structures in a specific state."""
        return [sid for sid, sdata in self.data['structures'].items() 
                if sdata['state'] == state]
    
    def get_running_count(self):
        """Count structures currently running (any stage)."""
        running_states = ['RELAX_RUNNING', 'SC_RUNNING', 'ELF_RUNNING', 'PARCHG_RUNNING']
        return sum(1 for s in self.data['structures'].values() 
                   if s['state'] in running_states)
    
    def get_stats(self):
        """Get overall statistics."""
        states = {}
        for s in self.data['structures'].values():
            state = s['state']
            states[state] = states.get(state, 0) + 1
        return {
            'total': len(self.data['structures']),
            'states': states,
            'running': self.get_running_count()
        }


class VASPWorkflowManager:
    """Manages VASP job submission and monitoring with batch control."""
    
    def __init__(self, db_path, max_concurrent=10, check_interval=60):
        self.db = WorkflowDatabase(db_path)
        self.max_concurrent = max_concurrent
        self.check_interval = check_interval
    
    def read_structures_from_zip(self, zip_path, max_structures=None):
        """Read CIF structures from zip file."""
        structures = []
        with zipfile.ZipFile(zip_path, 'r') as zf:
            cif_files = sorted([f for f in zf.namelist() if f.endswith('.cif')])
            if max_structures:
                cif_files = cif_files[:max_structures]
            
            for cif_file in cif_files:
                try:
                    with zf.open(cif_file) as f:
                        cif_content = f.read().decode('utf-8')
                        parser = CifParser(StringIO(cif_content))
                        structure = parser.parse_structures(primitive=True)[0]
                        structures.append(structure)
                except Exception as e:
                    print(f"  Warning: Could not parse {cif_file}: {e}")
        return structures
    
    def create_vasp_inputs(self, structure, job_dir, job_type='relax', parchg_settings=None):
        """
        Create VASP input files using pymatgen.
        
        Note: 
        - For relax jobs, structure is symmetrized using PyXtal before VASP input generation
        - For SC, ELF, and PARCHG jobs, POSCAR is deleted after generation
        - The relaxed structure (CONTCAR from Relax) will be copied by the SLURM script
        """
        job_dir = Path(job_dir)
        job_dir.mkdir(parents=True, exist_ok=True)
        
        if job_type == 'relax':
            # Symmetrize structure using PyXtal with progressive tolerance
            if PYXTAL_AVAILABLE:
                tolerances = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
                symmetrized = False
                for tol in tolerances:
                    try:
                        adaptor = AseAtomsAdaptor()
                        xtal = pyxtal()
                        xtal.from_seed(structure, tol=tol)
                        atoms = xtal.to_ase()
                        structure = adaptor.get_structure(atoms)
                        symmetrized = True
                        break
                    except Exception:
                        continue
                
                if not symmetrized:
                    print(f"    Warning: Could not symmetrize structure with tolerances {tolerances}")
                    print(f"    Proceeding with original structure...")
            
            vis = MPRelaxSet(structure, 
                user_incar_settings={
                    'PREC': 'Accurate',
                    'ALGO': 'Fast',
                    'ADDGRID': True,
                    'EDIFF': 1e-6,
                    'EDIFFG': -0.01,
                    'IBRION': 2,
                    'ISIF': 3,
                    'NELM': 60,
                    'NSW': 30,
                    'ISMEAR': 0,
                    'SIGMA': 0.05,
                    'ISPIN': 1,
                    'POTIM': 0.3,
                    'LWAVE': False,
                    'LCHARG': False,
                    'LAECHG': False,
                    'LASPH': False,
                    'LORBIT': 0,
                },
                user_kpoints_settings={'reciprocal_density': 64}
            )
        elif job_type == 'sc':
            vis = MPStaticSet(structure, 
                user_incar_settings={
                    'PREC': 'Accurate',
                    'ALGO': 'Normal',
                    'ADDGRID': True,
                    'EDIFF': 1e-6,
                    'IBRION': -1,
                    'NSW': 0,
                    'ISMEAR': 0,
                    'SIGMA': 0.05,
                    'ISPIN': 1,
                    'ISYM': 0,
                    'LWAVE': True,
                    'LCHARG': True,
                    'LAECHG': False,
                    'LASPH': False,
                    'LORBIT': 0,
                },
                user_kpoints_settings={'reciprocal_density': 64}
            )
        elif job_type == 'elf':
            vis = MPStaticSet(structure, 
                user_incar_settings={
                    'PREC': 'Accurate',
                    'ALGO': 'Normal',
                    'ADDGRID': True,
                    'EDIFF': 1e-6,
                    'IBRION': -1,
                    'NSW': 0,
                    'ISMEAR': -5,
                    'ISPIN': 1,
                    'ISYM': 0,
                    'LELF': True,
                    'ISTART': 1,
                    'ICHARG': 11,
                    'NEDOS': 1000,
                    'LWAVE': False,
                    'LCHARG': False,
                    'LAECHG': False,
                    'LASPH': False,
                    'LORBIT': 10,
                },
                user_kpoints_settings={'reciprocal_density': 64}
            )
        elif job_type == 'parchg':
            if not parchg_settings:
                raise ValueError("parchg_settings required for PARCHG job")
            
            base_settings = {
                'PREC': 'Accurate',
                'ALGO': 'Normal',
                'EDIFF': 1e-6,
                'IBRION': -1,
                'NSW': 0,
                'ISMEAR': 1,
                'SIGMA': 0.2,
                'ISPIN': 1,
                'ISYM': 0,
            }
            base_settings.update(parchg_settings)
            vis = MPStaticSet(structure, 
                user_incar_settings=base_settings,
                user_kpoints_settings={'reciprocal_density': 64}
            )
        else:
            raise ValueError(f"Unknown job_type: {job_type}")
        
        vis.write_input(job_dir)
        
        # For SC, ELF, and PARCHG, delete the generated POSCAR
        # We will copy CONTCAR from Relax in the SLURM script
        if job_type in ['sc', 'elf', 'parchg']:
            poscar_path = Path(job_dir) / 'POSCAR'
            if poscar_path.exists():
                poscar_path.unlink()
        
        return job_dir
    
    def create_slurm_script(self, job_dir, job_name, job_type='relax', prev_dir=None, relax_dir=None, parchg_label=None):
        """Create SLURM submission script."""
        job_dir = Path(job_dir).resolve()
        script_path = job_dir / 'job.sh'
        
        # Convert all paths to absolute paths
        if prev_dir:
            prev_dir = Path(prev_dir).resolve()
        if relax_dir:
            relax_dir = Path(relax_dir).resolve()
        
        script = f"""#!/bin/bash
#SBATCH --job-name={job_name}_{job_type}
#SBATCH --partition=Apus,Orion
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output={job_dir}/vasp_%j.out
#SBATCH --error={job_dir}/vasp_%j.err

# Load modules
module purge
module load intel/mkl/2024.0 intel/2024 intel-mpi/2021.11
ulimit -s unlimited

# Set environment
export OMP_NUM_THREADS=1
export PMG_VASP_PSP_DIR=$HOME/apps/PBE64

# Intel MPI settings for SLURM
if [ -e /opt/slurm/lib/libpmi.so ]; then
  export I_MPI_PMI_LIBRARY=/opt/slurm/lib/libpmi.so
else
  export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so.0
fi
export I_MPI_FABRICS=shm:ofi

# VASP executable (use srun for SLURM-native MPI launching)
VASP_CMD="srun $HOME/apps/vasp.6.2.1/bin/vasp_std"

# Change to job directory
cd {job_dir}

"""
        
        # For SC: copy CONTCAR from Relax (prev_dir is Relax)
        if prev_dir and job_type == 'sc':
            script += f"""
# Copy CONTCAR from relaxation as POSCAR
echo "Checking for relaxed structure..."
if [ ! -f "{prev_dir}/CONTCAR" ] || [ ! -s "{prev_dir}/CONTCAR" ]; then
    echo "ERROR: CONTCAR not found or empty in {prev_dir}"
    echo "Cannot proceed with SC calculation without relaxed structure"
    touch VASP_FAILED
    exit 1
fi

cp "{prev_dir}/CONTCAR" POSCAR
echo "Copied CONTCAR from {prev_dir} as POSCAR (size: $(du -h POSCAR | cut -f1))"

"""
        
        # For ELF: copy CONTCAR from Relax (relax_dir must be provided separately)
        if relax_dir and job_type == 'elf':
            script += f"""
# Copy CONTCAR from relaxation as POSCAR
echo "Checking for relaxed structure..."
if [ ! -f "{relax_dir}/CONTCAR" ] || [ ! -s "{relax_dir}/CONTCAR" ]; then
    echo "ERROR: CONTCAR not found or empty in {relax_dir}"
    echo "Cannot proceed with ELF calculation without relaxed structure"
    touch VASP_FAILED
    exit 1
fi

cp "{relax_dir}/CONTCAR" POSCAR
echo "Copied CONTCAR from {relax_dir} as POSCAR (size: $(du -h POSCAR | cut -f1))"

"""
        
        if prev_dir and job_type == 'elf':
            script += f"""
# Copy CHGCAR and WAVECAR from static calculation
echo "Checking for required files from SC calculation..."
if [ ! -f "{prev_dir}/CHGCAR" ] || [ ! -s "{prev_dir}/CHGCAR" ]; then
    echo "ERROR: CHGCAR not found or empty in {prev_dir}"
    echo "Cannot proceed with ELF calculation"
    touch VASP_FAILED
    exit 1
fi

if [ ! -f "{prev_dir}/WAVECAR" ] || [ ! -s "{prev_dir}/WAVECAR" ]; then
    echo "ERROR: WAVECAR not found or empty in {prev_dir}"
    echo "Cannot proceed with ELF calculation"
    touch VASP_FAILED
    exit 1
fi

cp "{prev_dir}/CHGCAR" .
echo "Copied CHGCAR from {prev_dir} (size: $(du -h {prev_dir}/CHGCAR | cut -f1))"

cp "{prev_dir}/WAVECAR" .
echo "Copied WAVECAR from {prev_dir} (size: $(du -h {prev_dir}/WAVECAR | cut -f1))"

# Verify copies succeeded
if [ ! -f "CHGCAR" ] || [ ! -s "CHGCAR" ]; then
    echo "ERROR: Failed to copy CHGCAR"
    touch VASP_FAILED
    exit 1
fi

if [ ! -f "WAVECAR" ] || [ ! -s "WAVECAR" ]; then
    echo "ERROR: Failed to copy WAVECAR"
    touch VASP_FAILED
    exit 1
fi

echo "Successfully verified CHGCAR and WAVECAR are present"

"""
        
        # For PARCHG: copy CONTCAR, CHGCAR, and WAVECAR from SC
        if job_type == 'parchg':
            script += f"""
# Copy CONTCAR from relaxation as POSCAR
echo "Checking for relaxed structure..."
if [ ! -f "{relax_dir}/CONTCAR" ] || [ ! -s "{relax_dir}/CONTCAR" ]; then
    echo "ERROR: CONTCAR not found or empty in {relax_dir}"
    touch VASP_FAILED
    exit 1
fi

cp "{relax_dir}/CONTCAR" POSCAR
echo "Copied CONTCAR from {relax_dir} as POSCAR"

# Copy CHGCAR and WAVECAR from SC calculation
echo "Checking for required files from SC calculation..."
if [ ! -f "{prev_dir}/CHGCAR" ] || [ ! -s "{prev_dir}/CHGCAR" ]; then
    echo "ERROR: CHGCAR not found or empty in {prev_dir}"
    touch VASP_FAILED
    exit 1
fi

if [ ! -f "{prev_dir}/WAVECAR" ] || [ ! -s "{prev_dir}/WAVECAR" ]; then
    echo "ERROR: WAVECAR not found or empty in {prev_dir}"
    touch VASP_FAILED
    exit 1
fi

cp "{prev_dir}/CHGCAR" .
cp "{prev_dir}/WAVECAR" .
echo "Copied CHGCAR and WAVECAR from {prev_dir}"

"""
        
        script += f"""
# Run VASP
echo "Starting VASP calculation: {job_type}"
echo "Working directory: $(pwd)"
echo "VASP command: $VASP_CMD"
echo "Start time: $(date)"

$VASP_CMD

EXIT_CODE=$?

echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"

# Check if successful
if [ $EXIT_CODE -eq 0 ]; then
"""
        
        # Add job-specific validation
        if job_type == 'sc':
            script += """
    # Verify critical files for SC calculation
    if [ -f "CHGCAR" ] && [ -s "CHGCAR" ] && [ -f "WAVECAR" ] && [ -s "WAVECAR" ]; then
        echo "VASP calculation completed successfully"
        echo "Verified CHGCAR and WAVECAR exist"
        
        # Keep CHGCAR and WAVECAR for downstream ELF/PARCHG jobs
        # They will be cleaned up after ELF calculation finishes
        rm -f CHG WFULL AECCAR* TMPCAR 2>/dev/null
        
        touch VASP_DONE
    else
        echo "VASP calculation failed: CHGCAR or WAVECAR missing/empty"
        [ ! -f "CHGCAR" ] && echo "  CHGCAR does not exist"
        [ -f "CHGCAR" ] && [ ! -s "CHGCAR" ] && echo "  CHGCAR is empty"
        [ ! -f "WAVECAR" ] && echo "  WAVECAR does not exist"
        [ -f "WAVECAR" ] && [ ! -s "WAVECAR" ] && echo "  WAVECAR is empty"
        touch VASP_FAILED
    fi
"""
        elif job_type == 'elf':
            # Generate relative paths to SC and PARCHG directories for cleanup
            elf_to_sc = os.path.relpath(prev_dir, job_dir) if prev_dir else None
            elf_to_parchg = "../PARCHG"
            
            script += f"""
    # Verify critical files for ELF calculation
    if [ -f "ELFCAR" ] && [ -s "ELFCAR" ]; then
        echo "VASP calculation completed successfully"
        echo "Verified ELFCAR exists"
        
        # ELF is typically the last job - clean up CHGCAR/WAVECAR from all directories
        # Clean up in current ELF directory
        rm -f CHGCAR CHG WAVECAR WFULL AECCAR* TMPCAR 2>/dev/null
"""
            
            if elf_to_sc:
                script += f"""        
        # Clean up in SC directory
        if [ -d "{elf_to_sc}" ]; then
            echo "  Cleaning {elf_to_sc}/"
            rm -f "{elf_to_sc}/CHGCAR" "{elf_to_sc}/WAVECAR" "{elf_to_sc}/CHG" "{elf_to_sc}/WFULL" 2>/dev/null
        fi
"""
            
            script += f"""        
        # Clean up in all PARCHG subdirectories
        if [ -d "{elf_to_parchg}" ]; then
            echo "  Cleaning {elf_to_parchg}/*/"
            find "{elf_to_parchg}" -maxdepth 2 -type f \\( -name "CHGCAR" -o -name "WAVECAR" -o -name "CHG" -o -name "WFULL" \\) -delete 2>/dev/null
        fi
        
        # Compress PARCHG files in ELF directory into a single tar.gz archive
        # Only compress if we have all 5 PARCHG files (band0, band1, e0025, e05, e10)
        echo "  Validating and compressing PARCHG files..."
        
        # Check if all 5 expected PARCHG files exist and are non-empty
        PARCHG_FILES="PARCHG-band0 PARCHG-band1 PARCHG-e0025 PARCHG-e05 PARCHG-e10"
        ALL_VALID=true
        
        for pf in $PARCHG_FILES; do
            if [ ! -f "$pf" ] || [ ! -s "$pf" ]; then
                echo "    Warning: $pf is missing or empty"
                ALL_VALID=false
            fi
        done
        
        if [ "$ALL_VALID" = true ]; then
            if [ ! -f "PARCHG.tar.gz" ]; then
                echo "    All 5 PARCHG files validated (non-empty)"
                if tar czf PARCHG.tar.gz $PARCHG_FILES 2>/dev/null; then
                    echo "    Created PARCHG.tar.gz in ELF directory"
                    
                    # Remove PARCHG-* files from ELF directory
                    rm -f PARCHG-*
                    echo "    Removed PARCHG-* files from ELF directory"
                    
                    # Remove PARCHG-* files from all PARCHG subdirectories
                    if [ -d "{elf_to_parchg}" ]; then
                        find "{elf_to_parchg}" -maxdepth 2 -type f -name "PARCHG-*" -delete 2>/dev/null
                        echo "    Removed PARCHG-* files from {elf_to_parchg}/*/ subdirectories"
                    fi
                    
                    echo "    Disk space saved: PARCHG files compressed into single archive"
                else
                    echo "    ERROR: Failed to create PARCHG.tar.gz"
                fi
            else
                echo "    PARCHG.tar.gz already exists"
            fi
        else
            echo "    ERROR: Not all PARCHG files are valid - skipping compression"
            echo "    This indicates PARCHG jobs may have failed"
        fi
        
        touch VASP_DONE
    else
        echo "VASP calculation failed: ELFCAR missing/empty"
        touch VASP_FAILED
    fi
"""
        elif job_type == 'relax':
            script += """
    # Verify critical files for Relax calculation
    if [ -f "CONTCAR" ] && [ -s "CONTCAR" ]; then
        echo "VASP calculation completed successfully"
        echo "Verified CONTCAR exists"
        
        # Clean up large unnecessary files to save disk space
        rm -f CHGCAR CHG WAVECAR WFULL AECCAR* TMPCAR 2>/dev/null
        
        touch VASP_DONE
    else
        echo "VASP calculation failed: CONTCAR missing/empty"
        touch VASP_FAILED
    fi
"""
        elif job_type == 'parchg':
            parchg_file = f"PARCHG-{parchg_label}" if parchg_label else "PARCHG"
            script += f"""
    # Verify PARCHG file and rename with label
    if [ -f "PARCHG" ] && [ -s "PARCHG" ]; then
        echo "VASP calculation completed successfully"
        echo "Verified PARCHG exists"
        mv PARCHG "{parchg_file}"
        echo "Renamed PARCHG to {parchg_file}"
        
        # Clean up large unnecessary files to save disk space
        rm -f CHGCAR CHG WAVECAR WFULL AECCAR* TMPCAR 2>/dev/null
        
        touch VASP_DONE
    else
        echo "VASP calculation failed: PARCHG missing/empty"
        touch VASP_FAILED
    fi
"""
        
        script += """
else
    echo "VASP calculation failed with exit code $EXIT_CODE"
    touch VASP_FAILED
fi
"""
        
        with open(script_path, 'w') as f:
            f.write(script)
        
        os.chmod(script_path, 0o755)
        return script_path
    
    def submit_job(self, script_path):
        """Submit a SLURM job and return job ID."""
        result = subprocess.run(
            ['sbatch', str(script_path)],
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            job_id = result.stdout.strip().split()[-1]
            return job_id
        else:
            raise RuntimeError(f"sbatch failed: {result.stderr}")
    
    def check_job_status(self, job_id):
        """Check SLURM job status. Returns: RUNNING, COMPLETED, FAILED, or NOTFOUND."""
        result = subprocess.run(
            ['squeue', '-j', job_id, '-h', '-o', '%T'],
            capture_output=True,
            text=True
        )
        
        if result.stdout.strip():
            slurm_state = result.stdout.strip()
            if slurm_state in ['RUNNING', 'PENDING', 'CONFIGURING']:
                return 'RUNNING'
            else:
                return 'RUNNING'  # Other states still in queue
        else:
            return 'NOTFOUND'
    
    def check_local_status(self, job_dir):
        """Check local directory for completion markers."""
        job_dir = Path(job_dir)
        if (job_dir / 'VASP_DONE').exists():
            return 'DONE'
        elif (job_dir / 'VASP_FAILED').exists():
            return 'FAILED'
        else:
            return 'UNKNOWN'
    
    def submit_relax(self, struct_id, structure):
        """Submit relaxation job for a structure."""
        sdata = self.db.get_structure(struct_id)
        if not sdata:
            return False
        
        relax_dir = Path(sdata['relax_dir'])
        job_name = struct_id
        
        print(f"  Submitting Relax: {struct_id}")
        
        try:
            self.create_vasp_inputs(structure, relax_dir, 'relax')
            script = self.create_slurm_script(relax_dir, job_name, 'relax')
            job_id = self.submit_job(script)
            
            self.db.update_state(struct_id, 'RELAX_RUNNING', relax_job_id=job_id)
            print(f"    Relax job ID: {job_id}")
            return True
        except Exception as e:
            print(f"    Error: {e}")
            self.db.update_state(struct_id, 'RELAX_FAILED', error=str(e))
            return False
    
    def submit_sc(self, struct_id, structure):
        """Submit SC job for a structure."""
        sdata = self.db.get_structure(struct_id)
        if not sdata:
            return False
        
        sc_dir = Path(sdata['sc_dir'])
        relax_dir = Path(sdata['relax_dir'])
        job_name = struct_id
        
        print(f"  Submitting SC: {struct_id}")
        
        try:
            self.create_vasp_inputs(structure, sc_dir, 'sc')
            script = self.create_slurm_script(sc_dir, job_name, 'sc', prev_dir=relax_dir)
            job_id = self.submit_job(script)
            
            self.db.update_state(struct_id, 'SC_RUNNING', sc_job_id=job_id)
            print(f"    SC job ID: {job_id}")
            return True
        except Exception as e:
            print(f"    Error: {e}")
            self.db.update_state(struct_id, 'SC_FAILED', error=str(e))
            return False
    
    def submit_elf(self, struct_id, structure):
        """Submit ELF job for a structure."""
        sdata = self.db.get_structure(struct_id)
        if not sdata:
            return False
        
        elf_dir = Path(sdata['elf_dir'])
        sc_dir = Path(sdata['sc_dir'])
        relax_dir = Path(sdata['relax_dir'])
        parchg_dir = Path(sdata['parchg_dir'])
        job_name = struct_id
        
        print(f"  Submitting ELF: {struct_id}")
        
        try:
            # Ensure ELF directory exists before copying PARCHG files
            elf_dir.mkdir(parents=True, exist_ok=True)
            
            # Copy PARCHG files from PARCHG/*/ subdirectories to ELF/ directory
            # This ensures Electride.py can find them for analysis
            parchg_labels = ['band0', 'band1', 'e0025', 'e05', 'e10']
            copied_files = []
            
            for label in parchg_labels:
                parchg_file = parchg_dir / label / f'PARCHG-{label}'
                if parchg_file.exists():
                    dest = elf_dir / f'PARCHG-{label}'
                    shutil.copy2(parchg_file, dest)
                    copied_files.append(label)
            
            if copied_files:
                print(f"    Copied PARCHG files to ELF dir: {', '.join(copied_files)}")
            else:
                print(f"    Warning: No PARCHG files found to copy (will use ELFCAR only)")
            
            self.create_vasp_inputs(structure, elf_dir, 'elf')
            script = self.create_slurm_script(elf_dir, job_name, 'elf', prev_dir=sc_dir, relax_dir=relax_dir)
            job_id = self.submit_job(script)
            
            self.db.update_state(struct_id, 'ELF_RUNNING', elf_job_id=job_id)
            print(f"    ELF job ID: {job_id}")
            return True
        except Exception as e:
            print(f"    Error: {e}")
            self.db.update_state(struct_id, 'ELF_FAILED', error=str(e))
            return False
    
    def submit_parchg(self, struct_id, structure):
        """
        Submit PARCHG jobs for semiconductors.
        Checks band gap from SC/vasprun.xml and submits band-specific PARCHG if needed.
        """
        sdata = self.db.get_structure(struct_id)
        if not sdata:
            return False
        
        sc_dir = Path(sdata['sc_dir'])
        parchg_dir = Path(sdata['parchg_dir'])
        relax_dir = Path(sdata['relax_dir'])
        
        vasprun_path = sc_dir / 'vasprun.xml'
        if not vasprun_path.exists():
            print(f"  {struct_id}: vasprun.xml not found in SC, skipping PARCHG")
            self.db.update_state(struct_id, 'PARCHG_SKIPPED', 
                               is_semiconductor=False, 
                               band_gap=None)
            return False
        
        print(f"  {struct_id}: Parsing band structure for PARCHG...")
        band_gap, is_semiconductor = parse_band_gap_from_vasprun(vasprun_path)
        
        self.db.update_state(struct_id, sdata['state'], 
                           band_gap=band_gap, 
                           is_semiconductor=is_semiconductor)
        
        if band_gap is None:
            print(f"    Warning: Could not parse band structure")
            print(f"    Will still attempt PARCHG (following HT-electride methodology)")
        else:
            material_type = "Semiconductor" if is_semiconductor else "Metal"
            print(f"    {material_type} (gap={band_gap:.4f} eV)")
        
        print(f"    Submitting PARCHG jobs (following HT-electride: PARCHG for ALL materials)...")
        
        try:
            vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=True)
            
            # Get highest occupied band (VBM) from NELECT
            nelect = int(vr.parameters['NELECT'])
            is_soc = vr.parameters.get('LSORBIT', False)
            fac = 1 if is_soc else 2  # SOC: 1 electron per band, no SOC: 2 electrons per band
            
            if nelect % 2 == 0:
                iband = int(nelect / fac)  # Highest occupied band (VBM)
            else:
                iband = int(nelect / fac) + 1  # Odd electrons
            
            grid = [vr.parameters['NGXF'], vr.parameters['NGYF'], vr.parameters['NGZF']]
            
            print(f"      NELECT={nelect}, VBM band index={iband}, NBANDS={vr.parameters['NBANDS']}")
            
            parchg_jobs = []
            
            # Band-specific PARCHG: VBM (band0) and VBM-1 (band1)
            for idx, band_offset in enumerate([0, -1]):
                band_idx = iband + band_offset
                label = f"band{idx}"
                job_dir = parchg_dir / label
                job_dir.mkdir(parents=True, exist_ok=True)
                
                parchg_settings = get_parchg_incar_settings_band(band_idx, grid)
                self.create_vasp_inputs(structure, job_dir, 'parchg', parchg_settings=parchg_settings)
                script = self.create_slurm_script(job_dir, f"{struct_id}_parchg", 'parchg', 
                                                prev_dir=sc_dir, relax_dir=relax_dir, 
                                                parchg_label=label)
                job_id = self.submit_job(script)
                parchg_jobs.append(job_id)
                print(f"      PARCHG {label}: job ID {job_id}")
            
            # Energy-window PARCHG: below Fermi level
            # EINT format: "E_min(relative to E_fermi) E_width"
            # Negative E_min means below Fermi level
            for e_val in [0.025, 0.5, 1.0]:
                label = f"e{str(e_val).replace('.', '')}"
                job_dir = parchg_dir / label
                job_dir.mkdir(parents=True, exist_ok=True)
                
                eint_str = f"{-e_val} 0.025"  # e.g., "-0.5 0.025" = 0.5 eV below E_fermi, width 25 meV
                parchg_settings = get_parchg_incar_settings_energy(eint_str, grid)
                self.create_vasp_inputs(structure, job_dir, 'parchg', parchg_settings=parchg_settings)
                script = self.create_slurm_script(job_dir, f"{struct_id}_parchg", 'parchg', 
                                                prev_dir=sc_dir, relax_dir=relax_dir, 
                                                parchg_label=label)
                job_id = self.submit_job(script)
                parchg_jobs.append(job_id)
                print(f"      PARCHG {label}: job ID {job_id}")
            
            self.db.update_state(struct_id, 'PARCHG_RUNNING', parchg_job_ids=parchg_jobs)
            return True
            
        except Exception as e:
            print(f"    Error submitting PARCHG: {e}")
            self.db.update_state(struct_id, 'PARCHG_FAILED', error=str(e))
            return False
    
    def update_structure_status(self, struct_id):
        """Check and update status of a structure."""
        sdata = self.db.get_structure(struct_id)
        if not sdata:
            return
        
        state = sdata['state']
        
        if state == 'RELAX_RUNNING':
            job_status = self.check_job_status(sdata['relax_job_id'])
            if job_status == 'NOTFOUND':
                local_status = self.check_local_status(sdata['relax_dir'])
                if local_status == 'DONE':
                    # Check convergence before marking as done
                    vasprun_path = Path(sdata['relax_dir']) / 'vasprun.xml'
                    if vasprun_path.exists():
                        try:
                            vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=False)
                            if not vr.converged_electronic:
                                self.db.update_state(struct_id, 'RELAX_FAILED', 
                                                   error='Electronic SCF not converged')
                                print(f"  {struct_id}: Relax FAILED (electronic not converged)")
                            else:
                                self.db.update_state(struct_id, 'RELAX_DONE')
                                print(f"  {struct_id}: Relax completed (electronic converged)")
                        except Exception as e:
                            self.db.update_state(struct_id, 'RELAX_FAILED', 
                                               error=f'Could not check convergence: {e}')
                            print(f"  {struct_id}: Relax FAILED (convergence check error)")
                    else:
                        self.db.update_state(struct_id, 'RELAX_FAILED', 
                                           error='vasprun.xml not found')
                        print(f"  {struct_id}: Relax FAILED (vasprun.xml missing)")
                elif local_status == 'FAILED':
                    self.db.update_state(struct_id, 'RELAX_FAILED')
                    print(f"  {struct_id}: Relax failed")
                else:
                    # Job not in queue and no completion marker - check if timed out
                    relax_dir = Path(sdata['relax_dir'])
                    err_files = list(relax_dir.glob('vasp_*.err'))
                    is_timeout = False
                    
                    if err_files:
                        err_file = max(err_files, key=lambda p: p.stat().st_mtime)
                        try:
                            with open(err_file, 'r') as f:
                                if 'DUE TO TIME LIMIT' in f.read():
                                    is_timeout = True
                        except Exception:
                            pass
                    
                    if is_timeout:
                        # Check if electronic converged (using OUTCAR) and CONTCAR exists
                        outcar_path = relax_dir / 'OUTCAR'
                        contcar_path = relax_dir / 'CONTCAR'
                        
                        if not outcar_path.exists():
                            self.db.update_state(struct_id, 'RELAX_FAILED',
                                               error='Job timed out, OUTCAR not found')
                            print(f"  {struct_id}: Relax FAILED (timeout, OUTCAR missing)")
                        elif not contcar_path.exists() or contcar_path.stat().st_size == 0:
                            self.db.update_state(struct_id, 'RELAX_FAILED',
                                               error='Job timed out, CONTCAR missing/empty')
                            print(f"  {struct_id}: Relax FAILED (timeout, CONTCAR missing)")
                        elif check_electronic_convergence_outcar(outcar_path):
                            self.db.update_state(struct_id, 'RELAX_TMOUT',
                                               error='Relaxation timed out but electronic converged')
                            print(f"  {struct_id}: Relax TMOUT (electronic converged, proceeding)")
                        else:
                            self.db.update_state(struct_id, 'RELAX_FAILED',
                                               error='Job timed out, electronic not converged')
                            print(f"  {struct_id}: Relax FAILED (timeout, electronic not converged)")
                    else:
                        self.db.update_state(struct_id, 'RELAX_FAILED', 
                                           error='Job terminated without completion marker (crash)')
                        print(f"  {struct_id}: Relax FAILED (crash)")
        
        elif state == 'SC_RUNNING':
            job_status = self.check_job_status(sdata['sc_job_id'])
            if job_status == 'NOTFOUND':
                local_status = self.check_local_status(sdata['sc_dir'])
                if local_status == 'DONE':
                    # Check convergence before marking as done
                    vasprun_path = Path(sdata['sc_dir']) / 'vasprun.xml'
                    if vasprun_path.exists():
                        try:
                            vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=False)
                            if not vr.converged_electronic:
                                self.db.update_state(struct_id, 'SC_FAILED', 
                                                   error='Electronic SCF not converged')
                                print(f"  {struct_id}: SC FAILED (electronic not converged)")
                            else:
                                self.db.update_state(struct_id, 'SC_DONE')
                                print(f"  {struct_id}: SC completed (electronic converged)")
                        except Exception as e:
                            self.db.update_state(struct_id, 'SC_FAILED', 
                                               error=f'Could not check convergence: {e}')
                            print(f"  {struct_id}: SC FAILED (convergence check error)")
                    else:
                        self.db.update_state(struct_id, 'SC_FAILED', 
                                           error='vasprun.xml not found')
                        print(f"  {struct_id}: SC FAILED (vasprun.xml missing)")
                elif local_status == 'FAILED':
                    self.db.update_state(struct_id, 'SC_FAILED')
                    print(f"  {struct_id}: SC failed")
                else:
                    # Job not in queue and no completion marker - likely timed out or crashed
                    self.db.update_state(struct_id, 'SC_FAILED', 
                                       error='Job terminated without completion marker (timeout or crash)')
                    print(f"  {struct_id}: SC FAILED (timeout or crash)")
        
        elif state == 'ELF_RUNNING':
            job_status = self.check_job_status(sdata['elf_job_id'])
            if job_status == 'NOTFOUND':
                local_status = self.check_local_status(sdata['elf_dir'])
                if local_status == 'DONE':
                    self.db.update_state(struct_id, 'ELF_DONE')
                    print(f"  {struct_id}: ELF completed")
                elif local_status == 'FAILED':
                    self.db.update_state(struct_id, 'ELF_FAILED')
                    print(f"  {struct_id}: ELF failed")
                else:
                    # Job not in queue and no completion marker - likely timed out or crashed
                    self.db.update_state(struct_id, 'ELF_FAILED', 
                                       error='Job terminated without completion marker (timeout or crash)')
                    print(f"  {struct_id}: ELF FAILED (timeout or crash)")
        
        elif state == 'PARCHG_RUNNING':
            all_done = True
            any_failed = False
            parchg_dir = Path(sdata['parchg_dir'])
            
            for job_id in sdata.get('parchg_job_ids', []):
                job_status = self.check_job_status(job_id)
                if job_status == 'RUNNING':
                    all_done = False
                    break
            
            if all_done:
                # Check local status for each PARCHG job (no convergence check needed)
                for subdir in ['band0', 'band1', 'e0025', 'e05', 'e10']:
                    job_dir = parchg_dir / subdir
                    if job_dir.exists():
                        local_status = self.check_local_status(job_dir)
                        if local_status == 'FAILED':
                            any_failed = True
                            break
                        elif local_status == 'UNKNOWN':
                            # Job terminated without completion marker (timeout or crash)
                            any_failed = True
                            print(f"    {subdir}: no completion marker (timeout or crash)")
                            break
                
                if any_failed:
                    self.db.update_state(struct_id, 'PARCHG_FAILED',
                                       error='One or more PARCHG jobs terminated without completion marker')
                    print(f"  {struct_id}: PARCHG failed")
                else:
                    self.db.update_state(struct_id, 'PARCHG_DONE')
                    print(f"  {struct_id}: PARCHG completed")
    
    def monitor_and_submit(self, structures_dict):
        """Main monitoring loop that checks status and submits new jobs."""
        print("\n" + "="*70)
        print("Starting workflow monitoring loop...")
        print(f"Max concurrent structures: {self.max_concurrent}")
        print(f"Check interval: {self.check_interval}s")
        print("="*70 + "\n")
        sys.stdout.flush()
        
        while True:
            print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Checking job status...")
            sys.stdout.flush()
            
            # Update status of all running jobs
            for struct_id in list(self.db.data['structures'].keys()):
                self.update_structure_status(struct_id)
            
            # Check if we can submit new jobs
            running_count = self.db.get_running_count()
            print(f"Currently running: {running_count}/{self.max_concurrent}")
            
            # Try to submit next stage for completed jobs
            for struct_id in list(self.db.data['structures'].keys()):
                if running_count >= self.max_concurrent:
                    break
                
                sdata = self.db.get_structure(struct_id)
                state = sdata['state']
                structure = structures_dict.get(struct_id)
                
                if not structure:
                    continue
                
                if state == 'PENDING':
                    if self.submit_relax(struct_id, structure):
                        running_count += 1
                
                elif state == 'RELAX_DONE' or state == 'RELAX_TMOUT':
                    if self.submit_sc(struct_id, structure):
                        running_count += 1
                
                elif state == 'SC_DONE':
                    # Always submit PARCHG (following HT-electride methodology)
                    if self.submit_parchg(struct_id, structure):
                        running_count += 1
                
                elif state == 'PARCHG_DONE' or state == 'PARCHG_SKIPPED':
                    if self.submit_elf(struct_id, structure):
                        running_count += 1
            
            # Print statistics
            stats = self.db.get_stats()
            print("\nStatistics:")
            for state, count in sorted(stats['states'].items()):
                print(f"  {state}: {count}")
            sys.stdout.flush()
            
            # Check if all done
            pending_count = len(self.db.get_by_state('PENDING'))
            if running_count == 0 and pending_count == 0:
                completed = len(self.db.get_by_state('ELF_DONE'))
                total = stats['total']
                failed_count = (len(self.db.get_by_state('ELF_FAILED')) + 
                               len(self.db.get_by_state('PARCHG_FAILED')) + 
                               len(self.db.get_by_state('RELAX_FAILED')) + 
                               len(self.db.get_by_state('SC_FAILED')))
                if completed + failed_count >= total:
                    print("\n" + "="*70)
                    print("All workflows completed!")
                    print("="*70)
                    sys.stdout.flush()
                    break
            
            print(f"\nSleeping for {self.check_interval}s...")
            sys.stdout.flush()
            time.sleep(self.check_interval)
    
    def initialize_structures(self, results_dir, output_dir, 
                             max_compositions=None, max_structures=5,
                             prescreen_results=None):
        """Scan results directory and initialize database."""
        results_dir = Path(results_dir)
        output_dir = Path(output_dir)
        
        print("="*70)
        print("Initializing VASP Workflow")
        print("="*70)
        print(f"Results directory: {results_dir}")
        print(f"Output directory: {output_dir}")
        print(f"Max concurrent: {self.max_concurrent}")
        print(f"Max compositions: {max_compositions or 'all'}")
        print(f"Max structures: {max_structures}")
        
        # Load pre-screening results if provided
        passed_structures = None
        if prescreen_results:
            prescreen_path = Path(prescreen_results)
            if prescreen_path.exists():
                print(f"Pre-screening results: {prescreen_path}")
                with open(prescreen_path, 'r') as f:
                    prescreen_data = json.load(f)
                
                passed_structures = set()
                for result in prescreen_data.get('results', []):
                    if result.get('passed_prescreening', False):
                        passed_structures.add(result['structure_id'])
                
                print(f"Structures passed pre-screening: {len(passed_structures)}")
                print(f"Energy threshold: {prescreen_data['summary']['hull_threshold']} eV/atom")
            else:
                print(f"Warning: Pre-screening file not found: {prescreen_path}")
                print("Will process all structures without filtering")
        else:
            print("No pre-screening filter (will process all structures)")
        
        print("="*70 + "\n")
        
        comp_dirs = sorted(results_dir.glob("*_structures"))
        if max_compositions:
            comp_dirs = comp_dirs[:max_compositions]
        
        structures_dict = {}
        
        for comp_dir in comp_dirs:
            comp_name = comp_dir.name.replace("_structures", "")
            zip_path = comp_dir / "generated_crystals_cif.zip"
            
            if not zip_path.exists():
                print(f"  Skipping {comp_name} (no ZIP file)")
                continue
            
            print(f"Scanning {comp_name}...")
            
            structures = self.read_structures_from_zip(zip_path, max_structures)
            if not structures:
                print(f"  No structures found")
                continue
            
            added_count = 0
            for idx, structure in enumerate(structures, 1):
                struct_id = f"{comp_name}_s{idx:03d}"
                
                # Skip if not in passed structures (when filtering is enabled)
                if passed_structures is not None and struct_id not in passed_structures:
                    continue
                
                structures_dict[struct_id] = structure
                
                if struct_id not in self.db.data['structures']:
                    elements = sorted([str(el) for el in structure.composition.elements])
                    chemsys = '-'.join(elements)
                    
                    self.db.add_structure(
                        struct_id, comp_name, idx,
                        output_dir / comp_name,
                        chemsys=chemsys
                    )
                    added_count += 1
            
            if added_count > 0:
                print(f"  Added {added_count} structures")
        
        # Load structures from database that aren't in structures_dict yet
        # This handles resume scenarios where structures exist in DB but weren't loaded from ZIP
        print("\nChecking database for additional structures...")
        loaded_from_contcar = 0
        skipped_count = 0
        
        for struct_id, sdata in self.db.data['structures'].items():
            if struct_id in structures_dict:
                continue  # Already loaded from ZIP
            
            # Try to load from Relax/CONTCAR for structures that have been processed
            if sdata['state'] not in ['PENDING', 'RELAX_RUNNING']:
                relax_dir = Path(sdata['relax_dir'])
                contcar_path = relax_dir / 'CONTCAR'
                
                if contcar_path.exists():
                    try:
                        structure = Structure.from_file(str(contcar_path))
                        structures_dict[struct_id] = structure
                        loaded_from_contcar += 1
                    except Exception as e:
                        print(f"  Warning: Could not load {struct_id} from CONTCAR: {e}")
                        skipped_count += 1
                else:
                    skipped_count += 1
        
        if loaded_from_contcar > 0:
            print(f"  Loaded {loaded_from_contcar} structures from CONTCAR files (resume)")
        if skipped_count > 0:
            print(f"  Skipped {skipped_count} structures (no CONTCAR available)")
        
        self.db.data['config'] = {
            'max_concurrent': self.max_concurrent,
            'results_dir': str(results_dir),
            'output_dir': str(output_dir),
            'max_structures': max_structures
        }
        self.db.save()
        
        print(f"\nTotal structures ready for workflow: {len(structures_dict)}")
        
        # Report structures in database but not in structures_dict
        missing_count = len(self.db.data['structures']) - len(structures_dict)
        if missing_count > 0:
            print(f"  Note: {missing_count} structures in database but not in structures_dict")
            print(f"        (These will be skipped during monitoring)")
        
        return structures_dict


def main():
    parser = argparse.ArgumentParser(
        description="VASP Workflow Manager - Batch submission with monitoring"
    )
    parser.add_argument(
        '--results-dir',
        type=str,
        required=True,
        help="MatterGen results directory"
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='/scratch/$USER/VASP_JOBS',
        help="Output directory for VASP jobs"
    )
    parser.add_argument(
        '--db',
        type=str,
        default='workflow.json',
        help="JSON database file path"
    )
    parser.add_argument(
        '--max-concurrent',
        type=int,
        default=10,
        help="Max concurrent structures running"
    )
    parser.add_argument(
        '--max-compositions',
        type=int,
        default=None,
        help="Max compositions to process"
    )
    parser.add_argument(
        '--max-structures',
        type=int,
        default=5,
        help="Max structures per composition"
    )
    parser.add_argument(
        '--check-interval',
        type=int,
        default=60,
        help="Status check interval in seconds"
    )
    parser.add_argument(
        '--init-only',
        action='store_true',
        help="Only initialize database, don't start monitoring"
    )
    parser.add_argument(
        '--prescreen-results',
        type=str,
        default=None,
        help="Path to prescreening_stability.json (filters structures by energy_above_hull)"
    )
    
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    db_path = Path(args.db).expanduser()
    
    if not db_path.is_absolute():
        db_path = output_dir / args.db
    
    # Create workflow manager
    manager = VASPWorkflowManager(
        db_path=db_path,
        max_concurrent=args.max_concurrent,
        check_interval=args.check_interval
    )
    
    # Initialize structures
    structures_dict = manager.initialize_structures(
        results_dir=results_dir,
        output_dir=output_dir,
        max_compositions=args.max_compositions,
        max_structures=args.max_structures,
        prescreen_results=args.prescreen_results
    )
    
    if args.init_only:
        print("\n" + "="*70)
        print("Initialization complete!")
        print(f"Database: {db_path}")
        print("="*70)
        return
    
    # Start monitoring and submission loop
    try:
        manager.monitor_and_submit(structures_dict)
    except KeyboardInterrupt:
        print("\n\nWorkflow interrupted by user.")
        print(f"Database saved to: {db_path}")
        print("Resume with same command to continue.")


if __name__ == '__main__':
    main()

