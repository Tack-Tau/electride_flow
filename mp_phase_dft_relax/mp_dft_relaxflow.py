#!/usr/bin/env python3
"""
VASP Workflow Manager for MP Phase DFT Relaxations

Manages VASP relaxation jobs for MP phases to validate DFT energies.
Simplified workflow - only performs relaxation (no SC/PARCHG/ELF).

Uses the same VASP settings as the main workflow for consistency.
"""

import os
import sys
import json
import time
import argparse
import subprocess
import warnings
from pathlib import Path
from datetime import datetime
from collections import defaultdict

from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet, BadInputSetWarning
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Suppress POTCAR warnings (they are informational only)
warnings.filterwarnings('ignore', category=BadInputSetWarning)

# MAGMOM overrides for elements with wrong pymatgen neutral-element defaults.
# pymatgen gives Co=0.6 (should be ~5 for magnetic intermetallics) and
# most rare earths get 0.6 (only Ce=5 and Eu=10 have neutral entries).
# RE values: Hund's rule spin-only moments (unpaired 4f electron count).
MAGMOM_OVERRIDE = {
    'Pr': 2.0, 'Nd': 3.0, 'Pm': 4.0, 'Sm': 5.0,
    'Gd': 7.0, 'Tb': 6.0, 'Dy': 5.0, 'Ho': 4.0,
    'Er': 3.0, 'Tm': 2.0, 'Yb': 1.0,
    'Co': 2.0,
}


def build_magmom(structure):
    """Build MAGMOM dict overriding wrong pymatgen defaults for intermetallics.
    
    When overriding via user_incar_settings, pymatgen does NOT merge with its
    config defaults, so all elements must be included in the returned dict.
    Returns dict (element -> value) if any override is needed,
    None otherwise (use pymatgen defaults).
    """
    elements = [str(el) for el in structure.composition.elements]
    if not any(el in MAGMOM_OVERRIDE for el in elements):
        return None
    pmg_defaults = {
        'Ce': 5, 'Eu': 10, 'Fe': 5, 'Ni': 5,
        'Mn': 5, 'Cr': 5, 'V': 5, 'Mo': 5, 'W': 5,
    }
    result = {}
    for el in elements:
        if el in MAGMOM_OVERRIDE:
            result[el] = MAGMOM_OVERRIDE[el]
        elif el in pmg_defaults:
            result[el] = pmg_defaults[el]
        else:
            result[el] = 0.6
    return result


def check_electronic_convergence_oszicar(relax_dir):
    """
    Check electronic convergence from OSZICAR and return energy.

    Searches backwards through ionic steps (F= lines) to find the most recent
    one with converged electronic SCF (iteration count < NELM). Returns the
    energy from that step for reliable timeout recovery.

    Args:
        relax_dir: Path to VASP relaxation directory (contains OSZICAR and INCAR)

    Returns:
        tuple: (converged: bool, total_energy: float or None)
    """
    relax_dir = Path(relax_dir)
    oszicar_path = relax_dir / 'OSZICAR'
    incar_path = relax_dir / 'INCAR'

    if not oszicar_path.exists():
        return False, None

    nelm = 60
    if incar_path.exists():
        try:
            with open(incar_path, 'r') as f:
                for line in f:
                    if 'NELM' in line and '=' in line:
                        val = line.split('=')[1].split()[0].strip()
                        nelm = int(val)
                        break
        except Exception:
            pass

    try:
        with open(oszicar_path, 'r') as f:
            lines = [l.rstrip() for l in f.readlines() if l.strip()]

        if len(lines) < 2:
            return False, None

        # Search backwards through F= lines for the most recent ionic step
        # with converged electronic SCF
        search_from = len(lines) - 1
        while search_from > 0:
            f_idx = None
            for i in range(search_from, -1, -1):
                if 'F=' in lines[i]:
                    f_idx = i
                    break

            if f_idx is None or f_idx < 1:
                return False, None

            # Extract total energy from F= line
            f_line = lines[f_idx]
            try:
                f_pos = f_line.index('F=')
                energy_str = f_line[f_pos + 2:].split()[0]
                total_energy = float(energy_str)
            except (ValueError, IndexError):
                search_from = f_idx - 1
                continue

            # Check electronic convergence: line before F= is the last
            # electronic SCF iteration for that ionic step
            scf_line = lines[f_idx - 1]
            parts = scf_line.split()
            if len(parts) >= 2:
                try:
                    e_step = int(parts[1])
                    if e_step < nelm:
                        return True, total_energy
                except ValueError:
                    pass

            # This step's electronic SCF didn't converge, try earlier step
            search_from = f_idx - 1

        return False, None

    except Exception:
        return False, None


def load_structures_from_cache(cache_dir):
    """
    Load structures from mp_cache_structs directory (POSCAR .vasp files).
    
    Returns:
        dict: {mp_id: {'structure': Structure, 'chemsys': str, 'formula': str}}
    """
    cache_dir = Path(cache_dir)
    structures = {}
    
    print("Loading structures from cache...")
    
    vasp_files = sorted(cache_dir.glob("mp-*.vasp"))
    
    if not vasp_files:
        print(f"  WARNING: No POSCAR (.vasp) files found in {cache_dir}")
        return structures
    
    print(f"  Found {len(vasp_files)} POSCAR files")
    
    for vasp_file in vasp_files:
        mp_id = vasp_file.stem  # mp-12345 from mp-12345.vasp
        
        try:
            structure = Structure.from_file(str(vasp_file))
            
            chemsys = '-'.join(sorted(set([el.symbol for el in structure.composition.elements])))
            
            structures[mp_id] = {
                'structure': structure,
                'chemsys': chemsys,
                'formula': structure.composition.reduced_formula
            }
            
        except Exception as e:
            print(f"  ERROR: Failed to load {mp_id}: {e}")
    
    print(f"Loaded {len(structures)} structures\n")
    return structures


def initialize_database(db_path, structures_dict):
    """
    Initialize or load workflow database.
    
    Database structure:
    {
        'structures': {
            'mp-12345': {
                'state': 'PENDING',
                'mp_id': 'mp-12345',
                'chemsys': 'B-Li-N',
                'formula': 'Li3BN2',
                'relax_dir': 'VASP_JOBS/mp_B-Li-N/mp-12345/Relax',
                'slurm_id': None,
                'submit_time': None,
                'update_time': None,
                'mp_energy_per_atom': -5.234,
                'vasp_energy_per_atom': None
            }
        }
    }
    """
    if db_path.exists():
        print(f"Loading existing database: {db_path}")
        with open(db_path, 'r') as f:
            db = json.load(f)
        print(f"  Found {len(db['structures'])} structures in database\n")
        return db
    
    print("Creating new database...")
    db = {'structures': {}}
    
    for mp_id, sdata in structures_dict.items():
        db['structures'][mp_id] = {
            'state': 'PENDING',
            'mp_id': mp_id,
            'chemsys': sdata['chemsys'],
            'formula': sdata['formula'],
            'relax_dir': None,
            'slurm_id': None,
            'submit_time': None,
            'update_time': None,
            'mp_energy_per_atom': None,  # Will be fetched by compare script
            'vasp_energy_per_atom': None
        }
    
    save_database(db_path, db)
    print(f"Initialized database with {len(db['structures'])} structures\n")
    
    return db


def save_database(db_path, db):
    """Save workflow database."""
    with open(db_path, 'w') as f:
        json.dump(db, f, indent=2)


def symmetrize_structure(structure, symprec_list=(0.01, 0.1)):
    """Symmetrize structure to clean up numerical noise from CIF conversion.

    Tries progressively looser tolerances until symmetry is detected.
    Prefers primitive cell, falls back to conventional if it matches original size.
    If neither matches, returns original structure with a WARNING.
    """
    n_orig = len(structure)
    for symprec in symprec_list:
        try:
            sga = SpacegroupAnalyzer(structure, symprec=symprec)
            prim = sga.get_primitive_standard_structure()
            spg = sga.get_space_group_symbol()
            spg_num = sga.get_space_group_number()
            if len(prim) == n_orig:
                print(f"    Symmetrized (prim): {spg} (#{spg_num}), "
                      f"symprec={symprec}, {len(prim)} atoms")
                return prim
            refined = sga.get_refined_structure()
            if len(refined) == n_orig:
                print(f"    Symmetrized (conv): {spg} (#{spg_num}), "
                      f"symprec={symprec}, {n_orig} atoms")
                return refined
            print(f"    WARNING: Symmetrization changes atom count "
                  f"({n_orig} -> prim {len(prim)}, conv {len(refined)}), "
                  f"symprec={symprec}, skipping")
        except Exception:
            continue
    return structure


def create_vasp_relax_inputs(structure, job_dir):
    """
    Create VASP input files for relaxation using same settings as bin_mag_flow.py.
    """
    job_dir = Path(job_dir)
    job_dir.mkdir(parents=True, exist_ok=True)

    structure = symmetrize_structure(structure)

    incar_settings = {
        'PREC': 'Normal',
        'ALGO': 'All',
        'ADDGRID': True,
        'EDIFF': 1e-4,
        'EDIFFG': -0.005,
        'IBRION': 1,
        'ISIF': 3,
        'NELM': 120,
        'NSW': 100,
        'ISMEAR': 1,
        'SIGMA': 0.05,
        'ISPIN': 2,
        'POTIM': 0.2,
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
        'LAECHG': False,
        'LASPH': True,
        'LORBIT': 11,
        'SYMPREC': 1e-4,
        'NCORE': 4,
    }
    
    # Override MAGMOM for RE and Co (pymatgen defaults are wrong for intermetallics)
    magmom = build_magmom(structure)
    if magmom is not None:
        incar_settings['MAGMOM'] = magmom
    
    vis = MPRelaxSet(structure,
        user_incar_settings=incar_settings,
        user_kpoints_settings={'reciprocal_density': 64}
    )
    
    vis.write_input(str(job_dir))


def create_slurm_script(job_dir, job_name):
    """Create SLURM submission script for VASP relaxation."""
    job_dir = Path(job_dir).resolve()
    
    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}_relax
#SBATCH --partition=Apus,Orion
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --exclude=str-c[85-97]
#SBATCH --output={job_dir}/vasp_%j.out
#SBATCH --error={job_dir}/vasp_%j.err

# Load modules
module purge
module load intel/mkl/2024.0 intel/2024 intel-mpi/2021.11
ulimit -s unlimited

# Set environment
export OMP_NUM_THREADS=1
export PMG_VASP_PSP_DIR=$HOME/apps/PBE52

# Intel MPI settings for SLURM
if [ -e /opt/slurm/lib/libpmi.so ]; then
  export I_MPI_PMI_LIBRARY=/opt/slurm/lib/libpmi.so
else
  export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so.0
fi
export I_MPI_FABRICS=shm:ofi

# VASP executable (use srun for SLURM-native MPI launching)
VASP_CMD="srun --mpi=pmi2 $HOME/apps/vasp.6.2.1/bin/vasp_std"

# Change to job directory
cd {job_dir}

# Run VASP
echo "Starting VASP coarse relaxation"
echo "Working directory: $(pwd)"
echo "VASP command: $VASP_CMD"
echo "Start time: $(date)"

$VASP_CMD

EXIT_CODE=$?

echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"

# Check if successful
if [ $EXIT_CODE -eq 0 ]; then
    # Verify critical files for Relax calculation
    if [ -f "CONTCAR" ] && [ -s "CONTCAR" ]; then
        echo "VASP calculation completed successfully"
        echo "Verified CONTCAR exists"
        
        # Clean up large unnecessary files to save disk space
        rm -f CHGCAR CHG WAVECAR WFULL AECCAR* TMPCAR PROCAR 2>/dev/null
        
        touch VASP_DONE
    else
        echo "VASP calculation failed: CONTCAR missing/empty"
        # Clean up large intermediate files to save disk space
        rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR PROCAR 2>/dev/null
        touch VASP_FAILED
    fi
else
    echo "VASP calculation failed with exit code $EXIT_CODE"
    # Clean up large intermediate files to save disk space
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR PROCAR 2>/dev/null
    touch VASP_FAILED
fi
"""
    
    script_path = job_dir / 'job.sh'
    with open(script_path, 'w') as f:
        f.write(script)
    
    script_path.chmod(0o755)
    return script_path


def submit_vasp_job(script_path):
    """
    Submit VASP job via SLURM.
    
    Returns:
        str: SLURM job ID or None if failed
    """
    try:
        result = subprocess.run(
            ['sbatch', str(script_path)],
            capture_output=True,
            text=True,
            check=True
        )
        
        # Parse job ID from "Submitted batch job 12345"
        output = result.stdout.strip()
        if 'Submitted batch job' in output:
            job_id = output.split()[-1]
            return job_id
        
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: SLURM submission failed: {e}")
    
    return None


def check_job_status(slurm_id):
    """
    Check SLURM job status.
    
    Returns:
        str: 'RUNNING'/'PENDING' (in queue), 'RELAX_DONE' (sacct COMPLETED),
             'TIMEOUT' (sacct TIMEOUT/CANCELLED), 'RELAX_FAILED' (other sacct),
             or None (sacct unavailable)
    """
    if slurm_id is None:
        return None
    
    try:
        result = subprocess.run(
            ['squeue', '-j', slurm_id, '-h', '-o', '%T'],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if result.returncode == 0 and result.stdout.strip():
            status = result.stdout.strip()
            return status
        
        # Job not in queue - check if completed or failed
        result = subprocess.run(
            ['sacct', '-j', slurm_id, '-n', '-o', 'State'],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if result.returncode == 0 and result.stdout.strip():
            status = result.stdout.strip().split()[0]
            if 'COMPLETED' in status:
                return 'RELAX_DONE'
            elif 'TIMEOUT' in status or 'CANCELLED' in status:
                return 'TIMEOUT'
            else:
                return 'RELAX_FAILED'
        
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        pass
    
    return None


def check_relax_convergence(relax_dir):
    """
    Check if VASP relaxation converged successfully.
    
    Returns:
        tuple: (converged: bool, energy_per_atom: float or None)
    """
    vasprun_path = Path(relax_dir) / 'vasprun.xml'
    
    if not vasprun_path.exists():
        return False, None
    
    try:
        vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=False)
        
        if not vr.converged_electronic:
            return False, None
        
        final_energy = vr.final_energy
        n_atoms = len(vr.final_structure)
        energy_per_atom = final_energy / n_atoms
        
        return True, energy_per_atom
        
    except Exception as e:
        print(f"  WARNING: Could not parse vasprun.xml: {e}")
        return False, None


def submit_pending_jobs(db, structures_dict, vasp_jobs_dir, max_concurrent):
    """
    Submit pending jobs up to max_concurrent limit.
    
    Returns:
        int: Number of jobs submitted
    """
    running_count = sum(1 for s in db['structures'].values() if s['state'] == 'RELAX_RUNNING')
    
    if running_count >= max_concurrent:
        return 0
    
    pending = [mp_id for mp_id, s in db['structures'].items() if s['state'] == 'PENDING']
    
    if not pending:
        return 0
    
    submitted = 0
    for mp_id in pending:
        if running_count + submitted >= max_concurrent:
            break
        
        if mp_id not in structures_dict:
            print(f"  WARNING: {mp_id} not in loaded structures, skipping...")
            db['structures'][mp_id]['state'] = 'SKIPPED'
            continue
        
        sdata = db['structures'][mp_id]
        structure = structures_dict[mp_id]['structure']
        chemsys = sdata['chemsys']
        
        relax_dir = Path(vasp_jobs_dir) / f"mp_{chemsys}" / mp_id / "Relax"
        
        print(f"Submitting {mp_id} ({sdata['formula']})...")
        
        try:
            create_vasp_relax_inputs(structure, relax_dir)
            script_path = create_slurm_script(relax_dir, f"{mp_id}")
            slurm_id = submit_vasp_job(script_path)
            
            if slurm_id:
                sdata['state'] = 'RELAX_RUNNING'
                sdata['slurm_id'] = slurm_id
                sdata['relax_dir'] = str(relax_dir)
                sdata['submit_time'] = datetime.now().isoformat()
                sdata['update_time'] = datetime.now().isoformat()
                
                print(f"    Submitted as job {slurm_id}")
                submitted += 1
            else:
                print(f"    Submission failed")
                
        except Exception as e:
            print(f"  ERROR: {e}")
    
    return submitted


def update_job_status(db):
    """
    Update status of running jobs.
    
    Returns:
        tuple: (completed_count, failed_count)
    """
    completed = 0
    failed = 0
    
    running_structures = [(mp_id, s) for mp_id, s in db['structures'].items() if s['state'] == 'RELAX_RUNNING']
    
    for mp_id, sdata in running_structures:
        slurm_id = sdata['slurm_id']
        slurm_status = check_job_status(slurm_id)
        
        if slurm_status == 'RELAX_DONE':
            # Check local marker first, then vasprun.xml convergence
            relax_dir = Path(sdata['relax_dir'])
            
            if (relax_dir / 'VASP_FAILED').exists():
                sdata['state'] = 'RELAX_FAILED'
                sdata['update_time'] = datetime.now().isoformat()
                print(f"    {mp_id}: RELAX_FAILED (VASP_FAILED marker)")
                failed += 1
                continue
            
            converged, energy_per_atom = check_relax_convergence(sdata['relax_dir'])
            
            if converged:
                sdata['state'] = 'RELAX_DONE'
                sdata['vasp_energy_per_atom'] = energy_per_atom
                sdata['update_time'] = datetime.now().isoformat()
                print(f"    {mp_id}: RELAX_DONE (E={energy_per_atom:.6f} eV/atom)")
                completed += 1
            else:
                sdata['state'] = 'RELAX_FAILED'
                sdata['update_time'] = datetime.now().isoformat()
                print(f"    {mp_id}: RELAX_FAILED (electronic not converged)")
                failed += 1
        
        elif slurm_status in ('RELAX_FAILED', 'TIMEOUT', None):
            relax_dir = Path(sdata['relax_dir'])
            
            if slurm_status is None:
                # Check local markers when sacct has no record
                if (relax_dir / 'VASP_DONE').exists():
                    converged, energy_per_atom = check_relax_convergence(sdata['relax_dir'])
                    if converged:
                        sdata['state'] = 'RELAX_DONE'
                        sdata['vasp_energy_per_atom'] = energy_per_atom
                        sdata['update_time'] = datetime.now().isoformat()
                        print(f"    {mp_id}: RELAX_DONE (E={energy_per_atom:.6f} eV/atom)")
                        completed += 1
                    else:
                        sdata['state'] = 'RELAX_FAILED'
                        sdata['update_time'] = datetime.now().isoformat()
                        print(f"    {mp_id}: RELAX_FAILED (VASP_DONE but vasprun not converged)")
                        failed += 1
                    continue
                
                if (relax_dir / 'VASP_FAILED').exists():
                    sdata['state'] = 'RELAX_FAILED'
                    sdata['update_time'] = datetime.now().isoformat()
                    print(f"    {mp_id}: RELAX_FAILED (VASP_FAILED marker)")
                    failed += 1
                    continue
            
            # Check for timeout: sacct may report TIMEOUT, CANCELLED, or
            # even FAILED depending on SLURM version/config -- always check
            # the .err file as a fallback.
            is_timeout = (slurm_status == 'TIMEOUT')
            
            if not is_timeout:
                err_files = list(relax_dir.glob('vasp_*.err'))
                if err_files:
                    err_file = max(err_files, key=lambda p: p.stat().st_mtime)
                    try:
                        with open(err_file, 'r') as f:
                            if 'DUE TO TIME LIMIT' in f.read():
                                is_timeout = True
                    except Exception:
                        pass
            
            if is_timeout:
                contcar_path = relax_dir / 'CONTCAR'
                
                if not contcar_path.exists() or contcar_path.stat().st_size == 0:
                    sdata['state'] = 'RELAX_FAILED'
                    sdata['update_time'] = datetime.now().isoformat()
                    print(f"    {mp_id}: RELAX_FAILED (timeout, CONTCAR missing/empty)")
                    failed += 1
                    continue

                converged, total_energy = check_electronic_convergence_oszicar(relax_dir)
                if converged and total_energy is not None:
                    try:
                        structure = Structure.from_file(str(contcar_path))
                        n_atoms = len(structure)
                        energy_per_atom = total_energy / n_atoms
                        
                        sdata['state'] = 'RELAX_TMOUT'
                        sdata['vasp_energy_per_atom'] = energy_per_atom
                        sdata['update_time'] = datetime.now().isoformat()
                        print(f"    {mp_id}: RELAX_TMOUT (E={energy_per_atom:.6f} eV/atom)")
                        completed += 1
                    except Exception as e:
                        sdata['state'] = 'RELAX_FAILED'
                        sdata['update_time'] = datetime.now().isoformat()
                        print(f"    {mp_id}: RELAX_FAILED (timeout, error reading CONTCAR: {e})")
                        failed += 1
                else:
                    sdata['state'] = 'RELAX_FAILED'
                    sdata['update_time'] = datetime.now().isoformat()
                    print(f"    {mp_id}: RELAX_FAILED (timeout, electronic not converged)")
                    failed += 1
            else:
                sdata['state'] = 'RELAX_FAILED'
                sdata['update_time'] = datetime.now().isoformat()
                print(f"    {mp_id}: RELAX_FAILED (SLURM: {slurm_status or 'unknown'})")
                failed += 1
    
    return completed, failed


def print_status_summary(db):
    """Print workflow status summary."""
    state_counts = defaultdict(int)
    for sdata in db['structures'].values():
        state_counts[sdata['state']] += 1
    
    print("\n" + "="*70)
    print("Workflow Status")
    print("="*70)
    for state in ['PENDING', 'RELAX_RUNNING', 'RELAX_DONE', 'RELAX_TMOUT', 'RELAX_FAILED', 'SKIPPED']:
        count = state_counts.get(state, 0)
        if count > 0 or state in ['PENDING', 'RELAX_RUNNING', 'RELAX_DONE', 'RELAX_FAILED']:
            print(f"  {state:14s}: {count:4d}")
    print("="*70 + "\n")
    sys.stdout.flush()


def main():
    parser = argparse.ArgumentParser(
        description="VASP Workflow Manager for MP Phase DFT Relaxations"
    )
    parser.add_argument(
        '--cache-dir',
        type=str,
        default='./mp_cache_structs',
        help="Directory with cached MP structures"
    )
    parser.add_argument(
        '--vasp-jobs-dir',
        type=str,
        default='./VASP_JOBS',
        help="Directory for VASP jobs"
    )
    parser.add_argument(
        '--db',
        type=str,
        default='./mp_relax_workflow.json',
        help="Workflow database file"
    )
    parser.add_argument(
        '--max-concurrent',
        type=int,
        default=10,
        help="Maximum concurrent VASP jobs"
    )
    parser.add_argument(
        '--check-interval',
        type=int,
        default=300,
        help="Job status check interval (seconds)"
    )
    parser.add_argument(
        '--chemsys',
        type=str,
        default=None,
        help="Process specific chemical system only (e.g., 'B-Li-N')"
    )
    
    args = parser.parse_args()
    
    cache_dir = Path(args.cache_dir).expanduser()
    vasp_jobs_dir = Path(args.vasp_jobs_dir).expanduser()
    db_path = Path(args.db).expanduser()
    
    vasp_jobs_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*70)
    print("MP Phase DFT Relaxation Workflow")
    print("="*70)
    print(f"Cache directory: {cache_dir}")
    print(f"VASP jobs directory: {vasp_jobs_dir}")
    print(f"Database: {db_path}")
    print(f"Max concurrent jobs: {args.max_concurrent}")
    print(f"Check interval: {args.check_interval}s")
    if args.chemsys:
        print(f"Chemical system filter: {args.chemsys}")
    print("="*70 + "\n")
    
    # Load structures
    structures_dict = load_structures_from_cache(cache_dir)
    
    if not structures_dict:
        print("ERROR: No structures found!")
        return 1
    
    # Filter by chemical system if requested
    if args.chemsys:
        structures_dict = {
            mp_id: sdata for mp_id, sdata in structures_dict.items()
            if sdata['chemsys'] == args.chemsys
        }
        print(f"Filtered to {len(structures_dict)} structures in {args.chemsys}\n")
    
    # Initialize database
    db = initialize_database(db_path, structures_dict)
    
    # If resuming with --chemsys filter, also filter the database
    if args.chemsys and db_path.exists():
        original_count = len(db['structures'])
        db['structures'] = {
            mp_id: sdata for mp_id, sdata in db['structures'].items()
            if sdata['chemsys'] == args.chemsys
        }
        filtered_count = len(db['structures'])
        if filtered_count < original_count:
            print(f"Filtered database to {filtered_count} structures in {args.chemsys}")
            print(f"  (removed {original_count - filtered_count} structures from other systems)\n")
    
    # Main workflow loop
    print("Starting workflow monitoring...\n")
    
    try:
        while True:
            # Update status of running jobs
            completed, failed = update_job_status(db)
            
            # Submit new jobs
            submitted = submit_pending_jobs(db, structures_dict, vasp_jobs_dir, args.max_concurrent)
            
            # Save database
            if submitted > 0 or completed > 0 or failed > 0:
                save_database(db_path, db)
            
            # Print status
            print_status_summary(db)
            
            # Check if done
            pending = sum(1 for s in db['structures'].values() if s['state'] == 'PENDING')
            running = sum(1 for s in db['structures'].values() if s['state'] == 'RELAX_RUNNING')
            
            if pending == 0 and running == 0:
                print("\n" + "="*70)
                print("All jobs completed! Workflow terminating.")
                print("="*70 + "\n")
                sys.stdout.flush()
                break
            
            # Show active jobs (for monitoring)
            if running > 0:
                running_jobs = [f"{mp_id} (job {s['slurm_id']})" 
                               for mp_id, s in db['structures'].items() 
                               if s['state'] == 'RELAX_RUNNING']
                print(f"Active jobs: {', '.join(running_jobs[:5])}")
                if len(running_jobs) > 5:
                    print(f"  ... and {len(running_jobs) - 5} more")
            
            # Wait before next check
            print(f"Waiting {args.check_interval}s before next check...")
            sys.stdout.flush()
            time.sleep(args.check_interval)
            
    except KeyboardInterrupt:
        print("\n\nInterrupted by user. Saving database...")
        save_database(db_path, db)
        print("Database saved. Exiting.")
        return 0
    
    # Final summary
    print("\n" + "="*70)
    print("Final Summary")
    print("="*70)
    
    done_count = sum(1 for s in db['structures'].values() if s['state'] == 'RELAX_DONE')
    tmout_count = sum(1 for s in db['structures'].values() if s['state'] == 'RELAX_TMOUT')
    failed_count = sum(1 for s in db['structures'].values() if s['state'] == 'RELAX_FAILED')
    
    print(f"Relaxed: {done_count}")
    print(f"Timed out (usable): {tmout_count}")
    print(f"Failed: {failed_count}")
    print(f"Total: {len(db['structures'])}")
    print("="*70 + "\n")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

