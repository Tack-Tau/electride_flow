#!/usr/bin/env python3
"""
MatterSim Pre-screening for VASPflow

Performs fast thermodynamic stability screening using MatterSim before expensive VASP calculations.
Outputs: prescreening_stability.json with structures passing energy_above_hull threshold.

IMPORTANT UPDATE (Modern MPRester):
- Uses mp_api.client.MPRester (modern API, not legacy pymatgen.ext.matproj)
- Queries ALL GGA/GGA+U entries via get_entries_in_chemsys() (not filtered by r2SCAN is_stable)
- Uses GGA-optimized structures (compatible with MatterSim relaxation)
- Matches methodology in compute_dft_e_hull.py and redo_mp_phase_relax.py
- Fixes false negative issue where r2SCAN geometries gave artificially high MatterSim energies
"""

import os
import sys
import json
import argparse
import zipfile
import warnings
import time
import random
import fcntl
import numpy as np
from pathlib import Path
from io import StringIO
from datetime import datetime

from pymatgen.core import Composition
from pymatgen.io.cif import CifParser
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
# Use modern MP API for correct GGA phase queries
try:
    from mp_api.client import MPRester
except ImportError:
    print("ERROR: mp-api package is required")
    print("Install with: pip install mp-api")
    sys.exit(1)
from pymatgen.io.ase import AseAtomsAdaptor

from ase.optimize import FIRE
from ase.filters import UnitCellFilter
from pyxtal import pyxtal
from pyxtal.db import database_topology

try:
    from mattersim.forcefield import MatterSimCalculator
    MATTERSIM_AVAILABLE = True
except ImportError:
    MATTERSIM_AVAILABLE = False
    print("ERROR: MatterSim not available!")
    sys.exit(1)

warnings.filterwarnings('ignore', category=UserWarning, message='.*POTCAR data with symbol.*')
warnings.filterwarnings('ignore', message='Using UFloat objects with std_dev==0')
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in det')


def validate_structure(pmg_struct):
    """
    Pre-validate structure to catch obviously malformed structures before MatterSim.
    Returns (is_valid, error_message).
    """
    try:
        # Check 1: Reasonable lattice parameters (not too small/large)
        lattice = pmg_struct.lattice
        abc = lattice.abc
        if any(a < 1.0 or a > 50.0 for a in abc):
            return False, f"Unrealistic lattice parameters: {abc}"
        
        # Check 2: Reasonable angles
        angles = lattice.angles
        if any(angle < 10.0 or angle > 170.0 for angle in angles):
            return False, f"Unrealistic angles: {angles}"
        
        # Check 3: Atoms too close (min distance check)
        dist_matrix = pmg_struct.distance_matrix
        min_dist = dist_matrix[dist_matrix > 0].min() if len(dist_matrix) > 1 else 1.0
        if min_dist < 0.5:  # Too close (< 0.5 Angstrom)
            return False, f"Atoms too close: min_dist={min_dist:.3f} Å"
        
        # Check 4: Try spglib symmetry analysis (catches spglib failures early)
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        try:
            sga = SpacegroupAnalyzer(pmg_struct, symprec=0.1)
            _ = sga.get_space_group_symbol()
        except Exception as e:
            # If spglib fails here, it will definitely fail later
            return False, f"Spglib symmetry detection failed: {str(e)[:100]}"
        
        return True, None
        
    except Exception as e:
        return False, f"Validation error: {str(e)[:100]}"


def relax_structure_mattersim(pmg_struct, calculator, structure_id=None, fmax=0.01, max_steps=500, logfile=None):
    """
    Relax structure using MatterSim + FIRE optimizer with symmetrization and cell relaxation.
    
    Args:
        pmg_struct: Pymatgen Structure object
        calculator: MatterSimCalculator instance (reused across batch)
        structure_id: Structure ID for error reporting (optional)
        fmax: Force convergence criterion (eV/Angstrom)
        max_steps: Maximum optimization steps
        logfile: None to suppress FIRE optimizer log, otherwise path to log file or set to '-' for stdout
    
    Returns:
        tuple: (relaxed_structure, energy_per_atom, error_message)
               Returns (None, None, error_msg) if relaxation fails
    
    Note:
        - Symmetrizes structure using PyXtal with progressive tolerance (1e-5 to 0.5)
        - Falls back to direct conversion if all tolerances fail
        - Uses UnitCellFilter to relax both cell and atomic positions
        - Calculator is reused for efficiency (batch processing)
        - Has extensive error handling for ASE calculator failures (NaN, inf, spglib crashes)
        - Never raises exceptions - returns None values on failure
    """
    
    sid_prefix = f"[{structure_id}] " if structure_id else ""
    
    # Symmetrize structure using PyXtal with progressive tolerance relaxation
    tolerances = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
    atoms = None
    
    for tol in tolerances:
        try:
            xtal = pyxtal()
            xtal.from_seed(pmg_struct, tol=tol)
            if len(xtal.check_short_distance(r=0.5)) > 0:
                continue
            atoms = xtal.to_ase()
            break
        except Exception:
            continue
    
    # If all tolerances failed, try direct conversion without symmetrization
    if atoms is None:
        try:
            adaptor = AseAtomsAdaptor()
            atoms = adaptor.get_atoms(pmg_struct)
        except Exception as e:
            return None, None, f"{sid_prefix}Structure conversion failed (PyXtal failed at all tolerances {tolerances}, direct conversion also failed): {e}"
    
    # Attach calculator (reused from batch)
    atoms.calc = calculator
    
    # Check initial energy calculation (catches NaN/inf early)
    try:
        initial_energy = atoms.get_potential_energy()
        if not np.isfinite(initial_energy):
            return None, None, f"{sid_prefix}Initial energy is not finite: {initial_energy}"
    except Exception as e:
        return None, None, f"{sid_prefix}Initial energy calculation failed: {e}"
    
    # Relax with error handling for ASE calculator failures
    try:
        ecf = UnitCellFilter(atoms)
        dyn = FIRE(ecf, a=0.1, logfile=logfile)
        dyn.run(fmax=fmax, steps=max_steps)
        
        # Get final energy
        energy = atoms.get_potential_energy()
        
        # Check for NaN/inf values
        if not np.isfinite(energy):
            # Clean up before returning
            try:
                del atoms, ecf, dyn
            except:
                pass
            return None, None, f"{sid_prefix}Final energy is not finite: {energy}"
        
        energy_per_atom = energy / len(atoms)
        
        # Convert back to pymatgen structure
        adaptor = AseAtomsAdaptor()
        relaxed_structure = adaptor.get_structure(atoms)
        
    except Exception as e:
        # Clean up before returning
        try:
            del atoms
            if 'ecf' in locals():
                del ecf
            if 'dyn' in locals():
                del dyn
        except:
            pass
        return None, None, f"{sid_prefix}Relaxation failed: {e}"
    
    # Clean up atoms object (but keep calculator for reuse)
    del atoms
    del ecf
    del dyn
    
    return relaxed_structure, energy_per_atom, None


def load_mp_cache_locked(cache_file):
    """
    Load MP cache with file locking for safe concurrent access.
    
    Returns:
        dict: Cache data keyed by (chemsys, entry_id)
    """
    cache_file = Path(cache_file)
    cached_data = {}
    
    if cache_file.exists():
        with open(cache_file, 'r') as f:
            # Acquire shared lock for reading
            fcntl.flock(f.fileno(), fcntl.LOCK_SH)
            try:
                cache_list = json.load(f)
                # Convert list to dict keyed by (chemsys, mp_id) for fast lookup
                for item in cache_list:
                    key = (item.get('chemsys', ''), item['entry_id'])
                    cached_data[key] = item
            finally:
                fcntl.flock(f.fileno(), fcntl.LOCK_UN)
    
    return cached_data


def save_mp_cache_locked(cache_file, new_entries_data):
    """
    Save MP cache entries with file locking for safe concurrent writes.
    
    Args:
        cache_file: Path to cache file
        new_entries_data: List of new cache entries to add
    """
    cache_file = Path(cache_file)
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Use a lock file to coordinate access
    lock_file = cache_file.with_suffix('.lock')
    
    with open(lock_file, 'w') as lock_f:
        # Acquire exclusive lock
        fcntl.flock(lock_f.fileno(), fcntl.LOCK_EX)
        try:
            # Load existing cache
            existing_cache = {}
            if cache_file.exists():
                with open(cache_file, 'r') as f:
                    cache_list = json.load(f)
                    for item in cache_list:
                        key = (item.get('chemsys', ''), item['entry_id'])
                        existing_cache[key] = item
            
            # Add new entries (avoid duplicates)
            for new_item in new_entries_data:
                key = (new_item['chemsys'], new_item['entry_id'])
                if key not in existing_cache:
                    existing_cache[key] = new_item
            
            # Write back to file
            all_cache_data = list(existing_cache.values())
            with open(cache_file, 'w') as f:
                json.dump(all_cache_data, f, indent=2)
            
            print(f"    Updated cache: {cache_file} ({len(all_cache_data)} total entries)")
        finally:
            fcntl.flock(lock_f.fileno(), fcntl.LOCK_UN)


def get_mp_stable_phases_mattersim(chemsys, mp_api_key, cache_file, calculator, batch_size=32):
    """
    Get MP GGA phases and relax with MatterSim for consistent energy reference.
    
    CRITICAL UPDATES (matches compute_dft_e_hull.py and redo_mp_phase_relax.py):
    - Uses modern mp_api.client.MPRester (not legacy pymatgen.ext.matproj)
    - Queries ALL GGA/GGA+U entries via get_entries_in_chemsys() (not filtered by r2SCAN is_stable)
    - Filters functionals: GGA/GGA+U accepted, r2SCAN/SCAN excluded
    - Deduplication: Prefers pure GGA over GGA+U for same formula
    - Uses GGA-optimized structures (compatible with MatterSim relaxation)
    
    This fixes the false negative issue where r2SCAN-stable + r2SCAN-geometry phases
    gave artificially high MatterSim energies due to geometry incompatibility.
    
    Key optimizations:
    - Uses single SHARED global cache file (mp_mattersim.json) across all parallel batches
    - Per-chemsys locking prevents duplicate MP queries across batches
    - Reuses calculator across batch for efficiency
    - Returns PDEntry objects for phase diagram analysis
    
    Args:
        chemsys: Chemical system string (e.g., 'B-Li-N')
        mp_api_key: Materials Project API key
        cache_file: Path to SHARED global cache file (mp_mattersim.json)
        calculator: MatterSimCalculator instance (reused for batch)
        batch_size: Batch size for processing (not used in this function but passed for consistency)
    
    Returns:
        List of PDEntry objects for GGA phases in this chemical system
    """
    cache_file = Path(cache_file)
    
    # Use per-chemical-system lock to prevent duplicate queries
    # This ensures only ONE batch queries MP for this chemsys at a time
    chemsys_lock_file = cache_file.parent / f"mp_cache_{chemsys.replace('-', '')}.lock"
    
    with open(chemsys_lock_file, 'w') as chemsys_lock:
        # Acquire exclusive lock for this chemical system
        fcntl.flock(chemsys_lock.fileno(), fcntl.LOCK_EX)
        
        try:
            # Double-check pattern: read cache AFTER acquiring lock
            # Another batch may have populated it while we waited for lock
            cached_data = load_mp_cache_locked(cache_file)
    
            # Extract entries for this chemsys from cache
            elements = chemsys.split('-')
            entries = []
            cached_count = 0
            
            for key, item in cached_data.items():
                item_chemsys = key[0]
                # Include if chemsys matches OR if it's a subset (e.g., 'B', 'Li' for 'B-Li-N')
                item_elements = sorted(item_chemsys.split('-'))
                if set(item_elements).issubset(set(elements)):
                    entry = PDEntry(
                        composition=Composition(item['composition']),
                        energy=item['energy'],
                        name=item['entry_id']
                    )
                    entries.append(entry)
                    cached_count += 1
            
            # Check if we have all required subsystems cached
            required_chemsys = set()
            for n in range(1, len(elements) + 1):
                from itertools import combinations
                for combo in combinations(elements, n):
                    required_chemsys.add('-'.join(sorted(combo)))
            
            cached_chemsys = set(item[0] for item in cached_data.keys())
            missing_chemsys = required_chemsys - cached_chemsys
            
            if cached_count > 0 and not missing_chemsys:
                print(f"    Using cached data: {cached_count} GGA phases for {chemsys}")
                return entries
    
            # Need to query MP for missing data
            print(f"    Querying MP for GGA phases in {chemsys}...")
            if missing_chemsys:
                print(f"      Missing subsystems: {sorted(missing_chemsys)}")
    
            try:
                with MPRester(mp_api_key) as mpr:
                    # Get ALL GGA entries in this chemical system (not filtered by r2SCAN stability)
                    # This ensures we get GGA data even for phases where r2SCAN is now preferred
                    computed_entries = mpr.get_entries_in_chemsys(
                        elements,
                        property_data=None,
                        conventional_unit_cell=False,
                        additional_criteria={}  # No is_stable filter - get all GGA/GGA+U entries
                    )
                    
                    print(f"    Retrieved {len(computed_entries)} entries from MP")
                    
                    # Filter for GGA/GGA+U, exclude r2SCAN/SCAN, and get structures
                    mp_phases = []
                    seen_formulas = {}
                    skipped_structure_retrieval = []
                    
                    for comp_entry in computed_entries:
                        entry_id = comp_entry.entry_id
                        
                        # Filter by functional based on entry_id suffix
                        has_U = False
                        has_SCAN = False
                        is_gga = False
                        
                        if '-GGA' in entry_id:
                            is_gga = True
                        elif '-GGA+U' in entry_id or '_GGA+U' in entry_id:
                            has_U = True
                            is_gga = True
                        elif '-r2SCAN' in entry_id or '_r2SCAN' in entry_id:
                            has_SCAN = True
                        elif '-SCAN' in entry_id or '_SCAN' in entry_id:
                            has_SCAN = True
                        else:
                            # No explicit suffix - likely legacy GGA
                            is_gga = True
                        
                        # Skip r2SCAN/SCAN
                        if has_SCAN or not is_gga:
                            continue
                        
                        # Get MP ID (base without functional suffix)
                        mp_id_parts = entry_id.split('-')
                        if len(mp_id_parts) >= 2:
                            mp_id = mp_id_parts[0] + '-' + mp_id_parts[1]
                        else:
                            mp_id = entry_id
                        
                        # Extract structure - try multiple methods
                        structure = None
                        
                        # Method 1: Try get_structure_by_material_id
                        try:
                            structure = mpr.get_structure_by_material_id(mp_id)
                        except Exception:
                            # Method 2: Try materials.summary.search
                            try:
                                summary_docs = mpr.materials.summary.search(
                                    material_ids=[mp_id],
                                    fields=["material_id", "structure"]
                                )
                                if summary_docs and len(summary_docs) > 0:
                                    structure = summary_docs[0].structure
                            except Exception:
                                skipped_structure_retrieval.append(f"{mp_id} ({comp_entry.composition.reduced_formula})")
                                continue
                        
                        if structure is None:
                            continue
                        
                        formula = structure.composition.reduced_formula
                        
                        # Handle duplicates: prefer pure GGA over GGA+U
                        if formula in seen_formulas:
                            existing_has_U = seen_formulas[formula]['has_U']
                            if not has_U and existing_has_U:
                                # Replace with pure GGA
                                mp_phases = [(mid, s, has_u) for mid, s, has_u in mp_phases 
                                            if s.composition.reduced_formula != formula]
                            else:
                                continue
                        
                        mp_phases.append((mp_id, structure, has_U))
                        seen_formulas[formula] = {'mp_id': mp_id, 'has_U': has_U}
                    
                    if skipped_structure_retrieval:
                        print(f"    WARNING: Could not retrieve structures for {len(skipped_structure_retrieval)} phases:")
                        for phase_info in skipped_structure_retrieval[:5]:
                            print(f"      {phase_info}")
                        if len(skipped_structure_retrieval) > 5:
                            print(f"      ... and {len(skipped_structure_retrieval) - 5} more")
                    
                    print(f"    Filtered to {len(mp_phases)} GGA/GGA+U phases with structures")
                    
                    # Verify we have terminal (elemental) phases
                    elements_found = set()
                    for mp_id, structure, has_U in mp_phases:
                        if len(structure.composition.elements) == 1:
                            elements_found.add(str(structure.composition.elements[0]))
                    
                    expected_elements = set(elements)
                    if elements_found != expected_elements:
                        missing = expected_elements - elements_found
                        print(f"    WARNING: Missing terminal phases for elements: {sorted(missing)}")
                        print(f"    Querying elemental phases directly with GGA filter...")
                        
                        for missing_elem in sorted(missing):
                            try:
                                print(f"      Querying {missing_elem} elemental phases...")
                                elem_docs = mpr.materials.summary.search(
                                    elements=[missing_elem],
                                    num_elements=(1, 1),
                                    fields=['material_id', 'structure', 'energy_per_atom']
                                )
                                
                                if not elem_docs:
                                    print(f"        ERROR: No phases found for {missing_elem}")
                                    continue
                                gga_docs = []
                                for doc in elem_docs:
                                    mp_id = doc.material_id
                                    if '-r2SCAN' in mp_id or '-SCAN' in mp_id or '_r2SCAN' in mp_id or '_SCAN' in mp_id:
                                        continue
                                    gga_docs.append(doc)
                                
                                if gga_docs:
                                    elem_doc = sorted(gga_docs, key=lambda d: d.energy_per_atom)[0]
                                    elem_mp_id = elem_doc.material_id
                                    elem_structure = elem_doc.structure
                                    print(f"        Found {elem_mp_id} (E={elem_doc.energy_per_atom:.4f} eV/atom)")
                                    mp_phases.append((elem_mp_id, elem_structure, False))
                                    elements_found.add(missing_elem)
                                else:
                                    print(f"        ERROR: No GGA phases found for {missing_elem}")
                            except Exception as e:
                                print(f"        ERROR querying {missing_elem}: {e}")
                    else:
                        print(f"      All terminal phases present: {sorted(elements_found)}")
                    
                    print(f"    Relaxing {len(mp_phases)} phases with MatterSim...")
                    
                    new_entries_data = []
                    success_count = 0
                    failed_count = 0
                    
                    for mp_id, structure, has_U in mp_phases:
                        # Determine chemsys for this phase
                        doc_elements = sorted([str(el) for el in structure.composition.elements])
                        doc_chemsys = '-'.join(doc_elements)
                        
                        # Skip if already cached
                        entry_id_cached = f"mp_mattersim_{mp_id}"
                        cache_key = (doc_chemsys, entry_id_cached)
                        if cache_key in cached_data:
                            continue
                        
                        # Relax GGA-optimized structure with MatterSim (using reused calculator)
                        relaxed_struct, energy_per_atom, error_msg = relax_structure_mattersim(
                            structure, calculator, structure_id=mp_id
                        )
                        
                        if relaxed_struct is None or energy_per_atom is None:
                            print(f"      Warning: Failed to relax {mp_id}: {error_msg}")
                            failed_count += 1
                            continue
                        
                        total_energy = energy_per_atom * relaxed_struct.composition.num_atoms
                        
                        # Create PDEntry
                        entry = PDEntry(
                            composition=relaxed_struct.composition,
                            energy=total_energy,
                            name=f"mp_mattersim_{mp_id}"
                        )
                        entries.append(entry)
                        
                        # Store for caching
                        new_entries_data.append({
                            'chemsys': doc_chemsys,
                            'composition': {str(el): float(amt) for el, amt in relaxed_struct.composition.items()},
                            'energy': total_energy,
                            'entry_id': f"mp_mattersim_{mp_id}",
                            'mp_id': mp_id
                        })
                        success_count += 1
                    
                    print(f"    Successfully relaxed {success_count} GGA phases")
                    if failed_count > 0:
                        print(f"    Failed to relax {failed_count} phases")
                    
                    # Update cache with new entries (with file locking)
                    if new_entries_data:
                        save_mp_cache_locked(cache_file, new_entries_data)
                    
                    return entries
                    
            except Exception as e:
                print(f"    Error querying MP for {chemsys}: {e}")
                import traceback
                traceback.print_exc()
                return entries
        
        finally:
            # Release the per-chemsys lock
            fcntl.flock(chemsys_lock.fileno(), fcntl.LOCK_UN)


def compute_energy_above_hull(structure, energy_per_atom, mp_entries):
    """
    Compute energy above hull for a structure using PDEntry.
    
    Args:
        structure: Pymatgen Structure
        energy_per_atom: Energy per atom (eV/atom)
        mp_entries: List of PDEntry objects for reference phases
    
    Returns:
        float: Energy above hull (eV/atom)
    """
    composition = structure.composition
    total_energy = energy_per_atom * composition.num_atoms
    
    entry = PDEntry(
        composition=composition,
        energy=total_energy,
        name='generated'
    )
    
    all_entries = list(mp_entries) + [entry]
    pd = PhaseDiagram(all_entries)
    decomp, e_above_hull = pd.get_decomp_and_e_above_hull(entry, allow_negative=True)
    
    return float(e_above_hull)


def read_structures_from_zip(zip_path, max_structures=None):
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


def main():
    parser = argparse.ArgumentParser(
        description="MatterSim Pre-screening for VASPflow"
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
        required=True,
        help="Output directory for results"
    )
    parser.add_argument(
        '--mp-api-key',
        type=str,
        default=None,
        help="Materials Project API key"
    )
    parser.add_argument(
        '--hull-threshold',
        type=float,
        default=0.1,
        help="Energy above hull threshold (eV/atom)"
    )
    parser.add_argument(
        '--device',
        type=str,
        default='cpu',
        choices=['cpu', 'cuda'],
        help="Device for MatterSim"
    )
    parser.add_argument(
        '--start-composition',
        type=int,
        default=0,
        help="Starting composition index for parallel batch processing (default: 0)"
    )
    parser.add_argument(
        '--max-compositions',
        type=int,
        default=None,
        help="Max compositions to process (count from start-composition)"
    )
    parser.add_argument(
        '--max-structures',
        type=int,
        default=0,
        help="Max structures per composition"
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=32,
        help="Batch size for structure processing (reuses MatterSimCalculator)"
    )
    parser.add_argument(
        '--mp-cache-dir',
        type=str,
        default=None,
        help="Directory to save MP cache file (default: same as output-dir)"
    )
    parser.add_argument(
        '--batch-id',
        type=str,
        default=None,
        help="Batch ID for parallel runs (creates unique checkpoint/output files)"
    )
    
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    
    mp_api_key = args.mp_api_key or os.environ.get('MP_API_KEY')
    if not mp_api_key:
        print("ERROR: Materials Project API key required!")
        print("  Use: --mp-api-key YOUR_KEY or export MP_API_KEY=YOUR_KEY")
        return 1
    
    # Use SHARED global cache file (with file locking for parallel safety)
    # All parallel batches share the same cache to avoid duplicate MP API calls
    mp_cache_dir = Path(args.mp_cache_dir).expanduser() if args.mp_cache_dir else output_dir
    mp_cache_file = mp_cache_dir / "mp_mattersim.json"
    
    print("="*70)
    print("MatterSim Pre-screening for VASPflow")
    print("="*70)
    print(f"Results directory: {results_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Hull threshold: {args.hull_threshold} eV/atom")
    print(f"Device: {args.device}")
    print(f"Batch size: {args.batch_size}")
    print(f"Start composition index: {args.start_composition}")
    print(f"Max compositions: {args.max_compositions or 'all'}")
    print(f"Max structures per composition: {args.max_structures}")
    print(f"MP cache file: {mp_cache_file}")
    print("="*70 + "\n")
    
    
    # Scan structures
    comp_dirs = sorted(results_dir.glob("*_structures"))
    
    # Apply start index for parallel batch processing
    if args.start_composition > 0:
        comp_dirs = comp_dirs[args.start_composition:]
    
    # Apply max compositions limit
    if args.max_compositions:
        comp_dirs = comp_dirs[:args.max_compositions]
    
    all_structures = []
    
    for comp_dir in comp_dirs:
        comp_name = comp_dir.name.replace("_structures", "")
        zip_path = comp_dir / "generated_crystals_cif.zip"
        
        if not zip_path.exists():
            continue
        
        print(f"Loading {comp_name}...")
        structures = read_structures_from_zip(zip_path, args.max_structures)
        
        for idx, structure in enumerate(structures, 1):
            struct_id = f"{comp_name}_s{idx:03d}"
            elements = sorted([str(el) for el in structure.composition.elements])
            chemsys = '-'.join(elements)
            
            all_structures.append({
                'id': struct_id,
                'composition': comp_name,
                'chemsys': chemsys,
                'structure': structure
            })
        
        print(f"  Added {len(structures)} structures")
    
    print(f"\nTotal structures loaded: {len(all_structures)}\n")
    
    # Pre-fetch MP stable phases for all chemical systems
    unique_chemsys = sorted(set(s['chemsys'] for s in all_structures))
    
    print("="*70)
    print("Pre-fetching MP GGA Phases (Modern MPRester)")
    print("="*70)
    print(f"Unique chemical systems: {len(unique_chemsys)}")
    print("Strategy: Query ALL GGA/GGA+U entries (not filtered by r2SCAN stability)")
    print("  This ensures GGA-optimized structures compatible with MatterSim")
    print("="*70 + "\n")
    
    # Create MatterSim calculator for MP phases (one-time)
    print(f"Creating MatterSimCalculator (device={args.device})...")
    try:
        import torch
        import gc
        calc_mp = MatterSimCalculator(
            load_path="MatterSim-v1.0.0-5M.pth",
            device=args.device
        )
        print("  Calculator created successfully\n")
    except Exception as e:
        print(f"ERROR: Failed to create MatterSimCalculator: {e}")
        return 1
    
    mp_entries_cache = {}
    for chemsys in unique_chemsys:
        print(f"Fetching MP stable phases for {chemsys}...")
        sys.stdout.flush()
        try:
            mp_entries = get_mp_stable_phases_mattersim(
                chemsys, mp_api_key, mp_cache_file, calc_mp, batch_size=args.batch_size
            )
            mp_entries_cache[chemsys] = mp_entries
            print(f"  → {len(mp_entries)} stable phases ready\n")
            sys.stdout.flush()
        except Exception as e:
            print(f"  → Error: {e}\n")
            sys.stdout.flush()
            mp_entries_cache[chemsys] = []
    
    # Clean up MP calculator
    del calc_mp
    if args.device == 'cuda':
        torch.cuda.empty_cache()
        gc.collect()
    
    print("="*70)
    print("MP stable phases pre-fetch complete")
    print("="*70 + "\n")
    
    # Run pre-screening with batch processing
    print("="*70)
    print("Running MatterSim Pre-screening (Batch Processing)")
    print("="*70)
    print(f"Processing {len(all_structures)} structures in batches of {args.batch_size}...\n")
    sys.stdout.flush()
    
    # Check for existing checkpoint
    batch_suffix = f"_batch{args.batch_id}" if args.batch_id else ""
    checkpoint_file = output_dir / f'prescreening_checkpoint{batch_suffix}.json'
    db_file = output_dir / f'prescreening_structures{batch_suffix}.db'
    
    if checkpoint_file.exists():
        print(f"Found checkpoint file, resuming from previous run...")
        with open(checkpoint_file, 'r') as f:
            checkpoint_data = json.load(f)
            results = checkpoint_data.get('results', [])
            passed = checkpoint_data.get('passed', 0)
            failed = checkpoint_data.get('failed', 0)
            processed_ids = {r['structure_id'] for r in results}
            currently_processing = checkpoint_data.get('currently_processing', None)
        
        print(f"Resuming: {len(results)}/{len(all_structures)} already processed")
        
        # Detect if we crashed during a structure (segfault detection)
        if currently_processing:
            print(f"\n   SEGFAULT DETECTED: Crashed while processing '{currently_processing}'")
            
            # Find the problematic structure
            problem_struct = None
            for s in all_structures:
                if s['id'] == currently_processing:
                    problem_struct = s
                    break
            
            if problem_struct:
                results.append({
                    'structure_id': currently_processing,
                    'composition': problem_struct['composition'],
                    'chemsys': problem_struct['chemsys'],
                    'mattersim_energy_per_atom': None,
                    'energy_above_hull': None,
                    'is_stable': False,
                    'passed_prescreening': False,
                    'error': 'Auto-skipped due to segfault (spglib crash)'
                })
                processed_ids.add(currently_processing)
                failed += 1
                print(f"     Added '{currently_processing}' to failed list\n")
            else:
                print(f"      Could not find structure in loaded data\n")
        
        # Open existing database for incremental updates
        if db_file.exists():
            print(f"Found existing database: {db_file}")
            db = database_topology(str(db_file))
        else:
            print(f"Creating new database: {db_file}")
            db = database_topology(str(db_file))
        print()
        sys.stdout.flush()
    else:
        results = []
        passed = 0
        failed = 0
        processed_ids = set()
        # Create new database for incremental updates
        db = database_topology(str(db_file))
        print(f"Created new database: {db_file}\n")
    
    # Filter unprocessed structures
    structures_to_process = [s for s in all_structures if s['id'] not in processed_ids]
    
    if not structures_to_process:
        print("All structures already processed!\n")
    else:
        print(f"Processing {len(structures_to_process)} remaining structures...\n")
        
        # Process in batches
        import torch
        import gc
        from tqdm import tqdm
        
        n_batches = (len(structures_to_process) + args.batch_size - 1) // args.batch_size
        
        for batch_idx in range(n_batches):
            batch_start = batch_idx * args.batch_size
            batch_end = min((batch_idx + 1) * args.batch_size, len(structures_to_process))
            batch = structures_to_process[batch_start:batch_end]
            
            print(f"\n{'='*70}")
            print(f"Batch {batch_idx + 1}/{n_batches}: Processing {len(batch)} structures")
            print(f"{'='*70}\n")
            
            # Mark the first structure in this batch as currently processing
            # This helps detect segfaults during calculator creation
            if batch:
                first_struct_id = batch[0]['id']
                checkpoint_data_temp = {
                    'results': [{k: v for k, v in r.items() if k != 'pmg'} for r in results],
                    'passed': passed,
                    'failed': failed,
                    'currently_processing': first_struct_id,
                    'timestamp': datetime.now().isoformat()
                }
                with open(checkpoint_file, 'w') as f:
                    json.dump(checkpoint_data_temp, f, indent=2)
            
            # Create calculator for this batch (with error handling)
            try:
                calc = MatterSimCalculator(
                    load_path="MatterSim-v1.0.0-5M.pth",
                    device=args.device
                )
            except Exception as e:
                print(f"ERROR creating calculator for batch {batch_idx + 1}: {e}")
                print(f"Skipping entire batch of {len(batch)} structures\n")
                sys.stdout.flush()
                
                # Mark all structures in this batch as failed
                for item in batch:
                    failed += 1
                    results.append({
                        'structure_id': item['id'],
                        'pmg': None,
                        'composition': item['composition'],
                        'chemsys': item['chemsys'],
                        'mattersim_energy_per_atom': None,
                        'energy_above_hull': None,
                        'is_stable': None,
                        'passed_prescreening': False,
                        'error': f'Calculator creation failed: {str(e)}'
                    })
                
                # Save checkpoint and continue to next batch
                results_json = [{k: v for k, v in r.items() if k != 'pmg'} for r in results]
                checkpoint_data = {
                    'results': results_json,
                    'passed': passed,
                    'failed': failed,
                    'currently_processing': None,
                    'timestamp': datetime.now().isoformat()
                }
                with open(checkpoint_file, 'w') as f:
                    json.dump(checkpoint_data, f, indent=2)
                
                continue
            
            # Process structures in batch
            for item in tqdm(batch, desc=f"Batch {batch_idx + 1}", unit="struct"):
                struct_id = item['id']
                structure = item['structure']
                chemsys = item['chemsys']
                
                # Mark this structure as currently processing (for segfault detection)
                # Save immediately to detect crashes
                checkpoint_data_temp = {
                    'results': [{k: v for k, v in r.items() if k != 'pmg'} for r in results],
                    'passed': passed,
                    'failed': failed,
                    'currently_processing': struct_id,  # Mark what we're processing
                    'timestamp': datetime.now().isoformat()
                }
                with open(checkpoint_file, 'w') as f:
                    json.dump(checkpoint_data_temp, f, indent=2)
                
                # Pre-validate structure to catch malformed structures early
                is_valid, validation_error = validate_structure(structure)
                if not is_valid:
                    # Skip MatterSim relaxation for invalid structures
                    tqdm.write(f"  SKIPPING {struct_id}: {validation_error}")
                    failed += 1
                    results.append({
                        'structure_id': struct_id,
                        'pmg': None,
                        'composition': item['composition'],
                        'chemsys': chemsys,
                        'mattersim_energy_per_atom': None,
                        'energy_above_hull': None,
                        'is_stable': None,
                        'passed_prescreening': False,
                        'error': f"Invalid structure: {validation_error}"
                    })
                    continue
                
                # Relax using batch calculator (handles errors internally)
                relaxed_struct, energy_per_atom, relax_error = relax_structure_mattersim(
                    structure, calc, structure_id=struct_id
                )
                
                # Check if relaxation failed
                if relaxed_struct is None or energy_per_atom is None:
                    tqdm.write(f"  ERROR on {struct_id}: {relax_error}")
                    failed += 1
                    results.append({
                        'structure_id': struct_id,
                        'pmg': None,
                        'composition': item['composition'],
                        'chemsys': chemsys,
                        'mattersim_energy_per_atom': None,
                        'energy_above_hull': None,
                        'is_stable': None,
                        'passed_prescreening': False,
                        'error': relax_error
                    })
                    continue
                
                # Relaxation succeeded - compute hull distance
                mp_entries = mp_entries_cache.get(chemsys, [])
                
                if not mp_entries:
                    e_hull = None
                    passed_prescreening = True
                else:
                    e_hull = compute_energy_above_hull(relaxed_struct, energy_per_atom, mp_entries)
                    passed_prescreening = e_hull < args.hull_threshold
                    
                    if passed_prescreening:
                        passed += 1
                    else:
                        failed += 1
                
                results.append({
                    'structure_id': struct_id,
                    'pmg': relaxed_struct,
                    'composition': item['composition'],
                    'chemsys': chemsys,
                    'mattersim_energy_per_atom': float(energy_per_atom),
                    'energy_above_hull': float(e_hull) if e_hull is not None else None,
                    'is_stable': e_hull < 0.001 if e_hull is not None else None,
                    'passed_prescreening': passed_prescreening
                })
            
            # Clean up calculator after batch
            del calc
            if args.device == 'cuda':
                torch.cuda.empty_cache()
                gc.collect()
            
            # Save checkpoint after each batch (clear currently_processing marker)
            results_json = [{k: v for k, v in r.items() if k != 'pmg'} for r in results]
            checkpoint_data = {
                'results': results_json,
                'passed': passed,
                'failed': failed,
                'currently_processing': None,  # Clear marker after batch completes
                'timestamp': datetime.now().isoformat()
            }
            with open(checkpoint_file, 'w') as f:
                json.dump(checkpoint_data, f, indent=2)
            
            # Save structures to database incrementally after each batch
            batch_db_count = 0
            batch_failed_count = 0
            batch_start_idx = max(0, len(results) - len(batch))
            tolerances_db = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
            
            for i in range(batch_start_idx, len(results)):
                res = results[i]
                
                # Find original structure for this result
                orig_struct = None
                for s in batch:
                    if s['id'] == res['structure_id']:
                        orig_struct = s['structure']
                        break
                
                if res['pmg'] is not None:
                    # Valid structure with successful relaxation - try progressive tolerances
                    added = False
                    for tol in tolerances_db:
                        try:
                            xtal = pyxtal()
                            xtal.from_seed(res['pmg'], tol=tol)
                            if len(xtal.check_short_distance(r=0.5)) > 0:
                                continue
                            db.add_xtal(
                                xtal,
                                kvp={
                                    'structure_id': res['structure_id'],
                                    'e_above_hull': res['energy_above_hull'],
                                    'e_mattersim': res['mattersim_energy_per_atom'],
                                    'composition': res['composition'],
                                    'chemsys': res['chemsys'],
                                    'passed_prescreening': res['passed_prescreening'],
                                    'status': 'valid',
                                    'symmetrized': True
                                }
                            )
                            batch_db_count += 1
                            added = True
                            break
                        except Exception:
                            continue
                    
                    if not added:
                        # PyXtal failed on relaxed structure - try original structure with progressive tolerances
                        tqdm.write(f"  Warning: PyXtal failed on relaxed structure {res['structure_id']}, trying original structure")
                        if orig_struct is not None:
                            orig_added = False
                            for tol in tolerances_db:
                                try:
                                    xtal_orig = pyxtal()
                                    xtal_orig.from_seed(orig_struct, tol=tol)
                                    db.add_xtal(
                                        xtal_orig,
                                        kvp={
                                            'structure_id': res['structure_id'],
                                            'e_above_hull': res['energy_above_hull'],
                                            'e_mattersim': res['mattersim_energy_per_atom'],
                                            'composition': res['composition'],
                                            'chemsys': res['chemsys'],
                                            'passed_prescreening': False,  # Never proceed to DFT (use relaxed structure failed)
                                            'status': 'failed_symmetrization_after_relax',
                                            'symmetrized': True,
                                            'note': 'original_structure_saved_due_to_relaxed_pyxtal_failure'
                                        }
                                    )
                                    batch_db_count += 1
                                    orig_added = True
                                    break
                                except Exception:
                                    continue
                            
                            if not orig_added:
                                tqdm.write(f"  Error: Could not save even original structure for {res['structure_id']} (PyXtal failed at all tolerances)")
                                batch_failed_count += 1
                        else:
                            batch_failed_count += 1
                else:
                    # Invalid/failed structure - try to save original structure for book-keeping
                    if orig_struct is not None:
                        orig_added = False
                        for tol in tolerances_db:
                            try:
                                xtal_failed = pyxtal()
                                xtal_failed.from_seed(orig_struct, tol=tol)
                                db.add_xtal(
                                    xtal_failed,
                                    kvp={
                                        'structure_id': res['structure_id'],
                                        'e_above_hull': None,
                                        'e_mattersim': None,
                                        'composition': res['composition'],
                                        'chemsys': res['chemsys'],
                                        'passed_prescreening': False,  # Never proceed to DFT
                                        'status': res.get('error', 'failed_prescreening')[:200],  # Truncate long errors
                                        'symmetrized': True,
                                        'note': 'original_structure_saved_for_bookkeeping'
                                    }
                                )
                                batch_db_count += 1
                                orig_added = True
                                break
                            except Exception:
                                continue
                        
                        if not orig_added:
                            # Even original structure failed PyXtal at all tolerances - skip database entry
                            # These are still tracked in checkpoint JSON with full error details
                            batch_failed_count += 1
                    else:
                        batch_failed_count += 1
            
            # Explicit commit after each batch
            if hasattr(db, 'db') and hasattr(db.db, 'commit'):
                db.db.commit()
            
            print(f"\nBatch {batch_idx + 1} complete.")
            print(f"  Checkpoint saved: {checkpoint_file.name}")
            print(f"  Database updated: {batch_db_count} structures added")
            if batch_failed_count > 0:
                print(f"  Could not save {batch_failed_count} structures to database (tracked in checkpoint JSON)")
            print(f"  Progress: {len(results)}/{len(all_structures)} structures\n")
    
    # Save results to JSON (exclude pmg Structure objects)
    results_json = [{k: v for k, v in r.items() if k != 'pmg'} for r in results]
    output_data = {
        'summary': {
            'total_structures': len(all_structures),
            'passed_prescreening': passed,
            'failed_prescreening': failed,
            'hull_threshold': args.hull_threshold,
            'energy_reference': 'MatterSim-v1.0.0-5M'
        },
        'results': results_json
    }
    
    output_file = output_dir / f'prescreening_stability{batch_suffix}.json'
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nJSON results saved to: {output_file}")
    
    # Database was already saved incrementally during batch processing (see lines 770-808)
    # Get final database statistics
    print("\nDatabase Statistics:")
    if db_file.exists():
        try:
            # Count structures in database
            db_check = database_topology(str(db_file))
            if hasattr(db_check, 'db') and hasattr(db_check.db, 'execute'):
                cursor = db_check.db.execute("SELECT COUNT(*) FROM systems")
                db_count = cursor.fetchone()[0]
                print(f"  PyXtal database: {db_file}")
                print(f"  Structures in database: {db_count}")
            else:
                print(f"  PyXtal database: {db_file} (count unavailable)")
        except Exception as e:
            print(f"  PyXtal database: {db_file} (statistics unavailable: {e})")
    else:
        print(f"  No database file found (all structures may have failed)")
    
    # Remove checkpoint file (no longer needed)
    try:
        if checkpoint_file.exists():
            checkpoint_file.unlink()
            print(f"  Checkpoint file removed")
    except Exception as e:
        print(f"  Warning: Could not remove checkpoint file: {e}")
    
    print("\n" + "="*70)
    print("Pre-screening Complete")
    print("="*70)
    print(f"Total structures: {len(all_structures)}")
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")
    print(f"JSON: {output_file}")
    print(f"Database: {db_file}")
    print("="*70)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

