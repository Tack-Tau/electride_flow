#!/usr/bin/env python3
"""
MatterSim Pre-screening for VASPflow

Performs fast thermodynamic stability screening using MatterSim before expensive VASP calculations.
Outputs: prescreening_stability.json with structures passing energy_above_hull threshold.
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
from pathlib import Path
from io import StringIO
from datetime import datetime

from pymatgen.core import Composition
from pymatgen.io.cif import CifParser
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from pymatgen.ext.matproj import MPRester
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


def relax_structure_mattersim(pmg_struct, calculator, sym_tol=1e-1, fmax=0.01, max_steps=500, logfile=None):
    """
    Relax structure using MatterSim + FIRE optimizer with symmetrization and cell relaxation.
    
    Args:
        pmg_struct: Pymatgen Structure object
        calculator: MatterSimCalculator instance (reused across batch)
        sym_tol: Symmetry tolerance for PyXtal symmetrization (default 1e-1)
        fmax: Force convergence criterion (eV/Angstrom)
        max_steps: Maximum optimization steps
        logfile: None to suppress FIRE optimizer log, otherwise path to log file or set to '-' for stdout
    
    Returns:
        tuple: (relaxed_structure, energy_per_atom)
    
    Note:
        - Symmetrizes structure using PyXtal (important for mattergen structures)
        - Uses UnitCellFilter to relax both cell and atomic positions
        - Calculator is reused for efficiency (batch processing)
    """
    
    # Symmetrize structure using PyXtal
    xtal = pyxtal()
    xtal.from_seed(pmg_struct, tol=sym_tol)
    atoms = xtal.to_ase()
    
    # Attach calculator (reused from batch)
    atoms.calc = calculator
    
    # Use UnitCellFilter to relax both cell and atomic positions
    ecf = UnitCellFilter(atoms)
    dyn = FIRE(ecf, a=0.1, logfile=logfile)
    dyn.run(fmax=fmax, steps=max_steps)
    
    energy = atoms.get_potential_energy()
    energy_per_atom = energy / len(atoms)
    
    # Convert back to pymatgen structure
    adaptor = AseAtomsAdaptor()
    relaxed_structure = adaptor.get_structure(atoms)
    
    # Clean up atoms object (but keep calculator for reuse)
    del atoms
    del ecf
    del dyn
    
    return relaxed_structure, energy_per_atom


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
    Get MP stable (on-hull) phases and relax with MatterSim for consistent energy reference.
    
    Key optimizations:
    - Queries ONLY stable phases (is_stable=True) from MP
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
        List of PDEntry objects for stable phases in this chemical system
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
                print(f"    Using cached data: {cached_count} stable phases for {chemsys}")
                return entries
    
            # Need to query MP for missing data
            print(f"    Querying MP for stable phases in {chemsys} and subsystems...")
            if missing_chemsys:
                print(f"      Missing: {sorted(missing_chemsys)}")
    
            try:
                with MPRester(mp_api_key) as mpr:
                    # Query ALL missing subsystems (including elementals and binaries)
                    # This ensures we have terminal entries for phase diagram construction
                    all_stable_docs = []
                    systems_to_query = missing_chemsys if missing_chemsys else [chemsys]
                    
                    for sys in sorted(systems_to_query):
                        # Query only STABLE phases (is_stable=True)
                        # This is MP's authoritative determination of thermodynamic stability
                        stable_docs = mpr.materials.summary.search(
                            chemsys=sys,
                            is_stable=True,  # Only on-hull phases
                            fields=["material_id", "formula_pretty", "energy_per_atom", "structure", "energy_above_hull"]
                        )
                        
                        if stable_docs:
                            all_stable_docs.extend(stable_docs)
                        
                        # Random sleep to avoid API rate limiting (0-500 ms)
                        time.sleep(random.uniform(0, 0.5))
                    
                    if not all_stable_docs:
                        print(f"    Warning: No stable MP phases found for {chemsys} and subsystems")
                        return entries
                    
                    print(f"    Found {len(all_stable_docs)} stable phases across all subsystems, relaxing with MatterSim...")
                    
                    new_entries_data = []
                    success_count = 0
                    
                    for doc in all_stable_docs:
                        mp_id = doc.material_id
                        
                        # Determine chemsys for this phase
                        doc_elements = sorted([str(el) for el in doc.structure.composition.elements])
                        doc_chemsys = '-'.join(doc_elements)
                        
                        # Skip if already cached
                        cache_key = (doc_chemsys, mp_id)
                        if cache_key in cached_data:
                            continue
                        
                        try:
                            # Relax with MatterSim using reused calculator
                            relaxed_struct, energy_per_atom = relax_structure_mattersim(
                                doc.structure, calculator
                            )
                            
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
                            
                        except Exception as e:
                            print(f"      Warning: Failed to relax {mp_id}: {e}")
                    
                    print(f"    Successfully relaxed {success_count} new phases")
                    
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
    print("Pre-fetching MP Stable Phases (on-hull only)")
    print("="*70)
    print(f"Unique chemical systems: {len(unique_chemsys)}")
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
    if checkpoint_file.exists():
        print(f"Found checkpoint file, resuming from previous run...")
        with open(checkpoint_file, 'r') as f:
            checkpoint_data = json.load(f)
            results = checkpoint_data.get('results', [])
            passed = checkpoint_data.get('passed', 0)
            failed = checkpoint_data.get('failed', 0)
            processed_ids = {r['structure_id'] for r in results}
        print(f"Resuming: {len(results)}/{len(all_structures)} already processed\n")
        sys.stdout.flush()
    else:
        results = []
        passed = 0
        failed = 0
        processed_ids = set()
    
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
            
            # Create calculator for this batch
            calc = MatterSimCalculator(
                load_path="MatterSim-v1.0.0-5M.pth",
                device=args.device
            )
            
            # Process structures in batch
            for item in tqdm(batch, desc=f"Batch {batch_idx + 1}", unit="struct"):
                struct_id = item['id']
                structure = item['structure']
                chemsys = item['chemsys']
                
                try:
                    # Relax using batch calculator
                    relaxed_struct, energy_per_atom = relax_structure_mattersim(
                        structure, calc
                    )
                    
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
                        'energy_per_atom': energy_per_atom,
                        'composition': item['composition'],
                        'chemsys': chemsys,
                        'mattersim_energy_per_atom': float(energy_per_atom),
                        'energy_above_hull': float(e_hull) if e_hull is not None else None,
                        'is_stable': e_hull < 0.01 if e_hull is not None else None,
                        'passed_prescreening': passed_prescreening
                    })
                    
                except Exception as e:
                    tqdm.write(f"  ERROR on {struct_id}: {e}")
                    failed += 1
                    results.append({
                        'structure_id': struct_id,
                        'pmg': None,
                        'energy_per_atom': None,
                        'composition': item['composition'],
                        'chemsys': chemsys,
                        'mattersim_energy_per_atom': None,
                        'energy_above_hull': None,
                        'is_stable': None,
                        'passed_prescreening': False,
                        'error': str(e)
                    })
            
            # Clean up calculator after batch
            del calc
            if args.device == 'cuda':
                torch.cuda.empty_cache()
                gc.collect()
            
            # Save checkpoint after each batch
            results_json = [{k: v for k, v in r.items() if k != 'pmg'} for r in results]
            checkpoint_data = {
                'results': results_json,
                'passed': passed,
                'failed': failed,
                'timestamp': datetime.now().isoformat()
            }
            with open(checkpoint_file, 'w') as f:
                json.dump(checkpoint_data, f, indent=2)
            
            print(f"\nBatch {batch_idx + 1} complete. Checkpoint saved.")
            print(f"Progress: {len(results)}/{len(all_structures)} structures\n")
    
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
    
    # Save PyXtal database with relaxed structures
    print("\nSaving PyXtal database...")
    db_file = output_dir / f'prescreening_structures{batch_suffix}.db'
    db = database_topology(str(db_file))
    
    successful_count = 0
    fallback_count = 0
    for res in results:
        if res['pmg'] is not None:
            try:
                # Try with default tolerance first
                xtal = pyxtal()
                xtal.from_seed(res['pmg'], tol=1e-1)
                db.add_xtal(
                    xtal, 
                    kvp={
                        'structure_id': res['structure_id'],
                        'e_above_hull': res['energy_above_hull'],
                        'e_mattersim': res['energy_per_atom'],
                        'composition': res['composition'],
                        'chemsys': res['chemsys'],
                        'passed_prescreening': res['passed_prescreening']
                    }
                )
                successful_count += 1
            except Exception as e:
                # Try again with more permissive tolerance
                try:
                    xtal = pyxtal()
                    xtal.from_seed(res['pmg'], tol=0.5)
                    db.add_xtal(
                        xtal, 
                        kvp={
                            'structure_id': res['structure_id'],
                            'e_above_hull': res['energy_above_hull'],
                            'e_mattersim': res['energy_per_atom'],
                            'composition': res['composition'],
                            'chemsys': res['chemsys'],
                            'passed_prescreening': res['passed_prescreening']
                        }
                    )
                    successful_count += 1
                    fallback_count += 1
                    print(f"  Note: {res['structure_id']} added with relaxed tolerance (tol=0.5)")
                except Exception as e2:
                    print(f"  Warning: Could not add {res['structure_id']} to database: {e2}")
    
    print(f"PyXtal database saved to: {db_file}")
    print(f"  Structures in database: {successful_count}")
    if fallback_count > 0:
        print(f"  Structures using relaxed tolerance: {fallback_count}")
    
    # Remove checkpoint file (no longer needed)
    if checkpoint_file.exists():
        checkpoint_file.unlink()
    
    print("\n" + "="*70)
    print("Pre-screening Complete")
    print("="*70)
    print(f"Total structures: {len(all_structures)}")
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")
    print(f"JSON: {output_file}")
    print(f"Database: {db_file} ({successful_count} structures)")
    print("="*70)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

