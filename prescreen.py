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
from pathlib import Path
from io import StringIO
from collections import defaultdict
from datetime import datetime

from pymatgen.core import Structure, Composition
from pymatgen.io.cif import CifParser
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram
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


def relax_structure_mattersim(pmg_struct, sym_tol=1e-1, device='cpu', fmax=0.01, max_steps=500, logfile=None):
    """
    Relax structure using MatterSim + FIRE optimizer with symmetrization and cell relaxation.
    
    Args:
        pmg_struct: Pymatgen Structure object
        sym_tol: Symmetry tolerance for PyXtal symmetrization (default 1e-1)
        device: 'cpu' or 'cuda'
        fmax: Force convergence criterion (eV/Angstrom)
        max_steps: Maximum optimization steps
        logfile: None to suppress FIRE optimizer log, otherwise path to log file or set to '-' for stdout
    
    Returns:
        tuple: (relaxed_structure, energy_per_atom)
    
    Note:
        - Symmetrizes structure using PyXtal (important for mattergen structures)
        - Uses UnitCellFilter to relax both cell and atomic positions
        - Creates a NEW calculator for each structure to avoid memory accumulation
    """
    
    # Symmetrize structure using PyXtal
    xtal = pyxtal()
    xtal.from_seed(pmg_struct, tol=sym_tol)
    atoms = xtal.to_ase()
    
    # Create fresh calculator for this structure (prevents memory issues)
    calc = MatterSimCalculator(
        load_path="MatterSim-v1.0.0-5M.pth",
        device=device
    )
    atoms.calc = calc
    
    # Use UnitCellFilter to relax both cell and atomic positions
    ecf = UnitCellFilter(atoms)
    dyn = FIRE(ecf, a=0.1, logfile=logfile)
    dyn.run(fmax=fmax, steps=max_steps)
    
    energy = atoms.get_potential_energy()
    energy_per_atom = energy / len(atoms)
    
    # Convert back to pymatgen structure
    adaptor = AseAtomsAdaptor()
    relaxed_structure = adaptor.get_structure(atoms)
    
    # Clean up
    del atoms
    del calc
    if device == 'cuda':
        import torch
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
    
    return relaxed_structure, energy_per_atom


def get_mp_competing_phases_mattersim(chemsys, mp_api_key, cache_dir, device='cpu'):
    """
    Get MP competing phases and re-relax with MatterSim for consistent energy reference.
    Results are cached per chemical system.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"mp_{chemsys}_mattersim.json"
    
    if cache_file.exists():
        print(f"    Using cached MP data for {chemsys}")
        with open(cache_file, 'r') as f:
            cached_data = json.load(f)
        
        entries = []
        for item in cached_data:
            entry = ComputedEntry(
                composition=item['composition'],
                energy=item['energy'],
                entry_id=item['entry_id']
            )
            entries.append(entry)
        
        # Check for missing elemental phases in cache
        elements = chemsys.split('-')
        present_elements = set()
        for entry in entries:
            comp = entry.composition
            if len(comp.elements) == 1:
                present_elements.add(str(comp.elements[0]))
        
        missing_elements = set(elements) - present_elements
        
        if missing_elements:
            print(f"    WARNING: Cache missing elemental phases: {sorted(missing_elements)}")
            print(f"    Cache will be regenerated...")
            # Force re-query by deleting cache and returning empty (will trigger full query)
            cache_file.unlink()
            print(f"    Deleted incomplete cache, will re-fetch from MP...")
            # Don't return here, let it fall through to the MP query below
        else:
            return entries
    
    print(f"    Querying MP for {chemsys} and all subsystems...")
    try:
        with MPRester(mp_api_key) as mpr:
            elements = chemsys.split('-')
            
            mp_entries = mpr.get_entries_in_chemsys(elements)
            
            if not mp_entries:
                print(f"    Warning: No MP entries found for {chemsys} and subsystems")
                return []
            
            phase_count = {}
            for entry in mp_entries:
                comp_els = tuple(sorted([str(el) for el in entry.composition.elements]))
                phase_count[comp_els] = phase_count.get(comp_els, 0) + 1
            
            for comp_els, count in sorted(phase_count.items()):
                sub_chemsys = '-'.join(comp_els)
                print(f"      {sub_chemsys}: {count} phases")
            
            print(f"    Total: {len(mp_entries)} MP phases, fetching structures and relaxing with MatterSim...")
            
            mattersim_entries = []
            fetch_errors = 0
            relax_errors = 0
            
            for mp_entry in mp_entries:
                try:
                    if hasattr(mp_entry, 'structure') and mp_entry.structure is not None:
                        mp_struct = mp_entry.structure
                    else:
                        mp_struct = mpr.get_structure_by_material_id(mp_entry.entry_id)
                        
                        if isinstance(mp_struct, list):
                            if len(mp_struct) == 0:
                                fetch_errors += 1
                                continue
                            mp_struct = mp_struct[0]
                    
                    if not hasattr(mp_struct, 'composition'):
                        fetch_errors += 1
                        continue
                    
                except Exception:
                    fetch_errors += 1
                    continue
                
                try:
                    relaxed_struct, energy_per_atom = relax_structure_mattersim(
                        mp_struct, device=device
                    )
                    
                    total_energy = energy_per_atom * relaxed_struct.composition.num_atoms
                    
                    entry = ComputedEntry(
                        composition=relaxed_struct.composition,
                        energy=total_energy,
                        entry_id=f"mp_mattersim_{mp_entry.entry_id}"
                    )
                    mattersim_entries.append(entry)
                except Exception:
                    relax_errors += 1
            
            print(f"    Successfully relaxed: {len(mattersim_entries)}/{len(mp_entries)} phases")
            if fetch_errors > 0:
                print(f"    Structure fetch errors: {fetch_errors}")
            if relax_errors > 0:
                print(f"    Relaxation errors: {relax_errors}")
            sys.stdout.flush()
            
            # Check for missing elemental phases and add fallback
            present_elements = set()
            for entry in mattersim_entries:
                comp = entry.composition
                if len(comp.elements) == 1:
                    present_elements.add(str(comp.elements[0]))
            
            missing_elements = set(elements) - present_elements
            
            if missing_elements:
                print(f"    WARNING: Missing elemental phases: {sorted(missing_elements)}")
                print(f"    Adding fallback elemental entries using MP reference energies...")
                
                for elem in sorted(missing_elements):
                    try:
                        # Query MP for elemental phases using summary.search
                        elem_docs = mpr.materials.summary.search(
                            elements=[elem],
                            num_elements=(1, 1),
                            fields=["material_id", "formula_pretty", "energy_per_atom", "structure"]
                        )
                        
                        if not elem_docs:
                            print(f"      ERROR: Could not find {elem} in MP!")
                            continue
                        
                        # Use the most stable elemental phase
                        elem_doc = sorted(elem_docs, key=lambda d: d.energy_per_atom)[0]
                        mp_id = elem_doc.material_id
                        mp_struct = elem_doc.structure
                        
                        print(f"      Found {elem} elemental phase: {mp_id} (MP E = {elem_doc.energy_per_atom:.4f} eV/atom)")
                        
                        # Relax with MatterSim
                        relaxed_struct, energy_per_atom = relax_structure_mattersim(
                            mp_struct, device=device
                        )
                        
                        total_energy = energy_per_atom * relaxed_struct.composition.num_atoms
                        
                        fallback_entry = ComputedEntry(
                            composition=relaxed_struct.composition,
                            energy=total_energy,
                            entry_id=f"mp_mattersim_{mp_id}_fallback"
                        )
                        mattersim_entries.append(fallback_entry)
                        print(f"      Added {elem} (MatterSim E = {energy_per_atom:.4f} eV/atom, fallback from {mp_id})")
                        
                    except Exception as e:
                        print(f"      ERROR: Failed to add {elem}: {e}")
                        import traceback
                        traceback.print_exc()
                
                sys.stdout.flush()
            
            if mattersim_entries:
                cached_data = []
                for entry in mattersim_entries:
                    cached_data.append({
                        'composition': entry.composition.as_dict(),
                        'energy': entry.energy,
                        'entry_id': entry.entry_id
                    })
                with open(cache_file, 'w') as f:
                    json.dump(cached_data, f, indent=2)
                print(f"    Cached {len(mattersim_entries)} MP phases for {chemsys}")
            
            return mattersim_entries
            
    except Exception as e:
        print(f"    Error querying MP for {chemsys}: {e}")
        return []


def compute_energy_above_hull(structure, energy_per_atom, mp_entries):
    """Compute energy above hull for a structure."""
    composition = structure.composition
    total_energy = energy_per_atom * composition.num_atoms
    
    entry = ComputedEntry(
        composition=composition,
        energy=total_energy,
        entry_id='generated'
    )
    
    all_entries = list(mp_entries) + [entry]
    pd = PhaseDiagram(all_entries)
    e_above_hull = pd.get_e_above_hull(entry)
    
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
        '--mp-cache-dir',
        type=str,
        default=None,
        help="Directory to cache MP data"
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
    
    mp_cache_dir = args.mp_cache_dir
    if mp_cache_dir is None:
        mp_cache_dir = output_dir / "mp_mattersim_cache"
    else:
        mp_cache_dir = Path(mp_cache_dir).expanduser()
    
    print("="*70)
    print("MatterSim Pre-screening for VASPflow")
    print("="*70)
    print(f"Results directory: {results_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Hull threshold: {args.hull_threshold} eV/atom")
    print(f"Device: {args.device}")
    print(f"Max compositions: {args.max_compositions or 'all'}")
    print(f"Max structures per composition: {args.max_structures}")
    print(f"MP cache directory: {mp_cache_dir}")
    print("="*70 + "\n")
    
    
    # Scan structures
    comp_dirs = sorted(results_dir.glob("*_structures"))
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
    
    # Pre-fetch MP data for all chemical systems
    unique_chemsys = sorted(set(s['chemsys'] for s in all_structures))
    
    print("="*70)
    print("Pre-fetching MP Competing Phases")
    print("="*70)
    print(f"Unique chemical systems: {len(unique_chemsys)}")
    print("="*70 + "\n")
    
    mp_entries_cache = {}
    for chemsys in unique_chemsys:
        print(f"Fetching MP data for {chemsys}...")
        sys.stdout.flush()  # Ensure output is visible immediately
        try:
            mp_entries = get_mp_competing_phases_mattersim(
                chemsys, mp_api_key, mp_cache_dir, device=args.device
            )
            mp_entries_cache[chemsys] = mp_entries
            print(f"  → Cached {len(mp_entries)} competing phases\n")
            sys.stdout.flush()
        except Exception as e:
            print(f"  → Error: {e}\n")
            sys.stdout.flush()
            mp_entries_cache[chemsys] = []
    
    print("="*70)
    print("MP data pre-fetch complete")
    print("="*70 + "\n")
    
    # Run pre-screening
    print("="*70)
    print("Running MatterSim Pre-screening")
    print("="*70)
    print(f"Processing {len(all_structures)} structures with progress tracking...\n")
    sys.stdout.flush()
    
    # Check for existing checkpoint
    checkpoint_file = output_dir / 'prescreening_checkpoint.json'
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
    
    # Process structures with progress bar
    from tqdm import tqdm
    for item in tqdm(all_structures, desc="Pre-screening", unit="struct"):
        struct_id = item['id']
        
        # Skip if already processed
        if struct_id in processed_ids:
            continue
        
        structure = item['structure']
        chemsys = item['chemsys']
        
        try:
            relaxed_struct, energy_per_atom = relax_structure_mattersim(
                structure, device=args.device
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
        
        # Save checkpoint every 50 structures
        if len(results) % 50 == 0:
            # Create JSON-safe version (exclude pmg Structure objects)
            results_json = [{k: v for k, v in r.items() if k != 'pmg'} for r in results]
            checkpoint_data = {
                'results': results_json,
                'passed': passed,
                'failed': failed,
                'timestamp': datetime.now().isoformat()
            }
            with open(checkpoint_file, 'w') as f:
                json.dump(checkpoint_data, f, indent=2)
    
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
    
    output_file = output_dir / 'prescreening_stability.json'
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nJSON results saved to: {output_file}")
    
    # Save PyXtal database with relaxed structures
    print("\nSaving PyXtal database...")
    db_file = output_dir / 'prescreening_structures.db'
    db = database_topology(str(db_file))
    
    successful_count = 0
    for res in results:
        if res['pmg'] is not None:
            try:
                xtal = pyxtal()
                xtal.from_seed(res['pmg'])
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
                print(f"  Warning: Could not add {res['structure_id']} to database: {e}")
    
    print(f"PyXtal database saved to: {db_file}")
    print(f"  Structures in database: {successful_count}")
    
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

