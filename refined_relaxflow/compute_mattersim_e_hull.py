#!/usr/bin/env python3
"""
Compute MatterSim energy_above_hull for refined VASP-relaxed structures.

This script:
1. Queries MP API for GGA/GGA+U reference phases in required chemical systems
2. Relaxes MP phases with MatterSim (fmax=0.001, max_steps=800 - matching refined VASP)
3. Loads high-precision VASP-relaxed structures from REFINE_VASP_JOBS
4. Relaxes refined structures with MatterSim (same tight convergence)
5. Computes energy_above_hull for refined structures using MatterSim-relaxed MP phases
6. Outputs mattersim_stability_results.json for direct DFT vs MatterSim comparison

Key differences from prescreen.py:
- Input: VASP-relaxed CONTCARs from refined 3-step relaxation (not CIFs from mattergen)
- Symmetrization: Structures are symmetrized via PyXtal
- Tighter convergence: fmax=0.001 vs 0.01, max_steps=800 vs 500
- MP phases: Query and relax with same tight convergence for consistency
- Output: Similar to compute_dft_e_hull.py for easy comparison

Usage:
    python3 compute_mattersim_e_hull.py \\
        --refine-jobs REFINE_VASP_JOBS/ \\
        --mp-api-key $MP_API_KEY \\
        --device cuda
"""

import os
import sys
import json
import argparse
import warnings
from pathlib import Path
from collections import defaultdict

from pymatgen.core import Structure, Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from pymatgen.io.ase import AseAtomsAdaptor

from ase.optimize import FIRE
from ase.filters import UnitCellFilter

try:
    from pyxtal import pyxtal
    PYXTAL_AVAILABLE = True
except ImportError:
    PYXTAL_AVAILABLE = False

try:
    from mattersim.forcefield import MatterSimCalculator
    MATTERSIM_AVAILABLE = True
except ImportError:
    MATTERSIM_AVAILABLE = False
    print("ERROR: MatterSim not available!")
    sys.exit(1)

try:
    from mp_api.client import MPRester
except ImportError:
    print("ERROR: mp-api package is required")
    print("Install with: pip install mp-api")
    sys.exit(1)

import numpy as np

warnings.filterwarnings('ignore', category=UserWarning, message='.*POTCAR data with symbol.*')
warnings.filterwarnings('ignore', message='Using UFloat objects with std_dev==0')
warnings.filterwarnings('ignore', category=DeprecationWarning, module='pkg_resources')
warnings.filterwarnings('ignore', category=UserWarning, message='.*pkg_resources is deprecated.*')


def relax_structure_mattersim(pmg_struct, calculator, structure_id=None, fmax=0.001, max_steps=800):
    """
    Relax structure using MatterSim + FIRE optimizer (tighter convergence than prescreen).
    
    Args:
        pmg_struct: Pymatgen Structure object
        calculator: MatterSimCalculator instance (reused)
        structure_id: Structure ID for error reporting
        fmax: Force convergence criterion (eV/Angstrom) - tighter than prescreen
        max_steps: Maximum optimization steps - more than prescreen
    
    Returns:
        tuple: (relaxed_structure, energy_per_atom, error_message)
    
    Note:
        - Symmetrizes structure using PyXtal with progressive tolerance (same as prescreen.py)
        - Falls back to direct conversion if all tolerances fail
        - Uses UnitCellFilter to relax both cell and atomic positions
    """
    sid_prefix = f"[{structure_id}] " if structure_id else ""
    
    # Symmetrize structure using PyXtal with progressive tolerance
    tolerances = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
    atoms = None
    
    for tol in tolerances:
        try:
            xtal = pyxtal()
            xtal.from_seed(pmg_struct, tol=tol)
            if not xtal.valid:
                continue
            if len(xtal.check_short_distances(r=0.5)) > 0:
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
            return None, None, f"{sid_prefix}Structure conversion failed: {e}"
    
    # Attach calculator
    atoms.calc = calculator
    
    # Check initial energy
    try:
        initial_energy = atoms.get_potential_energy()
        if not np.isfinite(initial_energy):
            return None, None, f"{sid_prefix}Initial energy not finite: {initial_energy}"
    except Exception as e:
        return None, None, f"{sid_prefix}Initial energy calculation failed: {e}"
    
    # Relax with error handling
    try:
        ecf = UnitCellFilter(atoms)
        dyn = FIRE(ecf, a=0.1, logfile=None)
        dyn.run(fmax=fmax, steps=max_steps)
        
        # Get final energy
        energy = atoms.get_potential_energy()
        
        # Check for NaN/inf
        if not np.isfinite(energy):
            try:
                del atoms, ecf, dyn
            except:
                pass
            return None, None, f"{sid_prefix}Final energy not finite: {energy}"
        
        energy_per_atom = energy / len(atoms)
        
        # Convert back to pymatgen
        adaptor = AseAtomsAdaptor()
        relaxed_structure = adaptor.get_structure(atoms)
        
    except Exception as e:
        try:
            del atoms
            if 'ecf' in locals():
                del ecf
            if 'dyn' in locals():
                del dyn
        except:
            pass
        return None, None, f"{sid_prefix}Relaxation failed: {e}"
    
    # Clean up
    del atoms, ecf, dyn
    
    return relaxed_structure, energy_per_atom, None


def load_workflow_database(db_path):
    """Load workflow.json database."""
    with open(db_path, 'r') as f:
        return json.load(f)


def load_mp_cache(cache_file):
    """
    Load MP cache from file.
    
    Returns:
        dict: Cache data keyed by (chemsys, entry_id)
    """
    cache_file = Path(cache_file)
    cached_data = {}
    
    if cache_file.exists():
        with open(cache_file, 'r') as f:
            cache_list = json.load(f)
            for item in cache_list:
                key = (item.get('chemsys', ''), item['entry_id'])
                cached_data[key] = item
    
    return cached_data


def save_mp_cache(cache_file, entries_data):
    """
    Save MP cache entries to file.
    
    Args:
        cache_file: Path to cache file
        entries_data: List of cache entries to save
    """
    cache_file = Path(cache_file)
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(cache_file, 'w') as f:
        json.dump(entries_data, f, indent=2)
    
    print(f"    Saved cache: {cache_file} ({len(entries_data)} total entries)")


def get_mp_stable_phases_mattersim(chemsys, mp_api_key, cached_data, calculator, fmax=0.001, max_steps=800):
    """
    Get MP GGA phases and relax with MatterSim using tight convergence.
    
    Uses same methodology as prescreen.py but with tighter convergence criteria
    to match refined VASP structures (fmax=0.001, max_steps=800).
    
    Args:
        chemsys: Chemical system string (e.g., 'B-Li-N')
        mp_api_key: Materials Project API key
        cached_data: Pre-loaded cache dict keyed by (chemsys, entry_id)
        calculator: MatterSimCalculator instance (reused)
        fmax: Force convergence criterion (default: 0.001, tighter than prescreen)
        max_steps: Maximum optimization steps (default: 800, more than prescreen)
    
    Returns:
        tuple: (List of PDEntry objects, List of new cache entries to add)
    """
            
    # Extract entries for this chemsys from cache
    elements = chemsys.split('-')
    entries = []
    cached_count = 0
    
    for key, item in cached_data.items():
        item_chemsys = key[0]
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
        return entries, []
    
    # Need to query MP for missing data
    print(f"    Querying MP for GGA phases in {chemsys}...")
    if missing_chemsys:
        print(f"      Missing subsystems: {sorted(missing_chemsys)}")
    
    try:
                with MPRester(mp_api_key) as mpr:
                    computed_entries = mpr.get_entries_in_chemsys(
                        elements,
                        property_data=None,
                        conventional_unit_cell=False,
                        additional_criteria={}
                    )
                    
                    print(f"    Retrieved {len(computed_entries)} entries from MP")
                    
                    # Filter for GGA/GGA+U, get structures
                    mp_phases = []
                    seen_formulas = {}
                    skipped_structure_retrieval = []
                    
                    for comp_entry in computed_entries:
                        entry_id = comp_entry.entry_id
                        
                        # Filter by functional
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
                            is_gga = True
                        
                        if has_SCAN or not is_gga:
                            continue
                        
                        # Get MP ID
                        mp_id_parts = entry_id.split('-')
                        if len(mp_id_parts) >= 2:
                            mp_id = mp_id_parts[0] + '-' + mp_id_parts[1]
                        else:
                            mp_id = entry_id
                        
                        # Extract structure
                        structure = None
                        try:
                            structure = mpr.get_structure_by_material_id(mp_id)
                        except Exception:
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
                                mp_phases = [(mid, s, has_u) for mid, s, has_u in mp_phases 
                                            if s.composition.reduced_formula != formula]
                            else:
                                continue
                        
                        mp_phases.append((mp_id, structure, has_U))
                        seen_formulas[formula] = {'mp_id': mp_id, 'has_U': has_U}
                    
                    if skipped_structure_retrieval:
                        print(f"    WARNING: Could not retrieve structures for {len(skipped_structure_retrieval)} phases")
                    
                    print(f"    Filtered to {len(mp_phases)} GGA/GGA+U phases with structures")
                    
                    # Verify terminal phases
                    elements_found = set()
                    for mp_id, structure, has_U in mp_phases:
                        if len(structure.composition.elements) == 1:
                            elements_found.add(str(structure.composition.elements[0]))
                    
                    expected_elements = set(elements)
                    if elements_found != expected_elements:
                        missing = expected_elements - elements_found
                        print(f"    WARNING: Missing terminal phases for elements: {sorted(missing)}")
                        print(f"    Querying elemental phases directly...")
                        
                        for missing_elem in sorted(missing):
                            try:
                                elem_docs = mpr.materials.summary.search(
                                    elements=[missing_elem],
                                    num_elements=(1, 1),
                                    fields=['material_id', 'structure', 'energy_per_atom']
                                )
                                
                                if not elem_docs:
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
                                    print(f"        Found {elem_mp_id} for {missing_elem}")
                                    mp_phases.append((elem_mp_id, elem_structure, False))
                                    elements_found.add(missing_elem)
                            except Exception as e:
                                print(f"        ERROR querying {missing_elem}: {e}")
                    else:
                        print(f"      All terminal phases present: {sorted(elements_found)}")
                    
                    print(f"    Relaxing {len(mp_phases)} phases with MatterSim (fmax={fmax}, max_steps={max_steps})...")
                    
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
                        
                        # Relax with MatterSim using tight convergence
                        relaxed_struct, energy_per_atom, error_msg = relax_structure_mattersim(
                            structure, calculator, structure_id=mp_id, fmax=fmax, max_steps=max_steps
                        )
                        
                        if relaxed_struct is None or energy_per_atom is None:
                            print(f"      Warning: Failed to relax {mp_id}: {error_msg}")
                            failed_count += 1
                            continue
                        
                        total_energy = energy_per_atom * relaxed_struct.composition.num_atoms
                        
                        entry = PDEntry(
                            composition=relaxed_struct.composition,
                            energy=total_energy,
                            name=f"mp_mattersim_{mp_id}"
                        )
                        entries.append(entry)
                        
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
                
                return entries, new_entries_data
                
    except Exception as e:
        print(f"    Error querying MP for {chemsys}: {e}")
        import traceback
        traceback.print_exc()
        return entries, []


def compute_energy_above_hull(structure, energy_per_atom, mp_entries):
    """
    Compute energy above hull using PDEntry.
    
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


def main():
    parser = argparse.ArgumentParser(
        description="Compute MatterSim energy_above_hull for refined VASP structures"
    )
    parser.add_argument(
        '--refine-jobs',
        type=str,
        default='./REFINE_VASP_JOBS',
        help="Refined VASP jobs directory (default: ./REFINE_VASP_JOBS)"
    )
    parser.add_argument(
        '--db',
        type=str,
        default='workflow.json',
        help="Workflow database (default: workflow.json in refine-jobs)"
    )
    parser.add_argument(
        '--mp-api-key',
        type=str,
        default=None,
        help="Materials Project API key (default: from MP_API_KEY env)"
    )
    parser.add_argument(
        '--device',
        type=str,
        default='cpu',
        choices=['cpu', 'cuda'],
        help="Device for MatterSim (default: cpu)"
    )
    parser.add_argument(
        '--output',
        type=str,
        default='mattersim_stability_results.json',
        help="Output JSON file (default: mattersim_stability_results.json)"
    )
    
    args = parser.parse_args()
    
    refine_jobs = Path(args.refine_jobs).expanduser()
    db_path = Path(args.db).expanduser()
    
    if not db_path.is_absolute():
        db_path = refine_jobs / args.db
    
    output_path = refine_jobs / args.output
    mp_cache_local = refine_jobs / "mp_mattersim.json"
    
    # Get MP API key
    mp_api_key = args.mp_api_key or os.environ.get('MP_API_KEY')
    if not mp_api_key:
        print("ERROR: Materials Project API key required!")
        print("  Use: --mp-api-key YOUR_KEY or export MP_API_KEY=YOUR_KEY")
        return 1
    
    print("="*70)
    print("MatterSim Energy Above Hull Calculation (Refined Structures)")
    print("="*70)
    print(f"Refined VASP jobs: {refine_jobs}")
    print(f"Workflow database: {db_path}")
    print(f"MP MatterSim cache: {mp_cache_local}")
    print(f"Output file: {output_path}")
    print(f"Device: {args.device}")
    print(f"Convergence: fmax=0.001 eV/Å, max_steps=800 (tighter than prescreen)")
    print("="*70 + "\n")
    
    # Load workflow database
    print("Loading workflow database...")
    if not db_path.exists():
        print(f"ERROR: Workflow database not found: {db_path}")
        return 1
    
    db = load_workflow_database(db_path)
    print(f"Found {len(db['structures'])} structures in database\n")
    
    # Scan for completed relaxations
    print("Scanning for RELAX_DONE structures...")
    completed_structures = []
    
    for struct_id, sdata in db['structures'].items():
        if sdata['state'] == 'RELAX_DONE':
            completed_structures.append(struct_id)
    
    print(f"Found {len(completed_structures)} completed refined relaxations\n")
    
    if not completed_structures:
        print("No completed structures found. Run refined workflow first.")
        return 0
    
    # Group by chemical system
    print("Grouping structures by chemical system...")
    structures_by_chemsys = defaultdict(list)
    
    for struct_id in completed_structures:
        sdata = db['structures'][struct_id]
        chemsys = sdata.get('chemsys')
        
        if not chemsys:
            comp = Composition(sdata['composition'])
            elements = sorted([str(el) for el in comp.elements])
            chemsys = '-'.join(elements)
        
        structures_by_chemsys[chemsys].append(struct_id)
    
    required_chemsys = set(structures_by_chemsys.keys())
    unique_chemsys = sorted(required_chemsys)
    print(f"Found {len(unique_chemsys)} unique chemical systems")
    print(f"  {unique_chemsys}\n")
    
    # Create MatterSim calculator (for both MP phases and structures)
    print("Creating MatterSimCalculator...")
    try:
        import torch
        import gc
        calc = MatterSimCalculator(
            load_path="MatterSim-v1.0.0-5M.pth",
            device=args.device
        )
        print("  Calculator created successfully\n")
    except Exception as e:
        print(f"ERROR: Failed to create MatterSimCalculator: {e}")
        return 1
    
    # Load existing MP cache
    print("Loading existing MP cache...")
    cached_data = load_mp_cache(mp_cache_local)
    print(f"  Loaded {len(cached_data)} cached entries\n")
    
    # Pre-fetch and relax MP stable phases with tight convergence
    print("="*70)
    print("Querying and Relaxing MP GGA Phases (Tight Convergence)")
    print("="*70)
    print(f"Chemical systems: {len(unique_chemsys)}")
    print(f"Convergence: fmax=0.001 eV/Å, max_steps=800 (matching refined structures)")
    print("Strategy: Query ALL GGA/GGA+U entries, relax with MatterSim")
    print("="*70 + "\n")
    
    mp_entries_cache = {}
    all_new_cache_entries = []
    
    for chemsys in unique_chemsys:
        print(f"Fetching MP stable phases for {chemsys}...")
        sys.stdout.flush()
        try:
            mp_entries, new_cache_entries = get_mp_stable_phases_mattersim(
                chemsys, mp_api_key, cached_data, calc, fmax=0.001, max_steps=800
            )
            mp_entries_cache[chemsys] = mp_entries
            
            # Collect new cache entries
            for new_entry in new_cache_entries:
                key = (new_entry['chemsys'], new_entry['entry_id'])
                if key not in cached_data:
                    all_new_cache_entries.append(new_entry)
                    cached_data[key] = new_entry
            
            print(f"  → {len(mp_entries)} stable phases ready\n")
            sys.stdout.flush()
        except Exception as e:
            print(f"  → Error: {e}\n")
            sys.stdout.flush()
            mp_entries_cache[chemsys] = []
    
    # Save updated cache
    if all_new_cache_entries:
        print(f"Saving updated MP cache ({len(all_new_cache_entries)} new entries)...")
        all_cache_list = list(cached_data.values())
        save_mp_cache(mp_cache_local, all_cache_list)
    
    print("="*70)
    print("MP stable phases query/relaxation complete")
    print("="*70 + "\n")
    
    # Process structures
    print("="*70)
    print("Relaxing Refined Structures with MatterSim")
    print("="*70)
    print(f"Processing {len(completed_structures)} structures...\n")
    
    results = []
    processed = 0
    failed = 0
    
    for chemsys, struct_ids in sorted(structures_by_chemsys.items()):
        print(f"\nProcessing {chemsys} ({len(struct_ids)} structures)...")
        
        mp_entries = mp_entries_cache.get(chemsys, [])
        print(f"  Using {len(mp_entries)} MP reference phases\n")
        
        for struct_id in struct_ids:
            sdata = db['structures'][struct_id]
            relax_dir = Path(sdata['relax_dir'])
            contcar_path = relax_dir / 'CONTCAR'
            
            print(f"  {struct_id}:")
            
            # Load VASP-relaxed structure
            if not contcar_path.exists():
                print(f"    FAILED: CONTCAR not found")
                failed += 1
                results.append({
                    'structure_id': struct_id,
                    'composition': sdata['composition'],
                    'chemsys': chemsys,
                    'mattersim_energy_per_atom': None,
                    'energy_above_hull': None,
                    'is_stable': None,
                    'passed_prescreening': False,
                    'error': 'CONTCAR not found'
                })
                continue
            
            try:
                structure = Structure.from_file(str(contcar_path))
            except Exception as e:
                print(f"    FAILED: Could not load CONTCAR: {e}")
                failed += 1
                results.append({
                    'structure_id': struct_id,
                    'composition': sdata['composition'],
                    'chemsys': chemsys,
                    'mattersim_energy_per_atom': None,
                    'energy_above_hull': None,
                    'is_stable': None,
                    'passed_prescreening': False,
                    'error': f'CONTCAR parsing failed: {e}'
                })
                continue
            
            # Relax with MatterSim (tighter convergence)
            relaxed_struct, energy_per_atom, error_msg = relax_structure_mattersim(
                structure, calc, structure_id=struct_id, fmax=0.001, max_steps=800
            )
            
            if relaxed_struct is None or energy_per_atom is None:
                print(f"    FAILED: {error_msg}")
                failed += 1
                results.append({
                    'structure_id': struct_id,
                    'composition': sdata['composition'],
                    'chemsys': chemsys,
                    'mattersim_energy_per_atom': None,
                    'energy_above_hull': None,
                    'is_stable': None,
                    'passed_prescreening': False,
                    'error': error_msg
                })
                continue
            
            print(f"    MatterSim energy: {energy_per_atom:.6f} eV/atom")
            
            # Compute hull
            if not mp_entries:
                e_hull = None
                print(f"    WARNING: No MP reference phases, skipping hull calculation")
            else:
                try:
                    e_hull = compute_energy_above_hull(relaxed_struct, energy_per_atom, mp_entries)
                    print(f"    E_hull: {e_hull:.6f} eV/atom")
                except Exception as e:
                    print(f"    FAILED: Hull calculation error: {e}")
                    e_hull = None
            
            is_stable = bool(e_hull < 0.001) if e_hull is not None else None
            processed += 1
            
            results.append({
                'structure_id': struct_id,
                'composition': sdata['composition'],
                'chemsys': chemsys,
                'mattersim_energy_per_atom': float(energy_per_atom),
                'energy_above_hull': float(e_hull) if e_hull is not None else None,
                'is_stable': is_stable,
                'passed_prescreening': True,
                'error': None
            })
    
    # Clean up calculator
    del calc
    if args.device == 'cuda':
        import torch
        torch.cuda.empty_cache()
        import gc
        gc.collect()
    
    # Save results
    print("\n" + "="*70)
    print("Saving MatterSim hull results...")
    
    output_data = {
        'summary': {
            'total_structures': len(completed_structures),
            'successfully_processed': processed,
            'failed': failed,
            'energy_reference': 'MatterSim-v1.0.0-5M',
            'convergence': 'fmax=0.001, max_steps=800',
            'mp_reference_source': 'MP API GGA/GGA+U phases relaxed with MatterSim (same convergence)',
            'mp_cache': str(mp_cache_local)
        },
        'results': results
    }
    
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"MatterSim results saved to: {output_path}")
    print("="*70)
    print(f"\nMatterSim Hull Summary:")
    print(f"  Total structures: {len(completed_structures)}")
    print(f"  Successfully processed: {processed}")
    print(f"  Failed: {failed}")
    print("="*70 + "\n")
    
    return 0


if __name__ == '__main__':
    exit(main())

