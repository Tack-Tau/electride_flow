#!/usr/bin/env python3
"""
Compute DFT energy_above_hull for VASP-relaxed structures.

Uses:
- VASP-PBE energies from completed relaxations (vasprun.xml)
- MP DFT-PBE UNCORRECTED energies for competing phases (queried via MP API)
- Cached MP entries from prescreening to minimize API calls

Important:
- Uses uncorrected_energy_per_atom from MP (raw DFT, no anion/composition corrections)
- This matches VASP energies which also have no MP-style corrections
- Filters out non-PBE functionals (R2SCAN, SCAN) when --pure-pbe is used

Output: dft_stability_results.json with DFT-level hull analysis
"""

import os
import sys
import json
import argparse
import time
import random
from pathlib import Path
from collections import defaultdict

from pymatgen.core import Composition
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from pymatgen.ext.matproj import MPRester

def load_workflow_database(db_path):
    """Load workflow.json database."""
    with open(db_path, 'r') as f:
        return json.load(f)


def get_vasp_energy_from_relax(relax_dir):
    """
    Extract final energy from VASP relaxation.
    
    Checks convergence before returning energy. If calculation is not converged,
    returns (None, None) which will be recorded as "VASP calculation not converged".
    
    Returns:
        tuple: (energy_per_atom, structure) or (None, None) if failed/not converged
    """
    vasprun_path = Path(relax_dir) / 'vasprun.xml'
    
    if not vasprun_path.exists():
        return None, None
    
    try:
        vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=False)
        
        # Check convergence first
        if not vr.converged:
            print(f"    Warning: VASP calculation not converged")
            return None, None
        
        final_energy = vr.final_energy  # Total energy in eV
        structure = vr.final_structure
        n_atoms = len(structure)
        energy_per_atom = final_energy / n_atoms
        
        return energy_per_atom, structure
        
    except Exception as e:
        print(f"    Error parsing vasprun.xml: {e}")
        return None, None


def save_dft_cache(new_entries_data, dft_cache_file):
    """
    Save DFT entries to single global cache file.
    
    Args:
        new_entries_data: List of dict entries to add to cache
        dft_cache_file: Path to global DFT cache file (mp_vaspdft.json)
    """
    dft_cache_file = Path(dft_cache_file)
    dft_cache_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Load existing cache
    existing_cache = {}
    if dft_cache_file.exists():
        with open(dft_cache_file, 'r') as f:
            cache_list = json.load(f)
            # Convert to dict for deduplication
            for item in cache_list:
                key = (item.get('chemsys', ''), item['entry_id'])
                existing_cache[key] = item
    
    # Add new entries (avoid duplicates)
    for new_item in new_entries_data:
        key = (new_item['chemsys'], new_item['entry_id'])
        existing_cache[key] = new_item
    
    # Save back to file
    all_cache_data = list(existing_cache.values())
    with open(dft_cache_file, 'w') as f:
        json.dump(all_cache_data, f, indent=2)
    
    print(f"    Updated DFT cache: {dft_cache_file} ({len(all_cache_data)} total entries)")


def load_dft_cache(chemsys, dft_cache_file):
    """
    Load DFT entries for a chemical system from single global cache file.
    
    Args:
        chemsys: Chemical system (e.g., 'Li-B-N')
        dft_cache_file: Path to global DFT cache file (mp_vaspdft.json)
    
    Returns:
        List of PDEntry objects for this chemsys, or None if cache doesn't exist
    """
    dft_cache_file = Path(dft_cache_file)
    
    if not dft_cache_file.exists():
        return None
    
    try:
        with open(dft_cache_file, 'r') as f:
            cache_list = json.load(f)
        
        # Extract entries for this chemsys and all subsystems
        elements = chemsys.split('-')
        entries = []
        
        for item in cache_list:
            item_chemsys = item.get('chemsys', '')
            item_elements = sorted(item_chemsys.split('-'))
            
            # Include if chemsys matches OR if it's a subset
            if set(item_elements).issubset(set(elements)):
                comp = Composition(item['composition'])
                entry = PDEntry(
                    composition=comp,
                    energy=item['energy'],
                    name=item['entry_id']
                )
                entries.append(entry)
        
        return entries if entries else None
    except Exception as e:
        print(f"    Warning: Could not load DFT cache: {e}")
        return None


def get_mp_stable_phases_dft(chemsys, mp_api_key, dft_cache_file, pure_pbe=False):
    """
    Get MP stable (on-hull) phases using UNCORRECTED DFT energies.
    
    Key optimizations:
    - Queries ONLY stable phases (is_stable=True) from MP
    - Uses single global cache file (mp_vaspdft.json) with chemsys field
    - Returns PDEntry objects for phase diagram analysis
    - Uses uncorrected_energy_per_atom (raw DFT, no MP corrections)
    
    Args:
        chemsys: Chemical system (e.g., 'Li-B-N')
        mp_api_key: Materials Project API key
        dft_cache_file: Path to global DFT cache file (mp_vaspdft.json)
        pure_pbe: If True, filter to GGA-PBE only (exclude PBE+U, R2SCAN, etc.)
                  If False (default), accept both PBE and PBE+U
    
    Returns:
        List of PDEntry objects for stable phases in this chemical system
    """
    # Check DFT cache first
    cached_entries = load_dft_cache(chemsys, dft_cache_file)
    if cached_entries is not None:
        print(f"    Loaded {len(cached_entries)} stable DFT entries from cache")
        return cached_entries
    
    # Need to query MP for stable phases and all subsystems
    print(f"    Querying MP for stable phases in {chemsys} and subsystems...")
    
    # Generate all required subsystems (elementals, binaries, ternaries, etc.)
    from itertools import combinations
    elements = chemsys.split('-')
    required_chemsys = set()
    for n in range(1, len(elements) + 1):
        for combo in combinations(elements, n):
            required_chemsys.add('-'.join(sorted(combo)))
    
    try:
        with MPRester(mp_api_key) as mpr:
            # Query ALL subsystems (including elementals and binaries)
            # This ensures we have terminal entries for phase diagram construction
            all_stable_docs = []
            
            for sys in sorted(required_chemsys):
                # Query only STABLE phases (is_stable=True)
                # This is MP's authoritative determination of thermodynamic stability
                stable_docs = mpr.materials.summary.search(
                    chemsys=sys,
                    is_stable=True,  # Only on-hull phases
                    fields=["material_id", "formula_pretty", "uncorrected_energy_per_atom", 
                            "energy_per_atom", "structure", "energy_above_hull"]
                )
                
                if stable_docs:
                    all_stable_docs.extend(stable_docs)
                
                # Random sleep to avoid API rate limiting (0-500 ms)
                time.sleep(random.uniform(0, 0.5))
            
            if not all_stable_docs:
                print(f"    Warning: No stable MP phases found for {chemsys} and subsystems")
                return []
            
            print(f"    Found {len(all_stable_docs)} stable phases across all subsystems")
            
            entries = []
            new_entries_data = []
            
            for doc in all_stable_docs:
                mp_id = doc.material_id
                
                if doc.structure is None:
                    continue
                
                # Determine chemsys for this phase
                doc_elements = sorted([str(el) for el in doc.structure.composition.elements])
                doc_chemsys = '-'.join(doc_elements)
                
                # Get uncorrected energy (raw DFT, no MP corrections)
                energy_per_atom = getattr(doc, 'uncorrected_energy_per_atom', None)
                if energy_per_atom is None:
                    # Fallback to energy_per_atom if uncorrected not available
                    energy_per_atom = getattr(doc, 'energy_per_atom', None)
                    if energy_per_atom is None:
                        print(f"      Warning: Skipping {mp_id} (missing energy)")
                        continue
                
                # Filter by functional if requested
                if pure_pbe:
                    try:
                        thermo_docs = mpr.materials.thermo.search(
                            material_ids=[mp_id],
                            fields=["entry_types"]
                        )
                        if thermo_docs and len(thermo_docs) > 0:
                            entry_types = getattr(thermo_docs[0], 'entry_types', [])
                            if entry_types:
                                et_str = ' '.join(str(et) for et in entry_types)
                                non_pbe_markers = ['+U', 'GGA_U', 'R2SCAN', 'SCAN', 'r2SCAN']
                                if any(marker in et_str for marker in non_pbe_markers):
                                    print(f"      Skipping {mp_id} (non-PBE: {et_str})")
                                    continue
                    except:
                        pass
                
                total_energy = energy_per_atom * doc.structure.composition.num_atoms
                
                # Create PDEntry
                entry = PDEntry(
                    composition=doc.structure.composition,
                    energy=total_energy,
                    name=mp_id
                )
                entries.append(entry)
                
                # Store for caching
                new_entries_data.append({
                    'chemsys': doc_chemsys,
                    'composition': {str(el): float(amt) for el, amt in doc.structure.composition.items()},
                    'energy': total_energy,
                    'entry_id': mp_id,
                    'mp_id': mp_id
                })
            
            print(f"    Retrieved {len(entries)} stable DFT phases")
            
            # Save to DFT cache
            if new_entries_data:
                save_dft_cache(new_entries_data, dft_cache_file)
            
            return entries
            
    except Exception as e:
        print(f"    Error querying MP for {chemsys}: {e}")
        import traceback
        traceback.print_exc()
        return []


def compute_dft_hull(target_entry, mp_entries):
    """
    Compute energy_above_hull using DFT energies with PDEntry.
    
    Args:
        target_entry: PDEntry for target structure (VASP total energy in eV)
        mp_entries: List of PDEntry from MP (DFT total energies in eV)
    
    Returns:
        tuple: (decomposition_products, energy_above_hull in eV/atom)
               decomposition_products is dict of {PDEntry: fraction}
    
    Note:
        PDEntry expects total energy (eV), not energy per atom.
        PhaseDiagram internally normalizes to per-atom basis for hull calculation.
        The returned e_hull is in eV/atom.
    """
    all_entries = [target_entry] + mp_entries
    
    try:
        pd = PhaseDiagram(all_entries)
        decomp, e_hull = pd.get_decomp_and_e_above_hull(target_entry, allow_negative=True)
        return decomp, e_hull
    except Exception as e:
        print(f"    Error computing phase diagram: {e}")
        return None, None


def main():
    parser = argparse.ArgumentParser(
        description="Compute DFT energy_above_hull for VASP-relaxed structures"
    )
    parser.add_argument(
        '--vasp-jobs',
        type=str,
        default='./VASP_JOBS',
        help="VASP jobs directory (default: ./VASP_JOBS)"
    )
    parser.add_argument(
        '--db',
        type=str,
        default='workflow.json',
        help="Workflow database (default: workflow.json in VASP_JOBS)"
    )
    parser.add_argument(
        '--mp-api-key',
        type=str,
        default=None,
        help="Materials Project API key (default: from MP_API_KEY env)"
    )
    parser.add_argument(
        '--output',
        type=str,
        default='dft_stability_results.json',
        help="Output JSON file (default: dft_stability_results.json)"
    )
    parser.add_argument(
        '--prescreen-results',
        type=str,
        default=None,
        help="Optional: prescreening_stability.json to filter structures"
    )
    parser.add_argument(
        '--pure-pbe',
        action='store_true',
        help="Filter MP entries to pure GGA-PBE only (exclude PBE+U). "
             "Default: accept both PBE and PBE+U for accurate phase diagrams. "
             "Use this flag only if your VASP calculations use pure PBE without +U corrections."
    )
    
    args = parser.parse_args()
    
    vasp_jobs = Path(args.vasp_jobs).expanduser()
    db_path = Path(args.db).expanduser()
    
    if not db_path.is_absolute():
        db_path = vasp_jobs / args.db
    
    output_path = vasp_jobs / args.output
    
    # Get MP API key
    mp_api_key = args.mp_api_key or os.environ.get('MP_API_KEY')
    if not mp_api_key:
        print("ERROR: MP_API_KEY not found in environment or arguments")
        print("Set it with: export MP_API_KEY=your_key")
        return 1
    
    if len(mp_api_key) != 32:
        print("ERROR: MP API key should be 32 characters (new API)")
        print("Get a new key from: https://next-gen.materialsproject.org/api")
        return 1
    
    # Set up cache paths (single global cache files)
    mp_mattersim_cache = vasp_jobs / "mp_mattersim.json"
    mp_dft_cache = vasp_jobs / "mp_vaspdft.json"
    
    print("="*70)
    print("DFT Energy Above Hull Calculation")
    print("="*70)
    print(f"VASP jobs directory: {vasp_jobs}")
    print(f"Workflow database: {db_path}")
    print(f"MP MatterSim cache: {mp_mattersim_cache}")
    print(f"MP DFT cache: {mp_dft_cache}")
    print(f"Output file: {output_path}")
    
    print(f"Energy source: MP stable phases ONLY (energy_above_hull = 0)")
    print(f"Energy corrections: NONE (using uncorrected DFT energies)")
    print("  VASP: raw DFT energies from vasprun.xml")
    print("  MP: uncorrected_energy_per_atom (no anion/composition corrections)")
    
    if args.pure_pbe:
        print(f"Functional filtering: Pure GGA-PBE only (PBE+U/R2SCAN/SCAN excluded)")
        print("  WARNING: This may reduce accuracy for transition metal systems")
    else:
        print(f"Functional filtering: Mixed PBE/PBE+U (MP recommended methodology)")
        print("  Note: Matches MP phase diagram methodology for best accuracy")
    
    print("="*70 + "\n")
    
    # Load workflow database
    print("Loading workflow database...")
    if not db_path.exists():
        print(f"ERROR: Workflow database not found: {db_path}")
        return 1
    
    db = load_workflow_database(db_path)
    print(f"Found {len(db['structures'])} structures in database\n")
    
    # Load pre-screening results if provided
    passed_structures = None
    if args.prescreen_results:
        prescreen_path = Path(args.prescreen_results).expanduser()
        
        # Handle relative paths: try as-is first, then try prepending vasp_jobs
        if not prescreen_path.is_absolute() and not prescreen_path.exists():
            # Try prepending vasp_jobs only if just a filename
            alt_path = vasp_jobs / args.prescreen_results
            if alt_path.exists():
                prescreen_path = alt_path
        
        if prescreen_path.exists():
            print(f"Loading pre-screening filter: {prescreen_path}")
            with open(prescreen_path, 'r') as f:
                prescreen_data = json.load(f)
            
            passed_structures = set()
            for result in prescreen_data.get('results', []):
                if result.get('passed_prescreening', False):
                    passed_structures.add(result['structure_id'])
            
            print(f"Will only process {len(passed_structures)} structures that passed pre-screening\n")
        else:
            print(f"Warning: Pre-screening file not found: {prescreen_path}")
            print("Processing all structures\n")
    
    # Scan for completed VASP relaxations
    print("Scanning for completed VASP relaxations...")
    completed_structures = []
    
    for struct_id, sdata in db['structures'].items():
        # Filter by pre-screening if provided
        if passed_structures is not None and struct_id not in passed_structures:
            continue
        
        state = sdata['state']
        
        # Consider RELAX_DONE and any downstream states as having completed relax
        completed_states = ['RELAX_DONE', 'SC_RUNNING', 'SC_DONE', 'SC_FAILED',
                          'PARCHG_RUNNING', 'PARCHG_DONE', 'PARCHG_FAILED', 'PARCHG_SKIPPED',
                          'ELF_RUNNING', 'ELF_DONE', 'ELF_FAILED']
        
        if state in completed_states:
            completed_structures.append(struct_id)
    
    print(f"Found {len(completed_structures)} completed structures\n")
    
    if not completed_structures:
        print("No completed structures found. Run VASP workflow first.")
        return 0
    
    # Group by chemical system for efficient MP queries
    print("Grouping structures by chemical system...")
    structures_by_chemsys = defaultdict(list)
    
    for struct_id in completed_structures:
        sdata = db['structures'][struct_id]
        chemsys = sdata.get('chemsys')
        
        if not chemsys:
            # Fallback: try to get from composition
            comp = Composition(sdata['composition'])
            elements = sorted([str(el) for el in comp.elements])
            chemsys = '-'.join(elements)
        
        structures_by_chemsys[chemsys].append(struct_id)
    
    print(f"Found {len(structures_by_chemsys)} unique chemical systems\n")
    
    # Process each chemical system
    results = []
    
    processed = 0
    failed = 0
    
    for chemsys, struct_ids in sorted(structures_by_chemsys.items()):
        print(f"\nProcessing {chemsys} ({len(struct_ids)} structures)...")
        
        # Get MP stable DFT entries once per chemical system
        print(f"  Querying MP for stable phases...")
        try:
            mp_entries = get_mp_stable_phases_dft(chemsys, mp_api_key, mp_dft_cache, pure_pbe=args.pure_pbe)
            
            if args.pure_pbe:
                print(f"  Retrieved {len(mp_entries)} stable phases from MP (pure GGA-PBE only)")
            else:
                print(f"  Retrieved {len(mp_entries)} stable phases from MP (mixed PBE/PBE+U)")
        except Exception as e:
            print(f"  ERROR: Could not get MP entries: {e}")
            failed += len(struct_ids)
            continue
        
        # Process each structure in this chemical system
        for struct_id in struct_ids:
            sdata = db['structures'][struct_id]
            relax_dir = Path(sdata['relax_dir'])
            
            print(f"\n  {struct_id}:")
            
            # Get VASP energy
            energy_per_atom, structure = get_vasp_energy_from_relax(relax_dir)
            
            if energy_per_atom is None:
                print(f"    FAILED: Could not extract VASP energy (not converged or parsing error)")
                failed += 1
                results.append({
                    'structure_id': struct_id,
                    'composition': sdata['composition'],
                    'chemsys': chemsys,
                    'vasp_energy_per_atom': None,
                    'dft_energy_above_hull': None,
                    'decomposition_products': None,
                    'error': 'VASP calculation not converged or vasprun.xml parsing failed'
                })
                continue
            
            print(f"    VASP energy: {energy_per_atom:.6f} eV/atom")
            
            # Create PDEntry for target structure
            composition = structure.composition
            total_energy = energy_per_atom * len(structure)
            target_entry = PDEntry(
                composition=composition,
                energy=total_energy,
                name=f"vasp_{struct_id}"
            )
            
            # Compute DFT hull
            decomp, e_hull = compute_dft_hull(target_entry, mp_entries)
            
            if e_hull is None:
                print(f"    FAILED: Could not compute phase diagram")
                failed += 1
                results.append({
                    'structure_id': struct_id,
                    'composition': sdata['composition'],
                    'chemsys': chemsys,
                    'vasp_energy_per_atom': energy_per_atom,
                    'dft_energy_above_hull': None,
                    'decomposition_products': None,
                    'error': 'Phase diagram calculation failed'
                })
                continue
            
            print(f"    DFT E_hull: {e_hull:.6f} eV/atom")
            
            # Format decomposition products
            decomp_str = None
            if decomp:
                decomp_parts = []
                for entry, frac in decomp.items():
                    formula = entry.composition.reduced_formula
                    decomp_parts.append(f"{formula} ({frac:.3f})")
                decomp_str = " + ".join(decomp_parts)
                print(f"    Decomposition: {decomp_str}")
            
            processed += 1
            
            results.append({
                'structure_id': struct_id,
                'composition': sdata['composition'],
                'chemsys': chemsys,
                'vasp_energy_per_atom': energy_per_atom,
                'dft_energy_above_hull': e_hull,
                'decomposition_products': decomp_str,
                'num_mp_stable_phases': len(mp_entries),
                'error': None
            })
    
    # Save results
    print("\n" + "="*70)
    print("Saving results...")
    
    output_data = {
        'summary': {
            'total_structures': len(completed_structures),
            'processed_successfully': processed,
            'failed': failed,
            'energy_reference': 'DFT uncorrected (VASP-PBE + MP stable phases only)' if not args.pure_pbe else 'DFT uncorrected (VASP-PBE + MP stable phases, pure PBE only)',
            'mp_phase_selection': 'stable_only (energy_above_hull = 0)',
            'mp_functional_filter': 'pure_pbe' if args.pure_pbe else 'mixed_pbe_pbeU',
            'mp_energy_corrections': 'none (using uncorrected_energy_per_atom)',
            'timestamp': str(Path(db_path).stat().st_mtime)
        },
        'results': results
    }
    
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"Results saved to: {output_path}")
    print("="*70)
    print(f"\nSummary:")
    print(f"  Total structures: {len(completed_structures)}")
    print(f"  Successfully processed: {processed}")
    print(f"  Failed: {failed}")
    print("="*70 + "\n")
    
    return 0


if __name__ == '__main__':
    exit(main())

