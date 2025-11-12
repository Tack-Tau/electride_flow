#!/usr/bin/env python3
"""
Compute DFT energy_above_hull for VASP-relaxed structures.

Uses:
- VASP-PBE energies from completed relaxations (vasprun.xml)
- MP DFT-PBE energies for competing phases (queried via MP API)
- Cached MP entries from prescreening to minimize API calls

Output: dft_stability_results.json with DFT-level hull analysis
"""

import os
import sys
import json
import argparse
from pathlib import Path
from collections import defaultdict

from pymatgen.core import Composition
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.ext.matproj import MPRester

def load_workflow_database(db_path):
    """Load workflow.json database."""
    with open(db_path, 'r') as f:
        return json.load(f)


def get_vasp_energy_from_relax(relax_dir):
    """
    Extract final energy from VASP relaxation.
    
    Returns:
        tuple: (energy_per_atom, structure) or (None, None) if failed
    """
    vasprun_path = Path(relax_dir) / 'vasprun.xml'
    
    if not vasprun_path.exists():
        return None, None
    
    try:
        vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=False)
        
        if not vr.converged:
            print(f"    Warning: Calculation not converged")
            return None, None
        
        final_energy = vr.final_energy  # Total energy in eV
        structure = vr.final_structure
        n_atoms = len(structure)
        energy_per_atom = final_energy / n_atoms
        
        return energy_per_atom, structure
        
    except Exception as e:
        print(f"    Error parsing vasprun.xml: {e}")
        return None, None


def extract_mp_ids_from_cache(cache_file):
    """
    Extract MP IDs from cached MatterSim file.
    
    Converts "mp_mattersim_mp-12345-GGA" to "mp-12345"
    """
    with open(cache_file, 'r') as f:
        cached_data = json.load(f)
    
    mp_ids = []
    for item in cached_data:
        entry_id = item['entry_id']
        # Extract mp-XXXXX from "mp_mattersim_mp-12345-GGA"
        if 'mp-' in entry_id:
            mp_id = entry_id.split('mp_mattersim_')[1].split('-GGA')[0]
            if mp_id not in mp_ids:
                mp_ids.append(mp_id)
    
    return mp_ids


def get_mp_dft_entries(chemsys, cache_dir, mp_api_key):
    """
    Get MP DFT entries for competing phases.
    
    Uses cached MP IDs from prescreening to minimize API calls.
    Returns original DFT energies from MP (not MatterSim).
    """
    cache_file = Path(cache_dir) / f"mp_{chemsys}_mattersim.json"
    
    if not cache_file.exists():
        print(f"    Warning: No cache file for {chemsys}, querying MP directly")
        # Fallback: query MP for all entries in chemical system
        with MPRester(mp_api_key) as mpr:
            docs = mpr.materials.summary.search(
                chemsys=chemsys,
                fields=["material_id", "energy_per_atom", "structure"]
            )
            entries = []
            for doc in docs:
                if doc.structure is None or doc.energy_per_atom is None:
                    continue
                
                total_energy = doc.energy_per_atom * doc.structure.composition.num_atoms
                entry = ComputedEntry(
                    composition=doc.structure.composition,
                    energy=total_energy,
                    entry_id=doc.material_id
                )
                entries.append(entry)
        return entries
    
    # Get MP IDs from cache
    mp_ids = extract_mp_ids_from_cache(cache_file)
    print(f"    Found {len(mp_ids)} MP entries in cache for {chemsys}")
    
    # Query MP for original DFT energies
    entries = []
    with MPRester(mp_api_key) as mpr:
        for mp_id in mp_ids:
            try:
                # Use new API to get energy and structure
                docs = mpr.materials.summary.search(
                    material_ids=[mp_id],
                    fields=["material_id", "energy_per_atom", "structure"]
                )
                
                if docs and len(docs) > 0:
                    doc = docs[0]
                    if doc.structure is None or doc.energy_per_atom is None:
                        print(f"    Warning: Skipping {mp_id} (missing structure or energy)")
                        continue
                    
                    total_energy = doc.energy_per_atom * doc.structure.composition.num_atoms
                    entry = ComputedEntry(
                        composition=doc.structure.composition,
                        energy=total_energy,
                        entry_id=doc.material_id
                    )
                    entries.append(entry)
            except Exception as e:
                print(f"    Warning: Could not get entry for {mp_id}: {e}")
    
    print(f"    Retrieved {len(entries)} DFT entries from MP")
    return entries


def compute_dft_hull(target_entry, mp_entries):
    """
    Compute energy_above_hull using DFT energies.
    
    Args:
        target_entry: ComputedEntry for target structure (VASP energy)
        mp_entries: List of ComputedEntry from MP (DFT energies)
    
    Returns:
        float: energy_above_hull in eV/atom
    """
    all_entries = [target_entry] + mp_entries
    
    try:
        pd = PhaseDiagram(all_entries)
        e_hull = pd.get_e_above_hull(target_entry)
        return e_hull
    except Exception as e:
        print(f"    Error computing phase diagram: {e}")
        return None


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
    
    print("="*70)
    print("DFT Energy Above Hull Calculation")
    print("="*70)
    print(f"VASP jobs directory: {vasp_jobs}")
    print(f"Workflow database: {db_path}")
    print(f"MP cache directory: {vasp_jobs / 'mp_mattersim_cache'}")
    print(f"Output file: {output_path}")
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
    cache_dir = vasp_jobs / 'mp_mattersim_cache'
    
    processed = 0
    failed = 0
    
    for chemsys, struct_ids in sorted(structures_by_chemsys.items()):
        print(f"\nProcessing {chemsys} ({len(struct_ids)} structures)...")
        
        # Get MP DFT entries once per chemical system
        print(f"  Querying MP for competing phases...")
        try:
            mp_entries = get_mp_dft_entries(chemsys, cache_dir, mp_api_key)
            print(f"  Retrieved {len(mp_entries)} competing phases from MP")
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
                print(f"    FAILED: Could not extract VASP energy")
                failed += 1
                results.append({
                    'structure_id': struct_id,
                    'composition': sdata['composition'],
                    'chemsys': chemsys,
                    'vasp_energy_per_atom': None,
                    'dft_energy_above_hull': None,
                    'error': 'Could not extract VASP energy from vasprun.xml'
                })
                continue
            
            print(f"    VASP energy: {energy_per_atom:.6f} eV/atom")
            
            # Create entry for target structure
            composition = structure.composition
            target_entry = ComputedEntry(
                composition=composition,
                energy=energy_per_atom * len(structure),
                entry_id=f"vasp_{struct_id}"
            )
            
            # Compute DFT hull
            e_hull = compute_dft_hull(target_entry, mp_entries)
            
            if e_hull is None:
                print(f"    FAILED: Could not compute phase diagram")
                failed += 1
                results.append({
                    'structure_id': struct_id,
                    'composition': sdata['composition'],
                    'chemsys': chemsys,
                    'vasp_energy_per_atom': energy_per_atom,
                    'dft_energy_above_hull': None,
                    'error': 'Phase diagram calculation failed'
                })
                continue
            
            print(f"    DFT E_hull: {e_hull:.6f} eV/atom")
            processed += 1
            
            results.append({
                'structure_id': struct_id,
                'composition': sdata['composition'],
                'chemsys': chemsys,
                'vasp_energy_per_atom': energy_per_atom,
                'dft_energy_above_hull': e_hull,
                'num_mp_competing_phases': len(mp_entries),
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
            'energy_reference': 'DFT (VASP-PBE + MP-PBE)',
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

