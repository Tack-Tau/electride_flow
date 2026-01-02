#!/usr/bin/env python3
"""
Extract MP structures using legacy MPRester for complete GGA coverage.

Reads chemical systems from mp_mattersim.json and queries MP directly
using legacy MPRester with get_entries_in_chemsys() to get GGA/GGA+U structures.

Features:
- Uses legacy pymatgen.ext.matproj.MPRester for complete GGA entry retrieval
- Queries by chemical system, not by MP ID
- Strict filtering: only accepts entries with '-GGA' or '-GGA+U' suffix
- Optional --pure-pbe flag to exclude PBE+U (use pure GGA-PBE only)
- Saves MP GGA-PBE uncorrected energies to mp_vaspdft.json
- Downloads structure files as CIF format
- Saves structures flat (no subdirectories per MP ID)

Output: 
- Structures saved to mp_phase_dft_relax/mp_cache_structs/mp-XXXXX.cif
- Energies saved to mp_phase_dft_relax/mp_cache_structs/mp_vaspdft.json
"""

import os
import sys
import json
import argparse
import time
import random
from pathlib import Path

from pymatgen.core import Structure

# Use legacy pymatgen MPRester for complete GGA entry retrieval
try:
    from pymatgen.ext.matproj import MPRester
except ImportError:
    print("ERROR: pymatgen package with MPRester is required")
    print("Install with: pip install pymatgen")
    sys.exit(1)


def extract_chemsys_from_cache(cache_data, chemsys_filter=None):
    """
    Extract unique chemical systems from cached MatterSim JSON data.
    
    Args:
        cache_data: List of cached MP entries
        chemsys_filter: Optional chemical system filter (e.g., 'B-Li-N')
    
    Returns:
        set: Unique chemical systems
    """
    unique_chemsys = set()
    
    for item in cache_data:
        item_chemsys = item.get('chemsys')
        
        if not item_chemsys:
            continue
        
        # Filter by chemsys if requested
        if chemsys_filter and item_chemsys != chemsys_filter:
            continue
        
        unique_chemsys.add(item_chemsys)
    
    return unique_chemsys


def query_and_download_structures(chemsys, output_dir, mp_api_key, skip_existing=True, pure_pbe=False):
    """
    Query MP for GGA/GGA+U structures and save as CIF files with MP energies.
    Uses legacy MPRester with get_entries_in_chemsys() for complete GGA coverage.
    
    Args:
        chemsys: Chemical system (e.g., 'B-Li-N')
        output_dir: Output directory for CIF files
        mp_api_key: Materials Project API key
        skip_existing: If True, skip structures that already have CIF files
        pure_pbe: If True, filter to GGA-PBE only (exclude PBE+U)
    
    Returns:
        dict: {mp_id: {'formula': str, 'num_sites': int, 'chemsys': str, 'mp_energy_per_atom': float}} for downloaded structures
    
    Note:
        - Uses legacy pymatgen.ext.matproj.MPRester for complete GGA coverage
        - Strict filtering: only accepts entries with '-GGA' or '-GGA+U' suffix
        - Saves MP GGA-PBE uncorrected energies (raw DFT, no composition corrections)
        - Deduplicates by entry_id: keeps all polymorphs (no formula-based dedup)
        - Saves flat: mp_cache_structs/mp-XXXXX.cif
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    elements = chemsys.split('-')
    downloaded = {}
    
    try:
        mpr = MPRester(mp_api_key)
        print(f"  Querying MP for GGA/GGA+U entries in {chemsys}...")
        
        # Get ALL entries (not filtered by stability)
        computed_entries = mpr.get_entries_in_chemsys(elements)
        
        print(f"    Retrieved {len(computed_entries)} entries from MP (filtering for GGA...)")
        
        # Process entries: filter for GGA only (entry_id ending with '-GGA' or '-GGA+U')
        # Same strict filtering logic as prescreen.py, compute_dft_e_hull.py, etc.
        mp_phases = []
        seen_entries = {}  # Track by entry_id to avoid duplicates (keep all polymorphs)
        skipped_structure_retrieval = []
        
        for comp_entry in computed_entries:
            entry_id = str(comp_entry.entry_id)
            
            # Skip if already seen
            if entry_id in seen_entries:
                continue
            
            # Only accept entries ending with '-GGA' or '-GGA+U' (strict filtering)
            is_pure_gga = entry_id.endswith('-GGA')
            is_gga_u = entry_id.endswith('-GGA+U')
            
            # Skip non-GGA entries (r2SCAN, SCAN, or no suffix)
            if not is_pure_gga and not is_gga_u:
                continue
            
            # Skip +U if pure_pbe requested
            if pure_pbe and is_gga_u:
                continue
            
            has_U = is_gga_u
            
            # Extract base MP ID (e.g., 'mp-540703' from 'mp-540703-GGA')
            parts = entry_id.split('-')
            if len(parts) >= 2:
                mp_id = parts[0] + '-' + parts[1]
            else:
                mp_id = entry_id
            
            # Get structure
            structure = None
            try:
                structure = mpr.get_structure_by_material_id(mp_id)
            except Exception as e:
                skipped_structure_retrieval.append(f"{mp_id} ({comp_entry.composition.reduced_formula})")
                continue
            
            if structure is None:
                continue
            
            # Get MP GGA-PBE uncorrected energy (raw DFT, no composition corrections)
            # This is the energy we want to compare with VASP
            mp_energy_per_atom = comp_entry.energy_per_atom  # From ComputedEntry
            
            mp_phases.append((mp_id, structure, has_U, entry_id, float(mp_energy_per_atom)))
            seen_entries[entry_id] = True
        
        if skipped_structure_retrieval:
            print(f"    WARNING: Could not retrieve structures for {len(skipped_structure_retrieval)} phases")
        
        print(f"    Filtered to {len(mp_phases)} GGA phases (strict '-GGA'/'-GGA+U' suffix)")
        
        # Download structures and save info
        skipped_count = 0
        
        for mp_id, structure, has_U, entry_id, mp_energy_per_atom in mp_phases:
            cif_file = output_dir / f"{mp_id}.cif"
            
            # Skip if already downloaded and skip_existing is True
            if skip_existing and cif_file.exists():
                skipped_count += 1
                # Still add to downloaded dict for energy tracking
                struct_elements = sorted([str(el) for el in structure.composition.elements])
                struct_chemsys = '-'.join(struct_elements)
                downloaded[mp_id] = {
                    'formula': structure.composition.reduced_formula,
                    'num_sites': len(structure),
                    'chemsys': struct_chemsys,
                    'entry_id': entry_id,
                    'mp_energy_per_atom': mp_energy_per_atom
                }
                continue
            
            try:
                # Save structure as CIF
                structure.to(filename=str(cif_file), fmt='cif')
                
                # Determine chemsys for this structure
                struct_elements = sorted([str(el) for el in structure.composition.elements])
                struct_chemsys = '-'.join(struct_elements)
                
                downloaded[mp_id] = {
                    'formula': structure.composition.reduced_formula,
                    'num_sites': len(structure),
                    'chemsys': struct_chemsys,
                    'entry_id': entry_id,
                    'mp_energy_per_atom': mp_energy_per_atom
                }
                
            except Exception as e:
                print(f"      WARNING: Failed to save structure for {mp_id}: {e}")
                continue
        
        # Verify we have terminal (elemental) phases
        elements_found = set()
        for mp_id, structure, has_U, entry_id, mp_energy in mp_phases:
            if len(structure.composition.elements) == 1:
                elements_found.add(str(structure.composition.elements[0]))
        
        expected_elements = set(elements)
        if elements_found != expected_elements:
            missing = expected_elements - elements_found
            print(f"    WARNING: Missing terminal phases for elements: {sorted(missing)}")
        else:
            print(f"      All terminal phases present: {sorted(elements_found)}")
        
        if skipped_count > 0:
            print(f"    Skipped {skipped_count} existing structures")
        print(f"    Downloaded {len(downloaded) - skipped_count} new structures")
            
    except Exception as e:
        print(f"    ERROR querying MP for {chemsys}: {e}")
        import traceback
        traceback.print_exc()
    
    return downloaded


def main():
    parser = argparse.ArgumentParser(
        description="Query and download MP structures using legacy MPRester for complete GGA coverage"
    )
    parser.add_argument(
        '--cache-file',
        type=str,
        default='./VASP_JOBS/mp_mattersim.json',
        help="Path to mp_mattersim.json cache file (for extracting chemical systems)"
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='./mp_cache_structs',
        help="Output directory for structures and energies"
    )
    parser.add_argument(
        '--mp-api-key',
        type=str,
        default=None,
        help="Materials Project API key"
    )
    parser.add_argument(
        '--chemsys',
        type=str,
        default=None,
        help="Process specific chemical system (e.g., 'B-Li-N'). If not provided, process all."
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help="Force re-download even if CIF files already exist"
    )
    parser.add_argument(
        '--pure-pbe',
        action='store_true',
        help="Filter MP entries to pure GGA-PBE only (exclude PBE+U). "
             "Default: accept both PBE and PBE+U for accurate phase diagrams. "
             "Use this flag to match DFT calculations using pure PBE without +U corrections."
    )
    
    args = parser.parse_args()
    
    # Get MP API key
    mp_api_key = args.mp_api_key or os.environ.get('MP_API_KEY')
    if not mp_api_key:
        print("ERROR: MP_API_KEY not found in environment or arguments")
        print("Set it with: export MP_API_KEY=your_key")
        return 1
    
    cache_file = Path(args.cache_file).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*70)
    print("MP Structure Download (Legacy MPRester)")
    print("="*70)
    print(f"Cache file: {cache_file}")
    print(f"Output directory: {output_dir}")
    if args.chemsys:
        print(f"Chemical system filter: {args.chemsys}")
    else:
        print(f"Chemical system filter: All")
    print(f"Force re-download: {args.force}")
    if args.pure_pbe:
        print(f"Functional filtering: Pure GGA-PBE only (PBE+U excluded)")
    else:
        print(f"Functional filtering: Mixed PBE/PBE+U (recommended for accuracy)")
    print()
    print("Strategy:")
    print("  1. Extract chemical systems from mp_mattersim.json")
    print("  2. Query MP using legacy MPRester get_entries_in_chemsys()")
    print("  3. Strict filtering: only '-GGA' or '-GGA+U' suffix")
    print("  4. Save MP GGA-PBE uncorrected energies (raw DFT)")
    print("  5. Download structures as CIF files")
    print("  6. Save energies to mp_vaspdft.json for VASP comparison")
    print("="*70)
    
    # Load cache file
    if not cache_file.exists():
        print(f"ERROR: Cache file not found: {cache_file}")
        return 1
    
    print(f"\nLoading cache file...")
    with open(cache_file, 'r') as f:
        cache_data = json.load(f)
    print(f"Loaded {len(cache_data)} cached MP entries\n")
    
    # Extract unique chemical systems (not MP IDs!)
    unique_chemsys = extract_chemsys_from_cache(cache_data, chemsys_filter=args.chemsys)
    
    if not unique_chemsys:
        print("ERROR: No chemical systems found in cache!")
        if args.chemsys:
            print(f"  (chemsys filter: {args.chemsys})")
        return 1
    
    print(f"Found {len(unique_chemsys)} chemical system(s):")
    for cs in sorted(unique_chemsys):
        print(f"  {cs}")
    print()
    
    # Query and download structures for each chemical system
    print("="*70)
    print("Querying MP and Downloading Structures")
    print("="*70)
    print()
    
    all_downloaded = {}
    
    try:
        for chemsys in sorted(unique_chemsys):
            print(f"Processing {chemsys}...")
            
            downloaded = query_and_download_structures(
                chemsys,
                output_dir, 
                mp_api_key, 
                skip_existing=not args.force,
                pure_pbe=args.pure_pbe
            )
            
            all_downloaded.update(downloaded)
            print()
            
            # Random sleep to avoid API rate limiting
            time.sleep(random.uniform(0.5, 1.0))
            
    except Exception as e:
        print(f"ERROR: Failed to query/download structures: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    # Save MP energies to JSON file for later comparison
    energy_json_file = output_dir / "mp_vaspdft.json"
    energy_data = []
    
    for mp_id, info in sorted(all_downloaded.items()):
        energy_data.append({
            'chemsys': info['chemsys'],
            'composition': {info['formula']: info['num_sites']},  # Simple representation
            'energy': info['mp_energy_per_atom'] * info['num_sites'],  # Total energy
            'energy_per_atom': info['mp_energy_per_atom'],
            'entry_id': info['entry_id'],
            'mp_id': mp_id
        })
    
    with open(energy_json_file, 'w') as f:
        json.dump(energy_data, f, indent=2)
    
    print(f"\nSaved MP energies to: {energy_json_file}")
    print(f"  ({len(energy_data)} structures with GGA-PBE uncorrected energies)")
    
    print("="*70)
    print("Summary")
    print("="*70)
    print(f"Chemical systems processed: {len(unique_chemsys)}")
    print(f"Total structures with energies: {len(all_downloaded)}")
    print(f"Output directory: {output_dir}")
    print(f"Structure files: {output_dir}/mp-*.cif")
    print(f"Energy file: {energy_json_file}")
    print("="*70)
    print()
    print("Note:")
    print("  - Uses legacy pymatgen.ext.matproj.MPRester for complete GGA coverage")
    print("  - Strict filtering: only entries with '-GGA' or '-GGA+U' suffix")
    if args.pure_pbe:
        print("  - Filtered to pure GGA-PBE only (PBE+U excluded)")
    else:
        print("  - Includes both PBE and PBE+U entries (recommended)")
    print("  - Saved MP GGA-PBE uncorrected energies (raw DFT, no corrections)")
    print("  - compare_mp_vasp_energies.py can load energies from mp_vaspdft.json")
    print("  - This ensures complete GGA coverage and avoids re-querying MP")
    print("="*70 + "\n")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

