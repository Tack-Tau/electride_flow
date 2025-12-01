#!/usr/bin/env python3
"""
Extract MP structures using modern MPRester (consistent with compute_dft_e_hull.py).

Reads chemical systems from mp_mattersim.json and queries MP directly
using modern MPRester with get_entries_in_chemsys() to get GGA/GGA+U structures.

Features:
- Uses modern mp_api.client.MPRester (same as compute_dft_e_hull.py)
- Queries by chemical system, not by MP ID
- Filters for GGA/GGA+U only (excludes r2SCAN/SCAN)
- Deduplicates by formula: prefers pure GGA over GGA+U (consistent with compute_dft_e_hull.py)
- Downloads structure files as CIF format
- Saves structures flat (no subdirectories per MP ID)

Output: Structures saved to mp_phase_dft_relax/mp_cache_structs/mp-XXXXX.cif
"""

import os
import sys
import json
import argparse
import time
import random
from pathlib import Path

from pymatgen.core import Structure

# Require modern mp-api client
try:
    from mp_api.client import MPRester
except ImportError:
    print("ERROR: mp-api package is required")
    print("Install with: pip install mp-api")
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


def query_and_download_structures(chemsys, output_dir, mp_api_key, skip_existing=True):
    """
    Query MP for GGA/GGA+U structures and save as CIF files.
    Uses modern MPRester with get_entries_in_chemsys() (same as compute_dft_e_hull.py).
    
    Args:
        chemsys: Chemical system (e.g., 'B-Li-N')
        output_dir: Output directory for CIF files
        mp_api_key: Materials Project API key
        skip_existing: If True, skip structures that already have CIF files
    
    Returns:
        dict: {mp_id: {'formula': str, 'num_sites': int, 'chemsys': str}} for downloaded structures
    
    Note:
        - Queries ALL GGA/GGA+U entries (no is_stable filter)
        - Filters out r2SCAN/SCAN by entry_id suffix
        - Deduplicates by formula: prefers pure GGA over GGA+U (consistent with compute_dft_e_hull.py)
        - Downloads one structure per formula (structure is same regardless of functional)
        - Saves flat: mp_cache_structs/mp-XXXXX.cif
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    elements = chemsys.split('-')
    downloaded = {}
    
    try:
        with MPRester(mp_api_key) as mpr:
            print(f"  Querying MP for GGA/GGA+U entries in {chemsys}...")
            
            # Get ALL entries (not filtered by stability, same as compute_dft_e_hull.py)
            computed_entries = mpr.get_entries_in_chemsys(
                elements,
                property_data=None,
                conventional_unit_cell=False,
                additional_criteria={}  # No is_stable filter
            )
            
            print(f"    Retrieved {len(computed_entries)} entries from MP")
            
            # Filter for GGA/GGA+U and extract structures
            # Use deduplication logic: prefer pure GGA over GGA+U for same formula
            gga_count = 0
            skipped_count = 0
            seen_formulas = {}  # Track duplicates by formula
            
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
                    is_gga = True
                
                # Skip r2SCAN/SCAN
                if has_SCAN:
                    continue
                
                # Skip if not GGA-based
                if not is_gga:
                    continue
                
                gga_count += 1
                
                # Extract base MP ID (without functional suffix)
                mp_id_parts = entry_id.split('-')
                if len(mp_id_parts) >= 2:
                    mp_id = mp_id_parts[0] + '-' + mp_id_parts[1]
                else:
                    mp_id = entry_id
                
                # Get structure from entry
                try:
                    structure = comp_entry.structure
                    formula = structure.composition.reduced_formula
                    
                    # Handle duplicates: prefer pure GGA over GGA+U for same formula
                    if formula in seen_formulas:
                        existing_has_U = seen_formulas[formula]['has_U']
                        current_has_U = has_U
                        
                        # Replace existing if current is pure GGA and existing is +U
                        if not current_has_U and existing_has_U:
                            # Remove old CIF file if it was downloaded
                            old_mp_id = seen_formulas[formula]['mp_id']
                            old_cif = output_dir / f"{old_mp_id}.cif"
                            if old_cif.exists():
                                old_cif.unlink()
                            # Remove from downloaded dict
                            if old_mp_id in downloaded:
                                del downloaded[old_mp_id]
                        else:
                            # Keep existing, skip current
                            continue
                    
                    # Skip if already downloaded and skip_existing is True
                    cif_file = output_dir / f"{mp_id}.cif"
                    if skip_existing and cif_file.exists() and formula not in seen_formulas:
                        skipped_count += 1
                        seen_formulas[formula] = {'mp_id': mp_id, 'has_U': has_U}
                        continue
                    
                    # Save structure as CIF
                    structure.to(filename=str(cif_file), fmt='cif')
                    
                    # Determine chemsys for this structure
                    struct_elements = sorted([str(el) for el in structure.composition.elements])
                    struct_chemsys = '-'.join(struct_elements)
                    
                    downloaded[mp_id] = {
                        'formula': formula,
                        'num_sites': len(structure),
                        'chemsys': struct_chemsys,
                        'entry_id': entry_id
                    }
                    
                    # Track this formula
                    seen_formulas[formula] = {'mp_id': mp_id, 'has_U': has_U}
                    
                except Exception as e:
                    print(f"      WARNING: Failed to extract structure for {mp_id}: {e}")
                    continue
            
            print(f"    Filtered to {gga_count} GGA/GGA+U entries")
            if skipped_count > 0:
                print(f"    Skipped {skipped_count} existing structures")
            print(f"    Downloaded {len(downloaded)} new structures")
            
    except Exception as e:
        print(f"    ERROR querying MP for {chemsys}: {e}")
        import traceback
        traceback.print_exc()
    
    return downloaded


def main():
    parser = argparse.ArgumentParser(
        description="Query and download MP structures using modern MPRester (consistent with compute_dft_e_hull.py)"
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
        help="Output directory for structures"
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
    print("MP Structure Download (Modern MPRester)")
    print("="*70)
    print(f"Cache file: {cache_file}")
    print(f"Output directory: {output_dir}")
    if args.chemsys:
        print(f"Chemical system filter: {args.chemsys}")
    else:
        print(f"Chemical system filter: All")
    print(f"Force re-download: {args.force}")
    print()
    print("Strategy:")
    print("  1. Extract chemical systems from mp_mattersim.json")
    print("  2. Query MP directly using get_entries_in_chemsys()")
    print("  3. Filter for GGA/GGA+U only (exclude r2SCAN/SCAN)")
    print("  4. Deduplicate: prefer pure GGA over GGA+U for same formula")
    print("  5. Download structures as CIF files")
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
                skip_existing=not args.force
            )
            
            all_downloaded.update(downloaded)
            print()
            
            # Random sleep to avoid API rate limiting
            time.sleep(random.uniform(0, 0.5))
            
    except Exception as e:
        print(f"ERROR: Failed to query/download structures: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    print("="*70)
    print("Summary")
    print("="*70)
    print(f"Chemical systems processed: {len(unique_chemsys)}")
    print(f"Total structures downloaded: {len(all_downloaded)}")
    print(f"Output directory: {output_dir}")
    print(f"Structure files: {output_dir}/mp-*.cif")
    print("="*70)
    print()
    print("Note:")
    print("  - Structures queried using modern MPRester (same as compute_dft_e_hull.py)")
    print("  - Deduplicated by formula: prefers pure GGA over GGA+U (same as compute_dft_e_hull.py)")
    print("  - MP energies are fetched separately by compare_mp_vasp_energies.py")
    print("  - This ensures GGA structure consistency across all scripts")
    print("="*70 + "\n")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

