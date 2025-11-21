#!/usr/bin/env python3
"""
Extract MP structures from cached JSON file and download them via MP API.

Reads mp_mattersim.json from VASP_JOBS/ directory
and downloads crystal structures from Materials Project.

Features:
- Downloads structure files only (CIF format)
- Saves structures flat (no subdirectories per MP ID)
- Lightweight and fast

Output: Structures saved to mp_phase_dft_relax/mp_cache_structs/mp-XXXXX.cif
"""

import os
import sys
import json
import argparse
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.ext.matproj import MPRester


def extract_mp_ids_from_cache(cache_data, chemsys=None):
    """
    Extract MP IDs from cached MatterSim JSON data.
    
    Args:
        cache_data: List of cached MP entries
        chemsys: Optional chemical system filter (e.g., 'B-Li-N')
    
    Returns:
        dict: {chemsys: [mp_id1, mp_id2, ...]}
    
    Converts "mp_mattersim_mp-12345" or "mp_mattersim_mp-12345_fallback" to "mp-12345"
    """
    by_chemsys = {}
    
    for item in cache_data:
        entry_id = item['entry_id']
        item_chemsys = item.get('chemsys')
        
        # Filter by chemsys if requested
        if chemsys and item_chemsys != chemsys:
            continue
        
        # Extract mp-XXXXX from various formats
        if 'mp-' in entry_id:
            # Split by 'mp_mattersim_' and take the part after it
            mp_id = entry_id.split('mp_mattersim_')[1]
            # Remove '_fallback' suffix if present
            mp_id = mp_id.replace('_fallback', '')
            # Remove any GGA suffix
            mp_id = mp_id.split('-GGA')[0]
            
            if item_chemsys not in by_chemsys:
                by_chemsys[item_chemsys] = []
            
            if mp_id not in by_chemsys[item_chemsys]:
                by_chemsys[item_chemsys].append(mp_id)
    
    return by_chemsys


def download_mp_structures(mp_ids, output_dir, mp_api_key, skip_existing=True):
    """
    Download structures from MP API and save as CIF files.
    
    Args:
        mp_ids: List of MP IDs to download
        output_dir: Output directory for CIF files
        mp_api_key: Materials Project API key
        skip_existing: If True, skip MP IDs that already have CIF files
    
    Returns:
        dict: {mp_id: {'formula': str, 'num_sites': int}} for successfully downloaded structures
    
    Note:
        - Downloads structure files only (no energy metadata)
        - Saves flat: mp_cache_structs/mp-XXXXX.cif
        - Energies are fetched separately by compare_mp_vasp_energies.py
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Separate new vs already-downloaded
    if skip_existing:
        new_mp_ids = [mp_id for mp_id in mp_ids if not (output_dir / f"{mp_id}.cif").exists()]
        existing_mp_ids = [mp_id for mp_id in mp_ids if (output_dir / f"{mp_id}.cif").exists()]
        print(f"  Total MP IDs: {len(mp_ids)}")
        print(f"  New (to download): {len(new_mp_ids)}")
        print(f"  Existing (skipping): {len(existing_mp_ids)}")
    else:
        new_mp_ids = mp_ids
        print(f"  Total MP IDs: {len(mp_ids)}")
    
    downloaded = {}
    
    with MPRester(mp_api_key) as mpr:
        for mp_id in new_mp_ids:
            try:
                # Query structure from MP (structure only, no energies)
                docs = mpr.materials.summary.search(
                    material_ids=[mp_id],
                    fields=["material_id", "structure", "formula_pretty"]
                )
                
                if not docs or len(docs) == 0:
                    print(f"  WARNING: {mp_id} not found in MP")
                    continue
                
                doc = docs[0]
                if doc.structure is None:
                    print(f"  WARNING: {mp_id} has no structure")
                    continue
                
                structure = doc.structure
                formula = doc.formula_pretty
                
                # Save structure as CIF (flat, no subdirectories)
                cif_file = output_dir / f"{mp_id}.cif"
                structure.to(filename=str(cif_file), fmt='cif')
                
                downloaded[mp_id] = {
                    'formula': formula,
                    'num_sites': len(structure)
                }
                
                print(f"    {mp_id}: {formula} ({len(structure)} atoms)")
                
            except Exception as e:
                print(f"  ERROR: Failed to download {mp_id}: {e}")
    
    print(f"\nSuccessfully downloaded: {len(downloaded)} new structures")
    if skip_existing:
        print(f"Existing structures: {len(existing_mp_ids)}")
    print(f"Total structures: {len(downloaded) + (len(existing_mp_ids) if skip_existing else 0)}")
    
    return downloaded


def main():
    parser = argparse.ArgumentParser(
        description="Extract and download MP structures from cached JSON file"
    )
    parser.add_argument(
        '--cache-file',
        type=str,
        default='./VASP_JOBS/mp_mattersim.json',
        help="Path to mp_mattersim.json cache file (default: ./VASP_JOBS/mp_mattersim.json)"
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
    print("MP Structure Download")
    print("="*70)
    print(f"Cache file: {cache_file}")
    print(f"Output directory: {output_dir}")
    if args.chemsys:
        print(f"Chemical system filter: {args.chemsys}")
    else:
        print(f"Chemical system filter: All")
    print(f"Force re-download: {args.force}")
    print("="*70)
    
    # Load cache file
    if not cache_file.exists():
        print(f"ERROR: Cache file not found: {cache_file}")
        return 1
    
    print(f"\nLoading cache file...")
    with open(cache_file, 'r') as f:
        cache_data = json.load(f)
    print(f"Loaded {len(cache_data)} cached MP entries\n")
    
    # Extract MP IDs grouped by chemical system
    mp_ids_by_chemsys = extract_mp_ids_from_cache(cache_data, chemsys=args.chemsys)
    
    if not mp_ids_by_chemsys:
        print("ERROR: No MP IDs found in cache!")
        if args.chemsys:
            print(f"  (chemsys filter: {args.chemsys})")
        return 1
    
    # Collect all unique MP IDs across all chemical systems
    all_mp_ids = set()
    for ids in mp_ids_by_chemsys.values():
        all_mp_ids.update(ids)
    
    print(f"Found {len(mp_ids_by_chemsys)} chemical system(s):")
    for cs, ids in sorted(mp_ids_by_chemsys.items()):
        print(f"  {cs}: {len(ids)} structures")
    print(f"\nTotal unique MP IDs: {len(all_mp_ids)}")
    print()
    
    # Download all structures (flat, no subdirectories)
    print("="*70)
    print("Downloading Structures")
    print("="*70)
    
    try:
        downloaded = download_mp_structures(
            list(all_mp_ids), 
            output_dir, 
            mp_api_key, 
            skip_existing=not args.force
        )
    except Exception as e:
        print(f"ERROR: Failed to download structures: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    print("\n" + "="*70)
    print("Summary")
    print("="*70)
    print(f"Total unique MP IDs: {len(all_mp_ids)}")
    print(f"Output directory: {output_dir}")
    print(f"Structure files: {output_dir}/mp-*.cif")
    print("="*70)
    print("\nNote: MP energies are fetched separately by compare_mp_vasp_energies.py")
    print("="*70 + "\n")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

