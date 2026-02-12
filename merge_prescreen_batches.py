#!/usr/bin/env python3
"""
Merge parallel prescreening batch results into a single file.

Usage:
    python3 merge_prescreen_batches.py --output-dir ./VASP_JOBS
    python3 merge_prescreen_batches.py --output-dir ./VASP_JOBS --clean-batches
    python3 merge_prescreen_batches.py --output-dir ./VASP_JOBS --skip-database
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Any
from collections import defaultdict

try:
    from pyxtal import pyxtal
    from pyxtal.db import database_topology
    PYXTAL_AVAILABLE = True
except ImportError:
    PYXTAL_AVAILABLE = False

try:
    from pymatgen.analysis.structure_matcher import StructureMatcher
    from pymatgen.io.ase import AseAtomsAdaptor
    PYMATGEN_AVAILABLE = True
except ImportError:
    PYMATGEN_AVAILABLE = False


def load_batch_results(output_dir: Path) -> List[Dict[str, Any]]:
    """
    Load all prescreening_stability_batch*.json files.
    
    Args:
        output_dir: Directory containing batch files
        
    Returns:
        List of batch result dictionaries
    """
    batch_files = sorted(output_dir.glob('prescreening_stability_batch*.json'))
    
    if not batch_files:
        print(f"ERROR: No batch files found in {output_dir}")
        print("Expected files matching: prescreening_stability_batch*.json")
        sys.exit(1)
    
    print(f"Found {len(batch_files)} batch files:")
    for bf in batch_files:
        print(f"  - {bf.name}")
    print()
    
    batch_results = []
    for bf in batch_files:
        try:
            with open(bf, 'r') as f:
                data = json.load(f)
                batch_results.append(data)
        except Exception as e:
            print(f"WARNING: Failed to load {bf.name}: {e}")
            continue
    
    return batch_results


def merge_batch_results(batch_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Merge batch results into a single result dictionary.
    
    Args:
        batch_results: List of batch result dictionaries
        
    Returns:
        Merged result dictionary
    """
    if not batch_results:
        print("ERROR: No valid batch results to merge")
        sys.exit(1)
    
    # Use first batch as template
    merged = {
        'summary': batch_results[0]['summary'].copy(),
        'results': []
    }
    
    # Track unique structures by composition + structure_id
    seen_structures = set()
    duplicate_count = 0
    
    for batch_data in batch_results:
        for result in batch_data.get('results', []):
            composition = result['composition']
            struct_id = result['structure_id']
            key = (composition, struct_id)
            
            if key in seen_structures:
                duplicate_count += 1
                continue
            
            seen_structures.add(key)
            merged['results'].append(result)
    
    # Sort by composition, then structure_id
    merged['results'].sort(key=lambda x: (x['composition'], x['structure_id']))
    
    # Aggregate summary statistics from all batches
    total_loaded = sum(b['summary'].get('total_structures_loaded', 0) for b in batch_results)
    total_duplicates = sum(b['summary'].get('duplicate_structures_removed', 0) for b in batch_results)
    total_unique = sum(b['summary'].get('unique_structures_processed', 0) for b in batch_results)
    
    # Count structures that passed prescreening
    passed_count = sum(1 for r in merged['results'] if r.get('passed_prescreening', False))
    failed_count = len(merged['results']) - passed_count
    
    # Count duplicates in merged results (structures with 'duplicate_of' field)
    duplicate_records_count = sum(1 for r in merged['results'] if 'duplicate_of' in r)
    
    # Update summary
    merged['summary']['total_batches'] = len(batch_results)
    merged['summary']['total_structures_loaded'] = total_loaded
    merged['summary']['duplicate_structures_removed'] = total_duplicates
    merged['summary']['unique_structures_processed'] = total_unique
    merged['summary']['total_structures'] = len(merged['results'])  # Keep for compatibility
    merged['summary']['passed_prescreening'] = passed_count
    merged['summary']['failed_prescreening'] = failed_count
    
    print(f"Merged {len(batch_results)} batches:")
    print(f"  Total structures loaded: {total_loaded}")
    print(f"  Pre-relaxation duplicates removed: {total_duplicates} ({duplicate_records_count} duplicate records in merged results)")
    print(f"  Unique structures processed: {total_unique}")
    print(f"  Total result records: {len(merged['results'])}")
    print(f"  Cross-batch duplicates skipped: {duplicate_count}")
    print(f"  Passed prescreening: {passed_count}")
    print(f"  Failed prescreening: {failed_count}")
    
    return merged


def check_mp_cache(output_dir: Path) -> bool:
    """
    Check if shared MP cache exists (no merging needed - batches use shared file).
    
    Args:
        output_dir: Directory containing MP cache
        
    Returns:
        bool: True if cache exists
    """
    mp_cache = output_dir / 'mp_mattersim.json'
    
    if not mp_cache.exists():
        print("\nNo MP cache file found (mp_mattersim.json)")
        print("  Note: Parallel batches share a single MP cache with file locking")
        return False
    
    print(f"\n{'='*70}")
    print("Checking Shared MP MatterSim Cache")
    print(f"{'='*70}")
    
    try:
        with open(mp_cache, 'r') as f:
            cache_data = json.load(f)
        
        print(f"MP cache: {mp_cache}")
        print(f"  Total entries: {len(cache_data)}")
        
        # Show breakdown by chemical system
        chemsys_counts = {}
        for entry in cache_data:
            chemsys = entry.get('chemsys', 'unknown')
            chemsys_counts[chemsys] = chemsys_counts.get(chemsys, 0) + 1
        
        print(f"  Breakdown by chemical system:")
        for chemsys in sorted(chemsys_counts.keys())[:10]:  # Show first 10
            print(f"    {chemsys}: {chemsys_counts[chemsys]} phases")
        if len(chemsys_counts) > 10:
            print(f"    ... and {len(chemsys_counts) - 10} more systems")
        
        return True
    except Exception as e:
        print(f"  ERROR: Could not read MP cache: {e}")
        return False


def merge_pyxtal_databases(output_dir: Path) -> bool:
    """
    Merge PyXtal database files from parallel batches using PyXtal/ASE API.
    
    Reads structures from each batch database using db.db.select() and writes
    them into a new merged database using db.add_xtal() or db.db.write().
    This avoids raw SQL INSERT OR IGNORE issues with auto-increment IDs.
    
    Args:
        output_dir: Directory containing batch database files
        
    Returns:
        bool: True if successful, False otherwise
    """
    if not PYXTAL_AVAILABLE:
        print("\nWARNING: PyXtal not available. Cannot merge database files.")
        return False
    
    batch_dbs = sorted(output_dir.glob('prescreening_structures_batch*.db'))
    
    if not batch_dbs:
        print("\nNo batch database files found to merge")
        return False
    
    print(f"\n{'='*70}")
    print("Merging PyXtal Database Files")
    print(f"{'='*70}")
    print(f"Found {len(batch_dbs)} database files:")
    for db_path in batch_dbs:
        print(f"  - {db_path.name}")
    print()
    
    merged_db_path = output_dir / 'prescreening_structures.db'
    
    # Remove existing merged database if it exists
    if merged_db_path.exists():
        merged_db_path.unlink()
        print(f"Removed existing merged database\n")
    
    try:
        # Create fresh merged database
        merged_db = database_topology(str(merged_db_path))
        
        # Track structure_ids to avoid duplicates across batches
        seen_structure_ids = set()
        total_added = 0
        total_skipped = 0
        tolerances = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
        
        for db_path in batch_dbs:
            print(f"Reading {db_path.name}...")
            
            try:
                batch_db = database_topology(str(db_path))
                batch_count = batch_db.db.count()
                print(f"  Contains {batch_count} structures")
            except Exception as e:
                print(f"  ERROR: Could not open {db_path.name}: {e}")
                continue
            
            batch_added = 0
            batch_skipped = 0
            
            for row in batch_db.db.select():
                # Get structure_id
                try:
                    struct_id = row.structure_id
                except AttributeError:
                    # Fallback: try to get from key-value pairs
                    struct_id = row.get('structure_id', f"unknown_{row.id}")
                
                # Skip if already seen (cross-batch duplicate)
                if struct_id in seen_structure_ids:
                    batch_skipped += 1
                    continue
                
                seen_structure_ids.add(struct_id)
                
                # Extract atoms and all key-value pairs from the row
                try:
                    atoms = row.toatoms()
                except Exception as e:
                    print(f"    Warning: Could not extract atoms for {struct_id}: {e}")
                    batch_skipped += 1
                    continue
                
                # Collect all key-value pairs stored in the row
                kvp = {}
                for key in row.__dict__:
                    if key.startswith('_') or key in ('id', 'ctime', 'mtime', 'user',
                                                       'numbers', 'positions', 'cell',
                                                       'pbc', 'calculator', 'calculator_parameters',
                                                       'key_value_pairs', 'data', 'unique_id',
                                                       'constraints', 'natoms', 'fmax',
                                                       'smax', 'mass', 'charge', 'dipole',
                                                       'magmom', 'magmoms', 'stress',
                                                       'stresses', 'forces', 'free_energy',
                                                       'energy', 'volume'):
                        continue
                    try:
                        val = getattr(row, key)
                        if val is not None:
                            kvp[key] = val
                    except Exception:
                        continue
                
                # Also get key_value_pairs dict if available
                if hasattr(row, 'key_value_pairs') and row.key_value_pairs:
                    for k, v in row.key_value_pairs.items():
                        if k not in kvp:
                            kvp[k] = v
                
                # Try to reconstruct PyXtal object and add via add_xtal
                added = False
                for tol in tolerances:
                    try:
                        xtal = pyxtal()
                        xtal.from_seed(atoms, tol=tol)
                        if not xtal.valid:
                            continue
                        if len(xtal.check_short_distances(r=0.5)) > 0:
                            continue
                        merged_db.add_xtal(xtal, kvp=kvp)
                        added = True
                        break
                    except Exception:
                        continue
                
                if not added:
                    # Fallback: write directly as ASE atoms with metadata
                    try:
                        merged_db.db.write(atoms, **kvp)
                        added = True
                    except Exception as e:
                        print(f"    Warning: Could not save {struct_id}: {e}")
                
                if added:
                    batch_added += 1
                else:
                    batch_skipped += 1
            
            total_added += batch_added
            total_skipped += batch_skipped
            print(f"  Added {batch_added} structures, skipped {batch_skipped}")
        
        print(f"\nDatabase merge complete!")
        print(f"  Total structures in merged database: {total_added}")
        print(f"  Skipped (duplicates/errors): {total_skipped}")
        print(f"  Merged database: {merged_db_path}")
        return True
        
    except Exception as e:
        print(f"\nERROR merging databases: {e}")
        import traceback
        traceback.print_exc()
        return False


def deduplicate_relaxed_structures(output_dir: Path, merged_json: Dict[str, Any]) -> Dict[str, Any]:
    """
    Post-relaxation deduplication: remove structures that became identical
    after MatterSim relaxation using StructureMatcher.
    
    Only compares structures within the same composition that passed prescreening.
    Updates both the merged JSON and the merged database.
    
    Args:
        output_dir: Directory containing the merged database
        merged_json: Merged JSON result dictionary
        
    Returns:
        Updated merged JSON with post-relaxation duplicates marked
    """
    if not PYMATGEN_AVAILABLE:
        print("\nWARNING: pymatgen not available. Skipping post-relaxation deduplication.")
        return merged_json
    
    if not PYXTAL_AVAILABLE:
        print("\nWARNING: PyXtal not available. Skipping post-relaxation deduplication.")
        return merged_json
    
    merged_db_path = output_dir / 'prescreening_structures.db'
    if not merged_db_path.exists():
        print("\nWARNING: Merged database not found. Skipping post-relaxation deduplication.")
        return merged_json
    
    print(f"\n{'='*70}")
    print("Post-Relaxation Deduplication (StructureMatcher)")
    print(f"{'='*70}")
    print(f"Tolerances: ltol=0.2, stol=0.2, angle_tol=5")
    print(f"Comparing relaxed structures within same composition")
    print()
    
    matcher = StructureMatcher(ltol=0.2, stol=0.2, angle_tol=5)
    adaptor = AseAtomsAdaptor()
    
    # Load structures from merged database, grouped by composition
    try:
        db = database_topology(str(merged_db_path))
    except Exception as e:
        print(f"ERROR: Could not open merged database: {e}")
        return merged_json
    
    # Build mapping: structure_id -> {atoms, composition, passed_prescreening, row_data}
    struct_data = {}
    for row in db.db.select():
        try:
            struct_id = row.structure_id
            atoms = row.toatoms()
            composition = getattr(row, 'composition', 'unknown')
            passed = getattr(row, 'passed_prescreening', False)
            struct_data[struct_id] = {
                'atoms': atoms,
                'composition': composition,
                'passed': passed
            }
        except Exception:
            continue
    
    print(f"Loaded {len(struct_data)} structures from merged database")
    
    # Group passed structures by composition
    by_composition = defaultdict(list)
    for struct_id, data in struct_data.items():
        if data['passed']:
            by_composition[data['composition']].append(struct_id)
    
    print(f"Compositions with passed structures: {len(by_composition)}")
    
    # Find duplicates within each composition
    post_relax_duplicates = []  # list of (duplicate_id, original_id, composition)
    
    for comp, struct_ids in sorted(by_composition.items()):
        if len(struct_ids) < 2:
            continue
        
        # Convert to pymatgen structures for comparison
        pmg_structures = []
        valid_ids = []
        for sid in struct_ids:
            try:
                pmg_struct = adaptor.get_structure(struct_data[sid]['atoms'])
                pmg_structures.append(pmg_struct)
                valid_ids.append(sid)
            except Exception:
                continue
        
        if len(pmg_structures) < 2:
            continue
        
        # Compare all pairs, mark later structures as duplicates of earlier ones
        is_duplicate = {}  # struct_id -> original_id
        for i in range(len(pmg_structures)):
            if valid_ids[i] in is_duplicate:
                continue
            for j in range(i + 1, len(pmg_structures)):
                if valid_ids[j] in is_duplicate:
                    continue
                try:
                    if matcher.fit(pmg_structures[i], pmg_structures[j]):
                        is_duplicate[valid_ids[j]] = valid_ids[i]
                        post_relax_duplicates.append((valid_ids[j], valid_ids[i], comp))
                except Exception:
                    continue
    
    if not post_relax_duplicates:
        print("\nNo post-relaxation duplicates found")
        return merged_json
    
    print(f"\nFound {len(post_relax_duplicates)} post-relaxation duplicates:")
    
    # Group by composition for reporting
    dup_by_comp = defaultdict(list)
    for dup_id, orig_id, comp in post_relax_duplicates:
        dup_by_comp[comp].append((dup_id, orig_id))
    
    for comp, dups in sorted(dup_by_comp.items()):
        print(f"  {comp}: {len(dups)} duplicates")
        for dup_id, orig_id in dups[:3]:
            print(f"    - {dup_id} (duplicate of {orig_id})")
        if len(dups) > 3:
            print(f"    ... and {len(dups) - 3} more")
    
    # Update merged JSON: mark post-relaxation duplicates
    dup_set = {dup_id for dup_id, _, _ in post_relax_duplicates}
    dup_map = {dup_id: orig_id for dup_id, orig_id, _ in post_relax_duplicates}
    
    passed_before = sum(1 for r in merged_json['results'] if r.get('passed_prescreening', False))
    
    for result in merged_json['results']:
        if result['structure_id'] in dup_set:
            result['passed_prescreening'] = False
            result['post_relax_duplicate_of'] = dup_map[result['structure_id']]
            if result.get('error') is None:
                result['error'] = f"Post-relaxation duplicate of {dup_map[result['structure_id']]}"
    
    passed_after = sum(1 for r in merged_json['results'] if r.get('passed_prescreening', False))
    failed_after = len(merged_json['results']) - passed_after
    
    # Update summary
    merged_json['summary']['post_relax_duplicates_removed'] = len(post_relax_duplicates)
    merged_json['summary']['passed_prescreening'] = passed_after
    merged_json['summary']['failed_prescreening'] = failed_after
    
    print(f"\nUpdated JSON results:")
    print(f"  Passed before dedup: {passed_before}")
    print(f"  Post-relaxation duplicates: {len(post_relax_duplicates)}")
    print(f"  Passed after dedup: {passed_after}")
    
    # Update merged database: mark duplicates as not passed
    try:
        merged_db = database_topology(str(merged_db_path))
        updated_count = 0
        for row in merged_db.db.select():
            try:
                struct_id = row.structure_id
                if struct_id in dup_set:
                    merged_db.db.update(row.id,
                                        passed_prescreening=False,
                                        post_relax_duplicate_of=dup_map[struct_id])
                    updated_count += 1
            except Exception:
                continue
        print(f"  Updated {updated_count} entries in merged database")
    except Exception as e:
        print(f"  WARNING: Could not update database: {e}")
    
    return merged_json


def main():
    parser = argparse.ArgumentParser(
        description="Merge parallel prescreening batch results"
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='./VASP_JOBS',
        help="Directory containing batch result files"
    )
    parser.add_argument(
        '--clean-batches',
        action='store_true',
        help="Delete batch files after successful merge (default: keep them)"
    )
    parser.add_argument(
        '--output-file',
        type=str,
        default=None,
        help="Output filename (default: prescreening_stability.json)"
    )
    parser.add_argument(
        '--skip-database',
        action='store_true',
        help="Skip merging PyXtal database files (only merge JSON)"
    )
    parser.add_argument(
        '--skip-dedup',
        action='store_true',
        help="Skip post-relaxation deduplication (StructureMatcher)"
    )
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir).expanduser()
    
    if not output_dir.exists():
        print(f"ERROR: Output directory not found: {output_dir}")
        sys.exit(1)
    
    print("="*70)
    print("Merging Prescreening Batch Results")
    print("="*70)
    print(f"Output directory: {output_dir}\n")
    
    # Load all batch results
    batch_results = load_batch_results(output_dir)
    
    # Merge results
    merged_data = merge_batch_results(batch_results)
    
    # Check shared MP cache (no merging needed - file is shared across batches)
    print()
    check_mp_cache(output_dir)
    
    # Merge PyXtal databases unless skipped
    db_success = False
    if not args.skip_database:
        db_success = merge_pyxtal_databases(output_dir)
        if not db_success:
            print("\nWARNING: Database merge failed or was skipped")
            print("  Batch database files will be kept for manual inspection")
    else:
        print("\nSkipping database merge (--skip-database)")
        db_success = True  # Consider it successful if skipped
    
    # Post-relaxation deduplication
    if not args.skip_dedup and db_success and not args.skip_database:
        merged_data = deduplicate_relaxed_structures(output_dir, merged_data)
    elif args.skip_dedup:
        print("\nSkipping post-relaxation deduplication (--skip-dedup)")
    
    # Write merged JSON file
    output_filename = args.output_file or 'prescreening_stability.json'
    output_file = output_dir / output_filename
    
    with open(output_file, 'w') as f:
        json.dump(merged_data, f, indent=2)
    
    print(f"\nMerged results saved to: {output_file}")
    
    # Cleanup batch files only if explicitly requested
    if args.clean_batches:
        if not db_success:
            print("\nBatch files KEPT due to database merge failure")
            print("  Fix the database issues first, then re-run with --clean-batches")
        else:
            print("\nCleaning up batch files (--clean-batches)...")
            batch_files = list(output_dir.glob('prescreening_stability_batch*.json'))
            checkpoint_files = list(output_dir.glob('prescreening_checkpoint_batch*.json'))
            db_files = list(output_dir.glob('prescreening_structures_batch*.db'))
            mp_lock_files = list(output_dir.glob('mp_mattersim.lock'))
            chemsys_lock_files = list(output_dir.glob('mp_cache_*.lock'))
            
            all_files = batch_files + checkpoint_files + db_files + mp_lock_files + chemsys_lock_files
            for bf in all_files:
                try:
                    bf.unlink()
                    print(f"  Deleted: {bf.name}")
                except Exception as e:
                    print(f"  WARNING: Failed to delete {bf.name}: {e}")
            
            print(f"Deleted {len(all_files)} batch and checkpoint files")
    else:
        print("\nBatch files kept (use --clean-batches to delete after verifying merge)")
    
    print()
    print("="*70)
    print("Merge complete!")
    print("="*70)


if __name__ == '__main__':
    main()
