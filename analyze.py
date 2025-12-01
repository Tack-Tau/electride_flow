#!/usr/bin/env python3
"""
Analyze completed ELF calculations for electride candidates.

Workflow:
1. Extract PARCHG-* files from PARCHG.tar.gz
2. Run Electride.py analysis
3. Clean up extracted PARCHG-* files
4. PARCHG.tar.gz is ALWAYS preserved for data traceability
5. Save electride candidates to PyXtal database

Features:
- Incremental analysis: skips structures already in CSV or database
- Adds e_above_hull from DFT calculations
- Adds spacegroup from PyXtal analysis of CONTCAR
- Saves to PyXtal database with adaptive tolerance

Note: This is much faster than re-compressing - we just extract, analyze, and delete extracted files.
"""

import sys
import os
import json
import tarfile
import shutil
import argparse
import traceback
import pandas as pd
from tabulate import tabulate
from pathlib import Path
from Electride import electride
from pyxtal import pyxtal
from pyxtal.db import database
from pymatgen.io.vasp.outputs import Poscar


def load_existing_results(csv_path):
    """Load existing CSV results if available."""
    if not os.path.exists(csv_path):
        return None, set()
    
    try:
        df = pd.read_csv(csv_path)
        analyzed_ids = set(df['formula'].values)
        print(f"Found existing results: {len(analyzed_ids)} structures already analyzed")
        return df, analyzed_ids
    except Exception as e:
        print(f"Warning: Could not load existing CSV: {e}")
        return None, set()


def load_existing_database(db_path):
    """Load existing PyXtal database if available."""
    if not os.path.exists(db_path):
        return None, set()
    
    try:
        db = database(db_path)
        db_ids = set()
        for entry in db.get_all():
            if 'structure_id' in entry:
                db_ids.add(entry['structure_id'])
        print(f"Found existing database: {len(db_ids)} structures already in database")
        return db, db_ids
    except Exception as e:
        print(f"Warning: Could not load existing database: {e}")
        return None, set()


def load_dft_stability(dft_json_path):
    """Load DFT stability results for e_above_hull values."""
    if not os.path.exists(dft_json_path):
        print(f"Warning: DFT stability file not found: {dft_json_path}")
        return {}
    
    try:
        with open(dft_json_path, 'r') as f:
            data = json.load(f)
        
        e_hull_map = {}
        for result in data.get('results', []):
            struct_id = result['structure_id']
            e_hull_map[struct_id] = result.get('energy_above_hull', None)
        
        print(f"Loaded DFT stability: {len(e_hull_map)} structures with e_above_hull")
        return e_hull_map
    except Exception as e:
        print(f"Warning: Could not load DFT stability file: {e}")
        return {}


def get_spacegroup_from_contcar(relax_dir):
    """Extract spacegroup from CONTCAR using PyXtal."""
    contcar_path = os.path.join(relax_dir, 'CONTCAR')
    if not os.path.exists(contcar_path):
        return None
    
    try:
        poscar = Poscar.from_file(contcar_path)
        struct = poscar.structure
        
        # Try PyXtal symmetrization with progressive tolerances
        tolerances = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.3, 0.5]
        for tol in tolerances:
            try:
                xtal = pyxtal()
                xtal.from_seed(struct, tol=tol)
                return xtal.group.number
            except:
                continue
        
        # If all symmetrization attempts failed, use P1 (no symmetry)
        return 1
    except Exception as e:
        return None


def load_workflow_db(db_path):
    """Load workflow database and extract ELF_DONE structures."""
    with open(db_path, 'r') as f:
        data = json.load(f)
    
    elf_done = []
    metals = []
    semiconductors = []
    
    for struct_id, sdata in data['structures'].items():
        if sdata['state'] == 'ELF_DONE' and sdata.get('elf_dir'):
            elf_dir = sdata['elf_dir']
            relax_dir = sdata.get('relax_dir', '')
            elf_done.append({
                'id': struct_id,
                'dir': elf_dir,
                'relax_dir': relax_dir,
                'composition': sdata.get('composition', ''),
                'is_semi': sdata.get('is_semiconductor', False)
            })
            if sdata.get('is_semiconductor'):
                semiconductors.append(struct_id)
            else:
                metals.append(struct_id)
    
    print(f"Total ELF_DONE: {len(elf_done)}")
    print(f"  Semiconductors (with PARCHG): {len(semiconductors)}")
    print(f"  Metals (ELFCAR only): {len(metals)}")
    print("")
    
    return elf_done


def extract_parchg_files(elf_dir):
    """
    Extract PARCHG files from PARCHG.tar.gz for analysis.
    
    Returns:
        bool: True if extraction was performed, False otherwise
    """
    elf_path = Path(elf_dir)
    tar_archive = elf_path / "PARCHG.tar.gz"
    
    # If no tar archive exists, check if files are already extracted
    if not tar_archive.exists():
        # Check if PARCHG-* files already exist
        if list(elf_path.glob("PARCHG-*")):
            return False  # Files already available, no extraction needed
        return False  # No archive and no files
    
    # Check if files are already extracted
    if list(elf_path.glob("PARCHG-*")):
        return False  # Already extracted
    
    # Extract tar archive
    try:
        with tarfile.open(tar_archive, 'r:gz') as tar:
            tar.extractall(path=elf_path)
        
        extracted_files = list(elf_path.glob("PARCHG-*"))
        print(f"  Extracted {len(extracted_files)} PARCHG files")
        return True  # Extraction was performed
    except Exception as e:
        print(f"  Warning: Could not extract {tar_archive.name}: {e}")
        return False


def cleanup_parchg_files(elf_dir, was_extracted):
    """
    Clean up extracted PARCHG-* files after analysis.
    
    IMPORTANT: Original PARCHG.tar.gz is NEVER deleted - only extracted PARCHG-* files are removed.
    This preserves the compressed archive for future data traceability.
    
    Args:
        elf_dir: ELF directory path
        was_extracted: If True, delete extracted files. If False, keep them.
    """
    if not was_extracted:
        return  # Don't delete if we didn't extract
    
    elf_path = Path(elf_dir)
    tar_archive = elf_path / "PARCHG.tar.gz"
    
    # Safety check: Verify tar.gz exists before deleting extracted files
    if not tar_archive.exists():
        print(f"  Warning: PARCHG.tar.gz not found, keeping extracted files for safety")
        return
    
    # Find all extracted PARCHG-* files (NOT the tar.gz!)
    parchg_files = sorted(elf_path.glob("PARCHG-*"))
    
    if not parchg_files:
        return  # Nothing to clean
    
    try:
        # Delete only extracted PARCHG-* files (tar.gz is preserved)
        for parchg_file in parchg_files:
            if parchg_file.is_file() and parchg_file.suffix != '.gz':  # Extra safety
                parchg_file.unlink()
        print(f"  Cleaned up {len(parchg_files)} extracted PARCHG files (tar.gz preserved)")
    except Exception as e:
        print(f"  Warning: Could not clean up PARCHG files: {e}")


def analyze_structure(struct_info, bader_cmd, threshold, keep_extracted=False):
    """
    Analyze a single structure for electride character.
    
    Args:
        struct_info: Dict with structure info (id, dir, relax_dir, composition, is_semi)
        bader_cmd: Bader executable command
        threshold: ELF threshold
        keep_extracted: If True, keep extracted PARCHG files. If False, delete after analysis.
    
    Returns:
        tuple: (struct_id, volumes, spacegroup, error)
    """
    struct_id = struct_info['id']
    elf_dir = struct_info['dir']
    relax_dir = struct_info.get('relax_dir', '')
    
    print(f"Analyzing: {struct_id} ({elf_dir})")
    
    # Check if ELFCAR exists
    elfcar_path = os.path.join(elf_dir, 'ELFCAR')
    if not os.path.exists(elfcar_path):
        print(f"  ERROR: ELFCAR not found at {elfcar_path}")
        return struct_id, [0, 0, 0, 0, 0], None, f"ELFCAR not found"
    
    # Extract PARCHG files from tar.gz archive
    was_extracted = extract_parchg_files(elf_dir)
    
    # Get spacegroup from CONTCAR
    spacegroup = None
    if relax_dir:
        spacegroup = get_spacegroup_from_contcar(relax_dir)
        if spacegroup:
            print(f"  Spacegroup: {spacegroup}")
    
    # Run electride analysis
    try:
        test = electride(elf_dir, cmd=bader_cmd, ELF_min=threshold)
        
        # Extract volumes for all 5 PARCHG types: e0025, e05, e10, band0, band1
        volumes = [float(test.volume[i]) for i in range(5)]
        
        # Clean up extracted files (delete them, tar.gz remains)
        if not keep_extracted:
            cleanup_parchg_files(elf_dir, was_extracted)
        
        # New electride criteria: ALL of e0025, e05, e10, band0 must be > 0
        is_electride = (volumes[0] > 0 and volumes[1] > 0 and 
                       volumes[2] > 0 and volumes[3] > 0)
        
        if test.error:
            print(f"  Warning: {test.error}")
            return struct_id, volumes, spacegroup, test.error
        elif is_electride:
            print(f"  *** ELECTRIDE CANDIDATE *** (e0025={volumes[0]:.2f}, e05={volumes[1]:.2f}, e10={volumes[2]:.2f}, band0={volumes[3]:.2f}, band1={volumes[4]:.2f})")
            return struct_id, volumes, spacegroup, None
        else:
            print(f"  Not electride (criteria: e0025 AND e05 AND e10 AND band0 > 0)")
            return struct_id, volumes, spacegroup, None
            
    except Exception as e:
        error_msg = str(e)
        print(f"  ERROR: {error_msg}")
        print(f"  Traceback:")
        traceback.print_exc()
        
        # Clean up on error
        if not keep_extracted:
            cleanup_parchg_files(elf_dir, was_extracted)
        
        return struct_id, [0, 0, 0, 0, 0], spacegroup, error_msg


def save_to_database(struct_id, relax_dir, composition, e_above_hull, is_electride, db, tolerances):
    """Save structure to PyXtal database with adaptive tolerance."""
    contcar_path = os.path.join(relax_dir, 'CONTCAR')
    if not os.path.exists(contcar_path):
        return False
    
    try:
        poscar = Poscar.from_file(contcar_path)
        struct = poscar.structure
        
        # Try progressive tolerances
        for tol in tolerances:
            try:
                xtal = pyxtal()
                xtal.from_seed(struct, tol=tol)
                db.add_xtal(
                    xtal,
                    kvp={
                        'structure_id': struct_id,
                        'e_above_hull': e_above_hull,
                        'composition': composition,
                        'symmetrized': True,
                        'is_electride': is_electride
                    }
                )
                return True
            except:
                continue
        
        # If all tolerances failed, use P1 (no symmetry) as final fallback
        try:
            xtal = pyxtal()
            xtal.from_seed(struct, tol=1.0)
            db.add_xtal(
                xtal,
                kvp={
                    'structure_id': struct_id,
                    'e_above_hull': e_above_hull,
                    'composition': composition,
                    'symmetrized': True,
                    'is_electride': is_electride,
                    'spacegroup': 1  # P1 fallback
                }
            )
            print(f"    Using P1 fallback for {struct_id}")
            return True
        except:
            pass
        
        return False
    except Exception as e:
        print(f"  Warning: Could not save {struct_id} to database: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Analyze completed ELF calculations for electride candidates (incremental mode supported)"
    )
    parser.add_argument(
        '--db',
        type=str,
        required=True,
        help="Path to workflow database (workflow.json)"
    )
    parser.add_argument(
        '--bader-exe',
        type=str,
        default='bader',
        help="Path to bader executable (default: bader)"
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=0.6,
        help="ELF threshold (default: 0.6)"
    )
    parser.add_argument(
        '--output',
        type=str,
        default='electride_analysis.csv',
        help="Output CSV file (default: electride_analysis.csv)"
    )
    parser.add_argument(
        '--pyxtal-db',
        type=str,
        default='electride_data.db',
        help="PyXtal database file (default: electride_data.db)"
    )
    parser.add_argument(
        '--dft-stability',
        type=str,
        default='./VASP_JOBS/dft_stability_results.json',
        help="DFT stability results JSON (default: ./VASP_JOBS/dft_stability_results.json)"
    )
    parser.add_argument(
        '--keep-extracted',
        action='store_true',
        help="Keep extracted PARCHG-* files after analysis (default: delete to save space). "
             "Note: PARCHG.tar.gz is always preserved regardless of this flag."
    )
    
    args = parser.parse_args()
    
    # Check if database exists
    if not os.path.exists(args.db):
        print(f"ERROR: Workflow database not found: {args.db}")
        sys.exit(1)
    
    print("="*70)
    print("Electride Analysis - Completed ELF Calculations")
    print("="*70)
    print(f"Workflow database: {args.db}")
    print(f"Bader executable: {args.bader_exe}")
    print(f"ELF threshold: {args.threshold}")
    print(f"Output CSV: {args.output}")
    print(f"PyXtal database: {args.pyxtal_db}")
    print(f"DFT stability: {args.dft_stability}")
    print()
    print("Electride Criteria: e0025 > 0 AND e05 > 0 AND e10 > 0 AND band0 > 0")
    print()
    print("Workflow: Extract → Analyze → Cleanup → Save to DB")
    print("  - PARCHG.tar.gz is always preserved")
    print("  - Only temporary extracted files are deleted")
    print("  - Incremental mode: skips already-analyzed structures")
    print("="*70)
    print("")
    
    # Load existing results (CSV and PyXtal database)
    existing_df, analyzed_csv_ids = load_existing_results(args.output)
    existing_pyxtal_db, analyzed_db_ids = load_existing_database(args.pyxtal_db)
    
    # Combined set of already-processed structures
    skip_ids = analyzed_csv_ids | analyzed_db_ids
    
    # Load DFT stability results
    e_hull_map = load_dft_stability(args.dft_stability)
    
    # Load workflow database
    elf_structures = load_workflow_db(args.db)
    
    if not elf_structures:
        print("No completed ELF calculations found.")
        sys.exit(0)
    
    # Filter out already-analyzed structures
    structures_to_analyze = [s for s in elf_structures if s['id'] not in skip_ids]
    
    print(f"Found {len(elf_structures)} ELF_DONE structures")
    print(f"  Already analyzed: {len(skip_ids)}")
    print(f"  New to analyze: {len(structures_to_analyze)}")
    print("")
    
    if len(structures_to_analyze) == 0:
        print("No new structures to analyze. All up to date!")
        if existing_df is not None:
            print(f"\nExisting results in {args.output}:")
            print(tabulate(existing_df.head(20), headers='keys', tablefmt='psql', showindex=False))
        sys.exit(0)
    
    # Initialize or open PyXtal database
    if existing_pyxtal_db:
        pyxtal_db = existing_pyxtal_db
    else:
        pyxtal_db = database(args.pyxtal_db)
    
    # Analyze new structures
    bader_cmd = args.bader_exe + ' '
    new_results = {
        'formula': [],
        'composition': [],
        'e0025': [],
        'e05': [],
        'e10': [],
        'band0': [],
        'band1': [],
        'spacegroup': [],
        'e_above_hull': []
    }
    
    tolerances = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.3, 0.5]
    db_saved_count = 0
    db_failed_count = 0
    
    for struct_info in structures_to_analyze:
        struct_id, volumes, spacegroup, error = analyze_structure(
            struct_info, bader_cmd, args.threshold, 
            keep_extracted=args.keep_extracted
        )
        
        # Get e_above_hull from DFT stability
        e_hull = e_hull_map.get(struct_id, None)
        
        # Determine if electride: ALL of e0025, e05, e10, band0 > 0
        is_electride = (volumes[0] > 0 and volumes[1] > 0 and 
                       volumes[2] > 0 and volumes[3] > 0)
        
        # Save to results
        new_results['formula'].append(struct_id)
        new_results['composition'].append(struct_info.get('composition', ''))
        new_results['e0025'].append(volumes[0])
        new_results['e05'].append(volumes[1])
        new_results['e10'].append(volumes[2])
        new_results['band0'].append(volumes[3])
        new_results['band1'].append(volumes[4])
        new_results['spacegroup'].append(spacegroup if spacegroup else 0)
        new_results['e_above_hull'].append(e_hull if e_hull is not None else float('inf'))
        
        # Save to PyXtal database
        relax_dir = struct_info.get('relax_dir', '')
        if relax_dir:
            saved = save_to_database(
                struct_id, relax_dir, 
                struct_info.get('composition', ''),
                e_hull, is_electride, 
                pyxtal_db, tolerances
            )
            if saved:
                db_saved_count += 1
                print(f"  Saved to database: {struct_id}")
            else:
                db_failed_count += 1
                print(f"  Failed to save to database: {struct_id}")
    
    # Commit database
    if hasattr(pyxtal_db, 'db') and hasattr(pyxtal_db.db, 'commit'):
        pyxtal_db.db.commit()
    
    # Create DataFrame for new results
    new_df = pd.DataFrame(new_results)
    
    # Merge with existing results if any
    if existing_df is not None:
        df = pd.concat([existing_df, new_df], ignore_index=True)
    else:
        df = new_df
    
    # Sort by e_above_hull (low to high)
    df = df.sort_values(by='e_above_hull', ascending=True).reset_index(drop=True)
    
    # Save updated CSV
    df.to_csv(args.output, index=False)
    
    print("")
    print("="*70)
    print("Summary Results")
    print("="*70)
    print(f"Total structures in CSV: {len(df)}")
    print(f"New structures analyzed: {len(new_df)}")
    print(f"PyXtal database: {db_saved_count} saved, {db_failed_count} failed")
    print("")
    print(f"Results saved to: {args.output}")
    print(f"Database saved to: {args.pyxtal_db}")
    print("")
    
    # Show top results (sorted by e_above_hull)
    print("Top 20 structures (sorted by e_above_hull):")
    print(tabulate(df.head(20), headers='keys', tablefmt='psql', showindex=False))
    
    # Summary statistics (electride: ALL of e0025, e05, e10, band0 > 0)
    electrides = df[(df['e0025'] > 0) & (df['e05'] > 0) & (df['e10'] > 0) & (df['band0'] > 0)]
    print("")
    print("="*70)
    print("Statistics")
    print("="*70)
    print(f"Total structures analyzed: {len(df)}")
    print(f"Electride candidates: {len(electrides)} ({100*len(electrides)/len(df):.1f}%)")
    print(f"  (Criteria: e0025 > 0 AND e05 > 0 AND e10 > 0 AND band0 > 0)")
    
    if len(electrides) > 0:
        print("")
        print("Top electride candidates (sorted by e_above_hull):")
        print(tabulate(
            electrides.head(10),
            headers='keys',
            tablefmt='psql',
            showindex=False
        ))
    
    print("="*70)


if __name__ == '__main__':
    main()

