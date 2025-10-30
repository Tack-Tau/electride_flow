#!/usr/bin/env python3
"""
Analyze completed ELF calculations for electride candidates.
Reads workflow database and calls Electride.py for batch analysis.
"""

import sys
import os
import json
import argparse
import traceback
import pandas as pd
from tabulate import tabulate
from pathlib import Path
from Electride import electride


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
            elf_done.append({
                'id': struct_id,
                'dir': elf_dir,
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


def analyze_structure(struct_info, bader_cmd, threshold):
    """Analyze a single structure for electride character."""
    struct_id = struct_info['id']
    elf_dir = struct_info['dir']
    
    print(f"Analyzing: {struct_id} ({elf_dir})")
    
    # Check if ELFCAR exists
    elfcar_path = os.path.join(elf_dir, 'ELFCAR')
    if not os.path.exists(elfcar_path):
        print(f"  ERROR: ELFCAR not found at {elfcar_path}")
        return struct_id, [0, 0, 0, 0, 0], f"ELFCAR not found"
    
    # Run electride analysis
    try:
        test = electride(elf_dir, cmd=bader_cmd, ELF_min=threshold)
        
        # Extract volumes for all 5 PARCHG types: e0025, e05, e10, band0, band1
        # Convert to float for clean output (removes np.float64 wrapping)
        volumes = [float(test.volume[i]) for i in range(5)]
        
        if test.error:
            print(f"  Warning: {test.error}")
            return struct_id, volumes, test.error
        elif any(v > 0 for v in volumes):
            # Pretty print with named values
            print(f"  *** ELECTRIDE CANDIDATE *** (e0025={volumes[0]:.2f}, e05={volumes[1]:.2f}, e10={volumes[2]:.2f}, band0={volumes[3]:.2f}, band1={volumes[4]:.2f})")
            return struct_id, volumes, None
        else:
            print(f"  Not electride (no interstitial electrons)")
            return struct_id, volumes, None
            
    except Exception as e:
        error_msg = str(e)
        print(f"  ERROR: {error_msg}")
        print(f"  Traceback:")
        traceback.print_exc()
        return struct_id, [0, 0, 0, 0, 0], error_msg


def main():
    parser = argparse.ArgumentParser(
        description="Analyze completed ELF calculations for electride candidates"
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
    
    args = parser.parse_args()
    
    # Check if database exists
    if not os.path.exists(args.db):
        print(f"ERROR: Workflow database not found: {args.db}")
        sys.exit(1)
    
    print("="*70)
    print("Electride Analysis - Completed ELF Calculations")
    print("="*70)
    print(f"Database: {args.db}")
    print(f"Bader executable: {args.bader_exe}")
    print(f"ELF threshold: {args.threshold}")
    print(f"Output file: {args.output}")
    print("="*70)
    print("")
    
    # Load workflow database
    elf_structures = load_workflow_db(args.db)
    
    if not elf_structures:
        print("No completed ELF calculations found.")
        sys.exit(0)
    
    # Analyze each structure
    bader_cmd = args.bader_exe + ' '
    results = {
        'formula': [],
        'e0025': [],
        'e05': [],
        'e10': [],
        'band0': [],
        'band1': []
    }
    
    for struct_info in elf_structures:
        struct_id, volumes, error = analyze_structure(
            struct_info, bader_cmd, args.threshold
        )
        
        results['formula'].append(struct_id)
        results['e0025'].append(volumes[0])
        results['e05'].append(volumes[1])
        results['e10'].append(volumes[2])
        results['band0'].append(volumes[3])
        results['band1'].append(volumes[4])
    
    # Create DataFrame and save
    df = pd.DataFrame(results)
    
    print("")
    print("="*70)
    print("Summary Results")
    print("="*70)
    print(tabulate(df, headers='keys', tablefmt='psql', showindex=True))
    
    df.to_csv(args.output, index=False)
    print("")
    print(f"Results saved to: {args.output}")
    
    # Summary statistics (electride if ANY volume > 0)
    electrides = df[(df['e0025'] > 0) | (df['e05'] > 0) | (df['e10'] > 0) | (df['band0'] > 0) | (df['band1'] > 0)]
    print("")
    print("="*70)
    print("Statistics")
    print("="*70)
    print(f"Total structures analyzed: {len(df)}")
    print(f"Electride candidates: {len(electrides)} ({100*len(electrides)/len(df):.1f}%)")
    
    if len(electrides) > 0:
        print("")
        print("Top candidates (sorted by e0025):")
        print(tabulate(
            electrides.sort_values('e0025', ascending=False).head(10),
            headers='keys',
            tablefmt='psql',
            showindex=False
        ))
    
    print("="*70)


if __name__ == '__main__':
    main()

