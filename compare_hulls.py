#!/usr/bin/env python3
"""
Compare MatterSim pre-screening results vs DFT hull calculations.

Analyzes correlation between MatterSim and DFT energy_above_hull values
to validate pre-screening accuracy.
"""

import json
import argparse
from pathlib import Path
import numpy as np

def load_results(mattersim_json, dft_json):
    """Load both MatterSim and DFT results."""
    with open(mattersim_json, 'r') as f:
        mattersim_data = json.load(f)
    
    with open(dft_json, 'r') as f:
        dft_data = json.load(f)
    
    return mattersim_data, dft_data


def match_structures(mattersim_results, dft_results):
    """Match structures between MatterSim and DFT results."""
    # Create lookup dict for DFT results
    dft_lookup = {}
    for result in dft_results:
        if result['dft_energy_above_hull'] is not None:
            dft_lookup[result['structure_id']] = result['dft_energy_above_hull']
    
    # Match with MatterSim results
    matched = []
    for result in mattersim_results:
        struct_id = result['structure_id']
        if struct_id in dft_lookup:
            matched.append({
                'structure_id': struct_id,
                'composition': result['composition'],
                'mattersim_e_hull': result['energy_above_hull'],
                'dft_e_hull': dft_lookup[struct_id],
                'passed_prescreen': result['passed_prescreening']
            })
    
    return matched


def analyze_correlation(matched_data, threshold=0.1):
    """Analyze correlation and pre-screening accuracy."""
    mattersim_vals = np.array([d['mattersim_e_hull'] for d in matched_data])
    dft_vals = np.array([d['dft_e_hull'] for d in matched_data])
    passed = np.array([d['passed_prescreen'] for d in matched_data])
    
    # Calculate correlation
    correlation = np.corrcoef(mattersim_vals, dft_vals)[0, 1]
    mae = np.mean(np.abs(mattersim_vals - dft_vals))
    rmse = np.sqrt(np.mean((mattersim_vals - dft_vals)**2))
    
    # Pre-screening accuracy
    dft_stable = dft_vals < threshold
    
    true_positives = np.sum(passed & dft_stable)
    false_positives = np.sum(passed & ~dft_stable)
    true_negatives = np.sum(~passed & ~dft_stable)
    false_negatives = np.sum(~passed & dft_stable)
    
    precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    return {
        'correlation': correlation,
        'mae': mae,
        'rmse': rmse,
        'true_positives': int(true_positives),
        'false_positives': int(false_positives),
        'true_negatives': int(true_negatives),
        'false_negatives': int(false_negatives),
        'precision': precision,
        'recall': recall,
        'f1_score': f1_score
    }


def main():
    parser = argparse.ArgumentParser(
        description="Compare MatterSim pre-screening vs DFT hull results"
    )
    parser.add_argument(
        '--mattersim-json',
        type=str,
        default='./VASP_JOBS/prescreening_stability.json',
        help="MatterSim pre-screening results JSON"
    )
    parser.add_argument(
        '--dft-json',
        type=str,
        default='./VASP_JOBS/dft_stability_results.json',
        help="DFT hull results JSON"
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=0.1,
        help="E_hull threshold for stability (eV/atom)"
    )
    parser.add_argument(
        '--output',
        type=str,
        default='hull_comparison.json',
        help="Output comparison results"
    )
    
    args = parser.parse_args()
    
    mattersim_json = Path(args.mattersim_json).expanduser()
    dft_json = Path(args.dft_json).expanduser()
    output_path = Path(args.output).expanduser()
    
    print("="*70)
    print("MatterSim vs DFT Hull Comparison")
    print("="*70)
    print(f"MatterSim results: {mattersim_json}")
    print(f"DFT results: {dft_json}")
    print(f"Stability threshold: {args.threshold} eV/atom")
    print("="*70 + "\n")
    
    # Load results
    print("Loading results...")
    mattersim_data, dft_data = load_results(mattersim_json, dft_json)
    
    print(f"MatterSim structures: {len(mattersim_data['results'])}")
    print(f"DFT structures: {len(dft_data['results'])}")
    
    # Match structures
    print("\nMatching structures...")
    matched = match_structures(mattersim_data['results'], dft_data['results'])
    print(f"Matched structures: {len(matched)}")
    
    if len(matched) == 0:
        print("ERROR: No matching structures found")
        return 1
    
    # Analyze
    print("\nAnalyzing correlation and accuracy...")
    stats = analyze_correlation(matched, threshold=args.threshold)
    
    # Print results
    print("\n" + "="*70)
    print("Results")
    print("="*70)
    print(f"\nCorrelation & Error:")
    print(f"  Pearson correlation: {stats['correlation']:.4f}")
    print(f"  Mean Absolute Error: {stats['mae']:.4f} eV/atom")
    print(f"  RMSE: {stats['rmse']:.4f} eV/atom")
    
    print(f"\nPre-screening Accuracy (threshold = {args.threshold} eV/atom):")
    print(f"  True Positives:  {stats['true_positives']:4d} (correct: stable → passed)")
    print(f"  False Positives: {stats['false_positives']:4d} (error: unstable → passed)")
    print(f"  True Negatives:  {stats['true_negatives']:4d} (correct: unstable → filtered)")
    print(f"  False Negatives: {stats['false_negatives']:4d} (error: stable → filtered)")
    
    print(f"\nMetrics:")
    print(f"  Precision: {stats['precision']:.2%} (of passed structures, how many are stable)")
    print(f"  Recall:    {stats['recall']:.2%} (of stable structures, how many passed)")
    print(f"  F1 Score:  {stats['f1_score']:.4f}")
    
    # Compute cost savings
    total = len(matched)
    passed = stats['true_positives'] + stats['false_positives']
    filtered = total - passed
    savings = filtered / total if total > 0 else 0
    
    print(f"\nComputational Savings:")
    print(f"  Total structures: {total}")
    print(f"  Passed pre-screening: {passed} ({passed/total:.1%})")
    print(f"  Filtered: {filtered} ({savings:.1%} VASP computations saved)")
    
    # Identify false negatives (stable but filtered)
    if stats['false_negatives'] > 0:
        print(f"\nWarning: {stats['false_negatives']} stable structures were filtered!")
        print("False negatives (DFT stable but filtered by MatterSim):")
        false_neg_list = [d for d in matched if not d['passed_prescreen'] and d['dft_e_hull'] < args.threshold]
        for d in false_neg_list[:10]:  # Show first 10
            print(f"  {d['structure_id']:<25} DFT: {d['dft_e_hull']:.4f}, MatterSim: {d['mattersim_e_hull']:.4f}")
        if len(false_neg_list) > 10:
            print(f"  ... and {len(false_neg_list) - 10} more")
    
    print("="*70 + "\n")
    
    # Save detailed comparison
    output_data = {
        'summary': {
            'total_matched': len(matched),
            'threshold': args.threshold,
            'correlation': stats['correlation'],
            'mae': stats['mae'],
            'rmse': stats['rmse'],
            'precision': stats['precision'],
            'recall': stats['recall'],
            'f1_score': stats['f1_score'],
            'computational_savings': savings
        },
        'confusion_matrix': {
            'true_positives': stats['true_positives'],
            'false_positives': stats['false_positives'],
            'true_negatives': stats['true_negatives'],
            'false_negatives': stats['false_negatives']
        },
        'matched_structures': matched
    }
    
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"Detailed comparison saved to: {output_path}\n")
    
    return 0


if __name__ == '__main__':
    exit(main())

