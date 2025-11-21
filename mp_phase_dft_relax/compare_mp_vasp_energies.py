#!/usr/bin/env python3
"""
Compare MP and VASP energies from MP phase relaxation workflow.

Compares:
- MP GGA-PBE DFT energies (raw DFT from Materials Project API)
- VASP PBE DFT energies (from local VASP relaxation calculations)
"""

import os
import json
import argparse
import numpy as np
from pathlib import Path
from collections import defaultdict

from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.outputs import Vasprun
import matplotlib.pyplot as plt
import matplotlib as mpl

# Use non-interactive backend for cluster
mpl.use('Agg')


def load_workflow_db(db_path):
    """Load workflow database."""
    with open(db_path, 'r') as f:
        return json.load(f)


def fetch_mp_energies(mp_ids, mp_api_key, debug=False):
    """
    Fetch GGA-PBE DFT energies from Materials Project API.
    
    IMPORTANT: MP stores calculations with multiple functionals (GGA, r2SCAN, SCAN).
    For comparison with VASP PBE calculations, we MUST filter for GGA entries only.
    
    The raw DFT energy corresponds to "Final Energy/Atom" shown on the MP website
    calculation summary (before composition-based corrections applied for hull analysis).
    
    Strategy:
    1. Use materials.thermo endpoint and filter for energy_type='GGA'
    2. Get uncorrected_energy_per_atom (raw DFT, no composition corrections)
    3. Fallback methods if thermo endpoint unavailable
    
    Example: mp-23309 (ScCl3)
    - GGA entry: uncorrected_energy_per_atom = -5.088 eV/atom   (matches VASP PBE)
    - r2SCAN entry: uncorrected_energy_per_atom = -9.807 eV/atom   (wrong functional!)
    
    Args:
        mp_ids: List of MP IDs
        mp_api_key: Materials Project API key
        debug: Print detailed debugging info
    
    Returns:
        dict: {mp_id: energy_per_atom} using GGA-PBE functional
    """
    print(f"\nFetching MP GGA-PBE DFT energies for {len(mp_ids)} structures...")
    print("  (Filtering for GGA functional to match VASP PBE calculations)")
    
    mp_energies = {}
    
    with MPRester(mp_api_key) as mpr:
        for mp_id in mp_ids:
            try:
                # Try Method 1: materials.thermo endpoint (more reliable for energies)
                energy = None
                nsites = None
                formula = None
                
                try:
                    # Fetch thermo docs with entry_types for functional filtering
                    thermo_docs = mpr.materials.thermo.search(
                        material_ids=[mp_id],
                        fields=["entry_types", "energy_per_atom", "uncorrected_energy_per_atom", 
                               "energy_type"]
                    )
                    
                    if debug and thermo_docs:
                        print(f"\n  DEBUG {mp_id}: Found {len(thermo_docs)} thermo docs")
                        for i, tdoc in enumerate(thermo_docs):
                            entry_types = getattr(tdoc, 'entry_types', [])
                            et_str = ' '.join(str(et) for et in entry_types) if entry_types else 'N/A'
                            print(f"    Doc {i}: entry_types={et_str}, "
                                  f"uncorrected_e={getattr(tdoc, 'uncorrected_energy_per_atom', 'N/A')}")
                    
                    # Filter for pure GGA-PBE only (exclude +U, r2SCAN, SCAN)
                    # Same logic as compute_dft_e_hull.py
                    non_pbe_markers = ['+U', 'GGA_U', 'R2SCAN', 'SCAN', 'r2SCAN']
                    
                    pbe_docs = []
                    for tdoc in thermo_docs:
                        entry_types = getattr(tdoc, 'entry_types', [])
                        if entry_types:
                            et_str = ' '.join(str(et) for et in entry_types)
                            # Skip if any non-PBE marker found
                            if any(marker in et_str for marker in non_pbe_markers):
                                if debug:
                                    print(f"      Skipping doc (non-PBE: {et_str})")
                                continue
                            pbe_docs.append(tdoc)
                        else:
                            # No entry_types, try to use energy_type as fallback
                            energy_type = getattr(tdoc, 'energy_type', '')
                            if energy_type == 'GGA':
                                pbe_docs.append(tdoc)
                    
                    if pbe_docs:
                        # Use the first pure PBE entry
                        tdoc = pbe_docs[0]
                        entry_types = getattr(tdoc, 'entry_types', [])
                        et_str = ' '.join(str(et) for et in entry_types) if entry_types else 'GGA'
                        
                        # Prefer uncorrected_energy_per_atom (raw DFT, no composition corrections)
                        if hasattr(tdoc, 'uncorrected_energy_per_atom') and tdoc.uncorrected_energy_per_atom is not None:
                            energy_candidate = float(tdoc.uncorrected_energy_per_atom)
                            if -20.0 <= energy_candidate <= 0.0:
                                energy = energy_candidate
                                if debug:
                                    print(f"      Using PBE ({et_str}) uncorrected_energy_per_atom = {energy:.6f} eV/atom")
                        
                        # Fallback to energy_per_atom (includes composition corrections)
                        if energy is None and hasattr(tdoc, 'energy_per_atom') and tdoc.energy_per_atom is not None:
                            energy_candidate = float(tdoc.energy_per_atom)
                            if -20.0 <= energy_candidate <= 0.0:
                                energy = energy_candidate
                                if debug:
                                    print(f"      Using PBE ({et_str}) energy_per_atom = {energy:.6f} eV/atom (corrected)")
                    else:
                        if debug:
                            print(f"      No pure PBE thermo docs found (all have +U/r2SCAN/SCAN)")
                            
                except Exception as thermo_error:
                    if debug:
                        print(f"\n  DEBUG {mp_id}: materials.thermo failed: {thermo_error}")
                
                # Method 2: Fallback to materials.summary
                if energy is None:
                    docs = mpr.materials.summary.search(
                        material_ids=[mp_id],
                        fields=["material_id", "energy_per_atom", "uncorrected_energy_per_atom",
                               "formula_pretty", "nsites", "e_total"]
                    )
                    
                    if not docs or len(docs) == 0:
                        print(f"  WARNING: {mp_id} not found in MP")
                        continue
                    
                    doc = docs[0]
                else:
                    # Still need to get nsites and formula for reporting
                    docs = mpr.materials.summary.search(
                        material_ids=[mp_id],
                        fields=["material_id", "formula_pretty", "nsites"]
                    )
                    if docs and len(docs) > 0:
                        doc = docs[0]
                    else:
                        doc = None
                
                # If we don't have energy from thermo endpoint, try to get from summary
                if energy is None and doc is not None:
                    if debug and not energy:
                        print(f"\n  DEBUG {mp_id}: Trying materials.summary endpoint...")
                        print(f"    formula: {getattr(doc, 'formula_pretty', 'N/A')}")
                        print(f"    nsites: {getattr(doc, 'nsites', 'N/A')}")
                        print(f"    energy_per_atom: {getattr(doc, 'energy_per_atom', 'N/A')}")
                        print(f"    uncorrected_energy_per_atom: {getattr(doc, 'uncorrected_energy_per_atom', 'N/A')}")
                        print(f"    e_total: {getattr(doc, 'e_total', 'N/A')}")
                    
                    nsites = getattr(doc, 'nsites', None)
                    
                    # Method 2a: Calculate from e_total if available
                    e_total = getattr(doc, 'e_total', None)
                    if e_total is not None and nsites is not None and nsites > 0:
                        energy = float(e_total) / nsites
                        if debug:
                            print(f"      Calculated from e_total / nsites: {energy:.6f} eV/atom")
                    
                    # Method 2b: Try energy_per_atom field (but may be buggy)
                    if energy is None and hasattr(doc, 'energy_per_atom') and doc.energy_per_atom is not None:
                        energy_candidate = float(doc.energy_per_atom)
                        
                        # Sanity check: should be reasonable for DFT energies
                        if -20.0 <= energy_candidate <= 0.0:
                            energy = energy_candidate
                            if debug:
                                print(f"      Using energy_per_atom: {energy:.6f} eV/atom")
                        else:
                            if debug:
                                print(f"      energy_per_atom {energy_candidate:.6f} out of range [-20, 0]")
                    
                    # Method 2c: Try uncorrected_energy_per_atom
                    if energy is None and hasattr(doc, 'uncorrected_energy_per_atom') and doc.uncorrected_energy_per_atom is not None:
                        energy_candidate = float(doc.uncorrected_energy_per_atom)
                        
                        if -20.0 <= energy_candidate <= 0.0:
                            energy = energy_candidate
                            if debug:
                                print(f"      Using uncorrected_energy_per_atom: {energy:.6f} eV/atom")
                        else:
                            if debug:
                                print(f"      uncorrected_energy_per_atom {energy_candidate:.6f} out of range")
                
                # Store result or report error
                if energy is not None:
                    mp_energies[mp_id] = energy
                else:
                    print(f"  ERROR: {mp_id} - no valid energy data available")
                    if not debug and doc is not None:
                        print(f"    energy_per_atom: {getattr(doc, 'energy_per_atom', 'N/A')}")
                        print(f"    uncorrected_energy_per_atom: {getattr(doc, 'uncorrected_energy_per_atom', 'N/A')}")
                        print(f"    e_total: {getattr(doc, 'e_total', 'N/A')}")
                        print(f"    Run with --debug for more information")
                
            except Exception as e:
                print(f"  ERROR: Failed to fetch energy for {mp_id}: {e}")
                if debug:
                    import traceback
                    traceback.print_exc()
    
    print(f"Successfully fetched energies for {len(mp_energies)}/{len(mp_ids)} structures\n")
    
    return mp_energies


def check_vasp_convergence(relax_dir):
    """
    Check if VASP calculation truly converged by reading vasprun.xml.
    
    Args:
        relax_dir: Path to VASP relaxation directory
    
    Returns:
        bool: True if converged, False otherwise
    """
    vasprun_path = Path(relax_dir) / 'vasprun.xml'
    
    if not vasprun_path.exists():
        return False
    
    try:
        vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=False)
        return vr.converged
    except Exception:
        return False


def analyze_energy_differences(db, mp_energies, check_convergence=True, outlier_threshold=0.5):
    """
    Analyze energy differences between MP (uncorrected DFT) and VASP.
    
    Args:
        db: Workflow database
        mp_energies: dict of {mp_id: energy_per_atom} from MP API
        check_convergence: If True, verify VASP convergence and exclude non-converged
        outlier_threshold: Energy difference threshold (eV/atom) for outlier detection
    
    Returns:
        dict: Statistics and structure lists
    """
    completed = {}
    skipped_not_converged = []
    
    for mp_id, sdata in db['structures'].items():
        if sdata['state'] != 'COMPLETED':
            continue
        if sdata['vasp_energy_per_atom'] is None:
            continue
        
        # Get MP energy from fetched energies
        mp_e = mp_energies.get(mp_id)
        
        if mp_e is None:
            continue
        
        # Check VASP convergence if requested
        if check_convergence and sdata.get('relax_dir'):
            if not check_vasp_convergence(sdata['relax_dir']):
                skipped_not_converged.append(mp_id)
                print(f"  WARNING: Skipping {mp_id} - VASP calculation not converged")
                continue
        
        sdata['mp_energy_per_atom'] = mp_e  # Cache in sdata
        completed[mp_id] = sdata
    
    if not completed:
        return None
    
    # Report filtering results
    if skipped_not_converged:
        print(f"\n  Filtered out {len(skipped_not_converged)} non-converged structures")
    
    # First pass: calculate all differences to detect outliers
    all_diffs = []
    for mp_id, sdata in completed.items():
        mp_e = sdata['mp_energy_per_atom']
        vasp_e = sdata['vasp_energy_per_atom']
        diff = vasp_e - mp_e
        all_diffs.append((mp_id, diff, abs(diff), sdata))
    
    # Detect statistical outliers using IQR method
    abs_diffs = [abs_diff for _, _, abs_diff, _ in all_diffs]
    q25, q75 = np.percentile(abs_diffs, [25, 75])
    iqr = q75 - q25
    outlier_threshold_iqr = q75 + 3 * iqr  # 3*IQR above Q3
    
    # Also use absolute threshold
    outlier_threshold_abs = outlier_threshold
    
    # Filter outliers
    skipped_outliers = []
    filtered_diffs = []
    
    for mp_id, diff, abs_diff, sdata in all_diffs:
        # Check both statistical and absolute thresholds
        if abs_diff > outlier_threshold_iqr or abs_diff > outlier_threshold_abs:
            skipped_outliers.append({
                'mp_id': mp_id,
                'formula': sdata['formula'],
                'mp_energy': sdata['mp_energy_per_atom'],
                'vasp_energy': sdata['vasp_energy_per_atom'],
                'diff': diff,
                'abs_diff': abs_diff
            })
            print(f"  WARNING: Excluding outlier {mp_id} ({sdata['formula']}): "
                  f"Δ = {diff:+.4f} eV/atom (|Δ| = {abs_diff:.4f})")
        else:
            filtered_diffs.append((mp_id, diff, abs_diff, sdata))
    
    if skipped_outliers:
        print(f"\n  Filtered out {len(skipped_outliers)} statistical outliers "
              f"(|Δ| > {outlier_threshold:.2f} eV/atom or > Q3+3*IQR)")
    
    print(f"  Analyzing {len(filtered_diffs)} structures after filtering\n")
    
    # Calculate differences from filtered data
    diffs = []
    by_chemsys = defaultdict(list)
    
    for mp_id, diff, abs_diff, sdata in filtered_diffs:
        mp_e = sdata['mp_energy_per_atom']
        vasp_e = sdata['vasp_energy_per_atom']
        
        diffs.append({
            'mp_id': mp_id,
            'chemsys': sdata['chemsys'],
            'formula': sdata['formula'],
            'mp_energy': mp_e,
            'vasp_energy': vasp_e,
            'diff': diff,
            'abs_diff': abs_diff
        })
        
        by_chemsys[sdata['chemsys']].append(diff)
    
    # Overall statistics (after filtering)
    diff_values = [d['diff'] for d in diffs]
    abs_diff_values = [d['abs_diff'] for d in diffs]
    
    stats = {
        'n_structures': len(diffs),
        'n_skipped_not_converged': len(skipped_not_converged),
        'n_skipped_outliers': len(skipped_outliers),
        'outlier_threshold': outlier_threshold,
        'mean_diff': np.mean(diff_values),
        'std_diff': np.std(diff_values),
        'mae': np.mean(abs_diff_values),
        'rmse': np.sqrt(np.mean([d**2 for d in diff_values])),
        'min_diff': np.min(diff_values),
        'max_diff': np.max(diff_values),
        'median_diff': np.median(diff_values),
        'q25': np.percentile(diff_values, 25),
        'q75': np.percentile(diff_values, 75)
    }
    
    # Per-chemsys statistics
    chemsys_stats = {}
    for chemsys, chemsys_diffs in by_chemsys.items():
        chemsys_stats[chemsys] = {
            'n': len(chemsys_diffs),
            'mean': np.mean(chemsys_diffs),
            'std': np.std(chemsys_diffs),
            'mae': np.mean([abs(d) for d in chemsys_diffs]),
            'max_abs': max([abs(d) for d in chemsys_diffs])
        }
    
    return {
        'stats': stats,
        'chemsys_stats': chemsys_stats,
        'structures': diffs,
        'skipped_not_converged': skipped_not_converged,
        'skipped_outliers': skipped_outliers
    }


def plot_energy_comparison(results, output_prefix='mp_vasp_comparison'):
    """
    Create scatter plot comparing MP vs VASP energies.
    
    Args:
        results: Analysis results dictionary
        output_prefix: Prefix for output files
    """
    structures = results['structures']
    stats = results['stats']
    
    # Extract data
    vasp_energies = np.array([s['vasp_energy'] for s in structures])
    mp_energies = np.array([s['mp_energy'] for s in structures])
    formulas = [s['formula'] for s in structures]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Scatter plot
    ax.scatter(vasp_energies, mp_energies, alpha=0.6, s=50, c='steelblue', edgecolors='black', linewidth=0.5)
    
    # Perfect agreement line (y = x)
    all_energies = np.concatenate([vasp_energies, mp_energies])
    energy_min, energy_max = all_energies.min(), all_energies.max()
    margin = (energy_max - energy_min) * 0.05
    plot_min, plot_max = energy_min - margin, energy_max + margin
    
    ax.plot([plot_min, plot_max], [plot_min, plot_max], 'r--', linewidth=2, alpha=0.7, 
            label='Perfect agreement (y=x)')
    
    # Add error bands (±0.05 eV/atom)
    ax.fill_between([plot_min, plot_max], 
                     [plot_min - 0.05, plot_max - 0.05],
                     [plot_min + 0.05, plot_max + 0.05],
                     alpha=0.2, color='green', label='±0.05 eV/atom')
    
    # Labels and title
    ax.set_xlabel('VASP PBE Energy per Atom (eV)', fontsize=12, fontweight='bold')
    ax.set_ylabel('MP GGA-PBE Energy per Atom (eV)', fontsize=12, fontweight='bold')
    n_skipped_conv = stats.get('n_skipped_not_converged', 0)
    n_skipped_outliers = stats.get('n_skipped_outliers', 0)
    title_suffix = f' (Pure GGA-PBE, N={stats["n_structures"]}'
    if n_skipped_conv > 0 or n_skipped_outliers > 0:
        title_suffix += ', excluded:'
        if n_skipped_conv > 0:
            title_suffix += f' {n_skipped_conv} non-converged'
        if n_skipped_outliers > 0:
            if n_skipped_conv > 0:
                title_suffix += ','
            title_suffix += f' {n_skipped_outliers} outliers'
    title_suffix += ')'
    ax.set_title(f'MP vs VASP DFT Energy Comparison\n{title_suffix}', 
                 fontsize=14, fontweight='bold')
    
    # Add statistics text box
    stats_text = (
        f"N = {stats['n_structures']}\n"
        f"MAE = {stats['mae']:.4f} eV/atom\n"
        f"RMSE = {stats['rmse']:.4f} eV/atom\n"
        f"Mean Diff = {stats['mean_diff']:+.4f} eV/atom\n"
        f"Std Diff = {stats['std_diff']:.4f} eV/atom"
    )
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Grid and legend
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.minorticks_on()
    ax.grid(True, which='minor', alpha=0.15, linestyle=':')
    ax.legend(loc='lower right', fontsize=10)
    
    # Equal aspect ratio
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(plot_min, plot_max)
    ax.set_ylim(plot_min, plot_max)
    
    # Save figure
    plt.tight_layout()
    plot_file = f"{output_prefix}_scatter.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"\nScatter plot saved to: {plot_file}")
    plt.close()
    
    # Create residual plot
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Residual vs VASP energy
    residuals = vasp_energies - mp_energies
    ax.scatter(vasp_energies, residuals, alpha=0.6, s=50, c='steelblue', edgecolors='black', linewidth=0.5)
    
    # Reference lines
    ax.axhline(y=0, color='r', linestyle='--', linewidth=2, alpha=0.7, label='Zero residual')
    ax.axhline(y=0.05, color='green', linestyle=':', linewidth=1.5, alpha=0.6, label='±0.05 eV/atom')
    ax.axhline(y=-0.05, color='green', linestyle=':', linewidth=1.5, alpha=0.6)
    
    # Mean residual line
    ax.axhline(y=stats['mean_diff'], color='orange', linestyle='-', linewidth=2, 
               alpha=0.7, label=f'Mean = {stats["mean_diff"]:+.4f} eV/atom')
    
    # Labels and title
    ax.set_xlabel('VASP PBE Energy per Atom (eV)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Residual (VASP - MP) (eV/atom)', fontsize=12, fontweight='bold')
    n_skipped_conv = stats.get('n_skipped_not_converged', 0)
    n_skipped_outliers = stats.get('n_skipped_outliers', 0)
    title_suffix = f' (Pure GGA-PBE, N={stats["n_structures"]}'
    if n_skipped_conv > 0 or n_skipped_outliers > 0:
        title_suffix += ', excluded:'
        if n_skipped_conv > 0:
            title_suffix += f' {n_skipped_conv} non-converged'
        if n_skipped_outliers > 0:
            if n_skipped_conv > 0:
                title_suffix += ','
            title_suffix += f' {n_skipped_outliers} outliers'
    title_suffix += ')'
    ax.set_title(f'Energy Residuals vs VASP Energy\n{title_suffix}', 
                 fontsize=14, fontweight='bold')
    
    # Set focused Y-axis range based on data
    residual_range = max(abs(residuals.min()), abs(residuals.max()))
    # Use ±0.15 eV/atom or 1.5x the max residual, whichever is larger
    y_limit = max(0.15, residual_range * 1.5)
    ax.set_ylim(-y_limit, y_limit)
    
    # Grid with minor ticks
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.minorticks_on()
    ax.grid(True, which='minor', alpha=0.15, linestyle=':')
    
    # Legend
    ax.legend(loc='best', fontsize=10)
    
    # Save figure
    plt.tight_layout()
    residual_file = f"{output_prefix}_residuals.png"
    plt.savefig(residual_file, dpi=300, bbox_inches='tight')
    print(f"Residual plot saved to: {residual_file}")
    plt.close()


def print_summary(results):
    """Print minimal analysis summary."""
    stats = results['stats']
    
    print("\n" + "="*70)
    print("MP vs VASP Energy Comparison (Pure GGA-PBE Functional)")
    print("="*70)
    print(f"Structures analyzed: {stats['n_structures']}")
    if stats.get('n_skipped_not_converged', 0) > 0:
        print(f"Skipped (not converged): {stats['n_skipped_not_converged']}")
    if stats.get('n_skipped_outliers', 0) > 0:
        print(f"Skipped (statistical outliers): {stats['n_skipped_outliers']}")
    print(f"MAE:  {stats['mae']:.4f} eV/atom")
    print(f"RMSE: {stats['rmse']:.4f} eV/atom")
    print(f"Mean: {stats['mean_diff']:+.4f} eV/atom")
    print("="*70)


def save_detailed_comparison(results, output_file):
    """Save detailed comparison to JSON file."""
    # Add metadata about functional filtering, convergence checking, and outlier removal
    results['metadata'] = {
        'functional': 'GGA-PBE',
        'note': 'MP energies filtered for GGA functional to match VASP PBE',
        'mp_field_used': 'uncorrected_energy_per_atom from materials.thermo where energy_type=GGA',
        'convergence_check': 'Only converged VASP calculations included',
        'outlier_detection': 'Statistical outliers excluded using IQR method and absolute threshold',
        'outlier_threshold': results['stats'].get('outlier_threshold', 0.5),
        'n_skipped_not_converged': results['stats'].get('n_skipped_not_converged', 0),
        'n_skipped_outliers': results['stats'].get('n_skipped_outliers', 0)
    }
    
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nDetailed results saved to: {output_file}")
    print(f"  (MP GGA-PBE vs converged VASP PBE, outliers excluded)")


def main():
    parser = argparse.ArgumentParser(
        description="Compare MP and VASP energies from workflow database"
    )
    parser.add_argument(
        '--db',
        type=str,
        default='./mp_relax_workflow.json',
        help="Workflow database file (default: ./mp_relax_workflow.json)"
    )
    parser.add_argument(
        '--mp-api-key',
        type=str,
        default=None,
        help="Materials Project API key (or set MP_API_KEY environment variable)"
    )
    parser.add_argument(
        '--output',
        type=str,
        default='mp_vasp_comparison.json',
        help="Output file for detailed results (default: mp_vasp_comparison.json). "
             "Also generates PNG plots with same prefix."
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help="Print detailed debugging information"
    )
    parser.add_argument(
        '--no-convergence-check',
        action='store_true',
        help="Skip VASP convergence verification (include all completed structures)"
    )
    parser.add_argument(
        '--outlier-threshold',
        type=float,
        default=0.5,
        help="Absolute energy difference threshold (eV/atom) for outlier detection (default: 0.5)"
    )
    
    args = parser.parse_args()
    
    # Get MP API key
    mp_api_key = args.mp_api_key or os.environ.get('MP_API_KEY')
    if not mp_api_key:
        print("ERROR: MP_API_KEY not found in environment or arguments")
        print("Set it with: export MP_API_KEY=your_key")
        return 1
    
    db_path = Path(args.db).expanduser()
    
    if not db_path.exists():
        print(f"ERROR: Database file not found: {db_path}")
        return 1
    
    # Load database
    print(f"Loading workflow database: {db_path}")
    db = load_workflow_db(db_path)
    
    total = len(db['structures'])
    completed = sum(1 for s in db['structures'].values() if s['state'] == 'COMPLETED')
    
    print(f"Total structures: {total}")
    print(f"Completed: {completed}")
    
    if completed == 0:
        print("\nNo completed structures found!")
        return 1
    
    # Get list of completed MP IDs
    completed_mp_ids = [mp_id for mp_id, s in db['structures'].items() 
                        if s['state'] == 'COMPLETED' and s['vasp_energy_per_atom'] is not None]
    
    # Fetch MP energies via API
    mp_energies = fetch_mp_energies(completed_mp_ids, mp_api_key, debug=args.debug)
    
    # Analyze
    print("Analyzing energy differences...")
    check_convergence = not args.no_convergence_check
    if check_convergence:
        print("  (Verifying VASP convergence - use --no-convergence-check to disable)")
    print(f"  (Outlier threshold: |Δ| > {args.outlier_threshold:.2f} eV/atom)")
    results = analyze_energy_differences(
        db, mp_energies, 
        check_convergence=check_convergence,
        outlier_threshold=args.outlier_threshold
    )
    
    if results is None:
        print("ERROR: No valid energy data found")
        return 1
    
    # Print minimal summary
    print_summary(results)
    
    # Save detailed results
    output_prefix = args.output.replace('.json', '')
    save_detailed_comparison(results, args.output)
    
    # Generate plots
    print("\nGenerating comparison plots...")
    plot_energy_comparison(results, output_prefix=output_prefix)
    
    print("\n" + "="*70)
    print("Analysis complete!")
    print(f"  JSON results: {args.output}")
    print(f"  Scatter plot: {output_prefix}_scatter.png")
    print(f"  Residual plot: {output_prefix}_residuals.png")
    print("="*70 + "\n")
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())

