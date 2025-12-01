#!/usr/bin/env python3
"""
Compare energies from different calculation methods:

1. MP On-Hull Phases:
   - MatterSim relaxed (mp_mattersim.json)
   - MP raw DFT energies (mp_vaspdft.json)

2. Generated Structures:
   - MatterSim relaxed (prescreening_stability.json)
   - VASP-DFT relaxed (workflow.json + Relax/vasprun.xml)

Purpose: Diagnose non-linear behavior in energy_above_hull calculations
"""

import os
import sys
import json
import argparse
import warnings
import numpy as np
from pathlib import Path
from collections import defaultdict

from pymatgen.core import Composition
from pymatgen.io.vasp.outputs import Vasprun, UnconvergedVASPWarning
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use('Agg')

# Suppress unconverged VASP warnings (already handled by checking vr.converged)
warnings.filterwarnings('ignore', category=UnconvergedVASPWarning)


def load_json_file(filepath):
    """Load JSON file safely."""
    filepath = Path(filepath)
    if not filepath.exists():
        print(f"ERROR: File not found: {filepath}")
        return None
    
    with open(filepath, 'r') as f:
        return json.load(f)


def extract_mp_id(entry_id):
    """
    Extract base MP ID from various formats.
    
    Examples:
        'mp_mattersim_mp-12345' -> 'mp-12345'
        'mp-12345-GGA' -> 'mp-12345'
        'mp-12345' -> 'mp-12345'
    """
    if 'mp_mattersim_' in entry_id:
        mp_id = entry_id.split('mp_mattersim_')[1]
    else:
        mp_id = entry_id
    
    # Remove functional suffixes
    mp_id = mp_id.split('-GGA')[0].split('-r2SCAN')[0].split('_fallback')[0]
    
    return mp_id


def compare_mp_phases(mattersim_cache, dft_cache):
    """
    Compare MP on-hull phases: MatterSim vs MP raw DFT energies.
    
    Returns:
        dict: Comparison results with statistics and matched entries
    """
    print("\n" + "="*70)
    print("Comparing MP On-Hull Phases: MatterSim vs MP Raw DFT")
    print("="*70)
    
    # Parse MatterSim cache
    ms_by_mpid = {}
    for item in mattersim_cache:
        mp_id = extract_mp_id(item['entry_id'])
        ms_by_mpid[mp_id] = {
            'energy': item['energy'],
            'composition': Composition(item['composition']),
            'chemsys': item['chemsys'],
            'entry_id': item['entry_id']
        }
    
    # Parse DFT cache
    dft_by_mpid = {}
    for item in dft_cache:
        mp_id = extract_mp_id(item['entry_id'])
        dft_by_mpid[mp_id] = {
            'energy': item['energy'],
            'composition': Composition(item['composition']),
            'chemsys': item['chemsys'],
            'entry_id': item['entry_id']
        }
    
    print(f"MatterSim cache: {len(ms_by_mpid)} unique MP IDs")
    print(f"DFT cache: {len(dft_by_mpid)} unique MP IDs")
    
    # Match entries
    common_ids = set(ms_by_mpid.keys()) & set(dft_by_mpid.keys())
    ms_only = set(ms_by_mpid.keys()) - set(dft_by_mpid.keys())
    dft_only = set(dft_by_mpid.keys()) - set(ms_by_mpid.keys())
    
    print(f"\nCache Overlap:")
    print(f"  Common to both: {len(common_ids)}")
    print(f"  MatterSim only: {len(ms_only)}")
    print(f"  DFT only: {len(dft_only)}")
    
    if len(ms_only) > 0:
        coverage = len(common_ids) / len(ms_by_mpid) * 100
        print(f"  Coverage: {coverage:.1f}% of MatterSim cache has DFT energies")
        if coverage < 50:
            print(f"    Note: Low coverage - DFT hull calculations may still be running")
    print()
    
    if not common_ids:
        print("ERROR: No common MP IDs found - cannot compare!\n")
        print("Possible reasons:")
        print("  - DFT hull calculations not started yet (mp_vaspdft.json empty)")
        print("  - Different chemical systems in the two caches")
        print("\nSuggestion: Wait for compute_dft_e_hull.py to process more structures\n")
        return None
    
    # Compare energies
    matched = []
    
    for mp_id in sorted(common_ids):
        ms_data = ms_by_mpid[mp_id]
        dft_data = dft_by_mpid[mp_id]
        
        # Get composition (should be same, but check)
        ms_comp = ms_data['composition']
        dft_comp = dft_data['composition']
        
        # Use reduced formula
        formula = ms_comp.reduced_formula
        
        # Total energies
        ms_energy = ms_data['energy']
        dft_energy = dft_data['energy']
        
        # Compute per-atom energies
        ms_natoms = ms_comp.num_atoms
        dft_natoms = dft_comp.num_atoms
        
        ms_energy_per_atom = ms_energy / ms_natoms
        dft_energy_per_atom = dft_energy / dft_natoms
        
        # Check composition consistency
        comp_match = ms_comp.reduced_formula == dft_comp.reduced_formula
        
        matched.append({
            'mp_id': mp_id,
            'formula': formula,
            'chemsys': ms_data['chemsys'],
            'mattersim_total_energy': ms_energy,
            'dft_total_energy': dft_energy,
            'mattersim_energy_per_atom': ms_energy_per_atom,
            'dft_energy_per_atom': dft_energy_per_atom,
            'natoms': ms_natoms,
            'composition_match': comp_match,
            'energy_diff_total': ms_energy - dft_energy,
            'energy_diff_per_atom': ms_energy_per_atom - dft_energy_per_atom
        })
    
    # Calculate statistics
    total_diffs = np.array([m['energy_diff_total'] for m in matched])
    per_atom_diffs = np.array([m['energy_diff_per_atom'] for m in matched])
    
    ms_total = np.array([m['mattersim_total_energy'] for m in matched])
    dft_total = np.array([m['dft_total_energy'] for m in matched])
    ms_per_atom = np.array([m['mattersim_energy_per_atom'] for m in matched])
    dft_per_atom = np.array([m['dft_energy_per_atom'] for m in matched])
    
    stats = {
        'n_structures': len(matched),
        'total_energy': {
            'mae': np.mean(np.abs(total_diffs)),
            'rmse': np.sqrt(np.mean(total_diffs**2)),
            'mean_diff': np.mean(total_diffs),
            'std_diff': np.std(total_diffs),
            'correlation': np.corrcoef(ms_total, dft_total)[0, 1]
        },
        'per_atom_energy': {
            'mae': np.mean(np.abs(per_atom_diffs)),
            'rmse': np.sqrt(np.mean(per_atom_diffs**2)),
            'mean_diff': np.mean(per_atom_diffs),
            'std_diff': np.std(per_atom_diffs),
            'correlation': np.corrcoef(ms_per_atom, dft_per_atom)[0, 1]
        }
    }
    
    # Calculate coverage percentage
    coverage = len(common_ids) / len(ms_by_mpid) * 100 if len(ms_by_mpid) > 0 else 0
    
    return {
        'matched': matched,
        'statistics': stats,
        'coverage': coverage,
        'n_mattersim': len(ms_by_mpid),
        'n_dft': len(dft_by_mpid),
        'n_common': len(common_ids)
    }


def compare_generated_structures(prescreen_json, workflow_json, vasp_jobs_dir):
    """
    Compare generated structures: MatterSim (prescreening) vs VASP-DFT (workflow).
    
    Returns:
        dict: Comparison results with statistics and matched structures
    """
    print("\n" + "="*70)
    print("Comparing Generated Structures: MatterSim vs VASP-DFT")
    print("="*70)
    
    # Load prescreening results
    prescreen_data = load_json_file(prescreen_json)
    if prescreen_data is None:
        return None
    
    # Load workflow database
    workflow_data = load_json_file(workflow_json)
    if workflow_data is None:
        return None
    
    # Extract MatterSim energies for passed structures
    ms_energies = {}
    for result in prescreen_data.get('results', []):
        if result.get('passed_prescreening', False) and result.get('mattersim_energy_per_atom') is not None:
            struct_id = result['structure_id']
            ms_energies[struct_id] = {
                'energy_per_atom': result['mattersim_energy_per_atom'],
                'composition': result['composition'],
                'chemsys': result['chemsys']
            }
    
    print(f"MatterSim prescreening: {len(ms_energies)} structures passed threshold")
    
    # Extract VASP energies from workflow
    vasp_energies = {}
    vasp_jobs_dir = Path(vasp_jobs_dir)
    vasp_not_converged = []
    vasp_no_vasprun = []
    
    for struct_id, sdata in workflow_data['structures'].items():
        # Check if VASP relaxation completed
        if sdata['state'] not in ['RELAX_DONE', 'SC_RUNNING', 'SC_DONE', 'SC_FAILED',
                                   'PARCHG_RUNNING', 'PARCHG_DONE', 'PARCHG_FAILED',
                                   'ELF_RUNNING', 'ELF_DONE', 'ELF_FAILED']:
            continue
        
        relax_dir = Path(sdata['relax_dir'])
        vasprun_path = relax_dir / 'vasprun.xml'
        
        if not vasprun_path.exists():
            # Track structures in workflow but no vasprun.xml yet
            if struct_id in ms_energies:
                vasp_no_vasprun.append(struct_id)
            continue
        
        try:
            vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=False)
            
            # Check convergence
            if not vr.converged:
                if struct_id in ms_energies:
                    vasp_not_converged.append(struct_id)
                continue
            
            final_energy = vr.final_energy
            structure = vr.final_structure
            energy_per_atom = final_energy / len(structure)
            
            vasp_energies[struct_id] = {
                'energy_per_atom': energy_per_atom,
                'composition': sdata['composition'],
                'chemsys': sdata.get('chemsys', ''),
                'converged': True
            }
            
        except Exception as e:
            # Track parsing errors
            if struct_id in ms_energies:
                vasp_no_vasprun.append(struct_id)
            continue
    
    print(f"VASP-DFT: {len(vasp_energies)} converged relaxations")
    
    # Match structures
    common_ids = set(ms_energies.keys()) & set(vasp_energies.keys())
    ms_only = set(ms_energies.keys()) - set(vasp_energies.keys())
    
    print(f"\nStructure Overlap:")
    print(f"  Passed prescreening: {len(ms_energies)}")
    print(f"  Completed VASP-DFT: {len(vasp_energies)}")
    print(f"  Common (for comparison): {len(common_ids)}")
    
    if len(ms_only) > 0:
        print(f"  Prescreened but no VASP-DFT: {len(ms_only)}")
        if len(vasp_not_converged) > 0:
            print(f"    - Not converged: {len(vasp_not_converged)}")
        if len(vasp_no_vasprun) > 0:
            print(f"    - No vasprun.xml / running: {len(vasp_no_vasprun)}")
        
        coverage = len(common_ids) / len(ms_energies) * 100
        print(f"  Coverage: {coverage:.1f}% of prescreened structures have VASP-DFT energies")
        
        if coverage < 30:
            print(f"    Note: Low coverage - VASP calculations may still be in progress")
        elif coverage < 70:
            print(f"  ℹ Note: Partial coverage - workflow still running")
    print()
    
    if not common_ids:
        print("ERROR: No common structure IDs found - cannot compare!\n")
        print("Possible reasons:")
        print("  - VASP workflow not started yet")
        print("  - All VASP calculations still running or failed to converge")
        print("  - Structures failed VASP relaxation")
        print("\nSuggestion: Wait for more VASP relaxations to complete\n")
        return None
    
    # Compare energies
    matched = []
    
    for struct_id in sorted(common_ids):
        ms_data = ms_energies[struct_id]
        vasp_data = vasp_energies[struct_id]
        
        ms_e = ms_data['energy_per_atom']
        vasp_e = vasp_data['energy_per_atom']
        
        matched.append({
            'structure_id': struct_id,
            'composition': ms_data['composition'],
            'chemsys': ms_data['chemsys'],
            'mattersim_energy_per_atom': ms_e,
            'vasp_energy_per_atom': vasp_e,
            'energy_diff': ms_e - vasp_e
        })
    
    # Calculate statistics
    diffs = np.array([m['energy_diff'] for m in matched])
    ms_vals = np.array([m['mattersim_energy_per_atom'] for m in matched])
    vasp_vals = np.array([m['vasp_energy_per_atom'] for m in matched])
    
    stats = {
        'n_structures': len(matched),
        'mae': np.mean(np.abs(diffs)),
        'rmse': np.sqrt(np.mean(diffs**2)),
        'mean_diff': np.mean(diffs),
        'std_diff': np.std(diffs),
        'correlation': np.corrcoef(ms_vals, vasp_vals)[0, 1],
        'min_diff': np.min(diffs),
        'max_diff': np.max(diffs)
    }
    
    # Calculate coverage percentage
    coverage = len(common_ids) / len(ms_energies) * 100 if len(ms_energies) > 0 else 0
    
    return {
        'matched': matched,
        'statistics': stats,
        'coverage': coverage,
        'n_prescreened': len(ms_energies),
        'n_vasp': len(vasp_energies),
        'n_common': len(common_ids)
    }


def plot_mp_phase_comparison(mp_results, output_prefix='mp_phases_comparison'):
    """Plot comparison of MP phases: MatterSim vs DFT."""
    if mp_results is None:
        return
    
    print("  Generating MP phases plots...")
    matched = mp_results['matched']
    stats = mp_results['statistics']
    
    # Extract data
    ms_per_atom = np.array([m['mattersim_energy_per_atom'] for m in matched])
    dft_per_atom = np.array([m['dft_energy_per_atom'] for m in matched])
    
    # Scatter plot
    fig, ax = plt.subplots(figsize=(12, 10))
    
    ax.scatter(dft_per_atom, ms_per_atom, alpha=0.6, s=50, 
              c='steelblue', edgecolors='black', linewidth=0.5)
    
    # Perfect agreement line
    all_vals = np.concatenate([ms_per_atom, dft_per_atom])
    val_min, val_max = all_vals.min(), all_vals.max()
    margin = (val_max - val_min) * 0.05
    plot_min, plot_max = val_min - margin, val_max + margin
    
    ax.plot([plot_min, plot_max], [plot_min, plot_max], 'r--', 
            linewidth=2, alpha=0.7, label='Perfect agreement (y=x)')
    
    # Labels
    ax.set_xlabel('MP Raw DFT Energy per Atom (eV)', fontsize=12, fontweight='bold')
    ax.set_ylabel('MatterSim Energy per Atom (eV)', fontsize=12, fontweight='bold')
    ax.set_title(f'MP On-Hull Phases: MatterSim vs MP Raw DFT\n(N={stats["n_structures"]})',
                fontsize=14, fontweight='bold')
    
    # Statistics text
    per_atom_stats = stats['per_atom_energy']
    stats_text = (
        f"N = {stats['n_structures']}\n"
        f"R = {per_atom_stats['correlation']:.4f}\n"
        f"MAE = {per_atom_stats['mae']:.4f} eV/atom\n"
        f"RMSE = {per_atom_stats['rmse']:.4f} eV/atom\n"
        f"Mean Diff = {per_atom_stats['mean_diff']:+.4f} eV/atom"
    )
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Grid and legend
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='lower right')
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(plot_min, plot_max)
    ax.set_ylim(plot_min, plot_max)
    
    plt.tight_layout()
    scatter_file = f"{output_prefix}_scatter.png"
    plt.savefig(scatter_file, dpi=300, bbox_inches='tight')
    print(f"      Scatter: {Path(scatter_file).name}")
    plt.close()
    
    # Residual plot
    fig, ax = plt.subplots(figsize=(12, 10))
    
    residuals = ms_per_atom - dft_per_atom
    
    ax.scatter(dft_per_atom, residuals, alpha=0.6, s=50,
              c='steelblue', edgecolors='black', linewidth=0.5)
    
    ax.axhline(y=0, color='r', linestyle='--', linewidth=2, alpha=0.7,
              label='Zero residual')
    ax.axhline(y=per_atom_stats['mean_diff'], color='orange', linestyle='-',
              linewidth=2, alpha=0.7, 
              label=f'Mean = {per_atom_stats["mean_diff"]:+.4f} eV/atom')
    
    ax.set_xlabel('MP Raw DFT Energy per Atom (eV)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Residual (MatterSim - MP DFT) (eV/atom)', fontsize=12, fontweight='bold')
    ax.set_title(f'MP On-Hull Phases: Energy Residuals\n(N={stats["n_structures"]})',
                fontsize=14, fontweight='bold')
    
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='best')
    
    plt.tight_layout()
    residual_file = f"{output_prefix}_residuals.png"
    plt.savefig(residual_file, dpi=300, bbox_inches='tight')
    print(f"      Residual: {Path(residual_file).name}")
    plt.close()


def plot_generated_structures_comparison(gen_results, output_prefix='generated_structures_comparison'):
    """Plot comparison of generated structures: MatterSim vs VASP-DFT."""
    if gen_results is None:
        return
    
    print("  Generating generated structures plots...")
    matched = gen_results['matched']
    stats = gen_results['statistics']
    
    # Extract data
    ms_vals = np.array([m['mattersim_energy_per_atom'] for m in matched])
    vasp_vals = np.array([m['vasp_energy_per_atom'] for m in matched])
    
    # Scatter plot
    fig, ax = plt.subplots(figsize=(12, 10))
    
    ax.scatter(vasp_vals, ms_vals, alpha=0.6, s=50,
              c='forestgreen', edgecolors='black', linewidth=0.5)
    
    # Perfect agreement line
    all_vals = np.concatenate([ms_vals, vasp_vals])
    val_min, val_max = all_vals.min(), all_vals.max()
    margin = (val_max - val_min) * 0.05
    plot_min, plot_max = val_min - margin, val_max + margin
    
    ax.plot([plot_min, plot_max], [plot_min, plot_max], 'r--',
            linewidth=2, alpha=0.7, label='Perfect agreement (y=x)')
    
    # Labels
    ax.set_xlabel('VASP-DFT Energy per Atom (eV)', fontsize=12, fontweight='bold')
    ax.set_ylabel('MatterSim Energy per Atom (eV)', fontsize=12, fontweight='bold')
    ax.set_title(f'Generated Structures: MatterSim vs VASP-DFT\n(N={stats["n_structures"]})',
                fontsize=14, fontweight='bold')
    
    # Statistics text
    stats_text = (
        f"N = {stats['n_structures']}\n"
        f"R = {stats['correlation']:.4f}\n"
        f"MAE = {stats['mae']:.4f} eV/atom\n"
        f"RMSE = {stats['rmse']:.4f} eV/atom\n"
        f"Mean Diff = {stats['mean_diff']:+.4f} eV/atom"
    )
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Grid and legend
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='lower right')
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(plot_min, plot_max)
    ax.set_ylim(plot_min, plot_max)
    
    plt.tight_layout()
    scatter_file = f"{output_prefix}_scatter.png"
    plt.savefig(scatter_file, dpi=300, bbox_inches='tight')
    print(f"      Scatter: {Path(scatter_file).name}")
    plt.close()
    
    # Residual plot
    fig, ax = plt.subplots(figsize=(12, 10))
    
    residuals = ms_vals - vasp_vals
    
    ax.scatter(vasp_vals, residuals, alpha=0.6, s=50,
              c='forestgreen', edgecolors='black', linewidth=0.5)
    
    ax.axhline(y=0, color='r', linestyle='--', linewidth=2, alpha=0.7,
              label='Zero residual')
    ax.axhline(y=stats['mean_diff'], color='orange', linestyle='-',
              linewidth=2, alpha=0.7,
              label=f'Mean = {stats["mean_diff"]:+.4f} eV/atom')
    
    ax.set_xlabel('VASP-DFT Energy per Atom (eV)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Residual (MatterSim - VASP) (eV/atom)', fontsize=12, fontweight='bold')
    ax.set_title(f'Generated Structures: Energy Residuals\n(N={stats["n_structures"]})',
                fontsize=14, fontweight='bold')
    
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='best')
    
    plt.tight_layout()
    residual_file = f"{output_prefix}_residuals.png"
    plt.savefig(residual_file, dpi=300, bbox_inches='tight')
    print(f"      Residual: {Path(residual_file).name}")
    plt.close()


def print_summary(mp_results, gen_results):
    """Print comprehensive summary of both comparisons."""
    print("\n" + "="*70)
    print("COMPREHENSIVE ENERGY COMPARISON SUMMARY")
    print("="*70)
    
    # Check if results are partial (workflow still running)
    partial_results = False
    partial_msg = []
    
    if mp_results:
        print("\n1. MP On-Hull Phases (MatterSim vs MP Raw DFT):")
        print("-" * 70)
        stats = mp_results['statistics']
        pa_stats = stats['per_atom_energy']
        
        print(f"  Structures compared: {stats['n_structures']}")
        print(f"\n  Per-Atom Energy:")
        print(f"    Correlation: {pa_stats['correlation']:.4f}")
        print(f"    MAE:  {pa_stats['mae']:.4f} eV/atom")
        print(f"    RMSE: {pa_stats['rmse']:.4f} eV/atom")
        print(f"    Mean Diff (MS - DFT): {pa_stats['mean_diff']:+.4f} eV/atom")
        print(f"    Std Dev: {pa_stats['std_diff']:.4f} eV/atom")
        
        # Check for systematic bias
        if abs(pa_stats['mean_diff']) > 0.05:
            print(f"\n    WARNING: Systematic bias detected!")
            print(f"    MatterSim is {'higher' if pa_stats['mean_diff'] > 0 else 'lower'} "
                  f"than MP DFT by {abs(pa_stats['mean_diff']):.4f} eV/atom on average")
    
    if gen_results:
        print("\n2. Generated Structures (MatterSim vs VASP-DFT):")
        print("-" * 70)
        stats = gen_results['statistics']
        
        print(f"  Structures compared: {stats['n_structures']}")
        print(f"\n  Energy per Atom:")
        print(f"    Correlation: {stats['correlation']:.4f}")
        print(f"    MAE:  {stats['mae']:.4f} eV/atom")
        print(f"    RMSE: {stats['rmse']:.4f} eV/atom")
        print(f"    Mean Diff (MS - VASP): {stats['mean_diff']:+.4f} eV/atom")
        print(f"    Std Dev: {stats['std_diff']:.4f} eV/atom")
        print(f"    Range: [{stats['min_diff']:+.4f}, {stats['max_diff']:+.4f}] eV/atom")
        
        # Check for systematic bias
        if abs(stats['mean_diff']) > 0.05:
            print(f"\n    WARNING: Systematic bias detected!")
            print(f"    MatterSim is {'higher' if stats['mean_diff'] > 0 else 'lower'} "
                  f"than VASP-DFT by {abs(stats['mean_diff']):.4f} eV/atom on average")
    
    # Check for partial results warning
    if mp_results and mp_results.get('coverage'):
        if mp_results['coverage'] < 50:
            partial_results = True
            partial_msg.append(f"MP phases: Only {mp_results['coverage']:.1f}% coverage (workflow running)")
    
    if gen_results and gen_results.get('coverage'):
        if gen_results['coverage'] < 50:
            partial_results = True
            partial_msg.append(f"Generated structures: Only {gen_results['coverage']:.1f}% coverage (VASP running)")
    
    if partial_results:
        print("\n" + "="*35)
        print("WARNING: PARTIAL RESULTS - Workflow Still Running")
        print("="*35)
        for msg in partial_msg:
            print(f"  {msg}")
        print("\nThese statistics are based on incomplete data.")
        print("For accurate diagnosis, wait for workflow to complete.")
        print("="*35)
    
    # Diagnosis for non-linear hull behavior
    print("\n" + "="*70)
    print("DIAGNOSIS: Non-Linear Hull Behavior")
    print("="*70)
    
    if mp_results and gen_results:
        mp_bias = mp_results['statistics']['per_atom_energy']['mean_diff']
        gen_bias = gen_results['statistics']['mean_diff']
        
        print("\nSystematic Energy Biases:")
        print(f"  MP on-hull phases:     {mp_bias:+.4f} eV/atom (MS - MP DFT)")
        print(f"  Generated structures:  {gen_bias:+.4f} eV/atom (MS - VASP DFT)")
        
        bias_diff = abs(mp_bias - gen_bias)
        if bias_diff > 0.02:
            print(f"\n⚠ CRITICAL: Different biases between reference phases and target structures!")
            print(f"  Bias difference: {bias_diff:.4f} eV/atom")
            print(f"\n  This can cause non-linear hull behavior because:")
            print(f"    - Reference phases (MP) have different MatterSim bias than")
            print(f"    - Target structures (generated)")
            print(f"    - Result: Inconsistent energy reference → incorrect E_hull")
            print(f"\n  Recommendation:")
            print(f"    Use consistent energy method for BOTH reference and target:")
            print(f"    - Option A: MP raw DFT + VASP DFT (compute_dft_e_hull.py)")
            print(f"    - Option B: All MatterSim (prescreening only)")
            print(f"    - AVOID: Mixing MatterSim reference with DFT target (or vice versa)")
        else:
            print(f"\n✓ Biases are consistent ({bias_diff:.4f} eV/atom difference)")
            print(f"  This suggests MatterSim has systematic offset from DFT,")
            print(f"  but it's consistent across different structure types.")
    
    print("\n" + "="*70)


def main():
    parser = argparse.ArgumentParser(
        description="Compare energies from different calculation methods"
    )
    parser.add_argument(
        '--vasp-jobs',
        type=str,
        default='./VASP_JOBS',
        help="VASP jobs directory (default: ./VASP_JOBS)"
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help="Output directory for plots (default: same as vasp-jobs)"
    )
    
    args = parser.parse_args()
    
    vasp_jobs = Path(args.vasp_jobs).expanduser()
    output_dir = Path(args.output_dir).expanduser() if args.output_dir else vasp_jobs
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # File paths
    mp_mattersim_cache = vasp_jobs / "mp_mattersim.json"
    mp_dft_cache = vasp_jobs / "mp_vaspdft.json"
    prescreen_json = vasp_jobs / "prescreening_stability.json"
    workflow_json = vasp_jobs / "workflow.json"
    
    print("="*70)
    print("Energy Method Comparison Script")
    print("="*70)
    print(f"VASP jobs directory: {vasp_jobs}")
    print(f"Output directory: {output_dir}")
    print("\nCache files:")
    print(f"  MatterSim MP cache: {mp_mattersim_cache}")
    print(f"  DFT MP cache: {mp_dft_cache}")
    print(f"  Prescreening results: {prescreen_json}")
    print(f"  Workflow database: {workflow_json}")
    print("="*70)
    
    # Load cache files
    print("\nLoading data files...")
    ms_cache = load_json_file(mp_mattersim_cache)
    dft_cache = load_json_file(mp_dft_cache)
    
    # Compare MP phases
    mp_results = None
    if ms_cache and dft_cache:
        mp_results = compare_mp_phases(ms_cache, dft_cache)
    else:
        print("\nWARNING: Could not load MP cache files - skipping MP phase comparison")
    
    # Compare generated structures
    gen_results = compare_generated_structures(prescreen_json, workflow_json, vasp_jobs)
    
    # Generate plots
    print("\n" + "="*70)
    print("Generating Plots")
    print("="*70)
    
    if mp_results:
        print("\nMP On-Hull Phases:")
        plot_mp_phase_comparison(mp_results, 
                                 output_prefix=str(output_dir / "mp_phases_comparison"))
    
    if gen_results:
        print("\nGenerated Structures:")
        plot_generated_structures_comparison(gen_results,
                                            output_prefix=str(output_dir / "generated_structures_comparison"))
    
    # Print summary
    print_summary(mp_results, gen_results)
    
    # Save detailed results to JSON
    output_json = output_dir / "energy_method_comparison.json"
    results_data = {
        'mp_phases': {
            'statistics': mp_results['statistics'] if mp_results else None,
            'n_matched': len(mp_results['matched']) if mp_results else 0,
            'coverage': mp_results.get('coverage', 0) if mp_results else 0
        },
        'generated_structures': {
            'statistics': gen_results['statistics'] if gen_results else None,
            'n_matched': len(gen_results['matched']) if gen_results else 0,
            'coverage': gen_results.get('coverage', 0) if gen_results else 0
        }
    }
    
    with open(output_json, 'w') as f:
        json.dump(results_data, f, indent=2)
    
    # Print final output summary
    print("\n" + "="*70)
    print("OUTPUT FILES SUMMARY")
    print("="*70)
    
    print("\n  Plots Generated:")
    plots_generated = []
    if mp_results:
        mp_scatter = output_dir / "mp_phases_comparison_scatter.png"
        mp_residual = output_dir / "mp_phases_comparison_residuals.png"
        if mp_scatter.exists():
            plots_generated.append(str(mp_scatter))
            print(f"    {mp_scatter.name}")
        if mp_residual.exists():
            plots_generated.append(str(mp_residual))
            print(f"    {mp_residual.name}")
    
    if gen_results:
        gen_scatter = output_dir / "generated_structures_comparison_scatter.png"
        gen_residual = output_dir / "generated_structures_comparison_residuals.png"
        if gen_scatter.exists():
            plots_generated.append(str(gen_scatter))
            print(f"    {gen_scatter.name}")
        if gen_residual.exists():
            plots_generated.append(str(gen_residual))
            print(f"    {gen_residual.name}")
    
    if not plots_generated:
        print("    No plots generated (insufficient data)")
    
    print(f"\n  Statistics JSON:")
    print(f"    {output_json.name}")
    
    print(f"\n  Output directory: {output_dir}")
    print("="*70)
    
    # Instructions for viewing plots
    if plots_generated:
        print("\n  To view plots:")
        print(f"   cd {output_dir}")
        if len(plots_generated) <= 4:
            for plot in plots_generated:
                print(f"   open {Path(plot).name}")
        else:
            print(f"   open *.png")
    
    print("="*70 + "\n")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

