#!/usr/bin/env python3
"""
Compare convex hulls: VASP-DFT vs MP-API GGA and MatterSim vs MP-API GGA.

For each binary and ternary chemical system, builds three independent convex
hulls using energies from:
  1. MP-API GGA (from mp_vaspdft.json)
  2. VASP-DFT (from mp_vasp_comparison.json -- structures + outliers)
  3. MatterSim (from mp_mattersim.json)

Then compares e_above_hull and on-hull status for each phase across sources.
Binary hulls are plotted as formation-energy vs composition overlays.
Ternary hulls are plotted as triangular phase diagrams with tie-lines.

Outputs:
  hull_vasp_vs_mp.pdf              - VASP-DFT binary hull overlaid on MP-GGA
  hull_mattersim_vs_mp.pdf         - MatterSim binary hull overlaid on MP-GGA
  hull_ternary_vasp_vs_mp.pdf      - VASP-DFT ternary hull overlaid on MP-GGA
  hull_ternary_mattersim_vs_mp.pdf - MatterSim ternary hull overlaid on MP-GGA
  form_energy_vasp_vs_mp.pdf       - VASP vs MP formation energy scatter
  form_energy_mattersim_vs_mp.pdf  - MatterSim vs MP formation energy scatter
  e_above_hull_vasp_vs_mp.pdf     - VASP vs MP energy above hull scatter
  e_above_hull_mattersim_vs_mp.pdf - MatterSim vs MP energy above hull scatter
  hull_comparison.json             - Full numerical comparison data
  hull_vasp_vs_mp.csv              - VASP vs MP tabular comparison
  hull_mattersim_vs_mp.csv         - MatterSim vs MP tabular comparison
"""

import json
import sys
import csv
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
from itertools import combinations

from pymatgen.core import Composition, Element
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry


def parse_comp_dict(comp_dict):
    """Parse composition dict.

    Handles two formats:
      - Element keys: {'Co': 9.0, 'Gd': 12.0} -> value is atom count per element.
      - Formula key:  {'GdCo5': 6} -> value is total atom count, NOT formula units.
    """
    if len(comp_dict) == 1:
        formula, num_atoms = next(iter(comp_dict.items()))
        base = Composition(formula)
        scale = num_atoms / base.num_atoms
        return base * scale
    total = Composition()
    for el, count in comp_dict.items():
        total += Composition(el) * count
    return total


def load_mp_entries(path):
    """Load MP-API GGA entries from mp_vaspdft.json (total energies)."""
    with open(path) as f:
        data = json.load(f)
    entries = []
    for item in data:
        comp = parse_comp_dict(item['composition'])
        entries.append(PDEntry(comp, item['energy'], name=item['mp_id']))
    return entries


def load_vasp_entries(path):
    """Load VASP entries from mp_vasp_comparison.json (per-atom energies * natoms).
    Includes all structures (both normal and outlier-flagged)."""
    with open(path) as f:
        data = json.load(f)
    entries = []
    seen = set()
    for item in data['structures']:
        mp_id = item['mp_id']
        if mp_id in seen:
            continue
        seen.add(mp_id)
        comp = Composition(item['formula'])
        entries.append(PDEntry(comp, item['vasp_energy'] * comp.num_atoms, name=mp_id))
    return entries


def load_mattersim_entries(path):
    """Load MatterSim entries from mp_mattersim.json (total energies)."""
    with open(path) as f:
        data = json.load(f)
    entries = []
    for item in data:
        comp = parse_comp_dict(item['composition'])
        entries.append(PDEntry(comp, item['energy'], name=item['mp_id']))
    return entries


def filter_for_chemsys(entries, chemsys):
    """Keep entries whose elements are a subset of chemsys elements."""
    elements = set(chemsys.split('-'))
    return [e for e in entries
            if set(str(el) for el in e.composition.elements).issubset(elements)]


def safe_build_pd(entries, chemsys):
    """Build PhaseDiagram; returns (pd, filtered_entries) or (None, [])."""
    sub = filter_for_chemsys(entries, chemsys)
    elements = set(chemsys.split('-'))
    found = set()
    for e in sub:
        for el in e.composition.elements:
            found.add(str(el))
    missing = elements - found
    if missing:
        print(f"  WARNING: missing elements {sorted(missing)} for {chemsys}")
        return None, sub
    try:
        pd = PhaseDiagram(sub)
        return pd, sub
    except Exception as exc:
        print(f"  ERROR building PD for {chemsys}: {exc}")
        return None, sub


def get_hull_info(pd, entries):
    """Compute e_above_hull for all entries. Returns {mp_id: info_dict}."""
    if pd is None:
        return {}
    result = {}
    for entry in entries:
        try:
            eah = pd.get_e_above_hull(entry)
            fe = pd.get_form_energy_per_atom(entry)
        except Exception:
            continue
        result[entry.name] = {
            'e_above_hull': round(float(eah), 6),
            'form_energy': round(float(fe), 6),
            'formula': entry.composition.reduced_formula,
            'on_hull': bool(eah < 1e-4),
        }
    return result


def detect_binary_systems(entries):
    """Auto-detect binary chemical systems from entry compositions."""
    chemsys_set = set()
    for e in entries:
        els = sorted(str(el) for el in e.composition.elements)
        if len(els) == 2:
            chemsys_set.add('-'.join(els))
    return sorted(chemsys_set)


def detect_ternary_systems(entries):
    """Auto-detect ternary chemical systems from entry compositions."""
    chemsys_set = set()
    for e in entries:
        els = sorted(str(el) for el in e.composition.elements)
        if len(els) == 3:
            chemsys_set.add('-'.join(els))
    return sorted(chemsys_set)


def print_comparison_table(chemsys, mp_info, alt_info, alt_name):
    """Print side-by-side hull comparison.

    Returns:
        changes: list of (mp_id, formula, change_type)
        rows: list of dicts for CSV output
    """
    all_ids = sorted(set(mp_info.keys()) | set(alt_info.keys()))

    w = 120
    print(f"\n{'='*w}")
    print(f"  {chemsys}: MP-GGA vs {alt_name}")
    print(f"{'='*w}")
    print(f"{'mp_id':<16} {'formula':<12} "
          f"{'MP Ef':>8} {'MP Ehull':>10} {'':>4} "
          f"{alt_name+' Ef':>12} {alt_name+' Ehull':>12} {'':>4}  change")
    print(f"{'':-<16} {'':-<12} "
          f"{'':-<8} {'':-<10} {'':-<4} "
          f"{'':-<12} {'':-<12} {'':-<4}  {'':-<8}")

    changes = []
    rows = []
    for mp_id in all_ids:
        mp = mp_info.get(mp_id)
        alt = alt_info.get(mp_id)
        formula = (mp or alt or {}).get('formula', '?')

        mp_fe = f"{mp['form_energy']:.4f}" if mp else "  ---"
        mp_eah = f"{mp['e_above_hull']:.4f}" if mp else "  ---"
        mp_tag = "ON" if mp and mp['on_hull'] else ("off" if mp else "---")
        alt_fe = f"{alt['form_energy']:.4f}" if alt else "  ---"
        alt_eah = f"{alt['e_above_hull']:.4f}" if alt else "  ---"
        alt_tag = "ON" if alt and alt['on_hull'] else ("off" if alt else "---")

        changed = ""
        if mp and alt:
            if mp['on_hull'] and not alt['on_hull']:
                changed = "LOST"
                changes.append((mp_id, formula, changed))
            elif not mp['on_hull'] and alt['on_hull']:
                changed = "GAINED"
                changes.append((mp_id, formula, changed))
        elif mp and not alt:
            changed = "MISSING"
        elif alt and not mp:
            changed = "NEW"

        print(f"{mp_id:<16} {formula:<12} "
              f"{mp_fe:>8} {mp_eah:>10} ({mp_tag:>3}) "
              f"{alt_fe:>12} {alt_eah:>12} ({alt_tag:>3})  {changed}")

        rows.append({
            'chemsys': chemsys,
            'mp_id': mp_id,
            'formula': formula,
            'mp_form_energy': mp['form_energy'] if mp else '',
            'mp_e_above_hull': mp['e_above_hull'] if mp else '',
            'mp_on_hull': mp['on_hull'] if mp else '',
            f'{alt_name}_form_energy': alt['form_energy'] if alt else '',
            f'{alt_name}_e_above_hull': alt['e_above_hull'] if alt else '',
            f'{alt_name}_on_hull': alt['on_hull'] if alt else '',
            'change': changed,
        })

    if changes:
        print(f"\nHull boundary changes ({len(changes)}):")
        for mp_id, formula, chg in changes:
            print(f"  {chg:>6}: {mp_id} ({formula})")
    else:
        print(f"\nNo hull boundary changes.")

    return changes, rows


def _draw_binary(ax, pd, entries, el_a, el_b, color, label, marker, ls, anno_y):
    """Draw one binary hull + off-hull points on axes."""
    if pd is None:
        return

    hull_pts, off_pts = [], []
    for entry in entries:
        entry_els = set(str(el) for el in entry.composition.elements)
        if not entry_els.issubset({el_a, el_b}):
            continue
        try:
            fe = pd.get_form_energy_per_atom(entry)
            eah = pd.get_e_above_hull(entry)
        except Exception:
            continue
        pt = {
            'x': entry.composition.get_atomic_fraction(Element(el_b)),
            'y': fe,
            'formula': entry.composition.reduced_formula,
        }
        if eah < 1e-4:
            hull_pts.append(pt)
        else:
            off_pts.append(pt)

    hull_pts.sort(key=lambda p: p['x'])

    if hull_pts:
        hx = [p['x'] for p in hull_pts]
        hy = [p['y'] for p in hull_pts]
        ax.plot(hx, hy, ls=ls, color=color, lw=2, label=label, zorder=3)
        ax.scatter(hx, hy, color=color, marker=marker, s=60, zorder=4,
                   edgecolors='black', linewidth=0.5)
        for p in hull_pts:
            if 0.01 < p['x'] < 0.99:
                ax.annotate(p['formula'], (p['x'], p['y']),
                            textcoords='offset points', xytext=(0, anno_y),
                            fontsize=6, ha='center', color=color)

    if off_pts:
        ax.scatter([p['x'] for p in off_pts], [p['y'] for p in off_pts],
                   color=color, marker=marker, s=20, alpha=0.3, zorder=2)


def plot_hull_overlay(mp_entries, alt_entries, alt_name, binary_systems,
                      output_path, alt_color):
    """Plot binary hull comparison as multi-page PDF (2 cols x 3 rows max)."""
    n = len(binary_systems)
    if n == 0:
        return
    ncols = 2
    nrows = 3
    per_page = ncols * nrows

    with PdfPages(output_path) as pdf:
        for page_start in range(0, n, per_page):
            page_systems = binary_systems[page_start:page_start + per_page]
            n_page = len(page_systems)
            rows_needed = (n_page + ncols - 1) // ncols
            fig, axes = plt.subplots(rows_needed, ncols,
                                     figsize=(6.5 * ncols, 5 * rows_needed),
                                     squeeze=False)

            for idx, chemsys in enumerate(page_systems):
                row, col = divmod(idx, ncols)
                ax = axes[row][col]
                els = sorted(chemsys.split('-'))
                el_a, el_b = els[0], els[1]

                mp_pd, mp_sub = safe_build_pd(mp_entries, chemsys)
                alt_pd, alt_sub = safe_build_pd(alt_entries, chemsys)

                _draw_binary(ax, mp_pd, mp_sub, el_a, el_b,
                             color='#1f77b4', label='MP-GGA', marker='s',
                             ls='--', anno_y=10)
                _draw_binary(ax, alt_pd, alt_sub, el_a, el_b,
                             color=alt_color, label=alt_name, marker='o',
                             ls='-', anno_y=-14)

                ax.set_xlabel(f'x({el_b})  in  {el_a}$_{{1-x}}${el_b}$_x$',
                              fontsize=10)
                ax.set_ylabel('Formation energy (eV/atom)', fontsize=10)
                ax.set_title(f'{el_a}\u2013{el_b}', fontsize=12,
                             fontweight='bold')
                ax.axhline(0, color='gray', lw=0.5)
                ax.legend(fontsize=8, loc='best')
                ax.set_xlim(-0.05, 1.05)

            for idx in range(n_page, rows_needed * ncols):
                row, col = divmod(idx, ncols)
                axes[row][col].set_visible(False)

            page_num = page_start // per_page + 1
            total_pages = (n + per_page - 1) // per_page
            fig.suptitle(f'Convex Hull: {alt_name} vs MP-GGA '
                         f'(page {page_num}/{total_pages})',
                         fontsize=14, fontweight='bold')
            fig.tight_layout(rect=[0, 0, 1, 0.96])
            pdf.savefig(fig, dpi=150, bbox_inches='tight')
            plt.close(fig)

    print(f"  Saved: {output_path} ({total_pages} pages)")


def _ternary_coords(comp, el_b, el_c):
    """Convert composition to ternary (x, y) in equilateral triangle.

    Vertices: el_a at (0,0), el_b at (1,0), el_c at (0.5, sqrt(3)/2).
    """
    fb = comp.get_atomic_fraction(Element(el_b))
    fc = comp.get_atomic_fraction(Element(el_c))
    return fb + 0.5 * fc, fc * np.sqrt(3) / 2


def _setup_ternary_axes(ax, el_a, el_b, el_c):
    """Draw composition triangle and element labels."""
    h = np.sqrt(3) / 2
    ax.plot([0, 1, 0.5, 0], [0, 0, h, 0], 'k-', lw=1.5, zorder=1)
    off = 0.06
    ax.text(0 - off, -off, el_a, ha='center', va='top',
            fontsize=11, fontweight='bold')
    ax.text(1 + off, -off, el_b, ha='center', va='top',
            fontsize=11, fontweight='bold')
    ax.text(0.5, h + off, el_c, ha='center', va='bottom',
            fontsize=11, fontweight='bold')
    ax.set_xlim(-0.15, 1.15)
    ax.set_ylim(-0.15, h + 0.15)
    ax.set_aspect('equal')
    ax.axis('off')


def _draw_ternary(ax, pd, entries, el_a, el_b, el_c,
                  color, label, marker, ls, anno_offset):
    """Draw one ternary hull (tie-lines + points) on axes."""
    if pd is None:
        return

    target_els = {el_a, el_b, el_c}
    hull_pts, off_pts = [], []
    for entry in entries:
        entry_els = set(str(el) for el in entry.composition.elements)
        if not entry_els.issubset(target_els):
            continue
        try:
            eah = pd.get_e_above_hull(entry)
        except Exception:
            continue
        x, y = _ternary_coords(entry.composition, el_b, el_c)
        pt = {'x': x, 'y': y,
              'formula': entry.composition.reduced_formula,
              'n_els': len(entry.composition.elements)}
        if eah < 1e-4:
            hull_pts.append(pt)
        else:
            off_pts.append(pt)

    # Tie-lines from hull facets
    for facet in pd.facets:
        coords = [_ternary_coords(pd.qhull_entries[i].composition,
                                  el_b, el_c) for i in facet]
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                ax.plot([coords[i][0], coords[j][0]],
                        [coords[i][1], coords[j][1]],
                        ls=ls, color=color, lw=1.5, alpha=0.5, zorder=2)

    if hull_pts:
        hx = [p['x'] for p in hull_pts]
        hy = [p['y'] for p in hull_pts]
        ax.scatter(hx, hy, color=color, marker=marker, s=60, zorder=4,
                   edgecolors='black', linewidth=0.5, label=label)
        for p in hull_pts:
            if p['n_els'] > 1:
                ax.annotate(p['formula'], (p['x'], p['y']),
                            textcoords='offset points', xytext=anno_offset,
                            fontsize=5, ha='center', color=color)

    if off_pts:
        ax.scatter([p['x'] for p in off_pts], [p['y'] for p in off_pts],
                   color=color, marker=marker, s=20, alpha=0.3, zorder=2)


def plot_ternary_overlay(mp_entries, alt_entries, alt_name, ternary_systems,
                         output_path, alt_color):
    """Plot ternary hull comparison as multi-page PDF (2 cols x 3 rows max)."""
    n = len(ternary_systems)
    if n == 0:
        return
    ncols = 2
    nrows = 3
    per_page = ncols * nrows

    with PdfPages(output_path) as pdf:
        for page_start in range(0, n, per_page):
            page_systems = ternary_systems[page_start:page_start + per_page]
            n_page = len(page_systems)
            rows_needed = (n_page + ncols - 1) // ncols
            fig, axes = plt.subplots(rows_needed, ncols,
                                     figsize=(7 * ncols, 6.5 * rows_needed),
                                     squeeze=False)

            for idx, chemsys in enumerate(page_systems):
                row, col = divmod(idx, ncols)
                ax = axes[row][col]
                els = sorted(chemsys.split('-'))
                el_a, el_b, el_c = els

                _setup_ternary_axes(ax, el_a, el_b, el_c)

                mp_pd, mp_sub = safe_build_pd(mp_entries, chemsys)
                alt_pd, alt_sub = safe_build_pd(alt_entries, chemsys)

                _draw_ternary(ax, mp_pd, mp_sub, el_a, el_b, el_c,
                              color='#1f77b4', label='MP-GGA', marker='s',
                              ls='--', anno_offset=(0, 8))
                _draw_ternary(ax, alt_pd, alt_sub, el_a, el_b, el_c,
                              color=alt_color, label=alt_name, marker='o',
                              ls='-', anno_offset=(0, -10))

                ax.set_title(f'{el_a}\u2013{el_b}\u2013{el_c}',
                             fontsize=12, fontweight='bold')
                ax.legend(fontsize=8, loc='upper right')

            for idx in range(n_page, rows_needed * ncols):
                row, col = divmod(idx, ncols)
                axes[row][col].set_visible(False)

            page_num = page_start // per_page + 1
            total_pages = (n + per_page - 1) // per_page
            fig.suptitle(f'Ternary Hull: {alt_name} vs MP-GGA '
                         f'(page {page_num}/{total_pages})',
                         fontsize=14, fontweight='bold')
            fig.tight_layout(rect=[0, 0, 1, 0.96])
            pdf.savefig(fig, dpi=150, bbox_inches='tight')
            plt.close(fig)

    print(f"  Saved: {output_path} ({total_pages} pages)")


def plot_hull_scatter(rows, alt_name, mp_key, alt_key, quantity_label,
                      output_path, alt_color):
    """Plot MP vs alt scatter with parity line and MAE/RMSE metrics."""
    mp_vals, alt_vals = [], []
    seen = set()
    for r in rows:
        mp_id = r.get('mp_id', '')
        if mp_id in seen:
            continue
        seen.add(mp_id)
        mp_val = r.get(mp_key, '')
        alt_val = r.get(alt_key, '')
        if mp_val == '' or alt_val == '':
            continue
        mp_vals.append(float(mp_val))
        alt_vals.append(float(alt_val))

    if len(mp_vals) < 2:
        return

    mp_vals = np.array(mp_vals)
    alt_vals = np.array(alt_vals)
    mae = np.mean(np.abs(mp_vals - alt_vals))
    rmse = np.sqrt(np.mean((mp_vals - alt_vals) ** 2))

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(mp_vals, alt_vals, s=20, alpha=0.6, color=alt_color,
               edgecolors='none')

    lo = min(mp_vals.min(), alt_vals.min())
    hi = max(mp_vals.max(), alt_vals.max())
    margin = 0.05 * max(hi - lo, 1e-6)
    ax.plot([lo - margin, hi + margin], [lo - margin, hi + margin],
            'k--', lw=1, zorder=0)
    ax.set_xlim(lo - margin, hi + margin)
    ax.set_ylim(lo - margin, hi + margin)

    ax.set_xlabel(f'MP-GGA {quantity_label} (eV/atom)', fontsize=11)
    ax.set_ylabel(f'{alt_name} {quantity_label} (eV/atom)', fontsize=11)
    ax.set_title(f'{quantity_label}: {alt_name} vs MP-GGA', fontsize=13,
                 fontweight='bold')
    ax.set_aspect('equal')

    textstr = (f'N = {len(mp_vals)}\n'
               f'MAE = {mae:.4f} eV/atom\n'
               f'RMSE = {rmse:.4f} eV/atom')
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=9,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Compare VASP / MatterSim / MP-GGA convex hulls"
    )
    parser.add_argument('--comparison-json', default='mp_vasp_comparison.json',
                        help="mp_vasp_comparison.json path")
    parser.add_argument('--mp-json', default='mp_vaspdft.json',
                        help="mp_vaspdft.json path (MP-API GGA entries)")
    parser.add_argument('--mattersim-json', default='mp_mattersim.json',
                        help="mp_mattersim.json path (MatterSim energies)")
    parser.add_argument('--output-dir', default='.',
                        help="Directory for output plots and JSON")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    mp_entries = load_mp_entries(args.mp_json)
    vasp_entries = load_vasp_entries(args.comparison_json)
    ms_entries = load_mattersim_entries(args.mattersim_json)
    print(f"  MP-GGA entries:   {len(mp_entries)}")
    print(f"  VASP entries:     {len(vasp_entries)}")
    print(f"  MatterSim entries:{len(ms_entries)}")

    binary_systems = detect_binary_systems(mp_entries)
    ternary_systems = detect_ternary_systems(mp_entries)
    all_systems = binary_systems + ternary_systems
    print(f"\nBinary systems detected:  {binary_systems}")
    print(f"Ternary systems detected: {ternary_systems}\n")

    # ------- VASP vs MP -------
    print("=" * 85)
    print("  PART 1: VASP-DFT vs MP-GGA")
    print("=" * 85)

    vasp_changes_all = []
    vasp_rows_all = []
    vasp_json = {}
    for chemsys in all_systems:
        mp_pd, mp_sub = safe_build_pd(mp_entries, chemsys)
        vasp_pd, vasp_sub = safe_build_pd(vasp_entries, chemsys)
        mp_info = get_hull_info(mp_pd, mp_sub)
        vasp_info = get_hull_info(vasp_pd, vasp_sub)
        changes, rows = print_comparison_table(chemsys, mp_info, vasp_info, "VASP")
        vasp_changes_all.extend([(chemsys, *c) for c in changes])
        vasp_rows_all.extend(rows)
        vasp_json[chemsys] = {'mp': mp_info, 'vasp': vasp_info}

    # ------- MatterSim vs MP -------
    print("\n" + "=" * 85)
    print("  PART 2: MatterSim vs MP-GGA")
    print("=" * 85)

    ms_changes_all = []
    ms_rows_all = []
    ms_json = {}
    for chemsys in all_systems:
        mp_pd, mp_sub = safe_build_pd(mp_entries, chemsys)
        ms_pd, ms_sub = safe_build_pd(ms_entries, chemsys)
        mp_info = get_hull_info(mp_pd, mp_sub)
        ms_info = get_hull_info(ms_pd, ms_sub)
        changes, rows = print_comparison_table(chemsys, mp_info, ms_info, "MatterSim")
        ms_changes_all.extend([(chemsys, *c) for c in changes])
        ms_rows_all.extend(rows)
        ms_json[chemsys] = {'mp': mp_info, 'mattersim': ms_info}

    # ------- Plots -------
    print("\n\nGenerating hull comparison plots...")
    plot_hull_overlay(mp_entries, vasp_entries, 'VASP-DFT', binary_systems,
                      out / 'hull_vasp_vs_mp.pdf', alt_color='#d62728')
    plot_hull_overlay(mp_entries, ms_entries, 'MatterSim', binary_systems,
                      out / 'hull_mattersim_vs_mp.pdf', alt_color='#2ca02c')
    plot_ternary_overlay(mp_entries, vasp_entries, 'VASP-DFT', ternary_systems,
                         out / 'hull_ternary_vasp_vs_mp.pdf', alt_color='#d62728')
    plot_ternary_overlay(mp_entries, ms_entries, 'MatterSim', ternary_systems,
                         out / 'hull_ternary_mattersim_vs_mp.pdf', alt_color='#2ca02c')

    # ------- Save JSON -------
    comparison_out = {
        'vasp_vs_mp': vasp_json,
        'mattersim_vs_mp': ms_json,
        'vasp_hull_changes': [
            {'chemsys': cs, 'mp_id': mid, 'formula': f, 'change': ch}
            for cs, mid, f, ch in vasp_changes_all
        ],
        'mattersim_hull_changes': [
            {'chemsys': cs, 'mp_id': mid, 'formula': f, 'change': ch}
            for cs, mid, f, ch in ms_changes_all
        ],
    }
    json_path = out / 'hull_comparison.json'
    with open(json_path, 'w') as f:
        json.dump(comparison_out, f, indent=2)
    print(f"  Saved: {json_path}")

    # ------- Save CSVs -------
    def write_csv(rows, path):
        if not rows:
            return
        with open(path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        print(f"  Saved: {path}")

    write_csv(vasp_rows_all, out / 'hull_vasp_vs_mp.csv')
    write_csv(ms_rows_all, out / 'hull_mattersim_vs_mp.csv')

    # ------- Scatter plots -------
    print("\nGenerating scatter plots...")
    plot_hull_scatter(vasp_rows_all, 'VASP',
                      'mp_form_energy', 'VASP_form_energy',
                      'Formation Energy',
                      out / 'form_energy_vasp_vs_mp.pdf', '#d62728')
    plot_hull_scatter(ms_rows_all, 'MatterSim',
                      'mp_form_energy', 'MatterSim_form_energy',
                      'Formation Energy',
                      out / 'form_energy_mattersim_vs_mp.pdf', '#2ca02c')
    plot_hull_scatter(vasp_rows_all, 'VASP',
                      'mp_e_above_hull', 'VASP_e_above_hull',
                      'Energy Above Hull',
                      out / 'e_above_hull_vasp_vs_mp.pdf', '#d62728')
    plot_hull_scatter(ms_rows_all, 'MatterSim',
                      'mp_e_above_hull', 'MatterSim_e_above_hull',
                      'Energy Above Hull',
                      out / 'e_above_hull_mattersim_vs_mp.pdf', '#2ca02c')

    # ------- Summary -------
    print("\n" + "=" * 85)
    print("  Summary")
    print("=" * 85)
    print(f"VASP hull boundary changes: {len(vasp_changes_all)}")
    for cs, mid, formula, chg in vasp_changes_all:
        print(f"  {cs}: {chg:>6} {mid} ({formula})")
    print(f"\nMatterSim hull boundary changes: {len(ms_changes_all)}")
    for cs, mid, formula, chg in ms_changes_all:
        print(f"  {cs}: {chg:>6} {mid} ({formula})")
    print(f"\nOutput directory: {out}")
    print("=" * 85)
    return 0


if __name__ == '__main__':
    sys.exit(main())
