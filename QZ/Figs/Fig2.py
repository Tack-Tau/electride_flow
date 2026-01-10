import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl

df1 = pd.read_csv('fig2-binary.csv')
df2 = pd.read_csv('fig2-ternary.csv')

fig = plt.figure(figsize=(18, 12))
gs = GridSpec(2, 2, figure=fig, width_ratios=[1, 1.25], hspace=0.5, wspace=0.4)
threshold = 0.1

for row, df in enumerate([df1, df2]):

    # LEFT plot data
    ax1 = fig.add_subplot(gs[row, 0])
    ax2 = fig.add_subplot(gs[row, 1])
    x1 = df['dft_vals_hull'].values
    y1 = df['mattersim_vals_hull'].values
    x2 = df['vasp_vals_gen'].values
    y2 = df['ms_vals_gen'].values

    all_vals1 = np.concatenate([y1, x1])
    val_min1, val_max1 = all_vals1.min(), all_vals1.max()
    margin1 = (val_max1 - val_min1) * 0.05
    plot_min1, plot_max1 = val_min1 - margin1, val_max1 + margin1

    all_vals2 = np.concatenate([y2, x2])
    val_min2, val_max2 = all_vals2.min(), all_vals2.max()
    margin2 = (val_max2 - val_min2) * 0.05
    plot_min2, plot_max2 = val_min2 - margin2, val_max2 + margin2

    # ===== LEFT: Generated Structures Scatter (E_abs) =====
    xy2 = np.vstack([x2, y2])
    z2 = gaussian_kde(xy2)(xy2)
    z2 = z2 * len(x2)

    colors_list = ['cyan', 'dodgerblue', 'black']
    cmap = LinearSegmentedColormap.from_list("density", colors_list, N=100)

    idx2 = z2.argsort()
    x2_sorted, y2_sorted, z2_sorted = x2[idx2], y2[idx2], z2[idx2]

    scatter1 = ax1.scatter(x2_sorted, y2_sorted, c=z2_sorted, cmap=cmap,
              norm=mpl.colors.LogNorm(), s=20, marker='s',
              edgecolors='none')

    ax1.plot([plot_min2, plot_max2], [plot_min2, plot_max2], 'r--',
        linewidth=2, alpha=0.7, label='Perfect agreement')
    ax1.axvline(x=threshold, color='green', linestyle=':', linewidth=2, alpha=0.6)

    correlation2 = np.corrcoef(y2, x2)[0, 1]
    mae2 = np.mean(np.abs(y2 - x2))

    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.minorticks_on()
    ax1.grid(True, which='minor', alpha=0.15, linestyle=':')
    ax1.legend(loc='lower right', fontsize=18, frameon=False)

    ax1.set_xlim(plot_min2, plot_max2)
    ax1.set_ylim(plot_min2, plot_max2)
    ax1.tick_params(axis='both', which='major', labelsize=18)

    # ===== RIGHT: Hull Comparison Scatter (E_hull) =====
    xy1 = np.vstack([x1, y1])
    z1 = gaussian_kde(xy1)(xy1)
    z1 = z1 * len(x1)

    idx1 = z1.argsort()
    x1_sorted, y1_sorted, z1_sorted = x1[idx1], y1[idx1], z1[idx1]

    scatter2 = ax2.scatter(x1_sorted, y1_sorted, c=z1_sorted, cmap=cmap,
              norm=mpl.colors.LogNorm(), s=20, marker='s',
              edgecolors='none')

    ax2.plot([plot_min1, plot_max1], [plot_min1, plot_max1], 'r--',
        linewidth=2, alpha=0.7, label='Perfect agreement')
    ax2.axvline(x=threshold, color='green', linestyle=':', linewidth=2, alpha=0.6)

    correlation1 = np.corrcoef(x1, y1)[0, 1]
    mae1 = np.mean(np.abs(y1 - x1))

    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.minorticks_on()
    ax2.grid(True, which='minor', alpha=0.15, linestyle=':')
    ax2.legend(loc='lower right', fontsize=18, frameon=False)

    ax2.set_xlim(plot_min1, plot_max1)
    ax2.set_ylim(plot_min1, threshold)
    ax2.tick_params(axis='both', which='major', labelsize=18)

    stats_text1 = (
        f"R2  = {correlation1:.4f}\n"
        f"MAE = {mae1:.4f} eV/atom"
    )

    stats_text2 = (
        f"R2  = {correlation2:.4f}\n"
        f"MAE = {mae2:.4f} eV/atom"
    )

    ax1.text(0.05, 0.8, stats_text2, transform=ax1.transAxes, fontsize=18)
    ax2.text(0.05, 0.8, stats_text1, transform=ax2.transAxes, fontsize=18)
    ax1.set_xlabel('$E_\\text{abs-DFT}$ (eV/atom)', fontsize=20)
    ax1.set_ylabel('$E_\\text{abs-MLP}$ (eV/atom)', fontsize=20)
    ax2.set_xlabel('$E_\\text{ref-hull-DFT}$ (eV/atom)', fontsize=20)
    ax2.set_ylabel('$E_\\text{ref-hull-MLP}$ (eV/atom)', fontsize=20)

# Add row title
fig.text(0.52, 0.92, '(a) Energy comparison for binaries', ha='center', va='top', fontsize=20, fontweight='bold')
fig.text(0.52, 0.46, '(b) Energy comparison for ternaries', ha='center', va='top', fontsize=20, fontweight='bold')
plt.tight_layout()
plt.savefig('fig2.pdf', dpi=300)#, bbox_inches='tight')
plt.close()


