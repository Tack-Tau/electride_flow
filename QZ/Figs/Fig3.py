import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FormatStrFormatter
try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except Exception:
    adjust_text = None
    HAS_ADJUSTTEXT = False

df1 = pd.read_csv('YN_full_hull_data.csv')
df2 = pd.read_csv('fig2-ternary.csv')

fig = plt.figure(figsize=(16, 8))
gs = GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.3)

def plot_hull(ax, df, e_max, 
              add_label=False):
    
    xs = df['x'].values
    ys = df['y'].values
    labels = df['label'].values
    on_hulls = df['on_hull'].values
    e_hulls = df['e_above_hull'].values
    electrides = df['is_electride'].values
    formulas = df['formula'].values

    # Plot MP hull
    mp_hulls = []
    for _x, _y, label, formula, on_hull in zip(xs, ys, labels, formulas, on_hulls):
        if label.startswith('mp') and on_hull == 'yes':
            mp_hulls.append((_x, _y, label, formula))
    print(f"{name} MP Hulls:", len(mp_hulls))

    # Plot refined hull
    new_hulls = []
    for _x, _y, label, on_hull in zip(xs, ys, labels, on_hulls):
        if on_hull == 'yes':
            new_hulls.append((_x, _y, label))
    print(f"{name} New Hulls:", len(new_hulls))

    # sort mp_hulls by x value
    mp_hulls = sorted(mp_hulls, key=lambda x: x[0])
    mp_x = [hull[0] for hull in mp_hulls]
    mp_y = [hull[1] for hull in mp_hulls]
    mp_formulas = [hull[3] for hull in mp_hulls]
    ax.plot(mp_x, mp_y, '--', label='MP Hull', 
             color='gray', linewidth=1)
    for mp_xi, mp_yi, mp_formula in zip(mp_x, mp_y, mp_formulas):
        if mp_xi < 0.01 or mp_xi > 0.99:  # only label pure elements
            mp_yi += 0.025 * -ys.min()  # small vertical offset for better visibility
            ha = 'center'
        else:
            mp_xi += 0.058  # small vertical offset for better visibility
            mp_yi -= 0.08
            ha = 'center'
        ax.text(mp_xi, mp_yi, mp_formula.split('_')[0], fontsize=11, 
                alpha=0.8, ha=ha, va='bottom')

    new_hulls = sorted(new_hulls, key=lambda x: x[0])
    new_x = [hull[0] for hull in new_hulls]
    new_y = [hull[1] for hull in new_hulls]
    ax.plot(new_x, new_y, '-', label='Refined Hull', color='black', linewidth=2)
    print(f"{name} Refined Hulls:", len(new_hulls))

    # Plot mp stable points (filled circle)
    pcs_mp = ax.scatter(mp_x, mp_y, c='gray', s=60, alpha=0.7, label='MP Stable',
                 marker='o', edgecolors='none')

    # Plot non-mp stable points (filled square)
    new_electrides = []
    for _x, _y, label, electrides, e, on_hull in zip(xs, ys, labels, electrides, e_hulls, on_hulls):
        if electrides == 'yes' and not label.startswith('mp') and e <= e_max:
            new_electrides.append((_x, _y, label, e, on_hull))
    new_electrides = sorted(new_electrides, key=lambda x: x[0])
    new_electrides_x = [point[0] for point in new_electrides]
    new_electrides_y = [point[1] for point in new_electrides]
    if len(new_electrides) > 0:
        pcs_elec = ax.scatter(new_electrides_x, new_electrides_y, c='none', s=60, alpha=0.7, 
                label='New Electrides', marker='s', edgecolors='red', linewidths=1.5)
        for _x, _y, label, e, on_hull in new_electrides:
            if on_hull == 'yes':
                ax.scatter(_x, _y, c='red', s=60, alpha=0.7, marker='s', edgecolors='none')
        if add_label:
            texts = []
            # dynamic vertical offset for clearer separation
            y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
            if len(new_electrides) < 3:
                y_offset = 0.02 * y_range
            else:
                y_offset = 0.1 * y_range 
            x_offset = 0.1  # small horizontal offset for better visibility
            
            # Calculate y-spacing to detect clusters
            y_coords = [point[1] for point in new_electrides]
            
            for i, (_x, _y, label, e, on_hull) in enumerate(new_electrides):
                label = label.split('_')[0]
                #"CaP1" should become "CaP"
                if label.endswith('1') and label[-2].isalpha():
                    label = label[:-1]
                
                y_pos = _y + y_offset
                if label == 'oS30-Y8N7':
                    y_pos -= 0.10 * y_range  # manual adjustment for this specific label
                    _x -= 0.01  # shift left for better visibility
                if label == 'oP15-Y8N7':
                    y_pos -= 0.10 * y_range  # manual adjustment for this specific label
                    _x -= 0.01  # shift left for better visibility
                if label == 'hR21-Y4N3':
                    _x += 0.1  # shift left for better visibility
                if label == 'hP16-Ca5P3':
                    y_pos -= 0.08 * y_range
                    _x += 0.24  # shift left for better visibility
                if label == 'hR51-Y9N8':
                    y_pos -= 0.08 * y_range
                    _x -= 0.02  # shift left for better visibility
                texts.append(ax.text(_x-0.01, y_pos, label, fontsize=11, 
                                     alpha=0.8,ha='right', va='bottom'))

    new_non_electrides = []
    for _x, _y, label, electrides, e, on_hull in zip(xs, ys, labels, electrides, e_hulls, on_hulls):
        if electrides == 'no' and not label.startswith('mp') and e <= e_max:
            new_non_electrides.append((_x, _y, label, e, on_hull))
    new_non_electrides = sorted(new_non_electrides, key=lambda x: x[0])
    new_non_electrides_x = [point[0] for point in new_non_electrides]
    new_non_electrides_y = [point[1] for point in new_non_electrides]
    if len(new_non_electrides) > 0:
        pcs_non_elec = ax.scatter(new_non_electrides_x, new_non_electrides_y, c='orange', s=60, alpha=0.7, 
                    label='New Non-Electrides', marker='s', edgecolors='none')
        for _x, _y, label, e, on_hull in new_non_electrides:
            if on_hull == 'yes':
                ax.scatter(_x, _y, c='orange', s=60, alpha=0.7, marker='s', edgecolors='none')


for col, name in enumerate(['Ca-P_0.1', 'Y-N_0.02']):
    system = name.split('_')[0]
    A=system.split('-')[0]
    B=system.split('-')[1]
    e_max = float(name.split('_')[1])
    df1 = pd.read_csv(f'{A}{B}_full_hull_data.csv')
    df2 = pd.read_csv(f'{A}{B}_zoom_hull_data.csv')


    
    ax1 = fig.add_subplot(gs[0, col])
    mark = '(a)' if col == 0 else '(b)'

    ax1.set_title(f'{mark} {A}-{B}', fontsize=16, fontweight='bold')
    ax1.text(0.5, 0.95, rf'$E_\text{{hull}} < {e_max:.3f}$ eV/atom', fontsize=14, 
             ha='center', va='top', transform=ax1.transAxes)

    label_formula = rf'$\mathrm{{{A}}}_{{1-x}}\mathrm{{{B}}}_{{x}}$'
    ax1.set_xlabel(f'Composition $x$ in {label_formula}', fontsize=15)
    ax1.set_ylabel('$E_\\text{formation}$ (eV/atom)', fontsize=15)
    ax1.tick_params(axis='both', labelsize=14)
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plot_hull(ax1, df1, e_max)
    ax1.set_ylim(min(df1['y'])*1.15, -min(df1['y'])*0.15)


    ax2 = fig.add_subplot(gs[1, col])
    label_formula = rf'$\mathrm{{{A}}}_{{1-x}}(\mathrm{{{A}{B}}})_{{x}}$'
    ax2.set_xlabel(f'Composition $x$ in {label_formula}', fontsize=15)
    ax2.set_ylabel('$E_\\text{formation}$ (eV/atom)', fontsize=15)
    ax2.tick_params(axis='both', labelsize=13)
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plot_hull(ax2, df2, e_max, True)
    ax2.set_ylim(min(df2['y'])*1.15, -min(df2['y'])*0.15)
    ax2.legend(loc=3, fontsize=12, frameon=False)

plt.tight_layout()
plt.savefig('fig3.pdf', dpi=300)
#plt.show()