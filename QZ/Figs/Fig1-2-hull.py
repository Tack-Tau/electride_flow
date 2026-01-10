import matplotlib.pyplot as plt
import numpy as np

# Create figure and axis
fig, ax = plt.subplots(figsize=(5, 3.5))


x_points = np.array([-2, 0, 0.5, 2])
y_points = np.array([4.5, 0.8, 1.3, 4.5])
ax.scatter(x_points, y_points, color='k', s=100)
# Create points for convex hull
points = np.column_stack([x_points, y_points])
# Plot the hull as connected line segments
ax.plot(x_points, y_points, 'b-',
        linewidth=2.5,
        label='Reference MLP Hull')

good_pt = [-0.5, 1.2]
bad_pt = [-.25, 2.2]

# Find intersection of vertical line at bad_pt[0] with the hull
x_bad = bad_pt[0]
# Find which line segment the x-coordinate falls on
for i in range(len(x_points) - 1):
    if x_points[i] <= x_bad <= x_points[i+1] or x_points[i+1] <= x_bad <= x_points[i]:
        # Linear interpolation to find y on the hull
        x1, y1 = x_points[i], y_points[i]
        x2, y2 = x_points[i+1], y_points[i+1]
        y_hull = y1 + (y2 - y1) * (x_bad - x1) / (x2 - x1)
        break

ax.scatter(*bad_pt, color='red', s=150, facecolors='none', edgecolors='red')
ax.text(bad_pt[0], bad_pt[1], 'x', fontsize=16, color='red', ha='center', va='center')
ax.plot([x_bad, x_bad], [bad_pt[1]-0.075, y_hull], color='gray', linestyle='--',
        linewidth=1.5)

x_good = good_pt[0]
for i in range(len(x_points) - 1):
    if x_points[i] <= x_good <= x_points[i+1] or x_points[i+1] <= x_good <= x_points[i]:
        x1, y1 = x_points[i], y_points[i]
        x2, y2 = x_points[i+1], y_points[i+1]
        y_hull = y1 + (y2 - y1) * (x_good - x1) / (x2 - x1)
        break
ax.scatter(*good_pt, color='green', s=150, facecolors='none', edgecolors='green')
ax.text(good_pt[0], good_pt[1], '✓', fontsize=16, color='green', ha='center', va='center')
ax.plot([good_pt[0], good_pt[0]], [good_pt[1]+0.075, y_hull], color='gray', linestyle='--',
        linewidth=1.5)
#ax.text(0, 3.2, '$E_\\text{ref-hull-mlp}$ < 0.1 eV/atom', fontsize=14,
#        ha='center', va='center', fontweight='bold')
# Set labels and title
ax.set_xlabel('Composition', fontsize=16, fontweight='bold')
ax.set_ylabel('Energy', fontsize=16, fontweight='bold')
ax.legend(fontsize=18, frameon=False)

ax.set_xticks([])
ax.set_yticks([])
# Adjust layout and save
plt.tight_layout()
plt.savefig('Fig1-2.pdf', dpi=300, bbox_inches='tight')
#plt.show()
