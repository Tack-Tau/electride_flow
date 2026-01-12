import matplotlib.pyplot as plt
import numpy as np


def get_hull_y(x_points, y_points, x):
    for i in range(len(x_points) - 1):
        if x_points[i] <= x <= x_points[i+1] or x_points[i+1] <= x <= x_points[i]:
            x1, y1 = x_points[i], y_points[i]
            x2, y2 = x_points[i+1], y_points[i+1]
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
    return None
    
# Create figure and axis
fig, ax = plt.subplots(figsize=(5, 3.5))

x_points = np.array([-2, 0, 0.5, 2])
y_points = np.array([4.5, 0.8, 1.3, 4.5])
good_pt = [-0.5, 1.2]
switch_pt = [-0.75, 2.4]
bad_pt = [-.25, 2.2]
dy = 0.025
color = 'plum'

x_points_refine = np.array([-2, -0.5, 0, 0.5, 2])
y_points_refine = np.array([4.5, 1.2, 0.8, 1.3, 4.5])

ymin, ymax = y_points.min()*0.8, y_points.max()*1.1
electride_region = plt.Rectangle((-1.0, ymin), 1, ymax*0.7, 
                                 color=color, alpha=0.5)
ax.add_patch(electride_region)

ax.scatter(x_points, y_points, color='k', s=100)
# Create points for convex hull
points = np.column_stack([x_points, y_points])
# Plot the hull as connected line segments
ax.plot(x_points, y_points, 'b-',
        linewidth=2.5,
        label='Reference MLP Hull')



x, y = bad_pt
ax.scatter(x, y, color='red', s=150, facecolors='none', edgecolors='red')
ax.text(x, y-0.01, 'x', fontsize=16, color='red', ha='center', va='center')
y_hull = get_hull_y(x_points, y_points, x)
ax.plot([x, x], [y-dy*4.5, y_hull], color='gray', ls='--', lw=1.5)
print(f"Bad point: {bad_pt}, Hull y-value at x={x}: {y_hull-y}")

x, y = good_pt
y_hull = get_hull_y(x_points, y_points, x)
ax.scatter(x, y, color='green', s=150, facecolors='none', 
           edgecolors='green', linewidths=2.5)
ax.text(x, y-dy*4, '✓', fontsize=16, color='green', ha='center')
ax.plot([x, x], [y+dy*3, y_hull], color='gray', ls='--', lw=1.5)
print(f"Good point: {good_pt}, Hull y-value at x={x}: {y_hull-y}")

x, y = switch_pt
y_hull = get_hull_y(x_points, y_points, x)
ax.scatter(x, y, color='green', s=150, facecolors='none', 
           edgecolors='green', linewidths=2.5)
ax.text(x, y-dy*4, '✓', fontsize=16, color='green', ha='center')
ax.plot([x, x], [y-dy*4, y_hull], color='gray', ls='--', lw=1.5)
print(f"Good point: {good_pt}, Hull y-value at x={x}: {y_hull-y}")

ax.set_xlabel('Composition', fontsize=16, fontweight='bold')
ax.set_ylabel('Energy', fontsize=16, fontweight='bold')
ax.legend(fontsize=18, frameon=False)
ax.set_ylim(ymin, ymax)

ax.set_xticks([])
ax.set_yticks([])
plt.tight_layout()
plt.savefig('Fig1-2.pdf', dpi=300, bbox_inches='tight')
plt.close()


fig, ax = plt.subplots(figsize=(5, 3.5))
electride_region = plt.Rectangle((-1.0, ymin), 1, ymax*0.7, 
                                 color=color, alpha=0.5)
ax.add_patch(electride_region)

ax.scatter(x_points, y_points, color='k', s=100)
points = np.column_stack([x_points_refine, y_points_refine])
ax.plot(points[:,0], points[:,1], 'b-', lw=2.5, label='Refined DFT Hull')

x, y = bad_pt
ax.scatter(x, y, color='red', s=150, facecolors='none', edgecolors='red')
ax.text(x, y-0.01, 'x', fontsize=16, color='red', ha='center', va='center')
y_hull = get_hull_y(x_points_refine, y_points_refine, x)
ax.plot([x, x], [y-dy*4.5, y_hull], color='gray', ls='--', lw=1.5)
print(f"Bad point: {bad_pt}, Hull y-value at x={x}: {y_hull-y}")

x, y = good_pt
y_hull = get_hull_y(x_points_refine, y_points_refine, x)
ax.scatter(x, y, color='green', s=150, facecolors='none', 
           edgecolors='green', linewidths=2.5)
ax.text(x, y-dy*4, '✓', fontsize=16, color='green', ha='center')
ax.plot([x, x], [y+dy*3, y_hull], color='gray', ls='--', lw=1.5)
print(f"Good point: {good_pt}, Hull y-value at x={x}: {y_hull-y}")

x, y = switch_pt
y_hull = get_hull_y(x_points_refine, y_points_refine, x)
ax.scatter(x, y, color='green', s=150, facecolors='none', edgecolors='red')
ax.text(x, y-dy*4, 'x', fontsize=16, color='red', ha='center')
ax.plot([x, x], [y-dy*3, y_hull], color='gray', ls='--', lw=1.5)
print(f"Good point: {good_pt}, Hull y-value at x={x}: {y_hull-y}") 


ax.set_xlabel('Composition', fontsize=16, fontweight='bold')
ax.set_ylabel('Energy', fontsize=16, fontweight='bold')
ax.legend(fontsize=18, frameon=False)
ax.set_ylim(ymin, ymax)

ax.set_xticks([])
ax.set_yticks([])
plt.tight_layout()
plt.savefig('Fig1-4.pdf', dpi=300, bbox_inches='tight')
plt.close()