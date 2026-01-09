import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from math import gcd
from functools import reduce

# Create a figure to have 3 nested for x, y, z loops and each from 0 to 21.
# Draw the intersection plane on (100), (010), (001).
# Show the region where 0 < (1*x - 3*y - 3*z) <= 4

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111, projection='3d')

# Set view angle for better visualization
ax.view_init(elev=20, azim=45)

# Define the tetrahedron vertices (x+y+z <= 20)
vertices = np.array([
    [0, 0, 0],
    [20, 0, 0],
    [0, 20, 0],
    [0, 0, 20]
])

# Draw tetrahedron edges
edges = [
    [vertices[0], vertices[1]],
    [vertices[0], vertices[2]],
    [vertices[0], vertices[3]],
    [vertices[1], vertices[2]],
    [vertices[1], vertices[3]],
    [vertices[2], vertices[3]]
]

for edge in edges:
    points = np.array(edge)
    ax.plot3D(points[:, 0], points[:, 1], points[:, 2], 'k-', linewidth=2, alpha=0.6)

# Draw the three coordinate planes with better styling
# Plane (100): x = 0
yy, zz = np.meshgrid(np.linspace(0, 20, 30), np.linspace(0, 20, 30))
xx = np.zeros_like(yy)
zz_clip = np.where(yy + zz <= 20, zz, np.nan)
ax.plot_surface(xx, yy, zz_clip, alpha=0.15, color='lightblue', edgecolor='none')

# Plane (010): y = 0
xx, zz = np.meshgrid(np.linspace(0, 20, 30), np.linspace(0, 20, 30))
yy = np.zeros_like(xx)
zz_clip = np.where(xx + zz <= 20, zz, np.nan)
ax.plot_surface(xx, yy, zz_clip, alpha=0.15, color='lightblue', edgecolor='none')

# Plane (001): z = 0
xx, yy = np.meshgrid(np.linspace(0, 20, 30), np.linspace(0, 20, 30))
zz = np.zeros_like(xx)
yy_clip = np.where(xx + yy <= 20, yy, np.nan)
ax.plot_surface(xx, yy_clip, zz, alpha=0.15, color='lightblue', edgecolor='none')

# Highlight the region where 0 < (x - 3*y - 3*z) <= 4
# Draw multiple slices to fill the region with smooth edges
for v in np.linspace(0.1, 4, 40):
    xx, yy = np.meshgrid(np.linspace(0, 20, 200), np.linspace(0, 20, 200))
    zz = (xx - v - 3*yy) / 3
    # Apply constraints: z in [0, 20], x + y + z <= 20
    # Use continuous masking for smooth edges
    mask = (zz >= 0) & (zz <= 20) & (xx + yy + zz <= 20)
    zz[~mask] = np.nan
    ax.plot_surface(xx, yy, zz, alpha=1.0, color='purple', edgecolor='none', shade=True, antialiased=True, rcount=100, ccount=100)

# Draw boundary planes for the highlighted region (v=0 and v=4)
for v, alpha_val in [(0, 0.4), (4, 0.4)]:
    xx, yy = np.meshgrid(np.linspace(0, 20, 200), np.linspace(0, 20, 200))
    zz = (xx - v - 3*yy) / 3
    mask = (zz >= 0) & (zz <= 20) & (xx + yy + zz <= 20)
    zz[~mask] = np.nan
    ax.plot_surface(xx, yy, zz, alpha=alpha_val, color='darkviolet', edgecolor='none', antialiased=True, rcount=100, ccount=100)

# Remove axis and grids
ax.set_axis_off()

# Set axis limits
ax.set_xlim(0, 20)
ax.set_ylim(0, 20)
ax.set_zlim(0, 20)

# Add labels at the vertices
ax.text(20.1, -3, -4, '$N_\\text{Cs}$', fontsize=15, fontweight='bold')
ax.text(0, 19, -3, '$N_\\text{Al}$', fontsize=15, fontweight='bold')
ax.text(0, -1.5, 21, '$N_\\text{P}$', fontsize=15, fontweight='bold')
ax.text(15, 0, -5, 'Electride region', fontsize=13, fontweight='bold', color='darkviolet')
#plt.tight_layout()

plt.savefig('Fig1-1-comp.pdf', dpi=300, bbox_inches='tight')
#plt.show()
