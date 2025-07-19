
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D

#### 1. cost function
same_c = np.loadtxt("./build/same.txt")
diff_c = np.loadtxt("./build/diff.txt")

same_y = np.linspace(0, 1, len(same_c))
diff_y = np.linspace(0, 1, len(diff_c))

# with open("./presentation/diff_c.csv","w") as file:
#     file.write("c,y\n")
#     for i in range(len(diff_c)):
#         if diff_c[i] >= 0 and i % 10 != 0:
#             continue
#         file.write(f"{diff_c[i]},{diff_y[i]}\n")

plt.plot(same_c, same_y, color='blue', marker='o')
plt.plot(diff_c, diff_y, color='red', marker='o')
plt.axhline(0, color='black', linewidth=2)  # Horizontal axis (y=0)
plt.axvline(0, color='black', linewidth=2)  # Vertical axis (x=0)
plt.grid(True)
plt.show()

### 2. sample points
# Load the CSV file
df_points = pd.read_csv('./build/points.csv')
df_planes = pd.read_csv('./build/planes.csv')

# Extract x, y, z columns
n = df_points['n']
p = df_points['p']
x = df_points['x']
y = df_points['y']
z = df_points['z']

p_plane = df_planes['p']
x_plane = df_planes['x']
y_plane = df_planes['y']
z_plane = df_planes['z']

# Normalize p for consistent coloring
norm = Normalize(vmin=np.min(p), vmax=np.max(p))
cmap = cm.viridis
colors = cmap(norm(p))

# Create a 3D plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Plot the planes
x_vals = np.linspace(x.min(), x.max(), 10)
y_vals = np.linspace(y.min(), y.max(), 10)
z_vals = np.linspace(z.min(), z.max(), 10)
ax.set_xlim(x.min(), x.max())
ax.set_ylim(y.min(), y.max())
ax.set_zlim(z.min(), z.max())
for i in range(len(df_planes)):
    if abs(z_plane[i]) < 1e-3:
        continue
    if abs(z_plane[i]) < 0.1:
        Y, Z = np.meshgrid(y_vals, z_vals)
        X = -(y_plane[i] * Y + z_plane[i] * Z) / x_plane[i]
    else:
        X, Y = np.meshgrid(x_vals, y_vals)
        Z = -(x_plane[i] * X + y_plane[i] * Y) / z_plane[i]
    ax.plot_surface(X, Y, Z, alpha=0.2, color=cmap(norm(p_plane[i])))

# Plot the points
ax.scatter(0, 0, 0, color='black', s=100, label='Origin')
ax.scatter(x, y, z, c=p, cmap=cmap, marker='o')

# Draw axis lines through the origin
x_range = [min(x.min(), 0), max(x.max(), 0)]
y_range = [min(y.min(), 0), max(y.max(), 0)]
z_range = [min(z.min(), 0), max(z.max(), 0)]

ax.plot(x_range, [0, 0], [0, 0], color='black')  # X-axis
ax.plot([0, 0], y_range, [0, 0], color='black')  # Y-axis
ax.plot([0, 0], [0, 0], z_range, color='black')  # Z-axis

# Label the axes
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')


plt.show()
