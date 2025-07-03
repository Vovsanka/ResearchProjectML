
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#### 1. cost function
same_c = np.loadtxt("./build/same.txt")
diff_c = np.loadtxt("./build/diff.txt")

same_y = np.array([i*(1.0/len(same_c)) for i in range(len(same_c))])
diff_y = np.array([i*(1.0/len(diff_c)) for i in range(len(diff_c))])

plt.plot(same_c, same_y, color='blue', marker='o')
plt.plot(diff_c, diff_y, color='red', marker='o')
plt.axhline(0, color='black', linewidth=2)  # Horizontal axis (y=0)
plt.axvline(0, color='black', linewidth=2)  # Vertical axis (x=0)
plt.grid(True)
plt.show()

### 2. sample points
# Load the CSV file
df = pd.read_csv('./build/points.csv')
print(df.columns)

# Extract x, y, z columns
p = df['p']
x = df['x']
y = df['y']
z = df['z']

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the points
ax.scatter(x, y, z, c=p, marker='o')


# Draw axis lines through the origin
x_range = [min(x.min(), 0), max(x.max(), 0)]
y_range = [min(y.min(), 0), max(y.max(), 0)]
z_range = [min(z.min(), 0), max(z.max(), 0)]

ax.plot(x_range, [0, 0], [0, 0], color='black')  # X-axis
ax.plot([0, 0], y_range, [0, 0], color='black')  # Y-axis
ax.plot([0, 0], [0, 0], z_range, color='black')  # Z-axis

# Label the axes
ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')

# Show the plot
plt.show()



