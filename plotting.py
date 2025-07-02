
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

# Label the axes
ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')

# Show the plot
plt.show()

####
same_data = np.loadtxt("./build/same.txt")
diff_data = np.loadtxt("./build/diff.txt")

same_x = np.array([i*(1.0/len(same_data)) for i in range(len(same_data))])
diff_x = np.array([i*(1.0/len(diff_data)) for i in range(len(diff_data))])

plt.plot(same_x, same_data, color='blue', marker='o')
plt.plot(diff_x, diff_data, color='red', marker='o')
plt.grid(True)
plt.show()

