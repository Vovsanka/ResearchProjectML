import numpy as np
import pandas as pd

df_points = pd.read_csv('./build/data.csv')

opt = np.sort(df_points['opt'])
acc = np.sort(df_points['acc'])
time = np.sort(df_points['time'])

print("min,q1,q2,q3,max")
print("opt:", opt[0], opt[3], opt[7], opt[11], opt[14])
print("acc:", acc[0], acc[3], acc[7], acc[11], acc[14])
print("time:", time[0], time[3], time[7], time[11], time[14])