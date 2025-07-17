import numpy as np
import pandas as pd

df_points = pd.read_csv('./build/data.csv')

opt = np.sort(df_points['opt'])
acc = np.sort(df_points['acc'])
time = np.sort(df_points['time'])

print("opt_q1,opt_med,opt_q3,acc_q1,acc_med,acc_q3,time_q1,time_med,time_q3")
# print(f"{opt[3]},{opt[7]},{opt[11]},{acc[3]},{acc[7]},{acc[11]},{time[3]},{time[7]},{time[11]}")
print(f"{opt[1]},{opt[3]},{opt[5]},{acc[1]},{acc[3]},{acc[5]},{time[1]},{time[3]},{time[5]}")
