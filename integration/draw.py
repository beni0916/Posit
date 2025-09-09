import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 載入所有資料
df_double = pd.read_csv('lorenz_output_double.csv')
df_mpfr = pd.read_csv('lorenz_output_mpfr.csv')
df_posit = pd.read_csv('lorenz_output_posit.csv')

# 計算軌跡間的歐幾里得距離
# 歐幾里得距離 = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
distance_double_mpfr = np.sqrt((df_double['x'] - df_mpfr['x'])**2 + 
                               (df_double['y'] - df_mpfr['y'])**2 + 
                               (df_double['z'] - df_mpfr['z'])**2)

distance_mpfr_posit = np.sqrt((df_mpfr['x'] - df_posit['x'])**2 + 
                              (df_mpfr['y'] - df_posit['y'])**2 + 
                              (df_mpfr['z'] - df_posit['z'])**2)

# 繪製距離隨時間的變化
fig, ax = plt.subplots(figsize=(12, 6))

ax.plot(distance_double_mpfr, label='Double vs MPFR', color='purple')
ax.plot(distance_mpfr_posit, label='MPFR vs Posit', color='orange')

ax.set_title('Distance Between Trajectories Over Time', fontsize=16)
ax.set_xlabel('Time Steps')
ax.set_ylabel('Euclidean Distance')
ax.set_yscale('log')  # 使用對數坐標可以更清楚地看到初始的微小差異
ax.legend()
ax.grid(True)
plt.show()