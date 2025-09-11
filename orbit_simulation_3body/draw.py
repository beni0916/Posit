import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation # 用於製作動畫

# --- 設定 ---
# 檔案名稱 (請替換成您的 CSV 檔案名稱)
# 您可以同時讀取 double 和 posit 的模擬結果
double_file = "double.csv"
# posit_file = "3body_simulation_posit.csv" # 如果您有 posit 的結果

# --- 讀取 CSV 檔案 ---
try:
    df_double = pd.read_csv(double_file)
    # 如果您有 Posit 的結果，也讀取它
    # df_posit = pd.read_csv(posit_file)
except FileNotFoundError:
    print(f"Error: One or more CSV files not found. Please ensure '{double_file}' (and '{posit_file}') exist.")
    exit()

# --- 繪製靜態軌跡圖 ---
plt.figure(figsize=(10, 8)) # 設定圖形大小

# 繪製太陽 (通常標示為一個點，因為它相對不動)
plt.plot(df_double['x_sun'], df_double['y_sun'], 'o', markersize=5, color='yellow', label='Sun')

# 繪製地球軌跡
plt.plot(df_double['x_earth'], df_double['y_earth'], '-', color='blue', label='Earth (Double)')
# 如果您有 Posit 的結果，也可以繪製出來
# plt.plot(df_posit['x_earth'], df_posit['y_earth'], '-', color='red', label='Earth (Posit)')

# 繪製月亮軌跡
plt.plot(df_double['x_moon'], df_double['y_moon'], '-', color='gray', label='Moon (Double)')
# 如果您有 Posit 的結果，也可以繪製出來
# plt.plot(df_posit['x_moon'], df_posit['y_moon'], '-', color='purple', label='Moon (Posit)')


# --- 設定圖表屬性 ---
plt.title('Orbital Simulation (100 Years)')
plt.xlabel('X Position (AU)')
plt.ylabel('Y Position (AU)')
plt.axis('equal') # 確保 X 和 Y 軸比例相同，軌道不會變形
plt.grid(True)
plt.legend()

# 顯示圖表
plt.show()

# --- (可選) 製作動畫 ---
# 如果您想看到物體如何移動，可以製作一個動畫
# 這裡只展示製作動畫的基本結構，需要更詳細的設定
# fig_anim, ax_anim = plt.subplots(figsize=(10, 8))
# ax_anim.set_title('Orbital Simulation Animation')
# ax_anim.set_xlabel('X Position (AU)')
# ax_anim.set_ylabel('Y Position (AU)')
# ax_anim.axis('equal')
# ax_anim.grid(True)

# # 繪製太陽 (靜態)
# ax_anim.plot(df_double['x_sun'].iloc[0], df_double['y_sun'].iloc[0], 'o', markersize=5, color='yellow', label='Sun')

# # 繪製各天體的初始點
# line_earth, = ax_anim.plot([], [], '-', color='blue', label='Earth (Double)')
# line_moon, = ax_anim.plot([], [], '-', color='gray', label='Moon (Double)')
# # 如果有 Posit 的數據，也加入 line_earth_posit, line_moon_posit

# # 設定動畫更新函數
# def update(frame):
#     line_earth.set_data(df_double['x_earth'][:frame], df_double['y_earth'][:frame])
#     line_moon.set_data(df_double['x_moon'][:frame], df_double['y_moon'][:frame])
#     # 如果有 Posit 的數據，也更新
#     return line_earth, line_moon, # return all artists that were modified

# # 創建動畫
# # interval 是每幀之間的延遲時間 (毫秒)，frames 是總幀數
# ani = animation.FuncAnimation(fig_anim, update, frames=len(df_double), interval=10, blit=True)

# plt.legend()
# plt.show()