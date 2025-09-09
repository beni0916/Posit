#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

// 定義一個表示2D向量的結構
struct Vector2D {
    double x;
    double y;
};

const double G_AU = 39.478440530965371; // 萬有引力常數 (單位：AU^3 yr^-2 M_sun^-1)
const double SOLAR_MASS = 1.0; // 太陽質量 (單位：M_sun)

// 重載向量運算子，以簡化程式碼
Vector2D operator+(const Vector2D& a, const Vector2D& b) {
    return {a.x + b.x, a.y + b.y};
}

Vector2D operator*(double scalar, const Vector2D& vec) {
    return {scalar * vec.x, scalar * vec.y};
}

Vector2D operator/(const Vector2D& vec, double scalar) {
    return {vec.x / scalar, vec.y / scalar};
}

// 根據當前位置計算加速度 (a = F/m)
Vector2D calculateAcceleration(const Vector2D& pos) {
    
    double distance = sqrt(pos.x * pos.x + pos.y * pos.y);
    double accMagnitude = (G_AU * SOLAR_MASS) / (distance * distance);

    Vector2D acc = {
        -accMagnitude * pos.x / distance,
        -accMagnitude * pos.y / distance
    };
    return acc;
}

// RK4 數值積分方法
void runSimulation_RK4(double dt, int steps, const char* filename) {
    Vector2D pos = {1.0, 0.0}; // 初始位置 (1 AU, x軸)
    double iniV = sqrt(G_AU);
    Vector2D vel = {0.0, iniV}; // 初始速度 (單位：AU yr^-1)

    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return;
    }
    outputFile << "steps,x,y" << std::endl;
    outputFile << 0 << "," << pos.x << "," << pos.y << std::endl;

    for (int i = 0; i < steps; ++i) {
        // k1
        Vector2D k1_vel = vel;
        Vector2D k1_acc = calculateAcceleration(pos);

        // k2
        Vector2D k2_pos_temp = pos + (dt / 2.0) * k1_vel;
        Vector2D k2_vel_temp = vel + (dt / 2.0) * k1_acc;
        Vector2D k2_acc = calculateAcceleration(k2_pos_temp);

        // k3
        Vector2D k3_pos_temp = pos + (dt / 2.0) * k2_vel_temp;
        Vector2D k3_vel_temp = vel + (dt / 2.0) * k2_acc;
        Vector2D k3_acc = calculateAcceleration(k3_pos_temp);
        
        // k4
        Vector2D k4_pos_temp = pos + dt * k3_vel_temp;
        Vector2D k4_vel_temp = vel + dt * k3_acc;
        Vector2D k4_acc = calculateAcceleration(k4_pos_temp);
        
        // 更新速度和位置
        vel = vel + (dt / 6.0) * (k1_acc + (2.0 * k2_acc) + (2.0 * k3_acc) + k4_acc);
        pos = pos + (dt / 6.0) * (k1_vel + (2.0 * k2_vel_temp) + (2.0 * k3_vel_temp) + k4_vel_temp);
        
        // 將結果寫入檔案
        outputFile << std::fixed << std::setprecision(30) << i + 1 << "," << pos.x << "," << pos.y << std::endl;
    }
    outputFile.close();
}

// 主程式
int main() {
    const double TIME_STEP_AU = 0.001; // 時間步長 (單位：yr)
    const int SIM_STEPS = 100000; // 模擬步數 (約10年)

    // 執行 RK4 模擬
    runSimulation_RK4(TIME_STEP_AU, SIM_STEPS, "double.csv");

    std::cout << "RK4 軌道模擬結果已儲存至 orbit_simulation_rk4.csv 檔案中。" << std::endl;
    std::cout << "請使用繪圖軟體（如 Python, Excel）視覺化結果。" << std::endl;

    return 0;
}