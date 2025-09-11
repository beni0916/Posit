#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "../Posit/myfdlibm.h"

// 定義一個表示2D向量的結構
struct Vector2D {
    Posit64 x;
    Posit64 y;
};

// 重載向量運算子，以簡化程式碼
Vector2D operator+(const Vector2D& a, const Vector2D& b) {
    return {a.x + b.x, a.y + b.y};
}

Vector2D operator-(const Vector2D& a, const Vector2D& b) {
    return {a.x - b.x, a.y - b.y};
}

Vector2D operator*(Posit64 scalar, const Vector2D& vec) {
    return {scalar * vec.x, scalar * vec.y};
}

Vector2D operator/(const Vector2D& vec, Posit64 scalar) {
    return {vec.x / scalar, vec.y / scalar};
}

// 三體問題中，計算每個物體所受的總加速度
std::vector<Vector2D> calculateAccelerations(
    const std::vector<Vector2D>& positions,
    const std::vector<Posit64>& masses) {
    
    std::vector<Vector2D> accelerations(3, {0.0, 0.0});
    const Posit64 G_AU = Posit64(39.47841760435743); // (AU^3 / yr^2 / M_sun)

    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            Vector2D r = positions[j] - positions[i];
            Posit64 dist2 = r.x * r.x + r.y * r.y;
            Posit64 dist3 = dist2 * Posit_sqrt(dist2);
            Vector2D force = (G_AU / dist3) * r;

            accelerations[i] = accelerations[i] + masses[j] * force;
            accelerations[j] = accelerations[j] - masses[i] * force;
        }
    }
    return accelerations;
}

// RK4 數值積分方法
void runSimulation_RK4_3Body(Posit64 dt, int steps, const char* filename) {
    // 初始質量 (太陽、地球、月亮)
    std::vector<Posit64> masses = {
        Posit64(1.0), // 太陽質量
        Posit64(3.00348959632e-6), // 地球質量
        Posit64(3.69460599661e-8) // 月亮質量
    };

    // 初始位置
    std::vector<Vector2D> positions = {
        {Posit64(0.0), Posit64(0.0)}, // 太陽
        {Posit64(1.0), Posit64(0.0)}, // 地球
        {Posit64(1.0 + 0.00257), Posit64(0.0)} // 月亮
    };

    // 初始速度
    std::vector<Vector2D> velocities = {
        {Posit64(0.0), Posit64(0.0)}, // 太陽
        {Posit64(0.0), Posit64(6.2831853)}, // 地球
        {Posit64(0.0), Posit64(6.2831853) + Posit64(0.2159)} // 月亮
    };

    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return;
    }
    outputFile << "steps,x_sun,y_sun,x_earth,y_earth,x_moon,y_moon" << std::endl;
    outputFile << std::fixed << std::setprecision(30) << 0 << "," << positions[0].x << "," << positions[0].y << ","
               << positions[1].x << "," << positions[1].y << ","
               << positions[2].x << "," << positions[2].y << std::endl;

    for (int i = 0; i < steps; ++i) {
        // RK4 核心邏輯
        std::vector<Vector2D> k1_vel = velocities;
        std::vector<Vector2D> k1_acc = calculateAccelerations(positions, masses);

        std::vector<Vector2D> k2_pos_temp(3), k2_vel_temp(3);
        for(int j = 0; j < 3; ++j) {
            k2_pos_temp[j] = positions[j] + (dt / 2.0) * k1_vel[j];
            k2_vel_temp[j] = velocities[j] + (dt / 2.0) * k1_acc[j];
        }
        std::vector<Vector2D> k2_acc = calculateAccelerations(k2_pos_temp, masses);

        std::vector<Vector2D> k3_pos_temp(3), k3_vel_temp(3);
        for(int j = 0; j < 3; ++j) {
            k3_pos_temp[j] = positions[j] + (dt / 2.0) * k2_vel_temp[j];
            k3_vel_temp[j] = velocities[j] + (dt / 2.0) * k2_acc[j];
        }
        std::vector<Vector2D> k3_acc = calculateAccelerations(k3_pos_temp, masses);

        std::vector<Vector2D> k4_pos_temp(3), k4_vel_temp(3);
        for(int j = 0; j < 3; ++j) {
            k4_pos_temp[j] = positions[j] + dt * k3_vel_temp[j];
            k4_vel_temp[j] = velocities[j] + dt * k3_acc[j];
        }
        std::vector<Vector2D> k4_acc = calculateAccelerations(k4_pos_temp, masses);

        for (int j = 0; j < 3; ++j) {
            velocities[j] = velocities[j] + (dt / 6.0) * (k1_acc[j] + (2.0 * k2_acc[j]) + (2.0 * k3_acc[j]) + k4_acc[j]);
            positions[j] = positions[j] + (dt / 6.0) * (k1_vel[j] + (2.0 * k2_vel_temp[j]) + (2.0 * k3_vel_temp[j]) + k4_vel_temp[j]);
        }
        
        outputFile << std::fixed << std::setprecision(30) << i + 1 << "," << positions[0].x << "," << positions[0].y << ","
                   << positions[1].x << "," << positions[1].y << ","
                   << positions[2].x << "," << positions[2].y << std::endl;
    }
    outputFile.close();
}

int main() {
    const Posit64 TIME_STEP_AU = 0.01; // 時間步長 (單位：yr)
    const int SIM_STEPS = 100000; // 模擬步數 (約10年)

    runSimulation_RK4_3Body(TIME_STEP_AU, SIM_STEPS, "Posit64.csv");

    std::cout << "三體問題模擬結果已儲存至 Posit64.csv 檔案中。" << std::endl;
    return 0;
}