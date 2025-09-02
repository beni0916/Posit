#include <iostream>
#include <fstream>
#include <iomanip>
#include "../Posit/myfdlibm.h"

// 洛倫茲方程組的參數
const Posit64 sigma = Posit64(10.0);
const Posit64 rho = Posit64(28.0);
const Posit64 beta = Posit64(8.0) / Posit64(3.0);

// 定義常微分方程組的函數
void lorenz_equations(Posit64 x, Posit64 y, Posit64 z, Posit64& dxdt, Posit64& dydt, Posit64& dzdt) {
    dxdt = sigma * (y - x);
    dydt = x * (rho - z) - y;
    dzdt = x * y - beta * z;
}

int main() {
    // 輸出檔案
    std::ofstream outputFile("lorenz_output_posit.csv");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return 1;
    }

    // 設置輸出格式：使用科學記號且精確到 15 位小數
    outputFile << std::fixed << std::setprecision(30);
    
    // 寫入 CSV 檔案的標頭
    outputFile << "steps,x,y,z" << std::endl;

    // 初始條件和步長
    Posit64 x = Posit64(0.0);
    Posit64 y = Posit64(1.0);
    Posit64 z = Posit64(1.0); 
    Posit64 h = Posit64(0.01); // 時間步長
    int steps = 10000; // 迭代次數，這裡增加到 10000 以獲得更完整的軌跡

    // 寫入初始值
    outputFile << 0 << "," << x << "," << y << "," << z << std::endl;
    std::cout << std::fixed << std::setprecision(30) << 0 << "," << x << "," << y << "," << z << std::endl;

    for (int i = 0; i < steps; ++i) {
        // 計算當前時間點的變化率
        Posit64 k1x, k1y, k1z;
        Posit64 k2x, k2y, k2z;
        Posit64 k3x, k3y, k3z;
        Posit64 k4x, k4y, k4z;

        // k1
        lorenz_equations(x, y, z, k1x, k1y, k1z);

        // k2
        lorenz_equations(x + 0.5 * h * k1x, y + 0.5 * h * k1y, z + 0.5 * h * k1z,
                         k2x, k2y, k2z);

        // k3
        lorenz_equations(x + 0.5 * h * k2x, y + 0.5 * h * k2y, z + 0.5 * h * k2z,
                         k3x, k3y, k3z);

        // k4
        lorenz_equations(x + h * k3x, y + h * k3y, z + h * k3z,
                         k4x, k4y, k4z);

        // 更新
        x += (h / 6.0) * (k1x + 2*k2x + 2*k3x + k4x);
        y += (h / 6.0) * (k1y + 2*k2y + 2*k3y + k4y);
        z += (h / 6.0) * (k1z + 2*k2z + 2*k3z + k4z);

        // 將結果寫入檔案
        outputFile << i+1 << "," << x << "," << y << "," << z << std::endl;
        std::cout << std::fixed << std::setprecision(30) << i+1 << "," << x << "," << y << "," << z << std::endl;
    }

    // 關閉檔案
    outputFile.close();
    std::cout << "Data has been successfully written to lorenz_output_posit.csv" << std::endl;

    return 0;
}