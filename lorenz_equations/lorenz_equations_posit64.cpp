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
    outputFile << "t,x,y,z" << std::endl;

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
        Posit64 dxdt, dydt, dzdt;
        lorenz_equations(x, y, z, dxdt, dydt, dzdt);
        
        // 歐拉法核心：同時更新 x, y, z 的值
        Posit64 nextX = x + h * dxdt;
        Posit64 nextY = y + h * dydt;
        Posit64 nextZ = z + h * dzdt;

        x = nextX;
        y = nextY;
        z = nextZ;

        // 將結果寫入檔案
        outputFile << (i + 1) * h << "," << x << "," << y << "," << z << std::endl;
        std::cout << std::fixed << std::setprecision(30) << (i + 1) * h << "," << x << "," << y << "," << z << std::endl;
    }

    // 關閉檔案
    outputFile.close();
    std::cout << "Data has been successfully written to lorenz_output_posit.csv" << std::endl;

    return 0;
}