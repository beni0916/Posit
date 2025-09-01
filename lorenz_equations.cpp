#include <iostream>
#include <fstream>
#include <iomanip>

// 洛倫茲方程組的參數
const double sigma = 10.0;
const double rho = 28.0;
const double beta = 8.0 / 3.0;

// 定義常微分方程組的函數
void lorenz_equations(double x, double y, double z, double& dxdt, double& dydt, double& dzdt) {
    dxdt = sigma * (y - x);
    dydt = x * (rho - z) - y;
    dzdt = x * y - beta * z;
}

int main() {
    // 輸出檔案
    std::ofstream outputFile("lorenz_output.csv");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return 1;
    }

    // 設置輸出格式：使用科學記號且精確到 15 位小數
    outputFile << std::fixed << std::setprecision(30);
    
    // 寫入 CSV 檔案的標頭
    outputFile << "t,x,y,z" << std::endl;

    // 初始條件和步長
    double x = 0.0;
    double y = 1.0;
    double z = 1.0; 
    double h = 0.01; // 時間步長
    int steps = 10000; // 迭代次數，這裡增加到 10000 以獲得更完整的軌跡

    // 寫入初始值
    outputFile << 0 << "," << x << "," << y << "," << z << std::endl;
    std::cout << std::fixed << std::setprecision(30) << 0 << "," << x << "," << y << "," << z << std::endl;

    for (int i = 0; i < steps; ++i) {
        // 計算當前時間點的變化率
        double dxdt, dydt, dzdt;
        lorenz_equations(x, y, z, dxdt, dydt, dzdt);
        
        // 歐拉法核心：同時更新 x, y, z 的值
        double nextX = x + h * dxdt;
        double nextY = y + h * dydt;
        double nextZ = z + h * dzdt;

        x = nextX;
        y = nextY;
        z = nextZ;

        // 將結果寫入檔案
        outputFile << (i + 1) * h << "," << x << "," << y << "," << z << std::endl;
        std::cout << std::fixed << std::setprecision(30) << (i + 1) * h << "," << x << "," << y << "," << z << std::endl;
    }

    // 關閉檔案
    outputFile.close();
    std::cout << "Data has been successfully written to lorenz_output.csv" << std::endl;

    return 0;
}