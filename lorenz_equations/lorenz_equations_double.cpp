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
    std::ofstream outputFile("lorenz_output_double.csv");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return 1;
    }

    outputFile << std::fixed << std::setprecision(30);
    outputFile << "steps,x,y,z" << std::endl;

    double x = 0.0, y = 1.0, z = 1.0;
    double h = 0.01;
    int steps = 10000;

    outputFile << 0 << "," << x << "," << y << "," << z << std::endl;
    std::cout << std::fixed << std::setprecision(30)
              << 0 << "," << x << "," << y << "," << z << std::endl;

    for (int i = 0; i < steps; i++) {
        double k1x, k1y, k1z;
        double k2x, k2y, k2z;
        double k3x, k3y, k3z;
        double k4x, k4y, k4z;

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

        outputFile << i+1 << "," << x << "," << y << "," << z << std::endl;
        std::cout << std::fixed << std::setprecision(30) << i+1 << "," << x << "," << y << "," << z << std::endl;
    }

    outputFile.close();
    std::cout << "Data has been successfully written to lorenz_output_double.csv" << std::endl;
    return 0;
}
