#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

// 定義我們要積分的函數
// 這裡使用 C++11 的函式指標 (function pointer) 或 lambda 表達式，
// 以保持函式的彈性。
double gaussianFunction(double x) {
    return exp(-x * x);
}

// 使用梯形法進行數值積分的函式
double numericalIntegrationTrapezoidal(double (*func)(double), double a, double b, int n) {
    std::ofstream outputFile("double.csv");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return 1;
    }
    // 每個小梯形的寬度
    double deltaX = (b - a) / n;

    // 根據梯形法公式計算
    // 兩端點只佔 0.5，中間點佔 1.0
    double integralResult = 0.5 * (func(a) + func(b));

    outputFile << std::fixed << std::setprecision(30);
    outputFile << "steps,integralResult" << std::endl;

    outputFile << 0 << "," << integralResult << std::endl;
    std::cout << std::fixed << std::setprecision(30)
              << 0 << "," << integralResult << std::endl;

    for (int i = 1; i < n; ++i) {
        // 計算每個子區間的 x 值
        double x = a + i * deltaX;
        integralResult += func(x);

        outputFile << i << "," << integralResult*deltaX << std::endl;
        std::cout << std::fixed << std::setprecision(30)
              << i << "," << integralResult*deltaX << std::endl;
    }

    integralResult *= deltaX;

    outputFile.close();
    return integralResult;
}

int main() {
    // 設定積分參數
    const double lowerBound = 0.0;
    const double upperBound = 1.0;
    // 增加 n 可以提高精度
    const int numberOfTrapezoids = 1000000;

    // 執行數值積分
    double result = numericalIntegrationTrapezoidal(gaussianFunction, lowerBound, upperBound, numberOfTrapezoids);

    // 設定輸出精度
    std::cout << std::fixed << std::setprecision(15);
    
    std::cout << "高斯函數 f(x) = e^(-x^2) 在區間 [" << lowerBound << ", " << upperBound << "] 上" << std::endl;
    std::cout << "使用 " << numberOfTrapezoids << " 個梯形進行數值積分的結果為: " << result << std::endl;

    return 0;
}