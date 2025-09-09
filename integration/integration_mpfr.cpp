#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <mpfr.h>

// 定義高精度下要積分的函數 f(x) = exp(-x^2)
void gaussianFunction(mpfr_t result, const mpfr_t x) {
    mpfr_t temp;
    mpfr_init2(temp, mpfr_get_prec(x));
    
    mpfr_sqr(temp, x, MPFR_RNDN);
    mpfr_neg(temp, temp, MPFR_RNDN);
    
    mpfr_exp(result, temp, MPFR_RNDN);

    mpfr_clear(temp);
}

// 使用梯形法進行高精度數值積分
void numericalIntegrationTrapezoidal_MPFR(void (*func)(mpfr_t, const mpfr_t), double a, double b, int n, const char* filename) {
    const mpfr_prec_t precision = 256;
    
    mpfr_t lowerBound, upperBound, integralSum, deltaX, x, temp, zero_point_five;

    mpfr_init2(lowerBound, precision);
    mpfr_init2(upperBound, precision);
    mpfr_init2(integralSum, precision);
    mpfr_init2(deltaX, precision);
    mpfr_init2(x, precision);
    mpfr_init2(temp, precision);
    mpfr_init2(zero_point_five, precision);

    mpfr_set_d(lowerBound, a, MPFR_RNDN);
    mpfr_set_d(upperBound, b, MPFR_RNDN);
    mpfr_set_d(zero_point_five, 0.5, MPFR_RNDN);

    mpfr_sub(deltaX, upperBound, lowerBound, MPFR_RNDN);
    mpfr_div_si(deltaX, deltaX, n, MPFR_RNDN);

    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return;
    }

    outputFile << "steps,integralResult_mpfr" << std::endl;

    // 先計算兩端點的函數值總和
    func(temp, lowerBound);
    mpfr_mul(integralSum, temp, zero_point_five, MPFR_RNDN);
    
    func(temp, upperBound);
    mpfr_mul(temp, temp, zero_point_five, MPFR_RNDN);
    mpfr_add(integralSum, integralSum, temp, MPFR_RNDN);

    // 輸出初始值
    mpfr_mul(temp, integralSum, deltaX, MPFR_RNDN);
    // mpfr_out_str 函式：
    // NULL: 使用預設的輸出緩衝區
    // 0: 輸出所有有效位數
    // 10: 十進位
    // MPFR_RNDN: 捨入模式
    char* result_str_init = mpfr_get_str(NULL, NULL, 10, 30, temp, MPFR_RNDN);
    outputFile << 0 << "," << result_str_init << std::endl;
    mpfr_free_str(result_str_init);

    for (int i = 1; i < n; ++i) {
        // 計算每個子區間的 x 值
        mpfr_set_d(temp, i, MPFR_RNDN);
        mpfr_mul(x, temp, deltaX, MPFR_RNDN);
        mpfr_add(x, x, lowerBound, MPFR_RNDN);

        // 將中間點的函數值加到總和中
        func(temp, x);
        mpfr_add(integralSum, integralSum, temp, MPFR_RNDN);
        
        // 輸出當前總和乘以 deltaX 的結果
        mpfr_mul(temp, integralSum, deltaX, MPFR_RNDN);
        char* result_str = mpfr_get_str(NULL, NULL, 10, 30, temp, MPFR_RNDN);
        outputFile << i << "," << result_str << std::endl;
        mpfr_free_str(result_str);
    }
    
    std::cout << "MPFR 256位元高斯函數在 [0, 1] 上的積分結果請參閱檔案 " << filename << std::endl;
    
    mpfr_clear(lowerBound);
    mpfr_clear(upperBound);
    mpfr_clear(integralSum);
    mpfr_clear(deltaX);
    mpfr_clear(x);
    mpfr_clear(temp);
    mpfr_clear(zero_point_five);

    outputFile.close();
}

int main() {
    const double lowerBound_d = 0.0;
    const double upperBound_d = 1.0;
    const int numberOfTrapezoids = 1000000;

    numericalIntegrationTrapezoidal_MPFR(gaussianFunction, lowerBound_d, upperBound_d, numberOfTrapezoids, "mpfr.csv");

    return 0;
}