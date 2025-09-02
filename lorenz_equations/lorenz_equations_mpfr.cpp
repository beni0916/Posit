#include <iostream>
#include <fstream>
#include <iomanip>
#include <mpfr.h>

// 洛倫茲方程組的參數
const int PRECISION = 256; // 運算精度 (bits)
mpfr_t sigma, rho, beta;

// 洛倫茲方程組的函數
void lorenzEquations(mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t dxdt, mpfr_t dydt, mpfr_t dzdt) {
    mpfr_t temp1, temp2, temp3;
    mpfr_inits2(PRECISION, temp1, temp2, temp3, (mpfr_t *)0);

    // dxdt = sigma * (y - x)
    mpfr_sub(temp1, y, x, MPFR_RNDN);         // temp1 = y - x
    mpfr_mul(dxdt, sigma, temp1, MPFR_RNDN);  // dxdt = sigma * temp1

    // dydt = x * (rho - z) - y
    mpfr_sub(temp2, rho, z, MPFR_RNDN);       // temp2 = rho - z
    mpfr_mul(temp3, x, temp2, MPFR_RNDN);     // temp3 = x * temp2
    mpfr_sub(dydt, temp3, y, MPFR_RNDN);      // dydt = temp3 - y

    // dzdt = x * y - beta * z
    mpfr_mul(temp1, x, y, MPFR_RNDN);         // temp1 = x * y
    mpfr_mul(temp2, beta, z, MPFR_RNDN);      // temp2 = beta * z
    mpfr_sub(dzdt, temp1, temp2, MPFR_RNDN);  // dzdt = temp1 - temp2

    mpfr_clears(temp1, temp2, temp3, (mpfr_t *)0);
}

int main() {
    // 初始化參數
    mpfr_inits2(PRECISION, sigma, rho, beta, (mpfr_t *)0);
    mpfr_set_d(sigma, 10.0, MPFR_RNDN);
    mpfr_set_d(rho, 28.0, MPFR_RNDN);
    mpfr_set_d(beta, 8.0 / 3.0, MPFR_RNDN);

    // 輸出檔案
    std::ofstream outputFile("lorenz_output_mpfr.csv");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return 1;
    }

    // 寫入 CSV 檔案的標頭
    outputFile << "t,x,y,z" << std::endl;

    // 初始條件和步長
    mpfr_t x, y, z, h, t;
    mpfr_inits2(PRECISION, x, y, z, h, t, (mpfr_t *)0);
    mpfr_set_d(x, 0.0, MPFR_RNDN);
    mpfr_set_d(y, 1.0, MPFR_RNDN);
    mpfr_set_d(z, 1.0, MPFR_RNDN);
    mpfr_set_d(h, 0.01, MPFR_RNDN);  // 時間步長
    mpfr_set_d(t, 0.0, MPFR_RNDN);   // 初始化時間

    int steps = 10000; // 迭代次數
    char buffer[512];

    // 寫入初始值
    mpfr_snprintf(buffer, sizeof(buffer), "%.30Rf", t);
    outputFile << buffer << ",";
    mpfr_snprintf(buffer, sizeof(buffer), "%.30Rf", x);
    outputFile << buffer << ",";
    mpfr_snprintf(buffer, sizeof(buffer), "%.30Rf", y);
    outputFile << buffer << ",";
    mpfr_snprintf(buffer, sizeof(buffer), "%.30Rf", z);
    outputFile << buffer << std::endl;
    std::cout << buffer << std::endl;

    for (int i = 0; i < steps; ++i) {
        // 計算當前時間點的變化率
        mpfr_t dxdt, dydt, dzdt;
        mpfr_inits2(PRECISION, dxdt, dydt, dzdt, (mpfr_t *)0);
        lorenzEquations(x, y, z, dxdt, dydt, dzdt);
        
        // 歐拉法更新
        mpfr_t nextX, nextY, nextZ;
        mpfr_inits2(PRECISION, nextX, nextY, nextZ, (mpfr_t *)0);

        mpfr_mul(nextX, h, dxdt, MPFR_RNDN);
        mpfr_add(nextX, x, nextX, MPFR_RNDN);

        mpfr_mul(nextY, h, dydt, MPFR_RNDN);
        mpfr_add(nextY, y, nextY, MPFR_RNDN);

        mpfr_mul(nextZ, h, dzdt, MPFR_RNDN);
        mpfr_add(nextZ, z, nextZ, MPFR_RNDN);

        mpfr_set(x, nextX, MPFR_RNDN);
        mpfr_set(y, nextY, MPFR_RNDN);
        mpfr_set(z, nextZ, MPFR_RNDN);
        mpfr_add(t, t, h, MPFR_RNDN);

        // 輸出到 CSV
        mpfr_snprintf(buffer, sizeof(buffer), "%.30Rf", t);
        outputFile << buffer << ",";
        mpfr_snprintf(buffer, sizeof(buffer), "%.30Rf", x);
        outputFile << buffer << ",";
        mpfr_snprintf(buffer, sizeof(buffer), "%.30Rf", y);
        outputFile << buffer << ",";
        mpfr_snprintf(buffer, sizeof(buffer), "%.30Rf", z);
        outputFile << buffer << std::endl;
        std::cout << buffer << std::endl;

        mpfr_clears(dxdt, dydt, dzdt, nextX, nextY, nextZ, (mpfr_t *)0);
    }

    // 清理 MPFR 變數
    mpfr_clears(sigma, rho, beta, x, y, z, h, t, (mpfr_t *)0);

    // 關閉檔案
    outputFile.close();
    std::cout << "Data has been successfully written to lorenz_output_mpfr.csv" << std::endl;
    std::cout << "Compile with: g++ lorenz.cpp -o lorenz -lmpfr -lgmp" << std::endl;

    return 0;
}
