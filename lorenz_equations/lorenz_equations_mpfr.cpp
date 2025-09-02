#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <mpfr.h>

// 固定使用最近似偶(MPFR_RNDN)
static const mpfr_rnd_t RND = MPFR_RNDN;

// 將 mpfr_t 以固定小數位輸出成 std::string（使用 MPFR 的 printf 轉字串）
std::string mpfrToString(mpfr_srcptr v, int fracDigits = 30) {
    // 先取得需要的長度（不含終止字元）
    int n = mpfr_snprintf(nullptr, 0, "%.*Rf", fracDigits, v);
    std::vector<char> buf(n + 1);
    mpfr_snprintf(buf.data(), buf.size(), "%.*Rf", fracDigits, v);
    return std::string(buf.data());
}

// 洛倫茲方程組：dxdt, dydt, dzdt
// dxdt = sigma * (y - x)
// dydt = x * (rho - z) - y
// dzdt = x * y - beta * z
void lorenzEquations(mpfr_srcptr x, mpfr_srcptr y, mpfr_srcptr z,
                     mpfr_ptr dxdt, mpfr_ptr dydt, mpfr_ptr dzdt,
                     mpfr_srcptr sigma, mpfr_srcptr rho, mpfr_srcptr beta) {
    mpfr_t t1, t2;
    mpfr_init2(t1, mpfr_get_prec(dxdt));
    mpfr_init2(t2, mpfr_get_prec(dxdt));

    // dxdt = sigma * (y - x)
    mpfr_sub(t1, y, x, RND);
    mpfr_mul(dxdt, sigma, t1, RND);

    // dydt = x * (rho - z) - y
    mpfr_sub(t1, rho, z, RND);
    mpfr_mul(t2, x, t1, RND);
    mpfr_sub(dydt, t2, y, RND);

    // dzdt = x * y - beta * z
    mpfr_mul(t1, x, y, RND);
    mpfr_mul(t2, beta, z, RND);
    mpfr_sub(dzdt, t1, t2, RND);

    mpfr_clear(t1);
    mpfr_clear(t2);
}

int main() {
    // ====== 參數設定 ======
    const mpfr_prec_t PREC = 256; // MPFR 256-bit precision
    const int printDigits = 30;   // 固定小數 30 位
    const int steps = 10000;      // 積分步數
    // 步長、初值、參數等都用十進位字串設定，避免二進位轉換誤差
    const char* hStr = "0.01";
    const char* x0Str = "0.0";
    const char* y0Str = "1.0";
    const char* z0Str = "1.0";
    const char* sigmaStr = "10.0";
    const char* rhoStr   = "28.0";
    const char* betaNumStr = "8.0";
    const char* betaDenStr = "3.0";

    // ====== 檔案輸出 ======
    std::ofstream outputFile("lorenz_output_mpfr.csv");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return 1;
    }
    outputFile << std::fixed; // 我們自己控制小數位數字串

    // ====== 初始化 MPFR 變數 ======
    mpfr_set_default_prec(PREC);
    mpfr_set_default_rounding_mode(RND);

    mpfr_t x, y, z, h;
    mpfr_inits2(PREC, x, y, z, h, (mpfr_ptr) 0);

    // 參數
    mpfr_t sigma, rho, beta, betaNum, betaDen;
    mpfr_inits2(PREC, sigma, rho, beta, betaNum, betaDen, (mpfr_ptr) 0);

    // 常數
    mpfr_t two, half, six, hDiv6, halfH;
    mpfr_inits2(PREC, two, half, six, hDiv6, halfH, (mpfr_ptr) 0);

    mpfr_set_str(h, hStr, 10, RND);
    mpfr_set_str(x, x0Str, 10, RND);
    mpfr_set_str(y, y0Str, 10, RND);
    mpfr_set_str(z, z0Str, 10, RND);

    mpfr_set_str(sigma, sigmaStr, 10, RND);
    mpfr_set_str(rho, rhoStr, 10, RND);
    mpfr_set_str(betaNum, betaNumStr, 10, RND);
    mpfr_set_str(betaDen, betaDenStr, 10, RND);
    mpfr_div(beta, betaNum, betaDen, RND); // beta = 8/3

    mpfr_set_ui(two, 2u, RND);
    mpfr_set_str(half, "0.5", 10, RND);
    mpfr_set_ui(six, 6u, RND);

    // 預先計算 h/6 與 0.5*h
    mpfr_div(hDiv6, h, six, RND);
    mpfr_mul(halfH, h, half, RND);

    // RK4 暫存量
    mpfr_t k1x, k1y, k1z, k2x, k2y, k2z, k3x, k3y, k3z, k4x, k4y, k4z;
    mpfr_inits2(PREC, k1x, k1y, k1z, k2x, k2y, k2z, k3x, k3y, k3z, k4x, k4y, k4z, (mpfr_ptr) 0);

    // 暫存 for 中間點與更新
    mpfr_t xt, yt, zt;               // x + α*h*ki*
    mpfr_t tmp1, tmp2, tmp3, tmp4;   // 各式暫存
    mpfr_inits2(PREC, xt, yt, zt, tmp1, tmp2, tmp3, tmp4, (mpfr_ptr) 0);

    // ====== 輸出表頭與初值 ======
    outputFile << "steps,x,y,z\n";
    std::cout << std::fixed;

    auto printLine = [&](int step) {
        std::string xs = mpfrToString(x, printDigits);
        std::string ys = mpfrToString(y, printDigits);
        std::string zs = mpfrToString(z, printDigits);
        outputFile << step << "," << xs << "," << ys << "," << zs << "\n";
        std::cout    << step << "," << xs << "," << ys << "," << zs << "\n";
    };

    printLine(0);

    // ====== 主迴圈：RK4 ======
    for (int i = 0; i < steps; ++i) {
        // k1 = f(x, y, z)
        lorenzEquations(x, y, z, k1x, k1y, k1z, sigma, rho, beta);

        // (x, y, z) + 0.5*h*k1
        mpfr_mul(tmp1, halfH, k1x, RND);
        mpfr_add(xt, x, tmp1, RND);
        mpfr_mul(tmp1, halfH, k1y, RND);
        mpfr_add(yt, y, tmp1, RND);
        mpfr_mul(tmp1, halfH, k1z, RND);
        mpfr_add(zt, z, tmp1, RND);

        // k2
        lorenzEquations(xt, yt, zt, k2x, k2y, k2z, sigma, rho, beta);

        // (x, y, z) + 0.5*h*k2
        mpfr_mul(tmp1, halfH, k2x, RND);
        mpfr_add(xt, x, tmp1, RND);
        mpfr_mul(tmp1, halfH, k2y, RND);
        mpfr_add(yt, y, tmp1, RND);
        mpfr_mul(tmp1, halfH, k2z, RND);
        mpfr_add(zt, z, tmp1, RND);

        // k3
        lorenzEquations(xt, yt, zt, k3x, k3y, k3z, sigma, rho, beta);

        // (x, y, z) + h*k3
        mpfr_mul(tmp1, h, k3x, RND);
        mpfr_add(xt, x, tmp1, RND);
        mpfr_mul(tmp1, h, k3y, RND);
        mpfr_add(yt, y, tmp1, RND);
        mpfr_mul(tmp1, h, k3z, RND);
        mpfr_add(zt, z, tmp1, RND);

        // k4
        lorenzEquations(xt, yt, zt, k4x, k4y, k4z, sigma, rho, beta);

        // x += (h/6) * (k1 + 2*k2 + 2*k3 + k4)
        // tmp1 = 2*k2x, tmp2 = 2*k3x, tmp3 = k1x + tmp1, tmp3 += tmp2, tmp3 += k4x
        mpfr_mul(tmp1, two, k2x, RND);
        mpfr_mul(tmp2, two, k3x, RND);
        mpfr_add(tmp3, k1x, tmp1, RND);
        mpfr_add(tmp3, tmp3, tmp2, RND);
        mpfr_add(tmp3, tmp3, k4x, RND);
        mpfr_mul(tmp4, hDiv6, tmp3, RND);
        mpfr_add(x, x, tmp4, RND);

        // y
        mpfr_mul(tmp1, two, k2y, RND);
        mpfr_mul(tmp2, two, k3y, RND);
        mpfr_add(tmp3, k1y, tmp1, RND);
        mpfr_add(tmp3, tmp3, tmp2, RND);
        mpfr_add(tmp3, tmp3, k4y, RND);
        mpfr_mul(tmp4, hDiv6, tmp3, RND);
        mpfr_add(y, y, tmp4, RND);

        // z
        mpfr_mul(tmp1, two, k2z, RND);
        mpfr_mul(tmp2, two, k3z, RND);
        mpfr_add(tmp3, k1z, tmp1, RND);
        mpfr_add(tmp3, tmp3, tmp2, RND);
        mpfr_add(tmp3, tmp3, k4z, RND);
        mpfr_mul(tmp4, hDiv6, tmp3, RND);
        mpfr_add(z, z, tmp4, RND);

        printLine(i + 1);
    }

    outputFile.close();
    std::cout << "Data has been successfully written to lorenz_output_mpfr256.csv\n";

    // 釋放資源
    mpfr_clears(x, y, z, h, sigma, rho, beta, betaNum, betaDen,
                two, half, six, hDiv6, halfH,
                k1x, k1y, k1z, k2x, k2y, k2z, k3x, k3y, k3z, k4x, k4y, k4z,
                xt, yt, zt, tmp1, tmp2, tmp3, tmp4,
                (mpfr_ptr) 0);
    mpfr_free_cache();

    return 0;
}
