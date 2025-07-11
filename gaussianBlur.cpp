#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#include "Posit/myfdlibm.h" // 假設你的 Posit 函式庫在這裡
#include <gmp.h>
#include <mpfr.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <filesystem>
#include <string>
#include <iomanip> // for std::setprecision
#include <bits/stdc++.h>
#include <fstream>

namespace fs = std::filesystem;
using namespace std;

// 函式：生成一維高斯核心
std::vector<double> generateGaussianKernel(int radius, float sigma) {
    int kernelSize = 2 * radius + 1;
    std::vector<double> kernel(kernelSize);
    double sum = 0.0f;

    for (int i = 0; i < kernelSize; ++i) {
        double x = i - radius;
        kernel[i] = std::exp(-(x * x) / (2 * sigma * sigma));
        sum += kernel[i];
    }

    // 歸一化核心
    for (int i = 0; i < kernelSize; ++i) {
        kernel[i] /= sum;
    }
    return kernel;
}

// 函式：應用高斯模糊 (水平方向)
void applyGaussianBlurHorizontal(int width, int height, const double* inputData, double* outputData, const std::vector<double>& kernel, int channels) {
    int radius = (kernel.size() - 1) / 2;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            for (int c = 0; c < channels; ++c) {
                double sum = 0.0f;
                for (int k = -radius; k <= radius; ++k) {
                    int pixelX = x + k;
                    if (pixelX < 0) pixelX = 0; // 邊界處理：複製邊界像素
                    if (pixelX >= width) pixelX = width - 1; // 邊界處理：複製邊界像素

                    sum += inputData[(y * width + pixelX) * channels + c] * kernel[k + radius];
                }
                outputData[(y * width + x) * channels + c] = sum;
            }
        }
    }
}

// 函式：應用高斯模糊 (垂直方向)
void applyGaussianBlurVertical(int width, int height, const double* inputData, double* outputData, const std::vector<double>& kernel, int channels) {
    int radius = (kernel.size() - 1) / 2;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            for (int c = 0; c < channels; ++c) {
                double sum = 0.0f;
                for (int k = -radius; k <= radius; ++k) {
                    int pixelY = y + k;
                    if (pixelY < 0) pixelY = 0; // 邊界處理：複製邊界像素
                    if (pixelY >= height) pixelY = height - 1; // 邊界處理：複製邊界像素

                    sum += inputData[(pixelY * width + x) * channels + c] * kernel[k + radius];
                }
                outputData[(y * width + x) * channels + c] = sum;
            }
        }
    }
}

// 高斯模糊主函式 (Float 版本) - 回傳浮點數結果
double* gaussianBlurFloatStb(int width, int height, const unsigned char* rgbData, double sigma, int radius) {
    int channels = 3; // RGB 影像有 3 個通道
    std::vector<double> kernel = generateGaussianKernel(radius, sigma);

    double* normalizedRgbData = new double[width * height * channels];
    for (int i = 0; i < width * height * channels; ++i) {
        normalizedRgbData[i] = static_cast<double>(rgbData[i]) / 255.0f;
    }

    double* tempHorizontal = new double[width * height * channels];
    // 直接在函式內部建立 tempVertical，並回傳
    double* tempVertical = new double[width * height * channels];

    // 先進行水平模糊
    applyGaussianBlurHorizontal(width, height, normalizedRgbData, tempHorizontal, kernel, channels);
    // 再進行垂直模糊
    applyGaussianBlurVertical(width, height, tempHorizontal, tempVertical, kernel, channels);

    delete[] normalizedRgbData;
    delete[] tempHorizontal;
    // 不在此處轉換回 unsigned char，而是直接回傳浮點數結果
    return tempVertical;
}

// 函式：生成一維高斯核心
std::vector<Posit64> generateGaussianKernelPosit(int radius, float sigma) {
    int kernelSize = 2 * radius + 1;
    std::vector<Posit64> kernel(kernelSize);
    Posit64 sum = Posit64(0.0);

    for (int i = 0; i < kernelSize; ++i) {
        Posit64 x = i - radius;
        kernel[i] = Posit_exp(-(x * x) / (2 * Posit64(sigma) * Posit64(sigma)));
        sum += kernel[i];
    }

    // 歸一化核心
    for (int i = 0; i < kernelSize; ++i) {
        kernel[i] /= sum;
    }
    return kernel;
}

// 函式：應用高斯模糊 (水平方向)
void applyGaussianBlurHorizontalPosit(int width, int height, const unsigned char* inputData, Posit64* outputData, const std::vector<Posit64>& kernel, int channels) {
    int radius = (kernel.size() - 1) / 2;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            for (int c = 0; c < channels; ++c) {
                Posit64 sum = Posit64(0.0);
                for (int k = -radius; k <= radius; ++k) {
                    int pixelX = x + k;
                    if (pixelX < 0) pixelX = 0; // 邊界處理：複製邊界像素
                    if (pixelX >= width) pixelX = width - 1; // 邊界處理：複製邊界像素

                    sum += (Posit64(inputData[(y * width + pixelX) * channels + c]) / Posit64(255.0)) * kernel[k + radius];
                }
                outputData[(y * width + x) * channels + c] = sum;
            }
        }
    }
}

// 函式：應用高斯模糊 (垂直方向)
void applyGaussianBlurVerticalPosit(int width, int height, const Posit64* inputData, Posit64* outputData, const std::vector<Posit64>& kernel, int channels) {
    int radius = (kernel.size() - 1) / 2;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            for (int c = 0; c < channels; ++c) {
                Posit64 sum = Posit64(0.0);
                for (int k = -radius; k <= radius; ++k) {
                    int pixelY = y + k;
                    if (pixelY < 0) pixelY = 0; // 邊界處理：複製邊界像素
                    if (pixelY >= height) pixelY = height - 1; // 邊界處理：複製邊界像素

                    sum += inputData[(pixelY * width + x) * channels + c] * kernel[k + radius];
                }
                outputData[(y * width + x) * channels + c] = sum;
            }
        }
    }
}

// 高斯模糊主函式 (Posit64 版本) - 回傳 Posit64 陣列結果
Posit64* gaussianBlurPosit64Stb(int width, int height, const unsigned char* rgbData, double sigma, int radius) {
    int channels = 3; // RGB 影像有 3 個通道
    std::vector<Posit64> kernel = generateGaussianKernelPosit(radius, sigma);

    Posit64* tempHorizontal = new Posit64[width * height * channels];
    // 直接在函式內部建立 tempVertical，並回傳
    Posit64* tempVertical = new Posit64[width * height * channels];

    // 先進行水平模糊
    applyGaussianBlurHorizontalPosit(width, height, rgbData, tempHorizontal, kernel, channels);
    // 再進行垂直模糊
    applyGaussianBlurVerticalPosit(width, height, tempHorizontal, tempVertical, kernel, channels);

    delete[] tempHorizontal;
    // 不在此處轉換回 unsigned char，而是直接回傳 Posit64 結果
    return tempVertical;
}

// 函式：生成一維高斯核心 (MPFR 版本)
std::vector<mpfr_t> generateGaussianKernelMpfr(int radius, mpfr_t sigma, mpfr_prec_t prec) {
    int kernelSize = 2 * radius + 1;
    std::vector<mpfr_t> kernel(kernelSize); // 儲存 mpfr_t 的 vector

    mpfr_t sum;
    mpfr_init2(sum, prec); // 初始化 sum
    mpfr_set_d(sum, 0.0, MPFR_RNDN); // sum 歸零

    mpfr_t x_val, term, sigma_sq, two_sigma_sq;
    mpfr_init2(x_val, prec);
    mpfr_init2(term, prec);
    mpfr_init2(sigma_sq, prec);
    mpfr_init2(two_sigma_sq, prec);

    mpfr_mul(sigma_sq, sigma, sigma, MPFR_RNDN); // sigma * sigma
    mpfr_mul_d(two_sigma_sq, sigma_sq, 2.0, MPFR_RNDN); // 2 * sigma * sigma

    for (int i = 0; i < kernelSize; ++i) {
        mpfr_init2(kernel[i], prec); // 初始化 vector 中的每個 mpfr_t
        mpfr_set_d(x_val, static_cast<double>(i - radius), MPFR_RNDN); // x = i - radius

        // 計算 -(x * x) / (2 * sigma * sigma)
        mpfr_mul(term, x_val, x_val, MPFR_RNDN); // x * x
        mpfr_neg(term, term, MPFR_RNDN); // -(x * x)
        mpfr_div(term, term, two_sigma_sq, MPFR_RNDN); // -(x * x) / (2 * sigma * sigma)

        mpfr_exp(kernel[i], term, MPFR_RNDN); // kernel[i] = exp(...)
        mpfr_add(sum, sum, kernel[i], MPFR_RNDN); // sum += kernel[i]
    }

    // 歸一化核心
    for (int i = 0; i < kernelSize; ++i) {
        mpfr_div(kernel[i], kernel[i], sum, MPFR_RNDN); // kernel[i] /= sum
    }

    // 清除臨時變數
    mpfr_clear(sum);
    mpfr_clear(x_val);
    mpfr_clear(term);
    mpfr_clear(sigma_sq);
    mpfr_clear(two_sigma_sq);

    return kernel;
}

// 函式：應用高斯模糊 (水平方向) (MPFR 版本)
void applyGaussianBlurHorizontalMpfr(int width, int height, const mpfr_t* inputData, mpfr_t* outputData, const std::vector<mpfr_t>& kernel, int channels, mpfr_prec_t prec) {
    int radius = (kernel.size() - 1) / 2;

    mpfr_t sum_pixel, temp_val; // 臨時變數
    mpfr_init2(sum_pixel, prec);
    mpfr_init2(temp_val, prec);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            for (int c = 0; c < channels; ++c) {
                mpfr_set_d(sum_pixel, 0.0, MPFR_RNDN); // 每個新像素的 sum 歸零

                for (int k = -radius; k <= radius; ++k) {
                    int pixelX = x + k;
                    // 邊界處理：複製邊界像素
                    if (pixelX < 0) pixelX = 0;
                    if (pixelX >= width) pixelX = width - 1;

                    // sum += inputData[...] * kernel[...]
                    // 注意 inputData 是已經正規化過的 mpfr_t 陣列
                    mpfr_mul(temp_val, inputData[(y * width + pixelX) * channels + c], kernel[k + radius], MPFR_RNDN);
                    mpfr_add(sum_pixel, sum_pixel, temp_val, MPFR_RNDN);
                }
                mpfr_set(outputData[(y * width + x) * channels + c], sum_pixel, MPFR_RNDN); // 將結果設置到 outputData
            }
        }
    }
    mpfr_clear(sum_pixel);
    mpfr_clear(temp_val);
}

// 函式：應用高斯模糊 (垂直方向) (MPFR 版本)
void applyGaussianBlurVerticalMpfr(int width, int height, const mpfr_t* inputData, mpfr_t* outputData, const std::vector<mpfr_t>& kernel, int channels, mpfr_prec_t prec) {
    int radius = (kernel.size() - 1) / 2;

    mpfr_t sum_pixel, temp_val; // 臨時變數
    mpfr_init2(sum_pixel, prec);
    mpfr_init2(temp_val, prec);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            for (int c = 0; c < channels; ++c) {
                mpfr_set_d(sum_pixel, 0.0, MPFR_RNDN); // 每個新像素的 sum 歸零

                for (int k = -radius; k <= radius; ++k) {
                    int pixelY = y + k;
                    // 邊界處理：複製邊界像素
                    if (pixelY < 0) pixelY = 0;
                    if (pixelY >= height) pixelY = height - 1;

                    // sum += inputData[...] * kernel[...]
                    mpfr_mul(temp_val, inputData[(pixelY * width + x) * channels + c], kernel[k + radius], MPFR_RNDN);
                    mpfr_add(sum_pixel, sum_pixel, temp_val, MPFR_RNDN);
                }
                mpfr_set(outputData[(y * width + x) * channels + c], sum_pixel, MPFR_RNDN); // 將結果設置到 outputData
            }
        }
    }
    mpfr_clear(sum_pixel);
    mpfr_clear(temp_val);
}

// 高斯模糊主函式 (MPFR 版本) - 回傳 MPFR 陣列結果
mpfr_t* gaussianBlurMpfrStb(int width, int height, const unsigned char* rgbData, double sigma_double, int radius) {
    mpfr_prec_t prec = 256;
    int channels = 3; // RGB 影像有 3 個通道

    // 1. 將輸入的 unsigned char 影像轉換為 MPFR 格式並正規化 (0-1)
    mpfr_t* inputMpfrData = new mpfr_t[width * height * channels];
    mpfr_t div255;
    mpfr_init2(div255, prec);
    mpfr_set_d(div255, 255.0, MPFR_RNDN);

    for (int i = 0; i < width * height * channels; ++i) {
        mpfr_init2(inputMpfrData[i], prec);
        mpfr_set_d(inputMpfrData[i], static_cast<double>(rgbData[i]), MPFR_RNDN);
        mpfr_div(inputMpfrData[i], inputMpfrData[i], div255, MPFR_RNDN); // 正規化
    }
    mpfr_clear(div255);

    // 準備 sigma 的 MPFR 版本
    mpfr_t sigma_mpfr;
    mpfr_init2(sigma_mpfr, prec);
    mpfr_set_d(sigma_mpfr, sigma_double, MPFR_RNDN);

    // 2. 生成高斯核心
    std::vector<mpfr_t> kernel = generateGaussianKernelMpfr(radius, sigma_mpfr, prec);

    // 3. 準備中間結果的 MPFR 陣列
    mpfr_t* tempHorizontalMpfr = new mpfr_t[width * height * channels];
    // 直接在函式內部建立 tempVerticalMpfr，並回傳
    mpfr_t* tempVerticalMpfr = new mpfr_t[width * height * channels];

    // 初始化中間結果的 mpfr_t 元素
    for (int i = 0; i < width * height * channels; ++i) {
        mpfr_init2(tempHorizontalMpfr[i], prec);
        mpfr_init2(tempVerticalMpfr[i], prec);
    }

    // 4. 先進行水平模糊
    applyGaussianBlurHorizontalMpfr(width, height, inputMpfrData, tempHorizontalMpfr, kernel, channels, prec);
    // 5. 再進行垂直模糊
    applyGaussianBlurVerticalMpfr(width, height, tempHorizontalMpfr, tempVerticalMpfr, kernel, channels, prec);

    // 7. 清除所有 MPFR 變數和陣列
    mpfr_clear(sigma_mpfr);

    for (int i = 0; i < width * height * channels; ++i) {
        mpfr_clear(inputMpfrData[i]);
        mpfr_clear(tempHorizontalMpfr[i]);
        // tempVerticalMpfr 將被回傳，在呼叫者處清除
    }
    delete[] inputMpfrData;
    delete[] tempHorizontalMpfr;

    // 清除 kernel 中的 mpfr_t
    for (size_t i = 0; i < kernel.size(); ++i) {
        mpfr_clear(kernel[i]);
    }
    
    // 回傳 MPFR 浮點數結果
    return tempVerticalMpfr;
}


// 取得檔案名稱（不含副檔名）的函數
std::string getFilenameWithoutExtension(const std::string& filename) {
    size_t dotPos = filename.find_last_of('.');
    if (dotPos != std::string::npos) {
        return filename.substr(0, dotPos);
    }
    return filename;
}

// 取得檔案副檔名的函數 (轉為小寫)
std::string getFileExtension(const std::string& filename) {
    size_t dotPos = filename.find_last_of('.');
    if (dotPos != std::string::npos) {
        std::string extension = filename.substr(dotPos + 1);
        for (char& c : extension) {
            c = std::tolower(c);
        }
        return extension;
    }
    return "";
}

double RMSE(const std::vector<double> &vec) {
    if (vec.empty()) {
        return 0.0; // 或者拋出例外，視你的需求而定
    }

    double sumOfSquares = std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0);
    double meanSquareError = sumOfSquares / vec.size();
    return std::sqrt(meanSquareError);
}

std::string Difference(std::string& num1, std::string& num2) {
    std::string result = "";
    
    int n1 = num1.length();
    int n2 = num2.length();
    if(n2 != n1) return "-1.0"; // 如果長度不同，表示有問題，返回錯誤值

    std::string str1 = num1;
    std::string str2 = num2;
    // 確保 str1 >= str2，方便減法
    if(str1.compare(str2) < 0) str1.swap(str2); // 使用 compare 而不是 < 比較整個字串

    // 處理小數點，如果有的話
    size_t dotPos1 = str1.find('.');
    size_t dotPos2 = str2.find('.');

    // 如果只有一個有小數點或小數點位置不同，需要進一步處理確保對齊
    if (dotPos1 != std::string::npos || dotPos2 != std::string::npos) {
        // 找到最長的小數部分
        size_t maxLenAfterDot = 0;
        if (dotPos1 != std::string::npos) maxLenAfterDot = std::max(maxLenAfterDot, str1.length() - 1 - dotPos1);
        if (dotPos2 != std::string::npos) maxLenAfterDot = std::max(maxLenAfterDot, str2.length() - 1 - dotPos2);

        // 補齊小數點後的零
        if (dotPos1 == std::string::npos) str1 += ".";
        if (dotPos2 == std::string::npos) str2 += ".";

        while (str1.length() - 1 - str1.find('.') < maxLenAfterDot) str1 += '0';
        while (str2.length() - 1 - str2.find('.') < maxLenAfterDot) str2 += '0';
    }

    // 補齊整數部分的零 (如果長度不一致)
    size_t maxLenBeforeDot = 0;
    if (dotPos1 != std::string::npos) maxLenBeforeDot = std::max(maxLenBeforeDot, dotPos1);
    else maxLenBeforeDot = std::max(maxLenBeforeDot, str1.length());
    
    if (dotPos2 != std::string::npos) maxLenBeforeDot = std::max(maxLenBeforeDot, dotPos2);
    else maxLenBeforeDot = std::max(maxLenBeforeDot, str2.length());

    std::string paddedStr1 = std::string(maxLenBeforeDot - (dotPos1 == std::string::npos ? str1.length() : dotPos1), '0') + str1;
    std::string paddedStr2 = std::string(maxLenBeforeDot - (dotPos2 == std::string::npos ? str2.length() : dotPos2), '0') + str2;

    // 現在確保字串總長度相同，並處理負數情況 (儘管我們已經保證 str1 >= str2)
    // 這裡只是為了一般化，但對於 RMSE 的絕對差值，可以直接用 str1 - str2
    // 再次檢查長度，確保減法時索引不會越界
    n1 = paddedStr1.length();
    n2 = paddedStr2.length();
    if(n1 != n2) { /* 這是個內部錯誤，表示補零邏輯有問題 */ return "-2.0"; }

    reverse(paddedStr1.begin(), paddedStr1.end());
    reverse(paddedStr2.begin(), paddedStr2.end());
    
    int carry = 0;
    for (int i = 0; i < n1; i++) { // 現在 n1 == n2
        if(paddedStr1[i] == '.'){
            result.push_back('.');
            continue;
        }

        int sub = ((paddedStr1[i] - '0') - (paddedStr2[i] - '0') - carry);

        if (sub < 0) {
            sub += 10;
            carry = 1;
        } else {
            carry = 0;
        }
        result.push_back(sub + '0');
    }

    // 移除前導零 (如果差值是小數，不移除小數點前的零)
    reverse(result.begin(), result.end());

    size_t first_digit = result.find_first_not_of('0');
    size_t dot_pos_result = result.find('.');

    if (first_digit == std::string::npos) { // All zeros
        return "0.0";
    } else if (first_digit > dot_pos_result && dot_pos_result != std::string::npos) { // e.g., 00.123
        return "0" + result.substr(dot_pos_result);
    } else if (first_digit == dot_pos_result && dot_pos_result != std::string::npos) { // e.g., .123 -> 0.123
        return "0" + result.substr(first_digit);
    } else {
        return result.substr(first_digit);
    }
}

std::string formatPositiveFloatString(const std::string& inputStr) {
    std::string result = inputStr;
    size_t dotPos = result.find('.');

    // --- 處理小數點前部分 ---
    std::string integerPart;
    if (dotPos == std::string::npos) {
        integerPart = result; // 如果沒有小數點，整個字串就是整數部分
    } else {
        integerPart = result.substr(0, dotPos); // 取得小數點前的部分
    }

    // 移除前導零，除非它就是 "0" (例如 "007" 變成 "7", "0" 保持 "0")
    size_t firstDigit = integerPart.find_first_not_of('0');
    if (firstDigit != std::string::npos) {
        integerPart = integerPart.substr(firstDigit);
    } else if (!integerPart.empty()) { // 處理全是零的情況，例如 "000"
        integerPart = "0";
    } else { // 處理空字串作為整數部分的情況 (例如輸入是 ".123" 或空字串)
        integerPart = "0";
    }

    // 補足小數點前到四位
    if (integerPart.length() < 4) {
        std::string padding(4 - integerPart.length(), '0');
        integerPart = padding + integerPart;
    }
    // 如果 integerPart.length() > 4，保持原樣，不截斷

    // --- 處理小數點後部分 ---
    std::string decimalPart;
    if (dotPos != std::string::npos) {
        decimalPart = result.substr(dotPos + 1); // 取得小數點後的部分
    }
    // 如果 dotPos == std::string::npos，decimalPart 保持為空字串

    // 補足或截斷小數點後到五十位
    if (decimalPart.length() < 50) {
        decimalPart.append(50 - decimalPart.length(), '0');
    } else if (decimalPart.length() > 50) {
        decimalPart.erase(50); // 從索引 50 開始移除，保留前 50 個字元
    }

    // --- 組合最終結果 ---
    return integerPart + "." + decimalPart;
}

// 更新 toString(Posit64) 以處理其字串表示，確保小數點後有五十位
string toString(Posit64 input) {
    ostringstream out;
    // 使用 long double 轉換以確保標準庫的格式化行為
    // 這移除了手動補零的邏輯，讓 ostringstream 和 setprecision 完全控制格式
    out << fixed << setprecision(50) << input; 

    std::string result = out.str();

    return formatPositiveFloatString(result);
}

// 更新 toString(float) 以處理其字串表示，確保小數點後有五十位
string toString(double input) {
    string numStr;
    ostringstream out;
    out << fixed << setprecision(50) << input; // 使用更高的精度
    numStr = out.str();

    // 根據需求，移除末尾零的邏輯將被移除
    // 'fixed' 和 'setprecision(50)' 會自動確保小數點後有 50 位，並補零

    return formatPositiveFloatString(numStr);
}


// MPFR 版本，修改為確保小數點後有五十位
string toString(mpfr_t input){
    // 使用足夠大的 buffer
    char buffer[200]; // 根據需要的精度調整大小，例如 50 位數字 + 小數點 + 符號 + null 終止符
    // mpfr_sprintf 使用 "%.50Rf" 會確保小數點後有 50 位，並補零
    mpfr_sprintf(buffer, "%.50Rf", input); 
    string num1(buffer);

    // 根據需求，移除末尾零的邏輯將被移除

    return formatPositiveFloatString(num1);
}

int main() {
    std::string folderPath = "testImg";
    double sigma = 1.0; // 高斯模糊的標準差
    int radius = 2;     // 高斯核心的半徑 (例如，半徑為 2 會生成 5x5 的核心)

    // 確保 output 資料夾存在
    if (!fs::exists("output")) {
        fs::create_directory("output");
    }

    for (const auto& entry : fs::directory_iterator(folderPath)) {
        if (fs::is_regular_file(entry)) {
            std::string filename = entry.path().string();
            std::string extension = getFileExtension(filename);

            if (extension == "jpg" || extension == "jpeg" || extension == "png" || extension == "bmp" || extension == "tif" || extension == "tiff") {
                int width, height, channels;
                unsigned char* rgbData = stbi_load(filename.c_str(), &width, &height, &channels, 3); // 強制讀取為 3 通道 (RGB)
                if (!rgbData) {
                    std::cerr << "無法讀取影像: " << filename << std::endl;
                    continue; // 讀取失敗，繼續處理下一個檔案
                }

                std::cout << "處理影像: " << filename << ", 寬度: " << width << ", 高度: " << height << std::endl;

                // ***** 獲取不同精度的高斯模糊浮點數結果 *****
                double* blurredFloatResult = gaussianBlurFloatStb(width, height, rgbData, sigma, radius);
                Posit64* blurredPosit64Result = gaussianBlurPosit64Stb(width, height, rgbData, sigma, radius);
                mpfr_t* blurredMpfrResult = gaussianBlurMpfrStb(width, height, rgbData, sigma, radius);

                // 將數值轉str並進行誤差運算 (在浮點數層級)
                vector<double> IEEE_RMSE_vals, POS_RMSE_vals;
                std::string filenameOnly = fs::path(filename).filename().string(); // 只取得檔案名稱
                std::string filenameWithoutExt = getFilenameWithoutExtension(filenameOnly); // 取得不含副檔名的檔案名稱
                // 創建一個新的文字檔案來儲存浮點數值
                std::string outputTextFilename = "output/" + filenameWithoutExt + "_fp_values.txt";
                std::ofstream outputFile(outputTextFilename);
                if (!outputFile.is_open()) {
                    std::cerr << "無法開啟浮點數輸出檔案: " << outputTextFilename << std::endl;
                }

                for (int i = 0; i < width * height * 3; i++) {
                    std::string ieeeStr = toString(blurredFloatResult[i]);
                    std::string positStr = toString(blurredPosit64Result[i]);
                    std::string mpfrStr = toString(blurredMpfrResult[i]);

                    // 將結果寫入檔案
                    outputFile << "Pixel " << i << ":\n";
                    outputFile << "  MPFR: " << mpfrStr << "\n";
                    outputFile << "  IEEE: " << ieeeStr << "\n";
                    outputFile << "  Posit: " << positStr << "\n";

                    double Posit_result_diff = stod(Difference(mpfrStr, positStr));
                    double IEEE754_result_diff = stod(Difference(mpfrStr, ieeeStr));

                    outputFile << "  IEEE d: " << IEEE754_result_diff << "\n";
                    outputFile << "  Posit d: " << Posit_result_diff << "\n";
                    outputFile << "--------------------\n";

                    IEEE_RMSE_vals.push_back(IEEE754_result_diff);
                    POS_RMSE_vals.push_back(Posit_result_diff);
                }
                outputFile.close(); // 關閉浮點數結果檔案
                std::cout << "浮點數結果已寫入: " << outputTextFilename << std::endl;

                // ***** 在將浮點數結果轉換回 0-255 之前，先計算 RMSE *****
                std::cout << fixed  << setprecision(50) << "IEEE RMSE (double vs MPFR):" << RMSE(IEEE_RMSE_vals) << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT RMSE (Posit64 vs MPFR):" << RMSE(POS_RMSE_vals) << std::endl;

                // ***** 將浮點數結果轉換回 0-255 的 unsigned char 陣列，以便儲存為圖像 *****
                unsigned char* blurredRgbDataFloat = new unsigned char[width * height * 3];
                unsigned char* blurredRgbDataPosit64 = new unsigned char[width * height * 3];
                unsigned char* blurredRgbDataMpfr = new unsigned char[width * height * 3]; // 用於儲存 MPFR 結果量化後的圖像

                for (int i = 0; i < width * height * 3; ++i) {
                    // 將 float 結果轉換回 unsigned char (0-255)
                    blurredRgbDataFloat[i] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, blurredFloatResult[i] * 255.0)));
                    // 將 Posit64 結果轉換回 unsigned char (0-255)
                    blurredRgbDataPosit64[i] = static_cast<unsigned char>((int)Posit_floor(blurredPosit64Result[i] * Posit64(255.0)));
                    // 將 MPFR 結果轉換回 unsigned char (0-255)
                    blurredRgbDataMpfr[i] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, mpfr_get_d(blurredMpfrResult[i], MPFR_RNDN) * 255.0)));
                }

                // 產生輸出影像檔案名稱並儲存
                std::string outputFilenameFloat = "output/" + filenameWithoutExt + "_blurred_float." + extension;
                stbi_write_png(outputFilenameFloat.c_str(), width, height, 3, blurredRgbDataFloat, width * 3);

                std::string outputFilenamePosit64 = "output/" + filenameWithoutExt + "_blurred_Posit64." + extension;
                stbi_write_png(outputFilenamePosit64.c_str(), width, height, 3, blurredRgbDataPosit64, width * 3);
                
                std::string outputFilenameMpfr = "output/" + filenameWithoutExt + "_blurred_Mpfr." + extension;
                stbi_write_png(outputFilenameMpfr.c_str(), width, height, 3, blurredRgbDataMpfr, width * 3);


                // ***** 記憶體釋放 *****
                stbi_image_free(rgbData);
                delete[] blurredFloatResult;
                delete[] blurredPosit64Result;
                // MPFR 陣列需要逐個清除 mpfr_t 元素，然後再刪除陣列本身
                for (int i = 0; i < width * height * 3; ++i) {
                    mpfr_clear(blurredMpfrResult[i]);
                }
                delete[] blurredMpfrResult;
                delete[] blurredRgbDataFloat;
                delete[] blurredRgbDataPosit64;
                delete[] blurredRgbDataMpfr;

                std::cout << "浮點數高斯模糊影像已儲存為: " << outputFilenameFloat << std::endl;
                std::cout << "Posit64 高斯模糊影像已儲存為: " << outputFilenamePosit64 << std::endl;
                std::cout << "MPFR 高斯模糊影像已儲存為: " << outputFilenameMpfr << std::endl;
            }
        }
    }

    return 0;
}