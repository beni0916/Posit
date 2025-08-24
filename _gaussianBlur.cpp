#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#include <gmp.h>
#include <mpfr.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <filesystem>
#include <string>
#include <iomanip>
#include <bits/stdc++.h>
#include <fstream>
#include "testTool.h"

namespace fs = std::filesystem;
using namespace std;
using Posit32 = sw::universal::posit<32, 2>;
using Posit16_1 = sw::universal::posit<16, 1>;
using Posit16_2 = sw::universal::posit<16, 2>;

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
void applyGaussianBlurHorizontal(int width, int height, const unsigned char* inputData, double* outputData, const std::vector<double>& kernel, int channels) {
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
double* gaussianBlur(int width, int height, const unsigned char* rgbData, double sigma, int radius) {
    int channels = 3; // RGB 影像有 3 個通道
    std::vector<double> kernel = generateGaussianKernel(radius, sigma);

    double* tempHorizontal = new double[width * height * channels];
    // 直接在函式內部建立 tempVertical，並回傳
    double* tempVertical = new double[width * height * channels];

    // 先進行水平模糊
    applyGaussianBlurHorizontal(width, height, rgbData, tempHorizontal, kernel, channels);
    // 再進行垂直模糊
    applyGaussianBlurVertical(width, height, tempHorizontal, tempVertical, kernel, channels);

    delete[] tempHorizontal;
    // 不在此處轉換回 unsigned char，而是直接回傳浮點數結果
    return tempVertical;
}

int main() {
    std::string folderPath = "testImg";
    double sigma = 1.0; // 高斯模糊的標準差
    int radius = 2;     // 高斯核心的半徑 (例如，半徑為 2 會生成 5x5 的核心)

    // 確保 output 資料夾存在
    if (!fs::exists("output")) {
        fs::create_directory("output");
    }

    // 開啟 RMSE 結果輸出檔案
    std::ofstream rmseOutputFile("output/rmse_results.txt", std::ios::app);
    if (!rmseOutputFile.is_open()) {
        std::cerr << "無法開啟 RMSE 輸出檔案: output/rmse_results.txt" << std::endl;
        return 1; // 開啟檔案失敗，終止程式
    }
    
    // 開啟 Exponent 分佈結果輸出檔案
    std::ofstream exponentOutputFile("output/exponent_distribution.txt", std::ios::app);
    if (!exponentOutputFile.is_open()) {
        std::cerr << "無法開啟指數分佈輸出檔案: output/exponent_distribution.txt" << std::endl;
        return 1; // 開啟檔案失敗，終止程式
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
                
                int pixelCount = width * height * 3;

                std::cout << "處理影像: " << filename << ", 寬度: " << width << ", 高度: " << height << std::endl;
                rmseOutputFile << "處理影像: " << filename << ", 寬度: " << width << ", 高度: " << height << std::endl;
                exponentOutputFile << "處理影像: " << filename << ", 寬度: " << width << ", 高度: " << height << std::endl;

                // ***** 獲取原始浮點數結果來進行分析 *****
                double* originalFloatData = new double[pixelCount];
                for(int i = 0; i < pixelCount; ++i) {
                     originalFloatData[i] = static_cast<double>(rgbData[i]) / 1.0;
                }
                analyzeFloatExponentDistribution("原始影像", originalFloatData, pixelCount, exponentOutputFile);

                // ***** 獲取不同精度的高斯模糊浮點數結果 *****
                double* blurredFloatResult = gaussianBlur(width, height, rgbData, sigma, radius);
                
                Posit64* blurredPosit64Result = new Posit64[pixelCount];
                Posit32* blurredPosit32Result = new Posit32[pixelCount];
                Posit16_1* blurredPosit16_1Result = new Posit16_1[pixelCount];
                Posit16_2* blurredPosit16_2Result = new Posit16_2[pixelCount];
                mpfr_t* blurredMpfrResult = new mpfr_t[pixelCount]; // mpfr_t 需要特別處理

                for (int i = 0; i < pixelCount; i++) {
                    blurredPosit64Result[i] = Posit64(blurredFloatResult[i]);
                    blurredPosit32Result[i] = Posit32(blurredFloatResult[i]);
                    blurredPosit16_1Result[i] = Posit16_1(blurredFloatResult[i]);
                    blurredPosit16_2Result[i] = Posit16_2(blurredFloatResult[i]);
                    mpfr_init2(blurredMpfrResult[i], 256);
                    mpfr_set_d(blurredMpfrResult[i], blurredFloatResult[i], MPFR_RNDN);
                }

                // 將數值轉str並進行誤差運算 (在浮點數層級)
                vector<double> ieeeRMSEVals, pos64RMSEVals, pos32RMSEVals, pos16_1RMSEVals, pos16_2RMSEVals;
                vector<double> MPFRVals;
                // 新增用於計算整數 RMSE 的向量
                vector<double> ieeeIntRMSEVals, pos64IntRMSEVals, pos32IntRMSEVals, pos16_1IntRMSEVals, pos16_2IntRMSEVals;
                
                unsigned char* blurredRgbDataMpfrRef = new unsigned char[width * height * 3];
                std::string filenameOnly = fs::path(filename).filename().string();
                std::string filenameWithoutExt = getFilenameWithoutExtension(filenameOnly);
                // 創建一個新的文字檔案來儲存浮點數值
                std::string outputTextFilename = "output/" + filenameWithoutExt + "_fp_values.txt";
                std::ofstream outputFile(outputTextFilename);
                if (!outputFile.is_open()) {
                    std::cerr << "無法開啟浮點數輸出檔案: " << outputTextFilename << std::endl;
                }

                for (int i = 0; i < pixelCount; i++) {
                    // 將 MPFR 的結果轉換為 double 作為參考值
                    string mpfrStr = toString(blurredMpfrResult[i]);
                    MPFRVals.push_back(stod(mpfrStr));

                    // 將數值轉str
                    std::string ieeeStr = toString(blurredFloatResult[i]);
                    std::string posit64Str = toString(blurredPosit64Result[i]);
                    std::string posit32Str = toString(blurredPosit32Result[i]);
                    std::string posit16_1Str = toString(blurredPosit16_1Result[i]);
                    std::string posit16_2Str = toString(blurredPosit16_2Result[i]);

                    // 將結果寫入檔案
                    outputFile << "Pixel " << i << ":\n";
                    outputFile << "  MPFR: " << mpfrStr << "\n";
                    outputFile << "  IEEE: " << ieeeStr << "\n";
                    outputFile << "  Posit64: " << posit64Str << "\n";
                    outputFile << "  Posit32: " << posit32Str << "\n";
                    outputFile << "  Posit16_1: " << posit16_1Str << "\n";
                    outputFile << "  Posit16_2: " << posit16_2Str << "\n";

                    // 計算浮點數層級的誤差
                    double ieee754ResultDiff = stod(difference(mpfrStr, ieeeStr));
                    double posit64ResultDiff = stod(difference(mpfrStr, posit64Str));
                    double posit32ResultDiff = stod(difference(mpfrStr, posit32Str));
                    double posit16_1ResultDiff = stod(difference(mpfrStr, posit16_1Str));
                    double posit16_2ResultDiff = stod(difference(mpfrStr, posit16_2Str));
                    
                    outputFile << "  IEEE d: " << ieee754ResultDiff << "\n";
                    outputFile << "  Posit64 d: " << posit64ResultDiff << "\n";
                    outputFile << "  Posit32 d: " << posit32ResultDiff << "\n";
                    outputFile << "  Posit16_1 d: " << posit16_1ResultDiff << "\n";
                    outputFile << "  Posit16_2 d: " << posit16_2ResultDiff << "\n";
                    outputFile << "--------------------\n";

                    ieeeRMSEVals.push_back(ieee754ResultDiff);
                    pos64RMSEVals.push_back(posit64ResultDiff);
                    pos32RMSEVals.push_back(posit32ResultDiff);
                    pos16_1RMSEVals.push_back(posit16_1ResultDiff);
                    pos16_2RMSEVals.push_back(posit16_2ResultDiff);
                    
                    // 計算整數層級的誤差
                    blurredRgbDataMpfrRef[i] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, mpfr_get_d(blurredMpfrResult[i], MPFR_RNDN) * 1.0)));
                    unsigned char blurredRgbDataFloatVal = static_cast<unsigned char>(std::min(255.0, std::max(0.0, blurredFloatResult[i] * 1.0)));
                    unsigned char blurredRgbDataPosit64Val = static_cast<unsigned char>((int)Posit_floor(blurredPosit64Result[i] * Posit64(1.0)));
                    unsigned char blurredRgbDataPosit32Val = static_cast<unsigned char>((int)Posit_floor(blurredPosit32Result[i] * Posit32(1.0)));
                    unsigned char blurredRgbDataPosit16_1Val = static_cast<unsigned char>((int)Posit_floor(blurredPosit16_1Result[i] * Posit16_1(1.0)));
                    unsigned char blurredRgbDataPosit16_2Val = static_cast<unsigned char>((int)Posit_floor(blurredPosit16_2Result[i] * Posit16_2(1.0)));
                    
                    ieeeIntRMSEVals.push_back(static_cast<double>(blurredRgbDataMpfrRef[i] - blurredRgbDataFloatVal));
                    pos64IntRMSEVals.push_back(static_cast<double>(blurredRgbDataMpfrRef[i] - blurredRgbDataPosit64Val));
                    pos32IntRMSEVals.push_back(static_cast<double>(blurredRgbDataMpfrRef[i] - blurredRgbDataPosit32Val));
                    pos16_1IntRMSEVals.push_back(static_cast<double>(blurredRgbDataMpfrRef[i] - blurredRgbDataPosit16_1Val));
                    pos16_2IntRMSEVals.push_back(static_cast<double>(blurredRgbDataMpfrRef[i] - blurredRgbDataPosit16_2Val));
                }
                outputFile.close(); // 關閉浮點數結果檔案
                std::cout << "浮點數結果已寫入: " << outputTextFilename << std::endl;
                
                // ***** 輸出浮點數層級 RMSE *****
                double ieeeFpRMSE = RMSE(ieeeRMSEVals);
                double posit64FpRMSE = RMSE(pos64RMSEVals);
                double posit32FpRMSE = RMSE(pos32RMSEVals);
                double posit16_1FpRMSE = RMSE(pos16_1RMSEVals);
                double posit16_2FpRMSE = RMSE(pos16_2RMSEVals);

                // ***** 浮點數層級 平均相對誤差 (MRE) *****
                double ieeeFpMRE = calculateMeanRelativeError(ieeeRMSEVals, MPFRVals);
                double posit64FpMRE = calculateMeanRelativeError(pos64RMSEVals, MPFRVals);
                double posit32FpMRE = calculateMeanRelativeError(pos32RMSEVals, MPFRVals);
                double posit16_1FpMRE = calculateMeanRelativeError(pos16_1RMSEVals, MPFRVals);
                double posit16_2FpMRE = calculateMeanRelativeError(pos16_2RMSEVals, MPFRVals);

                // ***** 浮點數層級 中位數誤差 *****
                double ieeeFpMedian = calculateMedian(ieeeRMSEVals);
                double posit64FpMedian = calculateMedian(pos64RMSEVals);
                double posit32FpMedian = calculateMedian(pos32RMSEVals);
                double posit16_1FpMedian = calculateMedian(pos16_1RMSEVals);
                double posit16_2FpMedian = calculateMedian(pos16_2RMSEVals);

                // ***** 浮點數層級 標準差誤差 *****
                double ieeeFpStdDev = calculateStandardDeviation(ieeeRMSEVals);
                double posit64FpStdDev = calculateStandardDeviation(pos64RMSEVals);
                double posit32FpStdDev = calculateStandardDeviation(pos32RMSEVals);
                double posit16_1FpStdDev = calculateStandardDeviation(pos16_1RMSEVals);
                double posit16_2FpStdDev = calculateStandardDeviation(pos16_2RMSEVals);
                
                // --- 輸出中位數 ---
                std::cout << fixed  << setprecision(50) << "--- 浮點數層級 中位數 (vs MPFR) ---" << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "--- 浮點數層級 中位數 (vs MPFR) ---" << std::endl;
                std::cout << fixed  << setprecision(50) << "IEEE FpMedian:" << ieeeFpMedian << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "IEEE FpMedian:" << ieeeFpMedian << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT64 FpMedian:" << posit64FpMedian << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT64 FpMedian:" << posit64FpMedian << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT32 FpMedian:" << posit32FpMedian << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT32 FpMedian:" << posit32FpMedian << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_1 FpMedian:" << posit16_1FpMedian << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_1 FpMedian:" << posit16_1FpMedian << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_2 FpMedian:" << posit16_2FpMedian << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_2 FpMedian:" << posit16_2FpMedian << std::endl;

                // --- 輸出標準差 ---
                std::cout << fixed  << setprecision(50) << "--- 浮點數層級 標準差 (vs MPFR) ---" << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "--- 浮點數層級 標準差 (vs MPFR) ---" << std::endl;
                std::cout << fixed  << setprecision(50) << "IEEE FpStdDev:" << ieeeFpStdDev << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "IEEE FpStdDev:" << ieeeFpStdDev << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT64 FpStdDev:" << posit64FpStdDev << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT64 FpStdDev:" << posit64FpStdDev << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT32 FpStdDev:" << posit32FpStdDev << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT32 FpStdDev:" << posit32FpStdDev << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_1 FpStdDev:" << posit16_1FpStdDev << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_1 FpStdDev:" << posit16_1FpStdDev << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_2 FpStdDev:" << posit16_2FpStdDev << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_2 FpStdDev:" << posit16_2FpStdDev << std::endl;

                // --- 輸出平均相對誤差（MRE） ---
                std::cout << fixed  << setprecision(50) << "--- 浮點數層級 MRE (vs MPFR) ---" << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "--- 浮點數層級 MRE (vs MPFR) ---" << std::endl;
                std::cout << fixed  << setprecision(50) << "IEEE FpMRE:" << ieeeFpMRE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "IEEE FpMRE:" << ieeeFpMRE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT64 FpMRE:" << posit64FpMRE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT64 FpMRE:" << posit64FpMRE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT32 FpMRE:" << posit32FpMRE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT32 FpMRE:" << posit32FpMRE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_1 FpMRE:" << posit16_1FpMRE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_1 FpMRE:" << posit16_1FpMRE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_2 FpMRE:" << posit16_2FpMRE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_2 FpMRE:" << posit16_2FpMRE << std::endl;

                std::cout << fixed  << setprecision(50) << "--- 浮點數層級 RMSE (vs MPFR) ---" << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "--- 浮點數層級 RMSE (vs MPFR) ---" << std::endl;
                std::cout << fixed  << setprecision(50) << "IEEE FpRMSE:" << ieeeFpRMSE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "IEEE FpRMSE:" << ieeeFpRMSE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT64 FpRMSE:" << posit64FpRMSE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT64 FpRMSE:" << posit64FpRMSE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT32 FpRMSE:" << posit32FpRMSE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT32 FpRMSE:" << posit32FpRMSE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_1 FpRMSE:" << posit16_1FpRMSE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_1 FpRMSE:" << posit16_1FpRMSE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_2 FpRMSE:" << posit16_2FpRMSE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_2 FpRMSE:" << posit16_2FpRMSE << std::endl;
                
                // ***** 輸出整數層級 RMSE *****
                double ieeeIntRMSE = RMSE(ieeeIntRMSEVals);
                double posit64IntRMSE = RMSE(pos64IntRMSEVals);
                double posit32IntRMSE = RMSE(pos32IntRMSEVals);
                double posit16_1IntRMSE = RMSE(pos16_1IntRMSEVals);
                double posit16_2IntRMSE = RMSE(pos16_2IntRMSEVals);
                
                std::cout << fixed  << setprecision(50) << "\n--- 整數層級 RMSE (vs MPFR 轉換後) ---" << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "\n--- 整數層級 RMSE (vs MPFR 轉換後) ---" << std::endl;
                std::cout << fixed  << setprecision(50) << "IEEE intRMSE:" << ieeeIntRMSE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "IEEE intRMSE:" << ieeeIntRMSE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT64 intRMSE:" << posit64IntRMSE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT64 intRMSE:" << posit64IntRMSE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT32 intRMSE:" << posit32IntRMSE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT32 intRMSE:" << posit32IntRMSE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_1 intRMSE:" << posit16_1IntRMSE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_1 intRMSE:" << posit16_1IntRMSE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_2 intRMSE:" << posit16_2IntRMSE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_2 intRMSE:" << posit16_2IntRMSE << std::endl;

                rmseOutputFile << "----------------------------------------\n";

                // ***** 將浮點數結果轉換回 0-255 的 unsigned char 陣列，以便儲存為圖像 *****
                unsigned char* blurredRgbDataFloat = new unsigned char[pixelCount];
                unsigned char* blurredRgbDataPosit64 = new unsigned char[pixelCount];
                unsigned char* blurredRgbDataPosit32 = new unsigned char[pixelCount];
                unsigned char* blurredRgbDataPosit16_1 = new unsigned char[pixelCount];
                unsigned char* blurredRgbDataPosit16_2 = new unsigned char[pixelCount];
                unsigned char* blurredRgbDataMpfr = new unsigned char[pixelCount];

                for (int i = 0; i < pixelCount; ++i) {
                    blurredRgbDataFloat[i] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, blurredFloatResult[i] * 1.0)));
                    blurredRgbDataPosit64[i] = static_cast<unsigned char>((int)Posit_floor(blurredPosit64Result[i] * Posit64(1.0)));
                    blurredRgbDataPosit32[i] = static_cast<unsigned char>((int)Posit_floor(blurredPosit32Result[i] * Posit32(1.0)));
                    blurredRgbDataPosit16_1[i] = static_cast<unsigned char>((int)Posit_floor(blurredPosit16_1Result[i] * Posit16_1(1.0)));
                    blurredRgbDataPosit16_2[i] = static_cast<unsigned char>((int)Posit_floor(blurredPosit16_2Result[i] * Posit16_2(1.0)));
                    blurredRgbDataMpfr[i] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, mpfr_get_d(blurredMpfrResult[i], MPFR_RNDN) * 1.0)));
                }

                // 產生輸出影像檔案名稱並儲存
                std::string outputFilenameFloat = "output/" + filenameWithoutExt + "_blurred_float." + extension;
                stbi_write_png(outputFilenameFloat.c_str(), width, height, 3, blurredRgbDataFloat, width * 3);

                std::string outputFilenamePosit64 = "output/" + filenameWithoutExt + "_blurred_Posit64." + extension;
                stbi_write_png(outputFilenamePosit64.c_str(), width, height, 3, blurredRgbDataPosit64, width * 3);
                
                std::string outputFilenamePosit32 = "output/" + filenameWithoutExt + "_blurred_Posit32." + extension;
                stbi_write_png(outputFilenamePosit32.c_str(), width, height, 3, blurredRgbDataPosit32, width * 3);

                std::string outputFilenamePosit16_1 = "output/" + filenameWithoutExt + "_blurred_Posit16_1." + extension;
                stbi_write_png(outputFilenamePosit16_1.c_str(), width, height, 3, blurredRgbDataPosit16_1, width * 3);

                std::string outputFilenamePosit16_2 = "output/" + filenameWithoutExt + "_blurred_Posit16_2." + extension;
                stbi_write_png(outputFilenamePosit16_2.c_str(), width, height, 3, blurredRgbDataPosit16_2, width * 3);
                
                std::string outputFilenameMpfr = "output/" + filenameWithoutExt + "_blurred_Mpfr." + extension;
                stbi_write_png(outputFilenameMpfr.c_str(), width, height, 3, blurredRgbDataMpfr, width * 3);


                // ***** 記憶體釋放 *****
                stbi_image_free(rgbData);
                delete[] originalFloatData;
                delete[] blurredFloatResult;
                delete[] blurredPosit64Result;
                delete[] blurredPosit32Result;
                delete[] blurredPosit16_1Result;
                delete[] blurredPosit16_2Result;
                // MPFR 陣列需要逐個清除 mpfr_t 元素，然後再刪除陣列本身
                for (int i = 0; i < width * height * 3; ++i) {
                    mpfr_clear(blurredMpfrResult[i]);
                }
                delete[] blurredMpfrResult;
                delete[] blurredRgbDataFloat;
                delete[] blurredRgbDataPosit64;
                delete[] blurredRgbDataPosit32;
                delete[] blurredRgbDataPosit16_1;
                delete[] blurredRgbDataPosit16_2;
                delete[] blurredRgbDataMpfr;
                delete[] blurredRgbDataMpfrRef;

                std::cout << "浮點數高斯模糊影像已儲存為: " << outputFilenameFloat << std::endl;
                std::cout << "Posit64 高斯模糊影像已儲存為: " << outputFilenamePosit64 << std::endl;
                std::cout << "Posit32 高斯模糊影像已儲存為: " << outputFilenamePosit32 << std::endl;
                std::cout << "Posit16_1 高斯模糊影像已儲存為: " << outputFilenamePosit16_1 << std::endl;
                std::cout << "Posit16_2 高斯模糊影像已儲存為: " << outputFilenamePosit16_2 << std::endl;
                std::cout << "MPFR 高斯模糊影像已儲存為: " << outputFilenameMpfr << std::endl;
            }
        }
    }
    
    // 關閉 RMSE 檔案 
    rmseOutputFile.close();
    // 關閉 exponent 檔案
    exponentOutputFile.close();

    return 0;
}
