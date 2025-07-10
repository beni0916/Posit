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

namespace fs = std::filesystem;
using namespace std;

// 函式：生成一維高斯核心
std::vector<float> generateGaussianKernel(int radius, float sigma) {
    int kernelSize = 2 * radius + 1;
    std::vector<float> kernel(kernelSize);
    float sum = 0.0f;

    for (int i = 0; i < kernelSize; ++i) {
        float x = i - radius;
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
void applyGaussianBlurHorizontal(int width, int height, const unsigned char* inputData, float* outputData, const std::vector<float>& kernel, int channels) {
    int radius = (kernel.size() - 1) / 2;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            for (int c = 0; c < channels; ++c) {
                float sum = 0.0f;
                for (int k = -radius; k <= radius; ++k) {
                    int pixelX = x + k;
                    if (pixelX < 0) pixelX = 0; // 邊界處理：複製邊界像素
                    if (pixelX >= width) pixelX = width - 1; // 邊界處理：複製邊界像素

                    sum += (static_cast<float>(inputData[(y * width + pixelX) * channels + c]) / 255.0f) * kernel[k + radius];
                }
                outputData[(y * width + x) * channels + c] = sum;
            }
        }
    }
}

// 函式：應用高斯模糊 (垂直方向)
void applyGaussianBlurVertical(int width, int height, const float* inputData, float* outputData, const std::vector<float>& kernel, int channels) {
    int radius = (kernel.size() - 1) / 2;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            for (int c = 0; c < channels; ++c) {
                float sum = 0.0f;
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

// 高斯模糊主函式 (Float 版本)
void gaussianBlurFloatStb(int width, int height, const unsigned char* rgbData, unsigned char* blurredRgbData, double sigma, int radius) {
    int channels = 3; // RGB 影像有 3 個通道
    std::vector<float> kernel = generateGaussianKernel(radius, sigma);

    float* tempHorizontal = new float[width * height * channels];
    float* tempVertical = new float[width * height * channels];

    // 先進行水平模糊
    applyGaussianBlurHorizontal(width, height, rgbData, tempHorizontal, kernel, channels);
    // 再進行垂直模糊
    applyGaussianBlurVertical(width, height, tempHorizontal, tempVertical, kernel, channels);

    // 將結果轉換回 unsigned char (0-255)
    for (int i = 0; i < width * height * channels; ++i) {
        blurredRgbData[i] = static_cast<unsigned char>(std::min(255.0f, std::max(0.0f, tempVertical[i] * 255.0f)));
    }

    delete[] tempHorizontal;
    delete[] tempVertical;
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

// 高斯模糊主函式 (Float 版本)
void gaussianBlurPosit64Stb(int width, int height, const unsigned char* rgbData, unsigned char* blurredRgbData, double sigma, int radius) {
    int channels = 3; // RGB 影像有 3 個通道
    std::vector<Posit64> kernel = generateGaussianKernelPosit(radius, sigma);

    Posit64* tempHorizontal = new Posit64[width * height * channels];
    Posit64* tempVertical = new Posit64[width * height * channels];

    // 先進行水平模糊
    applyGaussianBlurHorizontalPosit(width, height, rgbData, tempHorizontal, kernel, channels);
    // 再進行垂直模糊
    applyGaussianBlurVerticalPosit(width, height, tempHorizontal, tempVertical, kernel, channels);

    // 將結果轉換回 unsigned char (0-255)
    for (int i = 0; i < width * height * channels; ++i) {
        blurredRgbData[i] = static_cast<unsigned char>((int)Posit_floor(tempVertical[i] * Posit64(255.0)));
    }

    delete[] tempHorizontal;
    delete[] tempVertical;
}

#include <gmp.h>
#include <mpfr.h>
#include <iostream>
#include <vector>
#include <cmath> // For std::min, std::max if needed for conversions
#include <iomanip> // For std::setprecision, if needed for debugging output

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

// 高斯模糊主函式 (MPFR 版本)
// rgbData: 輸入原始的 unsigned char 影像資料
// blurredRgbData: 輸出模糊後的 unsigned char 影像資料
void gaussianBlurMpfrStb(int width, int height, const unsigned char* rgbData, unsigned char* blurredRgbData, double sigma_double, int radius) {
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

    // 6. 將結果轉換回 unsigned char (0-255)
    for (int i = 0; i < width * height * channels; ++i) {
        // 從 MPFR 轉換回 double，然後乘以 255
        double val = mpfr_get_d(tempVerticalMpfr[i], MPFR_RNDN) * 255.0;
        // 飽和處理並轉換為 unsigned char
        blurredRgbData[i] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, val)));
    }

    // 7. 清除所有 MPFR 變數和陣列
    mpfr_clear(sigma_mpfr);

    for (int i = 0; i < width * height * channels; ++i) {
        mpfr_clear(inputMpfrData[i]);
        mpfr_clear(tempHorizontalMpfr[i]);
        mpfr_clear(tempVerticalMpfr[i]);
    }
    delete[] inputMpfrData;
    delete[] tempHorizontalMpfr;
    delete[] tempVerticalMpfr;

    // 清除 kernel 中的 mpfr_t
    for (size_t i = 0; i < kernel.size(); ++i) {
        mpfr_clear(kernel[i]);
    }
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
    if(n2 != n1) return "-1.0";

    std::string str1 = num1;
    std::string str2 = num2;
    if(str1 < str2) str1.swap(str2);
    
    reverse(str1.begin(), str1.end());
    reverse(str2.begin(), str2.end());
    

    int carry = 0;

    for (int i = 0; i < n2; i++) {
        if(str1[i] == '.'){
            result.push_back('.');
            continue;
        }

        int sub = ((str1[i] - '0') - (str2[i] - '0') - carry);

        if (sub < 0) {
            sub += 10;
            carry = 1;
        } else {
            carry = 0;
        }

        result.push_back(sub + '0');
    }

    for (int i = n2; i < n1; i++) {
        int sub = ((str1[i] - '0') - carry);

        if (sub < 0) {
            sub += 10;
            carry = 1;
        } else {
            carry = 0;
        }

        result.push_back(sub + '0');
    }

    reverse(result.begin(), result.end());

    return result;
}

string toString(mpfr_t input){
    char buffer[20];

    mpfr_sprintf(buffer, "%.18Rf", input);
    string num1(buffer);

    return num1;
}

string toString(Posit64 input) {
    string num2;
    ostringstream out;

    out << fixed << setprecision(21) << input;
    num2 = out.str();
    int length = num2.length();
    size_t site = num2.find('.');

    if (site != string::npos) { // 檢查是否有小數點
        int prec = length - site - 1;
        while (prec < 18) {
            num2.push_back('0');
            prec++;
        }
        while (prec > 18) {
            num2.pop_back();
            prec--;
        }
    } else {
        num2 += ".";
        for (int i = 0; i < 18; ++i) {
            num2.push_back('0');
        }
    }

    out.str(""); // 清空 ostringstream (雖然在這裡可能不是必要的)

    return num2;
}

int main() {
    std::string folderPath = "testImg";
    double sigma = 1.0; // 高斯模糊的標準差
    int radius = 2;     // 高斯核心的半徑 (例如，半徑為 2 會生成 5x5 的核心)

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

                // 浮點數處理 - 高斯模糊
                unsigned char* blurredRgbDataFloat = new unsigned char[width * height * 3];
                gaussianBlurFloatStb(width, height, rgbData, blurredRgbDataFloat, sigma, radius);

                // Posit64 處理 
                unsigned char* blurredRgbDataPosit64 = new unsigned char[width * height * 3];
                gaussianBlurPosit64Stb(width, height, rgbData, blurredRgbDataPosit64, sigma, radius);

                // MPFR 處理
                unsigned char* blurredRgbDataMpfr = new unsigned char[width * height * 3];
                gaussianBlurMpfrStb(width, height, rgbData, blurredRgbDataMpfr, sigma, radius);

                //將數值轉str並進行誤差運算
                vector<double> IEEE, POS;
                for (int i = 0; i < width * height * 3; i++) {
                    std::string posit = toString(blurredRgbDataFloat[i]);
                    std::string ieee = toString(blurredRgbDataPosit64[i]);
                    std::string mpfr = toString(blurredRgbDataMpfr[i]);

                    double Posit_result = stod(Difference(mpfr, posit));
                    double IEEE754_result = stod(Difference(mpfr, ieee));
                    IEEE.push_back(IEEE754_result);
                    POS.push_back(Posit_result);
                }
                std::cout << fixed  << setprecision(50) << "IEEE RMSE:" << RMSE(IEEE) << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT RMSE:" << RMSE(POS) << std::endl;

                // 產生輸出檔案名稱
                std::string filenameOnly = fs::path(filename).filename().string(); // 只取得檔案名稱
                std::string filenameWithoutExt = getFilenameWithoutExtension(filenameOnly); // 取得不含副檔名的檔案名稱

                std::string outputFilenameFloat = "output/" + filenameWithoutExt + "_hsv_float." + extension;
                stbi_write_png(outputFilenameFloat.c_str(), width, height, 3, blurredRgbDataFloat, width * 3);

                std::string outputFilenamePosit64 = "output/" + filenameWithoutExt + "_hsv_Posit64." + extension;
                stbi_write_png(outputFilenamePosit64.c_str(), width, height, 3, blurredRgbDataPosit64, width * 3);

                stbi_image_free(rgbData);
                delete[] blurredRgbDataFloat;
                delete[] blurredRgbDataPosit64;
                delete[] blurredRgbDataMpfr;

                std::cout << "浮點數 HSV 影像已儲存為: " << outputFilenameFloat << std::endl;
                std::cout << "Posit64 HSV 影像已儲存為: " << outputFilenamePosit64 << std::endl;
            }
        }
    }

    return 0;
}
