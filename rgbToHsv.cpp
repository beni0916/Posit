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

// RGB to HSV 轉換函數 (浮點數版本)
void rgbToHsvFloatStb(int width, int height, const unsigned char* rgbData, float* hsvData) {
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int index = (i * width + j) * 3;
            float r = static_cast<float>(rgbData[index + 0]) / 255.0f;
            float g = static_cast<float>(rgbData[index + 1]) / 255.0f;
            float b = static_cast<float>(rgbData[index + 2]) / 255.0f;

            float h, s, v;
            float maxValue = std::max(r, std::max(g, b));
            float minValue = std::min(r, std::min(g, b));
            float delta = maxValue - minValue;

            v = maxValue;

            if (delta == 0) {
                h = 0;
                s = 0;
            } else {
                s = delta / maxValue;
                if (r == maxValue) {
                    h = (g - b) / delta;
                } else if (g == maxValue) {
                    h = 2 + (b - r) / delta;
                } else {
                    h = 4 + (r - g) / delta;
                }
                h *= 60;
                if (h < 0) {
                    h += 360;
                }
            }

            int hsvIndex = (i * width + j) * 3;
            hsvData[hsvIndex + 0] = h;
            hsvData[hsvIndex + 1] = s;
            hsvData[hsvIndex + 2] = v;
        }
    }
}

// RGB to HSV 轉換函數 (Posit64 版本)
void rgbToHsvPosit64Stb(int width, int height, const unsigned char* rgbData, Posit64* hsvData) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            int index = (i * width + j) * 3;
            Posit64 r = Posit64(rgbData[index + 0]) / Posit64(255.0);
            Posit64 g = Posit64(rgbData[index + 1]) / Posit64(255.0);
            Posit64 b = Posit64(rgbData[index + 2]) / Posit64(255.0);

            Posit64 h, s, v;
            Posit64 maxValue = std::max(r, std::max(g, b));
            Posit64 minValue = std::min(r, std::min(g, b));
            Posit64 delta = maxValue - minValue;

            v = maxValue;

            if (delta == Posit64(0)) {
                h = Posit64(0);
                s = Posit64(0);
            } else {
                s = delta / maxValue;
                if (r == maxValue) {
                    h = (g - b) / delta;
                } else if (g == maxValue) {
                    h = Posit64(2) + (b - r) / delta;
                } else {
                    h = Posit64(4) + (r - g) / delta;
                }
                h *= Posit64(60);
                if (h < Posit64(0)) {
                    h += Posit64(360);
                }
            }

            int hsvIndex = (i * width + j) * 3;
            hsvData[hsvIndex + 0] = h;
            hsvData[hsvIndex + 1] = s;
            hsvData[hsvIndex + 2] = v;
        }
    }
}

// MPFR 版本的 RGB to HSV 轉換函數 (高精度)
void rgbToHsvMpfrStb(int width, int height, const unsigned char* rgbData, mpfr_t* hsvData) {
    mpfr_prec_t prec = 256; // 设置 MPFR 精度
    mpfr_t r, g, b, maxVal, minVal, delta, temp1, temp2, temp60;

    // 初始化 MPFR 变量
    mpfr_init2(r, prec);
    mpfr_init2(g, prec);
    mpfr_init2(b, prec);
    mpfr_init2(maxVal, prec);
    mpfr_init2(minVal, prec);
    mpfr_init2(delta, prec);
    mpfr_init2(temp1, prec);
    mpfr_init2(temp2, prec);
    mpfr_init2(temp60, prec);
    mpfr_set_d(temp60, 60.0, MPFR_RNDN);

    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int index = (i * width + j) * 3;

            // 将 RGB 值转换为 MPFR 浮点数
            mpfr_set_d(r, static_cast<double>(rgbData[index + 0]) / 255.0, MPFR_RNDN);
            mpfr_set_d(g, static_cast<double>(rgbData[index + 1]) / 255.0, MPFR_RNDN);
            mpfr_set_d(b, static_cast<double>(rgbData[index + 2]) / 255.0, MPFR_RNDN);

            // 找到最大值、最小值和色差
            mpfr_max(maxVal, r, g, MPFR_RNDN);
            mpfr_max(maxVal, maxVal, b, MPFR_RNDN);
            mpfr_min(minVal, r, g, MPFR_RNDN);
            mpfr_min(minVal, minVal, b, MPFR_RNDN);
            mpfr_sub(delta, maxVal, minVal, MPFR_RNDN);

            // 计算 V
            mpfr_t& v = hsvData[(i * width + j) * 3 + 2];
            mpfr_init2(v, prec);
            mpfr_set(v, maxVal, MPFR_RNDN);

            // 如果 delta 是 0，则 S 和 H 都是 0
            if (mpfr_cmp_d(delta, 0.0) == 0) {
                mpfr_t& h = hsvData[(i * width + j) * 3 + 0];
                mpfr_init2(h, prec);
                mpfr_set_d(h, 0.0, MPFR_RNDN);

                mpfr_t& s = hsvData[(i * width + j) * 3 + 1];
                mpfr_init2(s, prec);
                mpfr_set_d(s, 0.0, MPFR_RNDN);

            } else {
                // 计算 S
                mpfr_t& s = hsvData[(i * width + j) * 3 + 1];
                mpfr_init2(s, prec);
                mpfr_div(s, delta, maxVal, MPFR_RNDN);

                // 计算 H
                mpfr_t& h = hsvData[(i * width + j) * 3 + 0];
                mpfr_init2(h, prec);
                if (mpfr_equal_p(maxVal, r)) {  // maxVal == r
                    mpfr_sub(temp1, g, b, MPFR_RNDN);
                    mpfr_div(h, temp1, delta, MPFR_RNDN);
                } else if (mpfr_equal_p(maxVal, g)) { // maxVal == g
                    mpfr_set_d(temp1, 2.0, MPFR_RNDN);
                    mpfr_sub(temp2, b, r, MPFR_RNDN);
                    mpfr_div(temp2, temp2, delta, MPFR_RNDN);
                    mpfr_add(h, temp1, temp2, MPFR_RNDN);
                } else { // maxVal == b
                    mpfr_set_d(temp1, 4.0, MPFR_RNDN);
                    mpfr_sub(temp2, r, g, MPFR_RNDN);
                    mpfr_div(temp2, temp2, delta, MPFR_RNDN);
                    mpfr_add(h, temp1, temp2, MPFR_RNDN);
                }
                mpfr_mul(h, h, temp60, MPFR_RNDN);
                if (mpfr_cmp_d(h, 0.0) < 0) {
                    mpfr_add_d(h, h, 360.0, MPFR_RNDN);
                }
            }
        }
    }

    mpfr_clear(r);
    mpfr_clear(g);
    mpfr_clear(b);
    mpfr_clear(maxVal);
    mpfr_clear(minVal);
    mpfr_clear(delta);
    mpfr_clear(temp1);
    mpfr_clear(temp2);
    mpfr_clear(temp60);
}

// HSV to RGB 轉換函數 (Posit64 版本)
void hsvToRgbPosit64Stb(int width, int height, const Posit64* hsvData, unsigned char* rgbData) {
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int hsvIndex = (i * width + j) * 3;
            Posit64 h = hsvData[hsvIndex + 0];
            Posit64 s = hsvData[hsvIndex + 1];
            Posit64 v = hsvData[hsvIndex + 2];

            Posit64 r = 0, g = 0, b = 0;
            if (s == Posit64(0)) {
                r = g = b = v; // achromatic
            } else {
                Posit64 hh = h / Posit64(60);
                int ii = static_cast<int>(Posit_floor(hh)); // 使用 posit64_to_double
                Posit64 ff = hh - Posit64(ii);
                Posit64 p = v * (Posit64(1) - s);
                Posit64 q = v * (Posit64(1) - s * ff);
                Posit64 t = v * (Posit64(1) - s * (Posit64(1) - ff));

                switch (ii) {
                    case 0:
                        r = v;
                        g = t;
                        b = p;
                        break;
                    case 1:
                        r = q;
                        g = v;
                        b = p;
                        break;
                    case 2:
                        r = p;
                        g = v;
                        b = t;
                        break;
                    case 3:
                        r = p;
                        g = q;
                        b = v;
                        break;
                    case 4:
                        r = t;
                        g = p;
                        b = v;
                        break;
                    default:
                        r = v;
                        g = p;
                        b = q;
                        break;
                }
            }

            rgbData[(i * width + j) * 3 + 0] = static_cast<unsigned char>((int)Posit_floor(r * Posit64(255.0)));
            rgbData[(i * width + j) * 3 + 1] = static_cast<unsigned char>((int)Posit_floor(g * Posit64(255.0)));
            rgbData[(i * width + j) * 3 + 2] = static_cast<unsigned char>((int)Posit_floor(b * Posit64(255.0)));
        }
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
string toString(float input) {
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

                // 浮點數處理
                float* hsvDataFloat = new float[width * height * 3];
                rgbToHsvFloatStb(width, height, rgbData, hsvDataFloat);

                // Posit64 處理 
                Posit64* hsvDataPosit64 = new Posit64[width * height * 3];
                rgbToHsvPosit64Stb(width, height, rgbData, hsvDataPosit64);

                // MPFR 處理
                mpfr_t* hsvDataMpfr = new mpfr_t[width * height * 3];
                rgbToHsvMpfrStb(width, height, rgbData, hsvDataMpfr);

                // 將 float 轉換為 Posit64
                Posit64* hsvDataFloatPosit64 = new Posit64[width * height * 3];
                for (int i = 0; i < width * height * 3; i++) {
                    hsvDataFloatPosit64[i] = Posit64(hsvDataFloat[i]);
                }

                //將數值轉str並進行誤差運算
                vector<double> IEEE, POS;
                std::string filenameOnly = fs::path(filename).filename().string(); // 只取得檔案名稱
                std::string filenameWithoutExt = getFilenameWithoutExtension(filenameOnly); // 取得不含副檔名的檔案名稱
                // 創建一個新的文字檔案來儲存浮點數值
                std::string outputTextFilename = "output/" + filenameWithoutExt + "_fp_values.txt";
                std::ofstream outputFile(outputTextFilename);
                if (!outputFile.is_open()) {
                    std::cerr << "無法開啟浮點數輸出檔案: " << outputTextFilename << std::endl;
                }
                for (int i = 0; i < width * height * 3; i++) {
                    outputFile << fixed  << setprecision(50) << "Pixel " << i << ":\n";
                    outputFile << "  Posit O: " << hsvDataPosit64[i] << "\n";
                    std::string posit = toString(hsvDataPosit64[i]);
                    std::string ieee = toString(hsvDataFloat[i]);
                    std::string mpfr = toString(hsvDataMpfr[i]);
                    outputFile << "  MPFR: " << mpfr << "\n";
                    outputFile << "  IEEE: " << ieee << "\n";
                    outputFile << "  Posit: " << posit << "\n";

                    double Posit_result = stod(Difference(mpfr, posit));
                    double IEEE754_result = stod(Difference(mpfr, ieee));
                    outputFile << "  IEEE d: " << IEEE754_result << "\n";
                    outputFile << "  Posit d: " << Posit_result << "\n";
                    if (Posit_result>IEEE754_result) {
                        outputFile << "FFFFF"<< "\n";
                    }else {
                        outputFile << "TTTTT"<< "\n";
                    }
                    outputFile << "--------------------\n";
                    IEEE.push_back(IEEE754_result);
                    POS.push_back(Posit_result);
                }
                std::cout << fixed  << setprecision(50) << "IEEE RMSE:" << RMSE(IEEE) << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT RMSE:" << RMSE(POS) << std::endl;

                // 使用 Posit64 版本的 hsvToRgb 進行轉換
                unsigned char* rgbDataOutFloat = new unsigned char[width * height * 3];
                hsvToRgbPosit64Stb(width, height, hsvDataFloatPosit64, rgbDataOutFloat);

                unsigned char* rgbDataOutPosit64 = new unsigned char[width * height * 3];
                hsvToRgbPosit64Stb(width, height, hsvDataPosit64, rgbDataOutPosit64);

                std::string outputFilenameFloat = "output/" + filenameWithoutExt + "_hsv_float." + extension;
                stbi_write_png(outputFilenameFloat.c_str(), width, height, 3, rgbDataOutFloat, width * 3);

                std::string outputFilenamePosit64 = "output/" + filenameWithoutExt + "_hsv_Posit64." + extension;
                stbi_write_png(outputFilenamePosit64.c_str(), width, height, 3, rgbDataOutPosit64, width * 3);

                stbi_image_free(rgbData);
                delete[] hsvDataFloat;
                delete[] hsvDataPosit64;
                delete[] hsvDataFloatPosit64;
                delete[] rgbDataOutFloat;
                delete[] rgbDataOutPosit64;

                std::cout << "浮點數 HSV 影像已儲存為: " << outputFilenameFloat << std::endl;
                std::cout << "Posit64 HSV 影像已儲存為: " << outputFilenamePosit64 << std::endl;
            }
        }
    }

    return 0;
}
