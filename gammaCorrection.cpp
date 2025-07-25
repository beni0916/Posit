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

// 函式：應用伽馬校正 (Float 版本)
void applyGammaCorrectionFloat(int width, int height, const unsigned char* inputData, double* outputData, double gamma, int channels) {
    for (int i = 0; i < width * height * channels; ++i) {
        double normalizedPixel = static_cast<double>(inputData[i]) / 255.0;
        outputData[i] = std::pow(normalizedPixel, gamma);
    }
}

// 函式：應用伽馬校正 (Posit64 版本)
void applyGammaCorrectionPosit64(int width, int height, const unsigned char* inputData, Posit64* outputData, Posit64 gamma, int channels) {
    for (int i = 0; i < width * height * channels; ++i) {
        Posit64 normalizedPixel = Posit64(inputData[i]) / Posit64(255.0);
        outputData[i] = Posit_pow(normalizedPixel, gamma);
    }
}

// 函式：應用伽馬校正 (MPFR 版本)
void applyGammaCorrectionMpfr(int width, int height, const unsigned char* inputData, mpfr_t* outputData, mpfr_t gamma, int channels, mpfr_prec_t prec) {
    mpfr_t normalizedPixel, tempVal, div255;
    mpfr_init2(normalizedPixel, prec);
    mpfr_init2(tempVal, prec);
    mpfr_init2(div255, prec);
    mpfr_set_d(div255, 255.0, MPFR_RNDN);

    for (int i = 0; i < width * height * channels; ++i) {
        mpfr_set_d(normalizedPixel, static_cast<double>(inputData[i]), MPFR_RNDN);
        mpfr_div(normalizedPixel, normalizedPixel, div255, MPFR_RNDN); // 正規化
        mpfr_pow(outputData[i], normalizedPixel, gamma, MPFR_RNDN); // 應用伽馬校正
    }
    mpfr_clear(normalizedPixel);
    mpfr_clear(tempVal);
    mpfr_clear(div255);
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

std::string difference(std::string& num1, std::string& num2) {
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

    size_t firstDigit = result.find_first_not_of('0');
    size_t dotPosResult = result.find('.');

    if (firstDigit == std::string::npos) { // All zeros
        return "0.0";
    } else if (firstDigit > dotPosResult && dotPosResult != std::string::npos) { // e.g., 00.123
        return "0" + result.substr(dotPosResult);
    } else if (firstDigit == dotPosResult && dotPosResult != std::string::npos) { // e.g., .123 -> 0.123
        return "0" + result.substr(firstDigit);
    } else {
        return result.substr(firstDigit);
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
std::string toString(Posit64 input) {
    std::ostringstream out;

    // 將 Posit64 強制轉換為 long double，並使用不同的變數名
    long double numericInput = static_cast<long double>(input); // <-- 這裡改名了

    // 對 ostringstream 應用 fixed 和 setprecision
    out << std::fixed << std::setprecision(25) << numericInput;
    
    // 從 ostringstream 獲取字串，並使用 'result' 作為變數名
    std::string result = out.str(); // <-- 這裡使用 'result'

    // 將字串結果傳遞給 formatPositiveFloatString 函數
    return formatPositiveFloatString(result); // <-- 這裡傳遞 'result'
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
    double gammaValue = 2.2; // 伽馬校正的伽馬值

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

                // ***** 獲取不同精度的伽馬校正浮點數結果 *****
                double* correctedFloatResult = new double[width * height * channels];
                applyGammaCorrectionFloat(width, height, rgbData, correctedFloatResult, gammaValue, channels);

                Posit64* correctedPosit64Result = new Posit64[width * height * channels];
                applyGammaCorrectionPosit64(width, height, rgbData, correctedPosit64Result, Posit64(gammaValue), channels);

                mpfr_t* correctedMpfrResult = new mpfr_t[width * height * channels];
                mpfr_prec_t prec = 256; // MPFR 精度
                mpfr_t gammaMpfr;
                mpfr_init2(gammaMpfr, prec);
                mpfr_set_d(gammaMpfr, gammaValue, MPFR_RNDN);
                for(int i = 0; i < width * height * channels; ++i) { // 初始化 MPFR 陣列元素
                    mpfr_init2(correctedMpfrResult[i], prec);
                }
                applyGammaCorrectionMpfr(width, height, rgbData, correctedMpfrResult, gammaMpfr, channels, prec);
                mpfr_clear(gammaMpfr); // 清理 gammaMpfr

                // 將數值轉str並進行誤差運算 (在浮點數層級)
                vector<double> ieeeRMSEVals, posRMSEVals;
                std::string filenameOnly = fs::path(filename).filename().string(); // 只取得檔案名稱
                std::string filenameWithoutExt = getFilenameWithoutExtension(filenameOnly); // 取得不含副檔名的檔案名稱
                // 創建一個新的文字檔案來儲存浮點數值
                std::string outputTextFilename = "output/" + filenameWithoutExt + "_gamma_fp_values.txt";
                std::ofstream outputFile(outputTextFilename);
                if (!outputFile.is_open()) {
                    std::cerr << "無法開啟浮點數輸出檔案: " << outputTextFilename << std::endl;
                }

                for (int i = 0; i < width * height * 3; i++) {
                    std::string ieeeStr = toString(correctedFloatResult[i]);
                    std::string positStr = toString(correctedPosit64Result[i]);
                    std::string mpfrStr = toString(correctedMpfrResult[i]);

                    // 將結果寫入檔案
                    outputFile << "Pixel " << i << ":\n";
                    outputFile << "  MPFR: " << mpfrStr << "\n";
                    outputFile << "  IEEE: " << ieeeStr << "\n";
                    outputFile << "  Posit: " << positStr << "\n";

                    double positResultDiff = stod(difference(mpfrStr, positStr));
                    double ieee754ResultDiff = stod(difference(mpfrStr, ieeeStr));

                    outputFile << "  IEEE d: " << ieee754ResultDiff << "\n";
                    outputFile << "  Posit d: " << positResultDiff << "\n";
                    outputFile << "--------------------\n";

                    ieeeRMSEVals.push_back(ieee754ResultDiff);
                    posRMSEVals.push_back(positResultDiff);
                }
                outputFile.close(); // 關閉浮點數結果檔案
                std::cout << "浮點數結果已寫入: " << outputTextFilename << std::endl;

                // ***** 在將浮點數結果轉換回 0-255 之前，先計算 RMSE *****
                std::cout << fixed  << setprecision(50) << "IEEE RMSE (double vs MPFR):" << RMSE(ieeeRMSEVals) << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT RMSE (Posit64 vs MPFR):" << RMSE(posRMSEVals) << std::endl;

                // ***** 將浮點數結果轉換回 0-255 的 unsigned char 陣列，以便儲存為圖像 *****
                unsigned char* correctedRgbDataFloat = new unsigned char[width * height * 3];
                unsigned char* correctedRgbDataPosit64 = new unsigned char[width * height * 3];
                unsigned char* correctedRgbDataMpfr = new unsigned char[width * height * 3]; // 用於儲存 MPFR 結果量化後的圖像

                for (int i = 0; i < width * height * 3; ++i) {
                    // 將 float 結果轉換回 unsigned char (0-255)
                    correctedRgbDataFloat[i] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, correctedFloatResult[i] * 255.0)));
                    // 將 Posit64 結果轉換回 unsigned char (0-255)
                    correctedRgbDataPosit64[i] = static_cast<unsigned char>((int)Posit_floor(correctedPosit64Result[i] * Posit64(255.0)));
                    // 將 MPFR 結果轉換回 unsigned char (0-255)
                    correctedRgbDataMpfr[i] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, mpfr_get_d(correctedMpfrResult[i], MPFR_RNDN) * 255.0)));
                }

                // 產生輸出影像檔案名稱並儲存
                std::string outputFilenameFloat = "output/" + filenameWithoutExt + "_gamma_corrected_float." + extension;
                stbi_write_png(outputFilenameFloat.c_str(), width, height, 3, correctedRgbDataFloat, width * 3);

                std::string outputFilenamePosit64 = "output/" + filenameWithoutExt + "_gamma_corrected_Posit64." + extension;
                stbi_write_png(outputFilenamePosit64.c_str(), width, height, 3, correctedRgbDataPosit64, width * 3);
                
                std::string outputFilenameMpfr = "output/" + filenameWithoutExt + "_gamma_corrected_Mpfr." + extension;
                stbi_write_png(outputFilenameMpfr.c_str(), width, height, 3, correctedRgbDataMpfr, width * 3);


                // ***** 記憶體釋放 *****
                stbi_image_free(rgbData);
                delete[] correctedFloatResult;
                delete[] correctedPosit64Result;
                // MPFR 陣列需要逐個清除 mpfr_t 元素，然後再刪除陣列本身
                for (int i = 0; i < width * height * 3; ++i) {
                    mpfr_clear(correctedMpfrResult[i]);
                }
                delete[] correctedMpfrResult;
                delete[] correctedRgbDataFloat;
                delete[] correctedRgbDataPosit64;
                delete[] correctedRgbDataMpfr;

                std::cout << "浮點數伽馬校正影像已儲存為: " << outputFilenameFloat << std::endl;
                std::cout << "Posit64 伽馬校正影像已儲存為: " << outputFilenamePosit64 << std::endl;
                std::cout << "MPFR 伽馬校正影像已儲存為: " << outputFilenameMpfr << std::endl;
            }
        }
    }

    return 0;
}