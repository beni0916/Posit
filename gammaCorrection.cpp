#include "testTool.h"

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