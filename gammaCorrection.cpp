#include "testTool.h"

namespace fs = std::filesystem;
using namespace std;
using Posit32 = sw::universal::posit<32, 2>;
using Posit16_1 = sw::universal::posit<16, 1>;
using Posit16_2 = sw::universal::posit<16, 2>;

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

// 函式：應用伽馬校正 (Posit32 版本)
void applyGammaCorrectionPosit32(int width, int height, const unsigned char* inputData, Posit32* outputData, Posit32 gamma, int channels) {
    for (int i = 0; i < width * height * channels; ++i) {
        Posit32 normalizedPixel = Posit32(inputData[i]) / Posit32(255.0);
        outputData[i] = Posit_pow(normalizedPixel, gamma);
    }
}

// 函式：應用伽馬校正 (Posit16_1 版本)
void applyGammaCorrectionPosit16_1(int width, int height, const unsigned char* inputData, Posit16_1* outputData, Posit16_1 gamma, int channels) {
    for (int i = 0; i < width * height * channels; ++i) {
        Posit16_1 normalizedPixel = Posit16_1(inputData[i]) / Posit16_1(255.0);
        outputData[i] = Posit_pow(normalizedPixel, gamma);
    }
}

// 函式：應用伽馬校正 (Posit16_2 版本)
void applyGammaCorrectionPosit16_2(int width, int height, const unsigned char* inputData, Posit16_2* outputData, Posit16_2 gamma, int channels) {
    for (int i = 0; i < width * height * channels; ++i) {
        Posit16_2 normalizedPixel = Posit16_2(inputData[i]) / Posit16_2(255.0);
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

    // 開啟 RMSE 結果輸出檔案
    std::ofstream rmseOutputFile("output/rmse_results.txt", std::ios::app);
    if (!rmseOutputFile.is_open()) {
        std::cerr << "無法開啟 RMSE 輸出檔案: output/rmse_results.txt" << std::endl;
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

                std::cout << "處理影像: " << filename << ", 寬度: " << width << ", 高度: " << height << std::endl;
                rmseOutputFile << "處理影像: " << filename << ", 寬度: " << width << ", 高度: " << height << std::endl;

                // ***** 獲取不同精度的伽馬校正浮點數結果 *****
                double* correctedFloatResult = new double[width * height * channels];
                applyGammaCorrectionFloat(width, height, rgbData, correctedFloatResult, gammaValue, channels);

                Posit64* correctedPosit64Result = new Posit64[width * height * channels];
                applyGammaCorrectionPosit64(width, height, rgbData, correctedPosit64Result, Posit64(gammaValue), channels);
                
                Posit32* correctedPosit32Result = new Posit32[width * height * channels];
                applyGammaCorrectionPosit32(width, height, rgbData, correctedPosit32Result, Posit32(gammaValue), channels);

                Posit16_1* correctedPosit16_1Result = new Posit16_1[width * height * channels];
                applyGammaCorrectionPosit16_1(width, height, rgbData, correctedPosit16_1Result, Posit16_1(gammaValue), channels);

                Posit16_2* correctedPosit16_2Result = new Posit16_2[width * height * channels];
                applyGammaCorrectionPosit16_2(width, height, rgbData, correctedPosit16_2Result, Posit16_2(gammaValue), channels);

                mpfr_t* correctedMpfrResult = new mpfr_t[width * height * channels];
                mpfr_prec_t prec = 256; // MPFR 精度
                mpfr_t gammaMpfr;
                mpfr_init2(gammaMpfr, prec);
                mpfr_set_d(gammaMpfr, gammaValue, MPFR_RNDN);
                for(int i = 0; i < width * height * channels; ++i) {
                    mpfr_init2(correctedMpfrResult[i], prec);
                }
                applyGammaCorrectionMpfr(width, height, rgbData, correctedMpfrResult, gammaMpfr, channels, prec);
                mpfr_clear(gammaMpfr);

                // 將數值轉str並進行誤差運算 (在浮點數層級)
                vector<double> ieeeRMSEVals, pos64RMSEVals, pos32RMSEVals, pos16_1RMSEVals, pos16_2RMSEVals;
                // 新增用於計算整數 RMSE 的向量
                vector<double> ieeeIntRMSEVals, pos64IntRMSEVals, pos32IntRMSEVals, pos16_1IntRMSEVals, pos16_2IntRMSEVals;

                // 新增 MPFR 轉換後的整數參考值
                unsigned char* correctedRgbDataMpfrRef = new unsigned char[width * height * 3];

                std::string filenameOnly = fs::path(filename).filename().string();
                std::string filenameWithoutExt = getFilenameWithoutExtension(filenameOnly);
                std::string outputTextFilename = "output/" + filenameWithoutExt + "_gamma_fp_values.txt";
                std::ofstream outputFile(outputTextFilename);
                if (!outputFile.is_open()) {
                    std::cerr << "無法開啟浮點數輸出檔案: " << outputTextFilename << std::endl;
                }

                for (int i = 0; i < width * height * 3; i++) {
                    std::string ieeeStr = toString(correctedFloatResult[i]);
                    std::string posit64Str = toString(correctedPosit64Result[i]);
                    std::string posit32Str = toString(correctedPosit32Result[i]);
                    std::string posit16_1Str = toString(correctedPosit16_1Result[i]);
                    std::string posit16_2Str = toString(correctedPosit16_2Result[i]);
                    std::string mpfrStr = toString(correctedMpfrResult[i]);

                    outputFile << "Pixel " << i << ":\n";
                    outputFile << "  MPFR: " << mpfrStr << "\n";
                    outputFile << "  IEEE: " << ieeeStr << "\n";
                    outputFile << "  Posit64: " << posit64Str << "\n";
                    outputFile << "  Posit32: " << posit32Str << "\n";
                    outputFile << "  Posit16_1: " << posit16_1Str << "\n";
                    outputFile << "  Posit16_2: " << posit16_2Str << "\n";

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
                    correctedRgbDataMpfrRef[i] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, mpfr_get_d(correctedMpfrResult[i], MPFR_RNDN) * 255.0)));
                    unsigned char correctedRgbDataFloatVal = static_cast<unsigned char>(std::min(255.0, std::max(0.0, correctedFloatResult[i] * 255.0)));
                    unsigned char correctedRgbDataPosit64Val = static_cast<unsigned char>((int)Posit_floor(correctedPosit64Result[i] * Posit64(255.0)));
                    unsigned char correctedRgbDataPosit32Val = static_cast<unsigned char>((int)Posit_floor(correctedPosit32Result[i] * Posit32(255.0)));
                    unsigned char correctedRgbDataPosit16_1Val = static_cast<unsigned char>((int)Posit_floor(correctedPosit16_1Result[i] * Posit16_1(255.0)));
                    unsigned char correctedRgbDataPosit16_2Val = static_cast<unsigned char>((int)Posit_floor(correctedPosit16_2Result[i] * Posit16_2(255.0)));
                    
                    ieeeIntRMSEVals.push_back(static_cast<double>(correctedRgbDataMpfrRef[i] - correctedRgbDataFloatVal));
                    pos64IntRMSEVals.push_back(static_cast<double>(correctedRgbDataMpfrRef[i] - correctedRgbDataPosit64Val));
                    pos32IntRMSEVals.push_back(static_cast<double>(correctedRgbDataMpfrRef[i] - correctedRgbDataPosit32Val));
                    pos16_1IntRMSEVals.push_back(static_cast<double>(correctedRgbDataMpfrRef[i] - correctedRgbDataPosit16_1Val));
                    pos16_2IntRMSEVals.push_back(static_cast<double>(correctedRgbDataMpfrRef[i] - correctedRgbDataPosit16_2Val));
                }
                outputFile.close();
                std::cout << "浮點數結果已寫入: " << outputTextFilename << std::endl;
                
                // ***** 輸出浮點數層級 RMSE *****
                double ieeeFpRMSE = RMSE(ieeeRMSEVals);
                double posit64FpRMSE = RMSE(pos64RMSEVals);
                double posit32FpRMSE = RMSE(pos32RMSEVals);
                double posit16_1FpRMSE = RMSE(pos16_1RMSEVals);
                double posit16_2FpRMSE = RMSE(pos16_2RMSEVals);
                
                std::cout << fixed  << setprecision(50) << "--- 浮點數層級 RMSE (vs MPFR) ---" << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "--- 浮點數層級 RMSE (vs MPFR) ---" << std::endl;
                std::cout << fixed  << setprecision(50) << "IEEE fpRMSE:" << ieeeFpRMSE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "IEEE fpRMSE:" << ieeeFpRMSE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT64 fpRMSE:" << posit64FpRMSE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT64 fpRMSE:" << posit64FpRMSE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT32 fpRMSE:" << posit32FpRMSE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT32 fpRMSE:" << posit32FpRMSE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_1 fpRMSE:" << posit16_1FpRMSE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_1 fpRMSE:" << posit16_1FpRMSE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_2 fpRMSE:" << posit16_2FpRMSE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_2 fpRMSE:" << posit16_2FpRMSE << std::endl;

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
                unsigned char* correctedRgbDataFloat = new unsigned char[width * height * 3];
                unsigned char* correctedRgbDataPosit64 = new unsigned char[width * height * 3];
                unsigned char* correctedRgbDataPosit32 = new unsigned char[width * height * 3];
                unsigned char* correctedRgbDataPosit16_1 = new unsigned char[width * height * 3];
                unsigned char* correctedRgbDataPosit16_2 = new unsigned char[width * height * 3];
                unsigned char* correctedRgbDataMpfr = new unsigned char[width * height * 3];

                for (int i = 0; i < width * height * 3; ++i) {
                    correctedRgbDataFloat[i] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, correctedFloatResult[i] * 255.0)));
                    correctedRgbDataPosit64[i] = static_cast<unsigned char>((int)Posit_floor(correctedPosit64Result[i] * Posit64(255.0)));
                    correctedRgbDataPosit32[i] = static_cast<unsigned char>((int)Posit_floor(correctedPosit32Result[i] * Posit32(255.0)));
                    correctedRgbDataPosit16_1[i] = static_cast<unsigned char>((int)Posit_floor(correctedPosit16_1Result[i] * Posit16_1(255.0)));
                    correctedRgbDataPosit16_2[i] = static_cast<unsigned char>((int)Posit_floor(correctedPosit16_2Result[i] * Posit16_2(255.0)));
                    correctedRgbDataMpfr[i] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, mpfr_get_d(correctedMpfrResult[i], MPFR_RNDN) * 255.0)));
                }

                // 產生輸出影像檔案名稱並儲存
                std::string outputFilenameFloat = "output/" + filenameWithoutExt + "_gamma_corrected_float." + extension;
                stbi_write_png(outputFilenameFloat.c_str(), width, height, 3, correctedRgbDataFloat, width * 3);

                std::string outputFilenamePosit64 = "output/" + filenameWithoutExt + "_gamma_corrected_Posit64." + extension;
                stbi_write_png(outputFilenamePosit64.c_str(), width, height, 3, correctedRgbDataPosit64, width * 3);
                
                std::string outputFilenamePosit32 = "output/" + filenameWithoutExt + "_gamma_corrected_Posit32." + extension;
                stbi_write_png(outputFilenamePosit32.c_str(), width, height, 3, correctedRgbDataPosit32, width * 3);

                std::string outputFilenamePosit16_1 = "output/" + filenameWithoutExt + "_gamma_corrected_Posit16_1." + extension;
                stbi_write_png(outputFilenamePosit16_1.c_str(), width, height, 3, correctedRgbDataPosit16_1, width * 3);

                std::string outputFilenamePosit16_2 = "output/" + filenameWithoutExt + "_gamma_corrected_Posit16_2." + extension;
                stbi_write_png(outputFilenamePosit16_2.c_str(), width, height, 3, correctedRgbDataPosit16_2, width * 3);

                std::string outputFilenameMpfr = "output/" + filenameWithoutExt + "_gamma_corrected_Mpfr." + extension;
                stbi_write_png(outputFilenameMpfr.c_str(), width, height, 3, correctedRgbDataMpfr, width * 3);

                // ***** 記憶體釋放 *****
                stbi_image_free(rgbData);
                delete[] correctedFloatResult;
                delete[] correctedPosit64Result;
                delete[] correctedPosit32Result;
                delete[] correctedPosit16_1Result;
                delete[] correctedPosit16_2Result;
                for (int i = 0; i < width * height * 3; ++i) {
                    mpfr_clear(correctedMpfrResult[i]);
                }
                delete[] correctedMpfrResult;
                delete[] correctedRgbDataFloat;
                delete[] correctedRgbDataPosit64;
                delete[] correctedRgbDataPosit32;
                delete[] correctedRgbDataPosit16_1;
                delete[] correctedRgbDataPosit16_2;
                delete[] correctedRgbDataMpfr;
                delete[] correctedRgbDataMpfrRef;


                std::cout << "浮點數伽馬校正影像已儲存為: " << outputFilenameFloat << std::endl;
                std::cout << "Posit64 伽馬校正影像已儲存為: " << outputFilenamePosit64 << std::endl;
                std::cout << "Posit32 伽馬校正影像已儲存為: " << outputFilenamePosit32 << std::endl;
                std::cout << "Posit16_1 伽馬校正影像已儲存為: " << outputFilenamePosit16_1 << std::endl;
                std::cout << "Posit16_2 伽馬校正影像已儲存為: " << outputFilenamePosit16_2 << std::endl;
                std::cout << "MPFR 伽馬校正影像已儲存為: " << outputFilenameMpfr << std::endl;
            }
        }
    }
    
    // 關閉 RMSE 檔案
    rmseOutputFile.close();

    return 0;
}