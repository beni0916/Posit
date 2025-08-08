#include "testTool.h"

namespace fs = std::filesystem;
using namespace std;

void analyzeFloatExponentDistribution(const std::string& analysisName, const double* data, int size, std::ofstream& outputFile) {
    std::cout << "--- " << analysisName << " 指數分佈 ---" << std::endl;
    outputFile << "--- " << analysisName << " 指數分佈 ---" << std::endl;
    std::map<int, int> exponentCount;
    for (int i = 0; i < size; ++i) {
        if (data[i] > 0.0) {
            int exponent;
            frexp(data[i], &exponent);
            exponentCount[exponent]++;
        }
    }
    for (const auto& pair : exponentCount) {
        std::cout << "Exponent " << pair.first << ": " << pair.second << " 次" << std::endl;
        outputFile << "Exponent " << pair.first << ": " << pair.second << " 次" << std::endl;
    }
    // 輸出直方圖
    std::cout << "\n--- " << analysisName << " 指數分佈 (直方圖) ---" << std::endl;
    outputFile << "\n--- " << analysisName << " 指數分佈 (直方圖) ---" << std::endl;
    int maxCount = 0;
    if (!exponentCount.empty()) {
        for (const auto& pair : exponentCount) {
            if (pair.second > maxCount) {
                maxCount = pair.second;
            }
        }
    }

    const int barWidth = 50; // 直方圖最寬的長度
    for (const auto& pair : exponentCount) {
        int exponent = pair.first;
        int count = pair.second;
        
        int barLength = 0;
        if (maxCount > 0) {
            barLength = static_cast<int>((static_cast<double>(count) / maxCount) * barWidth);
        }
        std::cout << "Exponent " << std::setw(3) << exponent << ": " << std::string(barLength, '#') << " (" << count << ")" << std::endl;
        outputFile << "Exponent " << std::setw(3) << exponent << ": " << std::string(barLength, '#') << " (" << count << ")" << std::endl;
    }
    std::cout << "--------------------------------\n";
    outputFile << "--------------------------------\n";
}

double calculateMeanRelativeError(const std::vector<double>& errors, const std::vector<double>& mpfrValues) {
    if (errors.empty() || mpfrValues.empty() || errors.size() != mpfrValues.size()) {
        return 0.0;
    }

    double sumOfRelativeErrors = 0.0;
    int validCount = errors.size(); // 計算所有像素的平均值

    for (size_t i = 0; i < errors.size(); ++i) {
        if (mpfrValues[i] == 0.0) {
            // 參考值為 0 的情況
            if (errors[i] == 0.0) {
                // 如果誤差也為 0 (Posit結果也是0)，則相對誤差為 0
                sumOfRelativeErrors += 0.0;
            } else {
                // 如果誤差不為 0 (Posit結果不為0)，代表發生極大誤差
                // 將其視為一個非常大的值，避免除以零且能反映誤差的嚴重性
                sumOfRelativeErrors += DBL_MAX;
            }
        } else {
            // 一般情況：參考值不為 0
            sumOfRelativeErrors += std::abs(errors[i] / mpfrValues[i]);
        }
    }

    // 計算平均值
    return (validCount > 0) ? sumOfRelativeErrors / validCount : 0.0;
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
    out << std::fixed << std::setprecision(50) << numericInput;
    
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