#include "testTool.h"

namespace fs = std::filesystem;
using namespace std;

double RMSE(const std::vector<double> &vec) {
    if (vec.empty()) {
        return 0.0; // 或者拋出例外，視你的需求而定
    }

    double sumOfSquares = std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0);
    double meanSquareError = sumOfSquares / vec.size();
    return std::sqrt(meanSquareError);
}

std::string difference(std::string num1, std::string num2) {
    std::string result = "";
    std::string str1 = num1;
    std::string str2 = num2;
    //std::cout << num1 << " " << num2 << endl;

    if (str1[0]=='-'&&str2[0]=='-') {
        str1 = str1.substr(1);
        str2 = str2.substr(1);
    } 
    else if (str1[0]=='-'||str2[0]=='-') {
        if (str1[0]=='-') {
            str1 = str1.substr(1);
        }
        else {
            str2 = str2.substr(1);
        }
        return add(str1,str2);
    }
    
    int n1 = str1.length();
    int n2 = str2.length();
    if(n2 != n1) return "-48763.0"; // 如果長度不同，表示有問題，返回錯誤值

    // 確保 str1 >= str2，方便減法
    if(str1.compare(str2) < 0) str1.swap(str2);

    reverse(str1.begin(), str1.end());
    reverse(str2.begin(), str2.end());
    
    int carry = 0;
    for (int i = 0; i < n1; i++) { // 現在 n1 == n2
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

    reverse(result.begin(), result.end());

    return result;
}

std::string add(std::string num1, std::string num2) {
    int carry = 0;
    int n1 = num1.length();
    std::string result = "";
    reverse(num1.begin(), num1.end());
    reverse(num2.begin(), num2.end());
    for (int i = 0; i < n1; i++) {
        if(num1[i] == '.'){
            result.push_back('.');
            continue;
        }

        int add = ((num1[i] - '0') + (num2[i] - '0') + carry);

        if (add > 9) {
            add -= 10;
            carry = 1;
        } else {
            carry = 0;
        }
        result.push_back(add + '0');
    }
    reverse(result.begin(), result.end());
    return result;
}

// String floatPoint 對齊
std::string formatPositiveFloatString(const std::string& inputStr, bool isNegative) {
    std::string result = inputStr;
    size_t dotPos = result.find('.');

    // --- 處理小數點前部分 ---
    std::string integerPart;
    if (dotPos == std::string::npos) {
        integerPart = result; // 如果沒有小數點，整個字串就是整數部分
    } 
    else {
        integerPart = result.substr(0, dotPos); // 取得小數點前的部分
    }

    // 移除前導零，除非它就是 "0" (例如 "007" 變成 "7", "0" 保持 "0")
    size_t firstDigit = integerPart.find_first_not_of('0');
    if (firstDigit != std::string::npos) {
        integerPart = integerPart.substr(firstDigit);
    } 
    else { // 處理空全是零和字串作為整數部分的情況 (例如輸入是 ".123" 或空字串)
        integerPart = "0";
    }

    // 補足小數點前到10位
    if (integerPart.length() < 10) {
        std::string padding(10 - integerPart.length(), '0');
        integerPart = padding + integerPart;
    }

    // --- 處理小數點後部分 ---
    std::string decimalPart;
    if (dotPos != std::string::npos) {
        decimalPart = result.substr(dotPos + 1); // 取得小數點後的部分
    }

    // 補足或截斷小數點後到70位
    if (decimalPart.length() < 70) {
        decimalPart.append(70 - decimalPart.length(), '0');
    } else if (decimalPart.length() > 70) {
        decimalPart.erase(70); // 從索引 70 開始移除，保留前 70 個字元
    }

    // --- 組合最終結果 ---
    string finalResult = integerPart + "." + decimalPart;
    if (isNegative) finalResult = "-" + finalResult;
    return finalResult;
}

std::string toString(Posit64 input) {
    bool isNegative = false;
    if (input < 0) {
        isNegative = true;
        input = -1*input;
    }
    long double numericInput = static_cast<long double>(input);
    std::ostringstream out;
    out << std::fixed << std::setprecision(70) << numericInput;
    std::string result = out.str();
    return formatPositiveFloatString(result,isNegative);
}

string toString(double input) {
    bool isNegative = false;
    if (input < 0) {
        isNegative = true;
        input = -1*input;
    }
    string numStr;
    ostringstream out;
    out << fixed << setprecision(70) << input;
    numStr = out.str();
    return formatPositiveFloatString(numStr,isNegative);
}


string toString(mpfr_t input){
    bool isNegative = (mpfr_sgn(input) < 0);
    mpfr_t abs_input;
    mpfr_init2(abs_input, 256);
    mpfr_abs(abs_input, input, MPFR_RNDN);

    char buffer[200]; 

    mpfr_sprintf(buffer, "%.70Rf", abs_input); 
    string num1(buffer);
    return formatPositiveFloatString(num1,isNegative);
}
