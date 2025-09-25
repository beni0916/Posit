#include "testTool.h"
using namespace std;

string MPFR(double input){
    mpfr_t x, result, one, m;
    char buffer[20];

    mpfr_init2(x, 256);
    mpfr_init2(result, 256);
    mpfr_init2(one, 256);
    mpfr_init2(m, 256);
    
    mpfr_set_d(x, input, MPFR_RNDN); 
    mpfr_set_d(one, 1.0, MPFR_RNDN);
    mpfr_exp(result, x, MPFR_RNDN);
    mpfr_sub(m, result, one, MPFR_RNDN);

    string num1 = toString(m);

    mpfr_clear(x);
    mpfr_clear(result);
    mpfr_clear(one);
    mpfr_clear(m);

    return num1;
}

string POSIT(double input){
    string num2;
    Posit64 ans = Posit_expMinus1(input);
    return toString(ans);
}

string IEEE754(double input){
    string num3;
    double ans = exp(input) - 1;
    return toString(ans);
}

int main() {
    ofstream ofs;
    random_device rd;
    mt19937 generator(rd());   
    mpfr_t result;
    double input;
    string num1, num2, num3;
    vector<double> IEEE, POS;
    ofs.open("output/Posit_expM1.csv");
    ofs <<  "index,input,MPFR,Posit64,Double,Posit64D,DoubleD,winner" << "\n";

    for(int i = 0; i < 1000; i++){                       
        uniform_real_distribution<double> range(0, 1);  
        input = range(generator);

        num1 = MPFR(input);
        num2 = POSIT(input);
        num3 = IEEE754(input);
        
        string Posit_string = difference(num1, num2);
        double Posit_result = stod(Posit_string);
        string IEEE754_string = difference(num1, num3);
        double IEEE754_result = stod(IEEE754_string);
        IEEE.push_back(IEEE754_result);
        POS.push_back(Posit_result);

        ofs << i << "," << input << "," << num1 << "," << num2 << "," << num3 << "," << Posit_string << "," << IEEE754_string << ",";
        if (abs(Posit_result)<abs(IEEE754_result)) {
            ofs << "Posit" << '\n';
        }
        else if (abs(Posit_result)>abs(IEEE754_result)) {
            ofs << "Double" << '\n';
        }
        else {
            ofs << "EQUAL" << '\n';
        }
    }
    ofs << "Posit64Rmse,DoubleRmse" << "\n";
    ofs << RMSE(POS) << "," << RMSE(IEEE) << "\n";
    ofs.close();
}