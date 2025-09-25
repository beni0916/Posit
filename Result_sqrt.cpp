#include "testTool.h"
using namespace std;

string MPFR(double input){
    mpfr_t x, result;
    char buffer[20];

    mpfr_init2(x, 256);
    mpfr_init2(result, 256);
    
    mpfr_set_d(x, input, MPFR_RNDN); 
    mpfr_sqrt(result, x, MPFR_RNDN);

    string num1 = toString(result);

    mpfr_clear(x);
    mpfr_clear(result);

    return num1;
}

string POSIT(double input){
    string num2;
    Posit64 ans = Posit_sqrt(input);
    return toString(ans);
}

string IEEE754(double input){
    string num3;
    double ans = sqrt(input);
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
    ofs.open("output/Posit_sqrt.csv");
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
