#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "myfdlibm.h"
// g++ -std=c++17 -o test_func test_func.cpp Posit_fabs.cpp Posit_sqrt.cpp Posit_exp2.cpp Posit_floor.cpp Posit_log2.cpp Posit_pow.cpp
using namespace std;

int main(void)
{
    Posit64 a{0.0000152587890625};
    __uint32_t ha, la;

    for(int i = 0; i < 32; i++)
    {
        GET_LOW_WORD(ha, a);
    
        cout << bitset<32>(ha) << endl;
        a = a * Posit64{2};
    }
    
    

    // double c = 8583691038.52;
    // double d = log2(c);
    // cout << d << endl;
    // __uint32_t c, d;
    // EXTRACT_WORDS(c, d, b);
    // cout << bitset<32>(c) << endl;
    // cout << bitset<32>(d) << endl;

    return 0;
}