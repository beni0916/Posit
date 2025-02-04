#include <iostream>
#include "myfdlibm.h"
using namespace std;

#ifdef __STDC__
	__int32_t Posit_remhalf(Posit64 x, Posit64 *y)
#else
	__int32_t Posit_remhalf(x)
	Posit64 x, y[];
#endif
{   
    Posit64 pi_over2{0.5};
    __int32_t n = (int)x * 2;
    Posit64 nx = x - pi_over2 * n;

    while(nx > pi_over2){
        nx = nx - pi_over2;
        n++;
    }

    while(nx < -pi_over2){
        nx = nx + pi_over2;
        n--;
    }

    if(nx < 0){
        nx = nx + pi_over2;
        n--;
    }

    y[0] = nx * PI;
    y[1] = 0;

    return n;
}

// int main(){
//     Posit64 value{-123};
//     Posit64 y[2] = {0, 0};
//     cout<< Posit_rempio2(value, y)<<endl;
//     cout<< y[0] <<endl;
// }
