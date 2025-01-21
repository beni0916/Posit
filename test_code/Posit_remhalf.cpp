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
    Posit64 pi_over2{0.5}, pi{3.141592025756836};
    __int32_t n = 0;

    while(x > pi_over2){
        x = x - pi_over2;
        n++;
    }

    while(x < -pi_over2){
        x = x + pi_over2;
        n--;
    }

    if(x < 0){
        x = x + pi_over2;
        n--;
    }

    y[0] = x * pi;
    y[1] = 0;

    return n;
}

// int main(){
//     Posit64 value{-123};
//     Posit64 y[2] = {0, 0};
//     cout<< Posit_rempio2(value, y)<<endl;
//     cout<< y[0] <<endl;
// }
