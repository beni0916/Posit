#include <iostream>
#include "myfdlibm.h"
using namespace std;

#ifdef __STDC__
	__int32_t Posit_rempio2_new(Posit64 x, Posit64 *y)
#else
	__int32_t Posit_rempio2_new(x)
	Posit64 x, y[];
#endif
{   
    Posit64 pi_over2{1.5707963267948966192313216916398};
    __int32_t n = 0;

    while(x > PI){
        x = x - PI;
    }

    while(x < 0){
        x = x + PI;
    }

    if(     x <= pi_over2                 && x > pi_over2 / 2)
    	n = 0;
    else if(x > pi_over2                  && x <= (pi_over2 + pi_over2 / 2))
    	n = 1;
    else if(x > (pi_over2 + pi_over2 / 2) && x <= PI)
    	n = 2;
    else
    	n = 3;
    
    y[0] = x;
    y[1] = 0;

    return n;
}

// int main(){
//     Posit64 value{-123};
//     Posit64 y[2] = {0, 0};
//     cout<< Posit_rempio2(value, y)<<endl;
//     cout<< y[0] <<endl;
// }
