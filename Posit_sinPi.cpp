#include <iostream>
#include "myfdlibm.h"
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Posit64 Posit_sinPi(Posit64 x)
#else
	Posit64 Posit_sinPi(x)
	Posit64 x;
#endif
{
	Posit64 y[2], z{0.0}, fabs_x = Posit_fabs(x);
	Posit64 nx = x * PI;
	__int32_t n, ix;

    // /* |x| ~< pi/4 */
	Posit64 limit{0.25};
	if(fabs_x <= limit) return __kernel_sin(nx, z, 0);

    // /* sin(Inf or NaN) is NaN */
	else if (x == NAR) return x-x;

    /* argument reduction needed */
	else{
		n = Posit_remhalf(x, y);
		switch(n & 3) {
			case 0: return  __kernel_sin(y[0], y[1], 1);
			case 1: return  __kernel_cos(y[0], y[1]);
			case 2: return -__kernel_sin(y[0], y[1], 1);
			default:
				return -__kernel_cos(y[0], y[1]);
		}
	}
	
}

#endif /* _DOUBLE_IS_32BITS */

// int main(){
//     Posit64 value{-12.3634};
// 	cout << "Input: " << value << "\n";
//     cout << fixed << setprecision(16) << "Result: " << sin(value) << "\n";
// }
