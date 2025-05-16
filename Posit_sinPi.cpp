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
	x = x * PI;
	Posit64 y[2], z{0.0}, fabs_x = Posit_fabs(x);
	__int32_t n, ix;

    // /* |x| ~< pi/4 */
	Posit64 limit{0.785398006439209};
	if(fabs_x <= limit) return __kernel_sin(x, z, 0);
	

    // /* sin(Inf or NaN) is NaN */
	else if (x == NAR) return x-x;

    /* argument reduction needed */
	else{
		Posit64 pi_over2{1.5707963267948966192313216916398};
		n = Posit_rempio2_sin(x, y);
		//cout << n << endl;
		switch(n)
		{
			case 0: 
				return  __kernel_cos(pi_over2 - y[0], y[1]);
			case 1: 
				return  __kernel_sin(y[0], y[1], 0);
			case 2: 
				return  -__kernel_sin(-y[0], y[1], 0);
			case 3:
				return -__kernel_cos(pi_over2 + y[0], y[1]);
			default:
				return  __kernel_sin(y[0], y[1], 0);
		}
	}
	
}

#endif /* _DOUBLE_IS_32BITS */
/*
 int main(){
     Posit64 value{12.3634};
 	cout << "Input: " << value << "\n";
     cout << fixed << setprecision(16) << "Result: " << Posit_sin(value) << "\n";
 }
 */
 //g++ -std=c++17 -o Result_sin Result_sin.cpp Posit_sin.cpp __kernel_sin.cpp __kernel_cos.cpp Posit_fabs.cpp Posit_rempio2.cpp Quire_qToP.cpp Quire_pToQ.cpp Quire_qMulAdd.cpp Quire_qSubQ.cpp Quire_qMulSub.cpp -lgmp -lmpfr
