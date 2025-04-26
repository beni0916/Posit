#include <iostream>
#include "myfdlibm.h"
using namespace std;

#ifdef __STDC__
	__int32_t Posit_rempio2_cos(Posit64 x, Posit64 *y)
#else
	__int32_t Posit_rempio2_cos(x)
	Posit64 x, y[];
#endif
{   
	Posit64 pi_over2{1.5707963267948966192313216916398};
	Posit64 d_pi = 2 * PI, zero{0};
	__int32_t n = 0;
	
	/*
	Quire64 qx,qpi_o2, qd_pi, qpi;
	qx	= Quire_pToQ(x);
	qpi_o2  = Quire_pToQ(pi_over2);
	qd_pi   = Quire_pToQ(d_pi);
	qpi     = Quire_pToQ(PI);
	*/

	// become circle
	while(x >= d_pi)
	{
		x = x - d_pi;
	}

	while(x < zero)
	{
		x = x + d_pi;
	}

	if(x >= PI)
		x = d_pi - x;

	// a
	if(pi_over2 / 2 > x && x >= 0)
		n = 0;
	else if(pi_over2 > x && x >= pi_over2 / 2)
		n = 1;
	else if((pi_over2 + pi_over2 / 2) > x && x >= pi_over2)
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
