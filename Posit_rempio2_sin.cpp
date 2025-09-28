#include <iostream>
#include "myfdlibm.h"
using namespace std;

#ifdef __STDC__
	__int32_t Posit_rempio2_sin(Posit64 x, Posit64 *y)
#else
	__int32_t Posit_rempio2_sin(x)
	Posit64 x, y[];
#endif
{   
	Posit64 pi_over2 = PI / Posit64{2};
	Posit64 d_pi = 2 * PI, zero{0};
	__int32_t n = 0;

	// become circle
	while(x >= d_pi)
	{
		x = x - d_pi;
	}

	while(x < zero)
	{
		x = x + d_pi;
	}


	
	//cout << x << endl;
	if(PI > x && x >= pi_over2)
		x = PI - x;
	else if((PI + pi_over2) > x && x >= PI)
		x = -(x - PI);
	else if(2 * PI > x && x >= PI)
		x = x - 2 * PI;
	
	// a
	if(x >= pi_over2 / 2)
		n = 0;
	else if(pi_over2 >x && x >= zero)
		n = 1;
	else if(zero > x && x >= (-pi_over2 / 2))
		n = 2;
	else
		n = 3;
	//cout << n << endl;
	//cout << x << endl;
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
