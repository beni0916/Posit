#include "myfdlibm.h"
#include <iostream>
#include <iomanip>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Posit64 Posit_tanPi(Posit64 x)
#else
	Posit64 Posit_tanPi(x)
	Posit64 x;
#endif
{
	Posit64 y[2], z=0.0, fabs_x = Posit_fabs(x);
	Posit64 nx = x * PI;
	__int32_t n,ix;

    /* High word of x. */
	GET_HIGH_WORD(ix,x);

    /* |x| ~< pi/4 */
	Posit64 limit{0.25};
	if(fabs_x <= limit) return __kernel_tan(nx, z, 1);

    /* tan(Inf or NaN) is NaN */
	//else if (ix>=0x7ff00000) return x-x;		/* NaN */

    /* argument reduction needed */
	else {
	    n = Posit_remhalf(x,y);
	    return __kernel_tan(y[0], y[1], 1-((n&1)<<1)); /*   1 -- n even
							-1 -- n odd */
	}
}

#endif /* _DOUBLE_IS_32BITS */

//int main(){
//   Posit64 x{3.7527};
//   cout << fixed << setprecision(8) << "Input: " << x << "\n";
//   cout << fixed << setprecision(20) << "Result: " << Posit_tanPi(x) << "\n";
//}
