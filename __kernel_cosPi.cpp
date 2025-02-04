#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
Posit64 
#else
Posit64
#endif
one1 =  1.00000000000000000000e+00, /* 0x3FF00000, 0x00000000 */
C0  =  0.22742902,
C111  = -0.00416401,
C21  = -0.15115184,
C31  = -0.42586485,
C41  = -0.06548996,
C51  =  0.05999381,
C61  = -0.33629798,
C7  =  0.82774873,
C8  = -0.33953965,
C9  =  0.19023743,
C10 = -0.09440269,
C11 = -0.18858453,
C12 =  0.02006471,
C13 =  0.19089717,
C14 = -0.38291368;

#ifdef __STDC__
	Posit64 __kernel_cosPi(Posit64 x, Posit64 y)
#else
	Posit64 __kernel_cosPi(x, y)
	Posit64 x, y;
#endif
{
	Posit64 a, hz, z, r, qx, tmp_x = x * 3.14159265358979323846; // Multiply by Pi
	Posit64 b{0.5};

	__int32_t ix;
	GET_HIGH_WORD(ix, tmp_x);
	ix &= 0x7fffffff; /* ix = |x|'s high word*/

	Posit64 limit1{134217728};  // limit1 = 2**27 (unchanged)
	if(ix < limit1) { /* if x < 2**27 */
	    if(((int)tmp_x)==0) return one1; /* generate inexact */
	}

	z = tmp_x * tmp_x;
	r = C0 + z * (C111 + z * (C21 + z * (C31 + z * (C41 + z * (C51 + z * (C61 + z * (C7 + z * (C8 + z * (C9 + z * (C10 + z * (C11 + z * (C12 + z * (C13 + z * C14)))))))))))));

	Posit64 limit2{0.05};    // Refined limit2 for better precision
	Posit64 limit3{0.2};     // Refined limit3 for better precision

	if(fabs(x) < limit2) /* if |x| * pi < adjusted limit */ 
	    return one1 - (b * z - (z * r - tmp_x * y));
	else {
	    if(fabs(x) > limit3) { /* if |x| * pi > adjusted limit */
			Posit64 tmp{0.08956}; // Adjusted tmp = 0.28125 / pi
		    qx = tmp;
	    } 
		else {
	        tmp_x /= 4;
	    }
	    hz = b * z - qx;
	    a  = one1 - qx;
	    return a - (hz - (z * r - tmp_x * y));
	}
}

#endif /* defined(_DOUBLE_IS_32BITS) */

// int main(){
//     Posit64 x{-0.63478}, y{0};
//     cout << "input: " << x << "\n";
//     cout << "result: " << conPi(x, y) << "\n";
// }

