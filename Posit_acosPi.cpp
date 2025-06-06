#include "myfdlibm.h"
#include <iostream>
#include <iomanip>
using namespace std; 

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
static const Posit64 
#else
static Posit64
#endif
one=  1.00000000000000000000e+00, /* 0x3FF00000, 0x00000000 */
pi =  3.14159265358979311600e+00, /* 0x400921FB, 0x54442D18 */
pio2_hi =  1.57079632679489655800e+00, /* 0x3FF921FB, 0x54442D18 */
pio2_lo =  6.12323399573676603587e-17, /* 0x3C91A626, 0x33145C07 */
pS0 =  1.66666666666666657415e-01, /* 0x3FC55555, 0x55555555 */
pS1 = -3.25565818622400915405e-01, /* 0xBFD4D612, 0x03EB6F7D */
pS2 =  2.01212532134862925881e-01, /* 0x3FC9C155, 0x0E884455 */
pS3 = -4.00555345006794114027e-02, /* 0xBFA48228, 0xB5688F3B */
pS4 =  7.91534994289814532176e-04, /* 0x3F49EFE0, 0x7501B288 */
pS5 =  3.47933107596021167570e-05, /* 0x3F023DE1, 0x0DFDF709 */
qS1 = -2.40339491173441421878e+00, /* 0xC0033A27, 0x1C8A2D4B */
qS2 =  2.02094576023350569471e+00, /* 0x40002AE5, 0x9C598AC8 */
qS3 = -6.88283971605453293030e-01, /* 0xBFE6066C, 0x1B8D0159 */
qS4 =  7.70381505559019352791e-02; /* 0x3FB3B8C5, 0xB12E9282 */

#ifdef __STDC__
	Posit64 Posit_acosPi(Posit64 x)
#else
	Posit64 Posit_acosPi(x)
	Posit64 x;
#endif
{
	if(x == NAR)
		return NAR;

	Posit64 z,p,q,r,w,s,c,df, fabs_x = Posit_fabs(x);
	__int32_t hx,ix;

	GET_HIGH_WORD(hx,x);
	ix = hx&0x7fffffff;

    Posit64 limit1{1}, limit2{0.5}, limit3{6.9388939e-18}, zero{0};
	if(fabs_x >= limit1) {	/* |x| >= 1 */
	    __uint32_t lx;
	    GET_LOW_WORD(lx,x);
	    if(fabs_x == limit1) {	/* |x|==1 */
		if(x > zero) return 0.0;		/* acos(1) = 0  */
		else return (pi + 2.0*pio2_lo) / PI;	/* acos(-1)= pi */
	    }
	    return (x-x)/(x-x);		/* acos(|x|>1) is NaN */
	}
	if(fabs_x < limit2) {	/* |x| < 0.5 */
	    if(ix <= limit3) return (pio2_hi + pio2_lo) / PI;/*if|x|<2**-57*/
	    z = x*x;
	    p = z*(pS0+z*(pS1+z*(pS2+z*(pS3+z*(pS4+z*pS5)))));
	    q = one+z*(qS1+z*(qS2+z*(qS3+z*qS4)));
	    r = p/q;
	    return (pio2_hi - (x - (pio2_lo-x*r))) / PI;
	} else  if (x < zero) {		/* x < -0.5 */
	    z = (one+x)*0.5;
	    p = z*(pS0+z*(pS1+z*(pS2+z*(pS3+z*(pS4+z*pS5)))));
	    q = one+z*(qS1+z*(qS2+z*(qS3+z*qS4)));
	    s = Posit_sqrt(z);
	    r = p/q;
	    w = r*s-pio2_lo;
	    return (pi - 2.0*(s+w)) / PI;
	} else {			/* x > 0.5 */
	    z = (one-x)*0.5;
	    s = Posit_sqrt(z);
	    df = s;
	    SET_LOW_WORD(df,0);
	    c  = (z-df*df)/(s+df);
	    p = z*(pS0+z*(pS1+z*(pS2+z*(pS3+z*(pS4+z*pS5)))));
	    q = one+z*(qS1+z*(qS2+z*(qS3+z*qS4)));
	    r = p/q;
	    w = r*s+c;
	    return (2.0*(df+w)) / PI;
	}
}

#endif /* defined(_DOUBLE_IS_32BITS) */

/*
#include <limits>
#include <cmath>

int main(void)
{
	double nan = std::nan("0");
	double inf = std::numeric_limits<double>::infinity();
	double max = 1;
	double min = 1e-15;
	
	cout << fixed << setprecision(16) << "nan: " << Posit_acosPi(nan) << "\n";
	cout << fixed << setprecision(16) << "inf: " << Posit_acosPi(inf) << "\n";
	cout << fixed << setprecision(16) << "max: " << Posit_acosPi(max) << "\n";
	cout << fixed << setprecision(16) << "min: " << Posit_acosPi(min) << "\n";
}
*/

// int main(){
//     Posit64 x{-0.6124};
//     cout << "Input: " << x << "\n";
//     cout << fixed << setprecision(25) << "Result: " << Posit_acos(x) << "\n";
// }
