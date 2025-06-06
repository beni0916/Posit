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
ln2 = 0.69314718055994530942,
two54 = 1.80143985094819840000e+16,	/* 43500000 00000000 */
Lg1 = 6.666666666666735130e-01,	/* 3FE55555 55555593 */
Lg2 = 3.999999999940941908e-01,	/* 3FD99999 9997FA04 */
Lg3 = 2.857142874366239149e-01,	/* 3FD24924 94229359 */
Lg4 = 2.222219843214978396e-01,	/* 3FCC71C5 1D8E78AF */
Lg5 = 1.818357216161805012e-01,	/* 3FC74664 96CB03DE */
Lg6 = 1.531383769920937332e-01,	/* 3FC39A09 D078C69F */
Lg7 = 1.479819860511658591e-01,	/* 3FC2F112 DF3E5244 */
zero = 0.0;

#ifdef __STDC__
	Posit64 Posit_log2Plus1(Posit64 x)
#else
	Posit64 Posit_log2Plus1(x)
	Posit64 x;
#endif
{
	Posit64 hfsq, f, s, z, R, w, t1, t2, dk;
	__int32_t k, hx, i, j;
	__uint32_t lx;
	x = x + Posit64{1};
	
	EXTRACT_WORDS(hx, lx, x);

	k = 0;
	if (hx < 0x00100000)
	{									/* x < 2**-1022  */
		if (((hx & 0x7fffffff) | lx) == 0)
			return -two54 / (x - x);	/* log(+-0)=-inf */
		if (hx < 0)
			return (x - x) / (x - x);	/* log(-#) = NaN */
		k -= 54;
		x *= two54;						/* subnormal number, scale up x */
		GET_HIGH_WORD(hx, x);
	}
	if (hx >= 0x7ff00000)
		return x + x;
	
	int data[5] = {0};
	getRegime(hx, data);
	getExponent(hx, data);

	int ex = 1;
	for(int i = 0; i < N; i++){
		ex *= 2;
	}

	k += (ex * data[0] + data[2]);
	hx &= (0x7fffffff >> (N + data[1])); 
	int offset = 12 - (N + 1 + data[1]);
	
	i = (hx + (0x95f64 << offset)) & (0x100000 << offset); 
	k += (i >> (20 + offset));
	Posit64 Posit_k{k};
	dk = Posit_k;
	
	Posit64 two{2};
    	while(x >= Posit_sqrt(two))
    	{
        	x /= 2;    
    	}

	while(x <= (Posit_sqrt(two) / two))
	{
        	x *= 2;    
    	}
	
	f = x - 1.0;
	
	
	if ((0x000fffff & (2 + hx)) < 3)
	{									/* |f| < 2**-20 */
		if (f == zero)
			return dk;
		
		R = f * f * (0.5 - 0.33333333333333333 * f);
		return dk - (R - f) / ln2;
	}
	s = f / (2.0 + f);
	z = s * s;
	i = hx - (0x6147a << offset);
	w = z * z;
	j = (0x6b851 << offset) - hx;
	t1 = w * (Lg2 + w * (Lg4 + w * Lg6));
	t2 = z * (Lg1 + w * (Lg3 + w * (Lg5 + w * Lg7)));
	i |= j;
	R = t2 + t1;
	if (i > 0)
	{
		hfsq = 0.5 * f * f;
		return dk - ((hfsq - (s * (hfsq + R))) - f) / ln2;
	} 
	else
	{
		return dk - ((s * (f - R)) - f) / ln2;
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
	double max = 1e10;
	double min = 1e-15;
	
	cout << fixed << setprecision(16) << "nan: " << Posit_log2Plus1(nan) << "\n";
	cout << fixed << setprecision(16) << "inf: " << Posit_log2Plus1(inf) << "\n";
	cout << fixed << setprecision(16) << "max: " << Posit_log2Plus1(max) << "\n";
	cout << fixed << setprecision(16) << "min: " << Posit_log2Plus1(min) << "\n";
}
*/

// int main(){
//     Posit64 x{45.2523};
//     cout << fixed << setprecision(20) << Posit_log10(x) << "\n";
// }

