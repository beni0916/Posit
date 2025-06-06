#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Posit64 Posit_hypot(Posit64 x, Posit64 y)
#else
	Posit64 Posit_hypot(x, y)
	Posit64 x, y;
#endif
{
	Posit64 a, b, t1, t2, y1, y2, w, fabs_x, fabs_y;
	__int32_t j, k, ha, hb;
	x = Posit_fabs(x);
	y = Posit_fabs(y);

	GET_HIGH_WORD(ha, x);
	//ha &= (0x7fffffff);
	GET_HIGH_WORD(hb, y);
	//hb &= (0x7fffffff);
	if (hb > ha)		// let a > b
	{
		a = y;
		b = x;
		j = ha;
		ha = hb;
		hb = j;
	} else
	{
		a = x;
		b = y;
	}
	SET_HIGH_WORD(a, ha);				/* a <- |a| */
	SET_HIGH_WORD(b, hb);				/* b <- |b| */
	k = 0;
	
	if (ha > (0x5f300000))
	{			
								/* a>2**500 */
		if (ha >= (0x7ff00000))
		{								/* Inf or NaN */
			__uint32_t low;

			w = a + b;					/* for sNaN */
			GET_LOW_WORD(low, a);
			if (((ha & (0xfffff)) | low) == 0)
				w = a;
			GET_LOW_WORD(low, b);
			if (((hb ^ (0x7ff00000)) | low) == 0)
				w = b;
			return w;
		}
		/* scale a and b by 2**-600 */
		ha -= (0x25800000);
		hb -= (0x25800000);
		k += (600);
		SET_HIGH_WORD(a, ha);
		SET_HIGH_WORD(b, hb);
	}
	if (hb < (0x20b00000))
	{		
								/* b < 2**-500 */
		if (hb <= (0x000fffff))
		{								/* subnormal b or 0 */
			__uint32_t low;

			GET_LOW_WORD(low, b);
			if ((hb | low) == 0)
				return a;
			t1 = 0;
			SET_HIGH_WORD(t1, 0x7fd00000);	/* t1=2^1022 */
			b *= t1;
			a *= t1;
			k -= 1022;
		} else
		{								/* scale a and b by 2^600 */
			ha += (0x25800000);		/* a *= 2^600 */
			hb += (0x25800000);		/* b *= 2^600 */
			k -= 600;
			SET_HIGH_WORD(a, ha);
			SET_HIGH_WORD(b, hb);
		}
	}
	
	/* medium size a and b */
	w = a - b;
	if (w > b)
	{
		t1 = Posit64{0};
		SET_HIGH_WORD(t1, ha);
		t2 = a - t1;
		w = Posit_sqrt(t1 * t1 - (b * (-b) - t2 * (a + t1)));
	} else
	{
		a = a + a;
		y1 = Posit64{0};
		SET_HIGH_WORD(y1, hb);
		y2 = b - y1;
		t1 = Posit64{0};
		SET_HIGH_WORD(t1, ha + (0x00100000));
		t2 = a - t1;
		w = Posit_sqrt(t1 * y1 - (w * (-w) - (t1 * y2 + t2 * b)));
	}

	return w;
}
#endif

/*
#include <limits>
#include <cmath>

int main(void)
{
	double nan = std::nan("0");
	double inf = std::numeric_limits<double>::infinity();
	double max = 1e9;
	double min = 1e-9;
	
	cout << fixed << setprecision(16) << "nan: " << Posit_asinPi(nan) << "\n";
	cout << fixed << setprecision(16) << "inf: " << Posit_asinPi(inf) << "\n";
	cout << fixed << setprecision(16) << "max: " << Posit_asinPi(max) << "\n";
	cout << fixed << setprecision(16) << "min: " << Posit_asinPi(min) << "\n";
}
*/

//int main(void)
//{
//	Posit64 a{-13}, b{14}, c{0};
//	c = Posit_hypot(a, b);
//	cout << "c = " << c << endl;
//}
