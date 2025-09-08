#include "myfdlibm.h"
#include <iostream>

using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
static const Posit64
#else
static Posit64
#endif
tiny = 1.0e-300,
zero = 0.0,
pi_o_4 = 7.8539816339744827900E-01,	/* 0x3FE921FB, 0x54442D18 */
pi_o_2 = 1.5707963267948965580E+00,	/* 0x3FF921FB, 0x54442D18 */
pi = 3.1415926535897931160E+00,		/* 0x400921FB, 0x54442D18 */
pi_lo = 1.2246467991473531772E-16;	/* 0x3CA1A626, 0x33145C07 */

#ifdef __STDC__
	Posit64 Posit_atan2(Posit64 y, Posit64 x)
#else
	Posit64 Posit_atan2(y, x)
	Posit64 y, x;
#endif
{
	Posit64 z;
	__int32_t k, m, hx, hy, ix, iy;
	__uint32_t lx, ly;
	
	EXTRACT_WORDS(hx, lx, x);
	ix = hx & (0x7fffffff);
	EXTRACT_WORDS(hy, ly, y);
	iy = hy & (0x7fffffff);
	
	//if (((ix | ((lx | -lx) >> 31)) > (0x7ff00000)) || ((iy | ((ly | -ly) >> 31)) > (0x7ff00000)))	/* x or y is NaN */
	//	return x + y;
	//if (((hx - (0x3ff00000)) | lx) == 0)
	if(x == Posit64{1.0})
		return Posit_atan(y);					/* x=1.0 */
	
	m = ((hy >> 31) & 1) | ((hx >> 30) & 2);	/* 2*sign(x)+sign(y) */

	/* when y = 0 */
	if ((iy | ly) == 0)
	{
		switch ((int)m)
		{
		case 0:
		case 1:
			return y;					/* atan(+-0,+anything)=+-0 */
		case 2:
			return pi + tiny;			/* atan(+0,-anything) = pi */
		case 3:
			return -pi - tiny;			/* atan(-0,-anything) =-pi */
		}
	}
	/* when x = 0 */
	if ((ix | lx) == 0)
	{
		return (hy < 0) ? -pi_o_2 - tiny : pi_o_2 + tiny;
	}

	/* when x is INF */
	if (ix == (0x7ff00000))
	{
		if (iy == (0x7ff00000))
		{
			switch ((int)m)
			{
			case 0:
				return pi_o_4 + tiny;	/* atan(+INF,+INF) */
			case 1:
				return -pi_o_4 - tiny;	/* atan(-INF,+INF) */
			case 2:
				return 3.0 * pi_o_4 + tiny;	/*atan(+INF,-INF) */
			case 3:
				return -3.0 * pi_o_4 - tiny;	/*atan(-INF,-INF) */
			}
		} else
		{
			switch ((int)m)
			{
			case 0:
				return zero;			/* atan(+...,+INF) */
			case 1:
				return -zero;			/* atan(-...,+INF) */
			case 2:
				return pi + tiny;		/* atan(+...,-INF) */
			case 3:
				return -pi - tiny;		/* atan(-...,-INF) */
			}
		}
	}
	/* when y is INF */
	if (iy == (0x7ff00000))
	{
		return (hy < 0) ? -pi_o_2 - tiny : pi_o_2 + tiny;
	}

	/* compute y/x */
	z = Posit_atan( Posit_fabs(y / x));			/* safe to do y/x */
	
	switch ((int)m)
	{
	case 0:
		return z;						/* atan(+,+) */
	case 1:
		{
			uint32_t zh;

			GET_HIGH_WORD(zh, z);
			SET_HIGH_WORD(z, zh ^ (0x80000000));
		}
		return z;						/* atan(-,+) */
	case 2:
		return pi - (z - pi_lo);		/* atan(+,-) */
	default:							/* case 3 */
		return (z - pi_lo) - pi;		/* atan(-,-) */
	}
}

#endif

//int main(void)
//{
//	Posit64 y{0.383032}, x{0.156005}, z;
//	z = Posit_atan2(y, x);
//	cout << "z = " << z << endl;		
//}
