#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
static Posit64
#else
static Posit64
#endif
MLN10 = 2.3025850929940456840179914546843642076011014886287729760333279;

#ifdef __STDC__
	Posit64 Posit_exp10Minus1(Posit64 x)
#else
	Posit64 Posit_exp10Minus1(x)
	Posit64 x;
#endif
{
	Posit64 tmp = Posit_exp(MLN10 * x), one{1.0};
	Posit64 res = tmp - one;
	return res;
}
#endif
