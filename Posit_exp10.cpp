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
	Posit64 Posit_exp10(Posit64 x)
#else
	Posit64 Posit_exp10(x)
	Posit64 x;
#endif
{
	return Posit_exp(MLN10 * x);
}
#endif
