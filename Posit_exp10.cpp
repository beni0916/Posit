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

/*
#include <limits>
#include <cmath>

int main(void)
{
	double nan = std::nan("0");
	double inf = std::numeric_limits<double>::infinity();
	double max = 1e2;
	double min = 1e-15;
	
	cout << fixed << setprecision(16) << "nan: " << Posit_exp10(nan) << "\n";
	cout << fixed << setprecision(16) << "inf: " << Posit_exp10(inf) << "\n";
	cout << fixed << setprecision(16) << "max: " << Posit_exp10(max) << "\n";
	cout << fixed << setprecision(16) << "min: " << Posit_exp10(min) << "\n";
}
*/

