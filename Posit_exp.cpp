#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
static Posit64
#else
static Posit64
#endif
invln2 = 1.4426950408889634;

#ifdef __STDC__
	Posit64 Posit_exp(Posit64 x)	/* default IEEE double exp */
#else
	Posit64 Posit_exp(x)	/* default IEEE double exp */
	Posit64 x;
#endif
{
	return Posit_exp2(x * invln2);
}

#endif /* defined(_DOUBLE_IS_32BITS) */

/*
int main()
{
	double x = 0.969413, i_res = 0;
	Posit64 y{x}, p_res{0};
	p_res = Posit_exp(y);
	i_res = exp(x);
	cout << "ieee  :" << bitset<64>(i_res) <<"\n";
	cout << "exp(x) : " << exp(x) << endl;
}
*/
