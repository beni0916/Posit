#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Posit64 Posit_abs(Posit64 x)
#else
	Posit64 Posit_abs(x)
	Posit64 x;
#endif
{
    Posit64 y{0.0};
	if(x<y) x *= -1;
    return x;
}
#endif /* _DOUBLE_IS_32BITS */

/*int main(){
    Posit64 x{-34552.26524};
    cout << x << Posit_abs(x) << "\n";
}*/