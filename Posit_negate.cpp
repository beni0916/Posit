#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Posit64 Posit_negate(Posit64 x)
#else
	Posit64 Posit_negate(x)
	Posit64 x;
#endif
{
    return -1*x;
}

#endif /* _DOUBLE_IS_32BITS */

/*int main(){
    Posit64 x{-34552.26524};
    cout << Posit_negate(x) << "\n";
}*/