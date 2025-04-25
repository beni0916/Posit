#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Posit64 Posit_ceil(Posit64 x)
#else
	Posit64 Posit_ceil(x)
	Posit64 x;
#endif
{
    Posit64 y = Posit_floor(x);
    if (x!=y) {
        y = y+1;
    }
    return y;
}

#endif 

/*int main(){
    Posit64 x{-34552.26524};
	Posit64 y{34552.26524};
    cout << x << y << Posit_ceil(x) << Posit_ceil(y) << "\n";
}*/
