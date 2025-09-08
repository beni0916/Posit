#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Posit64 Posit_floor(Posit64 x)
#else
	Posit64 Posit_floor(x)
	Posit64 x;
#endif
{	
	// if (isnan(x)) return NAR;
	Posit64 y = (int)x;
	if (x!=y) {
		if (y<Posit64{0.0}) {
			y = y-1;
		}
	}
	return y;
}

#endif /* _DOUBLE_IS_32BITS */

/*int main(){
    Posit64 x{-34552.26524};
	Posit64 y{34552.26524};
    cout << x << y << Posit_floor(x) << Posit_floor(y) << "\n";
}*/

