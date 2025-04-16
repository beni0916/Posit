#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Posit64 Posit_Nfloor(Posit64 x)
#else
	Posit64 Posit_Nfloor(x)
	Posit64 x;
#endif
{
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
    cout << x << y << Posit_Nfloor(x) << Posit_Nfloor(y) << "\n";
}*/

