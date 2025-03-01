#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Posit64 Posit_sign(Posit64 x)
#else
	Posit64 Posit_sign(x)
	Posit64 x;
#endif
{   
    Posit64 y{0.0};
    if (Posit_compareLess(x, y)) {
        return Posit64{-1.0};
    }
    else if (Posit_compareGreater(x, y)) {
        return Posit64{1.0};
    }
    else {
        return Posit64{0.0};
    }
}

#endif /* _DOUBLE_IS_32BITS */

/*int main(){
    Posit64 x{-34552.26524};
    cout << x << Posit_sign(x) << "\n";
}*/