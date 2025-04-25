#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
    bool Posit_compareGreater(Posit64 x, Posit64 y)
#else
    bool Posit_compareGreater(x, y)
	Posit64 x;
    Posit64 y;
#endif
{   
    if (isnan(x)||isnan(y)) return false;
    if (x>y) {
        return true;
    }
    else {
        return false;
    }
}

#endif /* _DOUBLE_IS_32BITS */

/*int main(){
    Posit64 x{-34552.26524};
    Posit64 y{-34552.26524};
    Posit64 z{-34552.26520};
    cout << x << y << z << Posit_compareGreater(x,y) << Posit_compareGreater(x,z) << "\n";
}*/