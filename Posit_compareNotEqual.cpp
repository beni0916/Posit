#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
    bool Posit_compareNotEqual(Posit64 x, Posit64 y)
#else
    bool Posit_compareNotEqual(x, y)
	Posit64 x;
    Posit64 y;
#endif
{   
    if (Posit_compareEqual(x, y)) {
        return false;
    }
    else {
        return true;
    }
}

#endif /* _DOUBLE_IS_32BITS */

/*int main(){
    Posit64 x{-34552.26524};
    Posit64 y{-34552.26524};
    Posit64 z{-34552.26520};
    cout << x << y << z << Posit_compareNotEqual(x,y) << Posit_compareNotEqual(x,z) << "\n";
}*/