#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Quire64 Quire_qSubP(Quire64 x, Posit64 y)
#else
    Quire64 Quire_qSubP(x, y)
    Quire64 x;
    Posit64 y;
#endif
{
    x-=y;
    return x;
}

#endif /* _DOUBLE_IS_32BITS */

/*int main(){
    Quire64 x{-34552.26524};
    Posit64 y{-34552}
    cout << Quire_qSubP(x,y) <<  x << "\n";
}*/