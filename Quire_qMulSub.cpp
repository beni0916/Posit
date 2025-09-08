#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Quire64 Quire_qMulSub(Quire64 x, Posit64 y, Posit64 z)
#else
    Quire64 Quire_qMulSub(x, y, z)
    Quire64 x;
    Posit64 y;
    Posit64 z;
#endif
{
    x-=quire_mul(y,z);
    return x;
}

#endif /* _DOUBLE_IS_32BITS */

/*int main(){
    Quire64 x{-8.26524};
    Posit64 y{-2}
    Posit64 z{4}
    cout << Quire_qMulSub(x,y,z) << x << "\n";
}*/