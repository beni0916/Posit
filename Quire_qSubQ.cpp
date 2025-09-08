#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Quire64 Quire_qSubQ(Quire64 x, Quire64 y)
#else
    Quire64 Quire_qSubQ(x, y)
    Quire64 x;
    Quire64 y;
#endif
{
    x-=y;
    return x;
}

#endif /* _DOUBLE_IS_32BITS */

/*int main(){
    Quire64 x{-34552.26524};
    Quire64 y{-34552}
    cout << Quire_qSubQ(x,y) <<  x << "\n";
}*/