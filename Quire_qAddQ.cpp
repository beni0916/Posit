#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Quire64 Quire_qAddQ(Quire64 x, Quire64 y)
#else
    Quire64 Quire_qAddQ(x, y)
    Quire64 x;
    Quire64 y;
#endif
{   
    return x+y;
}

#endif /* _DOUBLE_IS_32BITS */

/*int main(){
    Quire64 x{-34552.26524};
    Posit64 y{34552}
    cout << Quire_qAddQ(x,y) << x << "\n";
}*/