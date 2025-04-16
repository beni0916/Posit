#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
    Quire64  Quire_pToQ(Posit64 x)
#else
    Quire64  Quire_pToQ(x)
	Posit64 x;
#endif
{
    Quire64 y = x;
    return y;
}
#endif /* _DOUBLE_IS_32BITS */

/*int main(){
    Posit64 x{-8.26524};
    Quire64 y;
    y = pToQ(x);
    cout << x << y << "\n";
}*/

