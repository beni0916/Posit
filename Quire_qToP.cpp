#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
    Posit64  Quire_qToP(Quire64 x)
#else
    Posit64  Quire_qToP(x)
	Quire64 x;
#endif
{   
    Posit64 y;
    y = x.to_value();
    return y;
}
#endif /* _DOUBLE_IS_32BITS */

/*int main(){
    Quire64 x{-8.26524};
    Posit64 y;
    y = qToP(x);
    cout << x << y << "\n";
}*/
