#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Quire64 Quire_qAbs(Quire64 x)
#else
    Quire64 Quire_qAbs(x)
    Quire64 x;
#endif
{
    return abs(x);
}

#endif /* _DOUBLE_IS_32BITS */

/*int main(){
    Quire64 x{-34552.26524};
    cout << Quire_qAbs(x) << "\n";
}*/