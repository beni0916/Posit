#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Quire64 Quire_qNegate(Quire64 x)
#else
    Quire64 Quire_qNegate(x)
    Quire64 x;
#endif
{
    x.set_sign(!x.sign());
    return x;
}

#endif /* _DOUBLE_IS_32BITS */

/*int main(){
    Quire64 x{-34552.26524};
    cout << Quire_qNegate(x) << "\n";
}*/