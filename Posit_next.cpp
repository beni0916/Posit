#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Posit64 Posit_next(Posit64 x)
#else
	Posit64 Posit_next(x)
	Posit64 x;
#endif
{   
    Posit64 y;
    __uint32_t hx, lx;
    GET_HIGH_WORD(hx, x);
    GET_LOW_WORD(lx, x);
    lx += 1;
    if (lx==0) {
        hx += 1;
    }
    SET_HIGH_WORD(y, hx);
    SET_LOW_WORD(y, lx);
    return y;
}

#endif /* _DOUBLE_IS_32BITS */

/*int main(){
    Posit64 x{-34552.26524};
    __uint32_t hx, lx;
    GET_HIGH_WORD(hx, x);
    GET_LOW_WORD(lx, x);
    cout << hx << lx  << "\n";
    Posit64 y = Posit_next(x);
    GET_HIGH_WORD(hx, y);
    GET_LOW_WORD(lx, y);
    cout << hx << lx  << "\n";
}*/