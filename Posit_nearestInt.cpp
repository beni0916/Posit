#include "myfdlibm.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
	Posit64 Posit_nearestInt(Posit64 x)
#else
	Posit64 Posit_nearestInt(x)
	Posit64 x;
#endif
{   
    Posit64 floor = Posit_floor(x);
    Posit64 ceil = Posit_ceil(x);
    Posit64 floorD = x - floor;
    Posit64 ceilD = ceil - x;

    if (floorD > ceilD) {
        return ceil;
    }
    else if (floorD < ceilD) {
        return floor;
    }
    else {
        if (floorD == ceilD) {
            return floor;
        }
        else {
            if (((int)floor)%2==0) {
                return floor;
            }
            else {
                return ceil;
            }
        }
    }
}

#endif 

/*int main(){
    Posit64 x{-34552.26524};
	Posit64 y{-34552.76524};
    cout << x << y << Posit_nearestInt(x) << Posit_nearestInt(y) << "\n";
}*/
