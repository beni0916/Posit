#include <iostream>
#include "myfdlibm.h"
using namespace std;

#ifdef __STDC__
	__int32_t Posit_remhalf_hp(Posit64 x, Posit64 *y)
#else
	__int32_t Posit_remhalf_hp(x)
	Posit64 x, y[];
#endif
{
    static const Posit64 PI_high = 3.14159265358979323846264338327950288419716939937510;
    static const Posit64 PI_lo = 1.224646799147353207e-16;
    static const Posit64 half = 0.5;

    __int32_t n = (int)(x * 2); /* 快速計算 n = x / 0.5 */
    Posit64 high_correction = n * half;
    Posit64 low_correction = n * PI_lo;

    Posit64 nx = x - high_correction; /* 減去高位補償 */
    nx = nx - low_correction; /* 減去低位補償 */

    y[0] = nx * PI_high; /* 高精度結果 */
    y[1] = nx * PI_lo;   /* 保留低位 */
    return n;
}

