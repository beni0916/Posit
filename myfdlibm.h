#include "universal/number/posit/posit.hpp"
#ifndef N
    #define N 2
#endif
using Posit64 = sw::universal::posit<64, N>;

typedef union {
	Posit64 value{0};
	struct {
		__uint32_t lsw;
		__uint32_t msw;
	} parts;
} posit64_type;

#define GET_HIGH_WORD(i, d)					\
do {								\
    posit64_type gh_u;					\
    gh_u.value = (d);						\
    (i) = gh_u.parts.msw;						\
} while (0)

#define GET_LOW_WORD(i,d)					\
do {								\
    posit64_type gl_u;					\
    gl_u.value = (d);						\
    (i) = gl_u.parts.lsw;						\
} while (0)

#define SET_HIGH_WORD(d,v)					\
do {								\
    posit64_type sh_u;					\
    sh_u.value = (d);						\
    sh_u.parts.msw = (v);						\
    (d) = sh_u.value;						\
} while (0)

#define SET_LOW_WORD(d,v)					\
do {								\
    posit64_type sl_u;					\
    sl_u.value = (d);						\
    sl_u.parts.lsw = (v);						\
    (d) = sl_u.value;						\
} while (0)

#define EXTRACT_WORDS(ix0,ix1,d)				\
do {								\
    posit64_type ew_u;					\
    ew_u.value = (d);						\
    (ix0) = ew_u.parts.msw;					\
    (ix1) = ew_u.parts.lsw;					\
} while (0)

#define INSERT_WORDS(d,ix0,ix1)					\
do {								\
    posit64_type iw_u;					\
    iw_u.parts.msw = (ix0);					\
    iw_u.parts.lsw = (ix1);					\
    (d) = iw_u.value;						\
} while (0)

extern Posit64 Posit_sqrt (Posit64 x);
extern Posit64 Posit_log(Posit64 x);
extern Posit64 Posit_log10(Posit64 x);
extern Posit64 Posit_floor (Posit64 x);
extern Posit64 Posit_ceil (Posit64 x);
extern Posit64 Posit_fabs(Posit64 x);
extern Posit64 Posit_exp(Posit64 x);
extern Posit64 Posit_pow(Posit64 x, Posit64 y);
extern Posit64 __kernel_sin(Posit64 x, Posit64 y, int iy);
extern Posit64 __kernel_cos(Posit64 x, Posit64 y);
extern Posit64 __kernel_tan(Posit64 x, Posit64 y, int iy);
extern __int32_t Posit_rempio2(Posit64 x, Posit64 *y);
extern Posit64 Posit_sin(Posit64 x);
extern Posit64 Posit_cos(Posit64 x);
extern Posit64 Posit_tan(Posit64 x);
extern Posit64 Posit_sinh(Posit64 x);
extern Posit64 Posit_cosh(Posit64 x);
extern Posit64 Posit_tanh(Posit64 x);
extern Posit64 Posit_asin(Posit64 x);
extern Posit64 Posit_acos(Posit64 x);
extern Posit64 Posit_atan(Posit64 x);
extern Posit64 Posit_asinh(Posit64 x);
extern Posit64 Posit_acosh(Posit64 x);
extern Posit64 Posit_atanh(Posit64 x);

inline const Posit64 PI(3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233);
inline const Posit64 PI2(9.8696044010893586188344909998762);
extern Posit64 Posit_cosPi(Posit64 x);
extern Posit64 Posit_sinPi(Posit64 x);
extern Posit64 Posit_tanPi(Posit64 x);
extern __int32_t Posit_remhalf(Posit64 x, Posit64 *y);
extern Posit64 Posit_asinPi(Posit64 x);
extern Posit64 Posit_acosPi(Posit64 x);
extern Posit64 Posit_atanPi(Posit64 x);
extern Posit64 Posit_expMinus1(Posit64 x);
extern Posit64 Posit_exp2(Posit64 x);
extern Posit64 Posit_exp2Minus1(Posit64 x);
extern Posit64 Posit_exp10(Posit64 x);
extern Posit64 Posit_exp10Minus1(Posit64 x);
extern Posit64 Posit_rSqrt(Posit64 x);
extern Posit64 Posit_logPlus1(Posit64 x);
extern Posit64 Posit_log2(Posit64 x);
extern Posit64 Posit_log2Plus1(Posit64 x);
extern Posit64 Posit_log10Plus1(Posit64 x);
extern Posit64 Posit_negate(Posit64 x);
extern Posit64 Posit_abs(Posit64 x);
extern Posit64 Posit_sign(Posit64 x);
extern Posit64 Posit_nearestInt(Posit64 x);
extern Posit64 Posit_next(Posit64 x);
extern Posit64 Posit_prior(Posit64 x);
extern bool Posit_compareEqual(Posit64 x, Posit64 y);
extern bool Posit_compareNotEqual(Posit64 x, Posit64 y);
extern bool Posit_compareGreater(Posit64 x, Posit64 y);
extern bool Posit_compareGreaterEqual(Posit64 x, Posit64 y);
extern bool Posit_compareLess(Posit64 x, Posit64 y);
extern bool Posit_compareLessEqual(Posit64 x, Posit64 y);
extern Posit64 Posit_hypot(Posit64 x, Posit64 y);
extern Posit64 Posit_atan2(Posit64 y, Posit64 x);
extern Posit64 Posit_atan2Pi(Posit64 y, Posit64 x);


inline void getRegime(__int32_t ix, int *data){
    int r = 0, k = 0;                   //r紀錄regime bits代表的數值，k紀錄有幾個regime bits
    ix <<= 1;                           //清除sign bit
    if(ix > 0){                         //若regime bits第一個位元為0
        while(ix > 0){                  //此迴圈讓ix一直左移，直到ix < 0，代表遇到相反的位元(1)，即可停止
            k++;                        //regime bits數加1
            ix <<= 1;                   //左移ix
        }
        r = -k;                         //根據公式，若第一個位元為0，r = -k
    }
    else{                               //若regime bits第一個位元為1
        while(ix < 0){                  //此迴圈讓ix一直左移，直到ix > 0，代表遇到相反的位元(0)，即可停止
            k++;                        
            ix <<= 1;                   
        }
        r = k - 1;                      //根據公式，若第一個位元為0，r = k - 1
    }
    data[0] = r;                        
    data[1] = k + 1;
}

inline void getExponent(__int32_t ix, int *data){
    ix <<= data[1];                        //右移regime bits數量的位數 -> (0 or 1) (exponent bits) (fraction bits)
    ix &= 0x7fffffff;                      //若首位為1，之後右移會補1，所以先將首位設成0
    ix >>= (32 - N - 1);                   //最後右移正確位元數
    data[2] = ix;                          //得到正確的exponent
}

