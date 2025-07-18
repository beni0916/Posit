#include "myfdlibm.h"
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
Posit64	one{1.0}, tiny{1.0e-300};
#else
Posit64 one{1.0}, tiny{1.0e-300};
#endif

#ifdef __STDC__
	Posit64 Posit_sqrt(Posit64 x)
#else
	Posit64 Posit_sqrt(x)
	Posit64 x;
#endif
{
	return sqrt(x);
	static int T1[32]= {
        0,	1024,	3062,	5746,	9193,	13348,	18162,	23592,
        29598,	36145,	43202,	50740,	58733,	67158,	75992,	85215,
        83599,	71378,	60428,	50647,	41945,	34246,	27478,	21581,
        16499,	12183,	8588,	5674,	3403,	1742,	661,	130,};

	Posit64 y, temp; 
	__uint32_t ix1;
	__int32_t ix0, k, y0;

    EXTRACT_WORDS(ix0,ix1,x);						//取得x的高位與低位

    k  = (ix0>>1) + 0x1ff80000;						//用作計算初始平方根的估計
	y0 = k - T1[31&(k>>15)];						//使用查表 (T1) 進一步修正估計值
    INSERT_WORDS(y, y0, 0);							//將結果放回變數y的高位，0放入y的低位

	
	y = (y+x/y) / 2;								//開始牛頓迭代

	do{
		temp = y;
		y = (y+x/y) / 2;
	}while(y - temp <= -1);							//直到本次跟前一次迭代的結果誤差小於1

	y = (y+x/y) / 2;
	y = (y+x/y) / 2;
	y = y - (y-x/y) / 2;							//微調
	
	return y;
}
 
#endif /* defined(_DOUBLE_IS_32BITS) */

// int main(){
// 	double x = 314.322;
//  	Posit64 px{x};                  
// 	cout << fixed << setprecision(20)<< Posit_sqrt(px) << "\n";
// 	cout << fixed << setprecision(20)<< sqrt(x) << "\n";
//  }

/*
Other methods  (use floating-point arithmetic)
-------------
(This is a copy of a drafted paper by Prof W. Kahan 
and K.C. Ng, written in May, 1986)

	Two algorithms are given here to implement sqrt(x) 
	(IEEE double precision arithmetic) in software.
	Both supply sqrt(x) correctly rounded. The first algorithm (in
	Section A) uses newton iterations and involves four divisions.
	The second one uses reciproot iterations to avoid division, but
	requires more multiplications. Both algorithms need the ability
	to chop results of arithmetic operations instead of round them, 
	and the INEXACT flag to indicate when an arithmetic operation
	is executed exactly with no roundoff error, all part of the 
	standard (IEEE 754-1985). The ability to perform shift, add,
	subtract and logical AND operations upon 32-bit words is needed
	too, though not part of the standard.

A.  sqrt(x) by Newton Iteration

   (1)	Initial approximation

	Let x0 and x1 be the leading and the trailing 32-bit words of
	a floating point number x (in IEEE double format) respectively 

	    1    11		     52				  ...widths
	   ------------------------------------------------------
	x: |s|	  e     |	      f				|
	   ------------------------------------------------------
	      msb    lsb  msb				      lsb ...order

 
	     ------------------------  	     ------------------------
	x0:  |s|   e    |    f1     |	 x1: |          f2           |
	     ------------------------  	     ------------------------

	By performing shifts and subtracts on x0 and x1 (both regarded
	as integers), we obtain an 8-bit approximation of sqrt(x) as
	follows.

		k  := (x0>>1) + 0x1ff80000;
		y0 := k - T1[31&(k>>15)].	... y ~ sqrt(x) to 8 bits
	Here k is a 32-bit integer and T1[] is an integer array containing
	correction terms. Now magically the floating value of y (y's
	leading 32-bit word is y0, the value of its trailing word is 0)
	approximates sqrt(x) to almost 8-bit.

	Value of T1:
	static int T1[32]= {
	0,	1024,	3062,	5746,	9193,	13348,	18162,	23592,
	29598,	36145,	43202,	50740,	58733,	67158,	75992,	85215,
	83599,	71378,	60428,	50647,	41945,	34246,	27478,	21581,
	16499,	12183,	8588,	5674,	3403,	1742,	661,	130,};

    (2)	Iterative refinement

	Apply Heron's rule three times to y, we have y approximates 
	sqrt(x) to within 1 ulp (Unit in the Last Place):

		y := (y+x/y)/2		... almost 17 sig. bits
		y := (y+x/y)/2		... almost 35 sig. bits
		y := y-(y-x/y)/2	... within 1 ulp


	Remark 1.
	    Another way to improve y to within 1 ulp is:

		y := (y+x/y)		... almost 17 sig. bits to 2*sqrt(x)
		y := y - 0x00100006	... almost 18 sig. bits to sqrt(x)

				2
			    (x-y )*y
		y := y + 2* ----------	...within 1 ulp
			       2
			     3y  + x


	This formula has one division fewer than the one above; however,
	it requires more multiplications and additions. Also x must be
	scaled in advance to avoid spurious overflow in evaluating the
	expression 3y*y+x. Hence it is not recommended uless division
	is slow. If division is very slow, then one should use the 
	reciproot algorithm given in section B.

    (3) Final adjustment

	By twiddling y's last bit it is possible to force y to be 
	correctly rounded according to the prevailing rounding mode
	as follows. Let r and i be copies of the rounding mode and
	inexact flag before entering the square root program. Also we
	use the expression y+-ulp for the next representable floating
	numbers (up and down) of y. Note that y+-ulp = either fixed
	point y+-1, or multiply y by nextafter(1,+-inf) in chopped
	mode.

		I := FALSE;	... reset INEXACT flag I
		R := RZ;	... set rounding mode to round-toward-zero
		z := x/y;	... chopped quotient, possibly inexact
		If(not I) then {	... if the quotient is exact
		    if(z=y) {
		        I := i;	 ... restore inexact flag
		        R := r;  ... restore rounded mode
		        return sqrt(x):=y.
		    } else {
			z := z - ulp;	... special rounding
		    }
		}
		i := TRUE;		... sqrt(x) is inexact
		If (r=RN) then z=z+ulp	... rounded-to-nearest
		If (r=RP) then {	... round-toward-+inf
		    y = y+ulp; z=z+ulp;
		}
		y := y+z;		... chopped sum
		y0:=y0-0x00100000;	... y := y/2 is correctly rounded.
	        I := i;	 		... restore inexact flag
	        R := r;  		... restore rounded mode
	        return sqrt(x):=y.
		    
    (4)	Special cases

	Square root of +inf, +-0, or NaN is itself;
	Square root of a negative number is NaN with invalid signal.


B.  sqrt(x) by Reciproot Iteration

   (1)	Initial approximation

	Let x0 and x1 be the leading and the trailing 32-bit words of
	a floating point number x (in IEEE double format) respectively
	(see section A). By performing shifs and subtracts on x0 and y0,
	we obtain a 7.8-bit approximation of 1/sqrt(x) as follows.

	    k := 0x5fe80000 - (x0>>1);
	    y0:= k - T2[63&(k>>14)].	... y ~ 1/sqrt(x) to 7.8 bits

	Here k is a 32-bit integer and T2[] is an integer array 
	containing correction terms. Now magically the floating
	value of y (y's leading 32-bit word is y0, the value of
	its trailing word y1 is set to zero) approximates 1/sqrt(x)
	to almost 7.8-bit.

	Value of T2:
	static int T2[64]= {
	0x1500,	0x2ef8,	0x4d67,	0x6b02,	0x87be,	0xa395,	0xbe7a,	0xd866,
	0xf14a,	0x1091b,0x11fcd,0x13552,0x14999,0x15c98,0x16e34,0x17e5f,
	0x18d03,0x19a01,0x1a545,0x1ae8a,0x1b5c4,0x1bb01,0x1bfde,0x1c28d,
	0x1c2de,0x1c0db,0x1ba73,0x1b11c,0x1a4b5,0x1953d,0x18266,0x16be0,
	0x1683e,0x179d8,0x18a4d,0x19992,0x1a789,0x1b445,0x1bf61,0x1c989,
	0x1d16d,0x1d77b,0x1dddf,0x1e2ad,0x1e5bf,0x1e6e8,0x1e654,0x1e3cd,
	0x1df2a,0x1d635,0x1cb16,0x1be2c,0x1ae4e,0x19bde,0x1868e,0x16e2e,
	0x1527f,0x1334a,0x11051,0xe951,	0xbe01,	0x8e0d,	0x5924,	0x1edd,};

    (2)	Iterative refinement

	Apply Reciproot iteration three times to y and multiply the
	result by x to get an approximation z that matches sqrt(x)
	to about 1 ulp. To be exact, we will have 
		-1ulp < sqrt(x)-z<1.0625ulp.
	
	... set rounding mode to Round-to-nearest
	   y := y*(1.5-0.5*x*y*y)	... almost 15 sig. bits to 1/sqrt(x)
	   y := y*((1.5-2^-30)+0.5*x*y*y)... about 29 sig. bits to 1/sqrt(x)
	... special arrangement for better accuracy
	   z := x*y			... 29 bits to sqrt(x), with z*y<1
	   z := z + 0.5*z*(1-z*y)	... about 1 ulp to sqrt(x)

	Remark 2. The constant 1.5-2^-30 is chosen to bias the error so that
	(a) the term z*y in the final iteration is always less than 1; 
	(b) the error in the final result is biased upward so that
		-1 ulp < sqrt(x) - z < 1.0625 ulp
	    instead of |sqrt(x)-z|<1.03125ulp.

    (3)	Final adjustment

	By twiddling y's last bit it is possible to force y to be 
	correctly rounded according to the prevailing rounding mode
	as follows. Let r and i be copies of the rounding mode and
	inexact flag before entering the square root program. Also we
	use the expression y+-ulp for the next representable floating
	numbers (up and down) of y. Note that y+-ulp = either fixed
	point y+-1, or multiply y by nextafter(1,+-inf) in chopped
	mode.

	R := RZ;		... set rounding mode to round-toward-zero
	switch(r) {
	    case RN:		... round-to-nearest
	       if(x<= z*(z-ulp)...chopped) z = z - ulp; else
	       if(x<= z*(z+ulp)...chopped) z = z; else z = z+ulp;
	       break;
	    case RZ:case RM:	... round-to-zero or round-to--inf
	       R:=RP;		... reset rounding mod to round-to-+inf
	       if(x<z*z ... rounded up) z = z - ulp; else
	       if(x>=(z+ulp)*(z+ulp) ...rounded up) z = z+ulp;
	       break;
	    case RP:		... round-to-+inf
	       if(x>(z+ulp)*(z+ulp)...chopped) z = z+2*ulp; else
	       if(x>z*z ...chopped) z = z+ulp;
	       break;
	}

	Remark 3. The above comparisons can be done in fixed point. For
	example, to compare x and w=z*z chopped, it suffices to compare
	x1 and w1 (the trailing parts of x and w), regarding them as
	two's complement integers.

	...Is z an exact square root?
	To determine whether z is an exact square root of x, let z1 be the
	trailing part of z, and also let x0 and x1 be the leading and
	trailing parts of x.

	If ((z1&0x03ffffff)!=0)	... not exact if trailing 26 bits of z!=0
	    I := 1;		... Raise Inexact flag: z is not exact
	else {
	    j := 1 - [(x0>>20)&1]	... j = logb(x) mod 2
	    k := z1 >> 26;		... get z's 25-th and 26-th 
					    fraction bits
	    I := i or (k&j) or ((k&(j+j+1))!=(x1&3));
	}
	R:= r		... restore rounded mode
	return sqrt(x):=z.

	If multiplication is cheaper then the foregoing red tape, the 
	Inexact flag can be evaluated by

	    I := i;
	    I := (z*z!=x) or I.

	Note that z*z can overwrite I; this value must be sensed if it is 
	True.

	Remark 4. If z*z = x exactly, then bit 25 to bit 0 of z1 must be
	zero.

		    --------------------
		z1: |        f2        | 
		    --------------------
		bit 31		   bit 0

	Further more, bit 27 and 26 of z1, bit 0 and 1 of x1, and the odd
	or even of logb(x) have the following relations:

	-------------------------------------------------
	bit 27,26 of z1		bit 1,0 of x1	logb(x)
	-------------------------------------------------
	00			00		odd and even
	01			01		even
	10			10		odd
	10			00		even
	11			01		even
	-------------------------------------------------

    (4)	Special cases (see (4) of Section A).	
 
 */
