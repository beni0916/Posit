#include "myfdlibm.h"
#include "t_exp2.h"
#include <iostream>
using namespace std;

#ifndef _DOUBLE_IS_32BITS

#ifdef __STDC__
static Posit64
#else
static Posit64
#endif
DBL_MAX_EXP = 1024,
DBL_MIN_EXP = -1021,
DBL_MANT_DIG = 53,
DBL_EPSILON = 2.2204460492503131e-016,
huge_val = 1.0e+300;

#ifdef __STDC__
	Posit64 Posit_exp2(Posit64 x)
#else
	Posit64 Posit_exp2(x)
	Posit64 x;
#endif
{
    	static Posit64 himark = DBL_MAX_EXP;
    	static Posit64 lomark = DBL_MIN_EXP - DBL_MANT_DIG - 1;
    	
    	if(x < himark)
    	{
    		__int32_t tval, unsafe;
    		__uint32_t exp;
		Posit64 rx, x22, result, ex, t, x1;
		Posit64 ex2_u, scale_u, fabs_x = Posit_fabs(x);
		
		if(x < lomark)
			return 0;
		
		if(fabs_x < DBL_EPSILON / 4.0)
		{
			Posit64 ans;
			ans = Posit64{1.0} + x;
			return ans;
		}
		
		{
			/* 1. Argument reduction. */
			// x = ex + t / 512 + x1
			ex = Posit_floor(x); 
			rx = x - ex;
			t = Posit_floor(rx * Posit64{512.0});	
			x1 = rx - Posit64{(double)t} / Posit64{512.0};	// x1 = x - ex - t / 512
			tval = (__int32_t)ex * 512 + (__int32_t)t + 256;	// tval = ex * 512 + t + 256
			rx = rx + ex - x1;	// rx = ex + t / 512
			x = x1;		// x = x1
			
			/* 2. Adjust for accurate table entry. */
			x -= Posit64{exp2_deltatable[tval & 511]};
			
			/* 3. Compute ex2 = 2^(t/512+e+ex).  */
			ex2_u = Posit64{exp2_accuratetable[tval & 511]};
			
			tval >>= 9;
			unsafe = abs(tval) >= (double)-DBL_MIN_EXP - 56;
			GET_HIGH_WORD(exp, ex2_u);
			double tem = (double)ex2_u;
			
			int data[5] = {0};
			getRegime(exp, data);
			getExponent(exp, data);
			int offset = 12 - (N + 1 + data[1]);
			
			// 1.
			int pow_diff = data[0] * 8 + data[2]; // the total power
			int re = 0, ex = 0;	// (r, e)
			
			// 2.
			re = (tval + pow_diff) / 8;
			ex = (tval + pow_diff) % 8;
			if(ex < 0)
			{
				re = re - 1;
				ex = ex + 8;
			}		
			
			// 3.
			int base = 0;
			if(re >= 0)
			{
				base = 2;
				for(int i = 0; i < re; i++)
					base = base | (1 << (i + 2));
				base <<= N; 	// 110 000
				base += ex;	// 110 001
			}
			else
			{
				base = 0 + ex;
				base = base | (1 << N);
			}
	
			// 4. caculate shift by regime bit
			int ori_bit = data[0];
			int new_bit = re;
			
			if(ori_bit >= 0)
		 		ori_bit = ori_bit + 2;
			else
			 	ori_bit = abs(ori_bit) + 1;
			 
			if(new_bit >= 0)
			 	new_bit = new_bit + 2;
			else
			 	new_bit = abs(new_bit) + 1;
			 	
			int shift_bit = new_bit - ori_bit;
			int t_num = (0xffffffff >> (32 - shift_bit)) & exp;
			__uint32_t exp_l;
			GET_LOW_WORD(exp_l, ex2_u);
			if(shift_bit > 0)
			{
				exp_l >>= shift_bit;
				exp_l = (exp_l & (0xffffffff >> shift_bit)) | (0xfffffff & t_num) << (32 - shift_bit);
			}
		 	
		 	exp = exp >> shift_bit;
		 	base = base << (32 - (1 + new_bit + N));
	 		exp = (exp & (0xffffffff >> (1 + new_bit + N))) | base;
			INSERT_WORDS(ex2_u, exp, exp_l);
			
			getRegime(exp, data);
			getExponent(exp, data);
			exp = 0x3ff00000 + ((tval - (tval >> unsafe)) << (20 + offset));
			INSERT_WORDS(scale_u, exp, 0);
			
			/* 4 */
			x22 = (((Posit64{.0096181293647031180} * x + Posit64{.055504110254308625}) * x + Posit64{.240226506959100583}) * x + Posit64{.69314718055994495}) * ex2_u;
		}
		
		result = x22 * x + ex2_u;

		if (!unsafe)
			return result;
		return result * scale_u;
	}
		
	return huge_val;
}	

#endif

//int main(void)
//{
//	Posit64 num = Posit64{-0.00001};
//	cout<< fixed << setprecision(21)<<Posit_exp2(num)<<"\n";
//}

