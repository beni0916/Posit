#include "myfdlibm.h"
#include "t_exp2.h"
#include <iostream>
#include <iomanip>
using namespace std;
// g++ -std=c++17 -o Posit_exp2 Posit_exp2.cpp Posit_floor.cpp Posit_fabs.cpp
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
	Posit64 Posit_exp2Minus1(Posit64 x)
#else
	Posit64 Posit_exp2Minus1(x)
	Posit64 x;
#endif
{
	if(x == NAR)
		return NAR;
	if(N != 2)
		return exp2(x) - 1;

	static Posit64 himark = DBL_MAX_EXP;
	static Posit64 lomark = DBL_MIN_EXP - DBL_MANT_DIG - 1;

	// if(P_BIT == 32)
	// {
	// 	__uint32_t a, b = 0;
	// 	GET_LOW_WORD(a, x);
	// 	cout << bitset<32>(a) << endl;
	// 	INSERT_WORDS(x, a, b);
	// }
    	
	if(x < himark)
	{
		__int32_t tval, unsafe;
		__uint32_t exp;
		Posit64 rx, x22, result, ex, t, x1;
		Posit64 ex2_u, scale_u, fabs_x = Posit_fabs(x);
		
		if(x < lomark)
			return 0;
		
		// cout << "51 end" << endl;
		// cout << "fabs_x = " << fabs(x) << endl;

		if(fabs_x < DBL_EPSILON / 4.0)
		{
			Posit64 ans;
			ans = Posit64{1.0} + x;
			return ans;
		}
		
		// cout << "60 end" << endl; 

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
			//cout << "ex = " << ex << endl;
			//cout << "t  = " << t  << endl;
			//cout << "x1 = " << x1 << endl;
			//cout << "rx = " << ex + t / 512 << endl;
			//cout << "tval = " << tval << endl;
			//cout << "tval & 511 " << (tval & 511) << endl;
			
			/* 2. Adjust for accurate table entry. */
			x -= Posit64{exp2_deltatable[tval & 511]};
			
			/* 3. Compute ex2 = 2^(t/512+e+ex).  */
			ex2_u = Posit64{exp2_accuratetable[tval & 511]};
			
	//	cout << "ex2_u = " << ex2_u << endl;
			
			tval >>= 9;
			unsafe = abs(tval) >= (double)-DBL_MIN_EXP - 56;

			if(P_BIT == 64)
				GET_HIGH_WORD(exp, ex2_u);
			else
				GET_LOW_WORD(exp, ex2_u);
	//	 cout << "exp = " << bitset<32>(exp) << endl;

			double tem = (double)ex2_u;
			
			//cout << "3. tval = " << tval << endl;
			//cout << "unsafe = " << unsafe << endl;
			//cout << "tem = " << tem << endl;
			//cout << "after get high word exp = " << exp << endl;
			
			int data[5] = {0};
			getRegime(exp, data);
			getExponent(exp, data);
			int offset = 12 - (N + 1 + data[1]);
	//	cout << "data[0] = " << data[0] << endl;
	//	cout << "data[1] = " << data[1] << endl;

			//cout << "offset = " << offset << endl;
			//exp = ((exp & UC(0x7ff00000)) + ((tval >> unsafe) << IEEE754_DOUBLE_SHIFT)) | (exp & ~UC(0x7ff00000));
			//exp = ((exp & ((~(0x7fffffff >> (N + data[1]))) & 0x7fffffff)) + ((tval >> unsafe) << (20 + offset))) | (exp & ~((~(0x7fffffff >> (N + data[1]))) & 0x7fffffff));
			//cout << "ori exp = " << bitset<32>(exp) << endl;
			//cout << "exp = " << exp << endl;
			
			// 1.
			int pow_diff = data[0] * pow(2, N) + data[2]; // the total power
			int re = 0, ex = 0;	// (r, e)
	//	cout << "pow_diff = " << pow_diff << endl;
	//	cout << "tval = " << tval << endl;

			// 2.
			re = (tval + pow_diff) / pow(2, N);
	//	cout << "re = " << re << endl;
			ex = (tval + pow_diff) % ((int)pow(2, N));
			if(ex < 0)
			{
				re = re - 1;
				ex = ex + 2*2*N;
			}		
	//	cout << "re = " << re << endl;
			// 3.
			int base = 0;
			if(re >= 0)
			{
				base = N;
				for(int i = 0; i < re; i++)
					base = base | (1 << (i + N));
				base <<= N; 	// 110 000
				base += ex;	// 110 001
			}
			else
			{
				base = 0 + ex;
				base = base | (1 << N);
			}
			
			//cout << "re >= 0 " << re << endl;
			//cout << "3. base = " << bitset<32>(base) << endl;
			
			// 4. caculate shift by regime bit
			int ori_bit = data[0];
			int new_bit = re;
			
	//	cout << "ori_bit = " << ori_bit << endl;
	//	cout << "new_bit = " << new_bit << endl;
			
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
			
	//	cout << "shift bit = " << shift_bit << endl;
			
			__uint32_t exp_l;
			GET_LOW_WORD(exp_l, ex2_u);
		
			if(shift_bit > 0)
			{
				exp_l >>= shift_bit;
				exp_l = (exp_l & (0xffffffff >> shift_bit)) | (0xfffffff & t_num) << (32 - shift_bit);
			}
	//	cout << "1 exp = " << bitset<32>(exp) << endl;
		 	exp = exp >> shift_bit;
		 	if(N != 1)
				base = base << (32 - (1 + new_bit + N));
	//	cout << "base = " << bitset<32>(base) << endl;
	//	cout << "2 exp = " << bitset<32>(exp) << endl;
	 		exp = (exp & (0xffffffff >> (1 + new_bit + N))) | base;
	//	cout << "3 exp = " << bitset<32>(exp) << endl;
		
			__uint32_t exp_t = 0;
	//	 cout << "exp   = " << bitset<32>(exp) << endl;
	//	 cout << "exp_l = " << bitset<32>(exp_l) << endl;
			if(P_BIT == 64)
				INSERT_WORDS(ex2_u, exp, exp_l);
			else
				INSERT_WORDS(ex2_u, exp, exp);
			// cout << "ex2_u = " << ex2_u << endl;
			
			__uint32_t hh, ll;
			EXTRACT_WORDS(hh, ll, ex2_u);
			// cout << "hh = " << bitset<32>(hh) << endl;
			// cout << "ll = " << bitset<32>(ll) << endl;

			getRegime(exp, data);
			getExponent(exp, data);
			exp = 0x3ff00000 + ((tval - (tval >> unsafe)) << (20 + offset));
			if(P_BIT == 64)
				INSERT_WORDS(scale_u, exp, 0);
			else
				INSERT_WORDS(scale_u, exp, exp);
			// cout << "scale_u = " << scale_u << endl;

			/* 4 */
			x22 = (((Posit64{.0096181293647031180} * x + Posit64{.055504110254308625}) * x + Posit64{.240226506959100583}) * x + Posit64{.69314718055994495}) * ex2_u;
		}
		
		// cout << "ex2_u = " << ex2_u << endl;
		result = x22 * x + ex2_u;
		// cout << "result = " << result << endl;

		if(P_BIT == 32)
		{
			__uint32_t hi_wd,lo_wd,  tmp = 0;
			GET_HIGH_WORD(hi_wd, result);
			GET_LOW_WORD(lo_wd, result);
			// cout << "hi_wd = " << hi_wd << endl;
			// cout << "lo_wd = " << lo_wd << endl;
			INSERT_WORDS(result, lo_wd, lo_wd);
		}
		
		// cout << "result = " << result << endl;

		if (!unsafe)
			return result - 1;
		
		return result * scale_u - 1;
	}
		
	return huge_val - 1;
}	

#endif

/*
#include <limits>
#include <cmath>

int main(void)
{
	double nan = std::nan("0");
	double inf = std::numeric_limits<double>::infinity();
	double max = 1e2;
	double min = 1e-15;
	
	cout << fixed << setprecision(16) << "nan: " << Posit_exp2(nan) << "\n";
	cout << fixed << setprecision(16) << "inf: " << Posit_exp2(inf) << "\n";
	cout << fixed << setprecision(16) << "max: " << Posit_exp2(max) << "\n";
	cout << fixed << setprecision(16) << "min: " << Posit_exp2(min) << "\n";
}
*/

/*
 int main(void)
 {
 	cout << "-4" << endl;
 	Posit64 num = Posit64{-4};
 	//Posit64 ans = Posit_exp2(num);
 	//cout<< "ans = " << fixed << setprecision(21)<<ans<<"\n";

 	cout << "2.2" << endl;
 	Posit64 num1 = Posit64{2.2};
 	Posit64 ans = Posit_exp2(num1);
 	cout<< "ans = " << fixed << setprecision(21)<<ans<<"\n";
}
*/
