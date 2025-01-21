// test10.cpp
#include <iostream>
#include "myfdlibm.h"

using namespace std;

int main(void)
{
	Posit64 values[] = 
	{
		Posit64(0),
		Posit64(0.1001675),
		Posit64(-0.1001675),
		Posit64(1.17520119),
		Posit64(-1.17520119),
		Posit64(74.20321057),
		Posit64(-74.20321057),
		Posit64(4.8724e+09),
		Posit64(-4.8724e+09),
		Posit64(2.04587e+149),
		Posit64(-2.04587e+149)
	};
	
	for(auto x: values)
	{
		Posit64 res = Posit_asinh(x);
		cout << "asinh(" << x << ") = " << res << endl;
	}
	
	return 0;
}
