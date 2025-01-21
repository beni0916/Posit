// test11.cpp
#include <iostream>
#include "myfdlibm.h"

using namespace std;

int main(void)
{
	Posit64 values[] = 
	{
		Posit64(-0.16666666),
		Posit64(-0.25),
		Posit64(-0.33333333),
		Posit64(-0.5),
		Posit64(0.0),
		Posit64(0.16666666),
		Posit64(0.25),
		Posit64(0.33333333),
		Posit64(0.5),
		Posit64(0.66666666),
		Posit64(0.75),
		Posit64(0.83333333),
		Posit64(1.0),
	};
	
	for(auto x: values)
	{
		Posit64 res = Posit_sinPi(x);
		cout << "sinPi(" << x << ") = " << res << endl;
	}
	
	return 0;
}
