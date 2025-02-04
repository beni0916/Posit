// test9.cpp
#include <iostream>
#include "myfdlibm.h"

using namespace std;

int main(void)
{
	Posit64 values[] = 
	{
		Posit64(0.0),
		Posit64(0.1),
		Posit64(-0.1),
		Posit64(1.0),
		Posit64(-1.0),
		Posit64(5.0),
		Posit64(-5.0),
		Posit64(23.0),
		Posit64(-23.0),
		Posit64(710.0),
		Posit64(-710.0)
	};
	
	for(auto x: values)
	{
		Posit64 res = Posit_sinh(x);
		cout << "sinh(" << x << ") = " << res << endl;
	}
	
	return 0;
}
