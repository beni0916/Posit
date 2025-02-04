// test3.cpp
#include <iostream>
#include "myfdlibm.h"

using namespace std;

int main(void)
{
	Posit64 test_values[] = 
	{
		Posit64(123.0),
		Posit64(-123.0)
	};
	
	for(auto x: test_values)
	{
		Posit64 y[2] = {0, 0};
		Posit64 n = Posit_rempio2(x, y);
		cout << "Posit_rempio2(" << x << ") = " << y[0] << ", n = " << n << endl;
	}
	
	return 0;
}
