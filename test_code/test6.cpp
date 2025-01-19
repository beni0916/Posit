// test6.cpp
#include <iostream>
#include "myfdlibm.h"

using namespace std;

int main(void)
{
	Posit64 test_values[] = 
	{
		Posit64(-1.0),
		Posit64(0.0),
		Posit64(0.5),
		Posit64(0.540302),
		Posit64(0.707107),
		Posit64(0.866025), 
		Posit64(0.877583),
		Posit64(1.0),		
	};
	
	for(auto x: test_values)
	{
		Posit64 res = Posit_acos(x);
		cout << "acos(" << x << ") = " << res << endl;
	}
	
	return 0;
}

