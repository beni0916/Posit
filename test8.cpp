// test8.cpp
#include <iostream>
#include "myfdlibm.h"

using namespace std;

int main(void)
{
	Posit64 values[] = 
	{
		Posit64(NAR),
		Posit64(0.0),
		Posit64(0.809797),
		Posit64(1.31696),
		Posit64(1.76275)
	};
	
	for(auto x: values)
	{
		Posit64 res = Posit_cosh(x);
		cout << "cosh(" << x << ") = " << res << endl;
	}
	
	return 0;
}
