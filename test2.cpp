// test2.cpp
#include <iostream>
#include "myfdlibm.h"

using namespace std;

int main(void)
{
	Posit64 test_values[] = 
	{
		Posit64(0.0),
		Posit64(0.5),
		Posit64(-0.5),
		Posit64(1.0),
		Posit64(-1.0),
		Posit64(3.14159265 / 6), 	// 30 degree
		Posit64(3.14159265 / 4),	// 45 degree
		Posit64(3.14159265 / 3), 	// 60 degree
		Posit64(3.14159265 / 2)	// 90 degree
	};
	
	for(auto x: test_values)
	{
		Posit64 result = __kernel_tan(x, Posit64(0.0), 1);
		cout << "tan(" << x << ") = " << result << endl;
	} 
	
	return 0;
}
