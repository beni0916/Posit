#include <iostream>
#include "myfdlibm.h"

using namespace std;

int main(void)
{
	Posit64 values[] = 
	{
		Posit64(0.0),
		Posit64(1.0),
		Posit64(1.3462),
		Posit64(2.0),
		Posit64(3.0)
	};
	
	for(auto x: values)
	{
		Posit64 res = Posit_acosh(x);
		cout << "acosh(" << x << ") = " << res << endl;
	}
	
	return 0;
}

