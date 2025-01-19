// test5.cpp
#include <iostream>
#include "myfdlibm.h"

using namespace std;

int main(void)
{
/*
	Posit64 test_values[] = 
	{
		Posit64(0.0),
		Posit64(0.5),
		Posit64(-0.5),
		Posit64(1.0),
		Posit64(-1.0),
		Posit64(0.166666666666), 	// 30 degree
		Posit64(0.25),	// 45 degree
		Posit64(0.333333333333), 	// 60 degree
		Posit64(0.5)	// 90 degree
	};
	
	for(auto x: test_values)
	{
		Posit64 res = Posit_cosPi(x);
		cout << "cosPi(" << x << ") = " << res << endl;
	}
*/
	Posit64 n{0.25};
	cout << "Posit_cosPi" << endl;
	cout << fixed << setprecision(16) << "Result: " << Posit_cosPi(n) << "\n";
	
	Posit64 n1{0.5};
	cout << "Posit_cosPi" << endl;
	cout << fixed << setprecision(16) << "Result: " << Posit_cosPi(n1) << "\n";
	
	Posit64 value{0.785398006439209};
 	cout << "Posit_cos" << endl;
	cout << fixed << setprecision(16) << "Result: " << Posit_cos(value) << "\n";
	
	Posit64 v1{0.785398006439209 * 2.0};
 	cout << "Posit_cos" << endl;
	cout << fixed << setprecision(16) << "Result: " << Posit_cos(v1) << "\n";
	
	return 0;
}
