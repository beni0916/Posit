// test12.cpp
#include <iostream>
#include "myfdlibm.h"

using namespace std;

int main(void)
{
	Posit64 pi{3.141592025756836};

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
		Posit64 yP[2] = {Posit64(0.0), Posit64(0.0)};
		Posit64 nP = Posit_rempio2(x * pi, yP);
		
		Posit64 yH[2] = {Posit64(0.0), Posit64(0.0)};
		Posit64 nH = Posit_remhalf(x, yH);
		
		cout << "x = " << x << endl;
		cout << "nP = " << nP << "yP[0] = " << yP[0] << endl;
		cout << "nH = " << nH << "yH[0] = " << yH[0] << endl;
	}
	
	return 0;
}
