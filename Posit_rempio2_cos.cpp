#include <iostream>
#include "myfdlibm.h"
using namespace std;

#ifdef __STDC__
	__int32_t Posit_rempio2_cos(Posit64 x, Posit64 *y)
#else
	__int32_t Posit_rempio2_cos(x)
	Posit64 x, y[];
#endif
{   
	Posit64 d_pi = 2 * PI, zero{0}, two{2};
	Posit64 pi_over2 = PI / two;
	__int32_t n = 0;
	
	
	Quire64 qx,qpi_o2, qd_pi, qpi, qz{0}, qpi_o22;
	qx	= Quire_pToQ(x);
	qpi_o2  = Quire_pToQ(pi_over2);
	qd_pi   = Quire_pToQ(d_pi);
	qpi     = Quire_pToQ(PI);
	qpi_o22 = Quire_pToQ(pi_over2 / two);

	// become circle
/*
	while(x >= d_pi)
	{
		x = x - d_pi;
	}

	while(x < zero)
	{
		x = x + d_pi;
	}
*/

//======Quire become circle=====================================
	while(Quire_qToP(qx) >= Quire_qToP(qd_pi))
	{
		qx = Quire_qSubQ(qx, qd_pi);
	}
	
	while(qx < qz)
	{
		qx = Quire_qAddQ(qx, qd_pi);
	}
//==============================================================

/*
	if(x >= PI)
		x = d_pi - x;
*/

//==============================================================
	if((qx > qpi) || (Quire_qToP(qx) == Quire_qToP(qpi)))
		qx = Quire_qSubQ(qd_pi, qx);
//==============================================================

/*
	// a
	if(pi_over2 / 2 > x && x >= 0)
		n = 0;
	else if(pi_over2 > x && x >= pi_over2 / 2)
		n = 1;
	else if((pi_over2 + pi_over2 / 2) > x && x >= pi_over2)
		n = 2;
	else
		n = 3;
*/

//==============================================================
	if(qpi_o22 > qx && ((qx > qz) || (Quire_qToP(qx) == Quire_qToP(qz))))
		n = 0;
	else if(qpi_o2 > qx && ((qx > qpi_o22) || (Quire_qToP(qx) == Quire_qToP(qpi_o22))))
		n = 1;
	else if(Quire_qAddQ(qpi_o2, qpi_o22) > qx && ((qx > qpi_o2) || (Quire_qToP(qx) == Quire_qToP(qpi_o2))))
		n = 2;
	else
		n = 3;
//==============================================================



	y[0] = Quire_qToP(qx);
	y[1] = 0;

	return n;
}

// int main(){
//     Posit64 value{-123};
//     Posit64 y[2] = {0, 0};
//     cout<< Posit_rempio2(value, y)<<endl;
//     cout<< y[0] <<endl;
// }
