#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// gaussian distributed random 
// numbers - Koonin & Meredith 1986

//the rand() function should be replaced later (*already replaced*).

void GAUSS( double *gauss1, double *gauss2 )
{
	double x1, x2, twou, radius, theta;

	x1 = myran.doub();
	x2 = myran.doub();
	twou = -2*log( 1.0-x1 );
	radius = sqrt( twou);
	theta = 2*3.1415926535897932*x2;
	*gauss1 = radius*cos(theta);
	*gauss2 = radius*sin(theta);
}
