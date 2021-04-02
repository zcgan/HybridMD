#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// print acc & forces (debug)
void output_force()
{
	int i;
	double Ftotx, Ftoty, Ftotz, fx, fy, fz;

	Ftotx = 0.0;
	Ftoty = 0.0;
	Ftotz = 0.0;

	for ( i=0 ; i<N ; i++ )
	{
		fx =  mass[i]*acc_x[i];
		fy =  mass[i]*acc_y[i]; 
		fz =  mass[i]*acc_z[i];
		printf("Mass & acceleration:  %d %f %f %f %f -- %f %f %f\n", 
			i, mass[i], acc_x[i], acc_y[i], acc_z[i], fx, fy, fz );
		Ftotx = Ftotx + fx;
		Ftoty = Ftoty + fy;
		Ftotz = Ftotz + fz;
	}
	printf("Total force (x & y): %f %f %f\n", Ftotx, Ftoty, Ftotz );
}
