#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// output routine for position and velocity
void output_x_v()
{
	int i;
	for ( i=0 ; i<N ; i++ )
	{
		fprintf( stderr, "Position & velocity:  %d %f %f %f %f %f %f\n", 
			i, x[i], y[i], z[i], vx[i], vy[i], vz[i] ); 
	}
}
