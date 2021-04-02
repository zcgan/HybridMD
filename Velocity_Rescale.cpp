#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// rescale velocity on x, y, z direactions due to the termprature constrain.

void Velocity_Rescale( )
{
	int i;

	double factor_x, factor_y, factor_z;

	double kinetic_x=0.;

	double kinetic_y=0.;

	double kinetic_z=0.;


	for ( i=0 ; i<N; i++ )
	{
		kinetic_x = kinetic_x + 0.5 * mass[i] * vx[i]*vx[i];
		kinetic_y = kinetic_y + 0.5 * mass[i] * vy[i]*vy[i];
		kinetic_z = kinetic_z + 0.5 * mass[i] * vz[i]*vz[i];			
	}

	// velocity rescale.
	factor_x = sqrt(kinetic_x/(0.5*double(N)*KbT));
	factor_y = sqrt(kinetic_y/(0.5*double(N)*KbT));
	factor_z = sqrt(kinetic_z/(0.5*double(N)*KbT));

	for ( i=0 ; i<N ; i++ )
	{
		vx[i]= vx[i]/factor_x;
		vy[i]= vy[i]/factor_y;
		vz[i]= vz[i]/factor_z;
	}

}
