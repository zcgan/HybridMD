#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* force (acceleration) due to Coulomb potential
on each body caused by all other bodies. 
(Calculated using pairwise directsum)*/ 
                          
void Coulomb_accelerations_Directsum( int iprint )
{
	int i, j;
	double xij, yij, zij, rij, rij2, rij3;
	double Fijx, Fijy, Fijz;
	double factor;

	for ( i=0 ; i<N; i++ )
	{
		Fx[i] = 0.0;
		Fy[i] = 0.0;
		Fz[i] = 0.0;
	}

	for ( i=0 ; i<N-1 ; i++ )
	{
		for ( j= i+1 ; j<N ; j++ )
		{
			xij = x[i] - x[j];
			yij = y[i] - y[j];
			zij = z[i] - z[j];

			rij2 = xij*xij + yij*yij + zij*zij ;

			rij=sqrt(rij2);

			rij3=rij*rij2;

			factor = q[i]*q[j]/rij3; 

			Fijx = factor * xij;
			Fijy = factor * yij;
			Fijz = factor * zij;

			Fx[i] = Fx[i] + Fijx;
			Fy[i] = Fy[i] + Fijy;
			Fz[i] = Fz[i] + Fijz;
			Fx[j] = Fx[j] - Fijx;
			Fy[j] = Fy[j] - Fijy;
			Fz[j] = Fz[j] - Fijz;

		} 
	} 

	for ( i=0 ; i<N ; i++ )
	{
		acc_x[i] += Fx[i] / mass[i];
		acc_y[i] += Fy[i] / mass[i];
		acc_z[i] += Fz[i] / mass[i];
	}

	if ( iprint == 1 )
		output_force();
}
