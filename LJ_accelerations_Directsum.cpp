#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* force (acceleration) due to the LJ (6-12) potential
on each body caused by all other bodies. 
Note that, should care about acc+= or acc=*/                           
void LJ_accelerations_Directsum( int iprint )
{
	int i, j;
	double xij, yij, zij, rij,rij2, rij4, rij8, rij14;
	double c2, c4, c6, c12, delta_ij;
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
			delta_ij=r[i]+r[j]-c_lj;
		if(rij2<(delta_ij+factor_lj*c_lj)*(delta_ij+factor_lj*c_lj))
		{
			rij=sqrt(rij2)-delta_ij;
			rij2=rij*rij;
			rij4 = rij2*rij2;
			rij8 = rij4*rij4;
			rij14 = rij2*rij4*rij8;
			
			c2=c_lj*c_lj;
			c4=c2*c2;
			c6=c4*c2;
			c12=c6*c6;
			//epsion_LJ is set to be equal to KbT !!
			factor = 4.0 *KbT * ( 12.0*c12/rij14 - 6.0*c6/rij8);

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
	} 

	for ( i=0 ; i<N ; i++ )
	{
		acc_x[i] = Fx[i] / mass[i];
		acc_y[i] = Fy[i] / mass[i];
		acc_z[i] = Fz[i] / mass[i];
	}

	if ( iprint == 1 )
		output_force();
}
