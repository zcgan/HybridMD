#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// Velocity Verlet
void Velocity_Verlet_Langevin( int iprint )
{
	int i, j;
	double g1,g2,g3,na; //4 gauss random variables. the 4th one is not used.

	// Velocity Verlet
	for ( i=0 ; i<N ; i++ )
	{
		//random force.
		GAUSS(&g1,&g2);
		GAUSS(&g3,&na);
		
		Fran_x[i]=noise[i]*g1;
		Fran_y[i]=noise[i]*g2;
		Fran_z[i]=noise[i]*g3;
		
		vx[i]=vx[i]+0.5*(acc_x[i]+Fran_x[i]/mass[i]-FricCoef*vx[i])*dt;
		vy[i]=vy[i]+0.5*(acc_y[i]+Fran_y[i]/mass[i]-FricCoef*vy[i])*dt;
		vz[i]=vz[i]+0.5*(acc_z[i]+Fran_z[i]/mass[i]-FricCoef*vz[i])*dt;

		x[i]=x[i]+vx[i]*dt;
		y[i]=y[i]+vy[i]*dt;
		z[i]=z[i]+vz[i]*dt;
	}

		// accelerations at t+dt

#ifdef DIRECT_SUM		
		LJ_accelerations_Directsum( iprint);
#else
		Cell_list(iprint);
#endif
		//Coulomb_accelerations_Directsum( iprint );
		//Coulomb_accelerations_FMM( iprint );
		Coulomb_accelerations_Hybrid( iprint );
		LJ_Boundaryforce( iprint );
	//	Central_sph_force(iprint);

		for ( i=0 ; i<N ; i++ )
		{
		vx[i]=vx[i]+0.5*(acc_x[i]+Fran_x[i]/mass[i]-FricCoef*vx[i])*dt;
		vy[i]=vy[i]+0.5*(acc_y[i]+Fran_y[i]/mass[i]-FricCoef*vy[i])*dt;
		vz[i]=vz[i]+0.5*(acc_z[i]+Fran_z[i]/mass[i]-FricCoef*vz[i])*dt;
		}

	
}
