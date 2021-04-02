#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// Velocity Verlet
void Velocity_Verlet( int iprint )
{
	int i, j;

	// Velocity Verlet
	for ( i=0 ; i<N ; i++ )
	{
		
		vx[i]=vx[i]+0.5*(acc_x[i])*dt;
		vy[i]=vy[i]+0.5*(acc_y[i])*dt;
		vz[i]=vz[i]+0.5*(acc_z[i])*dt;

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
		vx[i]=vx[i]+0.5*(acc_x[i])*dt;
		vy[i]=vy[i]+0.5*(acc_y[i])*dt;
		vz[i]=vz[i]+0.5*(acc_z[i])*dt;
		}

	
}
