#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* force (acceleration) due to the LJ equipped boundary. 
Note that, should care about acc+= or acc=*/                           
void LJ_Boundaryforce( int iprint )
{
	int i, j;
	double  ri,ri2, rib,rib2, rib4, rib8, rib14;
	double factor;
	double c_boundary,c2,c4,c6,c12;
	double dmin; 
//dmin:truncated distance 2 the boundary that the LJ force need to be evaluated.


	for ( i=0 ; i<N; i++ )
	{
		Fx[i] = 0.0;
		Fy[i] = 0.0;
		Fz[i] = 0.0;
	}

	for ( i=0 ; i<N ; i++ )
	{

		ri2 = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
		dmin=factor_lj*r[i];
		if(ri2>(Rshell-dmin)*(Rshell-dmin))
		{
			c_boundary=r[i];
			c2=c_boundary*c_boundary;
			c4=c2*c2;
			c6=c4*c2;
			c12=c6*c6;
			ri = sqrt(ri2);	
			rib=Rshell-ri;
			rib2=rib*rib;
			rib4 = rib2*rib2;
			rib8 = rib4*rib4;
			rib14 = rib2*rib4*rib8;
			// epsilon_LJ is set to be equal to KbT !!
			factor = 4.0 * KbT * ( 12.0*c12/rib14 - 6.0*c6/rib8);
			/*the negative sign represents 
			the force is the inversion direction of the ri vector.*/
			Fx[i] = -1.0*factor * x[i]*rib/ri;
			Fy[i] = -1.0*factor * y[i]*rib/ri; 
			Fz[i] = -1.0*factor * z[i]*rib/ri;
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
