#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//#include<omp.h>

/* force (acceleration) due to pairwise Coulomb force
on each body caused by all other bodies. FMM is used to achieve O(N) (N large)*/              

extern "C" void lfmm3dpartself_(int *,int *,int *,double [][3], int *,\
complex [],int *,complex[],double [][3],int *,complex [],int *,complex[][3]);
             
void Coulomb_accelerations_FMM( int iprint )
{
	int i, j;			

	for ( i=0 ; i<N; i++ )
	{
		Fx[i] = 0.0;
		Fy[i] = 0.0;
		Fz[i] = 0.0;
		
		source[i][0]=x[i];
		source[i][1]=y[i];
		source[i][2]=z[i];
		charge[i].real=q[i];
		charge[i].imag=0.0;
	}
	int ier,iprec=3,ifcharge=1,ifdipole=0,ifpot=0,iffld=1;
	lfmm3dpartself_(&ier,&iprec,&N,source,&ifcharge,charge,&ifdipole,\
dipstr,dipvec,&ifpot,pot,&iffld,fld);
	
	for ( i=0 ; i<N ; i++ )
	{
		acc_x[i] += q[i]*fld[i][0].real/ mass[i];
		acc_y[i] += q[i]*fld[i][1].real/ mass[i];
		acc_z[i] += q[i]*fld[i][2].real/ mass[i];
	}

	if ( iprint == 1 )
		output_force();
}
