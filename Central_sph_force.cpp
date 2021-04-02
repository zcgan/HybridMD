#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double *im_int_x = new double[N*im];
double *im_int_y = new double[N*im];
double *im_int_z = new double[N*im];
double *im_int_q = new double[N*im];
double *wwwquad,*xxxquad;
//the Gauss-jacobi quadrature function
extern void cgqf ( int nt, int kind, double alpha, double beta, double a,
double b, double t[], double wts[] );

/* force (acceleration) due to the central spherical boundary. Note
that, should care about acc+= or acc=*/                           
void Central_sph_force( int iprint )
{
	int i, j;
	double  ri,ri2,ri3, rib,rib2, rib4, rib8, rib14;
	double factor;
	double c_boundary,c2,c4,c6,c12;
	double dmin; 
//dmin: the distance to the shell that the LJ force need to be evaluated.


	for ( i=0 ; i<N; i++ )
	{
		Fx[i] = 0.0;
		Fy[i] = 0.0;
		Fz[i] = 0.0;
	}

/*Generating the images inside the central spherical dielctric interface*/
/*Gan 10/23/2019: currently no longer used,i.e., no central sphere exists*/

	int imcount=0; //count the total number of image charges.
	int im=4;
	int ii, ind;
	double epsi_s=1.0;
	double sourcetol=4.0;
	double dx,dy,dz,dr,rk,rim;
	double start=0;
	double alpha=0.0;
	double end;
	double lamda=epsi_s/(epsi_s+epsi_M),beta=lamda-1.0,\
gamma=(epsi_M-epsi_s)/(epsi_s+epsi_M); 
//if we first assume epsi_M are all the same for different dielectrics.
	int kind=4;
	int order=im-1;

	if(iter_indicator==0)
			{
				wwwquad = new double[order];
				xxxquad = new double[order];
			}

	for (i=0;i<N;i++)
	{

		dx=x[i]-0.;
		dy=y[i]-0.;
		dz=z[i]-0.;
		dr=sqrt(dx*dx+dy*dy+dz*dz);

		if(dr/RM<sourcetol&&dr>RM)
		{						
			rk=RM*RM/dr;
			rim=rk/dr;
			im_int_x[imcount]=rim*dx+0.;
			im_int_y[imcount]=rim*dy+0.;
			im_int_z[imcount]=rim*dz+0.;
			im_int_q[imcount]=-gamma*RM*q[i]/dr; //the Kelvin image.

			//the (im-1) line image.
			end=rk;
			cgqf (order, kind, alpha, beta, start, end,\
			xxxquad, wwwquad);

			for(ii=0;ii<order;ii++)
			{
				ind=imcount+ii+1;
				im_int_x[ind]=0.+xxxquad[ii]*dx/dr;
				im_int_y[ind]=0.+xxxquad[ii]*dy/dr;
				im_int_z[ind]=0.+xxxquad[ii]*dz/dr;
				im_int_q[ind]=wwwquad[ii]*gamma*lamda\
				*pow(rk,1.0-lamda)*q[i]/RM;						
			}
			imcount=imcount+im;
		}					
	}

	for ( i=0 ; i<N ; i++ )
	{

		ri2 = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];

		/*Coulomb force*/
		ri=sqrt(ri2);
		ri3=ri*ri2;
		factor = q[i]*QM/ri3;
		Fx[i] = factor * x[i]; 
		Fy[i] = factor * y[i]; 
		Fz[i] = factor * z[i];

		/*LJ force*/
		dmin=factor_lj*r[i];
		if(ri2<(dmin+RM)*(dmin+RM))
		{

			c2=r[i]*r[i];
			c4=c2*c2;
			c6=c4*c2;
			c12=c6*c6;	
			rib=ri-RM;
			rib2=rib*rib;
			rib4 = rib2*rib2;
			rib8 = rib4*rib4;
			rib14 = rib2*rib4*rib8;

			factor = 4.0 * ( 12.0*c12/rib14 - 6.0*c6/rib8);
//the negative sign: the force is the inversion direction of the ri vector.
			Fx[i] = Fx[i]+1.0*factor * x[i]*rib/ri; 
			Fy[i] = Fy[i]+1.0*factor * y[i]*rib/ri; 
			Fz[i] = Fz[i]+1.0*factor * z[i]*rib/ri;
		}
	} 

	/*image-sph & ion interaction*/
	double temfx,temfy,temfz;	
	for(i=0;i<imcount;i++)
		for(ii=0;ii<N;ii++)
		{
			dx=-im_int_x[i]+x[ii];
			dy=-im_int_y[i]+y[ii];
			dz=-im_int_z[i]+z[ii];
			dr=sqrt(dx*dx+dy*dy+dz*dz);

			temfx=q[ii]*im_int_q[i]*dx/dr/dr/dr/epsi_s;
			temfy=q[ii]*im_int_q[i]*dy/dr/dr/dr/epsi_s;
			temfz=q[ii]*im_int_q[i]*dz/dr/dr/dr/epsi_s;

			Fx[i]+= temfx;
			Fy[i]+= temfy;
			Fz[i]+= temfz;
																																																						
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
