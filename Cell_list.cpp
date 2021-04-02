#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* Cell-list accelerated force (acceleration) due to 6-12 potential
on each body caused by all other bodies. 

Note that one should care about acc+= or acc=*/                           
void Cell_list( int iprint )
{
	int i,j,o;
	printed_lj_energy=0.;
	/*initialize the head and list arrow*/
	for(i=0;i<Cell_num_1D*Cell_num_1D*Cell_num_1D;i++)
	{
		head[i]=-1;
	}
    
/*Bug1 found by ziwei: run into infnite loop because of bad configuration. 
Reason(gan note): due to ions/cols goes outside the boundary*/
    double dr;
    for (i=0;i<N;i++)
    {
        dr=0.0;
        dr=x[i]*x[i]+y[i]*y[i]+z[i]*z[i];
        if (dr>=Rshell*Rshell)
        {   printf("particles run outside the boundary,"
" please check your input file, and make the equilibration stage longer!");
            exit(0);
        }
    }

	for(i=0;i<N;i++)
	{
		j=floor((x[i]+Rshell)*celli/2/Rshell)\
		 +floor((y[i]+Rshell)*celli/2/Rshell)*Cell_num_1D\
		 +floor((z[i]+Rshell)*celli/2/Rshell)*Cell_num_1D*Cell_num_1D;
		list[i]=head[j];
		head[j]=i;
	}

	/*initialize the acc array*/
	for ( i=0 ; i<N ; i++ )
	{
		acc_x[i] = 0.;
		acc_y[i] = 0.;
		acc_z[i] = 0.;
	}

	/*Searching the 27 neighbors of the target particle and
 evaluate the force on it by calling LJForce. */

	for(i=0;i<N;i++)
	{
		o=floor((x[i]+Rshell)*celli/2/Rshell)\
+floor((y[i]+Rshell)*celli/2/Rshell)*Cell_num_1D\
+floor((z[i]+Rshell)*celli/2/Rshell)*Cell_num_1D*Cell_num_1D;

		if(head[o]>=0)
			LJForce(i,head[o]);

		o=o+1;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}

		o=o-2;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}


		o=o+1+Cell_num_1D;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}



		o=o-2*Cell_num_1D;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}



		o=o+1+2*Cell_num_1D;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}


		o=o-2;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}


		o=o+2-2*Cell_num_1D;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}



		o=o-2;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}



		o=o+1+Cell_num_1D+Cell_num_1D*Cell_num_1D;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}



		o=o+1;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}


		o=o-2;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}

		o=o+1-2*Cell_num_1D*Cell_num_1D;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}

		o=o+1;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}


		o=o-2;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}

		o=o+1+Cell_num_1D;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}

		o=o+1;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}


		o=o-2;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}

		o=o+1-2*Cell_num_1D+2*Cell_num_1D*Cell_num_1D;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}


		o=o+1;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}


		o=o-2;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}

		o=o+1-2*Cell_num_1D*Cell_num_1D;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}

		o=o+1;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}


		o=o-2;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}

		o=o+1+2*Cell_num_1D+2*Cell_num_1D*Cell_num_1D;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}

		o=o+1;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}

		o=o-2;
		if(o>=0&&o<Cell_num_3D)
		{
			if(head[o]>=0)
				LJForce(i,head[o]);
		}


	}

	if ( iprint == 1 )
		output_force();
}

void LJForce(int indx, int site)
{
	/*evaluating the force according to the "list"*/
	double xo,yo,zo;
	double xio, yio, zio,rio,rio2, rio4,rio6,rio12;
	double r0,r02;
	double c2,c4,c6,c12,delta_ij;
	double Fiox, Fioy, Fioz;
	double factor;
	int i;
	double pot;

	if(site!=indx)
	{
		xo=x[site];
		yo=y[site];
		zo=z[site];
		xio=x[indx]-xo;
		yio=y[indx]-yo;
		zio=z[indx]-zo;
		r02=xio*xio + yio*yio + zio*zio;
		r0=sqrt(r02);
		delta_ij=r[site]+r[indx]-c_lj;
		if (r02<=delta_ij*delta_ij)
		  {
		    printf("Two particles (hard core) are overlapping!"
" Please decrease your timestep!\n");
		    exit(0);
		  }
	else if (r02<(delta_ij+factor_lj*c_lj)*(delta_ij+factor_lj*c_lj))
		{
		  rio=r0-delta_ij;
		  rio2=rio*rio;
		  rio4 = rio2*rio2;
		  rio6 = rio2*rio4;
		  rio12 = rio6*rio6;
		  c2=c_lj*c_lj;
		  c4=c2*c2;
		  c6=c4*c2;
		  c12=c6*c6;

			if(force_compute==1)
			{
		// Epsilon_LJ is set to be equal to KbT!!
			 factor = 4.0 * KbT*(12.0*c12/rio12-6.0*c6/rio6)/rio/r0;
			 Fiox = factor * xio;
			 Fioy = factor * yio;
			 Fioz = factor * zio;

			  acc_x[indx] += Fiox / mass[indx];
			  acc_y[indx] += Fioy / mass[indx];
			  acc_z[indx] += Fioz / mass[indx];
			 }
			if(energy_compute==1)
			{
			// Epsilon_LJ is set to be equal to KbT!!
			// Factor 0.5 is account for double-counting
			  pot = 0.5*(4.0 * (c12/rio12-c6/rio6)+1.0)*KbT;
			  printed_lj_energy = printed_lj_energy+ pot;
			}		
		}
	}

	i=list[site];
	if(i>=0)
	{
		LJForce(indx,i);
	}
}
