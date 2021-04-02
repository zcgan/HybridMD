#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// total energy & momentum 
// check conservation
void energy_momentum_check( int info, double time )
{
	double kinetic, potential, pot, etotal;

	double xij, yij, zij, rij,rij2, rij4, rij6, rij12;
	double dmin;
	double c_boundary,c2,c4,c6,c12,delta_ij;
	int i, j;

	kinetic = 0.0;

	for ( i=0 ; i<N; i++ )
	{
		kinetic = kinetic + 0.5 * mass[i] * 
			( vx[i]*vx[i] + vy[i]*vy[i] +  vz[i]*vz[i]);
	}
    
    	Cell_list(0);
	Coulomb_accelerations_Hybrid(0);
	
	//printf("ele=%f,lj=%f\n",printed_ele_energy,printed_lj_energy);

	potential =printed_ele_energy+printed_lj_energy;      

	// sum over pairs -no longer used
/*	for ( i=0 ; i<N; i++ )
	{
		xij = x[i] ;
		yij = y[i] ;
		zij = z[i] ;
		rij2 = xij*xij + yij*yij + zij*zij;
		dmin=factor_lj*r[i];
		
		for ( j=i+1 ; j<N ; j++ )
		{
			xij = x[i] - x[j];
			yij = y[i] - y[j];
			zij = z[i] - z[j];
			rij2 = xij*xij + yij*yij + zij*zij;
			rij=sqrt(rij2);
			pot=q[i]*q[j]/rij;  //the electrostatic potential
			potential = potential +pot;  
			delta_ij=r[i]+r[j]-c_lj;
			
			if(rij<(delta_ij+factor_lj*c_lj)) 
//if rij<rc, compute the shift-truncated lennard-jones potential.
			{
			rij=c_lj/(rij-delta_ij);
			rij2=rij*rij;	
			rij4 = rij2*rij2;
			rij6 = rij2*rij4;
			rij12 = rij6*rij6;
			pot = 4.0 * (rij12-rij6)+1.0;
			potential = potential + pot;
			
			}
		}
	}*/

	etotal = kinetic + potential;

	/*if ( info == 0 )
	{
		printf( "%f %f %f %f\n", time, kinetic, potential, etotal );
		return;
	}*/

	// center of mass momentum
	moment_x = 0.0;
	moment_y = 0.0;
	moment_z = 0.0;
	for ( i=0 ; i<N ; i++ )
	{
		moment_x = moment_x + mass[i] * vx[i];
		moment_y = moment_y + mass[i] * vy[i];
		moment_z = moment_z + mass[i] * vz[i];
	}
	moment_x = moment_x / tot_mass;
	moment_y = moment_y / tot_mass;
	moment_z = moment_z / tot_mass;

	if(info==1)
fprintf(output_energy, " %f %f %f %f %f %f %f\n",kinetic, printed_ele_energy,\
printed_lj_energy,kinetic*2/3/N,moment_x,moment_y,moment_z);

	if(info==0)
fprintf(output_energy_equilibration, "%f %f %f %f %f %f %f\n",kinetic, \
printed_ele_energy,printed_lj_energy,kinetic*2/3/N,moment_x,moment_y,moment_z);

}
