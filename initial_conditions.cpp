#include"MDpara.h"
//#include"ran.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*MD initializing steps:
1. initialize x, y,z randomly inside the spherical shell;
2. initialize vx, vy, vz as gaussian. 
3. adjust the velocity according to the zero total momentum constrain.*/ 
    
int Read_source_fromfile();  //one can either read initial config from file 
void Random_initial_config(); //or generate random initial config 
              
void initial_conditions( )
{
	int i,j,k,l;
	double g1, g2,g3, dumb,speed_scale, r1, r2, r3, mass_scale;
	double sq3, v;
	int dbgprint=0;
	

	double vec[3];
	double tempx, tempy, tempz,tempq,rij;
	double dmin; //depends on the ion and colloid size.

	//initialize radius array r[N]
	for(i=0;i<N;i++)
	{
	    if(i<N_ion)
		r[i]=r_ion;
	    else if(i<N_ion+N_col1)
		r[i]=r_col1;
	    else
	      r[i]=r_col2;
	}

	//initialize the cell-list parameters.
    if(Rc_lj1>Rc_lj2)
	Cell_num_1D=int(2*Rshell/Rc_lj1); 
//The number of cells in each dimension.
    else
        Cell_num_1D=int(2*Rshell/Rc_lj2);
    
	Cell_num_3D=Cell_num_1D*Cell_num_1D*Cell_num_1D; 
//The total number of cells in 3D.

	celli=double(Cell_num_1D);
	head= new int[Cell_num_1D*Cell_num_1D*Cell_num_1D]; 
//cell-list head array

	list= new int[N]; //cell-list particle list array.

	if(read_config==1)
	Read_source_fromfile();
	else if(read_config==0)
	Random_initial_config(); 
	
	// random velocities
	speed_scale = sqrt(KbT);

	for ( i=0 ; i<N ; i++ )
	{
	        mass_scale = 1.0/sqrt(mass[i]);
		GAUSS( &g1, &g2 );
		GAUSS( &g3, &dumb );
		r1 = g1;
		r2 = g2;
		r3 = g3;
		vx[i] = speed_scale * r1 * mass_scale;
		vy[i] = speed_scale * r2 * mass_scale;
		vz[i] = speed_scale * r3 * mass_scale;
	}

	// recale velocity to obtain 0 total momentum
	moment_x = 0.0;
	moment_y = 0.0;
	moment_z = 0.0;
	tot_mass = 0.0;
	for ( i=0 ; i<N; i++ )
	{
		tot_mass = tot_mass + mass[i];
		moment_x += mass[i] * vx[i];
		moment_y += mass[i] * vy[i];
		moment_z += mass[i] * vz[i];
	}
	for ( i=0 ; i<N ; i++ )
	{
		vx[i]= vx[i] - moment_x/tot_mass;
		vy[i]= vy[i] - moment_y/tot_mass;
		vz[i]= vz[i] - moment_z/tot_mass;

	}

	//obtaining the initial acc array
	//LJ_accelerations_Directsum(dbgprint );
	//Coulomb_accelerations_Directsum( dbgprint );
        Cell_list(dbgprint);
	Coulomb_accelerations_Hybrid( dbgprint );
	LJ_Boundaryforce (dbgprint);
	//Central_sph_force(dbgprint);
}

int Read_source_fromfile()
{
	int i,ind;
	double sizeratio;
	double *tempx, *tempy, *tempz;
	
	tempx=new double[N];
	tempy=new double[N];
	tempz=new double[N];

	FILE *jj;
	if((jj=fopen(config_name,"r"))==NULL)
	{
		printf("can not find the input configuration\
 file named %s, exiting the code!\n",config_name);
		exit (0);
	}
	i=0;
	while(!feof(jj))//make sure the input file format is correct.
	{ 
		fscanf(jj,"%d%lf%lf%lf",&ind,&tempx[i],&tempy[i],&tempz[i]); 
		i++;
	}
	fclose(jj);
	
	for (i=0;i<N;i++)
	{
		x[i]=tempx[i];
		y[i]=tempy[i];
		z[i]=tempz[i];
		if(i<N_ion1)
		  q[i]=q_ion1;
		else if(i<N_ion)
		  q[i]=q_ion2;
		else if(i<N_ion+N_col1)
		  q[i]=q_col1;
		else
		  q[i]=q_col2;
	}
	return 0;
}

void Random_initial_config()
{
/*random source charge locations 
(constrain 1: inside the box. 2: pairwise distance >dmin.)*/
	int i,j,k,l;
	double g1, g2,g3, dumb,speed_scale, r1, r2, r3;
	double sq3, v;
	int dbgprint=0;
	double vec[3];
	double tempx, tempy, tempz,tempq,rij;
	double dmin; //should depends on the ion and colloid size.
	double relax_factor=1.2; 
//the value 1.2 is empirical chosen. usually works if not too closely compacted
	
	for(j=0;j<N;j++)
	{
aaa:
		for(k=0;k<3;k++)
		{
			vec[k]=myran.doub();
		}

		x[j]=(vec[0]*2.0-1.0)*1.0*Rshell;

		y[j]=(vec[1]*2.0-1.0)*1.0*Rshell;

		z[j]=(vec[2]*2.0-1.0)*1.0*Rshell;


		tempx=x[j];
		tempy=y[j];
		tempz=z[j];

		rij=sqrt(tempx*tempx+tempy*tempy+tempz*tempz);
		dmin=Rshell-relax_factor*r_col;
		if(rij>dmin)
			goto aaa;

		for(l=0;l<j;l++)
		{
			rij=sqrt((tempx-x[l])*(tempx-x[l])+\
			(tempy-y[l])*(tempy-y[l])+(tempz-z[l])*(tempz-z[l]));

			dmin=relax_factor*(r[j]+r[l]);
			if(rij<dmin)
				goto aaa;
		}


		if(j<N_ion1)
			q[j]=q_ion1; //valences
        	else if (j<N_ion)
           		q[j]=q_ion2;
        	else if(j<N_ion+N_col1)
         	   	q[j]=q_col1;
       		else
            		q[j]=q_col2;
	}
	
}



