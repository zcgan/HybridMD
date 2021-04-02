#include "MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include<string.h>
void allocate_dynamic();
void set_precision(int prec);

int main( int argn, char * argv[] )
{
	int dbg_print, special, print_energy,i,j,k;
	int count_rescale=0;
	int count_sample=0;
	int count_steps=0;
   	int traj_ind=0;
	double t;
	int seed=1;
	Ran myran(seed);

//if there is inputfile name specified, use it;
//otherwise, try use the default name "para.txt"

	if(argv[1]!=NULL)
	{
		strcpy(paraname,argv[1]);
		printf("the input file name is: %s\n",paraname);
	}
	else
	{
		strcpy(paraname,"para.txt");
		printf("the input file name is: %s\n",paraname);
	}

	read_para(); //read the input parameters

	allocate_arrays(); //allocate the arrays

	dbg_print = 0; // extra information print for debug purpose if ==1.

	//record energy
	print_energy = 0;

	// recording trajectories in file
	if(Production_steps>=sample_cycle)
	{
		trajectory = fopen( "dump.lammpstrj", "w" );
		output_energy= fopen("energy_production", "w");
		fprintf(output_energy,"# time kinetic ele_pot"
		" LJ_pot temperature Px Py Pz\n");
	}
	
 	//record trajectory for the annealing process
	trajectory_annealing=fopen( "annealing_dump.lammpstrj", "w" );

	output_energy_equilibration= fopen("energy_equilibration", "w");
	output_finalconfig=fopen("finalconfig.txt","w");	
	fprintf(output_energy_equilibration,"# time kinetic"
	" ele_pot LJ_pot temperature Px Py Pz\n");
	
	//recording rdf in file
#ifdef RDF
	if(Production_steps>=sample_cycle)
	initialize_rdf();
#endif

	// mass
	mass_initialize( dbg_print );

	// time & noise setting	
	t = tmin;
	dt=dt_initial;
	for (i=0;i<N;i++)
	  {
	    noise[i]=sqrt(2.0*FricCoef*mass[i]*KbT/dt_final);
	  }
	
	// initial positions & velocities
	set_precision(prec_set);
    	force_compute=0;energy_compute=1;
	initial_conditions( );
	
	iter_indicator=1;
/*iter_indicator=1 means we can start use the result from last timestep as
the initial guess (usually a good guess) for the next time step.*/
    
   	 /*produce the initial output when t=0*/
    	fprintf(output_energy_equilibration, "%f ",t);
    	force_compute=0;energy_compute=1; //check energy, turn off force compute
    	energy_momentum_check( 0, t );
    	force_compute=1;energy_compute=0;//start MD, turn off energy compute

	/*****Starting the main simulation loops with 4 different stages*****/
	printf("#steps=%d, start running the simulation; dt=%f\n",\
	count_steps,dt_initial);

	printf("In the equilibration stage,"
	" we first set the precision to %d digits\n",prec_set);

	/*****Stage 1: the time and velocity rescale loop*****/
	for(i=0;i<time_stepratio+1;i++)
	{
		for(j=0;j<dt_inc_interval;j++)
		{
			count_rescale=count_rescale+1;
			// time step -- velocity-verlet
			Velocity_Verlet( dbg_print );

			t = t + dt;
			
		//rescale the velocity to maintain the temprature constrain.
			if(count_rescale%rescale_cycle==0)	
			{		
				Velocity_Rescale( );				
			}
			if(count_rescale%sample_cycle==0)	
			{		
				fprintf(output_energy_equilibration, "%f ",t);
               			printf("time=%f \n",t);
				force_compute=0;energy_compute=1;
				energy_momentum_check( 0, t );		
				force_compute=1;energy_compute=0;	
			}			
		}
			if(dt*2<=dt_final)
			{
				dt=dt*2;
				printf("#steps=%d, dt increased to %f\n",\
				count_rescale,dt);
			}
			else
			{
				dt=dt_final;
				if(i<time_stepratio)
				printf("#steps=%d, dt increased to %f\n",\
				count_rescale,dt);
			}
	}
		
		
	
	/*****Stage 2: the langevin themostat equilibration loop*****/
	printf("#steps=%d, now switching to Langevin thermostat"
	" with gamma=%f\n",count_rescale,FricCoef);

	printf("Now start running for %d timesteps with fixed step size"
	" dt=%f for thermostat equilibration\n",Langevin_equi_steps, dt_final);

	for(i=0;i<Langevin_equi_steps;i++)
	{
		count_steps=count_steps+1;
		Velocity_Verlet_Langevin( dbg_print );
		t = t + dt;
		if(count_steps%sample_cycle==0)	
		{		
			fprintf(output_energy_equilibration, "%f ",t);
          		printf("time=%f \n",t);
			force_compute=0;energy_compute=1;
			energy_momentum_check( 0, t );
			force_compute=1;energy_compute=0;
		}
	}
	
    /*****Stage 3: the annealing loop: swtiching between T and T_annealing*****/
    printf("To prevent trapped at a metastable state, we switch to anealling"
" between temperature %f and %f\n",temperature,Temprature_annealing);

    printf("Now start running for %d cycles with high temperature %d steps "
"followed with low temp %d steps for each cycle\n",\
N_annealing_cycle,N_hightemp_annealing, N_lowtemp_annealing);

    count_steps=0;

    for(i=0;i<N_annealing_cycle;i++)
    {
        /*update noise vector for high temperature*/
        
        for (k=0;k<N;k++)
        {
            noise[k]=sqrt(2.0*FricCoef*mass[k]*Temprature_annealing/dt_final);
        }
        
        for(j=0;j<N_hightemp_annealing;j++)
        {
            count_steps=count_steps+1;
            Velocity_Verlet_Langevin( dbg_print );
            t = t + dt;

        if(count_steps%sample_cycle==0)
        {
            fprintf(output_energy_equilibration, "%f ",t);
            printf("time=%f \n",t);
            force_compute=0;energy_compute=1;
            energy_momentum_check( 0, t );
            force_compute=1;energy_compute=0;
            fprintf(trajectory_annealing,\
 "ITEM: TIMESTEP\n%d \nITEM: NUMBER OF ATOMS \n%d \n", count_steps, N);

            fprintf(trajectory_annealing,\
 "ITEM: BOX BOUNDS pp pp pp \n%f %f \n%f %f \n%f %f \n",\
 -Rshell, Rshell, -Rshell, Rshell, -Rshell, Rshell);

            fprintf(trajectory_annealing,\
 "ITEM: ATOMS id type x y z \n");
            record_trajectories(traj_ind);
        }
        }
        
        /*swtich noise back to low temperature*/
        for (k=0;k<N;k++)
        {
            noise[k]=sqrt(2.0*FricCoef*mass[k]*KbT/dt_final);
        }
        
        for(j=0;j<N_lowtemp_annealing;j++)
        {
            count_steps=count_steps+1;
            Velocity_Verlet_Langevin( dbg_print );
            t = t + dt;

            if(count_steps%sample_cycle==0)
            {
                fprintf(output_energy_equilibration, "%f ",t);
                printf("time=%f \n",t);
                force_compute=0;energy_compute=1;
                energy_momentum_check( 0, t );
                force_compute=1;energy_compute=0;
                fprintf(trajectory_annealing,\
 "ITEM: TIMESTEP\n%d \nITEM: NUMBER OF ATOMS \n%d \n", count_steps, N);

                fprintf(trajectory_annealing,\
 "ITEM: BOX BOUNDS pp pp pp \n%f %f \n%f %f \n%f %f \n",\
 -Rshell, Rshell, -Rshell, Rshell, -Rshell, Rshell);

                fprintf(trajectory_annealing, "ITEM: ATOMS id type x y z \n");

                record_trajectories(traj_ind);
            }
        }
    }

    fclose( trajectory_annealing );
    fclose(output_energy_equilibration);

   	traj_ind=1;//store the production trajectory in another dump file
	iter_indicator=0;//for re-allocate the memory

	/*****Stage 4: the production loop*****/
	set_precision(prec_set);
	printf("#steps=%d, now we are starting the production time steps\n",\
	count_rescale+count_steps);

	printf("Reset the prec to the user-set value: %d digits\n",prec_set);

	printf("Running for %d timesteps with fixed step size dt=%f for data"
	" production\n",Production_steps, dt_final);

	for(i=0;i<Production_steps;i++)
	{		
		if(count_sample%1000==0) 
			printf("#production steps=%d, time=%f\n",\
			count_sample,dt*(count_sample));

		Velocity_Verlet_Langevin( dbg_print );

		iter_indicator=1;
		
		t = t + dt;
		
		count_sample=count_sample+1;

		if(count_sample%sample_cycle==0)
		{
				fprintf(output_energy, "%f ",(count_sample)*dt);
             		 	printf("time=%f \n",t);
				force_compute=0;energy_compute=1;
				energy_momentum_check( 1, t );
				force_compute=1;energy_compute=0;

				fprintf(trajectory,\
	 "ITEM: TIMESTEP\n%d \nITEM: NUMBER OF ATOMS \n%d \n", count_sample, N);

				fprintf(trajectory,\
 "ITEM: BOX BOUNDS pp pp pp \n%f %f \n%f %f \n%f %f \n",\
 -Rshell, Rshell, -Rshell, Rshell, -Rshell, Rshell);

				fprintf(trajectory,\
				 "ITEM: ATOMS id type x y z \n");

				record_trajectories(traj_ind);
#ifdef RDF
				if(Production_steps>=sample_cycle)
				record_rdf();
#endif				
		}		
	}
	
	/*****End of the MD 4-stage loops*****/

	
	/*****record the final configuration*/
	for ( i=0 ; i<N ; i++ )
	{
	fprintf(output_finalconfig, "%d %f %f %f\n",i, x[i], y[i], z[i]);
	}
	
	fprintf( stderr, "\nSimulation Finished!\n");
	
#ifdef RDF
	if(Production_steps>=sample_cycle)
	output_rdf();
#endif

	// close trajectories file
	if(Production_steps>=sample_cycle)
	{
	fclose( trajectory );
	fclose(output_energy);
	}
	fclose(output_finalconfig);

}

void set_precision(int prec)
{
	int i;
	i=prec-3;
	p=int(order_p[i]);
	im=int(im_num[i]);
	imm=int(imm_num[i]);
	sourcetol=source_tol[i];
	sphtol=sph_tol[i];
	gmrestol=1.0/pow(10.0,gmres_tol[i]);
	fmmtol=int(fmm_tol[i]);
	
	allocate_dynamic();
}

void allocate_dynamic()
{
    int i,j,k;
    M=N_col;
    nphi=2*p;
    ntheta=2*p;
    imnum_thresh=N_col*N*im;
    
    b=new double [M*(p+1)*(2*p+1)*2];
    Bknm1D=new double [M*(p+1)*(2*p+1)*2];
    
    imx = new double[imnum_thresh];
    imy = new double[imnum_thresh];
    imz = new double[imnum_thresh];
    imq = new double[imnum_thresh];
    imind = new int[imnum_thresh];
    
    srcPosarr=new double[3*N];
    srcDenarr=new double[(p+1)*(p+1)*N];
    pot_m=new double[3*N];
    
    sqrtk=new double[p+1];
    
    ionind= new int*[N];
    ionind[0] = new int[N*N];
    for (i=1; i<N; i++)
        ionind[i] = ionind[0] + i*N;
    
    sph_imind= new int*[M];
    sph_imind[0] = new int[imnum_thresh];
    for (i=1; i<M; i++)
        sph_imind[i] = sph_imind[0] + i*imnum_thresh/M;
    
    weights= new double*[2*p];
    weights[0] = new double[4*p*p];
    for (i=1; i<2*p; i++)
        weights[i] = weights[0] + i*2*p;
    
    rnodes=new double**[2*p];
    rnodes[0]=new double*[2*p*2*p];
    for (i=1; i<2*p; i++)
        rnodes[i] = rnodes[0] + i*2*p;
    rnodes[0][0]=new double[4*p*p*3];
    for (i=0; i<2*p; i++)
        for(j=0;j<2*p;j++)
            rnodes[i][j] = rnodes[0][0] + i*2*p*3+j*3;
    
    ynm=new complex**[2*p];
    ynm[0]=new complex*[2*p*(2*p+1)];
    for (i=1; i<2*p; i++)
        ynm[i] = ynm[0] + i*(2*p+1);
    ynm[0][0]=new complex[2*p*(2*p+1)*(p+1)];
    for (i=0; i<2*p; i++)
        for(j=0;j<2*p+1;j++)
            ynm[i][j] = ynm[0][0] + i*(2*p+1)*(p+1)+j*(p+1);
    
    Bknm=new complex**[N_col];
    Bknm[0]=new complex*[N_col*(2*p+1)];
    for (i=1; i<N_col; i++)
        Bknm[i] = Bknm[0] + i*(2*p+1);
    Bknm[0][0]=new complex[N_col*(2*p+1)*(p+1)];
    for (i=0; i<N_col; i++)
        for(j=0;j<2*p+1;j++)
            Bknm[i][j] = Bknm[0][0] + i*(2*p+1)*(p+1)+j*(p+1);
    
    BknmCopy=new complex**[N_col];
    BknmCopy[0]=new complex*[N_col*(2*p+1)];
    for (i=1; i<N_col; i++)
        BknmCopy[i] = BknmCopy[0] + i*(2*p+1);
    BknmCopy[0][0]=new complex[N_col*(2*p+1)*(p+1)];
    for (i=0; i<N_col; i++)
        for(j=0;j<2*p+1;j++)
            BknmCopy[i][j] = BknmCopy[0][0] + i*(2*p+1)*(p+1)+j*(p+1);
    
    Bcs=new double**[N_col];
    Bcs[0]=new double*[N_col*(2*p+1)];
    for (i=1; i<N_col; i++)
        Bcs[i] = Bcs[0] + i*(2*p+1);
    Bcs[0][0]=new double[N_col*(2*p+1)*(p+1)];
    for (i=0; i<N_col; i++)
        for(j=0;j<2*p+1;j++)
            Bcs[i][j] = Bcs[0][0] + i*(2*p+1)*(p+1)+j*(p+1);
    
    target= new double*[M*4*p*p];
    target[0] = new double[M*12*p*p];
    for (i=1; i<M*4*p*p; i++)
        target[i] = target[0] + i*3;
    
    fldtarg= new complex*[M*4*p*p];
    fldtarg[0] = new complex[M*12*p*p];
    for (i=1; i<M*4*p*p; i++)
        fldtarg[i] = fldtarg[0] + i*3;
    
    pottarg=new complex[M*4*p*p];
    
    fgrid=new complex**[N_col];
    fgrid[0]=new complex*[N_col*2*p];
    for (i=1; i<N_col; i++)
        fgrid[i] = fgrid[0] + i*2*p;
    fgrid[0][0]=new complex[N_col*4*p*p];
    for (i=0; i<N_col; i++)
        for(j=0;j<2*p;j++)
            fgrid[i][j] = fgrid[0][0] + i*(2*p)*(2*p)+j*(p*2);
    
    ycoef=new complex**[N_col];
    ycoef[0]=new complex*[N_col*(2*p+1)];
    for (i=1; i<N_col; i++)
        ycoef[i] = ycoef[0] + i*(2*p+1);
    ycoef[0][0]=new complex[N_col*(2*p+1)*(p+1)];
    for (i=0; i<N_col; i++)
        for(j=0;j<2*p+1;j++)
            ycoef[i][j] = ycoef[0][0] + i*(2*p+1)*(p+1)+j*(p+1);
    
    Mnodes=new double***[N_col];
    Mnodes[0]=new double**[N_col*2*p];
    for (i=1; i<N_col; i++)
        Mnodes[i] = Mnodes[0] + i*2*p;
    Mnodes[0][0]=new double*[N_col*4*p*p];
    for (i=0; i<N_col; i++)
        for(j=0;j<2*p;j++)
            Mnodes[i][j] = Mnodes[0][0] + i*(2*p)*(2*p)+j*(p*2);
    Mnodes[0][0][0]=new double[N_col*4*p*p*3];
    for (i=0; i<N_col; i++)
        for(j=0;j<2*p;j++)
            for(k=0;k<2*p;k++)
                Mnodes[i][j][k] = Mnodes[0][0][0]+i*(2*p)*(2*p)*3+j*(p*2)*3+k*3;
    
    mpole=new complex***[N_col];
    mpole[0]=new complex**[N_col*N_col];
    for (i=1; i<N_col; i++)
        mpole[i] = mpole[0] + i*N_col;
    mpole[0][0]=new complex*[N_col*N_col*(2*p+1)];
    for (i=0; i<N_col; i++)
        for(j=0;j<N_col;j++)
            mpole[i][j] = mpole[0][0] + i*N_col*(2*p+1)+j*(p*2+1);
    mpole[0][0][0]=new complex[N_col*N_col*(2*p+1)*(p+1)];
    for (i=0; i<N_col; i++)
        for(j=0;j<N_col;j++)
            for(k=0;k<2*p+1;k++)
                mpole[i][j][k]=\
		mpole[0][0][0]+i*M*(2*p+1)*(p+1)+j*(p*2+1)*(p+1)+k*(p+1);
    
    local=new complex***[N_col];
    local[0]=new complex**[N_col*N_col];
    for (i=1; i<N_col; i++)
        local[i] = local[0] + i*N_col;
    local[0][0]=new complex*[N_col*N_col*(2*p+1)];
    for (i=0; i<N_col; i++)
        for(j=0;j<N_col;j++)
            local[i][j] = local[0][0] + i*N_col*(2*p+1)+j*(p*2+1);
    local[0][0][0]=new complex[N_col*N_col*(2*p+1)*(p+1)];
    for (i=0; i<N_col; i++)
        for(j=0;j<N_col;j++)
            for(k=0;k<2*p+1;k++)
                local[i][j][k] =\
		local[0][0][0] + i*M*(2*p+1)*(p+1)+j*(p*2+1)*(p+1)+k*(p+1);
    
    source= new double*[imnum_thresh+N];
    source[0] = new double[(imnum_thresh+N)*3];
    for (i=1; i<imnum_thresh+N; i++)
        source[i] = source[0] + i*3;
    
    mpsource= new double*[N-N_ion];
    mpsource[0] = new double[(N-N_ion)*3];
    for (i=1; i<N-N_ion; i++)
        mpsource[i] = mpsource[0] + i*3;
    
    dipvec= new double*[imnum_thresh+N];
    dipvec[0] = new double[(imnum_thresh+N)*3];
    for (i=1; i<imnum_thresh+N; i++)
        dipvec[i] = dipvec[0] + i*3;
    
    mpdipvec= new double*[N_col];
    mpdipvec[0] = new double[(N_col)*3];
    for (i=1; i<N_col; i++)
        mpdipvec[i] = mpdipvec[0] + i*3;
    
    
    fld= new complex*[imnum_thresh+N];
    fld[0] = new complex[(imnum_thresh+N)*3];
    for (i=1; i<imnum_thresh+N; i++)
        fld[i] = fld[0] + i*3;
    
    charge=new complex[imnum_thresh+N];
    dipstr=new complex[imnum_thresh+N];
    mpdipstr=new complex[N_col];
    mptarget=new double[3];
    pot=new complex[imnum_thresh+N];
    scarray=new double[10*(p+2)*(p+2)];
    wlege=new double[10*(p+2)*(p+2)];
}


