#ifndef _MDPARA_H_
#define _MDPARA_H_

#include"ran.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*Gan note 10/24/2019: for more detailed comments, look at MDpara.cpp instead*/

extern int N;    //number of particles.

extern int Ntype; 
//default to be 2, ions and colloids, 
//they are treated differently because colloids can have dielectric.

extern int N_ion; //number of ions

extern int N_col,M, N_col1, N_col2, N_ion1, N_ion2; //number of colloids

// locations and charges of the sources
extern double *x;
extern double *y;
extern double *z;
extern double *q;

//charge and radius of the central spherical interface
extern double QM;
extern double RM;


//radius of the sources
extern double *r;

extern double r_ion; 
extern double r_col,r_col1,r_col2;
extern double q_ion;
extern double q_col, q_col1, q_col2, q_ion1, q_ion2;
extern double m_col1, m_col2, m_ion;

// velocities, of the sources
extern double *vx ;
extern double *vy ;
extern double *vz ;

// accelerations of the sources
extern double *acc_x;
extern double *acc_y;
extern double *acc_z;

// Forces of the sources
extern double *Fx;
extern double *Fy;
extern double *Fz;

//KbT
extern double KbT;

//noise for random force
extern double *noise;

//the friction coeff for langevin thermostat.
extern double FricCoef; 

//random force for langevin thermostat.
extern double *Fran_x;
extern double *Fran_y;
extern double *Fran_z;

//mass of the sources
extern double *mass;

//the total mass
extern double tot_mass; 

// System momentum
extern double moment_x, moment_y, moment_z;

//file storing the trajectory
extern FILE *trajectory;
extern FILE *trajectory_annealing;
//file for output the simulation process
extern FILE *output_equilibration;
extern FILE *output_production;
extern FILE *output_energy;
extern FILE *output_energy_equilibration;
extern FILE *output_finalconfig;

//file storing the rdf
extern FILE *rdf_file;

//arrays for rdf recording statistics 
extern double numsample_rdf;
extern double *rdf;
extern double bin_size;
extern int bin_num;
extern double *bin_ii;
extern double *bin_ic;
extern double *bin_cc;
extern double *bin_ic1;
extern double *bin_ic2;
extern double *bin_c1c1;
extern double *bin_c1c2;
extern double *bin_c2c2;

extern int sample_cycle;

//the boundary shell radius
extern double Rshell;

// time grid
extern double dt, tmin, tmax; 

extern int step_max, step_rescale, step_burnin; 

extern int time_stepratio;

extern double dt_initial,dt_final;

extern int dt_inc_interval,Langevin_equi_steps,Production_steps,Sample_interval;

extern int N_hightemp_annealing, N_lowtemp_annealing, N_annealing_cycle;

extern double Temprature_annealing;

extern double t_rescale, t_burnin;

extern double size_scaling;

extern double temperature;

//in the beginning of a MD run, every rescale_cycle time step, 
//we do a velocity rescaling.
extern int rescale_cycle;

/*variables for the cell-list algorithm*/
extern int Cell_num_1D; //The number of cells in each dimension.
extern int Cell_num_3D; //The number of cells in 3D.
extern double celli;
extern int *head; //cell-list head array
extern int *list; //cell-list particle list array.

/*variables for the shift-truncated Lennard-Jones potential & force*/
extern double c_lj; 
extern double Rc_lj,Rc_lj1,Rc_lj2; 
extern double factor_lj; //2^(1/6)

/*complex number structure*/
struct complex {double real,imag;}; 

extern int iter_indicator;
extern int equi_indicator;

extern char outputname[100];
extern char paraname[100];
extern char config_name[100];

extern double ei,ei1,ei2;
extern double epsi_ion;
extern double epsi_M;

/*indicator for energy or force computing*/
extern int energy_compute, force_compute;
extern double printed_ele_energy, printed_lj_energy;

/*random seed*/
extern Ran myran;

/*precision setting parameters*/
extern double precision[4],order_p[4],im_num[4],\
imm_num[4],source_tol[4], sph_tol[4],gmres_tol[4],fmm_tol[4];
extern int prec_set, read_config;


/*----------------hybrid method parameters and arrays------------------------*/
extern double *srcPosarr;//[3*Nmax];
extern double *srcDenarr;//[4*Nmax];
extern double *pot_m;//[3*Nmax];


/*global parameters for the system*/

//int M=200; //the number of dielectric spheres in the environment.
extern int p; //the truncated order of the infinite sum moment series.
extern int im; 
//the number of image charges we use for the near interface charges.
extern int imm; 
//the number of discrete multipoles we use for the image multipole integral.
//double a=1.0; //the sphere radius
//double tau=0.2;

//the diectric coefficients for the environment and for each dielectric sphere.
extern double epsi_s;
extern double *epsi_i;
extern double *lamda_i;
extern double *beta_i;
extern double *gamma_i;

//allocate memory for the locations and radius for the dielectric spheres.
extern double *ox;
extern double *oy;
extern double *oz;
extern double *orad;
extern double *osigma;

//paramters for the hybrid method criterion.

// allocate memory for the locations and charges of the image charges.
extern double *imx; //= new double[Nimax*Ncmax*IMmax];//N_ion*im*M
extern double *imy; //= new double[Nimax*Ncmax*IMmax];
extern double *imz; //= new double[Nimax*Ncmax*IMmax];
extern double *imq; //= new double[Nimax*Ncmax*IMmax];
extern int *imind;//=new int[Nimax*Ncmax*IMmax];
extern int **ionind;//[Nmax][Nmax];
extern int **sph_imind;//[Nc][20*im]

extern double **rotangle;

extern double sourcetol;
// if for a source charge, dr/rad<sourcetol, then generate its image charges.

extern double sphtol; 
extern double m2ltr;//64.0;

extern double gmrestol;
extern int fmmtol;

//precomputed factorials for the coefficient recurrence relations.
extern double **ncoef1;
extern double ***ncoef2;
extern double ****ncoef3;

extern double **ocoef1;
extern double ***ocoef2;
extern double ****ocoef3;

//the array storing all the temporary coefficients
extern double ******tempcoef;
extern double ******tempcoefimag;

//parameters for the Gauss-type quadrature
extern double *wquad,*xquad;//storing the quadrature weights and abscissas.

//precompute the powers
extern double *powerd;// = new double[2*p+2];
extern double **powerxquad;
extern double ***newpowerd;//[Ncmax][Ncmax][100];


//arrays in the matvec function
extern complex ****mpole;//[Ncmax][Ncmax][2*Pmax+1][Pmax+1];
//complex mpolerot[Ncmax][Ncmax][2*Pmax+1][Pmax+1];
extern complex ****local;//[Ncmax][Ncmax][2*Pmax+1][Pmax+1];

/*complex immpoletemp[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletemp_i[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletemprot[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletemprot_i[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletrans[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletrans_i[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletransmpmp[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
*/

extern double ***rnodes;
extern double **weights;
extern complex ***ynm;
extern double ****Mnodes;//static array Mnodes[M][ntheta][nphi][3].
//evaluate the potentials on each nodes using FMM methods.
extern complex ***fgrid; //fgrid[M][ntheta][nphi]
extern complex ***ycoef;//ycoef[M][2p+1][p+1]

/*arrays for FMM3dlib calls*/
extern double **source, **dipvec, **mpsource,\
 **mpdipvec,*mptarget,*scarray,*wlege;//[nim+N][3],dipvec[nim+N][3];
extern complex *charge,*dipstr,*pot,**fld, *mpdipstr;
extern double **target; //Nc*2p*2p,3
extern complex *pottarg;//Nc*2p*2p
extern complex **fldtarg;//Nc*2p*2p,3
extern complex ***Bknm;//[Ncmax][2*Pmax+1][Pmax+1]; /
extern complex ***BknmCopy;//[Ncmax][2*Pmax+1][Pmax+1];
extern double ***Bcs;

extern double *wwquad,*xxquad; 
extern double *sqrtk;//=new double[p+1];
extern double *b;
extern double *Bknm1D;

//calling the sshgrid function to obtain 
//the x,y,z coordinates and their weights on a unit sphere.
extern int nphi;//=2*p;
extern int ntheta;//=2*p;
extern int nnodes;

extern int imnum_thresh;
extern int FMM_thresh;
extern int *count_nim;
/*-----------------------------------------------------------------------*/

/*protocals*/
extern void allocate_arrays();
extern void read_para( );
extern void energy_momentum_check( int info, double time );
extern void GAUSS( double *gauss1, double *gauss2 );
extern void initial_conditions( );
extern void LJ_accelerations_Directsum( int iprint );
extern void Cell_list( int iprint );
extern void LJForce(int indx, int site);
extern void Coulomb_accelerations_Directsum( int iprint );
extern void Coulomb_accelerations_FMM( int iprint );
extern void Coulomb_accelerations_Hybrid( int iprint );
extern void LJ_Boundaryforce( int iprint );
extern void Central_sph_force( int iprint );
extern void Velocity_Rescale( );
extern void mass_initialize( int iprint);
extern void output_force();
extern void output_x_v();
extern void record_trajectories(int traj_ind);
extern void initialize_rdf();
extern void record_rdf();
extern void output_rdf();
extern void Velocity_Verlet( int iprint );
extern void Velocity_Verlet_Langevin( int iprint );
extern void Velocity_Rescale_ion( );
extern void Velocity_Verlet_ion( int iprint );
extern void Velocity_Verlet_Langevin_ion( int iprint );
#endif
