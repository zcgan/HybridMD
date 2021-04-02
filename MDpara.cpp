#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int N;    //Total number of particles.

int Ntype; 
//default to be 2, ions and colloids, they are treated differently because colloids can have dielectric.

int N_ion,N_ion1,N_ion2; //number of ions

int N_col,M,N_col1,N_col2; //number of colloids

// locations and charges of the sources
double *x;
double *y;
double *z;
double *q;

//charge and radius of the central spherical interface
double QM;
double RM;

//radius of the sources
double *r;

/*now we first assume ion and colloid have their own size. 
However, we can further add that ions or colloids have different sizes.*/
/*Gan note 10/24/2019: already added two type of colloids.*/
double r_ion;
double r_col,r_col1,r_col2;
double q_ion,q_ion1,q_ion2;
double q_col,q_col1,q_col2;

// mass of ions and colloids
double m_col1,m_col2,m_ion;

// velocities, of the sources
double *vx ;
double *vy ;
double *vz ;

// accelerations of the sources
double *acc_x;
double *acc_y;
double *acc_z;

// Forces of the sources
double *Fx;
double *Fy;
double *Fz;

//KbT=temperature
double KbT;

//noise for random force
double *noise;

//the friction coeff for langevin thermostat.
double FricCoef; 

//random force for langevin thermostat.
double *Fran_x;
double *Fran_y;
double *Fran_z;

//mass of the sources
double *mass;

//the total mass
double tot_mass; 

// System momentum
double moment_x, moment_y, moment_z;

//file storing the trajectory
FILE *trajectory;
FILE *trajectory_annealing;

//file for output the simulation process
FILE *output_equilibration;
FILE *output_production;
FILE *output_energy;
FILE *output_energy_equilibration;
FILE *output_finalconfig;

//file storing the rdf
FILE *rdf_file;

//arrays for rdf recording statistics 
double numsample_rdf=0.0;
double *rdf;
double bin_size;
int bin_num;
double *bin_ii;
double *bin_ic;
double *bin_cc;
double *bin_ic1;
double *bin_ic2;
double *bin_c1c1;
double *bin_c1c2;
double *bin_c2c2;

int sample_cycle;

//the boundary shell radius
double Rshell;

// time grid
double dt, tmin, tmax;
int step_max, step_rescale, step_burnin; 
int time_stepratio;

double dt_initial,dt_final;

int dt_inc_interval,Langevin_equi_steps,Production_steps,Sample_interval;

int N_hightemp_annealing, N_lowtemp_annealing, N_annealing_cycle;

double Temprature_annealing;

double t_rescale, t_burnin;

double size_scaling=1.0;

//temperature unit(infact energy unit kbt)
double temperature;

//in the 1-st stage, every rescale_cycle time step, we do a velocity rescaling.
//in order to "calm down" the system
int rescale_cycle;

/*variables for the cell-list algorithm*/
int Cell_num_1D; //The number of cells in each dimension.
int Cell_num_3D; //The number of cells in 3D.
double celli;
int *head; //cell-list head array
int *list; //cell-list particle list array.

/*variables for the shift-truncated Lennard-Jones potential & force*/

double c_lj; 
/*distance unit constant. 
By default, choose to be the diameter of the smallest particle.*/


double Rc_lj, Rc_lj1, Rc_lj2; 
/*the truncated distance for deciding the cell-list cell size. 
By default, choose to be the Rc for the largest particle.*/

double factor_lj=1.122462048309373; //2^(1/6)

int iter_indicator=0;
int equi_indicator=0;

char outputname[100];
char paraname[100];
char config_name[100];

double ei, ei1,ei2;
double epsi_ion;
double epsi_M;

/*indicator for energy or force computing*/
int energy_compute=0, force_compute=1;
double printed_ele_energy, printed_lj_energy;


/*random seed*/
Ran myran(1);

/*precision setting parameters*/
double precision[4],order_p[4],im_num[4],imm_num[4],\
source_tol[4], sph_tol[4],gmres_tol[4],fmm_tol[4];

int prec_set, read_config;


/*----------------hybrid method parameters and arrays------------------------*/
double *srcPosarr;//[3*Nmax];
double *srcDenarr;//[4*Nmax];
double *pot_m;//[3*Nmax];


/*global parameters for the system*/

int p; //the truncated order of the infinite sum moment series.
int im; //the number of image charges we use for the near interface charges.
int imm; //# of discrete multipoles we use for the image multipole integral.
//int M=200; //the number of dielectric spheres in the environment.
//double a=1.0; //the sphere radius
//double tau=0.2;

//the diectric coefficients for the environment and for each dielectric sphere.
double epsi_s;
double *epsi_i;
double *lamda_i;
double *beta_i;
double *gamma_i;


//allocate memory for the locations and radius for the dielectric spheres.
double *ox;
double *oy;
double *oz;
double *orad;
double *osigma;

//paramters for the hybrid method criterion.

// allocate memory for the locations and charges of the image charges.
double *imx; //= new double[Nimax*Ncmax*IMmax];//N_ion*im*M
double *imy; //= new double[Nimax*Ncmax*IMmax];
double *imz; //= new double[Nimax*Ncmax*IMmax];
double *imq; //= new double[Nimax*Ncmax*IMmax];
int *imind;//=new int[Nimax*Ncmax*IMmax];
int **ionind;//[Nmax][Nmax];
int **sph_imind;//[Nc][20*im]

double **rotangle;

double sourcetol;
// if for a source charge, dr/rad<sourcetol, then generate its image charges.

//double corrpara=1.4; 
/*if the ion distance with A and B <corrpara simutaniously,
add the force corrections.
Gan note 10/24/2019: no longer used. Current force calculation is exact.
NO CORRECTION IS NEEDED*/

double sphtol; 
//if the distance between 2 sphs d<sphtol, generate the image multipole coef

double m2ltr;//64.0 is a good choice in many case;

double gmrestol;
int fmmtol;

//precomputed factorials for the coefficient recurrence relations.
double **ncoef1;
double ***ncoef2;
double ****ncoef3;

double **ocoef1;
double ***ocoef2;
double ****ocoef3;

//the array storing all the temporary coefficients
double ******tempcoef;
double ******tempcoefimag;

//parameters for the Gauss-type quadrature
double *wquad,*xquad;//storing the quadrature weights and abscissas.

//precompute the powers, it is not necessary to do that in each time-step
double *powerd;// = new double[2*p+2];
double **powerxquad;
double ***newpowerd;//[Ncmax][Ncmax][100];


//arrays in the matvec function
complex ****mpole;//[Ncmax][Ncmax][2*Pmax+1][Pmax+1];
//complex mpolerot[Ncmax][Ncmax][2*Pmax+1][Pmax+1];
complex ****local;//[Ncmax][Ncmax][2*Pmax+1][Pmax+1];

/*
complex immpoletemp[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletemp_i[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletemprot[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletemprot_i[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletrans[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletrans_i[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletransmpmp[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
double oim[Ncmax][Ncmax][IMMmax][3]; 
//oim is the image multipole center inside each sphere due to all other spheres. 
thus oim[M][M][imm][3]
double oim_i[Ncmax][Ncmax][IMMmax][3]; 
//oim_i is the image multipole center outside each sphere due to all other
 spheres. thus oim[M][M][imm][3]
*/

double ***rnodes;
double **weights;
complex ***ynm;

//double rnodes[2*Pmax][2*Pmax][3]; 
//the fortran code uses static arrays, 
//so we define here: rnodes[ntheta][nphi][3], weights[ntheta][nphi].
//double weights[2*Pmax][2*Pmax];
//complex ynm[2*Pmax][2*Pmax+1][Pmax+1]; //ynm[ntheta][2*p+1][p+1]

double ****Mnodes;//static array Mnodes[M][ntheta][nphi][3].

//evaluate the potentials on each nodes using FMM methods.
complex ***fgrid; //fgrid[M][ntheta][nphi]
complex ***ycoef;//ycoef[M][2p+1][p+1]

/*arrays for FMM3dlib calls*/
double **source, **dipvec, **mpsource,\
**mpdipvec,*mptarget,*scarray,*wlege;//[nim+N][3],dipvec[nim+N][3];

complex *charge,*dipstr,*pot,**fld,*mpdipstr;
//[nim+N],dipstr[nim+N],pot[nim+N],fld[nim+N][3];

double **target; //Nc*2p*2p,3
complex *pottarg;//Nc*2p*2p
complex **fldtarg;//Nc*2p*2p,3
complex ***Bknm;//[Ncmax][2*Pmax+1][Pmax+1]; 
/*the multipole expansion coefficients Bknm[M][2p+1][p+1] for all the spheres, 
which is infact the unknown vector x.*/
complex ***BknmCopy;
double ***Bcs;

double *wwquad,*xxquad;
/*storing the quadrature weights and abscissas. 
since the image line charge integral is singular, 
we use Gauss-Jacobi quadrature to increase the accuracy.*/

double *sqrtk;//=new double[p+1];
double *b;
double *Bknm1D;

/*calling the sshgrid function to obtain the x,y,z coordinates
 and their weights on a unit sphere.*/

int nphi;//=2*p;
int ntheta;//=2*p;
int nnodes;

int imnum_thresh;
int FMM_thresh;
int *count_nim;
/*------------------------------------------------------------------------*/


/*protocals*/
void allocate_arrays();
void read_para( );
void energy_momentum_check( int info, double time );
void GAUSS( double *gauss1, double *gauss2 );
void initial_conditions( );
void LJ_accelerations_Directsum( int iprint );
void Cell_list( int iprint );
void LJForce(int indx, int site);
void Coulomb_accelerations_Directsum( int iprint );
void Coulomb_accelerations_FMM( int iprint );
void Coulomb_accelerations_Hybrid( int iprint );
void LJ_Boundaryforce( int iprint );
void Central_sph_force( int iprint );
void mass_initialize( int iprint);
void output_force();
void output_x_v();
void record_trajectories(int traj_ind);
void initialize_rdf();
void record_rdf();
void output_rdf();
void Velocity_Rescale( );
void Velocity_Verlet( int iprint );
void Velocity_Verlet_Langevin( int iprint );
void Velocity_Rescale_ion( );
void Velocity_Verlet_ion( int iprint );
void Velocity_Verlet_Langevin_ion( int iprint );
