//Hybrid method for electrostatics force field with spherical dielectric objects
#include"Numupbound.h"
#include"MDpara.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "itlin.h"
//#include<omp.h>
/*parameters and functions for the iterative GMRES method.*/
void matvec(int n, double *vecx, double *b); //the matrix vector product function.
void preconr(int n, double *x, double *b); // the preconditioner function.
void initialization();
//double gmres_norm2(int n, double *v); // the evaluation of the gmres residue norm.

extern "C" void gmres(int n, double *y, MATVEC *matvec,PRECON *preconr, PRECON *preconl, double *b, struct ITLIN_OPT *opt, struct ITLIN_INFO *info);

//the Gauss-jacobi quadrature function
extern void cgqf ( int nt, int kind, double alpha, double beta, double a, double b,double t[], double wts[]);

//factorial function
double factorial(int n);

//generate the constant coefficients for the reccurence relation of the multipole expansion coefficients
void generate_coeff();

//generate the image multipoles for a given set of multipole coefficients
void generate_imagecoef(int n, double *vecx);
void generate_power();

/*functions and structure involved with the spherical harmonics (in fortran) */

//struct complex {double real,imag;};  //define a complex structure to make the  C code consistent with fortran codes.

extern "C" void sshgrid_(int *,int *,double [][2*Pmax][3],double [][2*Pmax],int *);
extern "C" void sshgfun_(int *,int *,complex [][2*Pmax+1][Pmax+1]);
extern "C" void sshcoef_(complex[][2*Pmax],int *,int *,int *,complex [][2*Pmax+1][Pmax+1],complex [][Pmax+1]);
extern "C" void sshevalg_(complex [][Pmax+1],int *,int *,int*,complex[][2*Pmax+1][Pmax+1],complex[][2*Pmax]);
extern "C" void ssheval_(complex [][Pmax+1],int *,double [3], complex *);
extern "C" void l3dmplocquadu_(double *, double [],complex[][Pmax+1],int *,double *, double [], complex[][Pmax+1], int *, int *);
extern "C" void l3dmpmpquadu_(double *, double [],complex[][Pmax+1],int *,double *, double [], complex[][Pmax+1], int *, int *);
extern "C" void lfmm3dparttarg_(int *,int *,int *,double [][3], int *,complex [],int *,complex[],double [][3],int *,complex [],int *,complex[][3],int *,double [][3],int *,complex[],int *,complex[][3]);
extern "C" void lfmm3dpartself_(int *,int *,int *,double [][3], int *,complex [],int *,complex[],double [][3],int *,complex [],int *,complex[][3]);
extern "C" void rotviarecur3f90_(double *,int *,int *,int *,complex[][Pmax+1],int *,complex[][Pmax+1],int *);


double diag = 1.0; //the right-handside preconditioner factor.
double pi=4.0*atan(1.0),npi=-1.0*pi;

double srcPosarr[3*Nmax];
double srcDenarr[4*Nmax];
double pot_m[3*Nmax];


/*global parameters for the system*/

//int M=200; //the number of dielectric spheres in the environment.
int p=1; //the truncated order of the infinite sum moment series.
int im=2; //the number of image charges we use for the near interface charges.
int imm=5; //the number of discrete multipoles we use for the image multipole integral.
//double a=1.0; //the sphere radius
//double tau=0.2;

//the diectric coefficients for the environment and for each dielectric sphere.
double epsi_s;
double *epsi_i;

//allocate memory for the locations and radius for the dielectric spheres.
double *ox;
double *oy;
double *oz;
double *orad;
double *osigma;

//paramters for the hybrid method criterion.

// allocate memory for the locations and charges of the image charges.
double *imx = new double[Nimax*Ncmax*IMmax];//N_ion*im*M
double *imy = new double[Nimax*Ncmax*IMmax];
double *imz = new double[Nimax*Ncmax*IMmax];
double *imq = new double[Nimax*Ncmax*IMmax];
int *imind=new int[Nimax*Ncmax*IMmax];
int ionind[Nmax][Nmax];

double **rotangle;

double sourcetol=2.0;// if for a source charge, dr/rad<sourcetol, then generate its image charges.
//double corrpara=1.4; //if the ion distance with A and B <corrpara simutaniously, then add the force corrections.
double sphtol=100.0; //if the distance between two spheres d<sphtol, then generate the image multipole coefficients.
double m2ltr=16;//64.0;
//double D=10.0; //the contact distance for the spheres.

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

//precompute the powers, so that it is not necessary to do that again and again in GMRES iterations.
double *powerd = new double[2*p+2];
double **powerxquad;
double newpowerd[Ncmax][Ncmax][100];


//arrays in the matvec function
complex mpole[Ncmax][Ncmax][2*Pmax+1][Pmax+1];
complex mpolerot[Ncmax][Ncmax][2*Pmax+1][Pmax+1];
complex local[Ncmax][Ncmax][2*Pmax+1][Pmax+1];

complex immpoletemp[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletemp_i[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletemprot[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletemprot_i[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletrans[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletrans_i[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
complex immpoletransmpmp[Ncmax][Ncmax][IMMmax][2*Pmax+1][Pmax+1];
double oim[Ncmax][Ncmax][IMMmax][3]; //oim is the image multipole center inside each sphere due to all other spheres. thus oim[M][M][imm][3]
double oim_i[Ncmax][Ncmax][IMMmax][3]; //oim_i is the image multipole center outside each sphere due to all other spheres. thus oim[M][M][imm][3]

double rnodes[2*Pmax][2*Pmax][3]; //the fortran code uses static arrays, so we define here: rnodes[ntheta][nphi][3], weights[ntheta][nphi].
//double **weights;
double weights[2*Pmax][2*Pmax];
complex ynm[2*Pmax][2*Pmax+1][Pmax+1]; //ynm[ntheta][2*p+1][p+1]
double Mnodes[Ncmax][2*Pmax][2*Pmax][3]; //static array Mnodes[M][ntheta][nphi][3].
//evaluate the potentials on each nodes using FMM methods.
complex fgrid[Ncmax][2*Pmax][2*Pmax]; //fgrid[M][ntheta][nphi]
complex ycoef[Ncmax][2*Pmax+1][Pmax+1]; //ycoef[M][2p+1][p+1]


double target[Ncmax*Nimax*IMmax][3]; //here (7) is infact imcount=M*N*im.
complex pottarg[Ncmax*Nimax*IMmax],fldtarg[Ncmax*Nimax*IMmax][3];
complex Bknm[Ncmax][2*Pmax+1][Pmax+1]; //the multipole expansion coefficients Bknm[M][2p+1][p+1] for all the spheres, which is infact the unknown vector x.
struct ITLIN_OPT *opt;
struct ITLIN_INFO *info;
int  maxiter = 1000;
int rcode;

double *wwquad,*xxquad;//storing the quadrature weights and abscissas. since the image line charge integral is singular, we use Gauss-Jacobi quadrature to increase the accuracy.
double *sqrtk=new double[p+1];
double *b;
double *Bknm1D;

//calling the sshgrid function to obtain the x,y,z coordinates and their weights on a unit sphere.
int nphi=2*p;
int ntheta=2*p;
int nnodes;

void Coulomb_accelerations_Hybrid( int iprint )
{	

	int i,j,k,l,ii,jj,kk,ind;
	double tempx,tempy,tempz,rji,powr;
	double aa,bb,cc;
	double vec[5];
	M=N_col;
	
	if(iter_indicator==0)
		initialization();

	/*allocate the GMRES solver options*/

	if(iter_indicator==0)
	{
		opt= (struct ITLIN_OPT*)malloc(sizeof(struct ITLIN_OPT));
		info= (struct ITLIN_INFO*)malloc(sizeof(struct ITLIN_INFO));
		TERM_CHECK termcheck = CheckEachIter;
	//	fprintf(stdout,"\n Start of of GMRES\n\n");

		opt->tol = 1.0e-4; //iteration tolerence
		opt->i_max = 10; //the maximum number of iteration for restart.
		opt->termcheck = termcheck;
		opt->maxiter = maxiter;
		opt->errorlevel = Verbose;
		opt->monitorlevel = None;
		opt->datalevel = None;
		opt->errorfile = stdout;
		opt->monitorfile = NULL;
		opt->datafile =NULL;// fopen("test_gmres.data","w");
		//if (!opt->datafile) fprintf(stdout,"\n open of test_gmres.data failed!\n");
		opt->iterfile =NULL;//fopen("test_gmres_iter.data","w");
		opt->resfile  =NULL;// fopen("test_gmres_res.data","w");
		opt->miscfile =NULL;// fopen("test_gmres_misc.data","w");
		
		b=new double [M*(p+1)*(2*p+1)*2];
		Bknm1D=new double [M*(p+1)*(2*p+1)*2];
		
		epsi_i= new double[M];

		//allocate memory for the locations and radius for the dielectric spheres.
		ox = new double[M];
		oy = new double[M];
		oz = new double[M];
		orad = new double[M];
		osigma = new double [M];
	}
	//allocate all the spherical harmonic coefficient arrays

	//spheres location
	for(i=0;i<M;i++)
	{
		ind=i+N_ion;
		ox[i]=x[ind];
		oy[i]=y[ind];
		oz[i]=z[ind];
		orad[i]=1.0;
		osigma[i]=q[ind];
	}

	//set the dielectric constants.
	if(iter_indicator==0)
	{
		epsi_s=epsi_ion;

		for(i=0;i<M;i++)
		{epsi_i[i]=ei;}
	}

	if(iter_indicator==0)
	{
		for(k=0;k<p+1;k++)
		{

			sqrtk[k]=sqrt(2.0*k+1.0);
		} //rescale.
	}


	if(iter_indicator==0)
		sshgrid_(&nphi,&ntheta,rnodes,weights,&nnodes);
	

	//calculate all the ssh functions we will use in the ssh transform.
	if(iter_indicator==0)
		sshgfun_(&p,&ntheta,ynm);


	/*After all the initialization, starts to compute the computation time.*/
	clock_t tstart_total, tfinish_total;
	clock_t tstart_fmm, tfinish_fmm;
	clock_t tstart_sht, tfinish_sht;
	clock_t tstart_gmres, tfinish_gmres;
	clock_t tstart_direct, tfinish_direct;

	double duration;
	//tstart_total = clock();

	//translate all the coordinates onto each spheres.


	for (i=0;i<M;i++)
		for(j=0;j<ntheta;j++)
			for(k=0;k<nphi;k++)
			{
				Mnodes[i][j][k][0]=rnodes[j][k][0]*orad[i]+ox[i];
				Mnodes[i][j][k][1]=rnodes[j][k][1]*orad[i]+oy[i];
				Mnodes[i][j][k][2]=rnodes[j][k][2]*orad[i]+oz[i];
				//printf("rnodes=%f\n",rnodes[i][j][1]);
			}
		//	getchar();




			//First find out whether a source charge approaches the interface, if so, add it's first level images to u1, which we compute in the following.
			//	tstart_total = clock();
			int imcount=0; //count the total number of image charges.
			double dx,dy,dz,dr,rk,rim;
			double start=0;
			double alpha=0.0;
			double end;
			double lamda=epsi_s/(epsi_s+epsi_i[0]),beta=lamda-1.0,gamma=(epsi_i[0]-epsi_s)/(epsi_s+epsi_i[0]); //if we first assume epsi_i are all the same for different dielectrics.
			int kind=4;
			int order=im-1;

			if(iter_indicator==0)
			{
				wwquad = new double[order];
				xxquad = new double[order];
			}


			for (i=0;i<N;i++)
				for(j=0;j<M;j++)
				{

					dx=x[i]-ox[j];
					dy=y[i]-oy[j];
					dz=z[i]-oz[j];
					dr=sqrt(dx*dx+dy*dy+dz*dz);

					if(dr/orad[j]<sourcetol&&dr>orad[j])
					{
						ionind[i][j]=1;//indicating there is image for i inside j
						rk=orad[j]*orad[j]/dr;
						rim=rk/dr;
						imx[imcount]=rim*dx+ox[j];
						imy[imcount]=rim*dy+oy[j];
						imz[imcount]=rim*dz+oz[j];
						imq[imcount]=-gamma*orad[j]*q[i]/dr;  //the Kelvin image.
						imind[imcount]=j+N_ion;
						//the (im-1) line image.
						end=rk;
						cgqf (order, kind, alpha, beta, start, end, xxquad, wwquad);

						for(ii=0;ii<order;ii++)
						{
							ind=imcount+ii+1;
							imx[ind]=ox[j]+xxquad[ii]*dx/dr;
							imy[ind]=oy[j]+xxquad[ii]*dy/dr;
							imz[ind]=oz[j]+xxquad[ii]*dz/dr;
							imq[ind]=wwquad[ii]*gamma*lamda*pow(rk,1.0-lamda)*q[i]/orad[j];
							imind[ind]=j+N_ion;
						}

						imcount=imcount+im;

					}
					else
					{ionind[i][j]=-1;}
				}

				//	tfinish_total = clock();
				/*direct sum test*/



				/*		if(N<12)
				{
				for (i=0;i<M;i++)
				for(j=0;j<ntheta;j++)
				for(k=0;k<nphi;k++)
				{
				dx=x[0]-Mnodes[i][j][k][0];
				dy=y[0]-Mnodes[i][j][k][1];
				dz=z[0]-Mnodes[i][j][k][2];
				dr=sqrt(dx*dx+dy*dy+dz*dz);
				fgrid[i][j][k].real=q[0]/dr/epsi_s;
				fgrid[i][j][k].imag=0.0;
				}
				}*/

				int Ntot=imcount, ier,iprec=3,ifcharge=1,ifdipole=0,ifpot=0,iffld=0,ntarget=M*ntheta*nphi,ifpottarg=1, iffldtarg=0; 
				//here iprec defines the accuracy for FMM. (1, 3digit, 2, 6digit, 3, 9digit.)



				/*	for (i=0;i<imcount;i++)
				{
				source[i][0]=imx[i];
				source[i][1]=imy[i];
				source[i][2]=imz[i];
				charge[i].real=imq[i];
				charge[i].imag=0.0;
				}*/
				//	tstart_fmm= clock();
				for (i=0;i<M;i++)
					for(j=0;j<ntheta;j++)
						for(k=0;k<nphi;k++)
						{
							//pottarg[i*ntheta*nphi+j*nphi+k].real=0.0;
							//pottarg[i*ntheta*nphi+j*nphi+k].imag=0.0;
							fgrid[i][j][k].real=0.0;
							fgrid[i][j][k].imag=0.0;
						}



						//tstart_sht =0;//omp_get_wtime();

						for (i=0;i<M;i++)
							for(j=0;j<ntheta;j++)
								for(k=0;k<nphi;k++)
								{
									//	#pragma omp parallel for reduction(+:pottarg[i*ntheta*nphi+j*nphi+k].real)
									for(ii=0;ii<imcount;ii++)				
									{				
										if(imind[ii]!=i+N_ion)
										{
											dx=imx[ii]-Mnodes[i][j][k][0];
											dy=imy[ii]-Mnodes[i][j][k][1];
											dz=imz[ii]-Mnodes[i][j][k][2];
											dr=sqrt(dx*dx+dy*dy+dz*dz);
											fgrid[i][j][k].real=fgrid[i][j][k].real+imq[ii]/dr;
										}

									}

									for(ii=0;ii<N;ii++)
									{
										if(ionind[ii][i]<0&&ii!=i+N_ion)
										{
											dx=x[ii]-Mnodes[i][j][k][0];
											dy=y[ii]-Mnodes[i][j][k][1];
											dz=z[ii]-Mnodes[i][j][k][2];
											dr=sqrt(dx*dx+dy*dy+dz*dz);
											fgrid[i][j][k].real=fgrid[i][j][k].real+q[ii]/dr;
										}
									}
								}


								//	tfinish_fmm = clock();

								/*for (i=0;i<M;i++)
								for(j=0;j<ntheta;j++)
								for(k=0;k<nphi;k++)
								{
								target[i*ntheta*nphi+j*nphi+k][0]=Mnodes[i][j][k][0];
								target[i*ntheta*nphi+j*nphi+k][1]=Mnodes[i][j][k][1];
								target[i*ntheta*nphi+j*nphi+k][2]=Mnodes[i][j][k][2];
								}

								lfmm3dparttarg_(&ier,&iprec,&Ntot,source,&ifcharge,charge,&ifdipole,dipstr,dipvec,&ifpot,pot,&iffld,fld,&ntarget,target,&ifpottarg,pottarg,&iffldtarg,fldtarg);
								*/

								//	#pragma omp parallel for
								/*	for (i=0;i<M;i++)
								for(j=0;j<ntheta;j++)
								for(k=0;k<nphi;k++)
								{
								fgrid[i][j][k].real=pottarg[i*ntheta*nphi+j*nphi+k].real;
								for(ii=0;ii<N_ion;ii++)
								for(jj=0;jj<im;jj++)
								{
								dx=imx[ii*M*im+i*im+jj]-Mnodes[i][j][k][0];
								dy=imy[ii*M*im+i*im+jj]-Mnodes[i][j][k][1];
								dz=imz[ii*M*im+i*im+jj]-Mnodes[i][j][k][2];
								dr=sqrt(dx*dx+dy*dy+dz*dz);
								fgrid[i][j][k].real=fgrid[i][j][k].real-imq[ii*M*im+i*im+jj]/dr;
								}
								fgrid[i][j][k].real=fgrid[i][j][k].real/epsi_s;
								fgrid[i][j][k].imag=0.0;
								//	printf("wei!%d %d %d %.13f\n",i,j,k,fgrid[i][j][k].real);
								}*/


								//	getchar();

								//calling sshcoef to obtain all the harmonic coefficients around each spheres.
								//	#pragma omp parallel for

								//	tstart_sht = clock();
								for(i=0;i<M;i++)
									sshcoef_(fgrid[i],&p,&nphi,&ntheta,ynm,ycoef[i]);
								//	tfinish_sht = clock();
								//tfinish_sht =0;//omp_get_wtime();

								/*for(i=0;i<M;i++)
								for(j=0;j<2*p+1;j++)
								for(k=0;k<p+1;k++)
								printf("%d %d %d %.27f %.27f\n",i, k, j-p,ycoef[i][j][k].real, ycoef[i][j][k].imag);*/  //check the output of the sshcoef_.


								//using sshevalg to obtain the scalar function fgridappr based on the truncated coefficients and compare with the fgrid.
								//complex fgridappr[2][12][12]; //fgridappr[M][ntheta][nphi]
								//for(i=0;i<M;i++)
								//	sshevalg_(ycoef[i],&p,&nphi,&ntheta,ynm,fgridappr[i]);

								/*	for(i=0;i<M;i++)
								for(j=0;j<2*p;j++)
								for(k=0;k<2*p;k++)
								printf("%d %d %d %f %f %f\n",i, j, k,fgrid[i][j][k].real,fgridappr[i][j][k].real,pottarg[i*ntheta*nphi+j*nphi+k].real/epsi_s); */ //check the difference between fgrid and fgridappr.

								//finally, rearrange the coefficients ycoef, we obtain the right-hand side vector b.
								//tstart_gmres= clock();
								for(i=0;i<M;i++)
									for(j=0;j<2*p+1;j++)
										for(k=0;k<p+1;k++)
										{

											//printf("test real imag %d %d %d %.26f %.26f\n",i, j-p, k,ycoef[i][j][k].real,ycoef[i][j][k].imag);
											if(k==0)
											{
												if(j==p)
													b[i*(2*p+1)*(p+1)+j*(p+1)+k]=0.0; //the scaling factor related to orad need to be checked!
												else
													b[i*(2*p+1)*(p+1)+j*(p+1)+k]=0.0;

												b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]=0.0;
											}
											else
											{
												b[i*(2*p+1)*(p+1)+j*(p+1)+k]=ycoef[i][j][k].real*(epsi_s-epsi_i[i]); //the scaling factor related to orad need to be checked!
												b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]=ycoef[i][j][k].imag*(epsi_s-epsi_i[i]);
											}
										}


										/*we are essentially solving a Ax=b system. till now, the vector b is obtained. In the following, we are doing the matrix vector product, i.e., A*x. Then we can use GMRES to obtain the solution.*/

										//initialize Bknm if it is the first MD iteration.
										if(iter_indicator==0)
										{
											for(i=0;i<M;i++)
												for(j=0;j<2*p+1;j++)
													for(k=0;k<p+1;k++)
													{

														Bknm[i][j][k].real=b[i*(2*p+1)*(p+1)+j*(p+1)+k]; //the initial guess vector is chosen to be the same as vector b. (better guesses?)
														Bknm[i][j][k].imag=b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k];

													}
										}

										//rearrange Bknm.real to obtain the initial vector Bknm1D.

										for(i=0;i<M;i++)
											for(j=0;j<2*p+1;j++)
												for(k=0;k<p+1;k++)
												{
													Bknm1D[i*(2*p+1)*(p+1)+j*(p+1)+k]=Bknm[i][j][k].real;
													Bknm1D[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]=Bknm[i][j][k].imag;
												}
												//perform the GMRES iteration.
												//first generate all the coefficients.
												//generate_power();
												//		generate_coeff();
												//		tstart_gmres = 0;//omp_get_wtime();
												int matsize=2*M*(p+1)*(2*p+1); //the size of the matrix.
												gmres(matsize,Bknm1D,&matvec,NULL,NULL,b,opt,info);
												//			tfinish_gmres = clock();
												//		tfinish_gmres =0;//omp_get_wtime();
												//	generate_imagecoef(matsize,Bknm1D);

												//end of the computation time. We get all the coefficients and then we can compute what ever we want.
												//		tfinish_total = clock();


												/*	for(j=0;j<matsize;j++)
												{
												printf("the solution: x%d=%.7f\n",j,Bknm1D[j]);   //check the solution.
												}*/

												/* close output files and print statistics */
												/*fclose(opt->datafile);
												fclose(opt->iterfile);
												fclose(opt->miscfile);
												fprintf(stdout,"\n Return code:                  %5i\n",info->rcode);
												fprintf(stdout," Iterations:                   %5i\n",info->iter);
												fprintf(stdout," Matvec calls:                 %5i\n",info->nomatvec);
												fprintf(stdout," Right Precon calls:           %5i\n",info->noprecr);
												fprintf(stdout," Left Precon calls:            %5i\n",info->noprecl);
												fprintf(stdout," Estimated residual reduction: %e\n",info->precision);*/


												//Check the accuracy of the obtained harmonic coeeficients Bknm1D by computing the total polarization potential energy.

												//	tstart_direct= clock();

												double energy=0.0;
												complex fval;
												double loc[3];
												double r;
												double aaaaaa=0.0;
												for(i=0;i<M;i++)
													for(j=0;j<2*p+1;j++)
														for(k=0;k<p+1;k++)
														{
															if(abs(j-p)<=k)
															{
																Bknm[i][j][k].real=Bknm1D[i*(2*p+1)*(p+1)+j*(p+1)+k];
																Bknm[i][j][k].imag=Bknm1D[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k];	
																//	printf("%d %d %d %f %f\n",i,k,j-p,Bknm[i][j][k].real,Bknm[i][j][k].imag);	
															}															
														}
														//getchar();
														/*----------------------Compute the force----------------------------*/

														/*testing the force using cartesian multipole expansion*/
														/*First of all, transform Bknm into Cartesian coefficients*/
														if(force_compute==1)
														{

															for(i=0;i<M;i++)
																for(j=0;j<2*p+1;j++)
																	for(k=0;k<p+1;k++)
																	{

																		Bknm[i][j][k].real=Bknm[i][j][k].real*sqrtk[k];
																		Bknm[i][j][k].imag=Bknm[i][j][k].imag*sqrtk[k];
																	} //rescale.

																	double Bcs[Ncmax][2*Pmax+1][Pmax+1];
																	double temper;
																	double sign;
																	for(i=0;i<M;i++)
																		for(j=p;j<2*p+1;j++)
																			for(k=0;k<p+1;k++)
																			{
																				if(abs(j-p)<=k)
																				{
																					if(abs(j-p)==0)
																					{
																						Bcs[i][j][k]=Bknm[i][j][k].real;
																					}
																					else if(abs(j-p)>0)
																					{
																						if(abs(j-p)%2==0)
																							sign=1.41421356;
																						else
																							sign=-1.41421356;

																						temper=(Bknm[i][j][k].real+Bknm[i][2*p-j][k].real)/sign;
																						Bcs[i][j][k]=temper;

																						temper=(-Bknm[i][j][k].imag+Bknm[i][2*p-j][k].imag)/sign;
																						Bcs[i][2*p-j][k]=temper;

																					}

																				}

																			} 




																			//printf("N=%d\n",N);getchar();
																			for(i=0;i<N;i++)
																			{
																				srcPosarr[3*i+0]=x[i];
																				srcPosarr[3*i+1]=y[i];
																				srcPosarr[3*i+2]=z[i];

																			}				
																			//output_x_v();getchar();														


																			for(i=0;i<4*N;i++)
																				srcDenarr[i]=0.0;

																			//output_x_v();getchar();

																			for(i=0;i<N_ion;i++)
																				srcDenarr[4*i]=q[i];  //here the source or image point source have to divide by eps_s, because eps_s is already hidden in the harmonic coefficients.
																			for(ii=0;ii<M;ii++)
																			{
																				ind=4*(N_ion+ii);
																				for(k=0;k<p+1;k++)
																					for(j=p;j<2*p+1;j++)
																					{
																						if(abs(j-p)<=k) 
																						{
																							if(abs(j-p)==0)
																							{

																								srcDenarr[ind]=Bcs[ii][j][k];
																								ind=ind+1;
																							}
																							else if (abs(j-p)>0)
																							{
																								srcDenarr[ind]=Bcs[ii][j][k]; 
																								srcDenarr[ind+1]=Bcs[ii][2*p-j][k];
																								ind=ind+2;
																							}
																						}
																					}

																			}	

																			for(i=N_ion;i<N;i++)
																				srcDenarr[4*i]=q[i];
																			/*	for(i=0;i<M;i++)
																			for(j=0;j<2*p+1;j++)
																			for(k=0;k<p+1;k++)
																			{
																			if(abs(j-p)<=k)
																			{	
																			printf("%d %d %d %f \n",i,k,j-p,Bcs[i][j][k]);	
																			}															
																			}*/

																			/*	for(i=N_ion;i<N;i++)
																			printf("%d %f %f %f %f\n",i,srcDenarr[4*i],srcDenarr[4*i+1],srcDenarr[4*i+2],srcDenarr[4*i+3]);
																			getchar();*/

																			/*then direct compute the potential and force using Cartesian harmonic coefficients*/


																			for(i=0;i<3*N;i++)
																				pot_m[i]=0.0;


																			for (int t=0; t<N; t++)
																			{
																				double p[3]={0,0,0};
																				double tx=srcPosarr[3*t];
																				double ty=srcPosarr[3*t+1];       
																				double tz=srcPosarr[3*t+2];
																				//printf("tx=%.15f\n",tx);
																				double scale=-1; 
																				/*----------------------If up to p=6-------------------------------*/

																				int dims=4; //(p+1)*(p+1)
																				for (int s=N_ion; s<N; s++)
																				{
																					double dX_reg=srcPosarr[3*s+0]-tx;
																					double dY_reg=srcPosarr[3*s+1]-ty;
																					double dZ_reg=srcPosarr[3*s+2]-tz;
																					dX_reg=dX_reg/scale;
																					dY_reg=dY_reg/scale;
																					dZ_reg=dZ_reg/scale;
																					//double x=dX_reg;
																					//double y=dY_reg;
																					//double z=dZ_reg;
																					double x2=dX_reg*dX_reg;
																					double y2=dY_reg*dY_reg;
																					double z2=dZ_reg*dZ_reg;
																					double invR = (x2+y2+z2);
																					if (invR!=0)
																					{
																						invR = 1.0/sqrt(invR);
																						dr=1.0/invR;
																						double pm1c=dX_reg, pm1s=dY_reg, pm2c=x2-y2, pm2s=dX_reg*dY_reg,pm3c=dX_reg*(x2-3.0*y2),pm3s=dY_reg*(3.*x2-y2),pm4c=x2*x2+y2*y2-6.0*x2*y2,pm4s=pm2s*pm2c;
																						double pm5c=dX_reg*(x2*x2-10.0*x2*y2+5.0*y2*y2),pm5s=dY_reg*(y2*y2-10.0*x2*y2+5.0*x2*x2);
																						double pm6c=x2*x2*x2-15.0*x2*x2*y2+15.0*x2*y2*y2-y2*y2*y2;
																						double pm6s=pm2s*(6.0*x2*x2-20.0*x2*y2+6.0*y2*y2);

																						double cl2=1.732050807568877;  //sqrt(3)
																						double cl31=2.449489742783178; //sqrt(6)
																						double cl32=3.872983346207417; //sqrt(15)
																						double cl33=3.162277660168380; //sqrt(10)
																						double cl41=cl33;
																						double cl42=2.236067977499790; //sqrt(5)
																						double cl43=8.366600265340756; //sqrt(70)
																						double cl44=5.916079783099616; //sqrt(35)
																						double cl51=cl32;
																						double cl52=10.246950765959598;
																						double cl53=4.183300132670378;
																						double cl54=5.916079783099616;
																						double cl55=1.870828693386971;
																						double cl61=4.582575694955840;
																						double cl62=7.245688373094719;
																						double cl63=cl62;
																						double cl64=2.645751311064591;
																						double cl65=6.204836822995429;
																						double cl66=10.747092630102339;

																						double invR2=invR*invR;
																						double invR3=invR2*invR;
																						double invR5=invR3*invR2;
																						double invR7=invR5*invR2;
																						double invR9=invR7*invR2;
																						double invR11=invR9*invR2;
																						double invR13=invR11*invR2;
																						double invR15=invR13*invR2;

																						double phi = srcDenarr[dims*s]*invR;

																						phi+=srcDenarr[dims*s+1]*dZ_reg*invR3;
																						phi+=srcDenarr[dims*s+2]*pm1c*invR3;
																						phi+=srcDenarr[dims*s+3]*pm1s*invR3;	

																						/*	phi+=srcDenarr[dims*s+4]*0.5*(3.0*z2-1.0/invR2)*invR5;
																						phi+=srcDenarr[dims*s+5]*cl2*pm1c*dZ_reg*invR5;
																						phi+=srcDenarr[dims*s+6]*cl2*pm1s*dZ_reg*invR5;
																						phi+=srcDenarr[dims*s+7]*0.5*cl2*pm2c*invR5;
																						phi+=srcDenarr[dims*s+8]*cl2*pm2s*invR5;

																						phi+=srcDenarr[dims*s+9]*0.5*(5.0*z2*dZ_reg-3.0*dZ_reg/invR2)*invR7;
																						phi+=srcDenarr[dims*s+10]*0.25*cl31*pm1c*(5.0*z2-1.0/invR2)*invR7;
																						phi+=srcDenarr[dims*s+11]*0.25*cl31*pm1s*(5.0*z2-1.0/invR2)*invR7;
																						phi+=srcDenarr[dims*s+12]*0.5*cl32*dZ_reg*pm2c*invR7;
																						phi+=srcDenarr[dims*s+13]*cl32*dZ_reg*pm2s*invR7;
																						phi+=srcDenarr[dims*s+14]*0.25*cl33*pm3c*invR7;
																						phi+=srcDenarr[dims*s+15]*0.25*cl33*pm3s*invR7;

																						phi+=srcDenarr[dims*s+16]*0.125*(8.0*z2*z2-24.*(x2+y2)*z2+3.0*(x2*x2+2.0*x2*y2+y2*y2))*invR9;
																						phi+=srcDenarr[dims*s+17]*0.25*cl41*(4.0*dX_reg*dZ_reg*z2-3.0*dX_reg*dZ_reg*(x2+y2))*invR9;
																						phi+=srcDenarr[dims*s+18]*0.25*cl41*(4.0*dY_reg*dZ_reg*z2-3.0*dY_reg*dZ_reg*(x2+y2))*invR9;
																						phi+=srcDenarr[dims*s+19]*0.25*cl42*pm2c*(6.0*z2-x2-y2)*invR9;
																						phi+=srcDenarr[dims*s+20]*0.5*cl42*pm2s*(6.0*z2-x2-y2)*invR9;
																						phi+=srcDenarr[dims*s+21]*0.25*cl43*dZ_reg*pm3c*invR9;
																						phi+=srcDenarr[dims*s+22]*0.25*cl43*dZ_reg*pm3s*invR9;
																						phi+=srcDenarr[dims*s+23]*0.125*cl44*pm4c*invR9;
																						phi+=srcDenarr[dims*s+24]*0.5*cl44*pm4s*invR9;

																						phi+=srcDenarr[dims*s+25]*0.125*(63.0*z2*z2*dZ_reg-70.0*z2*dZ_reg/invR2+15.0*dZ_reg/invR3/invR)*invR11;
																						phi+=srcDenarr[dims*s+26]*0.125*cl51*pm1c*(1.0/invR3/invR-14.0*z2/invR2+21.0*z2*z2)*invR11;
																						phi+=srcDenarr[dims*s+27]*0.125*cl51*pm1s*(1.0/invR3/invR-14.0*z2/invR2+21.0*z2*z2)*invR11;
																						phi+=srcDenarr[dims*s+28]*0.125*cl52*pm2c*(6.0*z2*dZ_reg-2.0*dZ_reg/invR2)*invR11;
																						phi+=srcDenarr[dims*s+29]*0.25*cl52*pm2s*(6.0*z2*dZ_reg-2.0*dZ_reg/invR2)*invR11;
																						phi+=srcDenarr[dims*s+30]*0.0625*cl53*pm3c*(18.0*z2-2.0/invR2)*invR11;
																						phi+=srcDenarr[dims*s+31]*0.0625*cl53*pm3s*(18.0*z2-2.0/invR2)*invR11;
																						phi+=srcDenarr[dims*s+32]*0.375*cl54*pm4c*dZ_reg*invR11;
																						phi+=srcDenarr[dims*s+33]*1.5*cl54*pm4s*dZ_reg*invR11;
																						phi+=srcDenarr[dims*s+34]*0.375*cl55*pm5c*invR11;
																						phi+=srcDenarr[dims*s+35]*0.375*cl55*pm5s*invR11;

																						phi+=srcDenarr[dims*s+36]*0.0625*(231.0*z2*z2*z2-315.0*z2*z2/invR2+105.0*z2/invR3/invR-5.0/invR3/invR3)*invR13;
																						phi+=srcDenarr[dims*s+37]*0.125*cl61*pm1c*dZ_reg*(33.0*z2*z2-30.0*z2/invR2+5.0/invR2/invR2)*invR13;
																						phi+=srcDenarr[dims*s+38]*0.125*cl61*pm1s*dZ_reg*(33.0*z2*z2-30.0*z2/invR2+5.0/invR2/invR2)*invR13;
																						phi+=srcDenarr[dims*s+39]*0.0625*cl62*pm2c*(33.0*z2*z2-18.0*z2/invR2+1.0/invR2/invR2)*invR13;
																						phi+=srcDenarr[dims*s+40]*0.125*cl62*pm2s*(33.0*z2*z2-18.0*z2/invR2+1.0/invR2/invR2)*invR13;
																						phi+=srcDenarr[dims*s+41]*0.0625*cl63*pm3c*dZ_reg*(22.0*z2-6.0/invR2)*invR13;
																						phi+=srcDenarr[dims*s+42]*0.0625*cl63*pm3s*dZ_reg*(22.0*z2-6.0/invR2)*invR13;
																						phi+=srcDenarr[dims*s+43]*0.09375*cl64*pm4c*(22.0*z2-2.0/invR2)*invR13;
																						phi+=srcDenarr[dims*s+44]*0.375*cl64*pm4s*(22.0*z2-2.0/invR2)*invR13;
																						phi+=srcDenarr[dims*s+45]*0.375*cl65*pm5c*dZ_reg*invR13;
																						phi+=srcDenarr[dims*s+46]*0.375*cl65*pm5s*dZ_reg*invR13;
																						phi+=srcDenarr[dims*s+47]*0.0625*cl66*pm6c*invR13;
																						phi+=srcDenarr[dims*s+48]*0.0625*cl66*pm6s*invR13;*/

																						double fx =srcDenarr[dims*s]*dX_reg*invR3; 
																						fx+=srcDenarr[dims*s+1]*3.0*dX_reg*dZ_reg*invR5;
																						fx+=srcDenarr[dims*s+2]*(3.0*x2-1.0/invR2)*invR5;
																						fx+=srcDenarr[dims*s+3]*3.0*dX_reg*dY_reg*invR5;

																						/*fx+=srcDenarr[dims*s+4]*1.5*dX_reg*(5.*z2-1.0/invR2)*invR7;
																						fx+=srcDenarr[dims*s+5]*cl2*dZ_reg*(5.0*x2-1.0/invR2)*invR7;
																						fx+=srcDenarr[dims*s+6]*cl2*5.0*pm1c*pm1s*dZ_reg*invR7;
																						fx+=srcDenarr[dims*s+7]*0.5*cl2*dX_reg*(5.0*pm2c-2.0/invR2)*invR7;
																						fx+=srcDenarr[dims*s+8]*cl2*pm1s*(5.0*x2-1.0/invR2)*invR7;

																						fx+=srcDenarr[dims*s+9]*dX_reg*dZ_reg*(10.0*z2-7.5*x2-7.5*y2)*invR9;
																						fx+=srcDenarr[dims*s+10]*0.25*cl31*(-4*x2*x2 +y2*y2 - 3*y2*z2 - 4*z2*z2 - 3*x2*y2 - 9*z2)*invR9;
																						fx+=srcDenarr[dims*s+11]*0.25*cl31*(-5*pm2s*(x2 + y2 - 6*z2))*invR9;
																						fx+=srcDenarr[dims*s+12]*0.5*cl32*(dX_reg*dZ_reg*(5*x2 - 9*y2 - 2*z2))*invR9;
																						fx+=srcDenarr[dims*s+13]*cl32*(-1.0*(dY_reg*dZ_reg*(-6*x2 + y2 + z2)))*invR9;
																						fx+=srcDenarr[dims*s+14]*0.25*cl33*(4*x2*x2 + 3*y2*(y2 + z2) - 
																						3*x2*(7*y2 + z2))*invR9;
																						fx+=srcDenarr[dims*s+15]*0.25*cl33*(pm2s*(15*x2 - 13*y2 - 6*z2))*invR9;

																						fx+=srcDenarr[dims*s+16]*(dX_reg*(1.875*x2*x2 + 1.875*y2*y2 - 22.5*y2*z2 + 
																						15.*z2*z2 + x2*(3.75*y2 - 22.5*z2)))*invR11;
																						fx+=srcDenarr[dims*s+17]*0.25*cl41*(-1.0*z*(18*x2*x2 - 3*y2*y2 + y2*z2 + 
																						4*z2*z2 + x2*(15*y2 - 41*z2)))*invR11;
																						fx+=srcDenarr[dims*s+18]*0.25*cl41*(-21*x*y*z*(x2 + y2 - 2*z2))*invR11;
																						fx+=srcDenarr[dims*s+19]*0.25*cl42*(x*(-5*x2*x2 + 9*y2*y2 - 66*y2*z2 - 
																						12*z2*z2 + x2*(4*y2 + 46*z2)))*invR11;
																						fx+=srcDenarr[dims*s+20]*0.5*cl42*(y*(-6*x2*x2 + y2*y2 - 5*y2*z2 - 
																						6*z2*z2 + x2*(-5*y2 + 51*z2)))*invR11;
																						fx+=srcDenarr[dims*s+21]*0.25*cl43*(3*(2*x2*x2 + y2*(y2 + z2) - 
																						x2*(9*y2 + z2)))*invR11;
																						fx+=srcDenarr[dims*s+22]*0.25*cl43*(3*x*y*(7*x2 - 5*y2 - 2*z2))*invR11;
																						fx+=srcDenarr[dims*s+23]*0.125*cl44*(x*(5*x2*x2 - 2*x2*(23*y2 + 2*z2) + 
																						3*y2*(7*y2 + 4*z2)))*invR11;
																						fx+=srcDenarr[dims*s+24]*0.5*cl44*(y*(6*x2*x2 + y2*(y2 + z2) - 
																						x2*(11*y2 + 3*z2)))*invR11;

																						fx+=srcDenarr[dims*s+25]*(x*z*(13.125*x2*x2 + 13.125*y2*y2 - 
																						52.5*y2*z2 + 21.*z2*z2 + 
																						x2*(26.25*y2 - 52.5*z2)))*invR13;
																						fx+=srcDenarr[dims*s+26]*0.125*cl51*(6*x2*x2*x2 - y2*y2*y2 + 11*y2*y2*z2 + 
																						4*y2*z2*z2 - 8*z2*z2*z2 + 
																						x2*x2*(11*y2 - 101*z2) + 
																						2*x2*(2*y2*y2 - 45*y2*z2 + 
																						58*z2*z2))*invR13;
																						fx+=srcDenarr[dims*s+27]*0.125*cl51*(7*x*y*(x2*x2 + y2*y2 - 16*y2*z2 + 
																						16*z2*z2 + 2*x2*(y2 - 8*z2)))*invR13;
																						fx+=srcDenarr[dims*s+28]*0.125*cl52*(2*x*z*(-7*x2*x2 + 11*y2*y2 - 26*y2*z2 - 
																						4*z2*z2 + x2*(4*y2 + 22*z2)))*invR13;
																						fx+=srcDenarr[dims*s+29]*0.25*cl52*(-2*y*z*(8*x2*x2 - y2*y2 + y2*z2 + 
																						2*z2*z2 + x2*(7*y2 - 23*z2)))*invR13;
																						fx+=srcDenarr[dims*s+30]*0.0625*cl53*(-6*(2*x2*x2*x2 + y2*y2*y2 - 7*y2*y2*z2 - 
																						8*y2*z2*z2 - 
																						x2*x2*(7*y2 + 23*z2) + 
																						x2*(-8*y2*y2 + 90*y2*z2 + 
																						8*z2*z2)))*invR13;
																						fx+=srcDenarr[dims*s+31]*0.0625*cl53*(-6*x*y*(7*x2*x2 - 5*y2*y2 + 44*y2*z2 + 
																						16*z2*z2 + 2*x2*(y2 - 38*z2)))*invR13;
																						fx+=srcDenarr[dims*s+32]*0.375*cl54*(x*z*(7*x2*x2 + 23*y2*y2 + 12*y2*z2 - 
																						2*x2*(29*y2 + 2*z2)))*invR13;
																						fx+=srcDenarr[dims*s+33]*1.5*cl54*(y*z*(8*x2*x2 + y2*(y2 + z2) - 
																						x2*(13*y2 + 3*z2)))*invR13;
																						fx+=srcDenarr[dims*s+34]*0.375*cl55*(6*x2*x2*x2 - 5*y2*y2*(y2 + z2) - 
																						5*x2*x2*(17*y2 + z2) + 
																						10*x2*(8*y2*y2 + 3*y2*z2))*invR13;
																						fx+=srcDenarr[dims*s+35]*0.375*cl55*(x*y*(35*x2*x2 + 31*y2*y2 + 20*y2*z2 - 
																						10*x2*(11*y2 + 2*z2)))*invR13;

																						fx+=srcDenarr[dims*s+36]*(x*(-2.1875*x2*x2*x2 - 2.1875*y2*y2*y2 + 
																						52.5*y2*y2*z2 - 105.*y2*z2*z2 + 
																						28.*z2*z2*z2 + x2*x2*
																						(-6.5625*y2 + 52.5*z2) + 
																						x2*(-6.5625*y2*y2 + 105.*y2*z2 - 
																						105.*z2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+37]*0.125*cl61*(z*(40*x2*x2*x2 - 5*y2*y2*y2 + 15*y2*y2*z2 + 
																						12*y2*z2*z2 - 8*z2*z2*z2 + 
																						75*x2*x2*(y2 - 3*z2) + 
																						6*x2*(5*y2*y2 - 35*y2*z2 + 
																						26*z2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+38]*0.125*cl61*(3*x*y*z*(15*x2*x2 + 15*y2*y2 - 80*y2*z2 + 
																						48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;
																						fx+=srcDenarr[dims*s+39]*0.0625*cl62*(x*(7*x2*x2*x2 - 11*y2*y2*y2 + 210*y2*y2*z2 - 
																						240*y2*z2*z2 - 32*z2*z2*z2 + 
																						3*x2*x2*(y2 - 50*z2) - 
																						15*x2*(y2*y2 - 4*y2*z2 - 
																						16*z2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+40]*0.125*cl62*(y*(8*x2*x2*x2 - y2*y2*y2 + 15*y2*y2*z2 - 
																						16*z2*z2*z2 + 15*x2*x2*(y2 - 11*z2) + 
																						6*x2*(y2*y2 - 25*y2*z2 + 
																						40*z2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+41]*0.0625*cl63*(-2*z*(24*x2*x2*x2 - 5*x2*x2*(15*y2 + 19*z2) + 
																						3*y2*(3*y2*y2 - 5*y2*z2 - 
																						8*z2*z2) + x2*
																						(-90*y2*y2 + 330*y2*z2 + 24*z2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+42]*0.0625*cl63*(-2*x*y*z*(81*x2*x2 - 51*y2*y2 + 140*y2*z2 + 
																						48*z2*z2 + 30*x2*(y2 - 10*z2)))*invR15;
																						fx+=srcDenarr[dims*s+43]*0.09375*cl64*(-2*x*(7*x2*x2*x2 + 23*y2*y2*y2 - 240*y2*y2*z2 - 
																						120*y2*z2*z2 - 
																						3*x2*x2*(17*y2 + 32*z2) + 
																						x2*(-35*y2*y2 + 720*y2*z2 + 
																						40*z2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+44]*0.375*cl64*(-2*y*(8*x2*x2*x2 + y2*y2*y2 - 9*y2*y2*z2 - 
																						10*y2*z2*z2 - 
																						5*x2*x2*(y2 + 21*z2) - 
																						6*x2*(2*y2*y2 - 25*y2*z2 - 
																						5*z2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+45]*0.375*cl65*(z*(8*x2*x2*x2 - 5*y2*y2*(y2 + z2) + 
																						30*x2*y2*(3*y2 + z2) - 
																						5*x2*x2*(21*y2 + z2)))*invR15;
																						fx+=srcDenarr[dims*s+46]*0.375*cl65*(x*y*z*(45*x2*x2 + 33*y2*y2 + 20*y2*z2 - 
																						10*x2*(13*y2 + 2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+47]*0.0625*cl66*(x*(7*x2*x2*x2 - 43*y2*y2*y2 - 30*y2*y2*z2 - 
																						3*x2*x2*(47*y2 + 2*z2) + 
																						15*x2*(15*y2*y2 + 4*y2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+48]*0.0625*cl66*(2*y*(24*x2*x2*x2 - 3*y2*y2*(y2 + z2) - 
																						5*x2*x2*(23*y2 + 3*z2) + 
																						6*x2*(11*y2*y2 + 5*y2*z2)))*invR15;*/

																						double fy =srcDenarr[dims*s]*dY_reg*invR3;
																						fy+=srcDenarr[dims*s+1]*3.0*dY_reg*dZ_reg*invR5;
																						fy+=srcDenarr[dims*s+2]*3.0*dX_reg*dY_reg*invR5;
																						fy+=srcDenarr[dims*s+3]*(3.0*y2-1.0/invR2)*invR5;

																						/*	fy+=srcDenarr[dims*s+4]*1.5*dY_reg*(5.*z2-1.0/invR2)*invR7;
																						fy+=srcDenarr[dims*s+5]*cl2*5.0*pm1c*pm1s*dZ_reg*invR7;
																						fy+=srcDenarr[dims*s+6]*cl2*dZ_reg*(5.0*y2-1.0/invR2)*invR7;
																						fy+=srcDenarr[dims*s+7]*0.5*cl2*dY_reg*(5.0*pm2c-2.0/invR2)*invR7;
																						fy+=srcDenarr[dims*s+8]*cl2*pm1c*(5.0*y2-1.0/invR2)*invR7;

																						fy+=srcDenarr[dims*s+9]*dY_reg*dZ_reg*(10.0*z2-7.5*x2-7.5*y2)*invR9;
																						fy+=srcDenarr[dims*s+10]*0.25*cl31*(-5.0*pm2s*(x2 + y2 - 6*z2))*invR9;
																						fy+=srcDenarr[dims*s+11]*0.25*cl31*(x2*x2 - 4*y2*y2 + 27*y2*z2 - 4*z2*z2 - 
																						3*x2*(y2 + z2))*invR9;
																						fy+=srcDenarr[dims*s+12]*0.5*cl32*(dY_reg*dZ_reg*(9*x2 - 5*y2 + 2*z2))*invR9;
																						fy+=srcDenarr[dims*s+13]*cl32*(-(x*z*(x2 - 6*y2 + z2)))*invR9;
																						fy+=srcDenarr[dims*s+14]*0.25*cl33*(pm2s*(13*x2 - 15*y2 + 6*z2))*invR9;
																						fy+=srcDenarr[dims*s+15]*0.25*cl33*(-3*x2*x2 - 4*y2*y2 + 3*y2*z2 + 
																						3*x2*(7*y2 - z2))*invR9;

																						fy+=srcDenarr[dims*s+16]*(dY_reg*(1.875*x2*x2 + 1.875*y2*y2 - 22.5*y2*z2 + 
																						15.*z2*z2 + x2*(3.75*y2 - 22.5*z2)))*invR11;
																						fy+=srcDenarr[dims*s+17]*0.25*cl41*(-21*x*y*z*(x2 + y2 - 2*z2))*invR11;
																						fy+=srcDenarr[dims*s+18]*0.25*cl41*(-1.0*z*(-3*x2*x2 + 18*y2*y2 - 41*y2*z2 + 
																						4*z2*z2 + x2*(15*y2 + z2)))*invR11;
																						fy+=srcDenarr[dims*s+19]*0.25*cl42*(y*(-9*x2*x2 + 5*y2*y2 - 46*y2*z2 + 
																						12*z2*z2 + x2*(-4*y2 + 66*z2)))*invR11;
																						fy+=srcDenarr[dims*s+20]*0.5*cl42*(x*(x2*x2 - 6*y2*y2 + 51*y2*z2 - 
																						6*z2*z2 - 5*x2*(y2 + z2)))*invR11;
																						fy+=srcDenarr[dims*s+21]*0.25*cl43*(3*x*y*(5*x2 - 7*y2 + 2*z2))*invR11;
																						fy+=srcDenarr[dims*s+22]*0.25*cl43*(-3*(x2*x2 + 2*y2*y2 - y2*z2 + 
																						x2*(-9*y2 + z2)))*invR11;
																						fy+=srcDenarr[dims*s+23]*0.125*cl44*(y*(21*x2*x2 + 5*y2*y2 - 4*y2*z2 + 
																						x2*(-46*y2 + 12*z2)))*invR11;
																						fy+=srcDenarr[dims*s+24]*0.5*cl44*(-1.0*x*(x2*x2 + 6*y2*y2 - 3*y2*z2 + 
																						x2*(-11*y2 + z2)))*invR11;

																						fy+=srcDenarr[dims*s+25]*(y*z*(13.125*x2*x2 + 13.125*y2*y2 - 
																						52.5*y2*z2 + 21.*z2*z2 + 
																						x2*(26.25*y2 - 52.5*z2)))*invR13;
																						fy+=srcDenarr[dims*s+26]*0.125*cl51*(7*x*y*(x2*x2 + y2*y2 - 16*y2*z2 + 
																						16*z2*z2 + 2*x2*(y2 - 8*z2)))*invR13;
																						fy+=srcDenarr[dims*s+27]*0.125*cl51*(-x2*x2*x2 + 6*y2*y2*y2 - 101*y2*y2*z2 + 
																						116*y2*z2*z2 - 8*z2*z2*z2 + 
																						x2*x2*(4*y2 + 11*z2) + 
																						x2*(11*y2*y2 - 90*y2*z2 + 4*z2*z2))*invR13;
																						fy+=srcDenarr[dims*s+28]*0.125*cl52*(y*z*(-11*x2*x2 + 7*y2*y2 - 22*y2*z2 + 
																						4*z2*z2 + x2*(-4*y2 + 26*z2)))*invR13;
																						fy+=srcDenarr[dims*s+29]*0.25*cl52*(2*x*z*(x2*x2 - 8*y2*y2 + 23*y2*z2 - 
																						2*z2*z2 - x2*(7*y2 + z2)))*invR13;
																						fy+=srcDenarr[dims*s+30]*0.0625*cl53*(-6*x*y*(5*x2*x2 - 7*y2*y2 + 76*y2*z2 - 
																						16*z2*z2 - 2*x2*(y2 + 22*z2)))*invR13;
																						fy+=srcDenarr[dims*s+31]*0.0625*cl53*(6*(x2*x2*x2 + 2*y2*y2*y2 - 23*y2*y2*z2 + 
																						8*y2*z2*z2 - 
																						x2*x2*(8*y2 + 7*z2) + 
																						x2*(-7*y2*y2 + 90*y2*z2 - 
																						8*z2*z2)))*invR13;
																						fy+=srcDenarr[dims*s+32]*0.375*cl54*(y*z*(23*x2*x2 + 7*y2*y2 - 4*y2*z2 + 
																						x2*(-58*y2 + 12*z2)))*invR13;
																						fy+=srcDenarr[dims*s+33]*1.5*cl54*(-1.0*x*z*(x2*x2 + 8*y2*y2 - 3*y2*z2 + 
																						x2*(-13*y2 + z2)))*invR13;
																						fy+=srcDenarr[dims*s+34]*0.375*cl55*(x*y*(31*x2*x2 + 5*y2*(7*y2 - 4*z2) + 
																						x2*(-110*y2 + 20*z2)))*invR13;
																						fy+=srcDenarr[dims*s+35]*0.375*cl55*(-5*x2*x2*x2 + 6*y2*y2*y2 - 5*y2*y2*z2 + 
																						x2*x2*(80*y2 - 5*z2) + 
																						x2*(-85*y2*y2 + 30*y2*z2))*invR13;

																						fy+=srcDenarr[dims*s+36]*(y*(-2.1875*x2*x2*x2 - 2.1875*y2*y2*y2 + 
																						52.5*y2*y2*z2 - 105.*y2*z2*z2 + 
																						28.*z2*z2*z2 + x2*x2*
																						(-6.5625*y2 + 52.5*z2) + 
																						x2*(-6.5625*y2*y2 + 105.*y2*z2 - 
																						105.*z2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+37]*0.125*cl61*(3*x*y*z*(15*x2*x2 + 15*y2*y2 - 80*y2*z2 + 
																						48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;
																						fy+=srcDenarr[dims*s+38]*0.125*cl61*(z*(-5*x2*x2*x2 + 40*y2*y2*y2 - 225*y2*y2*z2 + 
																						156*y2*z2*z2 - 8*z2*z2*z2 + 
																						15*x2*x2*(2*y2 + z2) + 
																						3*x2*(25*y2*y2 - 70*y2*z2 + 
																						4*z2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+39]*0.0625*cl62*(y*(11*x2*x2*x2 - 7*y2*y2*y2 + 150*y2*y2*z2 - 
																						240*y2*z2*z2 + 32*z2*z2*z2 + 
																						15*x2*x2*(y2 - 14*z2) - 
																						3*x2*(y2*y2 + 20*y2*z2 - 
																						80*z2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+40]*0.125*cl62*(x*(-x2*x2*x2 + 8*y2*y2*y2 - 165*y2*y2*z2 + 
																						240*y2*z2*z2 - 16*z2*z2*z2 + 
																						3*x2*x2*(2*y2 + 5*z2) + 
																						15*x2*(y2*y2 - 10*y2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+41]*0.0625*cl63*(2*x*y*z*(-51*x2*x2 + 81*y2*y2 - 300*y2*z2 + 
																						48*z2*z2 + 10*x2*(3*y2 + 14*z2)))*invR15;
																						fy+=srcDenarr[dims*s+42]*0.0625*cl63*(2*z*(9*x2*x2*x2 + 24*y2*y2*y2 - 95*y2*y2*z2 + 
																						24*y2*z2*z2 - 
																						15*x2*x2*(6*y2 + z2) - 
																						3*x2*(25*y2*y2 - 110*y2*z2 + 
																						8*z2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+43]*0.09375*cl64*(-2*y*(23*x2*x2*x2 + 7*y2*y2*y2 - 96*y2*y2*z2 + 
																						40*y2*z2*z2 - 
																						5*x2*x2*(7*y2 + 48*z2) - 
																						3*x2*(17*y2*y2 - 240*y2*z2 + 
																						40*z2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+44]*0.375*cl64*(2*x*(x2*x2*x2 + 8*y2*y2*y2 - 105*y2*y2*z2 + 
																						30*y2*z2*z2 - 
																						3*x2*x2*(4*y2 + 3*z2) - 
																						5*x2*(y2*y2 - 30*y2*z2 + 
																						2*z2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+45]*0.375*cl65*(x*y*z*(33*x2*x2 + 5*y2*(9*y2 - 4*z2) + 
																						x2*(-130*y2 + 20*z2)))*invR15;
																						fy+=srcDenarr[dims*s+46]*0.375*cl65*(z*(-5*x2*x2*x2 + 8*y2*y2*y2 - 5*y2*y2*z2 + 
																						x2*x2*(90*y2 - 5*z2) - 
																						15*x2*(7*y2*y2 - 2*y2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+47]*0.0625*cl66*(y*(43*x2*x2*x2 - 7*y2*y2*y2 + 6*y2*y2*z2 + 
																						x2*x2*(-225*y2 + 30*z2) + 
																						3*x2*(47*y2*y2 - 20*y2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+48]*0.0625*cl66*(-2*x*(3*x2*x2*x2 - 24*y2*y2*y2 + 15*y2*y2*z2 + 
																						x2*x2*(-66*y2 + 3*z2) + 
																						5*x2*(23*y2*y2 - 6*y2*z2)))*invR15;*/

																						double fz =srcDenarr[dims*s]*dZ_reg*invR3;
																						fz+=srcDenarr[dims*s+1]*(3.0*z2-1.0/invR2)*invR5;
																						fz+=srcDenarr[dims*s+2]*3.0*dX_reg*dZ_reg*invR5;
																						fz+=srcDenarr[dims*s+3]*3.0*dY_reg*dZ_reg*invR5;

																						/*	fz+=srcDenarr[dims*s+4]*1.5*dZ_reg*(5.*z2-3.0/invR2)*invR7;
																						fz+=srcDenarr[dims*s+5]*cl2*dX_reg*(5.0*z2-1.0/invR2)*invR7;
																						fz+=srcDenarr[dims*s+6]*cl2*dY_reg*(5.0*z2-1.0/invR2)*invR7;
																						fz+=srcDenarr[dims*s+7]*0.5*cl2*5.0*dZ_reg*pm2c*invR7;
																						fz+=srcDenarr[dims*s+8]*cl2*5.0*pm1c*pm1s*dZ_reg*invR7;

																						fz+=srcDenarr[dims*s+9]*(1.5*x2*x2+1.5*y2*y2-12.*y2*z2+4.*z2*z2+z2*(3.*y2-12.*z2))*invR9;
																						fz+=srcDenarr[dims*s+10]*0.25*cl31*(-5.0*dX_reg*dZ_reg*(3*x2 + 3*y2 - 4*x2))*invR9;
																						fz+=srcDenarr[dims*s+11]*0.25*cl31*(-5*dY_reg*dZ_reg*(3*x2 + 3*y2 - 4*z2))*invR9;
																						fz+=srcDenarr[dims*s+12]*0.5*cl32*((-x2 + y2)*(x2 + y2 - 6*z2))*invR9;
																						fz+=srcDenarr[dims*s+13]*cl32*(-(x*y*(x2 + y2 - 6*z2)))*invR9;
																						fz+=srcDenarr[dims*s+14]*0.25*cl33*(7*dX_reg*dZ_reg*(x2 - 3*y2))*invR9;
																						fz+=srcDenarr[dims*s+15]*0.25*cl33*(7*dY_reg*dZ_reg*(3*x2 - y2))*invR9;

																						fz+=srcDenarr[dims*s+16]*(dZ_reg*(9.375*x2*x2 + 9.375*y2*y2 - 25.*y2*z2 + 
																						5.*z2*z2 + x2*(18.75*y2 - 25.*z2)))*invR11;
																						fz+=srcDenarr[dims*s+17]*0.25*cl41*(3*x*(x2*x2 + y2*y2 - 12*y2*z2 + 
																						8*z2*z2 + 2*x2*(y2 - 6*z2)))*invR11;
																						fz+=srcDenarr[dims*s+18]*0.25*cl41*(3*y*(x2*x2 + y2*y2 - 12*y2*z2 + 
																						8*z2*z2 + 2*x2*(y2 - 6*z2)))*invR11;
																						fz+=srcDenarr[dims*s+19]*0.25*cl42*(-21*(x2 - y2)*z*(x2 + y2 - 2*z2))*invR11;
																						fz+=srcDenarr[dims*s+20]*0.5*cl42*(-21*x*y*z*(x2 + y2 - 2*z2))*invR11;
																						fz+=srcDenarr[dims*s+21]*0.25*cl43*dZ_reg*(9*(x2 - 3*y2)*z*x)*invR11;
																						fz+=srcDenarr[dims*s+22]*0.25*cl43*dZ_reg*(9*(3*x2 - y2)*z*y)*invR11;
																						fz+=srcDenarr[dims*s+23]*0.125*cl44*(9*(x2*x2 - 6*x2*y2 + y2*y2)*z)*invR11;
																						fz+=srcDenarr[dims*s+24]*0.5*cl44*(9*x*y*(x2 - y2)*z)*invR11;

																						fz+=srcDenarr[dims*s+25]*(-1.875*x2*x2*x2 - 1.875*y2*y2*y2 + 33.75*y2*y2*z2 - 
																						45.*y2*z2*z2 + 6.*z2*z2*z2 + 
																						x2*x2*(-5.625*y2 + 33.75*z2) + 
																						x2*(-5.625*y2*y2 + 67.5*y2*z2 - 
																						45.*z2*z2))*invR13;
																						fz+=srcDenarr[dims*s+26]*0.125*cl51*(7*x*z*(5*x2*x2 + 5*y2*y2 - 20*y2*z2 + 
																						8*z2*z2 + 10*x2*(y2 - 2*z2)))*invR13;
																						fz+=srcDenarr[dims*s+27]*0.125*cl51*(7*y*z*(5*x2*x2 + 5*y2*y2 - 20*y2*z2 + 
																						8*z2*z2 + 10*x2*(y2 - 2*z2)))*invR13;
																						fz+=srcDenarr[dims*s+28]*0.125*cl52*(2*(x2 - y2)*(x2*x2 + y2*y2 - 
																						16*y2*z2 + 16*z2*z2 + 
																						2*x2*(y2 - 8*z2)))*invR13;
																						fz+=srcDenarr[dims*s+29]*0.25*cl52*(2*x*y*(x2*x2 + y2*y2 - 16*y2*z2 + 
																						16*z2*z2 + 2*x2*(y2 - 8*z2)))*invR13;
																						fz+=srcDenarr[dims*s+30]*0.0625*cl53*(-18*x*(x2 - 3*y2)*z*
																						(3*x2 + 3*y2 - 8*z2))*invR13;
																						fz+=srcDenarr[dims*s+31]*0.0625*cl53*(18*y*(-3*x2 + y2)*z*
																						(3*x2 + 3*y2 - 8*z2))*invR13;
																						fz+=srcDenarr[dims*s+32]*0.375*cl54*(-1.0*(x2*x2 - 6*x2*y2 + y2*y2)*
																						(x2 + y2 - 10*z2))*invR13;
																						fz+=srcDenarr[dims*s+33]*1.5*cl54*(-1.0*x*y*(x2 - y2)*(x2 + y2 - 10*z2))*invR13;
																						fz+=srcDenarr[dims*s+34]*0.375*cl55*(11*x*(x2*x2 - 10*x2*y2 + 5*y2*y2)*z)*invR13;
																						fz+=srcDenarr[dims*s+35]*0.375*cl55*(-5*x2*x2*x2 + 6*y2*y2*y2 - 5*y2*y2*z2 + 
																						x2*x2*(80*y2 - 5*z2) + 
																						x2*(-85*y2*y2 + 30*y2*z2))*invR13;

																						fz+=srcDenarr[dims*s+36]*(z*(-15.3125*x2*x2*x2 - 15.3125*y2*y2*y2 + 
																						91.875*y2*y2*z2 - 73.5*y2*z2*z2 + 
																						7.*z2*z2*z2 + x2*x2*
																						(-45.9375*y2 + 91.875*z2) + 
																						x2*(-45.9375*y2*y2 + 183.75*y2*z2 - 
																						73.5*z2*z2)))*invR15;
																						fz+=srcDenarr[dims*s+37]*0.125*cl61*(3*x*y*z*(15*x2*x2 + 15*y2*y2 - 80*y2*z2 + 
																						48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;
																						fz+=srcDenarr[dims*s+38]*0.125*cl61*(-1.0*y*(5*x2*x2*x2 + 5*y2*y2*y2 - 120*y2*y2*z2 + 
																						240*y2*z2*z2 - 64*z2*z2*z2 + 
																						15*x2*x2*(y2 - 8*z2) + 
																						15*x2*(y2*y2 - 16*y2*z2 + 
																						16*z2*z2)))*invR15;
																						fz+=srcDenarr[dims*s+39]*0.0625*cl62*(3*(x2 - y2)*z*
																						(15*x2*x2 + 15*y2*y2 - 80*y2*z2 + 
																						48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;
																						fz+=srcDenarr[dims*s+40]*0.125*cl62*(3*x*y*z*(15*x2*x2 + 15*y2*y2 - 80*y2*z2 + 
																						48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;
																						fz+=srcDenarr[dims*s+41]*0.0625*cl63*(2*x*(x2 - 3*y2)*
																						(3*x2*x2 + 3*y2*y2 - 60*y2*z2 + 
																						80*z2*z2 + 6*x2*(y2 - 10*z2)))*invR15;
																						fz+=srcDenarr[dims*s+42]*0.0625*cl63*(-2*y*(-3*x2 + y2)*
																						(3*x2*x2 + 3*y2*y2 - 60*y2*z2 + 
																						80*z2*z2 + 6*x2*(y2 - 10*z2)))*invR15;
																						fz+=srcDenarr[dims*s+43]*0.09375*cl64*(-22*(x2*x2 - 6*x2*y2 + y2*y2)*z*
																						(3*x2 + 3*y2 - 10*z2))*invR15;
																						fz+=srcDenarr[dims*s+44]*0.375*cl64*(-22*x*y*(x2 - y2)*z*
																						(3*x2 + 3*y2 - 10*z2))*invR15;
																						fz+=srcDenarr[dims*s+45]*0.375*cl65*(-1.0*x*(x2*x2 - 10*x2*y2 + 5*y2*y2)*
																						(x2 + y2 - 12*z2))*invR15;
																						fz+=srcDenarr[dims*s+46]*0.375*cl65*(-1.0*y*(5*x2*x2 - 10*x2*y2 + y2*y2)*
																						(x2 + y2 - 12*z2))*invR15;
																						fz+=srcDenarr[dims*s+47]*0.0625*cl66*(13*(x2*x2*x2 - 15*x2*x2*y2 + 15*x2*y2*y2 - 
																						y2*y2*y2)*z)*invR15;
																						fz+=srcDenarr[dims*s+48]*0.0625*cl66*(26*x*y*(3*x2*x2 - 10*x2*y2 + 3*y2*y2)*z)*invR15;*/

																						//p[0] += phi;
																						if(t<N_ion)
																						{
																							p[0] += fx*srcDenarr[dims*t]; 
																							p[1] += fy*srcDenarr[dims*t];
																							p[2] += fz*srcDenarr[dims*t];

																							//	pot_m[4*s] -= phi;
																							pot_m[3*s+0] -= fx*srcDenarr[dims*t]; 
																							pot_m[3*s+1] -= fy*srcDenarr[dims*t];
																							pot_m[3*s+2] -= fz*srcDenarr[dims*t];

																							/*add correction*/
																							/*	if(dr<corrpara)
																							{
																							for(int c=N_ion;c<N;c++)
																							{
																							dx=x[c]-x[t];
																							dy=y[c]-y[t];
																							dz=z[c]-z[t];

																							double dr2=dx*dx+dy*dy+dz*dz;
																							if(dr2<corrpara*corrpara&&c!=s)
																							{
																							pot_m[3*c+0] =pot_m[3*c+0]+ fx*srcDenarr[dims*t]-srcDenarr[dims*t]*srcDenarr[dims*s]*(srcPosarr[3*t+0]-srcPosarr[3*s+0])/dr/dr/dr; 
																							pot_m[3*c+1] =pot_m[3*c+0]+ fy*srcDenarr[dims*t]-srcDenarr[dims*t]*srcDenarr[dims*s]*(srcPosarr[3*t+1]-srcPosarr[3*s+1])/dr/dr/dr; 
																							pot_m[3*c+2] =pot_m[3*c+0]+ fz*srcDenarr[dims*t]-srcDenarr[dims*t]*srcDenarr[dims*s]*(srcPosarr[3*t+2]-srcPosarr[3*s+2])/dr/dr/dr;;

																							}

																							}
																							}*/


																						}
																						else
																						{
																							p[0] += 0.5*fx*srcDenarr[dims*t]; 
																							p[1] += 0.5*fy*srcDenarr[dims*t];
																							p[2] += 0.5*fz*srcDenarr[dims*t];

																							//	pot_m[4*s] -= phi;
																							pot_m[3*s+0] -= 0.5*fx*srcDenarr[dims*t]; 
																							pot_m[3*s+1] -= 0.5*fy*srcDenarr[dims*t];
																							pot_m[3*s+2] -= 0.5*fz*srcDenarr[dims*t];
																						}



																					}
																				}
																				//	pot_m[4*t] += p[0];
																				pot_m[3*t+0] += p[0]; 
																				pot_m[3*t+1] += p[1];
																				pot_m[3*t+2] += p[2];

																			}



																			double temfx,temfy,temfz;
																			/*image-ion interaction*/
																			for(i=0;i<imcount;i++)
																				for(ii=0;ii<N_ion;ii++)
																				{
																					dx=-imx[i]+x[ii];
																					dy=-imy[i]+y[ii];
																					dz=-imz[i]+z[ii];
																					dr=sqrt(dx*dx+dy*dy+dz*dz);

																					temfx=q[ii]*imq[i]*dx/dr/dr/dr/epsi_s;
																					temfy=q[ii]*imq[i]*dy/dr/dr/dr/epsi_s;
																					temfz=q[ii]*imq[i]*dz/dr/dr/dr/epsi_s;


																					pot_m[3*ii+0]+= temfx;
																					pot_m[3*ii+1]+= temfy;
																					pot_m[3*ii+2]+= temfz;


																					pot_m[3*imind[i]+0]-=temfx;
																					pot_m[3*imind[i]+1]-=temfy;
																					pot_m[3*imind[i]+2]-=temfz;																																																						
																				}





																				/*image -sphere interaction*/



																		/*		for(i=0;i<imcount;i++)
																					for(ii=N_ion;ii<N;ii++)
																					{
																						dx=-imx[i]+x[ii];
																						dy=-imy[i]+y[ii];
																						dz=-imz[i]+z[ii];
																						dr=sqrt(dx*dx+dy*dy+dz*dz);

																						if(dr>1.0)
																						{
																							if(imind[i]<N_ion)
																							{
																								temfx=q[ii]*imq[i]*dx/dr/dr/dr/epsi_s;
																								temfy=q[ii]*imq[i]*dy/dr/dr/dr/epsi_s;
																								temfz=q[ii]*imq[i]*dz/dr/dr/dr/epsi_s;
																							}
																							else
																							{
																								temfx=q[ii]*imq[i]*dx/dr/dr/dr/epsi_s;
																								temfy=q[ii]*imq[i]*dy/dr/dr/dr/epsi_s;
																								temfz=q[ii]*imq[i]*dz/dr/dr/dr/epsi_s;
																							}

																							pot_m[3*ii+0]+=temfx;
																							pot_m[3*ii+1]+=temfy;
																							pot_m[3*ii+2]+=temfz;

																							pot_m[3*imind[i]+0]-=temfx;
																							pot_m[3*imind[i]+1]-=temfy;
																							pot_m[3*imind[i]+2]-=temfz;
																						}		
																																			
																					}*/

																					/*image-image interaction*/
																					for(i=0;i<imcount;i++)
																						for(ii=0;ii<imcount;ii++)
																						{
																							dx=-imx[i]+imx[ii];
																							dy=-imy[i]+imy[ii];
																							dz=-imz[i]+imz[ii];
																							dr=sqrt(dx*dx+dy*dy+dz*dz);

																							
																								if(imind[i]!=imind[ii])
																								{
																									temfx=imq[ii]*imq[i]*dx/dr/dr/dr/epsi_s;
																									temfy=imq[ii]*imq[i]*dy/dr/dr/dr/epsi_s;
																									temfz=imq[ii]*imq[i]*dz/dr/dr/dr/epsi_s;

																									pot_m[3*imind[ii]+0]+=temfx;
																									pot_m[3*imind[ii]+1]+=temfy;
																									pot_m[3*imind[ii]+2]+=temfz;

																									pot_m[3*imind[i]+0]-=temfx;
																									pot_m[3*imind[i]+1]-=temfy;
																									pot_m[3*imind[i]+2]-=temfz;
																								}													
																							}		
																						
																						/*image-multipole & image-sphere  interaction*/
																						for (int t=0; t<imcount; t++)
																			{
																				double p[3]={0,0,0};
																				double tx=imx[t];
																				double ty=imy[t];       
																				double tz=imz[t];
																				//printf("tx=%.15f\n",tx);
																				double scale=-1; 
																				/*----------------------If up to p=6-------------------------------*/

																				int dims=4;
																				for (int s=N_ion; s<N; s++)
																				{
																					double dX_reg=srcPosarr[3*s+0]-tx;
																					double dY_reg=srcPosarr[3*s+1]-ty;
																					double dZ_reg=srcPosarr[3*s+2]-tz;
																					dX_reg=dX_reg/scale;
																					dY_reg=dY_reg/scale;
																					dZ_reg=dZ_reg/scale;
																					//double x=dX_reg;
																					//double y=dY_reg;
																					//double z=dZ_reg;
																					double x2=dX_reg*dX_reg;
																					double y2=dY_reg*dY_reg;
																					double z2=dZ_reg*dZ_reg;
																					double invR = (x2+y2+z2);
																					if (invR!=0.0&&imind[t]!=s)
																					{
																						invR = 1.0/sqrt(invR);
																						dr=1.0/invR;
																						double pm1c=dX_reg, pm1s=dY_reg, pm2c=x2-y2, pm2s=dX_reg*dY_reg,pm3c=dX_reg*(x2-3.0*y2),pm3s=dY_reg*(3.*x2-y2),pm4c=x2*x2+y2*y2-6.0*x2*y2,pm4s=pm2s*pm2c;
																						double pm5c=dX_reg*(x2*x2-10.0*x2*y2+5.0*y2*y2),pm5s=dY_reg*(y2*y2-10.0*x2*y2+5.0*x2*x2);
																						double pm6c=x2*x2*x2-15.0*x2*x2*y2+15.0*x2*y2*y2-y2*y2*y2;
																						double pm6s=pm2s*(6.0*x2*x2-20.0*x2*y2+6.0*y2*y2);

																						double cl2=1.732050807568877;  //sqrt(3)
																						double cl31=2.449489742783178; //sqrt(6)
																						double cl32=3.872983346207417; //sqrt(15)
																						double cl33=3.162277660168380; //sqrt(10)
																						double cl41=cl33;
																						double cl42=2.236067977499790; //sqrt(5)
																						double cl43=8.366600265340756; //sqrt(70)
																						double cl44=5.916079783099616; //sqrt(35)
																						double cl51=cl32;
																						double cl52=10.246950765959598;
																						double cl53=4.183300132670378;
																						double cl54=5.916079783099616;
																						double cl55=1.870828693386971;
																						double cl61=4.582575694955840;
																						double cl62=7.245688373094719;
																						double cl63=cl62;
																						double cl64=2.645751311064591;
																						double cl65=6.204836822995429;
																						double cl66=10.747092630102339;

																						double invR2=invR*invR;
																						double invR3=invR2*invR;
																						double invR5=invR3*invR2;
																						double invR7=invR5*invR2;
																						double invR9=invR7*invR2;
																						double invR11=invR9*invR2;
																						double invR13=invR11*invR2;
																						double invR15=invR13*invR2;

																						double phi = srcDenarr[dims*s]*invR;

																						phi+=srcDenarr[dims*s+1]*dZ_reg*invR3;
																						phi+=srcDenarr[dims*s+2]*pm1c*invR3;
																						phi+=srcDenarr[dims*s+3]*pm1s*invR3;	

																						/*	phi+=srcDenarr[dims*s+4]*0.5*(3.0*z2-1.0/invR2)*invR5;
																						phi+=srcDenarr[dims*s+5]*cl2*pm1c*dZ_reg*invR5;
																						phi+=srcDenarr[dims*s+6]*cl2*pm1s*dZ_reg*invR5;
																						phi+=srcDenarr[dims*s+7]*0.5*cl2*pm2c*invR5;
																						phi+=srcDenarr[dims*s+8]*cl2*pm2s*invR5;

																						phi+=srcDenarr[dims*s+9]*0.5*(5.0*z2*dZ_reg-3.0*dZ_reg/invR2)*invR7;
																						phi+=srcDenarr[dims*s+10]*0.25*cl31*pm1c*(5.0*z2-1.0/invR2)*invR7;
																						phi+=srcDenarr[dims*s+11]*0.25*cl31*pm1s*(5.0*z2-1.0/invR2)*invR7;
																						phi+=srcDenarr[dims*s+12]*0.5*cl32*dZ_reg*pm2c*invR7;
																						phi+=srcDenarr[dims*s+13]*cl32*dZ_reg*pm2s*invR7;
																						phi+=srcDenarr[dims*s+14]*0.25*cl33*pm3c*invR7;
																						phi+=srcDenarr[dims*s+15]*0.25*cl33*pm3s*invR7;

																						phi+=srcDenarr[dims*s+16]*0.125*(8.0*z2*z2-24.*(x2+y2)*z2+3.0*(x2*x2+2.0*x2*y2+y2*y2))*invR9;
																						phi+=srcDenarr[dims*s+17]*0.25*cl41*(4.0*dX_reg*dZ_reg*z2-3.0*dX_reg*dZ_reg*(x2+y2))*invR9;
																						phi+=srcDenarr[dims*s+18]*0.25*cl41*(4.0*dY_reg*dZ_reg*z2-3.0*dY_reg*dZ_reg*(x2+y2))*invR9;
																						phi+=srcDenarr[dims*s+19]*0.25*cl42*pm2c*(6.0*z2-x2-y2)*invR9;
																						phi+=srcDenarr[dims*s+20]*0.5*cl42*pm2s*(6.0*z2-x2-y2)*invR9;
																						phi+=srcDenarr[dims*s+21]*0.25*cl43*dZ_reg*pm3c*invR9;
																						phi+=srcDenarr[dims*s+22]*0.25*cl43*dZ_reg*pm3s*invR9;
																						phi+=srcDenarr[dims*s+23]*0.125*cl44*pm4c*invR9;
																						phi+=srcDenarr[dims*s+24]*0.5*cl44*pm4s*invR9;

																						phi+=srcDenarr[dims*s+25]*0.125*(63.0*z2*z2*dZ_reg-70.0*z2*dZ_reg/invR2+15.0*dZ_reg/invR3/invR)*invR11;
																						phi+=srcDenarr[dims*s+26]*0.125*cl51*pm1c*(1.0/invR3/invR-14.0*z2/invR2+21.0*z2*z2)*invR11;
																						phi+=srcDenarr[dims*s+27]*0.125*cl51*pm1s*(1.0/invR3/invR-14.0*z2/invR2+21.0*z2*z2)*invR11;
																						phi+=srcDenarr[dims*s+28]*0.125*cl52*pm2c*(6.0*z2*dZ_reg-2.0*dZ_reg/invR2)*invR11;
																						phi+=srcDenarr[dims*s+29]*0.25*cl52*pm2s*(6.0*z2*dZ_reg-2.0*dZ_reg/invR2)*invR11;
																						phi+=srcDenarr[dims*s+30]*0.0625*cl53*pm3c*(18.0*z2-2.0/invR2)*invR11;
																						phi+=srcDenarr[dims*s+31]*0.0625*cl53*pm3s*(18.0*z2-2.0/invR2)*invR11;
																						phi+=srcDenarr[dims*s+32]*0.375*cl54*pm4c*dZ_reg*invR11;
																						phi+=srcDenarr[dims*s+33]*1.5*cl54*pm4s*dZ_reg*invR11;
																						phi+=srcDenarr[dims*s+34]*0.375*cl55*pm5c*invR11;
																						phi+=srcDenarr[dims*s+35]*0.375*cl55*pm5s*invR11;

																						phi+=srcDenarr[dims*s+36]*0.0625*(231.0*z2*z2*z2-315.0*z2*z2/invR2+105.0*z2/invR3/invR-5.0/invR3/invR3)*invR13;
																						phi+=srcDenarr[dims*s+37]*0.125*cl61*pm1c*dZ_reg*(33.0*z2*z2-30.0*z2/invR2+5.0/invR2/invR2)*invR13;
																						phi+=srcDenarr[dims*s+38]*0.125*cl61*pm1s*dZ_reg*(33.0*z2*z2-30.0*z2/invR2+5.0/invR2/invR2)*invR13;
																						phi+=srcDenarr[dims*s+39]*0.0625*cl62*pm2c*(33.0*z2*z2-18.0*z2/invR2+1.0/invR2/invR2)*invR13;
																						phi+=srcDenarr[dims*s+40]*0.125*cl62*pm2s*(33.0*z2*z2-18.0*z2/invR2+1.0/invR2/invR2)*invR13;
																						phi+=srcDenarr[dims*s+41]*0.0625*cl63*pm3c*dZ_reg*(22.0*z2-6.0/invR2)*invR13;
																						phi+=srcDenarr[dims*s+42]*0.0625*cl63*pm3s*dZ_reg*(22.0*z2-6.0/invR2)*invR13;
																						phi+=srcDenarr[dims*s+43]*0.09375*cl64*pm4c*(22.0*z2-2.0/invR2)*invR13;
																						phi+=srcDenarr[dims*s+44]*0.375*cl64*pm4s*(22.0*z2-2.0/invR2)*invR13;
																						phi+=srcDenarr[dims*s+45]*0.375*cl65*pm5c*dZ_reg*invR13;
																						phi+=srcDenarr[dims*s+46]*0.375*cl65*pm5s*dZ_reg*invR13;
																						phi+=srcDenarr[dims*s+47]*0.0625*cl66*pm6c*invR13;
																						phi+=srcDenarr[dims*s+48]*0.0625*cl66*pm6s*invR13;*/

																						double fx =srcDenarr[dims*s]*dX_reg*invR3; 
																						fx+=srcDenarr[dims*s+1]*3.0*dX_reg*dZ_reg*invR5;
																						fx+=srcDenarr[dims*s+2]*(3.0*x2-1.0/invR2)*invR5;
																						fx+=srcDenarr[dims*s+3]*3.0*dX_reg*dY_reg*invR5;

																						/*fx+=srcDenarr[dims*s+4]*1.5*dX_reg*(5.*z2-1.0/invR2)*invR7;
																						fx+=srcDenarr[dims*s+5]*cl2*dZ_reg*(5.0*x2-1.0/invR2)*invR7;
																						fx+=srcDenarr[dims*s+6]*cl2*5.0*pm1c*pm1s*dZ_reg*invR7;
																						fx+=srcDenarr[dims*s+7]*0.5*cl2*dX_reg*(5.0*pm2c-2.0/invR2)*invR7;
																						fx+=srcDenarr[dims*s+8]*cl2*pm1s*(5.0*x2-1.0/invR2)*invR7;

																						fx+=srcDenarr[dims*s+9]*dX_reg*dZ_reg*(10.0*z2-7.5*x2-7.5*y2)*invR9;
																						fx+=srcDenarr[dims*s+10]*0.25*cl31*(-4*x2*x2 +y2*y2 - 3*y2*z2 - 4*z2*z2 - 3*x2*y2 - 9*z2)*invR9;
																						fx+=srcDenarr[dims*s+11]*0.25*cl31*(-5*pm2s*(x2 + y2 - 6*z2))*invR9;
																						fx+=srcDenarr[dims*s+12]*0.5*cl32*(dX_reg*dZ_reg*(5*x2 - 9*y2 - 2*z2))*invR9;
																						fx+=srcDenarr[dims*s+13]*cl32*(-1.0*(dY_reg*dZ_reg*(-6*x2 + y2 + z2)))*invR9;
																						fx+=srcDenarr[dims*s+14]*0.25*cl33*(4*x2*x2 + 3*y2*(y2 + z2) - 
																						3*x2*(7*y2 + z2))*invR9;
																						fx+=srcDenarr[dims*s+15]*0.25*cl33*(pm2s*(15*x2 - 13*y2 - 6*z2))*invR9;

																						fx+=srcDenarr[dims*s+16]*(dX_reg*(1.875*x2*x2 + 1.875*y2*y2 - 22.5*y2*z2 + 
																						15.*z2*z2 + x2*(3.75*y2 - 22.5*z2)))*invR11;
																						fx+=srcDenarr[dims*s+17]*0.25*cl41*(-1.0*z*(18*x2*x2 - 3*y2*y2 + y2*z2 + 
																						4*z2*z2 + x2*(15*y2 - 41*z2)))*invR11;
																						fx+=srcDenarr[dims*s+18]*0.25*cl41*(-21*x*y*z*(x2 + y2 - 2*z2))*invR11;
																						fx+=srcDenarr[dims*s+19]*0.25*cl42*(x*(-5*x2*x2 + 9*y2*y2 - 66*y2*z2 - 
																						12*z2*z2 + x2*(4*y2 + 46*z2)))*invR11;
																						fx+=srcDenarr[dims*s+20]*0.5*cl42*(y*(-6*x2*x2 + y2*y2 - 5*y2*z2 - 
																						6*z2*z2 + x2*(-5*y2 + 51*z2)))*invR11;
																						fx+=srcDenarr[dims*s+21]*0.25*cl43*(3*(2*x2*x2 + y2*(y2 + z2) - 
																						x2*(9*y2 + z2)))*invR11;
																						fx+=srcDenarr[dims*s+22]*0.25*cl43*(3*x*y*(7*x2 - 5*y2 - 2*z2))*invR11;
																						fx+=srcDenarr[dims*s+23]*0.125*cl44*(x*(5*x2*x2 - 2*x2*(23*y2 + 2*z2) + 
																						3*y2*(7*y2 + 4*z2)))*invR11;
																						fx+=srcDenarr[dims*s+24]*0.5*cl44*(y*(6*x2*x2 + y2*(y2 + z2) - 
																						x2*(11*y2 + 3*z2)))*invR11;

																						fx+=srcDenarr[dims*s+25]*(x*z*(13.125*x2*x2 + 13.125*y2*y2 - 
																						52.5*y2*z2 + 21.*z2*z2 + 
																						x2*(26.25*y2 - 52.5*z2)))*invR13;
																						fx+=srcDenarr[dims*s+26]*0.125*cl51*(6*x2*x2*x2 - y2*y2*y2 + 11*y2*y2*z2 + 
																						4*y2*z2*z2 - 8*z2*z2*z2 + 
																						x2*x2*(11*y2 - 101*z2) + 
																						2*x2*(2*y2*y2 - 45*y2*z2 + 
																						58*z2*z2))*invR13;
																						fx+=srcDenarr[dims*s+27]*0.125*cl51*(7*x*y*(x2*x2 + y2*y2 - 16*y2*z2 + 
																						16*z2*z2 + 2*x2*(y2 - 8*z2)))*invR13;
																						fx+=srcDenarr[dims*s+28]*0.125*cl52*(2*x*z*(-7*x2*x2 + 11*y2*y2 - 26*y2*z2 - 
																						4*z2*z2 + x2*(4*y2 + 22*z2)))*invR13;
																						fx+=srcDenarr[dims*s+29]*0.25*cl52*(-2*y*z*(8*x2*x2 - y2*y2 + y2*z2 + 
																						2*z2*z2 + x2*(7*y2 - 23*z2)))*invR13;
																						fx+=srcDenarr[dims*s+30]*0.0625*cl53*(-6*(2*x2*x2*x2 + y2*y2*y2 - 7*y2*y2*z2 - 
																						8*y2*z2*z2 - 
																						x2*x2*(7*y2 + 23*z2) + 
																						x2*(-8*y2*y2 + 90*y2*z2 + 
																						8*z2*z2)))*invR13;
																						fx+=srcDenarr[dims*s+31]*0.0625*cl53*(-6*x*y*(7*x2*x2 - 5*y2*y2 + 44*y2*z2 + 
																						16*z2*z2 + 2*x2*(y2 - 38*z2)))*invR13;
																						fx+=srcDenarr[dims*s+32]*0.375*cl54*(x*z*(7*x2*x2 + 23*y2*y2 + 12*y2*z2 - 
																						2*x2*(29*y2 + 2*z2)))*invR13;
																						fx+=srcDenarr[dims*s+33]*1.5*cl54*(y*z*(8*x2*x2 + y2*(y2 + z2) - 
																						x2*(13*y2 + 3*z2)))*invR13;
																						fx+=srcDenarr[dims*s+34]*0.375*cl55*(6*x2*x2*x2 - 5*y2*y2*(y2 + z2) - 
																						5*x2*x2*(17*y2 + z2) + 
																						10*x2*(8*y2*y2 + 3*y2*z2))*invR13;
																						fx+=srcDenarr[dims*s+35]*0.375*cl55*(x*y*(35*x2*x2 + 31*y2*y2 + 20*y2*z2 - 
																						10*x2*(11*y2 + 2*z2)))*invR13;

																						fx+=srcDenarr[dims*s+36]*(x*(-2.1875*x2*x2*x2 - 2.1875*y2*y2*y2 + 
																						52.5*y2*y2*z2 - 105.*y2*z2*z2 + 
																						28.*z2*z2*z2 + x2*x2*
																						(-6.5625*y2 + 52.5*z2) + 
																						x2*(-6.5625*y2*y2 + 105.*y2*z2 - 
																						105.*z2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+37]*0.125*cl61*(z*(40*x2*x2*x2 - 5*y2*y2*y2 + 15*y2*y2*z2 + 
																						12*y2*z2*z2 - 8*z2*z2*z2 + 
																						75*x2*x2*(y2 - 3*z2) + 
																						6*x2*(5*y2*y2 - 35*y2*z2 + 
																						26*z2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+38]*0.125*cl61*(3*x*y*z*(15*x2*x2 + 15*y2*y2 - 80*y2*z2 + 
																						48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;
																						fx+=srcDenarr[dims*s+39]*0.0625*cl62*(x*(7*x2*x2*x2 - 11*y2*y2*y2 + 210*y2*y2*z2 - 
																						240*y2*z2*z2 - 32*z2*z2*z2 + 
																						3*x2*x2*(y2 - 50*z2) - 
																						15*x2*(y2*y2 - 4*y2*z2 - 
																						16*z2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+40]*0.125*cl62*(y*(8*x2*x2*x2 - y2*y2*y2 + 15*y2*y2*z2 - 
																						16*z2*z2*z2 + 15*x2*x2*(y2 - 11*z2) + 
																						6*x2*(y2*y2 - 25*y2*z2 + 
																						40*z2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+41]*0.0625*cl63*(-2*z*(24*x2*x2*x2 - 5*x2*x2*(15*y2 + 19*z2) + 
																						3*y2*(3*y2*y2 - 5*y2*z2 - 
																						8*z2*z2) + x2*
																						(-90*y2*y2 + 330*y2*z2 + 24*z2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+42]*0.0625*cl63*(-2*x*y*z*(81*x2*x2 - 51*y2*y2 + 140*y2*z2 + 
																						48*z2*z2 + 30*x2*(y2 - 10*z2)))*invR15;
																						fx+=srcDenarr[dims*s+43]*0.09375*cl64*(-2*x*(7*x2*x2*x2 + 23*y2*y2*y2 - 240*y2*y2*z2 - 
																						120*y2*z2*z2 - 
																						3*x2*x2*(17*y2 + 32*z2) + 
																						x2*(-35*y2*y2 + 720*y2*z2 + 
																						40*z2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+44]*0.375*cl64*(-2*y*(8*x2*x2*x2 + y2*y2*y2 - 9*y2*y2*z2 - 
																						10*y2*z2*z2 - 
																						5*x2*x2*(y2 + 21*z2) - 
																						6*x2*(2*y2*y2 - 25*y2*z2 - 
																						5*z2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+45]*0.375*cl65*(z*(8*x2*x2*x2 - 5*y2*y2*(y2 + z2) + 
																						30*x2*y2*(3*y2 + z2) - 
																						5*x2*x2*(21*y2 + z2)))*invR15;
																						fx+=srcDenarr[dims*s+46]*0.375*cl65*(x*y*z*(45*x2*x2 + 33*y2*y2 + 20*y2*z2 - 
																						10*x2*(13*y2 + 2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+47]*0.0625*cl66*(x*(7*x2*x2*x2 - 43*y2*y2*y2 - 30*y2*y2*z2 - 
																						3*x2*x2*(47*y2 + 2*z2) + 
																						15*x2*(15*y2*y2 + 4*y2*z2)))*invR15;
																						fx+=srcDenarr[dims*s+48]*0.0625*cl66*(2*y*(24*x2*x2*x2 - 3*y2*y2*(y2 + z2) - 
																						5*x2*x2*(23*y2 + 3*z2) + 
																						6*x2*(11*y2*y2 + 5*y2*z2)))*invR15;*/

																						double fy =srcDenarr[dims*s]*dY_reg*invR3;
																						fy+=srcDenarr[dims*s+1]*3.0*dY_reg*dZ_reg*invR5;
																						fy+=srcDenarr[dims*s+2]*3.0*dX_reg*dY_reg*invR5;
																						fy+=srcDenarr[dims*s+3]*(3.0*y2-1.0/invR2)*invR5;

																						/*	fy+=srcDenarr[dims*s+4]*1.5*dY_reg*(5.*z2-1.0/invR2)*invR7;
																						fy+=srcDenarr[dims*s+5]*cl2*5.0*pm1c*pm1s*dZ_reg*invR7;
																						fy+=srcDenarr[dims*s+6]*cl2*dZ_reg*(5.0*y2-1.0/invR2)*invR7;
																						fy+=srcDenarr[dims*s+7]*0.5*cl2*dY_reg*(5.0*pm2c-2.0/invR2)*invR7;
																						fy+=srcDenarr[dims*s+8]*cl2*pm1c*(5.0*y2-1.0/invR2)*invR7;

																						fy+=srcDenarr[dims*s+9]*dY_reg*dZ_reg*(10.0*z2-7.5*x2-7.5*y2)*invR9;
																						fy+=srcDenarr[dims*s+10]*0.25*cl31*(-5.0*pm2s*(x2 + y2 - 6*z2))*invR9;
																						fy+=srcDenarr[dims*s+11]*0.25*cl31*(x2*x2 - 4*y2*y2 + 27*y2*z2 - 4*z2*z2 - 
																						3*x2*(y2 + z2))*invR9;
																						fy+=srcDenarr[dims*s+12]*0.5*cl32*(dY_reg*dZ_reg*(9*x2 - 5*y2 + 2*z2))*invR9;
																						fy+=srcDenarr[dims*s+13]*cl32*(-(x*z*(x2 - 6*y2 + z2)))*invR9;
																						fy+=srcDenarr[dims*s+14]*0.25*cl33*(pm2s*(13*x2 - 15*y2 + 6*z2))*invR9;
																						fy+=srcDenarr[dims*s+15]*0.25*cl33*(-3*x2*x2 - 4*y2*y2 + 3*y2*z2 + 
																						3*x2*(7*y2 - z2))*invR9;

																						fy+=srcDenarr[dims*s+16]*(dY_reg*(1.875*x2*x2 + 1.875*y2*y2 - 22.5*y2*z2 + 
																						15.*z2*z2 + x2*(3.75*y2 - 22.5*z2)))*invR11;
																						fy+=srcDenarr[dims*s+17]*0.25*cl41*(-21*x*y*z*(x2 + y2 - 2*z2))*invR11;
																						fy+=srcDenarr[dims*s+18]*0.25*cl41*(-1.0*z*(-3*x2*x2 + 18*y2*y2 - 41*y2*z2 + 
																						4*z2*z2 + x2*(15*y2 + z2)))*invR11;
																						fy+=srcDenarr[dims*s+19]*0.25*cl42*(y*(-9*x2*x2 + 5*y2*y2 - 46*y2*z2 + 
																						12*z2*z2 + x2*(-4*y2 + 66*z2)))*invR11;
																						fy+=srcDenarr[dims*s+20]*0.5*cl42*(x*(x2*x2 - 6*y2*y2 + 51*y2*z2 - 
																						6*z2*z2 - 5*x2*(y2 + z2)))*invR11;
																						fy+=srcDenarr[dims*s+21]*0.25*cl43*(3*x*y*(5*x2 - 7*y2 + 2*z2))*invR11;
																						fy+=srcDenarr[dims*s+22]*0.25*cl43*(-3*(x2*x2 + 2*y2*y2 - y2*z2 + 
																						x2*(-9*y2 + z2)))*invR11;
																						fy+=srcDenarr[dims*s+23]*0.125*cl44*(y*(21*x2*x2 + 5*y2*y2 - 4*y2*z2 + 
																						x2*(-46*y2 + 12*z2)))*invR11;
																						fy+=srcDenarr[dims*s+24]*0.5*cl44*(-1.0*x*(x2*x2 + 6*y2*y2 - 3*y2*z2 + 
																						x2*(-11*y2 + z2)))*invR11;

																						fy+=srcDenarr[dims*s+25]*(y*z*(13.125*x2*x2 + 13.125*y2*y2 - 
																						52.5*y2*z2 + 21.*z2*z2 + 
																						x2*(26.25*y2 - 52.5*z2)))*invR13;
																						fy+=srcDenarr[dims*s+26]*0.125*cl51*(7*x*y*(x2*x2 + y2*y2 - 16*y2*z2 + 
																						16*z2*z2 + 2*x2*(y2 - 8*z2)))*invR13;
																						fy+=srcDenarr[dims*s+27]*0.125*cl51*(-x2*x2*x2 + 6*y2*y2*y2 - 101*y2*y2*z2 + 
																						116*y2*z2*z2 - 8*z2*z2*z2 + 
																						x2*x2*(4*y2 + 11*z2) + 
																						x2*(11*y2*y2 - 90*y2*z2 + 4*z2*z2))*invR13;
																						fy+=srcDenarr[dims*s+28]*0.125*cl52*(y*z*(-11*x2*x2 + 7*y2*y2 - 22*y2*z2 + 
																						4*z2*z2 + x2*(-4*y2 + 26*z2)))*invR13;
																						fy+=srcDenarr[dims*s+29]*0.25*cl52*(2*x*z*(x2*x2 - 8*y2*y2 + 23*y2*z2 - 
																						2*z2*z2 - x2*(7*y2 + z2)))*invR13;
																						fy+=srcDenarr[dims*s+30]*0.0625*cl53*(-6*x*y*(5*x2*x2 - 7*y2*y2 + 76*y2*z2 - 
																						16*z2*z2 - 2*x2*(y2 + 22*z2)))*invR13;
																						fy+=srcDenarr[dims*s+31]*0.0625*cl53*(6*(x2*x2*x2 + 2*y2*y2*y2 - 23*y2*y2*z2 + 
																						8*y2*z2*z2 - 
																						x2*x2*(8*y2 + 7*z2) + 
																						x2*(-7*y2*y2 + 90*y2*z2 - 
																						8*z2*z2)))*invR13;
																						fy+=srcDenarr[dims*s+32]*0.375*cl54*(y*z*(23*x2*x2 + 7*y2*y2 - 4*y2*z2 + 
																						x2*(-58*y2 + 12*z2)))*invR13;
																						fy+=srcDenarr[dims*s+33]*1.5*cl54*(-1.0*x*z*(x2*x2 + 8*y2*y2 - 3*y2*z2 + 
																						x2*(-13*y2 + z2)))*invR13;
																						fy+=srcDenarr[dims*s+34]*0.375*cl55*(x*y*(31*x2*x2 + 5*y2*(7*y2 - 4*z2) + 
																						x2*(-110*y2 + 20*z2)))*invR13;
																						fy+=srcDenarr[dims*s+35]*0.375*cl55*(-5*x2*x2*x2 + 6*y2*y2*y2 - 5*y2*y2*z2 + 
																						x2*x2*(80*y2 - 5*z2) + 
																						x2*(-85*y2*y2 + 30*y2*z2))*invR13;

																						fy+=srcDenarr[dims*s+36]*(y*(-2.1875*x2*x2*x2 - 2.1875*y2*y2*y2 + 
																						52.5*y2*y2*z2 - 105.*y2*z2*z2 + 
																						28.*z2*z2*z2 + x2*x2*
																						(-6.5625*y2 + 52.5*z2) + 
																						x2*(-6.5625*y2*y2 + 105.*y2*z2 - 
																						105.*z2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+37]*0.125*cl61*(3*x*y*z*(15*x2*x2 + 15*y2*y2 - 80*y2*z2 + 
																						48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;
																						fy+=srcDenarr[dims*s+38]*0.125*cl61*(z*(-5*x2*x2*x2 + 40*y2*y2*y2 - 225*y2*y2*z2 + 
																						156*y2*z2*z2 - 8*z2*z2*z2 + 
																						15*x2*x2*(2*y2 + z2) + 
																						3*x2*(25*y2*y2 - 70*y2*z2 + 
																						4*z2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+39]*0.0625*cl62*(y*(11*x2*x2*x2 - 7*y2*y2*y2 + 150*y2*y2*z2 - 
																						240*y2*z2*z2 + 32*z2*z2*z2 + 
																						15*x2*x2*(y2 - 14*z2) - 
																						3*x2*(y2*y2 + 20*y2*z2 - 
																						80*z2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+40]*0.125*cl62*(x*(-x2*x2*x2 + 8*y2*y2*y2 - 165*y2*y2*z2 + 
																						240*y2*z2*z2 - 16*z2*z2*z2 + 
																						3*x2*x2*(2*y2 + 5*z2) + 
																						15*x2*(y2*y2 - 10*y2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+41]*0.0625*cl63*(2*x*y*z*(-51*x2*x2 + 81*y2*y2 - 300*y2*z2 + 
																						48*z2*z2 + 10*x2*(3*y2 + 14*z2)))*invR15;
																						fy+=srcDenarr[dims*s+42]*0.0625*cl63*(2*z*(9*x2*x2*x2 + 24*y2*y2*y2 - 95*y2*y2*z2 + 
																						24*y2*z2*z2 - 
																						15*x2*x2*(6*y2 + z2) - 
																						3*x2*(25*y2*y2 - 110*y2*z2 + 
																						8*z2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+43]*0.09375*cl64*(-2*y*(23*x2*x2*x2 + 7*y2*y2*y2 - 96*y2*y2*z2 + 
																						40*y2*z2*z2 - 
																						5*x2*x2*(7*y2 + 48*z2) - 
																						3*x2*(17*y2*y2 - 240*y2*z2 + 
																						40*z2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+44]*0.375*cl64*(2*x*(x2*x2*x2 + 8*y2*y2*y2 - 105*y2*y2*z2 + 
																						30*y2*z2*z2 - 
																						3*x2*x2*(4*y2 + 3*z2) - 
																						5*x2*(y2*y2 - 30*y2*z2 + 
																						2*z2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+45]*0.375*cl65*(x*y*z*(33*x2*x2 + 5*y2*(9*y2 - 4*z2) + 
																						x2*(-130*y2 + 20*z2)))*invR15;
																						fy+=srcDenarr[dims*s+46]*0.375*cl65*(z*(-5*x2*x2*x2 + 8*y2*y2*y2 - 5*y2*y2*z2 + 
																						x2*x2*(90*y2 - 5*z2) - 
																						15*x2*(7*y2*y2 - 2*y2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+47]*0.0625*cl66*(y*(43*x2*x2*x2 - 7*y2*y2*y2 + 6*y2*y2*z2 + 
																						x2*x2*(-225*y2 + 30*z2) + 
																						3*x2*(47*y2*y2 - 20*y2*z2)))*invR15;
																						fy+=srcDenarr[dims*s+48]*0.0625*cl66*(-2*x*(3*x2*x2*x2 - 24*y2*y2*y2 + 15*y2*y2*z2 + 
																						x2*x2*(-66*y2 + 3*z2) + 
																						5*x2*(23*y2*y2 - 6*y2*z2)))*invR15;*/

																						double fz =srcDenarr[dims*s]*dZ_reg*invR3;
																						fz+=srcDenarr[dims*s+1]*(3.0*z2-1.0/invR2)*invR5;
																						fz+=srcDenarr[dims*s+2]*3.0*dX_reg*dZ_reg*invR5;
																						fz+=srcDenarr[dims*s+3]*3.0*dY_reg*dZ_reg*invR5;

																						/*	fz+=srcDenarr[dims*s+4]*1.5*dZ_reg*(5.*z2-3.0/invR2)*invR7;
																						fz+=srcDenarr[dims*s+5]*cl2*dX_reg*(5.0*z2-1.0/invR2)*invR7;
																						fz+=srcDenarr[dims*s+6]*cl2*dY_reg*(5.0*z2-1.0/invR2)*invR7;
																						fz+=srcDenarr[dims*s+7]*0.5*cl2*5.0*dZ_reg*pm2c*invR7;
																						fz+=srcDenarr[dims*s+8]*cl2*5.0*pm1c*pm1s*dZ_reg*invR7;

																						fz+=srcDenarr[dims*s+9]*(1.5*x2*x2+1.5*y2*y2-12.*y2*z2+4.*z2*z2+z2*(3.*y2-12.*z2))*invR9;
																						fz+=srcDenarr[dims*s+10]*0.25*cl31*(-5.0*dX_reg*dZ_reg*(3*x2 + 3*y2 - 4*x2))*invR9;
																						fz+=srcDenarr[dims*s+11]*0.25*cl31*(-5*dY_reg*dZ_reg*(3*x2 + 3*y2 - 4*z2))*invR9;
																						fz+=srcDenarr[dims*s+12]*0.5*cl32*((-x2 + y2)*(x2 + y2 - 6*z2))*invR9;
																						fz+=srcDenarr[dims*s+13]*cl32*(-(x*y*(x2 + y2 - 6*z2)))*invR9;
																						fz+=srcDenarr[dims*s+14]*0.25*cl33*(7*dX_reg*dZ_reg*(x2 - 3*y2))*invR9;
																						fz+=srcDenarr[dims*s+15]*0.25*cl33*(7*dY_reg*dZ_reg*(3*x2 - y2))*invR9;

																						fz+=srcDenarr[dims*s+16]*(dZ_reg*(9.375*x2*x2 + 9.375*y2*y2 - 25.*y2*z2 + 
																						5.*z2*z2 + x2*(18.75*y2 - 25.*z2)))*invR11;
																						fz+=srcDenarr[dims*s+17]*0.25*cl41*(3*x*(x2*x2 + y2*y2 - 12*y2*z2 + 
																						8*z2*z2 + 2*x2*(y2 - 6*z2)))*invR11;
																						fz+=srcDenarr[dims*s+18]*0.25*cl41*(3*y*(x2*x2 + y2*y2 - 12*y2*z2 + 
																						8*z2*z2 + 2*x2*(y2 - 6*z2)))*invR11;
																						fz+=srcDenarr[dims*s+19]*0.25*cl42*(-21*(x2 - y2)*z*(x2 + y2 - 2*z2))*invR11;
																						fz+=srcDenarr[dims*s+20]*0.5*cl42*(-21*x*y*z*(x2 + y2 - 2*z2))*invR11;
																						fz+=srcDenarr[dims*s+21]*0.25*cl43*dZ_reg*(9*(x2 - 3*y2)*z*x)*invR11;
																						fz+=srcDenarr[dims*s+22]*0.25*cl43*dZ_reg*(9*(3*x2 - y2)*z*y)*invR11;
																						fz+=srcDenarr[dims*s+23]*0.125*cl44*(9*(x2*x2 - 6*x2*y2 + y2*y2)*z)*invR11;
																						fz+=srcDenarr[dims*s+24]*0.5*cl44*(9*x*y*(x2 - y2)*z)*invR11;

																						fz+=srcDenarr[dims*s+25]*(-1.875*x2*x2*x2 - 1.875*y2*y2*y2 + 33.75*y2*y2*z2 - 
																						45.*y2*z2*z2 + 6.*z2*z2*z2 + 
																						x2*x2*(-5.625*y2 + 33.75*z2) + 
																						x2*(-5.625*y2*y2 + 67.5*y2*z2 - 
																						45.*z2*z2))*invR13;
																						fz+=srcDenarr[dims*s+26]*0.125*cl51*(7*x*z*(5*x2*x2 + 5*y2*y2 - 20*y2*z2 + 
																						8*z2*z2 + 10*x2*(y2 - 2*z2)))*invR13;
																						fz+=srcDenarr[dims*s+27]*0.125*cl51*(7*y*z*(5*x2*x2 + 5*y2*y2 - 20*y2*z2 + 
																						8*z2*z2 + 10*x2*(y2 - 2*z2)))*invR13;
																						fz+=srcDenarr[dims*s+28]*0.125*cl52*(2*(x2 - y2)*(x2*x2 + y2*y2 - 
																						16*y2*z2 + 16*z2*z2 + 
																						2*x2*(y2 - 8*z2)))*invR13;
																						fz+=srcDenarr[dims*s+29]*0.25*cl52*(2*x*y*(x2*x2 + y2*y2 - 16*y2*z2 + 
																						16*z2*z2 + 2*x2*(y2 - 8*z2)))*invR13;
																						fz+=srcDenarr[dims*s+30]*0.0625*cl53*(-18*x*(x2 - 3*y2)*z*
																						(3*x2 + 3*y2 - 8*z2))*invR13;
																						fz+=srcDenarr[dims*s+31]*0.0625*cl53*(18*y*(-3*x2 + y2)*z*
																						(3*x2 + 3*y2 - 8*z2))*invR13;
																						fz+=srcDenarr[dims*s+32]*0.375*cl54*(-1.0*(x2*x2 - 6*x2*y2 + y2*y2)*
																						(x2 + y2 - 10*z2))*invR13;
																						fz+=srcDenarr[dims*s+33]*1.5*cl54*(-1.0*x*y*(x2 - y2)*(x2 + y2 - 10*z2))*invR13;
																						fz+=srcDenarr[dims*s+34]*0.375*cl55*(11*x*(x2*x2 - 10*x2*y2 + 5*y2*y2)*z)*invR13;
																						fz+=srcDenarr[dims*s+35]*0.375*cl55*(-5*x2*x2*x2 + 6*y2*y2*y2 - 5*y2*y2*z2 + 
																						x2*x2*(80*y2 - 5*z2) + 
																						x2*(-85*y2*y2 + 30*y2*z2))*invR13;

																						fz+=srcDenarr[dims*s+36]*(z*(-15.3125*x2*x2*x2 - 15.3125*y2*y2*y2 + 
																						91.875*y2*y2*z2 - 73.5*y2*z2*z2 + 
																						7.*z2*z2*z2 + x2*x2*
																						(-45.9375*y2 + 91.875*z2) + 
																						x2*(-45.9375*y2*y2 + 183.75*y2*z2 - 
																						73.5*z2*z2)))*invR15;
																						fz+=srcDenarr[dims*s+37]*0.125*cl61*(3*x*y*z*(15*x2*x2 + 15*y2*y2 - 80*y2*z2 + 
																						48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;
																						fz+=srcDenarr[dims*s+38]*0.125*cl61*(-1.0*y*(5*x2*x2*x2 + 5*y2*y2*y2 - 120*y2*y2*z2 + 
																						240*y2*z2*z2 - 64*z2*z2*z2 + 
																						15*x2*x2*(y2 - 8*z2) + 
																						15*x2*(y2*y2 - 16*y2*z2 + 
																						16*z2*z2)))*invR15;
																						fz+=srcDenarr[dims*s+39]*0.0625*cl62*(3*(x2 - y2)*z*
																						(15*x2*x2 + 15*y2*y2 - 80*y2*z2 + 
																						48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;
																						fz+=srcDenarr[dims*s+40]*0.125*cl62*(3*x*y*z*(15*x2*x2 + 15*y2*y2 - 80*y2*z2 + 
																						48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;
																						fz+=srcDenarr[dims*s+41]*0.0625*cl63*(2*x*(x2 - 3*y2)*
																						(3*x2*x2 + 3*y2*y2 - 60*y2*z2 + 
																						80*z2*z2 + 6*x2*(y2 - 10*z2)))*invR15;
																						fz+=srcDenarr[dims*s+42]*0.0625*cl63*(-2*y*(-3*x2 + y2)*
																						(3*x2*x2 + 3*y2*y2 - 60*y2*z2 + 
																						80*z2*z2 + 6*x2*(y2 - 10*z2)))*invR15;
																						fz+=srcDenarr[dims*s+43]*0.09375*cl64*(-22*(x2*x2 - 6*x2*y2 + y2*y2)*z*
																						(3*x2 + 3*y2 - 10*z2))*invR15;
																						fz+=srcDenarr[dims*s+44]*0.375*cl64*(-22*x*y*(x2 - y2)*z*
																						(3*x2 + 3*y2 - 10*z2))*invR15;
																						fz+=srcDenarr[dims*s+45]*0.375*cl65*(-1.0*x*(x2*x2 - 10*x2*y2 + 5*y2*y2)*
																						(x2 + y2 - 12*z2))*invR15;
																						fz+=srcDenarr[dims*s+46]*0.375*cl65*(-1.0*y*(5*x2*x2 - 10*x2*y2 + y2*y2)*
																						(x2 + y2 - 12*z2))*invR15;
																						fz+=srcDenarr[dims*s+47]*0.0625*cl66*(13*(x2*x2*x2 - 15*x2*x2*y2 + 15*x2*y2*y2 - 
																						y2*y2*y2)*z)*invR15;
																						fz+=srcDenarr[dims*s+48]*0.0625*cl66*(26*x*y*(3*x2*x2 - 10*x2*y2 + 3*y2*y2)*z)*invR15;*/

																						//p[0] += phi;
																						
																							pot_m[3*imind[t]+0] += fx*imq[t]; 
																							pot_m[3*imind[t]+1] += fy*imq[t]; 
																							pot_m[3*imind[t]+2] += fz*imq[t]; 

																							
																							pot_m[3*s+0] -= fx*imq[t]; 
																							pot_m[3*s+1] -= fy*imq[t]; 
																							pot_m[3*s+2] -= fz*imq[t]; 																	
																						
																						
																					}
																				}																																								
																			}

																						/*ion-ion interaction*/

																						for(i=0;i<N_ion-1;i++)
																							for(ii=i+1;ii<N_ion;ii++)
																							{
																								dx=-x[i]+x[ii];
																								dy=-y[i]+y[ii];
																								dz=-z[i]+z[ii];
																								dr=sqrt(dx*dx+dy*dy+dz*dz);

																								temfx=q[ii]*q[i]*dx/dr/dr/dr/epsi_s;
																								temfy=q[ii]*q[i]*dy/dr/dr/dr/epsi_s;
																								temfz=q[ii]*q[i]*dz/dr/dr/dr/epsi_s;

																								pot_m[3*ii+0]+=temfx;
																								pot_m[3*ii+1]+=temfy;
																								pot_m[3*ii+2]+=temfz;

																								pot_m[3*i+0]-=temfx;
																								pot_m[3*i+1]-=temfy;
																								pot_m[3*i+2]-=temfz;																																						
																							}


																							if(iter_indicator==0)
																							{
																								for ( i=0 ; i<N ; i++ )
																								{
																									acc_x[i] += pot_m[3*i+0]/ mass[i];
																									acc_y[i] += pot_m[3*i+1]/ mass[i];
																									acc_z[i] += pot_m[3*i+2]/ mass[i];
																									//	fprintf(force,"%d %f %f %f\n",i+1,pot_m[3*i+0],pot_m[3*i+1],pot_m[3*i+2]);
																								}
																							}
																							else
																							{
																								for ( i=0 ; i<N ; i++ )
																								{
																									acc_x[i] += pot_m[3*i+0]/ mass[i];
																									acc_y[i] += pot_m[3*i+1]/ mass[i];
																									acc_z[i] += pot_m[3*i+2]/ mass[i];
																									//fprintf(force,"%d %f %f %f\n",i,acc_x[i],acc_y[i],acc_z[i]);
																									//	fprintf(force,"%d %f %f %f\n",i,acc_x[i],acc_y[i],acc_z[i]);
																								}
																							}



																							/*rescale back bknm for next iteration*/
																							for(i=0;i<M;i++)
																								for(j=0;j<2*p+1;j++)
																									for(k=0;k<p+1;k++)
																									{
																										Bknm[i][j][k].real=Bknm[i][j][k].real/sqrtk[k];
																										Bknm[i][j][k].imag=Bknm[i][j][k].imag/sqrtk[k];
																									} 

																									/*----------------------Compute the force----------------------------*/
														}

														/*--------------compute the potential energy------------------------*/	
														if(energy_compute==1)
														{
															for(i=0;i<M;i++)
																for(ii=0;ii<N_ion;ii++)
																{ 
																	loc[0]=x[ii]-ox[i];
																	loc[1]=y[ii]-oy[i];
																	loc[2]=z[ii]-oz[i];
																	r=sqrt(loc[0]*loc[0]+loc[1]*loc[1]+loc[2]*loc[2]);

																	powr=1.0;
																	for(k=0;k<p+1;k++)
																	{
																		powr=powr*r;
																		for(j=0;j<2*p+1;j++)

																		{

																			Bknm[i][j][k].real=Bknm[i][j][k].real/powr;
																			Bknm[i][j][k].imag=Bknm[i][j][k].imag/powr;

																		} //rescale.
																	}
																	//getchar();
																	ssheval_(Bknm[i],&p,loc,&fval);

																	energy=energy+0.5*q[ii]*fval.real;

																	powr=1.0;
																	for(k=0;k<p+1;k++)
																	{
																		powr=powr*r;
																		for(j=0;j<2*p+1;j++)

																		{
																			Bknm[i][j][k].real=Bknm[i][j][k].real*powr;
																			Bknm[i][j][k].imag=Bknm[i][j][k].imag*powr;
																		} //rescale back.
																	}

																}

																//printf("total self energy  (harmonic-ion)= %.15f\n",energy);

																for(i=0;i<M;i++)
																	for(ii=0;ii<M;ii++)
																	{

																		if(i!=ii)
																		{
																			loc[0]=ox[ii]-ox[i];
																			loc[1]=oy[ii]-oy[i];
																			loc[2]=oz[ii]-oz[i];
																			r=sqrt(loc[0]*loc[0]+loc[1]*loc[1]+loc[2]*loc[2]);

																			powr=1.0;
																			for(k=0;k<p+1;k++)
																			{
																				powr=powr*r;
																				for(j=0;j<2*p+1;j++)
																				{
																					Bknm[i][j][k].real=Bknm[i][j][k].real/powr;
																					Bknm[i][j][k].imag=Bknm[i][j][k].imag/powr;
																				} //rescale.
																			}

																			ssheval_(Bknm[i],&p,loc,&fval);

																			energy=energy+0.5*osigma[ii]*fval.real;
																			powr=1.0;
																			for(k=0;k<p+1;k++)
																			{
																				powr=powr*r;
																				for(j=0;j<2*p+1;j++)
																				{
																					Bknm[i][j][k].real=Bknm[i][j][k].real*powr;
																					Bknm[i][j][k].imag=Bknm[i][j][k].imag*powr;
																				} //rescale back.
																			}
																		}
																		/*	else
																		{
																		for(j=0;j<ntheta;j++)
																		for(k=0;k<nphi;k++)
																		{
																		loc[0]=-ox[i]+Mnodes[ii][j][k][0];
																		loc[1]=-oy[i]+Mnodes[ii][j][k][1];
																		loc[2]=-oz[i]+Mnodes[ii][j][k][2];


																		ssheval_(Bknm[i],&p,loc,&fval);

																		energy=energy+0.5*weights[j][k]*osigma[ii]*fval.real;





																		}	
																		}*/

																	}

																	//		printf("total self energy  (harmonic-sph)= %.15f\n",energy);



																	double tttt[1000];
																	if(N_ion>2000)
																	{
																		tstart_direct =0;//omp_get_wtime();

																		for(i=0;i<100*im;i++)																			
																		{
																			//#pragma omp parallel for reduction(+:energy) //default(shared) private(dx,dy,dz,dr,tttt) 
																			for(j=0;j<100;j++)
																			{

																				dx=imx[i]-x[j];
																				dy=imy[i]-y[j];
																				dz=imz[i]-z[j];
																				dr=sqrt(dx*dx+dy*dy+dz*dz);
																				if(i>=j*im*M&&i<(j+1)*im*M)
																				{
																					tttt[j]=0.5*q[j]*imq[i]/dr/epsi_s;
																					energy=energy+tttt[j];
																				}
																				else
																				{
																					tttt[j]=q[j]*imq[i]/dr/epsi_s;
																					energy=energy+tttt[j];
																				}
																				//printf("imq[i]=%d %.13f %.13f %.13f %.13f\n",i,imx[i],imy[i],imz[i],imq[i]);
																			}
																		}

																		tfinish_direct = 0;//omp_get_wtime();
																	}

																	else if(N_ion<=2000)
																	{
																		for(i=0;i<imcount;i++)
																			for(ii=0;ii<N_ion;ii++)
																			{
																				dx=imx[i]-x[ii];
																				dy=imy[i]-y[ii];
																				dz=imz[i]-z[ii];
																				dr=sqrt(dx*dx+dy*dy+dz*dz);
																				//if(i>=ii*M*im&&i<(ii*M*im+M*im))
																				energy=energy+0.5*q[ii]*imq[i]/dr;
																				//else
																				//energy=energy+0.5*q[ii]*imq[i]/dr/epsi_s;

																			}

																	}
																	//	printf("total energy (image-source)= %.15f\n",energy*epsi_s);

																	for(i=0;i<imcount;i++)
																		for(ii=N_ion;ii<N;ii++)
																		{
																			dx=imx[i]-x[ii];
																			dy=imy[i]-y[ii];
																			dz=imz[i]-z[ii];
																			dr=sqrt(dx*dx+dy*dy+dz*dz);
																			if(dr>1.0)
																				energy=energy+0.5*q[ii]*imq[i]/dr;																																						
																		}

																		printf("total energy (image-sph)= %.15f\n",energy*epsi_s);

																		for(i=0;i<N_ion;i++)
																			for(ii=i+1;ii<N_ion;ii++)
																			{
																				loc[0]=x[ii]-x[i];
																				loc[1]=y[ii]-y[i];
																				loc[2]=z[ii]-z[i];
																				r=sqrt(loc[0]*loc[0]+loc[1]*loc[1]+loc[2]*loc[2]);
																				energy=energy+q[i]*q[ii]/r;
																			}

																			printf("total energy (ion-ion)= %.15f\n",energy*epsi_s);

																			for(i=0;i<N_ion;i++)
																				for(ii=N_ion;ii<N;ii++)
																				{
																					loc[0]=x[ii]-x[i];
																					loc[1]=y[ii]-y[i];
																					loc[2]=z[ii]-z[i];
																					r=sqrt(loc[0]*loc[0]+loc[1]*loc[1]+loc[2]*loc[2]);
																					energy=energy+q[i]*q[ii]/r;

																					//printf("%f %f %f\n",r, q[i],q[ii]);
																				}
																				printf("total energy(ion-sph)= %.15f\n",energy*epsi_s);

																				for(i=N_ion;i<N;i++)
																					for(ii=i+1;ii<N;ii++)
																					{
																						loc[0]=x[ii]-x[i];
																						loc[1]=y[ii]-y[i];
																						loc[2]=z[ii]-z[i];
																						r=sqrt(loc[0]*loc[0]+loc[1]*loc[1]+loc[2]*loc[2]);
																						energy=energy+q[i]*q[ii]/r;
																					}
																					printf("total energy(sph-sph)= %.15f\n",energy*epsi_s);

														}


														//iprint=1;				
														if ( iprint == 1 )
															output_force();
														//	getchar();
}

void matvec(int n, double *vecx, double *b)
{
	int i,j,k,l;
	int ii,jj,kk,ll,iii;
	double sc1=1.0, sc2=1.0;
	int ier=0;
	double x0y0z0[3];
	double xnynzn[3];
	double d,dx,dy,dz,dr,rk,rim;

	//generate_imagecoef(n, vecx);

	for(i=0;i<M;i++)
		for(ii=0;ii<M;ii++)
			for(j=0;j<2*p+1;j++)
				for(k=0;k<p+1;k++)
				{
					if(i!=ii)
					{

						mpole[i][ii][j][k].real=vecx[i*(2*p+1)*(p+1)+j*(p+1)+k];
						mpole[i][ii][j][k].imag=vecx[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k];
						local[i][ii][j][k].real=0.0;
						local[i][ii][j][k].imag=0.0;
					}
				}

				for(i=0;i<M;i++)
					for(ii=0;ii<M;ii++)
					{
						x0y0z0[0]=ox[i];
						x0y0z0[1]=oy[i];
						x0y0z0[2]=oz[i];

						xnynzn[0]=ox[ii];
						xnynzn[1]=oy[ii];
						xnynzn[2]=oz[ii];
						if(i!=ii)
						{
							dx=ox[i]-ox[ii];
							dy=oy[i]-oy[ii];
							dz=oz[i]-oz[ii];
							dr=dx*dx+dy*dy+dz*dz;

							if(dr<m2ltr)
							{
								for(j=0;j<2*p+1;j++)
									for(k=0;k<p+1;k++)
									{


										mpole[i][ii][j][k].real=mpole[i][ii][j][k].real*sqrtk[k];
										mpole[i][ii][j][k].imag=mpole[i][ii][j][k].imag*sqrtk[k];

									}
									l3dmplocquadu_(&sc1,x0y0z0,mpole[i][ii],&p,&sc2,xnynzn,local[i][ii],&p,&ier);

							}
						}
					}

					for(i=0;i<M;i++)
						for(ii=0;ii<M;ii++)
						{
							if(i!=ii)
							{
								for(j=0;j<2*p+1;j++)
									for(k=0;k<p+1;k++)
									{
										if(abs(j-p)<=k)
										{
											local[i][ii][j][k].real=local[i][ii][j][k].real/sqrtk[k];
											local[i][ii][j][k].imag=local[i][ii][j][k].imag/sqrtk[k];
											//printf("lpole=%d %d %d %d %f %f\n",i, ii,k,j-p, local[i][ii][j][k].real,local[i][ii][j][k].imag);
										}
									}
							}
						}



						for(i=0;i<M;i++)
							for(j=0;j<2*p+1;j++)
								for(k=0;k<p+1;k++)
								{
									b[i*(2*p+1)*(p+1)+j*(p+1)+k]=0.0;
									b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]=0.0;
								}



								for(i=0;i<M;i++)
									for(ii=0;ii<M;ii++)
										for(j=0;j<2*p+1;j++)
											for(k=0;k<p+1;k++)						
											{
												if(abs(j-p)<=k)
												{
													if(ii==i)
													{
														if(k==0)
														{

															b[i*(2*p+1)*(p+1)+j*(p+1)+k]=b[i*(2*p+1)*(p+1)+j*(p+1)+k]+vecx[i*(2*p+1)*(p+1)+j*(p+1)+k]*((k+1)*epsi_s+k*epsi_i[i]);
															b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]=b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]+vecx[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]*((k+1)*epsi_s+k*epsi_i[i]);
														}
														else
														{
															b[i*(2*p+1)*(p+1)+j*(p+1)+k]=b[i*(2*p+1)*(p+1)+j*(p+1)+k]+vecx[i*(2*p+1)*(p+1)+j*(p+1)+k]*((k+1)*epsi_s/k+epsi_i[i]);
															b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]=b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]+vecx[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]*((k+1)*epsi_s/k+epsi_i[i]);	
														}


													}
													else
													{
														if(k==0)
														{
															b[i*(2*p+1)*(p+1)+j*(p+1)+k]=b[i*(2*p+1)*(p+1)+j*(p+1)+k]+local[ii][i][j][k].real*k*(epsi_i[i]-epsi_s);
															b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]=b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]+local[ii][i][j][k].imag*k*(epsi_i[i]-epsi_s);
														}
														else
														{
															b[i*(2*p+1)*(p+1)+j*(p+1)+k]=b[i*(2*p+1)*(p+1)+j*(p+1)+k]+local[ii][i][j][k].real*(epsi_i[i]-epsi_s);
															b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]=b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]+local[ii][i][j][k].imag*(epsi_i[i]-epsi_s);
														}
													}

												}
											}

											/*	for(i=0;i<M;i++)
											for(j=0;j<2*p+1;j++)
											for(k=0;k<p+1;k++)
											{
											if(abs(j-p)<=k)
											printf("b=%d %d %d %f %f\n",i, k,j-p, b[i*(2*p+1)*(p+1)+j*(p+1)+k],b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]);
											}
											getchar();*/


}

void preconr(int n, double *x, double *b)
{ int j; double div = ( diag != 0.0 ? diag : 1.0);
for (j=0;j<n;j++) b[j] = x[j]/div;
}

double factorial(int n)
{
	double factn=1.0;
	double counti=0.0;
	int i;

	if(n==0)
	{
		return(factn);
	}
	else if(n>0)
	{
		for(i=1;i<n+1;i++)
		{
			counti=counti+1.0;
			factn=factn*counti;
		}

	}
	return(factn);


}

void generate_coeff( )
{
	int i,j,k,l,ii,jj,kk;
	double gamma=(epsi_i[0]-epsi_s)/(epsi_s+epsi_i[0]);

	ncoef1 = new double *[p+1];
	for(i=0;i<p+1;i++)
		ncoef1[i] = new double [p+1];

	ncoef2=new double**[p+1];
	for(i=0;i<p+1;i++)
	{ncoef2[i]=new double*[p+1];}
	for(i=0;i<p+1;i++)
		for(j=0;j<p+1;j++)
		{ncoef2[i][j]=new double[p+1];}

		ncoef3=new double***[p+1];
		for(i=0;i<p+1;i++)
		{ncoef3[i]=new double**[p+1];}
		for(i=0;i<p+1;i++)
			for(j=0;j<p+1;j++)
			{ncoef3[i][j]=new double*[p+1];}
			for(i=0;i<p+1;i++)
				for(j=0;j<p+1;j++)
					for(k=0;k<p+1;k++)
					{ncoef3[i][j][k]=new double[p+1];}

					ocoef1 = new double *[p+1];
					for(i=0;i<p+1;i++)
						ocoef1[i] = new double [p+1];

					ocoef2=new double**[p+1];
					for(i=0;i<p+1;i++)
					{ocoef2[i]=new double*[p+1];}
					for(i=0;i<p+1;i++)
						for(j=0;j<p+1;j++)
						{ocoef2[i][j]=new double[p+1];}

						ocoef3=new double***[p+1];
						for(i=0;i<p+1;i++)
						{ocoef3[i]=new double**[p+1];}
						for(i=0;i<p+1;i++)
							for(j=0;j<p+1;j++)
							{ocoef3[i][j]=new double*[p+1];}
							for(i=0;i<p+1;i++)
								for(j=0;j<p+1;j++)
									for(k=0;k<p+1;k++)
									{ocoef3[i][j][k]=new double[p+1];}

									for(l=0;l<p+1;l++)
										for(k=0;k<=l;k++)
										{
											ncoef1[l][k]=pow(-1.0,double(l+1))*sqrt(factorial(l+k)/(factorial(l-k)*factorial(2*k)))*gamma;
											//ocoef1[l][k]=0.4*pow(-1.0,double(l+2))*(sqrt(double(factorial(2*k))/(double(factorial(l-k)*factorial(l+k)))))*double(factorial(l))*/double(factorial(k));
											//printf("coefs1= %d %d %.7f\n",l,k,ncoef1[l][k]);
										}

										for(l=0;l<p+1;l++)
											for(k=0;k<=l;k++)
												for(j=k+1;j<=l;j++)
												{
													ncoef2[l][k][j]=pow(-1.0,double(l+1))*sqrt(factorial(l+j)/(factorial(l+k)*factorial(j-k)))*sqrt(factorial(l+j)/(factorial(l-k)*factorial(j+k)))*gamma;
													//ocoef2[l][k][j]=-0.4*double(factorial(k+l-j+1))*double(2*k+1)/double(factorial(2*k+l-j+1))/double(k+1);
													//	printf("coefs2=%d %d %d %.7f\n",l,k,j,ncoef2[l][k][j]);
												}


												for(l=0;l<p+1;l++)
													for(k=0;k<=l;k++)
														for(j=k+1;j<=l;j++)
															for(i=1;i<=j-k;i++)
															{
																ncoef3[l][k][j][i]=	sqrt(factorial(j+k)/(factorial(i)*factorial(j+k-i)))*sqrt(factorial(j-k)/(factorial(i)*factorial(j-k-i)));
																//	ocoef3[l][k][j][i]=	ncoef3[l][k][j][i];
																//	printf("coefs3=%d %d %d  %d %.7f\n",l,k,j,i,ncoef3[l][k][j][i]);
															}



}



void generate_imagecoef(int n, double *vecx)
{
	int i,j,k,l;
	int ii,jj,kk,ll,iii;
	double sc1=1.0, sc2=1.0;
	int ier=0;
	double x0y0z0[3];
	double xnynzn[3];
	double d,dx,dy,dz,dr,rk,rim;

	double lamda=epsi_s/(epsi_s+epsi_i[0]),gamma=(epsi_i[0]-epsi_s)/(epsi_s+epsi_i[0]);



	for(i=0;i<M;i++)
		for(ii=0;ii<M;ii++)
			for(j=0;j<2*p+1;j++)
				for(k=0;k<p+1;k++)
				{
					if(i!=ii)
					{

						mpole[i][ii][j][k].real=vecx[i*(2*p+1)*(p+1)+j*(p+1)+k];
						mpole[i][ii][j][k].imag=vecx[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k];
						local[i][ii][j][k].real=0.0;
						local[i][ii][j][k].imag=0.0;


						mpolerot[i][ii][j][k].real=0.0;
						mpolerot[i][ii][j][k].imag=0.0;
					}
				}

				for(i=0;i<M;i++)
					for(ii=0;ii<M;ii++)
					{
						x0y0z0[0]=ox[i];
						x0y0z0[1]=oy[i];
						x0y0z0[2]=oz[i];

						xnynzn[0]=ox[ii];
						xnynzn[1]=oy[ii];
						xnynzn[2]=oz[ii];
						if(i!=ii)
						{
							for(j=0;j<2*p+1;j++)
								for(k=0;k<p+1;k++)
								{


									mpole[i][ii][j][k].real=mpole[i][ii][j][k].real*sqrt(2*k+1);
									mpole[i][ii][j][k].imag=mpole[i][ii][j][k].imag*sqrt(2*k+1);

								}
								l3dmplocquadu_(&sc1,x0y0z0,mpole[i][ii],&p,&sc2,xnynzn,local[i][ii],&p,&ier);
						}
					}


					//find the close compact sphere pairs, and then for each close pair, add the image of multipole expansion. (should add the rotation matrix later.)

					/*first find all the oim, i.e., all the image multipole centers*/
					for(i=0;i<M;i++)
						for(ii=0;ii<M;ii++)
						{

							if(i!=ii&&newpowerd[i][ii][1]<orad[i]+orad[ii]+sphtol) //the creterion can be modified further.
							{
								dx=ox[i]-ox[ii];
								dy=oy[i]-oy[ii];
								dz=oz[i]-oz[ii];
								dr=sqrt(dx*dx+dy*dy+dz*dz);
								rk=orad[ii]*orad[ii]/dr;
								rim=rk/dr;
								oim[i][ii][0][0]=rim*dx+ox[ii];
								oim[i][ii][0][1]=rim*dy+oy[ii];
								oim[i][ii][0][2]=rim*dz+oz[ii]; //the Kelvin point.

								//the (im-1) legendre points.


								for(j=0;j<imm-1;j++)
								{
									oim[i][ii][j+1][0]=ox[ii]+xquad[j]*rim*dx;
									oim[i][ii][j+1][1]=oy[ii]+xquad[j]*rim*dy;
									oim[i][ii][j+1][2]=oz[ii]+xquad[j]*rim*dz;
								}

							}
						}





						/*first find all the oim_i, i.e., all the image multipole centers for the internal field*/
						/*	for(i=0;i<M;i++)
						for(ii=0;ii<M;ii++)
						{

						if(i!=ii&&powerd[1]<orad[i]+orad[ii]+sphtol) //the creterion can be modified further.
						{
						dx=ox[i]-ox[ii];
						dy=oy[i]-oy[ii];
						dz=oz[i]-oz[ii];
						dr=sqrt(dx*dx+dy*dy+dz*dz);
						rk=orad[ii]*orad[ii]/dr;
						rim=rk/dr;
						oim_i[i][ii][0][0]=dx+ox[ii];
						oim_i[i][ii][0][1]=dy+oy[ii];
						oim_i[i][ii][0][2]=dz+oz[ii]; //the Kelvin point.

						//the (im-1) legendre points.


						for(j=0;j<imm-1;j++)
						{
						oim_i[i][ii][j+1][0]=ox[ii]+dx/xquad[j];
						oim_i[i][ii][j+1][1]=oy[ii]+dy/xquad[j];
						oim_i[i][ii][j+1][2]=oz[ii]+dz/xquad[j];
						}

						}
						}*/

						for(i=0;i<M;i++)
							for(ii=0;ii<M;ii++)
							{
								if(i!=ii) //the criterion can be modified further.
								{

									dx=ox[i]-ox[ii];
									dy=oy[i]-oy[ii];
									dz=oz[i]-oz[ii];
									rotangle[i][ii]=-acos(dz/sqrt(dx*dx+dy*dy+dz*dz));

									rotviarecur3f90_(&rotangle[i][ii],&p,&p,&p,mpole[i][ii],&p,mpolerot[i][ii],&p);
								}
							}


							/*second, find all the image multipole coefficients.*/
							//----------first, for each Mlk, find Nkk--------------------//
							for(i=0;i<M;i++)
								for(ii=0;ii<M;ii++)
								{

									if(i!=ii&&newpowerd[i][ii][1]<orad[i]+orad[ii]+sphtol) //the criterion can be modified further.
									{



										for(jj=0;jj<imm;jj++)
											for(l=0;l<=p;l++)
												for(k=0;k<=l;k++)

												{

													if(jj==0)
													{
														if(i!=ii)
														{
															tempcoef[i][ii][jj][l][k][k]=mpolerot[i][ii][k+p][l].real*ncoef1[l][k]/newpowerd[i][ii][l+k+1];
															tempcoefimag[i][ii][jj][l][k][k]=mpolerot[i][ii][k+p][l].imag*ncoef1[l][k]/newpowerd[i][ii][l+k+1];
														}
														else
														{
															tempcoef[i][ii][jj][l][k][k]=mpolerot[i][ii][k+p][l].real*ncoef1[l][k]/newpowerd[i][ii][l+k+1];
															tempcoefimag[i][ii][jj][l][k][k]=mpolerot[i][ii][k+p][l].imag*ncoef1[l][k]/newpowerd[i][ii][l+k+1];
														}
													}
													else
													{
														tempcoef[i][ii][jj][l][k][k]=tempcoef[i][ii][0][l][k][k]*(-0.5)*powerxquad[jj-1][k];
														tempcoefimag[i][ii][jj][l][k][k]=tempcoefimag[i][ii][0][l][k][k]*(-0.5)*powerxquad[jj-1][k];
													}






												}


									}
								}

								//--------------for Mlk, fnd Njk, j!=k--------------//
								for(i=0;i<M;i++)
									for(ii=0;ii<M;ii++)
									{

										if(i!=ii&&newpowerd[i][ii][1]<orad[i]+orad[ii]+sphtol) //the criterion can be modified further.
										{
											//	d=sqrt((ox[i]-ox[ii])*(ox[i]-ox[ii])+(oy[i]-oy[ii])*(oy[i]-oy[ii])+(oz[i]-oz[ii])*(oz[i]-oz[ii]));

											for(jj=0;jj<imm;jj++)
												for(l=0;l<=p;l++)
													for(k=0;k<=l;k++)
														for(j=k+1;j<=l;j++)
														{



															if(jj==0)
															{
																if(i!=ii)
																{
																	tempcoef[i][ii][jj][l][k][j]=mpolerot[i][ii][k+p][l].real*ncoef2[l][k][j]/newpowerd[i][ii][l+j+1];
																	for(iii=1;iii<=j-k;iii++)
																		tempcoef[i][ii][jj][l][k][j]=tempcoef[i][ii][jj][l][k][j]-ncoef3[l][k][j][iii]*tempcoef[i][ii][jj][l][k][j-iii]/newpowerd[i][ii][iii];

																	tempcoefimag[i][ii][jj][l][k][j]=mpolerot[i][ii][k+p][l].imag*ncoef2[l][k][j]/newpowerd[i][ii][l+j+1];;
																	for(iii=1;iii<=j-k;iii++)
																		tempcoefimag[i][ii][jj][l][k][j]=tempcoefimag[i][ii][jj][l][k][j]-ncoef3[l][k][j][iii]*tempcoefimag[i][ii][jj][l][k][j-iii]/newpowerd[i][ii][iii];
																}
																else
																{
																	tempcoef[i][ii][jj][l][k][j]=mpolerot[i][ii][k+p][l].real*ncoef2[l][k][j]/newpowerd[i][ii][l+j+1];
																	for(iii=1;iii<=j-k;iii++)
																		tempcoef[i][ii][jj][l][k][j]=tempcoef[i][ii][jj][l][k][j]-ncoef3[l][k][j][iii]*tempcoef[i][ii][jj][l][k][j-iii]/newpowerd[i][ii][iii];
																	tempcoefimag[i][ii][jj][l][k][j]=mpolerot[i][ii][k+p][l].imag*ncoef2[l][k][j]/newpowerd[i][ii][l+j+1];
																	for(iii=1;iii<=j-k;iii++)
																		tempcoefimag[i][ii][jj][l][k][j]=tempcoefimag[i][ii][jj][l][k][j]-ncoef3[l][k][j][iii]*tempcoefimag[i][ii][jj][l][k][j-iii]/newpowerd[i][ii][iii];
																}
															}

															else
															{

																tempcoef[i][ii][jj][l][k][j]=tempcoef[i][ii][0][l][k][j]*(-0.5)*powerxquad[jj-1][j];

																tempcoefimag[i][ii][jj][l][k][j]=tempcoefimag[i][ii][0][l][k][j]*(-0.5)*powerxquad[jj-1][j];

															}



														}
										}


									}




									//adding all the coefficients of the same order into immpoletemp[M][M][im][2p+1][p+1];
									for(i=0;i<M;i++)
										for(ii=0;ii<M;ii++)
										{
											if(i!=ii)
											{
												for(jj=0;jj<imm;jj++)
													for(k=0;k<2*p+1;k++)
														for(j=0;j<p+1;j++)
														{
															immpoletemp[i][ii][jj][k][j].real=0.0;
															immpoletemp[i][ii][jj][k][j].imag=0.0;
															immpoletemprot[i][ii][jj][k][j].real=0.0;
															immpoletemprot[i][ii][jj][k][j].imag=0.0;
														}

														for(jj=0;jj<imm;jj++)
															for(l=0;l<=p;l++)
																for(k=0;k<=p;k++)
																	for(j=k;j<=p;j++)
																	{
																		if(jj==0)
																		{
																			immpoletemp[i][ii][jj][k+p][j].real=immpoletemp[i][ii][jj][k+p][j].real+tempcoef[i][ii][jj][l][k][j];
																			immpoletemp[i][ii][jj][k+p][j].imag=immpoletemp[i][ii][jj][k+p][j].imag+tempcoefimag[i][ii][jj][l][k][j];
																		}
																		else
																		{
																			immpoletemp[i][ii][jj][k+p][j].real=immpoletemp[i][ii][jj][k+p][j].real+tempcoef[i][ii][jj][l][k][j]*wquad[jj-1];
																			immpoletemp[i][ii][jj][k+p][j].imag=immpoletemp[i][ii][jj][k+p][j].imag+tempcoefimag[i][ii][jj][l][k][j]*wquad[jj-1];
																		}
																	}

																	for(jj=0;jj<imm;jj++)
																		for(j=0;j<=p;j++)
																			for(k=1;k<=j;k++)
																			{
																				immpoletemp[i][ii][jj][p-k][j].real=immpoletemp[i][ii][jj][k+p][j].real;
																				immpoletemp[i][ii][jj][p-k][j].imag=-1.0*immpoletemp[i][ii][jj][k+p][j].imag;
																			}

																			for(jj=0;jj<imm;jj++)
																				for(j=0;j<=p;j++)

																				{
																					//immpoletemp[i][ii][jj][p][j].real=2.0*immpoletemp[i][ii][jj][p][j].real;
																					immpoletemp[i][ii][jj][p][j].imag=0.0;
																				}






											}
										}

										//obtain the image coeffcients for the field inside the sphere.
										/*		for(i=0;i<M;i++)
										for(ii=0;ii<M;ii++)
										{

										if(i!=ii&&powerd[1]<orad[i]+orad[ii]+sphtol) //the criterion can be modified further.
										{
										//	d=sqrt((ox[i]-ox[ii])*(ox[i]-ox[ii])+(oy[i]-oy[ii])*(oy[i]-oy[ii])+(oz[i]-oz[ii])*(oz[i]-oz[ii]));

										for(jj=0;jj<imm;jj++)
										for(k=p;k<2*p+1;k++)
										for(j=0;j<p+1;j++)
										{



										if(jj==0)
										{
										if(oz[i]>oz[ii])
										{
										immpoletemp_i[i][ii][jj][k][j].real=mpole[i][ii][k][j].real*(1-gamma);
										immpoletemp_i[i][ii][jj][k][j].imag=mpole[i][ii][k][j].imag*(1-gamma);
										}
										else
										{
										immpoletemp_i[i][ii][jj][k][j].real=mpolerot[i][ii][k][j].real*(1-gamma);
										immpoletemp_i[i][ii][jj][k][j].imag=mpolerot[i][ii][k][j].imag*(1-gamma);
										}
										}

										else
										{

										if(oz[i]>oz[ii])
										{
										immpoletemp_i[i][ii][jj][k][j].real=mpole[i][ii][k][j].real*gamma*0.7*pow(xquad[jj-1],-2*k-1+2*p+j)*wquad[jj-1];
										immpoletemp_i[i][ii][jj][k][j].imag=mpole[i][ii][k][j].imag*gamma*0.7*pow(xquad[jj-1],-2*k-1+2*p+j)*wquad[jj-1];
										}
										else
										{
										immpoletemp_i[i][ii][jj][k][j].real=mpolerot[i][ii][k][j].real*gamma*0.7*pow(xquad[jj-1],-2*k-1+2*p+j)*wquad[jj-1];
										immpoletemp_i[i][ii][jj][k][j].imag=mpolerot[i][ii][k][j].imag*gamma*0.7*pow(xquad[jj-1],-2*k-1+2*p+j)*wquad[jj-1];
										}

										}

										}

										for(jj=0;jj<imm;jj++)
										for(k=1;k<=j;k++)
										for(j=0;j<p+1;j++)
										{
										immpoletemp_i[i][ii][jj][p-k][j].real=immpoletemp_i[i][ii][jj][k+p][j].real;
										immpoletemp_i[i][ii][jj][p-k][j].imag=-1.0*immpoletemp_i[i][ii][jj][k+p][j].imag;
										}

										for(jj=0;jj<imm;jj++)
										for(j=0;j<p+1;j++)
										{
										//immpoletemp_i[i][ii][jj][p][j].real=immpoletemp_i[i][ii][jj][k+p][j].real;
										immpoletemp_i[i][ii][jj][p][j].imag=0.0;
										}

										}


										}		*/





										double neangle;

										//--------rotate back all the image coefficients----------//
										for(i=0;i<M;i++)
											for(ii=0;ii<M;ii++)
												for(jj=0;jj<imm;jj++)
												{
													if(i!=ii)
													{
														neangle=-1.0*rotangle[i][ii];
														rotviarecur3f90_(&neangle,&p,&p,&p,immpoletemp[i][ii][jj],&p,immpoletemprot[i][ii][jj],&p);
														//rotviarecur3f90_(&npi,&p,&p,&p,immpoletemp_i[i][ii][jj],&p,immpoletemprot_i[i][ii][jj],&p);
													}
												}




												//-----------------transform the image to local expansion-------
												for(i=0;i<M;i++)
													for(ii=0;ii<M;ii++)
														for(jj=0;jj<imm;jj++)
														{
															if(i!=ii)
															{

																x0y0z0[0]=oim[ii][i][jj][0];
																x0y0z0[1]=oim[ii][i][jj][1];
																x0y0z0[2]=oim[ii][i][jj][2];

																xnynzn[0]=ox[ii];
																xnynzn[1]=oy[ii];
																xnynzn[2]=oz[ii];

																//		if(oz[ii]>oz[i])
																//		l3dmplocquadu_(&sc1,x0y0z0,immpoletemp[ii][i][jj],&p,&sc2,xnynzn,immpoletrans[ii][i][jj],&p,&ier);
																//else
																l3dmplocquadu_(&sc1,x0y0z0,immpoletemprot[ii][i][jj],&p,&sc2,xnynzn,immpoletrans[ii][i][jj],&p,&ier);

																/*		x0y0z0[0]=oim_i[ii][i][jj][0];
																x0y0z0[1]=oim_i[ii][i][jj][1];
																x0y0z0[2]=oim_i[ii][i][jj][2];

																xnynzn[0]=ox[i];
																xnynzn[1]=oy[i];
																xnynzn[2]=oz[i];

																if(oz[ii]>oz[i])
																l3dmplocquadu_(&sc1,x0y0z0,immpoletemp_i[ii][i][jj],&p,&sc2,xnynzn,immpoletrans_i[ii][i][jj],&p,&ier);
																else
																l3dmplocquadu_(&sc1,x0y0z0,immpoletemprot_i[ii][i][jj],&p,&sc2,xnynzn,immpoletrans_i[ii][i][jj],&p,&ier);*/

																//transform the image multipole in the external region into a multipole at the center
																x0y0z0[0]=oim[ii][i][jj][0];
																x0y0z0[1]=oim[ii][i][jj][1];
																x0y0z0[2]=oim[ii][i][jj][2];

																xnynzn[0]=ox[i];
																xnynzn[1]=oy[i];
																xnynzn[2]=oz[i];

																//		if(oz[ii]>oz[i])
																//		l3dmpmpquadu_(&sc1,x0y0z0,immpoletemp[ii][i][jj],&p,&sc2,xnynzn,immpoletransmpmp[ii][i][jj],&p,&ier);
																//else
																l3dmpmpquadu_(&sc1,x0y0z0,immpoletemprot[ii][i][jj],&p,&sc2,xnynzn,immpoletransmpmp[ii][i][jj],&p,&ier);
															}
														}

														//rescale immpoletrans--------------
														for(i=0;i<M;i++)
															for(ii=0;ii<M;ii++)
																for(jj=0;jj<imm;jj++)
																{
																	if(i!=ii)
																	{
																		for(j=0;j<2*p+1;j++)
																			for(k=0;k<p+1;k++)
																			{

																				immpoletrans[ii][i][jj][j][k].real=immpoletrans[ii][i][jj][j][k].real/sqrt(2*k+1);
																				immpoletrans[ii][i][jj][j][k].imag=immpoletrans[ii][i][jj][j][k].imag/sqrt(2*k+1);

																				//	immpoletrans_i[ii][i][jj][j][k].real=immpoletrans_i[ii][i][jj][j][k].real/sqrt(2*k+1);
																				//	immpoletrans_i[ii][i][jj][j][k].imag=immpoletrans_i[ii][i][jj][j][k].imag/sqrt(2*k+1);

																				immpoletransmpmp[ii][i][jj][j][k].real=immpoletransmpmp[ii][i][jj][j][k].real/sqrt(2*k+1);
																				immpoletransmpmp[ii][i][jj][j][k].imag=immpoletransmpmp[ii][i][jj][j][k].imag/sqrt(2*k+1);

																			}
																	}
																}


																//---------------rescale back immpole temp and immpoletemprot
																for(i=0;i<M;i++)
																	for(ii=0;ii<M;ii++)
																		for(jj=0;jj<imm;jj++)
																		{
																			if(i!=ii)
																			{
																				for(j=0;j<2*p+1;j++)
																					for(k=0;k<p+1;k++)
																					{


																						immpoletemprot[ii][i][jj][j][k].real=immpoletemprot[ii][i][jj][j][k].real/sqrt(2*k+1);
																						immpoletemprot[ii][i][jj][j][k].imag=immpoletemprot[ii][i][jj][j][k].imag/sqrt(2*k+1);


																					}
																			}

																			else
																			{
																				for(j=0;j<2*p+1;j++)
																					for(k=0;k<p+1;k++)
																					{



																						immpoletemprot[ii][i][jj][j][k].real=immpoletemprot[ii][i][jj][j][k].real/sqrt(2*k+1);
																						immpoletemprot[ii][i][jj][j][k].imag=immpoletemprot[ii][i][jj][j][k].imag/sqrt(2*k+1);

																					}
																			}




																		}

																		for(i=0;i<M;i++)
																			for(ii=0;ii<M;ii++)
																			{
																				if(i!=ii)
																				{
																					for(j=0;j<2*p+1;j++)
																						for(k=0;k<p+1;k++)
																						{

																							mpole[i][ii][j][k].real=mpole[i][ii][j][k].real/sqrt(2*k+1);
																							mpole[i][ii][j][k].imag=mpole[i][ii][j][k].imag/sqrt(2*k+1);
																							local[i][ii][j][k].real=local[i][ii][j][k].real/sqrt(2*k+1);
																							local[i][ii][j][k].imag=local[i][ii][j][k].imag/sqrt(2*k+1);

																						}
																				}
																			}

}

void initialization()
{

	int i,j,k,l,ii,jj,kk;

	tempcoef=new double*****[M];
	for(i=0;i<M;i++)
	{tempcoef[i]=new double****[M];}
	for(i=0;i<M;i++)
		for(j=0;j<M;j++)
		{tempcoef[i][j]=new double***[imm];}
		for(i=0;i<M;i++)
			for(j=0;j<M;j++)
				for(k=0;k<imm;k++)
				{tempcoef[i][j][k]=new double**[p+1];}
				for(i=0;i<M;i++)
					for(j=0;j<M;j++)
						for(k=0;k<imm;k++)
							for(l=0;l<p+1;l++)
							{tempcoef[i][j][k][l]=new double*[p+1];}
							for(i=0;i<M;i++)
								for(j=0;j<M;j++)
									for(k=0;k<imm;k++)
										for(l=0;l<p+1;l++)
											for(ii=0;ii<p+1;ii++)
											{tempcoef[i][j][k][l][ii]=new double[p+1];}

											tempcoefimag=new double*****[M];
											for(i=0;i<M;i++)
											{tempcoefimag[i]=new double****[M];}
											for(i=0;i<M;i++)
												for(j=0;j<M;j++)
												{tempcoefimag[i][j]=new double***[imm];}
												for(i=0;i<M;i++)
													for(j=0;j<M;j++)
														for(k=0;k<imm;k++)
														{tempcoefimag[i][j][k]=new double**[p+1];}
														for(i=0;i<M;i++)
															for(j=0;j<M;j++)
																for(k=0;k<imm;k++)
																	for(l=0;l<p+1;l++)
																	{tempcoefimag[i][j][k][l]=new double*[p+1];}
																	for(i=0;i<M;i++)
																		for(j=0;j<M;j++)
																			for(k=0;k<imm;k++)
																				for(l=0;l<p+1;l++)
																					for(ii=0;ii<p+1;ii++)
																					{tempcoefimag[i][j][k][l][ii]=new double[p+1];}
																					rotangle= new double *[M];
																					for(i=0;i<M;i++)
																						rotangle[i] = new double [M];
}

void generate_power()
{
	//generate the quadrature nodes for the image multipole.
	double lamda=epsi_s/(epsi_s+epsi_i[0]),beta=lamda-1.0,gamma=(epsi_i[0]-epsi_s)/(epsi_s+epsi_i[0]);
	double dx,dy,dz,dr;
	int i,j,k,l,ii,jj,kk;
	wquad = new double[imm-1];
	xquad = new double[imm-1];

	cgqf (imm-1, 1, 0, 0, -1, 1, xquad, wquad);
	for(i=0;i<imm-1;i++)
	{

		xquad[i]=pow((1.0-xquad[i])/2.0,1.0/lamda);

	}

	//given a distribution of spheres, pre-compute the power of the distance between spheres so that we do not need to compute that again and again in the GMRES iteration.
	double d;
	double a=1.0, D=2.0;
	//i=0;
	//ii=1;

	d=2*a+D;

	powerd[0]=1.0;
	for(i=0;i<2*p+1;i++)
	{
		powerd[i+1]=powerd[i]*d;
	}

	for(i=0;i<M;i++)
		for(ii=0;ii<M;ii++)
		{
			if(i!=ii)
			{
				newpowerd[i][ii][0]=1.0;
				dx=ox[i]-ox[ii];
				dy=oy[i]-oy[ii];
				dz=oz[i]-oz[ii];
				dr=sqrt(dx*dx+dy*dy+dz*dz);
				for(j=0;j<2*p+1;j++)
				{
					newpowerd[i][ii][j+1]=newpowerd[i][ii][j]*dr;
				}

			}

		}

		powerxquad = new double *[imm-1];
		for(i=0;i<imm-1;i++)
			powerxquad[i] = new double [p+2];

		for(i=0;i<imm-1;i++)
		{
			powerxquad[i][0]=1.0;
			for(j=0;j<p+1;j++)
			{
				powerxquad[i][j+1]=powerxquad[i][j]*xquad[i];
			}
		}
}



