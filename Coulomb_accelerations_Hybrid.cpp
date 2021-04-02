//Hybrid method for electrostatics force field with spherical dielectric objects
//#include"Numupbound.h"
#include"MDpara.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "itlin.h"
//#include<omp.h>

/*parameters and functions for the iterative GMRES method.*/

void matvec(int n, double *vecx, double *b); //the matrix vector product
void preconr(int n, double *x, double *b); // the preconditioner
void initialization();
extern "C" void gmres(int n, double *y, MATVEC *matvec,PRECON *preconr,\
 PRECON *preconl, double *b, struct ITLIN_OPT *opt, struct ITLIN_INFO *info);

//the Gauss-jacobi quadrature function
extern void cgqf ( int nt, int kind, double alpha, double beta, double a,\
		 double b,double t[], double wts[]);

//factorial function
double factorial(int n);

//precompute some of the constant coefficients for the
//reccurence relation of the multipole expansion coefficients
void generate_coeff();

//generate the image multipoles for a given set of multipole coefficients
void generate_imagecoef(int n, double *vecx);
void generate_power();

/*functions and structure involved with the spherical harmonics (in fortran) */

//struct complex {double real,imag;};
//define a complex structure to make the  C code consistent with fortran codes.

extern "C" void sshgrid_(int *,int *,double *,double *,int *);
extern "C" void sshgfun_(int *,int *,complex *);
extern "C" void sshcoef_(complex *,int *,int *,int *,complex *,complex *);
//extern "C" void sshevalg_(complex [][Pmax+1],int *,int *,int*,\
		complex[][2*Pmax+1][Pmax+1],complex[][2*Pmax]);
extern "C" void ssheval_(complex *,int *,double [3], complex *);
extern "C" void l3dmplocquadu_(double *, double [],complex *,int *,double *,\
				 double [], complex *, int *, int *);
//extern "C" void l3dmpmpquadu_(double *, double [],complex[][Pmax+1],int *,\
			double *, double [], complex[][Pmax+1], int *, int *);
extern "C" void lfmm3dparttarg_(int *,int *,int *,double *, int *,complex [],\
		int *,complex[],double *,int *,complex [],int *,complex *,\
		int *,double *,int *,complex[],int *,complex *);

extern "C" void lfmm3dpartself_(int *,int *,int *,double *, int *,\
complex [],int *,complex[],double *,int *,complex [],int *,complex *);

//extern "C" void rotviarecur3f90_(double *,int *,int *,int *,\
		complex[][Pmax+1],int *,complex[][Pmax+1],int *);

extern "C" void lpotfld3dallhess_dp_(int *, int *, double *, complex *,\
		double *, int *, double *, complex *, complex *, complex *);

extern "C" void l3dmpevalhessd_trunc_(double *,double [], complex *, int *,\
 double *, complex *, int *, complex *,\
 int *, complex *, double *, double *, int *);

extern "C" void l3dmpevalhessdini_(int *, double *);

double diag = 1.0; //the right-handside preconditioner diag-factor.

double pi=4.0*atan(1.0),npi=-1.0*pi;

struct ITLIN_OPT *opt;
struct ITLIN_INFO *info;

int  maxiter = 500;
int rcode;

void Coulomb_accelerations_Hybrid( int iprint )
{

	int i,j,k,l,ii,jj,kk,ind;
	double tempx,tempy,tempz,rji,powr;
	double aa,bb,cc;
	double vec[5];

	FMM_thresh=2000;
	m2ltr=sphtol*sphtol;

    	double cutoff=4.0;//for short range interaction. not used now

/*before MD starts, first initialize*/

	if(iter_indicator==0)
		initialization();

	/*allocate the GMRES solver options*/

	if(iter_indicator==0)
	{
		opt= (struct ITLIN_OPT*)malloc(sizeof(struct ITLIN_OPT));
		info= (struct ITLIN_INFO*)malloc(sizeof(struct ITLIN_INFO));
		TERM_CHECK termcheck = CheckEachIter;
//	fprintf(stdout,"\n Start of of GMRES\n\n");

		opt->tol = gmrestol; //iteration tolerence
		opt->i_max = 10; //the maximum number of iteration for restart.
		opt->termcheck = termcheck;
		opt->maxiter = maxiter;
		opt->errorlevel = Verbose;
		opt->monitorlevel = None;
		opt->datalevel = None;
		opt->errorfile = stdout;
		opt->monitorfile = NULL;
		opt->datafile =NULL;// fopen("test_gmres.data","w");
		//if (!opt->datafile) fprintf(stdout,"\n open"
			//" of test_gmres.data failed!\n");
		opt->iterfile =NULL;//fopen("test_gmres_iter.data","w");
		opt->resfile  =NULL;// fopen("test_gmres_res.data","w");
		opt->miscfile =NULL;// fopen("test_gmres_misc.data","w");

		epsi_i= new double[M];
      	        lamda_i= new double[M];
                beta_i= new double[M];
    	        gamma_i= new double[M];

		//allocate for locations and radius for the dielectric spheres.
		ox = new double[M];
		oy = new double[M];
		oz = new double[M];
		orad = new double[M];
		osigma = new double [M];
		count_nim= new int [M];

       	        l3dmpevalhessdini_(&p, scarray);//initialize the ssh hesssd
	}

	//allocate all the spherical harmonic coefficient arrays

	//spheres location
	for(i=0;i<M;i++)
	{
		ind=i+N_ion;
		ox[i]=x[ind];
		oy[i]=y[ind];
		oz[i]=z[ind];
		osigma[i]=q[ind];
//gan: does the scaled q rely on radius? //2017/10/17 answer: no.2017/12/19
      	        orad[i]=r[ind]-c_lj/2.0;
	}

	//set the dielectric constants.
	if(iter_indicator==0)
	{
		epsi_s=epsi_ion;

		for(i=0;i<N_col1;i++)
		{
            epsi_i[i]=ei1;
        }
        for(i=N_col1;i<N_col;i++)
        {
            epsi_i[i]=ei2;
        }
	}

	if(iter_indicator==0)
	{
		for(k=0;k<p+1;k++)
		{

			sqrtk[k]=sqrt(2.0*k+1.0);
		} //rescale.
		for(k=0;k<M;k++)
			count_nim[k]=0;
	}

	if(iter_indicator==0)
		sshgrid_(&nphi,&ntheta,rnodes[0][0],weights[0],&nnodes);

	//calculate all the ssh functions we will use in the ssh transform.
	if(iter_indicator==0)
		sshgfun_(&p,&ntheta,ynm[0][0]);


	/*After all the initialization, start to compute the computation time.*/
	/*clock_t tstart_total, tfinish_total;
	clock_t tstart_fmm, tfinish_fmm;
	clock_t tstart_sht, tfinish_sht;
	clock_t tstart_gmres, tfinish_gmres;
	clock_t tstart_img_generation, tfinish_img_generation;
    	clock_t tstart_energy, tfinish_energy;
    	clock_t tstart_force, tfinish_force;*/

	//double duration;
	//tstart_total = clock();

	//translate all the coordinates onto each spheres.

	//printf("after fortran call\n");
	for (i=0;i<M;i++)
		for(j=0;j<ntheta;j++)
			for(k=0;k<nphi;k++)
			{
				Mnodes[i][j][k][0]=rnodes[j][k][0]*orad[i]+ox[i];
				Mnodes[i][j][k][1]=rnodes[j][k][1]*orad[i]+oy[i];
				Mnodes[i][j][k][2]=rnodes[j][k][2]*orad[i]+oz[i];
			//	printf("rnodes=%f\n",rnodes[i][j][1]);
			}

/*First find out whether a source charge approaches the interface,
if so, add it's first level images to u1, which we compute in the following.*/

			//tstart_img_generation = clock();
			int imcount=0; //count the total number of image charges.
			double dx,dy,dz,dr,rk,rim;
			double start=0;
			double alpha=0.0;
			double end;
/*	double lamda=epsi_s/(epsi_s+epsi_i[0]),beta=lamda-1.0,\
gamma=(epsi_i[0]-epsi_s)/(epsi_s+epsi_i[0]);
if we first assume epsi_i are all the same for different dielectrics.*/

			int kind=4;
			int order=im-1;

			if(iter_indicator==0)
			{
				wwquad = new double[order];
				xxquad = new double[order];
			}

   			for (i=0;i<M;i++)
    			{
    			    lamda_i[i]=epsi_s/(epsi_s+epsi_i[i]);
     			    beta_i[i]=lamda_i[i]-1.0;
        		    gamma_i[i]=(epsi_i[i]-epsi_s)/(epsi_s+epsi_i[i]);
    			}


    for (i=0;i<N;i++)
	for(j=0;j<M;j++)
	{
		dx=x[i]-ox[j];
		dy=y[i]-oy[j];
		dz=z[i]-oz[j];
		dr=dx*dx+dy*dy+dz*dz;

		if(dr/orad[j]/orad[j]<sourcetol*sourcetol&&dr>orad[j]*orad[j])
		{
			dr=sqrt(dr);
			ionind[i][j]=1;//indicating there is image for i inside j
			rk=orad[j]*orad[j]/dr;
			rim=rk/dr;
			imx[imcount]=rim*dx+ox[j];
			imy[imcount]=rim*dy+oy[j];
			imz[imcount]=rim*dz+oz[j];
			imq[imcount]=-gamma_i[j]*orad[j]*q[i]/dr;  //the Kelvin image.
			imind[imcount]=j+N_ion;
			//sph_imind[j][count_nim[j]]=imcount;
			//count_nim[j]=count_nim[j]+1;

			/*the (im-1) line image.*/
			end=rk;
		cgqf (order, kind, alpha, beta_i[j], start, end, xxquad, wwquad);

			for(ii=0;ii<order;ii++)
			{
				ind=imcount+ii+1;
				imx[ind]=ox[j]+xxquad[ii]*dx/dr;
				imy[ind]=oy[j]+xxquad[ii]*dy/dr;
				imz[ind]=oz[j]+xxquad[ii]*dz/dr;
				imq[ind]=wwquad[ii]*gamma_i[j]*lamda_i[j]\
				*pow(rk,1.0-lamda_i[j])*q[i]/orad[j];
				imind[ind]=j+N_ion;

			//	sph_imind[j][count_nim[j]]=ind;
			//	count_nim[j]=count_nim[j]+1;
			}

			imcount=imcount+im;

			if(imcount>imnum_thresh)
				{
					printf("Sorry! The configuration"
" generates image charges more than the preset threshold number! Exit!\n");
					exit (0);
				}

		}
		else
		{ionind[i][j]=-1;}
	}

//printf("img num=%d\n",imcount);
//   tfinish_img_generation = clock();


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

			for (i=0;i<imcount;i++)
			{
			source[i][0]=imx[i];
			source[i][1]=imy[i];
			source[i][2]=imz[i];
			charge[i].real=imq[i];
			charge[i].imag=0.0;
			}

			for(i=imcount;i<imcount+N;i++)
			{
				source[i][0]=x[i-imcount];
				source[i][1]=y[i-imcount];
				source[i][2]=z[i-imcount];
				charge[i].real=q[i-imcount];
				charge[i].imag=0.0;
			}
 //   tstart_fmm = clock();

		/*If Charge number is reasonably small, do direct summation*/
			if(imcount+N<FMM_thresh||M*ntheta*nphi<FMM_thresh)
			{
			for (i=0;i<M;i++)
			for(j=0;j<ntheta;j++)
			for(k=0;k<nphi;k++)
			{
			for(ii=0;ii<imcount+N;ii++)
				{

				dx=source[ii][0]-Mnodes[i][j][k][0];
				dy=source[ii][1]-Mnodes[i][j][k][1];
				dz=source[ii][2]-Mnodes[i][j][k][2];
				dr=sqrt(dx*dx+dy*dy+dz*dz);
				fgrid[i][j][k].real=\
					fgrid[i][j][k].real+charge[ii].real/dr;

				}
			}

			for (i=0;i<M;i++)
			for(j=0;j<ntheta;j++)
			for(k=0;k<nphi;k++)
			{
			for(ii=0;ii<imcount;ii++)
				{
				if(imind[ii]==i+N_ion)
				{
					dx=imx[ii]-Mnodes[i][j][k][0];
					dy=imy[ii]-Mnodes[i][j][k][1];
					dz=imz[ii]-Mnodes[i][j][k][2];
					dr=sqrt(dx*dx+dy*dy+dz*dz);
					fgrid[i][j][k].real=\
						fgrid[i][j][k].real-imq[ii]/dr;
				}
			}

			for(ii=0;ii<N;ii++)
			{
				if(ionind[ii][i]>0||ii==i+N_ion)
				{
					dx=x[ii]-Mnodes[i][j][k][0];
					dy=y[ii]-Mnodes[i][j][k][1];
					dz=z[ii]-Mnodes[i][j][k][2];
					dr=sqrt(dx*dx+dy*dy+dz*dz);
					fgrid[i][j][k].real=\
						fgrid[i][j][k].real-q[ii]/dr;
				}
			}
			}
			}

/*-----------------another way to do direct summation------------------------*/

/*	for (i=0;i<M;i++)
		for(j=0;j<ntheta;j++)
			for(k=0;k<nphi;k++)
			{
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
						fgrid[i][j][k].real=\
					fgrid[i][j][k].real+q[ii]/dr;
					}
				}
			}*/


/*---------------If Charge number is large, use FMM------------------------*/

else
{
int Ntot=imcount+N, ier,iprec=fmmtol,ifcharge=1,\
ifdipole=0,ifpot=0,iffld=0,ntarget=M*ntheta*nphi,ifpottarg=1, iffldtarg=0;
//here iprec defines the accuracy for FMM. (1, 3digit, 2, 6digit, 3, 9digit.)

	for (i=0;i<M;i++)
	for(j=0;j<ntheta;j++)
	for(k=0;k<nphi;k++)
	{
	target[i*ntheta*nphi+j*nphi+k][0]=Mnodes[i][j][k][0];
	target[i*ntheta*nphi+j*nphi+k][1]=Mnodes[i][j][k][1];
	target[i*ntheta*nphi+j*nphi+k][2]=Mnodes[i][j][k][2];
	}

lfmm3dparttarg_(&ier,&iprec,&Ntot,source[0],&ifcharge,charge,&ifdipole,\
dipstr,dipvec[0],&ifpot,pot,&iffld,fld[0],&ntarget,target[0],&ifpottarg,\
pottarg,&iffldtarg,fldtarg[0]);

	for (i=0;i<M;i++)
	for(j=0;j<ntheta;j++)
	for(k=0;k<nphi;k++)
	{
	fgrid[i][j][k].real=pottarg[i*ntheta*nphi+j*nphi+k].real;
	}

	for (i=0;i<M;i++)
	for(j=0;j<ntheta;j++)
	for(k=0;k<nphi;k++)
	{
		for(ii=0;ii<imcount;ii++)
		{
			if(imind[ii]==i+N_ion)
			{
				dx=imx[ii]-Mnodes[i][j][k][0];
				dy=imy[ii]-Mnodes[i][j][k][1];
				dz=imz[ii]-Mnodes[i][j][k][2];
				dr=sqrt(dx*dx+dy*dy+dz*dz);
				fgrid[i][j][k].real=\
					fgrid[i][j][k].real-imq[ii]/dr;
			}
		}

	for(ii=0;ii<N;ii++)
	{
		if(ionind[ii][i]>0||ii==i+N_ion)
		{
			dx=x[ii]-Mnodes[i][j][k][0];
			dy=y[ii]-Mnodes[i][j][k][1];
			dz=z[ii]-Mnodes[i][j][k][2];
			dr=sqrt(dx*dx+dy*dy+dz*dz);
			fgrid[i][j][k].real=fgrid[i][j][k].real-q[ii]/dr;
		}
	}
}
}

     //          tfinish_fmm = clock();

//calling sshcoef to obtain all the harmonic coefficients around each spheres.
//		tstart_sht = clock();
	for(i=0;i<M;i++)
	sshcoef_(fgrid[i][0],&p,&nphi,&ntheta,ynm[0][0],ycoef[i][0]);
	//tfinish_sht = clock();

//finally, rearrange the coefficients ycoef, we obtain the RHS vector b.
	//tstart_gmres= clock();
	for(i=0;i<M;i++)
	for(j=0;j<2*p+1;j++)
	for(k=0;k<p+1;k++)
	{

//printf("test real imag %d %d %d %.26f %.26f\n",i, j-p, k,ycoef[i][j][k].real,ycoef[i][j][k].imag);
		if(k==0&&j==p)
		{
                b[i*(2*p+1)*(p+1)+j*(p+1)+k]=\
				ycoef[i][j][k].real*k*(epsi_s-epsi_i[i]);//+osigma[i];

                b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]=\
				ycoef[i][j][k].imag*k*(epsi_s-epsi_i[i]);
		}
		else
		{
                b[i*(2*p+1)*(p+1)+j*(p+1)+k]=\
			ycoef[i][j][k].real*(epsi_s-epsi_i[i])*pow(orad[i],k+1);

 		b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]=\
			ycoef[i][j][k].imag*(epsi_s-epsi_i[i])*pow(orad[i],k+1);
		}
 	}


/*we are essentially solving a Ax=b system. till now, the vector b is obtained.
In the following, we are doing the matrix vector product, i.e., A*x. Then we
can use GMRES to obtain the solution.*/


//initialize Bknm if it is the first MD iteration.
if(iter_indicator==0)
	{
		for(i=0;i<M;i++)
			for(j=0;j<2*p+1;j++)
			for(k=0;k<p+1;k++)
			{

				Bknm[i][j][k].real=b[i*(2*p+1)*(p+1)+j*(p+1)+k];
//the initial guess vector is chosen to be the same as vector b. (better guesses?)
				Bknm[i][j][k].imag=\
				b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k];

                        }
	}

//rearrange Bknm.real to obtain the initial vector Bknm1D.

for(i=0;i<M;i++)
	for(j=0;j<2*p+1;j++)
	for(k=0;k<p+1;k++)
	{
		Bknm1D[i*(2*p+1)*(p+1)+j*(p+1)+k]=Bknm[i][j][k].real;
		Bknm1D[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]=\
							Bknm[i][j][k].imag;
	}
	//perform the GMRES iteration.

	//tstart_gmres = 0;//omp_get_wtime();
	int matsize=2*M*(p+1)*(2*p+1); //the size of the matrix.
	gmres(matsize,Bknm1D,&matvec,NULL,NULL,b,opt,info);
	//		tfinish_gmres = clock();
	//		tfinish_gmres =0;//omp_get_wtime();
	//	generate_imagecoef(matsize,Bknm1D);

	//end of the computation time. We get all the coefficients and then we can compute what ever we want.


	//Check the accuracy of the obtained harmonic coeeficients Bknm1D
	// by computing the total polarization potential energy.

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
			Bknm[i][j][k].imag=\
			Bknm1D[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k];
			}
		}



    /*use BknmCopy to store all the p>=1 multipole coefficients*/
    if(p>=1)
    {
        for(i=0;i<M;i++)
            for(j=0;j<2*p+1;j++)
                for(k=0;k<p+1;k++)
                {
                    if(k>=1)
                    {
                        BknmCopy[i][j][k].real=Bknm[i][j][k].real;
                        BknmCopy[i][j][k].imag=Bknm[i][j][k].imag;
                    }
                    else
                    {
                        BknmCopy[i][j][k].real=0.0;
                        BknmCopy[i][j][k].imag=0.0;
                    }

                }
    }


/*First of all, transform Bknm into Cartesian coefficients, which
can be conveniently used in the FMM calls both in force or energy calculations*/
for(i=0;i<M;i++)
for(j=0;j<2*p+1;j++)
for(k=0;k<p+1;k++)
{

Bknm[i][j][k].real=Bknm[i][j][k].real*sqrtk[k];
Bknm[i][j][k].imag=Bknm[i][j][k].imag*sqrtk[k];
} //rescale.


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


	for(i=0;i<N;i++)
	{
		srcPosarr[3*i+0]=x[i];
		srcPosarr[3*i+1]=y[i];
		srcPosarr[3*i+2]=z[i];
	}

    	for(i=0;i<(p+1)*(p+1)*N;i++)
        	srcDenarr[i]=0.0;


	for(i=0;i<N_ion;i++)
        	srcDenarr[(p+1)*(p+1)*i]=q[i];

/*here the source or image point source have to divide by eps_s,
because eps_s is already hidden in the harmonic coefficients.*/

for(ii=0;ii<M;ii++)
{
	ind=(p+1)*(p+1)*(N_ion+ii);
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



    /*rearange the poles to form dipstr and dipvec*/

    for(i=0;i<imcount+N_ion;i++)
    {
        dipvec[i][0]=1.0;
        dipvec[i][1]=0.0;
        dipvec[i][2]=0.0;
        dipstr[i].real=0.0;
        dipstr[i].imag=0.0;
    }
    for(i=imcount+N_ion;i<imcount+N;i++)
    {
        double rescaler=sqrt(srcDenarr[4*(i-imcount)+1]*\
srcDenarr[4*(i-imcount)+1]+srcDenarr[4*(i-imcount)+2]*srcDenarr[4*(i-imcount)+2]\
+srcDenarr[4*(i-imcount)+3]*srcDenarr[4*(i-imcount)+3]);

        if(rescaler>0.000001)
        {
            dipvec[i][2]=srcDenarr[4*(i-imcount)+1]/rescaler;
            dipvec[i][0]=srcDenarr[4*(i-imcount)+2]/rescaler;
            dipvec[i][1]=srcDenarr[4*(i-imcount)+3]/rescaler;
            dipstr[i].real=rescaler;
            dipstr[i].imag=0.0;

            mpdipvec[i-imcount-N_ion][0]=dipvec[i][0];
            mpdipvec[i-imcount-N_ion][1]=dipvec[i][1];
            mpdipvec[i-imcount-N_ion][2]=dipvec[i][2];
            mpdipstr[i-imcount-N_ion].real=rescaler;
            mpdipstr[i-imcount-N_ion].imag=0.0;


        }
        else
        {
            dipvec[i][0]=1.0;
            dipvec[i][1]=0.0;
            dipvec[i][2]=0.0;
            dipstr[i].real=0.0;
            dipstr[i].imag=0.0;

            mpdipvec[i-imcount-N_ion][0]=1.;
            mpdipvec[i-imcount-N_ion][1]=0.;
            mpdipvec[i-imcount-N_ion][2]=0.;
            mpdipstr[i-imcount-N_ion].real=0.;
            mpdipstr[i-imcount-N_ion].imag=0.0;

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

/*Calculating the force using cartesian multipole expansion*/

if(force_compute==1)
{
 /*	for(i=N_ion;i<N;i++)
	srcDenarr[(p+1)*(p+1)*i]=q[i];
*/

/*If the particle number is small, directly calculate the forces*/

if(imcount+N<imnum_thresh+N)
{

for(i=0;i<3*N;i++)

pot_m[i]=0.0;  //here the pot_m means the force field, thus 3*N.

/*multipole-charge (sph+ion) interaction*/

for (int t=0; t<N; t++)
{
double p_temp[3]={0,0,0};
double tx=srcPosarr[3*t];
double ty=srcPosarr[3*t+1];
double tz=srcPosarr[3*t+2];
double scale=-1;
/*----------------------If up to p=6-------------------------------*/

int dims=(p+1)*(p+1);

for (int s=N_ion; s<N; s++)
{
double phi, fx=0.0,fy=0.0,fz=0.0;
double dX_reg=srcPosarr[3*s+0]-tx;
double dY_reg=srcPosarr[3*s+1]-ty;
double dZ_reg=srcPosarr[3*s+2]-tz;
dX_reg=dX_reg/scale;
dY_reg=dY_reg/scale;
dZ_reg=dZ_reg/scale;
double x=dX_reg;
double y=dY_reg;
double z=dZ_reg;
double x2=dX_reg*dX_reg;
double y2=dY_reg*dY_reg;
double z2=dZ_reg*dZ_reg;
double invR = (x2+y2+z2);
if (invR!=0)
{
invR = 1.0/sqrt(invR);
dr=1.0/invR;
double pm1c=dX_reg, pm1s=dY_reg, pm2c=x2-y2, pm2s=dX_reg*dY_reg,pm3c=\
dX_reg*(x2-3.0*y2),pm3s=dY_reg*(3.*x2-y2),\
pm4c=x2*x2+y2*y2-6.0*x2*y2,pm4s=pm2s*pm2c;
double pm5c=dX_reg*(x2*x2-10.0*x2*y2+5.0*y2*y2),\
pm5s=dY_reg*(y2*y2-10.0*x2*y2+5.0*x2*x2);
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
/*      if(p>=0)
{
	phi = srcDenarr[dims*s]*invR;
}
if(p>=1)
{
	phi+=srcDenarr[dims*s+1]*dZ_reg*invR3;
	phi+=srcDenarr[dims*s+2]*pm1c*invR3;
	phi+=srcDenarr[dims*s+3]*pm1s*invR3;
}

if(p>=2)
{
	phi+=srcDenarr[dims*s+4]*0.5*(3.0*z2-1.0/invR2)*invR5;
	phi+=srcDenarr[dims*s+5]*cl2*pm1c*dZ_reg*invR5;
	phi+=srcDenarr[dims*s+6]*cl2*pm1s*dZ_reg*invR5;
	phi+=srcDenarr[dims*s+7]*0.5*cl2*pm2c*invR5;
	phi+=srcDenarr[dims*s+8]*cl2*pm2s*invR5;
}

if(p>=3)
{
	phi+=srcDenarr[dims*s+9]*0.5*(5.0*z2*dZ_reg-3.0*dZ_reg/invR2)*invR7;
	phi+=srcDenarr[dims*s+10]*0.25*cl31*pm1c*(5.0*z2-1.0/invR2)*invR7;
	phi+=srcDenarr[dims*s+11]*0.25*cl31*pm1s*(5.0*z2-1.0/invR2)*invR7;
	phi+=srcDenarr[dims*s+12]*0.5*cl32*dZ_reg*pm2c*invR7;
	phi+=srcDenarr[dims*s+13]*cl32*dZ_reg*pm2s*invR7;
	phi+=srcDenarr[dims*s+14]*0.25*cl33*pm3c*invR7;
	phi+=srcDenarr[dims*s+15]*0.25*cl33*pm3s*invR7;
}

if(p>=4)
{
	phi+=srcDenarr[dims*s+16]*0.125*(8.0*z2*z2-24.*(x2+y2)*z2+3.0*(x2*x2+2.0*x2*y2+y2*y2))*invR9;
	phi+=srcDenarr[dims*s+17]*0.25*cl41*(4.0*dX_reg*dZ_reg*z2-3.0*dX_reg*dZ_reg*(x2+y2))*invR9;
	phi+=srcDenarr[dims*s+18]*0.25*cl41*(4.0*dY_reg*dZ_reg*z2-3.0*dY_reg*dZ_reg*(x2+y2))*invR9;
	phi+=srcDenarr[dims*s+19]*0.25*cl42*pm2c*(6.0*z2-x2-y2)*invR9;
	phi+=srcDenarr[dims*s+20]*0.5*cl42*pm2s*(6.0*z2-x2-y2)*invR9;
	phi+=srcDenarr[dims*s+21]*0.25*cl43*dZ_reg*pm3c*invR9;
	phi+=srcDenarr[dims*s+22]*0.25*cl43*dZ_reg*pm3s*invR9;
	phi+=srcDenarr[dims*s+23]*0.125*cl44*pm4c*invR9;
	phi+=srcDenarr[dims*s+24]*0.5*cl44*pm4s*invR9;
}

if(p>=5)
{
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
}

if(p>=6)
{
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
	phi+=srcDenarr[dims*s+48]*0.0625*cl66*pm6s*invR13;
}*/

/* if(p>=0)
{
	fx =srcDenarr[dims*s]*dX_reg*invR3;
}*/

if(p>=1)
{
	fx+=srcDenarr[dims*s+1]*3.0*dX_reg*dZ_reg*invR5;
	fx+=srcDenarr[dims*s+2]*(3.0*x2-1.0/invR2)*invR5;
	fx+=srcDenarr[dims*s+3]*3.0*dX_reg*dY_reg*invR5;
}

if(p>=2)
{
	fx+=srcDenarr[dims*s+4]*1.5*dX_reg*(5.*z2-1.0/invR2)*invR7;
	fx+=srcDenarr[dims*s+5]*cl2*dZ_reg*(5.0*x2-1.0/invR2)*invR7;
	fx+=srcDenarr[dims*s+6]*cl2*5.0*pm1c*pm1s*dZ_reg*invR7;
	fx+=srcDenarr[dims*s+7]*0.5*cl2*dX_reg*(5.0*pm2c-2.0/invR2)*invR7;
	fx+=srcDenarr[dims*s+8]*cl2*pm1s*(5.0*x2-1.0/invR2)*invR7;
}

if(p>=3)
{
	fx+=srcDenarr[dims*s+9]*dX_reg*dZ_reg*(10.0*z2-7.5*x2-7.5*y2)*invR9;
	fx+=srcDenarr[dims*s+10]*0.25*cl31*(-4*x2*x2 +y2*y2\
				- 3*y2*z2 - 4*z2*z2 - 3*x2*y2 - 9*z2)*invR9;
	fx+=srcDenarr[dims*s+11]*0.25*cl31*(-5*pm2s*(x2 + y2 - 6*z2))*invR9;
	fx+=srcDenarr[dims*s+12]*0.5*cl32*(dX_reg*dZ_reg*(5*x2-9*y2-2*z2))*invR9;
	fx+=srcDenarr[dims*s+13]*cl32*(-1.0*(dY_reg*dZ_reg*(-6*x2+y2+z2)))*invR9;
	fx+=srcDenarr[dims*s+14]*0.25*cl33*(4*x2*x2 + 3*y2*(y2 + z2)-\
							3*x2*(7*y2 + z2))*invR9;
	fx+=srcDenarr[dims*s+15]*0.25*cl33*(pm2s*(15*x2 - 13*y2 - 6*z2))*invR9;
}

if(p>=4)
{
	fx+=srcDenarr[dims*s+16]*(dX_reg*(1.875*x2*x2+1.875*y2*y2-22.5*y2*z2+\
				15.*z2*z2 + x2*(3.75*y2 - 22.5*z2)))*invR11;
	fx+=srcDenarr[dims*s+17]*0.25*cl41*(-1.0*z*(18*x2*x2-3.*y2*y2+y2*z2+\
				4*z2*z2 + x2*(15*y2 - 41*z2)))*invR11;
	fx+=srcDenarr[dims*s+18]*0.25*cl41*(-21*x*y*z*(x2 + y2-2*z2))*invR11;
	fx+=srcDenarr[dims*s+19]*0.25*cl42*(x*(-5*x2*x2 + 9*y2*y2 - 66*y2*z2-\
					12*z2*z2 + x2*(4*y2 + 46*z2)))*invR11;
	fx+=srcDenarr[dims*s+20]*0.5*cl42*(y*(-6*x2*x2+y2*y2-5*y2*z2-6*z2*z2+\
						x2*(-5*y2 + 51*z2)))*invR11;
	fx+=srcDenarr[dims*s+21]*0.25*cl43*(3*(2*x2*x2 + y2*(y2 + z2)-\
					x2*(9*y2 + z2)))*invR11;
	fx+=srcDenarr[dims*s+22]*0.25*cl43*(3*x*y*(7*x2 - 5*y2 - 2*z2))*invR11;
	fx+=srcDenarr[dims*s+23]*0.125*cl44*(x*(5*x2*x2 - 2*x2*(23*y2 + 2*z2)+\
				3*y2*(7*y2 + 4*z2)))*invR11;
	fx+=srcDenarr[dims*s+24]*0.5*cl44*(y*(6*x2*x2 + y2*(y2 + z2)-\
						x2*(11*y2 + 3*z2)))*invR11;
}

if(p>=5)
{
	fx+=srcDenarr[dims*s+25]*(x*z*(13.125*x2*x2 + 13.125*y2*y2 -52.5*y2*z2 \
				+ 21.*z2*z2 +x2*(26.25*y2 - 52.5*z2)))*invR13;
	fx+=srcDenarr[dims*s+26]*0.125*cl51*(6*x2*x2*x2 - y2*y2*y2 + \
		11*y2*y2*z2 +4*y2*z2*z2 - 8*z2*z2*z2 +x2*x2*(11*y2 - 101*z2) +\
				2*x2*(2*y2*y2 - 45*y2*z2 +58*z2*z2))*invR13;
	fx+=srcDenarr[dims*s+27]*0.125*cl51*(7*x*y*(x2*x2 + y2*y2 - 16*y2*z2 +\
	16*z2*z2 + 2*x2*(y2 - 8*z2)))*invR13;
	fx+=srcDenarr[dims*s+28]*0.125*cl52*(2*x*z*(-7*x2*x2+11*y2*y2-26*y2*z2-\
	4*z2*z2 + x2*(4*y2 + 22*z2)))*invR13;

	fx+=srcDenarr[dims*s+29]*0.25*cl52*(-2*y*z*(8*x2*x2-y2*y2+y2*z2+2*z2*z2+\
	x2*(7*y2 - 23*z2)))*invR13;

	fx+=srcDenarr[dims*s+30]*0.0625*cl53*\
	(-6*(2*x2*x2*x2+y2*y2*y2-7*y2*y2*z2-8*y2*z2*z2-x2*x2*(7*y2 + 23*z2)\
	+x2*(-8*y2*y2 + 90*y2*z2+8*z2*z2)))*invR13;

	fx+=srcDenarr[dims*s+31]*0.0625*cl53*(-6*x*y*(7*x2*x2-5*y2*y2+44*y2*z2+\
	16*z2*z2 + 2*x2*(y2 - 38*z2)))*invR13;
	fx+=srcDenarr[dims*s+32]*0.375*cl54*(x*z*(7*x2*x2 + 23*y2*y2+12*y2*z2-\
	2*x2*(29*y2 + 2*z2)))*invR13;
	fx+=srcDenarr[dims*s+33]*1.5*cl54*(y*z*(8*x2*x2 + y2*(y2 + z2)-\
	x2*(13*y2 + 3*z2)))*invR13;
	fx+=srcDenarr[dims*s+34]*0.375*cl55*(6*x2*x2*x2 - 5*y2*y2*(y2 + z2)-\
	5*x2*x2*(17*y2 + z2) +10*x2*(8*y2*y2 + 3*y2*z2))*invR13;
	fx+=srcDenarr[dims*s+35]*0.375*cl55*(x*y*(35*x2*x2 + 31*y2*y2+20*y2*z2-\
	10*x2*(11*y2 + 2*z2)))*invR13;
}

if(p>=6)
{
	fx+=srcDenarr[dims*s+36]*(x*(-2.1875*x2*x2*x2-2.1875*y2*y2*y2+\
52.5*y2*y2*z2 - 105.*y2*z2*z2 +28.*z2*z2*z2 + x2*x2*(-6.5625*y2 + 52.5*z2)+\
x2*(-6.5625*y2*y2 + 105.*y2*z2 -105.*z2*z2)))*invR15;

	fx+=srcDenarr[dims*s+37]*0.125*cl61*(z*(40*x2*x2*x2 - 5*y2*y2*y2 + \
	15*y2*y2*z2 +12*y2*z2*z2 - 8*z2*z2*z2 +75*x2*x2*(y2 - 3*z2) +6*x2*(5*y2*y2 - \
	35*y2*z2 +26*z2*z2)))*invR15;

	fx+=srcDenarr[dims*s+38]*0.125*cl61*(3*x*y*z*(15*x2*x2 + 15*y2*y2 - \
	80*y2*z2 +48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;
	fx+=srcDenarr[dims*s+39]*0.0625*cl62*(x*(7*x2*x2*x2 - 11*y2*y2*y2 + \
	210*y2*y2*z2 -240*y2*z2*z2 - 32*z2*z2*z2 +3*x2*x2*(y2 - 50*z2) -\
	15*x2*(y2*y2 - 4*y2*z2 -16*z2*z2)))*invR15;

	fx+=srcDenarr[dims*s+40]*0.125*cl62*(y*(8*x2*x2*x2-y2*y2*y2+15*y2*y2*z2-\
	16*z2*z2*z2 + 15*x2*x2*(y2 - 11*z2)+6*x2*(y2*y2-25*y2*z2+40*z2*z2)))*invR15;

	fx+=srcDenarr[dims*s+41]*0.0625*cl63*(-2*z*(24*x2*x2*x2-5*x2*x2*(15*y2+\
	19*z2)+3*y2*(3*y2*y2 - 5*y2*z2 -8*z2*z2) +\
	x2*(-90*y2*y2 + 330*y2*z2 + 24*z2*z2)))*invR15;

	fx+=srcDenarr[dims*s+42]*0.0625*cl63*(-2*x*y*z*(81*x2*x2 - 51*y2*y2 +\
	140*y2*z2 +48*z2*z2 + 30*x2*(y2 - 10*z2)))*invR15;

	fx+=srcDenarr[dims*s+43]*0.09375*cl64*(-2*x*(7*x2*x2*x2 + 23*y2*y2*y2-\
	240*y2*y2*z2 -120*y2*z2*z2 -3*x2*x2*(17*y2 + 32*z2)+\
	x2*(-35*y2*y2 + 720*y2*z2 +40*z2*z2)))*invR15;

	fx+=srcDenarr[dims*s+44]*0.375*cl64*(-2*y*(8*x2*x2*x2 + y2*y2*y2-\
	9*y2*y2*z2 -10*y2*z2*z2 -5*x2*x2*(y2 + 21*z2) -6*x2*(2*y2*y2 -\
	25*y2*z2-5*z2*z2)))*invR15;

	fx+=srcDenarr[dims*s+45]*0.375*cl65*(z*(8*x2*x2*x2-5*y2*y2*(y2+z2)+\
	30*x2*y2*(3*y2+z2)-5*x2*x2*(21*y2 + z2)))*invR15;

	fx+=srcDenarr[dims*s+46]*0.375*cl65*(x*y*z*(45*x2*x2+33*y2*y2+20*y2*z2-\
	10*x2*(13*y2 + 2*z2)))*invR15;

	fx+=srcDenarr[dims*s+47]*0.0625*cl66*(x*(7*x2*x2*x2-43*y2*y2*y2-\
	30*y2*y2*z2-3*x2*x2*(47*y2 + 2*z2)+15*x2*(15*y2*y2 + 4*y2*z2)))*invR15;

	fx+=srcDenarr[dims*s+48]*0.0625*cl66*(2*y*(24*x2*x2*x2 -\
3*y2*y2*(y2 + z2) -5*x2*x2*(23*y2 + 3*z2) +6*x2*(11*y2*y2 + 5*y2*z2)))*invR15;
}

/*   if(p>=0)
{
	fy =srcDenarr[dims*s]*dY_reg*invR3;
}*/

if(p>=1)
{
	fy+=srcDenarr[dims*s+1]*3.0*dY_reg*dZ_reg*invR5;
	fy+=srcDenarr[dims*s+2]*3.0*dX_reg*dY_reg*invR5;
	fy+=srcDenarr[dims*s+3]*(3.0*y2-1.0/invR2)*invR5;
}

if(p>=2)
{
	fy+=srcDenarr[dims*s+4]*1.5*dY_reg*(5.*z2-1.0/invR2)*invR7;
	fy+=srcDenarr[dims*s+5]*cl2*5.0*pm1c*pm1s*dZ_reg*invR7;
	fy+=srcDenarr[dims*s+6]*cl2*dZ_reg*(5.0*y2-1.0/invR2)*invR7;
	fy+=srcDenarr[dims*s+7]*0.5*cl2*dY_reg*(5.0*pm2c-2.0/invR2)*invR7;
	fy+=srcDenarr[dims*s+8]*cl2*pm1c*(5.0*y2-1.0/invR2)*invR7;
}

if(p>=3)
{
	fy+=srcDenarr[dims*s+9]*dY_reg*dZ_reg*(10.0*z2-7.5*x2-7.5*y2)*invR9;
	fy+=srcDenarr[dims*s+10]*0.25*cl31*(-5.0*pm2s*(x2 + y2 - 6*z2))*invR9;
	fy+=srcDenarr[dims*s+11]*0.25*cl31*(x2*x2-4*y2*y2+27*y2*z2-4*z2*z2 -\
	3*x2*(y2 + z2))*invR9;
	fy+=srcDenarr[dims*s+12]*0.5*cl32*(dY_reg*dZ_reg*(9*x2-5*y2+2*z2))*invR9;
	fy+=srcDenarr[dims*s+13]*cl32*(-(x*z*(x2 - 6*y2 + z2)))*invR9;
	fy+=srcDenarr[dims*s+14]*0.25*cl33*(pm2s*(13*x2 - 15*y2 + 6*z2))*invR9;
	fy+=srcDenarr[dims*s+15]*0.25*cl33*(-3*x2*x2 - 4*y2*y2+3*y2*z2+\
	3*x2*(7*y2 - z2))*invR9;
}

if(p>=4)
{
	fy+=srcDenarr[dims*s+16]*(dY_reg*(1.875*x2*x2 + 1.875*y2*y2-\
	22.5*y2*z2+15.*z2*z2 + x2*(3.75*y2 - 22.5*z2)))*invR11;
	fy+=srcDenarr[dims*s+17]*0.25*cl41*(-21*x*y*z*(x2 + y2 - 2*z2))*invR11;
	fy+=srcDenarr[dims*s+18]*0.25*cl41*(-1.0*z*(-3*x2*x2+18*y2*y2-41*y2*z2+\
	4*z2*z2 + x2*(15*y2 + z2)))*invR11;
	fy+=srcDenarr[dims*s+19]*0.25*cl42*(y*(-9*x2*x2 + 5*y2*y2 - 46*y2*z2 +\
	12*z2*z2 + x2*(-4*y2 + 66*z2)))*invR11;
	fy+=srcDenarr[dims*s+20]*0.5*cl42*(x*(x2*x2 - 6*y2*y2 + 51*y2*z2 -\
	6*z2*z2 - 5*x2*(y2 + z2)))*invR11;
	fy+=srcDenarr[dims*s+21]*0.25*cl43*(3*x*y*(5*x2 - 7*y2 + 2*z2))*invR11;
	fy+=srcDenarr[dims*s+22]*0.25*cl43*(-3*(x2*x2 + 2*y2*y2 - y2*z2 +\
	x2*(-9*y2 + z2)))*invR11;
	fy+=srcDenarr[dims*s+23]*0.125*cl44*(y*(21*x2*x2 + 5*y2*y2 - 4*y2*z2 +\
	x2*(-46*y2 + 12*z2)))*invR11;
	fy+=srcDenarr[dims*s+24]*0.5*cl44*(-1.0*x*(x2*x2 + 6*y2*y2 - 3*y2*z2 +\
	x2*(-11*y2 + z2)))*invR11;
}

if(p>=5)
{
	fy+=srcDenarr[dims*s+25]*(y*z*(13.125*x2*x2 + 13.125*y2*y2 -52.5*y2*z2+\
	21.*z2*z2 +x2*(26.25*y2 - 52.5*z2)))*invR13;
	fy+=srcDenarr[dims*s+26]*0.125*cl51*(7*x*y*(x2*x2 + y2*y2 - 16*y2*z2 +\
	16*z2*z2 + 2*x2*(y2 - 8*z2)))*invR13;
	fy+=srcDenarr[dims*s+27]*0.125*cl51*(-x2*x2*x2 + 6*y2*y2*y2 - \
	101*y2*y2*z2 +116*y2*z2*z2 - 8*z2*z2*z2 +x2*x2*(4*y2 + 11*z2) +\
	x2*(11*y2*y2 - 90*y2*z2 + 4*z2*z2))*invR13;

	fy+=srcDenarr[dims*s+28]*0.125*cl52*(y*z*(-11*x2*x2 + 7*y2*y2 -\
	22*y2*z2+4*z2*z2 + x2*(-4*y2 + 26*z2)))*invR13;
	fy+=srcDenarr[dims*s+29]*0.25*cl52*(2*x*z*(x2*x2 - 8*y2*y2 + 23*y2*z2 -\
	2*z2*z2 - x2*(7*y2 + z2)))*invR13;
	fy+=srcDenarr[dims*s+30]*0.0625*cl53*(-6*x*y*(5*x2*x2 - 7*y2*y2 + \
	76*y2*z2 -16*z2*z2 - 2*x2*(y2 + 22*z2)))*invR13;
	fy+=srcDenarr[dims*s+31]*0.0625*cl53*(6*(x2*x2*x2 + 2*y2*y2*y2 -\
	23*y2*y2*z2 +8*y2*z2*z2 -x2*x2*(8*y2 + 7*z2) +x2*(-7*y2*y2 + \
	90*y2*z2 -8*z2*z2)))*invR13;

	fy+=srcDenarr[dims*s+32]*0.375*cl54*(y*z*(23*x2*x2 + 7*y2*y2 - 4*y2*z2+\
	x2*(-58*y2 + 12*z2)))*invR13;
	fy+=srcDenarr[dims*s+33]*1.5*cl54*(-1.0*x*z*(x2*x2 + 8*y2*y2 - 3*y2*z2+\
	x2*(-13*y2 + z2)))*invR13;
	fy+=srcDenarr[dims*s+34]*0.375*cl55*(x*y*(31*x2*x2 + 5*y2*(7*y2-4*z2)+\
	x2*(-110*y2 + 20*z2)))*invR13;
	fy+=srcDenarr[dims*s+35]*0.375*cl55*(-5*x2*x2*x2 + 6*y2*y2*y2 - \
	5*y2*y2*z2+x2*x2*(80*y2 - 5*z2) +\
	x2*(-85*y2*y2 + 30*y2*z2))*invR13;
}

if(p>=6)
{
	fy+=srcDenarr[dims*s+36]*(y*(-2.1875*x2*x2*x2 - 2.1875*y2*y2*y2 +\
	52.5*y2*y2*z2 - 105.*y2*z2*z2 +28.*z2*z2*z2 +\
	x2*x2*(-6.5625*y2 + 52.5*z2) +\
	x2*(-6.5625*y2*y2 + 105.*y2*z2 -105.*z2*z2)))*invR15;

	fy+=srcDenarr[dims*s+37]*0.125*cl61*(3*x*y*z*(15*x2*x2 + 15*y2*y2 - \
	80*y2*z2 +48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;

	fy+=srcDenarr[dims*s+38]*0.125*cl61*(z*(-5*x2*x2*x2 + 40*y2*y2*y2 - \
	225*y2*y2*z2 +156*y2*z2*z2 - 8*z2*z2*z2 +15*x2*x2*(2*y2 + z2) +\
	3*x2*(25*y2*y2 - 70*y2*z2 +4*z2*z2)))*invR15;

	fy+=srcDenarr[dims*s+39]*0.0625*cl62*(y*(11*x2*x2*x2 - 7*y2*y2*y2 +\
	150*y2*y2*z2 -240*y2*z2*z2 + 32*z2*z2*z2 +15*x2*x2*(y2 - 14*z2) -\
	3*x2*(y2*y2 + 20*y2*z2 -80*z2*z2)))*invR15;

	fy+=srcDenarr[dims*s+40]*0.125*cl62*(x*(-x2*x2*x2 + 8*y2*y2*y2 -\
	165*y2*y2*z2 +240*y2*z2*z2 - 16*z2*z2*z2 +3*x2*x2*(2*y2 + 5*z2) +\
	15*x2*(y2*y2 - 10*y2*z2)))*invR15;

	fy+=srcDenarr[dims*s+41]*0.0625*cl63*(2*x*y*z*(-51*x2*x2 + 81*y2*y2 -\
	300*y2*z2 +48*z2*z2 + 10*x2*(3*y2 + 14*z2)))*invR15;

	fy+=srcDenarr[dims*s+42]*0.0625*cl63*(2*z*(9*x2*x2*x2 + 24*y2*y2*y2 -\
	95*y2*y2*z2 +24*y2*z2*z2 -15*x2*x2*(6*y2 + z2) -3*x2*(25*y2*y2 -\
	110*y2*z2 +8*z2*z2)))*invR15;

	fy+=srcDenarr[dims*s+43]*0.09375*cl64*(-2*y*(23*x2*x2*x2 + 7*y2*y2*y2-\
	96*y2*y2*z2 +40*y2*z2*z2 -5*x2*x2*(7*y2 + 48*z2) -\
	3*x2*(17*y2*y2 - 240*y2*z2 +40*z2*z2)))*invR15;

	fy+=srcDenarr[dims*s+44]*0.375*cl64*(2*x*(x2*x2*x2 + 8*y2*y2*y2 - \
	105*y2*y2*z2 +30*y2*z2*z2 -3*x2*x2*(4*y2 + 3*z2) -\
	5*x2*(y2*y2 - 30*y2*z2 +2*z2*z2)))*invR15;

	fy+=srcDenarr[dims*s+45]*0.375*cl65*(x*y*z*(33*x2*x2 + 5*y2*(9*y2 -\
	4*z2) +x2*(-130*y2 + 20*z2)))*invR15;

	fy+=srcDenarr[dims*s+46]*0.375*cl65*(z*(-5*x2*x2*x2 + 8*y2*y2*y2 -\
	5*y2*y2*z2 +x2*x2*(90*y2 - 5*z2) -15*x2*(7*y2*y2 - 2*y2*z2)))*invR15;

	fy+=srcDenarr[dims*s+47]*0.0625*cl66*(y*(43*x2*x2*x2 - 7*y2*y2*y2 +\
	6*y2*y2*z2 +x2*x2*(-225*y2 + 30*z2)+3*x2*(47*y2*y2 - 20*y2*z2)))*invR15;

	fy+=srcDenarr[dims*s+48]*0.0625*cl66*(-2*x*(3*x2*x2*x2 - 24*y2*y2*y2 +\
	15*y2*y2*z2 +x2*x2*(-66*y2 + 3*z2) +5*x2*(23*y2*y2 - 6*y2*z2)))*invR15;
}


/* if(p>=0)
{
	fz =srcDenarr[dims*s]*dZ_reg*invR3;
}*/

if(p>=1)
{
	fz+=srcDenarr[dims*s+1]*(3.0*z2-1.0/invR2)*invR5;
	fz+=srcDenarr[dims*s+2]*3.0*dX_reg*dZ_reg*invR5;
	fz+=srcDenarr[dims*s+3]*3.0*dY_reg*dZ_reg*invR5;
}

if(p>=2)
{
	fz+=srcDenarr[dims*s+4]*1.5*dZ_reg*(5.*z2-3.0/invR2)*invR7;
	fz+=srcDenarr[dims*s+5]*cl2*dX_reg*(5.0*z2-1.0/invR2)*invR7;
	fz+=srcDenarr[dims*s+6]*cl2*dY_reg*(5.0*z2-1.0/invR2)*invR7;
	fz+=srcDenarr[dims*s+7]*0.5*cl2*5.0*dZ_reg*pm2c*invR7;
	fz+=srcDenarr[dims*s+8]*cl2*5.0*pm1c*pm1s*dZ_reg*invR7;
}

if(p>=3)
{
	fz+=srcDenarr[dims*s+9]*(1.5*x2*x2+1.5*y2*y2-12.*y2*z2+4.*z2*z2+\
	z2*(3.*y2-12.*z2))*invR9;
	fz+=srcDenarr[dims*s+10]*0.25*cl31*(-5.0*dX_reg*dZ_reg*(3*x2 +\
	3*y2 - 4*x2))*invR9;
	fz+=srcDenarr[dims*s+11]*0.25*cl31*(-5*dY_reg*dZ_reg*(3*x2 + 3*y2 -\
	4*z2))*invR9;
	fz+=srcDenarr[dims*s+12]*0.5*cl32*((-x2 + y2)*(x2 + y2 -6*z2))*invR9;
	fz+=srcDenarr[dims*s+13]*cl32*(-(x*y*(x2 + y2 - 6*z2)))*invR9;
	fz+=srcDenarr[dims*s+14]*0.25*cl33*(7*dX_reg*dZ_reg*(x2 - 3*y2))*invR9;
	fz+=srcDenarr[dims*s+15]*0.25*cl33*(7*dY_reg*dZ_reg*(3*x2 - y2))*invR9;
}


if(p>=4)
{
	fz+=srcDenarr[dims*s+16]*(dZ_reg*(9.375*x2*x2+9.375*y2*y2-25.*y2*z2+\
	5.*z2*z2 + x2*(18.75*y2 - 25.*z2)))*invR11;

	fz+=srcDenarr[dims*s+17]*0.25*cl41*(3*x*(x2*x2 + y2*y2 - 12*y2*z2 +\
	8*z2*z2 + 2*x2*(y2 - 6*z2)))*invR11;

	fz+=srcDenarr[dims*s+18]*0.25*cl41*(3*y*(x2*x2 + y2*y2 - 12*y2*z2 +\
	8*z2*z2 + 2*x2*(y2 - 6*z2)))*invR11;
	fz+=srcDenarr[dims*s+19]*0.25*cl42*(-21*(x2-y2)*z*(x2+y2-2*z2))*invR11;
	fz+=srcDenarr[dims*s+20]*0.5*cl42*(-21*x*y*z*(x2 + y2 - 2*z2))*invR11;
	fz+=srcDenarr[dims*s+21]*0.25*cl43*dZ_reg*(9*(x2 - 3*y2)*z*x)*invR11;
	fz+=srcDenarr[dims*s+22]*0.25*cl43*dZ_reg*(9*(3*x2 - y2)*z*y)*invR11;
	fz+=srcDenarr[dims*s+23]*0.125*cl44*(9*(x2*x2-6*x2*y2+y2*y2)*z)*invR11;
	fz+=srcDenarr[dims*s+24]*0.5*cl44*(9*x*y*(x2-y2)*z)*invR11;
}

if(p>=5)
{
	fz+=srcDenarr[dims*s+25]*(-1.875*x2*x2*x2-1.875*y2*y2*y2+、
	33.75*y2*y2*z2-45.*y2*z2*z2 + 6.*z2*z2*z2+x2*x2*(-5.625*y2 + 33.75*z2)+\
	x2*(-5.625*y2*y2 + 67.5*y2*z2 -45.*z2*z2))*invR13;

	fz+=srcDenarr[dims*s+26]*0.125*cl51*(7*x*z*(5*x2*x2 + 5*y2*y2 - \
	20*y2*z2 +8*z2*z2 + 10*x2*(y2 - 2*z2)))*invR13;

	fz+=srcDenarr[dims*s+27]*0.125*cl51*(7*y*z*(5*x2*x2+5*y2*y2-20*y2*z2+\
	8*z2*z2 + 10*x2*(y2 - 2*z2)))*invR13;
	fz+=srcDenarr[dims*s+28]*0.125*cl52*(2*(x2 - y2)*(x2*x2 + y2*y2 -\
	16*y2*z2 + 16*z2*z2 +2*x2*(y2 - 8*z2)))*invR13;
	fz+=srcDenarr[dims*s+29]*0.25*cl52*(2*x*y*(x2*x2 + y2*y2 - 16*y2*z2 +\
	16*z2*z2 + 2*x2*(y2 - 8*z2)))*invR13;
	fz+=srcDenarr[dims*s+30]*0.0625*cl53*(-18*x*(x2 - 3*y2)*z*\
	(3*x2 + 3*y2 - 8*z2))*invR13;
	fz+=srcDenarr[dims*s+31]*0.0625*cl53*(18*y*(-3*x2 + y2)*z*\
	(3*x2 + 3*y2 - 8*z2))*invR13;
	fz+=srcDenarr[dims*s+32]*0.375*cl54*(-1.0*(x2*x2 - 6*x2*y2 + y2*y2)*\
	(x2 + y2 - 10*z2))*invR13;
	fz+=srcDenarr[dims*s+33]*1.5*cl54*(-1.0*x*y*(x2 - y2)*\
	(x2 + y2 - 10*z2))*invR13;
	fz+=srcDenarr[dims*s+34]*0.375*cl55*(11*x*(x2*x2 - 10*x2*y2 +\
	5*y2*y2)*z)*invR13;
	fz+=srcDenarr[dims*s+35]*0.375*cl55*(-5*x2*x2*x2 + 6*y2*y2*y2 -\
	5*y2*y2*z2+x2*x2*(80*y2 - 5*z2) +x2*(-85*y2*y2 + 30*y2*z2))*invR13;
}


if(p>=6)
{
	fz+=srcDenarr[dims*s+36]*(z*(-15.3125*x2*x2*x2 - 15.3125*y2*y2*y2 +\
	91.875*y2*y2*z2-73.5*y2*z2*z2+7.*z2*z2*z2+x2*x2*(-45.9375*y2 +\
	91.875*z2) +x2*(-45.9375*y2*y2 + 183.75*y2*z2 -73.5*z2*z2)))*invR15;

	fz+=srcDenarr[dims*s+37]*0.125*cl61*(3*x*y*z*(15*x2*x2 + 15*y2*y2 -\
	80*y2*z2 +48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;

	fz+=srcDenarr[dims*s+38]*0.125*cl61*(-1.0*y*(5*x2*x2*x2 + 5*y2*y2*y2 -\
	120*y2*y2*z2 +240*y2*z2*z2 - 64*z2*z2*z2 +15*x2*x2*(y2 - 8*z2) +\
	15*x2*(y2*y2 - 16*y2*z2 +16*z2*z2)))*invR15;

	fz+=srcDenarr[dims*s+39]*0.0625*cl62*(3*(x2 - y2)*z*(15*x2*x2 +\
	15*y2*y2 - 80*y2*z2 +48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;

	fz+=srcDenarr[dims*s+40]*0.125*cl62*(3*x*y*z*(15*x2*x2 + 15*y2*y2 -\
	80*y2*z2 +48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;

	fz+=srcDenarr[dims*s+41]*0.0625*cl63*(2*x*(x2 - 3*y2)*\
	(3*x2*x2 + 3*y2*y2 - 60*y2*z2 +80*z2*z2 + 6*x2*(y2 - 10*z2)))*invR15;

	fz+=srcDenarr[dims*s+42]*0.0625*cl63*(-2*y*(-3*x2 + y2)*\
	(3*x2*x2 + 3*y2*y2 - 60*y2*z2 +80*z2*z2 + 6*x2*(y2 - 10*z2)))*invR15;

	fz+=srcDenarr[dims*s+43]*0.09375*cl64*(-22*(x2*x2 - 6*x2*y2 + y2*y2)*z*\
	(3*x2 + 3*y2 - 10*z2))*invR15;

	fz+=srcDenarr[dims*s+44]*0.375*cl64*(-22*x*y*(x2 - y2)*z*\
	(3*x2 + 3*y2 - 10*z2))*invR15;

	fz+=srcDenarr[dims*s+45]*0.375*cl65*(-1.0*x*(x2*x2 - 10*x2*y2+5*y2*y2)*\
	(x2 + y2 - 12*z2))*invR15;

	fz+=srcDenarr[dims*s+46]*0.375*cl65*(-1.0*y*(5*x2*x2 - 10*x2*y2+y2*y2)*\
	(x2 + y2 - 12*z2))*invR15;

	fz+=srcDenarr[dims*s+47]*0.0625*cl66*(13*(x2*x2*x2 - 15*x2*x2*y2 +\
	15*x2*y2*y2-y2*y2*y2)*z)*invR15;
	fz+=srcDenarr[dims*s+48]*0.0625*cl66*(26*x*y*(3*x2*x2 - 10*x2*y2 +\
	3*y2*y2)*z)*invR15;
}

			//p[0] += phi;

		if(t<N_ion)
		{
			p_temp[0] += fx*srcDenarr[dims*t];
			p_temp[1] += fy*srcDenarr[dims*t];
			p_temp[2] += fz*srcDenarr[dims*t];

			//	pot_m[4*s] -= phi;
			pot_m[3*s+0] -= fx*srcDenarr[dims*t];
			pot_m[3*s+1] -= fy*srcDenarr[dims*t];
			pot_m[3*s+2] -= fz*srcDenarr[dims*t];
		}
		else
		{
			p_temp[0] += fx*srcDenarr[dims*t];
			p_temp[1] += fy*srcDenarr[dims*t];
			p_temp[2] += fz*srcDenarr[dims*t];

			//	pot_m[4*s] -= phi;
			pot_m[3*s+0] -= fx*srcDenarr[dims*t];
			pot_m[3*s+1] -= fy*srcDenarr[dims*t];
			pot_m[3*s+2] -= fz*srcDenarr[dims*t];
		}



		}
	}
	//	pot_m[4*t] += p[0];
	pot_m[3*t+0] += p_temp[0];
	pot_m[3*t+1] += p_temp[1];
	pot_m[3*t+2] += p_temp[2];

}







/*image-multipole & (image-sphere)  interaction*/
for (int t=0; t<imcount; t++)
{
double p_temp[3]={0,0,0};
double tx=imx[t];
double ty=imy[t];
double tz=imz[t];

double scale=-1;
/*----------------------If up to p=6-------------------------------*/

int dims=(p+1)*(p+1);
for (int s=N_ion; s<N; s++)
{
double phi, fx=0.0, fy=0.0,fz=0.0;
double dX_reg=srcPosarr[3*s+0]-tx;
double dY_reg=srcPosarr[3*s+1]-ty;
double dZ_reg=srcPosarr[3*s+2]-tz;
dX_reg=dX_reg/scale;
dY_reg=dY_reg/scale;
dZ_reg=dZ_reg/scale;
double x=dX_reg;
double y=dY_reg;
double z=dZ_reg;
double x2=dX_reg*dX_reg;
double y2=dY_reg*dY_reg;
double z2=dZ_reg*dZ_reg;
double invR = (x2+y2+z2);
if (invR!=0.0&&imind[t]!=s)
{
invR = 1.0/sqrt(invR);
dr=1.0/invR;
double pm1c=dX_reg, pm1s=dY_reg, pm2c=x2-y2, pm2s=dX_reg*dY_reg,\
pm3c=dX_reg*(x2-3.0*y2),pm3s=dY_reg*(3.*x2-y2),\
pm4c=x2*x2+y2*y2-6.0*x2*y2,pm4s=pm2s*pm2c;
double pm5c=dX_reg*(x2*x2-10.0*x2*y2+5.0*y2*y2),\
pm5s=dY_reg*(y2*y2-10.0*x2*y2+5.0*x2*x2);
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

/*     if(p>=0)
{
phi = srcDenarr[dims*s]*invR;
}
if(p>=1)
{   //order: z-》x-》y
phi+=srcDenarr[dims*s+1]*dZ_reg*invR3;
phi+=srcDenarr[dims*s+2]*pm1c*invR3;
phi+=srcDenarr[dims*s+3]*pm1s*invR3;
}

if(p>=2)
{
phi+=srcDenarr[dims*s+4]*0.5*(3.0*z2-1.0/invR2)*invR5;
phi+=srcDenarr[dims*s+5]*cl2*pm1c*dZ_reg*invR5;
phi+=srcDenarr[dims*s+6]*cl2*pm1s*dZ_reg*invR5;
phi+=srcDenarr[dims*s+7]*0.5*cl2*pm2c*invR5;
phi+=srcDenarr[dims*s+8]*cl2*pm2s*invR5;
}

if(p>=3)
{
phi+=srcDenarr[dims*s+9]*0.5*(5.0*z2*dZ_reg-3.0*dZ_reg/invR2)*invR7;
phi+=srcDenarr[dims*s+10]*0.25*cl31*pm1c*(5.0*z2-1.0/invR2)*invR7;
phi+=srcDenarr[dims*s+11]*0.25*cl31*pm1s*(5.0*z2-1.0/invR2)*invR7;
phi+=srcDenarr[dims*s+12]*0.5*cl32*dZ_reg*pm2c*invR7;
phi+=srcDenarr[dims*s+13]*cl32*dZ_reg*pm2s*invR7;
phi+=srcDenarr[dims*s+14]*0.25*cl33*pm3c*invR7;
phi+=srcDenarr[dims*s+15]*0.25*cl33*pm3s*invR7;
}

if(p>=4)
{
phi+=srcDenarr[dims*s+16]*0.125*(8.0*z2*z2-24.*(x2+y2)*z2+3.0*(x2*x2+2.0*x2*y2+y2*y2))*invR9;
phi+=srcDenarr[dims*s+17]*0.25*cl41*(4.0*dX_reg*dZ_reg*z2-3.0*dX_reg*dZ_reg*(x2+y2))*invR9;
phi+=srcDenarr[dims*s+18]*0.25*cl41*(4.0*dY_reg*dZ_reg*z2-3.0*dY_reg*dZ_reg*(x2+y2))*invR9;
phi+=srcDenarr[dims*s+19]*0.25*cl42*pm2c*(6.0*z2-x2-y2)*invR9;
phi+=srcDenarr[dims*s+20]*0.5*cl42*pm2s*(6.0*z2-x2-y2)*invR9;
phi+=srcDenarr[dims*s+21]*0.25*cl43*dZ_reg*pm3c*invR9;
phi+=srcDenarr[dims*s+22]*0.25*cl43*dZ_reg*pm3s*invR9;
phi+=srcDenarr[dims*s+23]*0.125*cl44*pm4c*invR9;
phi+=srcDenarr[dims*s+24]*0.5*cl44*pm4s*invR9;
}

if(p>=5)
{
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
}

if(p>=6)
{
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
phi+=srcDenarr[dims*s+48]*0.0625*cl66*pm6s*invR13;
}*/

/*  if(p>=0)
{
fx =srcDenarr[dims*s]*dX_reg*invR3;
}*/

if(p>=1)
{
fx+=srcDenarr[dims*s+1]*3.0*dX_reg*dZ_reg*invR5;
fx+=srcDenarr[dims*s+2]*(3.0*x2-1.0/invR2)*invR5;
fx+=srcDenarr[dims*s+3]*3.0*dX_reg*dY_reg*invR5;
}

if(p>=2)
{
fx+=srcDenarr[dims*s+4]*1.5*dX_reg*(5.*z2-1.0/invR2)*invR7;
fx+=srcDenarr[dims*s+5]*cl2*dZ_reg*(5.0*x2-1.0/invR2)*invR7;
fx+=srcDenarr[dims*s+6]*cl2*5.0*pm1c*pm1s*dZ_reg*invR7;
fx+=srcDenarr[dims*s+7]*0.5*cl2*dX_reg*(5.0*pm2c-2.0/invR2)*invR7;
fx+=srcDenarr[dims*s+8]*cl2*pm1s*(5.0*x2-1.0/invR2)*invR7;
}

if(p>=3)
{
fx+=srcDenarr[dims*s+9]*dX_reg*dZ_reg*(10.0*z2-7.5*x2-7.5*y2)*invR9;
fx+=srcDenarr[dims*s+10]*0.25*cl31*(-4*x2*x2 +y2*y2-3*y2*z2-4*z2*z2-\
3*x2*y2-9*z2)*invR9;
fx+=srcDenarr[dims*s+11]*0.25*cl31*(-5*pm2s*(x2 + y2 - 6*z2))*invR9;
fx+=srcDenarr[dims*s+12]*0.5*cl32*(dX_reg*dZ_reg*(5*x2 - 9*y2 - 2*z2))*invR9;
fx+=srcDenarr[dims*s+13]*cl32*(-1.0*(dY_reg*dZ_reg*(-6*x2 + y2 + z2)))*invR9;
fx+=srcDenarr[dims*s+14]*0.25*cl33*(4*x2*x2+3*y2*(y2+z2)-3*x2*(7*y2+z2))*invR9;
fx+=srcDenarr[dims*s+15]*0.25*cl33*(pm2s*(15*x2-13*y2-6*z2))*invR9;
}

if(p>=4)
{
fx+=srcDenarr[dims*s+16]*(dX_reg*(1.875*x2*x2 + 1.875*y2*y2 - 22.5*y2*z2 +\
15.*z2*z2 + x2*(3.75*y2 - 22.5*z2)))*invR11;
fx+=srcDenarr[dims*s+17]*0.25*cl41*(-1.0*z*(18*x2*x2 - 3*y2*y2 + y2*z2 +\
4*z2*z2 + x2*(15*y2 - 41*z2)))*invR11;
fx+=srcDenarr[dims*s+18]*0.25*cl41*(-21*x*y*z*(x2 + y2 - 2*z2))*invR11;
fx+=srcDenarr[dims*s+19]*0.25*cl42*(x*(-5*x2*x2 + 9*y2*y2 - 66*y2*z2 -\
12*z2*z2 + x2*(4*y2 + 46*z2)))*invR11;
fx+=srcDenarr[dims*s+20]*0.5*cl42*(y*(-6*x2*x2 + y2*y2 - 5*y2*z2 -\
6*z2*z2 + x2*(-5*y2 + 51*z2)))*invR11;
fx+=srcDenarr[dims*s+21]*0.25*cl43*(3*(2*x2*x2 + y2*(y2 + z2) -\
x2*(9*y2 + z2)))*invR11;
fx+=srcDenarr[dims*s+22]*0.25*cl43*(3*x*y*(7*x2 - 5*y2 - 2*z2))*invR11;
fx+=srcDenarr[dims*s+23]*0.125*cl44*(x*(5*x2*x2 - 2*x2*(23*y2 + 2*z2) +\
3*y2*(7*y2 + 4*z2)))*invR11;
fx+=srcDenarr[dims*s+24]*0.5*cl44*(y*(6*x2*x2 + y2*(y2 + z2) -\
x2*(11*y2 + 3*z2)))*invR11;
}

if(p>=5)
{
fx+=srcDenarr[dims*s+25]*(x*z*(13.125*x2*x2 + 13.125*y2*y2 -\
52.5*y2*z2 + 21.*z2*z2 +x2*(26.25*y2 - 52.5*z2)))*invR13;
fx+=srcDenarr[dims*s+26]*0.125*cl51*(6*x2*x2*x2 - y2*y2*y2 + 11*y2*y2*z2 +\
4*y2*z2*z2 - 8*z2*z2*z2 +x2*x2*(11*y2 - 101*z2) +2*x2*(2*y2*y2 - 45*y2*z2 +\
58*z2*z2))*invR13;
fx+=srcDenarr[dims*s+27]*0.125*cl51*(7*x*y*(x2*x2 + y2*y2 - 16*y2*z2 +\
16*z2*z2 + 2*x2*(y2 - 8*z2)))*invR13;
fx+=srcDenarr[dims*s+28]*0.125*cl52*(2*x*z*(-7*x2*x2 + 11*y2*y2 - 26*y2*z2 -\
4*z2*z2 + x2*(4*y2 + 22*z2)))*invR13;
fx+=srcDenarr[dims*s+29]*0.25*cl52*(-2*y*z*(8*x2*x2 - y2*y2 + y2*z2 +\
2*z2*z2 + x2*(7*y2 - 23*z2)))*invR13;
fx+=srcDenarr[dims*s+30]*0.0625*cl53*(-6*(2*x2*x2*x2 + y2*y2*y2 - 7*y2*y2*z2 -\
8*y2*z2*z2 -x2*x2*(7*y2 + 23*z2) +x2*(-8*y2*y2 + 90*y2*z2 +8*z2*z2)))*invR13;

fx+=srcDenarr[dims*s+31]*0.0625*cl53*(-6*x*y*(7*x2*x2 - 5*y2*y2 + 44*y2*z2 +\
16*z2*z2 + 2*x2*(y2 - 38*z2)))*invR13;

fx+=srcDenarr[dims*s+32]*0.375*cl54*(x*z*(7*x2*x2 + 23*y2*y2 + 12*y2*z2 -\
2*x2*(29*y2 + 2*z2)))*invR13;

fx+=srcDenarr[dims*s+33]*1.5*cl54*(y*z*(8*x2*x2 + y2*(y2 + z2) -\
x2*(13*y2 + 3*z2)))*invR13;
fx+=srcDenarr[dims*s+34]*0.375*cl55*(6*x2*x2*x2 - 5*y2*y2*(y2 + z2) -\
5*x2*x2*(17*y2 + z2) +10*x2*(8*y2*y2 + 3*y2*z2))*invR13;

fx+=srcDenarr[dims*s+35]*0.375*cl55*(x*y*(35*x2*x2 + 31*y2*y2 + 20*y2*z2 -\
10*x2*(11*y2 + 2*z2)))*invR13;
}

if(p>=6)
{
fx+=srcDenarr[dims*s+36]*(x*(-2.1875*x2*x2*x2 - 2.1875*y2*y2*y2 +\
52.5*y2*y2*z2 - 105.*y2*z2*z2 +28.*z2*z2*z2 + x2*x2*(-6.5625*y2 + 52.5*z2) +\
x2*(-6.5625*y2*y2 + 105.*y2*z2 -105.*z2*z2)))*invR15;

fx+=srcDenarr[dims*s+37]*0.125*cl61*(z*(40*x2*x2*x2 - 5*y2*y2*y2 + 15*y2*y2*z2 +\
12*y2*z2*z2 - 8*z2*z2*z2 +75*x2*x2*(y2 - 3*z2)+6*x2*(5*y2*y2 - 35*y2*z2 +\
26*z2*z2)))*invR15;

fx+=srcDenarr[dims*s+38]*0.125*cl61*(3*x*y*z*(15*x2*x2 + 15*y2*y2 - 80*y2*z2 +\
48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;

fx+=srcDenarr[dims*s+39]*0.0625*cl62*(x*(7*x2*x2*x2-11*y2*y2*y2+210*y2*y2*z2-\
240*y2*z2*z2 - 32*z2*z2*z2 +3*x2*x2*(y2 - 50*z2)-15*x2*(y2*y2 - 4*y2*z2 -\
16*z2*z2)))*invR15;

fx+=srcDenarr[dims*s+40]*0.125*cl62*(y*(8*x2*x2*x2 - y2*y2*y2 + 15*y2*y2*z2 -\
16*z2*z2*z2 + 15*x2*x2*(y2 - 11*z2) +6*x2*(y2*y2-25*y2*z2+40*z2*z2)))*invR15;

fx+=srcDenarr[dims*s+41]*0.0625*cl63*(-2*z*(24*x2*x2*x2-5*x2*x2*(15*y2+19*z2)+\
3*y2*(3*y2*y2 - 5*y2*z2-8*z2*z2)+x2*(-90*y2*y2 + 330*y2*z2 + 24*z2*z2)))*invR15;

fx+=srcDenarr[dims*s+42]*0.0625*cl63*(-2*x*y*z*(81*x2*x2 - 51*y2*y2 +140*y2*z2+\
48*z2*z2 + 30*x2*(y2 - 10*z2)))*invR15;

fx+=srcDenarr[dims*s+43]*0.09375*cl64*(-2*x*(7*x2*x2*x2 + 23*y2*y2*y2 - \
240*y2*y2*z2-120*y2*z2*z2-3*x2*x2*(17*y2+32*z2)+x2*(-35*y2*y2 + 720*y2*z2 +\
40*z2*z2)))*invR15;

fx+=srcDenarr[dims*s+44]*0.375*cl64*(-2*y*(8*x2*x2*x2 + y2*y2*y2 - 9*y2*y2*z2 -\
10*y2*z2*z2 -5*x2*x2*(y2 + 21*z2)-6*x2*(2*y2*y2-25*y2*z2-5*z2*z2)))*invR15;

fx+=srcDenarr[dims*s+45]*0.375*cl65*(z*(8*x2*x2*x2 - 5*y2*y2*(y2 + z2) +\
30*x2*y2*(3*y2 + z2) -5*x2*x2*(21*y2 + z2)))*invR15;

fx+=srcDenarr[dims*s+46]*0.375*cl65*(x*y*z*(45*x2*x2 + 33*y2*y2 + 20*y2*z2 -\
10*x2*(13*y2 + 2*z2)))*invR15;

fx+=srcDenarr[dims*s+47]*0.0625*cl66*(x*(7*x2*x2*x2 - 43*y2*y2*y2 -30*y2*y2*z2-\
3*x2*x2*(47*y2 + 2*z2) +15*x2*(15*y2*y2 + 4*y2*z2)))*invR15;

fx+=srcDenarr[dims*s+48]*0.0625*cl66*(2*y*(24*x2*x2*x2 - 3*y2*y2*(y2 + z2) -\
5*x2*x2*(23*y2 + 3*z2) +6*x2*(11*y2*y2 + 5*y2*z2)))*invR15;
}

/*     if(p>=0)
{
fy =srcDenarr[dims*s]*dY_reg*invR3;
}*/

if(p>=1)
{
fy+=srcDenarr[dims*s+1]*3.0*dY_reg*dZ_reg*invR5;
fy+=srcDenarr[dims*s+2]*3.0*dX_reg*dY_reg*invR5;
fy+=srcDenarr[dims*s+3]*(3.0*y2-1.0/invR2)*invR5;
}

if(p>=2)
{
fy+=srcDenarr[dims*s+4]*1.5*dY_reg*(5.*z2-1.0/invR2)*invR7;
fy+=srcDenarr[dims*s+5]*cl2*5.0*pm1c*pm1s*dZ_reg*invR7;
fy+=srcDenarr[dims*s+6]*cl2*dZ_reg*(5.0*y2-1.0/invR2)*invR7;
fy+=srcDenarr[dims*s+7]*0.5*cl2*dY_reg*(5.0*pm2c-2.0/invR2)*invR7;
fy+=srcDenarr[dims*s+8]*cl2*pm1c*(5.0*y2-1.0/invR2)*invR7;
}

if(p>=3)
{
fy+=srcDenarr[dims*s+9]*dY_reg*dZ_reg*(10.0*z2-7.5*x2-7.5*y2)*invR9;
fy+=srcDenarr[dims*s+10]*0.25*cl31*(-5.0*pm2s*(x2 + y2 - 6*z2))*invR9;
fy+=srcDenarr[dims*s+11]*0.25*cl31*(x2*x2 - 4*y2*y2 + 27*y2*z2 - 4*z2*z2 -\
3*x2*(y2 + z2))*invR9;
fy+=srcDenarr[dims*s+12]*0.5*cl32*(dY_reg*dZ_reg*(9*x2 - 5*y2 + 2*z2))*invR9;
fy+=srcDenarr[dims*s+13]*cl32*(-(x*z*(x2 - 6*y2 + z2)))*invR9;
fy+=srcDenarr[dims*s+14]*0.25*cl33*(pm2s*(13*x2 - 15*y2 + 6*z2))*invR9;
fy+=srcDenarr[dims*s+15]*0.25*cl33*(-3*x2*x2 - 4*y2*y2 + 3*y2*z2 +\
3*x2*(7*y2 - z2))*invR9;
}

if(p>=4)
{
fy+=srcDenarr[dims*s+16]*(dY_reg*(1.875*x2*x2 + 1.875*y2*y2 - 22.5*y2*z2 +\
15.*z2*z2 + x2*(3.75*y2 - 22.5*z2)))*invR11;
fy+=srcDenarr[dims*s+17]*0.25*cl41*(-21*x*y*z*(x2 + y2 - 2*z2))*invR11;
fy+=srcDenarr[dims*s+18]*0.25*cl41*(-1.0*z*(-3*x2*x2 + 18*y2*y2 - 41*y2*z2 +\
4*z2*z2 + x2*(15*y2 + z2)))*invR11;
fy+=srcDenarr[dims*s+19]*0.25*cl42*(y*(-9*x2*x2 + 5*y2*y2 - 46*y2*z2 +\
12*z2*z2 + x2*(-4*y2 + 66*z2)))*invR11;
fy+=srcDenarr[dims*s+20]*0.5*cl42*(x*(x2*x2 - 6*y2*y2 + 51*y2*z2 -\
6*z2*z2 - 5*x2*(y2 + z2)))*invR11;
fy+=srcDenarr[dims*s+21]*0.25*cl43*(3*x*y*(5*x2 - 7*y2 + 2*z2))*invR11;
fy+=srcDenarr[dims*s+22]*0.25*cl43*(-3*(x2*x2 + 2*y2*y2 - y2*z2 +\
x2*(-9*y2 + z2)))*invR11;
fy+=srcDenarr[dims*s+23]*0.125*cl44*(y*(21*x2*x2 + 5*y2*y2 - 4*y2*z2 +\
x2*(-46*y2 + 12*z2)))*invR11;
fy+=srcDenarr[dims*s+24]*0.5*cl44*(-1.0*x*(x2*x2 + 6*y2*y2 - 3*y2*z2 +\
x2*(-11*y2 + z2)))*invR11;
}

if(p>=5)
{
fy+=srcDenarr[dims*s+25]*(y*z*(13.125*x2*x2 + 13.125*y2*y2 -52.5*y2*z2 + \
21.*z2*z2 +x2*(26.25*y2 - 52.5*z2)))*invR13;

fy+=srcDenarr[dims*s+26]*0.125*cl51*(7*x*y*(x2*x2 + y2*y2 - 16*y2*z2 +\
16*z2*z2 + 2*x2*(y2 - 8*z2)))*invR13;

fy+=srcDenarr[dims*s+27]*0.125*cl51*(-x2*x2*x2 + 6*y2*y2*y2 - 101*y2*y2*z2 +\
116*y2*z2*z2 - 8*z2*z2*z2 +x2*x2*(4*y2 + 11*z2)+\
x2*(11*y2*y2 - 90*y2*z2 + 4*z2*z2))*invR13;

fy+=srcDenarr[dims*s+28]*0.125*cl52*(y*z*(-11*x2*x2 + 7*y2*y2 - 22*y2*z2 +\
4*z2*z2 + x2*(-4*y2 + 26*z2)))*invR13;

fy+=srcDenarr[dims*s+29]*0.25*cl52*(2*x*z*(x2*x2 - 8*y2*y2 + 23*y2*z2 -\
2*z2*z2 - x2*(7*y2 + z2)))*invR13;

fy+=srcDenarr[dims*s+30]*0.0625*cl53*(-6*x*y*(5*x2*x2 - 7*y2*y2 + 76*y2*z2 -\
16*z2*z2 - 2*x2*(y2 + 22*z2)))*invR13;

fy+=srcDenarr[dims*s+31]*0.0625*cl53*(6*(x2*x2*x2 + 2*y2*y2*y2 - 23*y2*y2*z2 +\
8*y2*z2*z2 -x2*x2*(8*y2 + 7*z2) +x2*(-7*y2*y2 + 90*y2*z2 -8*z2*z2)))*invR13;

fy+=srcDenarr[dims*s+32]*0.375*cl54*(y*z*(23*x2*x2 + 7*y2*y2 - 4*y2*z2 +\
x2*(-58*y2 + 12*z2)))*invR13;

fy+=srcDenarr[dims*s+33]*1.5*cl54*(-1.0*x*z*(x2*x2 + 8*y2*y2 - 3*y2*z2 +\
x2*(-13*y2 + z2)))*invR13;

fy+=srcDenarr[dims*s+34]*0.375*cl55*(x*y*(31*x2*x2 + 5*y2*(7*y2 - 4*z2) +\
x2*(-110*y2 + 20*z2)))*invR13;

fy+=srcDenarr[dims*s+35]*0.375*cl55*(-5*x2*x2*x2 + 6*y2*y2*y2 - 5*y2*y2*z2 +\
x2*x2*(80*y2 - 5*z2) +x2*(-85*y2*y2 + 30*y2*z2))*invR13;
}

if(p>=6)
{
fy+=srcDenarr[dims*s+36]*(y*(-2.1875*x2*x2*x2 - 2.1875*y2*y2*y2 +52.5*y2*y2*z2-\
105.*y2*z2*z2 +28.*z2*z2*z2 + x2*x2*(-6.5625*y2 + 52.5*z2) +x2*(-6.5625*y2*y2 +\
105.*y2*z2 -105.*z2*z2)))*invR15;

fy+=srcDenarr[dims*s+37]*0.125*cl61*(3*x*y*z*(15*x2*x2 + 15*y2*y2 - 80*y2*z2 +\
48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;

fy+=srcDenarr[dims*s+38]*0.125*cl61*(z*(-5*x2*x2*x2+40*y2*y2*y2-225*y2*y2*z2+\
156*y2*z2*z2 - 8*z2*z2*z2 +15*x2*x2*(2*y2 + z2) +3*x2*(25*y2*y2 - 70*y2*z2 +\
4*z2*z2)))*invR15;

fy+=srcDenarr[dims*s+39]*0.0625*cl62*(y*(11*x2*x2*x2 - 7*y2*y2*y2 +\
150*y2*y2*z2 -240*y2*z2*z2 + 32*z2*z2*z2 +15*x2*x2*(y2 - 14*z2) -3*x2*(y2*y2 +\
20*y2*z2 -80*z2*z2)))*invR15;

fy+=srcDenarr[dims*s+40]*0.125*cl62*(x*(-x2*x2*x2 + 8*y2*y2*y2 - 165*y2*y2*z2 +\
240*y2*z2*z2-16*z2*z2*z2+3*x2*x2*(2*y2+5*z2)+15*x2*(y2*y2-10*y2*z2)))*invR15;

fy+=srcDenarr[dims*s+41]*0.0625*cl63*(2*x*y*z*(-51*x2*x2+81*y2*y2-300*y2*z2 +\
48*z2*z2 + 10*x2*(3*y2 + 14*z2)))*invR15;

fy+=srcDenarr[dims*s+42]*0.0625*cl63*(2*z*(9*x2*x2*x2+24*y2*y2*y2-95*y2*y2*z2+\
24*y2*z2*z2 -15*x2*x2*(6*y2 + z2)-3*x2*(25*y2*y2 - 110*y2*z2+8*z2*z2)))*invR15;

fy+=srcDenarr[dims*s+43]*0.09375*cl64*(-2*y*(23*x2*x2*x2 + 7*y2*y2*y2 -\
96*y2*y2*z2 +40*y2*z2*z2 -5*x2*x2*(7*y2 + 48*z2) -3*x2*(17*y2*y2 - 240*y2*z2 +\
40*z2*z2)))*invR15;

fy+=srcDenarr[dims*s+44]*0.375*cl64*(2*x*(x2*x2*x2 + 8*y2*y2*y2 - 105*y2*y2*z2+\
30*y2*z2*z2 -3*x2*x2*(4*y2 + 3*z2) -5*x2*(y2*y2 - 30*y2*z2 +2*z2*z2)))*invR15;

fy+=srcDenarr[dims*s+45]*0.375*cl65*(x*y*z*(33*x2*x2 + 5*y2*(9*y2 - 4*z2) +\
x2*(-130*y2 + 20*z2)))*invR15;

fy+=srcDenarr[dims*s+46]*0.375*cl65*(z*(-5*x2*x2*x2 + 8*y2*y2*y2 - 5*y2*y2*z2+\
x2*x2*(90*y2 - 5*z2)-15*x2*(7*y2*y2 - 2*y2*z2)))*invR15;

fy+=srcDenarr[dims*s+47]*0.0625*cl66*(y*(43*x2*x2*x2 - 7*y2*y2*y2 + 6*y2*y2*z2+\
x2*x2*(-225*y2 + 30*z2) +3*x2*(47*y2*y2 - 20*y2*z2)))*invR15;

fy+=srcDenarr[dims*s+48]*0.0625*cl66*(-2*x*(3*x2*x2*x2 - 24*y2*y2*y2 +\
15*y2*y2*z2 +x2*x2*(-66*y2 + 3*z2) +5*x2*(23*y2*y2 - 6*y2*z2)))*invR15;
}


/*if(p>=0)
{
fz =srcDenarr[dims*s]*dZ_reg*invR3;
}*/

if(p>=1)
{
fz+=srcDenarr[dims*s+1]*(3.0*z2-1.0/invR2)*invR5;
fz+=srcDenarr[dims*s+2]*3.0*dX_reg*dZ_reg*invR5;
fz+=srcDenarr[dims*s+3]*3.0*dY_reg*dZ_reg*invR5;
}

if(p>=2)
{
fz+=srcDenarr[dims*s+4]*1.5*dZ_reg*(5.*z2-3.0/invR2)*invR7;
fz+=srcDenarr[dims*s+5]*cl2*dX_reg*(5.0*z2-1.0/invR2)*invR7;
fz+=srcDenarr[dims*s+6]*cl2*dY_reg*(5.0*z2-1.0/invR2)*invR7;
fz+=srcDenarr[dims*s+7]*0.5*cl2*5.0*dZ_reg*pm2c*invR7;
fz+=srcDenarr[dims*s+8]*cl2*5.0*pm1c*pm1s*dZ_reg*invR7;
}

if(p>=3)
{
fz+=srcDenarr[dims*s+9]*(1.5*x2*x2+1.5*y2*y2-12.*y2*z2+4.*z2*z2+\
z2*(3.*y2-12.*z2))*invR9;
fz+=srcDenarr[dims*s+10]*0.25*cl31*(-5.0*dX_reg*dZ_reg*(3*x2+3*y2-4*x2))*invR9;
fz+=srcDenarr[dims*s+11]*0.25*cl31*(-5*dY_reg*dZ_reg*(3*x2+3*y2-4*z2))*invR9;
fz+=srcDenarr[dims*s+12]*0.5*cl32*((-x2 + y2)*(x2 + y2 - 6*z2))*invR9;
fz+=srcDenarr[dims*s+13]*cl32*(-(x*y*(x2 + y2 - 6*z2)))*invR9;
fz+=srcDenarr[dims*s+14]*0.25*cl33*(7*dX_reg*dZ_reg*(x2 - 3*y2))*invR9;
fz+=srcDenarr[dims*s+15]*0.25*cl33*(7*dY_reg*dZ_reg*(3*x2 - y2))*invR9;
}


if(p>=4)
{
fz+=srcDenarr[dims*s+16]*(dZ_reg*(9.375*x2*x2 + 9.375*y2*y2 - 25.*y2*z2 +\
5.*z2*z2 + x2*(18.75*y2 - 25.*z2)))*invR11;
fz+=srcDenarr[dims*s+17]*0.25*cl41*(3*x*(x2*x2 + y2*y2 - 12*y2*z2 +\
8*z2*z2 + 2*x2*(y2 - 6*z2)))*invR11;
fz+=srcDenarr[dims*s+18]*0.25*cl41*(3*y*(x2*x2 + y2*y2 - 12*y2*z2 +\
8*z2*z2 + 2*x2*(y2 - 6*z2)))*invR11;
fz+=srcDenarr[dims*s+19]*0.25*cl42*(-21*(x2 - y2)*z*(x2 + y2 - 2*z2))*invR11;
fz+=srcDenarr[dims*s+20]*0.5*cl42*(-21*x*y*z*(x2 + y2 - 2*z2))*invR11;
fz+=srcDenarr[dims*s+21]*0.25*cl43*dZ_reg*(9*(x2 - 3*y2)*z*x)*invR11;
fz+=srcDenarr[dims*s+22]*0.25*cl43*dZ_reg*(9*(3*x2 - y2)*z*y)*invR11;
fz+=srcDenarr[dims*s+23]*0.125*cl44*(9*(x2*x2 - 6*x2*y2 + y2*y2)*z)*invR11;
fz+=srcDenarr[dims*s+24]*0.5*cl44*(9*x*y*(x2 - y2)*z)*invR11;
}

if(p>=5)
{
fz+=srcDenarr[dims*s+25]*(-1.875*x2*x2*x2 - 1.875*y2*y2*y2 + 33.75*y2*y2*z2 -\
45.*y2*z2*z2 + 6.*z2*z2*z2 +x2*x2*(-5.625*y2 + 33.75*z2) +x2*(-5.625*y2*y2 +\
67.5*y2*z2 -45.*z2*z2))*invR13;

fz+=srcDenarr[dims*s+26]*0.125*cl51*(7*x*z*(5*x2*x2 + 5*y2*y2 - 20*y2*z2 +\
8*z2*z2 + 10*x2*(y2 - 2*z2)))*invR13;

fz+=srcDenarr[dims*s+27]*0.125*cl51*(7*y*z*(5*x2*x2 + 5*y2*y2 - 20*y2*z2 +\
8*z2*z2 + 10*x2*(y2 - 2*z2)))*invR13;

fz+=srcDenarr[dims*s+28]*0.125*cl52*(2*(x2 - y2)*(x2*x2 + y2*y2 -\
16*y2*z2 + 16*z2*z2 +2*x2*(y2 - 8*z2)))*invR13;

fz+=srcDenarr[dims*s+29]*0.25*cl52*(2*x*y*(x2*x2 + y2*y2 - 16*y2*z2 +16*z2*z2 +\
2*x2*(y2 - 8*z2)))*invR13;

fz+=srcDenarr[dims*s+30]*0.0625*cl53*(-18*x*(x2-3*y2)*z*(3*x2+3*y2-8*z2))*invR13;
fz+=srcDenarr[dims*s+31]*0.0625*cl53*(18*y*(-3*x2 + y2)*z*\
(3*x2 + 3*y2 - 8*z2))*invR13;

fz+=srcDenarr[dims*s+32]*0.375*cl54*(-1.0*(x2*x2 - 6*x2*y2 + y2*y2)*\
(x2 + y2 - 10*z2))*invR13;

fz+=srcDenarr[dims*s+33]*1.5*cl54*(-1.0*x*y*(x2 - y2)*(x2 + y2 - 10*z2))*invR13;
fz+=srcDenarr[dims*s+34]*0.375*cl55*(11*x*(x2*x2 - 10*x2*y2 + 5*y2*y2)*z)*invR13;
fz+=srcDenarr[dims*s+35]*0.375*cl55*(-5*x2*x2*x2 + 6*y2*y2*y2 - 5*y2*y2*z2 +\
x2*x2*(80*y2 - 5*z2) +x2*(-85*y2*y2 + 30*y2*z2))*invR13;
}


if(p>=6)
{
fz+=srcDenarr[dims*s+36]*(z*(-15.3125*x2*x2*x2 - 15.3125*y2*y2*y2+\
91.875*y2*y2*z2 - 73.5*y2*z2*z2 +7.*z2*z2*z2 + x2*x2*(-45.9375*y2 + 91.875*z2)+\
x2*(-45.9375*y2*y2 + 183.75*y2*z2 -73.5*z2*z2)))*invR15;

fz+=srcDenarr[dims*s+37]*0.125*cl61*(3*x*y*z*(15*x2*x2 + 15*y2*y2 - 80*y2*z2 +\
48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;

fz+=srcDenarr[dims*s+38]*0.125*cl61*(-1.0*y*(5*x2*x2*x2 + 5*y2*y2*y2 -\
120*y2*y2*z2 +240*y2*z2*z2 - 64*z2*z2*z2 +15*x2*x2*(y2 - 8*z2)+\
15*x2*(y2*y2 - 16*y2*z2 +16*z2*z2)))*invR15;

fz+=srcDenarr[dims*s+39]*0.0625*cl62*(3*(x2 - y2)*z*(15*x2*x2 + 15*y2*y2-\
80*y2*z2 +48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;

fz+=srcDenarr[dims*s+40]*0.125*cl62*(3*x*y*z*(15*x2*x2 + 15*y2*y2 - 80*y2*z2+\
48*z2*z2 + 10*x2*(3*y2 - 8*z2)))*invR15;

fz+=srcDenarr[dims*s+41]*0.0625*cl63*(2*x*(x2 - 3*y2)*(3*x2*x2 + 3*y2*y2 -\
60*y2*z2 +80*z2*z2 + 6*x2*(y2 - 10*z2)))*invR15;

fz+=srcDenarr[dims*s+42]*0.0625*cl63*(-2*y*(-3*x2 + y2)*(3*x2*x2 + 3*y2*y2 -\
60*y2*z2 +80*z2*z2 + 6*x2*(y2 - 10*z2)))*invR15;

fz+=srcDenarr[dims*s+43]*0.09375*cl64*(-22*(x2*x2 - 6*x2*y2 + y2*y2)*z*\
(3*x2 + 3*y2 - 10*z2))*invR15;

fz+=srcDenarr[dims*s+44]*0.375*cl64*(-22*x*y*(x2 - y2)*z*\
(3*x2 + 3*y2 - 10*z2))*invR15;

fz+=srcDenarr[dims*s+45]*0.375*cl65*(-1.0*x*(x2*x2 - 10*x2*y2 + 5*y2*y2)*\
(x2 + y2 - 12*z2))*invR15;

fz+=srcDenarr[dims*s+46]*0.375*cl65*(-1.0*y*(5*x2*x2 - 10*x2*y2 + y2*y2)*\
(x2 + y2 - 12*z2))*invR15;

fz+=srcDenarr[dims*s+47]*0.0625*cl66*(13*(x2*x2*x2 - 15*x2*x2*y2 + 15*x2*y2*y2-\
y2*y2*y2)*z)*invR15;

fz+=srcDenarr[dims*s+48]*0.0625*cl66*(26*x*y*(3*x2*x2-\
10*x2*y2+3*y2*y2)*z)*invR15;

}

pot_m[3*imind[t]+0] += fx*imq[t];
pot_m[3*imind[t]+1] += fy*imq[t];
pot_m[3*imind[t]+2] += fz*imq[t];


pot_m[3*s+0] -= fx*imq[t];
pot_m[3*s+1] -= fy*imq[t];
pot_m[3*s+2] -= fz*imq[t];


}
}
}
//  printf("%f %f %f\n", pot_m[3*(N-1)+0],pot_m[3*(N-1)+1],pot_m[3*(N-1)+2]);

/*multipole-dipole interaction (up to dipole-dipole)*/

//Gan note 10/28/2019: removed. wehave new exact method using hess.
                                                                               
 /*---------New exact method: use hess for dipole-dipole interaction----------*/
int Ntot=N-N_ion, ier,iprec=fmmtol,ifcharge=1,\
	ifdipole=0,ifpot=0,iffld=0, ifhess=1;


/*ultimate version: use hess to calculate the multipole-dipole force*/
double rscale=1.0;
int nlege=p;

 for(i=0;i<N-N_ion; i++)
      for(j=0;j<N-N_ion; j++)
      {
          if(j!=i)
          {
              double center[3];

              center[0]=ox[i];
              center[1]=oy[i];
              center[2]=oz[i];

              mptarget[0]=ox[j];
              mptarget[1]=oy[j];
              mptarget[2]=oz[j];

              complex realpot;

              l3dmpevalhessd_trunc_(&rscale,center,BknmCopy[i][0],&p,mptarget,\
		&realpot,&iffld,pot,&ifhess,fld[0],scarray,wlege,&nlege);

              double fx,fy,fz,dx,dy,dz,pxx,pyy,pzz,pxy,pxz,pyz;
              dx=mpdipstr[j].real*mpdipvec[j][0];
              dy=mpdipstr[j].real*mpdipvec[j][1];
              dz=mpdipstr[j].real*mpdipvec[j][2];
              pxx=fld[0][0].real;
              pyy=fld[0][1].real;
              pzz=fld[0][2].real;
              pxy=fld[1][0].real;
              pxz=fld[1][1].real;
              pyz=fld[1][2].real;

              fx=-1.0*(dx*pxx+dy*pxy+dz*pxz);
              fy=-1.0*(dx*pxy+dy*pyy+dz*pyz);
              fz=-1.0*(dx*pxz+dy*pyz+dz*pzz);

              pot_m[3*(j+N_ion)+0]+=fx;
              pot_m[3*(j+N_ion)+1]+=fy;
              pot_m[3*(j+N_ion)+2]+=fz;

          }
      }


/*use the fmm for coulomb interaction between all point charges 
(dipoles are embedded in the multipoles, so not counted here)*/

Ntot=imcount+N;
iffld=1;
//redefine the Ntot and iffld flags

if(Ntot<FMM_thresh)
    iprec=6; //if the number is not big, do direct sum.

lfmm3dpartself_(&ier,&iprec,&Ntot,source[0],&ifcharge,charge,&ifdipole,dipstr,\
		dipvec[0],&ifpot,pot,&iffld,fld[0]);

iprec=fmmtol;


for(i=0;i<N;i++)
{
    pot_m[3*i+0]+=fld[i+imcount][0].real*q[i];
    pot_m[3*i+1]+=fld[i+imcount][1].real*q[i];
    pot_m[3*i+2]+=fld[i+imcount][2].real*q[i];
}

for(i=0;i<imcount;i++)
{
    pot_m[3*imind[i]+0]+=fld[i][0].real*imq[i];
    pot_m[3*imind[i]+1]+=fld[i][1].real*imq[i];
    pot_m[3*imind[i]+2]+=fld[i][2].real*imq[i];
}


//  tfinish_force = clock();


/*-------------------------------------------------------------*/


for ( i=0 ; i<N ; i++ )
{
	acc_x[i] += pot_m[3*i+0]/ mass[i];
	acc_y[i] += pot_m[3*i+1]/ mass[i];
	acc_z[i] += pot_m[3*i+2]/ mass[i];
}
								}

/*---------if n_im+N is large, use FMM to calculate force----------*/
else
{
int Ntot=imcount+N, ier,iprec=fmmtol,ifcharge=1,ifdipole=1,ifpot=0,iffld=1;
//here iprec defines the accuracy for FMM. (1, 3digit, 2, 6digit, 3,digit.)


lfmm3dpartself_(&ier,&iprec,&Ntot,source[0],&ifcharge,charge,&ifdipole,dipstr,\
		dipvec[0],&ifpot,pot,&iffld,fld[0]);


for(i=imcount;i<imcount+N;i++)
{
	fld[i][0].real*=q[i-imcount];
	fld[i][1].real*=q[i-imcount];
	fld[i][2].real*=q[i-imcount];
}

/*dipole-source interaction*/
for(i=imcount+N_ion;i<imcount+N;i++)
{
	for(j=0;j<imcount+N;j++)
	{
		if(j!=i)
		{
	double testx=source[i][0]-source[j][0];
	double testy=source[i][1]-source[j][1];
	double testz=source[i][2]-source[j][2];
	double testr=sqrt(testx*testx+testy*testy+testz*testz);
	double r3=testr*testr*testr;
	double pdotr=dipstr[i].real*(dipvec[i][0]*testx+\
		dipvec[i][1]*testy+dipvec[i][2]*testz)/testr;
	fld[i][0].real-=charge[j].real*(3*pdotr*testx/testr-\
			dipstr[i].real*dipvec[i][0])/r3;
	fld[i][1].real-=charge[j].real*(3*pdotr*testy/testr-\
			dipstr[i].real*dipvec[i][1])/r3;
	fld[i][2].real-=charge[j].real*(3*pdotr*testz/testr-\
			dipstr[i].real*dipvec[i][2])/r3;
		}
	}

}


/*dipole-multipole interaction (up to dipole)*/

for(i=N_ion;i<N;i++)
{           
//i is the target dipole, j is the source dipole.
    for (j=i+1;j<N;j++)
    {
        double fddx=0.,fddy=0.,fddz=0.;
        int dims=(p+1)*(p+1);

        double dX=srcPosarr[3*i+0]-srcPosarr[3*j+0];
        double dY=srcPosarr[3*i+1]-srcPosarr[3*j+1];
        double dZ=srcPosarr[3*i+2]-srcPosarr[3*j+2];

        double x2=dX*dX;
        double y2=dY*dY;
        double z2=dZ*dZ;
        double invR =(x2+y2+z2);
        invR=1.0/sqrt(invR);
        double invR5=invR*invR*invR*invR*invR;
        double invR7=invR5*invR*invR;

        double innerprodc1, innerprodc2, innerprodc3;
        innerprodc1=srcDenarr[dims*i+1]*srcDenarr[dims*j+1]+srcDenarr[dims*i+2]*\
	srcDenarr[dims*j+2]+srcDenarr[dims*i+3]*srcDenarr[dims*j+3]; 
	//between source and target

        innerprodc2=dZ*srcDenarr[dims*j+1]+dX*srcDenarr[dims*j+2]+\
			dY*srcDenarr[dims*j+3];//between source and dr

        innerprodc3=dZ*srcDenarr[dims*i+1]+dX*srcDenarr[dims*i+2]+\
			dY*srcDenarr[dims*i+3];//between target and dr

        //dp-dp term1
        fddx+=innerprodc1*3*dX*invR5;
        fddy+=innerprodc1*3*dY*invR5;
        fddz+=innerprodc1*3*dZ*invR5;



        //dp-dp term2
        fddx+=-15.0*dX*innerprodc2*innerprodc3*invR7+3.0*srcDenarr[dims*i+2]*\
	innerprodc2*invR5+3.0*srcDenarr[dims*j+2]*innerprodc3*invR5;
        
	fddy+=-15.0*dY*innerprodc2*innerprodc3*invR7+3.0*srcDenarr[dims*i+3]*\
	innerprodc2*invR5+3.0*srcDenarr[dims*j+3]*innerprodc3*invR5;
        
	fddz+=-15.0*dZ*innerprodc2*innerprodc3*invR7+3.0*srcDenarr[dims*i+1]*\
	innerprodc2*invR5+3.0*srcDenarr[dims*j+1]*innerprodc3*invR5;


        fld[i+imcount][0].real += fddx;
        fld[i+imcount][1].real += fddy;
        fld[i+imcount][2].real += fddz;


        fld[j+imcount][0].real -= fddx;
        fld[j+imcount][1].real -= fddy;
        fld[j+imcount][2].real -= fddz;


    }
}


	for(i=0;i<imcount;i++)
	{
		fld[imind[i]+imcount][0].real+=fld[i][0].real*imq[i];
		fld[imind[i]+imcount][1].real+=fld[i][1].real*imq[i];
		fld[imind[i]+imcount][2].real+=fld[i][2].real*imq[i];
	}


	for ( i=0 ; i<N ; i++ )
	{
		acc_x[i] += fld[i+imcount][0].real/ mass[i];
		acc_y[i] += fld[i+imcount][1].real/ mass[i];
		acc_z[i] += fld[i+imcount][2].real/ mass[i];
	}
}
}
/*-----------------------------------------------------------------*/

/*--------------compute the potential energy------------------------*/
if(energy_compute==1)
{
//   tstart_energy = clock();

if(imcount<imnum_thresh)
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
        ssheval_(Bknm[i][0],&p,loc,&fval);

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

            ssheval_(Bknm[i][0],&p,loc,&fval);

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
    }
//    printf("total self energy  (harmonic-col)= %.15f\n",energy);

for(i=0;i<imcount;i++)
    for(ii=0;ii<N_ion;ii++)
    {
        dx=imx[i]-x[ii];
        dy=imy[i]-y[ii];
        dz=imz[i]-z[ii];
        dr=sqrt(dx*dx+dy*dy+dz*dz);
        energy=energy+0.5*q[ii]*imq[i]/dr;
    }


//	printf("total self energy  (image-ion)= %.15f\n",energy);

for(i=0;i<imcount;i++)
    for(ii=N_ion;ii<N;ii++)
    {
        dx=imx[i]-x[ii];
        dy=imy[i]-y[ii];
        dz=imz[i]-z[ii];
        dr=sqrt(dx*dx+dy*dy+dz*dz);
        if(dr>orad[ii-N_ion])
            energy=energy+0.5*q[ii]*imq[i]/dr;
    }

//	printf("total self energy  (image-col)= %.15f\n",energy);

for(i=0;i<N_ion;i++)
    for(ii=i+1;ii<N_ion;ii++)
    {
        loc[0]=x[ii]-x[i];
        loc[1]=y[ii]-y[i];
        loc[2]=z[ii]-z[i];
        r=sqrt(loc[0]*loc[0]+loc[1]*loc[1]+loc[2]*loc[2]);
        energy=energy+q[i]*q[ii]/r;
    }

//printf("total self energy  (ion-ion)= %.15f\n",energy);

for(i=N_ion;i<N;i++)
    for(ii=i+1;ii<N;ii++)
    {
        loc[0]=x[ii]-x[i];
        loc[1]=y[ii]-y[i];
        loc[2]=z[ii]-z[i];
        r=sqrt(loc[0]*loc[0]+loc[1]*loc[1]+loc[2]*loc[2]);
        energy=energy+q[i]*q[ii]/r;
    }

//printf("total self energy  (col-col)= %.15f\n",energy);

for(i=0;i<N_ion;i++)
    for(ii=N_ion;ii<N;ii++)
    {
        loc[0]=x[ii]-x[i];
        loc[1]=y[ii]-y[i];
        loc[2]=z[ii]-z[i];
        r=sqrt(loc[0]*loc[0]+loc[1]*loc[1]+loc[2]*loc[2]);
        energy=energy+q[i]*q[ii]/r; 
	//if sigma treated in the boundary cond, add 0.5
    }

}
else
{
int Ntot=imcount+N, ier,iprec=fmmtol,ifcharge=1,ifdipole=1,ifpot=1,iffld=0;
//here iprec defines the accuracy for FMM. (1, 3digit, 2, 6digit, 3,digit.)


lfmm3dpartself_(&ier,&iprec,&Ntot,source[0],&ifcharge,charge,&ifdipole,dipstr,\
		dipvec[0],&ifpot,pot,&iffld,fld[0]);

energy=0.0;

for(i=0;i<N;i++)
{
	energy+=0.5*q[i]*pot[i+imcount].real;
}

//substract the image contribution from the own sphere 
//(which should not be counted)

for (i=0;i<M;i++)
    for(ii=0;ii<imcount;ii++)
        {
            if(imind[ii]==i+N_ion)
            {
                dx=imx[ii]-ox[i];
                dy=imy[ii]-oy[i];
                dz=imz[ii]-oz[i];
                dr=sqrt(dx*dx+dy*dy+dz*dz);
                energy-=0.5*osigma[i]*imq[ii]/dr;
            }
        }

//FMM only computes the potential upto dipole. now we sum the contrubtions 
//from higher order multipoles with direct sum and truncation (short-ranged)
if(p>1)
{
for(i=0;i<M;i++)
for(ii=0;ii<N_ion;ii++)
{
    loc[0]=x[ii]-ox[i];
    loc[1]=y[ii]-oy[i];
    loc[2]=z[ii]-oz[i];
    r=sqrt(loc[0]*loc[0]+loc[1]*loc[1]+loc[2]*loc[2]);

    if(r<cutoff*orad[i])
    {
    powr=1.0;
    for(k=0;k<p+1;k++)
    {
        powr=powr*r;
        for(j=0;j<2*p+1;j++)
        {

            BknmCopy[i][j][k].real=BknmCopy[i][j][k].real/powr;
            BknmCopy[i][j][k].imag=BknmCopy[i][j][k].imag/powr;

        } //rescale.
    }

    ssheval_(BknmCopy[i][0],&p,loc,&fval);

    energy=energy+0.5*q[ii]*fval.real;


	powr=1.0;
	for(k=0;k<p+1;k++)
	{
	powr=powr*r;
	for(j=0;j<2*p+1;j++)

	{
	    BknmCopy[i][j][k].real=BknmCopy[i][j][k].real*powr;
	    BknmCopy[i][j][k].imag=BknmCopy[i][j][k].imag*powr;
	} //rescale back.
	}

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
        if(r<cutoff*orad[i])
        {

        powr=1.0;
        for(k=0;k<p+1;k++)
        {
            powr=powr*r;
            for(j=0;j<2*p+1;j++)
            {
                BknmCopy[i][j][k].real=BknmCopy[i][j][k].real/powr;
                BknmCopy[i][j][k].imag=BknmCopy[i][j][k].imag/powr;
            } //rescale.
        }

        ssheval_(BknmCopy[i][0],&p,loc,&fval);

        energy=energy+0.5*osigma[ii]*fval.real;
        powr=1.0;
        for(k=0;k<p+1;k++)
        {
            powr=powr*r;
            for(j=0;j<2*p+1;j++)
            {
                BknmCopy[i][j][k].real=BknmCopy[i][j][k].real*powr;
                BknmCopy[i][j][k].imag=BknmCopy[i][j][k].imag*powr;
            } //rescale back.
        }
    }


}
}

}
}

printf("total electrostatic energy= %.15f\n",energy);

printed_ele_energy=energy;
// tfinish_energy = clock();
}
//iprint=1;
if ( iprint == 1 )
	output_force();
    	//tfinish_total = clock();
/*
clock_t tstart_total, tfinish_total;
clock_t tstart_fmm, tfinish_fmm;
clock_t tstart_sht, tfinish_sht;
clock_t tstart_gmres, tfinish_gmres;
clock_t tstart_img_generation, tfinish_img_generation;
clock_t tstart_energy, tfinish_energy;
clock_t tstart_force, tfinish_force;*/

/*double total_time = (double)(tfinish_total-tstart_total) / CLOCKS_PER_SEC;
printf( "total time is %f seconds\n", total_time);

double fmm_time = (double)(tfinish_fmm-tstart_fmm) / CLOCKS_PER_SEC;
printf( "fmm time is %f seconds\n", fmm_time);

double sht_time = (double)(tfinish_sht-tstart_sht) / CLOCKS_PER_SEC;
printf( "sht time is %f seconds\n", sht_time);

double gmres_time = (double)(tfinish_gmres-tstart_gmres) / CLOCKS_PER_SEC;
printf( "gmres time is %f seconds\n", gmres_time);

double img_time = (double)(tfinish_img_generation-tstart_img_generation) / CLOCKS_PER_SEC;
printf( "img time is %f seconds\n", img_time);

double energy_time = (double)(tfinish_energy-tstart_energy) / CLOCKS_PER_SEC;
printf( "energy time is %f seconds\n", energy_time);

double force_time = (double)(tfinish_force-tstart_force) / CLOCKS_PER_SEC;
printf( "force time is %f seconds\n", force_time);*/

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
			mpole[i][ii][j][k].imag=vecx[M*(p+1)*(2*p+1)+\
			i*(2*p+1)*(p+1)+j*(p+1)+k];
			
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

			if(dr<m2ltr*orad[ii]*orad[ii])
			{
				for(j=0;j<2*p+1;j++)
					for(k=0;k<p+1;k++)
					{
						mpole[i][ii][j][k].real=\
						mpole[i][ii][j][k].real*sqrtk[k];
						mpole[i][ii][j][k].imag=\
						mpole[i][ii][j][k].imag*sqrtk[k];
					}
			l3dmplocquadu_(&sc1,x0y0z0,mpole[i][ii][0],&p,&sc2,\
					xnynzn,local[i][ii][0],&p,&ier);

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
			//printf("lpole=%d %d %d %d %f %f\n",i, \
			ii,k,j-p, local[i][ii][j][k].real,local[i][ii][j][k].imag);
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

				b[i*(2*p+1)*(p+1)+j*(p+1)+k]=\
			b[i*(2*p+1)*(p+1)+j*(p+1)+k]+vecx[i*(2*p+1)*(p+1)+\
			j*(p+1)+k]*((k+1)*epsi_s+k*epsi_i[i]);

				b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]=\
			b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]+\
			vecx[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+\
			j*(p+1)+k]*((k+1)*epsi_s+k*epsi_i[i]);
			}
			else
			{
				b[i*(2*p+1)*(p+1)+j*(p+1)+k]=\
			b[i*(2*p+1)*(p+1)+j*(p+1)+k]+vecx[i*(2*p+1)*(p+1)+\
			j*(p+1)+k]*((k+1)*epsi_s/k+epsi_i[i]);

				b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]=\
			b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]+\
			vecx[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+\
			j*(p+1)+k]*((k+1)*epsi_s/k+epsi_i[i]);
			}
		}
		else
		{
			if(k==0)
			{
				b[i*(2*p+1)*(p+1)+j*(p+1)+k]=\
			b[i*(2*p+1)*(p+1)+j*(p+1)+k]+\
			local[ii][i][j][k].real*k*(epsi_i[i]-epsi_s);

				b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]=\
			b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]+\
			local[ii][i][j][k].imag*k*(epsi_i[i]-epsi_s);
			}
			else
			{
				b[i*(2*p+1)*(p+1)+j*(p+1)+k]=\
			b[i*(2*p+1)*(p+1)+j*(p+1)+k]+\
		local[ii][i][j][k].real*(epsi_i[i]-epsi_s)*pow(orad[i],2*k+1);

				b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]=\
			b[M*(p+1)*(2*p+1)+i*(2*p+1)*(p+1)+j*(p+1)+k]+\
		local[ii][i][j][k].imag*(epsi_i[i]-epsi_s)*pow(orad[i],2*k+1);
			}
		}

	}
}

}

void preconr(int n, double *x, double *b)
{ 
int j; 
double div = ( diag != 0.0 ? diag : 1.0);
for (j=0;j<n;j++) 
b[j] = x[j]/div;
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
		ncoef1[l][k]=pow(-1.0,double(l+1))*\
		sqrt(factorial(l+k)/(factorial(l-k)*factorial(2*k)))*gamma;
			//ocoef1[l][k]=0.4*pow(-1.0,double(l+2))*\
	(sqrt(double(factorial(2*k))/(double(factorial(l-k)*factorial(l+k)))))*\
	double(factorial(l))*/double(factorial(k));
	}

	for(l=0;l<p+1;l++)
	for(k=0;k<=l;k++)
	for(j=k+1;j<=l;j++)
	{
	ncoef2[l][k][j]=pow(-1.0,double(l+1))*sqrt(factorial(l+j)/(factorial(l+k)\
*factorial(j-k)))*sqrt(factorial(l+j)/(factorial(l-k)*factorial(j+k)))*gamma;
	//ocoef2[l][k][j]=-0.4*double(factorial(k+l-j+1))*double(2*k+1)\
		/double(factorial(2*k+l-j+1))/double(k+1);
	}


	for(l=0;l<p+1;l++)
	for(k=0;k<=l;k++)
	for(j=k+1;j<=l;j++)
	for(i=1;i<=j-k;i++)
	{
	ncoef3[l][k][j][i]=sqrt(factorial(j+k)/(factorial(i)*factorial(j+k-i)))*\
			sqrt(factorial(j-k)/(factorial(i)*factorial(j-k-i)));
	//	ocoef3[l][k][j][i]=	ncoef3[l][k][j][i];
	//	printf("coefs3=%d %d %d  %d %.7f\n",l,k,j,i,ncoef3[l][k][j][i]);
	}
}


/*start the gemerate image function*/
//Gan note: removed, not practical for MD purpose.
/*end*/

void initialization()
{

	int i,j,k,l,ii,jj,kk;

	powerd= new double[2*p+2];

	newpowerd=new double**[N_col];
	for(i=0;i<N_col;i++)
	{newpowerd[i]=new double*[N_col];}
	for(i=0;i<N_col;i++)
		for(j=0;j<N_col;j++)
		{newpowerd[i][j]=new double[100];}

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
	double lamda=epsi_s/(epsi_s+epsi_i[0]),beta=lamda-1.0,\
	gamma=(epsi_i[0]-epsi_s)/(epsi_s+epsi_i[0]);
	
	double dx,dy,dz,dr;
	int i,j,k,l,ii,jj,kk;
	
	wquad = new double[imm-1];
	xquad = new double[imm-1];

	cgqf (imm-1, 1, 0, 0, -1, 1, xquad, wquad);
	for(i=0;i<imm-1;i++)
	{

		xquad[i]=pow((1.0-xquad[i])/2.0,1.0/lamda);

	}

//given a distribution of spheres, pre-compute the power of the distance between
//spheres so that we do not need to compute that again and again in the GMRES iteration.
	double d;
	double a=1.0, D=2.0;

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

