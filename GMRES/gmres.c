/*
        GMRES - RESiduum minimization based iterative linear solver

 *  Written by        L. Weimann 
 *  Purpose           Iterative solution of large linear systems
 *  Category          ???. - Linear systems
 *  Keywords          large linear system, iterative solver
 *  Version           1.0
 *  Revision          June 2006
 *  Latest Change     June 2006
 *  Library           NewtonLib
 *  Code              C, Double Precision
 *  Environment       Standard C environment on PC's,
                      workstations and hosts.
 *  Copyright     (c) Konrad-Zuse-Zentrum fuer
                      Informationstechnik Berlin (ZIB)
                      Takustrasse 7, D-14195 Berlin-Dahlem
                      phone : + 49/30/84185-0
                      fax   : + 49/30/84185-125
 *  Contact           Lutz Weimann
                      ZIB, Division Scientific Computing, 
                           Department Numerical Analysis and Modelling
                      phone : + 49/30/84185-185
                      fax   : + 49/30/84185-107
                      e-mail: weimann@zib.de

*    References:

         P. Deuflhard:
         Newton Methods for Nonlinear Problems. -
         Affine Invariance and Adaptive Algorithms.
         Series Computational Mathematics 35, Springer (2004)

   ---------------------------------------------------------------
 
 * Licence
     You may use or modify this code for your own non commercial
     purposes for an unlimited time. 
     In any case you should not deliver this code without a special 
     permission of ZIB.
     In case you intend to use the code commercially, we oblige you
     to sign an according licence agreement with ZIB.
 
 * Warranty 
     This code has been tested up to a certain level. Defects and
     weaknesses, which may be included in the code, do not establish
     any warranties by ZIB. ZIB does not take over any liabilities
     which may follow from acquisition or application of this code.
 
 * Software status 
     This code is under care of ZIB and belongs to ZIB software class 2.
 
      ------------------------------------------------------------
 
 *    Parameters description
      ======================
 
      The calling interface looks as follows:

      extern void gmres(int n, double *y, MATVEC *matvec,
                        PRECON *preconr, PRECON *preconl, double *b,
                        struct ITLIN_OPT *opt, struct ITLIN_INFO *info)

      The structures used within the parameter list are defined
      as follows:
      ---
      struct ITLIN_OPT
      {
         double tol, rho;
         int i_max, maxiter;
         TERM_CHECK termcheck;   /* GMRES only * /
         CONV_CHECK convcheck;   /* PCG only   * /
         LOGICAL rescale;        /* GBIT only  * /
         PRINT_LEVEL errorlevel, monitorlevel, datalevel;
         FILE *errorfile, *monitorfile, *datafile,
              *iterfile, *resfile, *miscfile;
         double *scale;
      };
      
      where the applicable types used within this structure are defined by
      typedef enum {None=0, Minimum=1, Verbose=2, Debug=3} PRINT_LEVEL;
      typedef enum { CheckOnRestart=0, CheckEachIter=1 } TERM_CHECK ;
      ---
	  struct ITLIN_INFO
	  {
		 double precision, normdx, residuum;
		 int iter, rcode, subcode, nomatvec, noprecon, noprecl, noprecr;
	  };
      ---
      
      A detailed description of the parameters follows: 
      
      int n :
      The number of equations and unknown variables of the linear system.
      
      double *y :
      A pointer to an array of double values of size n.
      The array must contain on input an initial guess of the linear system
      solution, which is used as the start-vector of the iteration.
      On output, the pointed array contains an approximate solution vector y*,
      which fits the preconditioned residuum reduction condition
      ||B_l*r^i||/||B_l*r^0|| <= opt->tol,
      where ||w|| denotes the Euclidian norm of w, and B_l denotes the matrix
      of the left preconditioner.
      
      void *matvec(int n, double *y, double *z);
      A pointer to the matrix times vector multiplication user routine.
      This routine is required - no default routine is supplied.
      The parameters are:
        int     n     input  Number of vector components.
        double *y     input  Vector of unknowns, of size n .
        double *z     output Vector which holds the matrix-vector product A*y.
 
      void *preconr(int n, double *z, double *w);
      A pointer to the right preconditioner user routine, which computes 
      w=B_r*z, where B_r should be an approximation of the inverse of the
      matrix A. If a null pointer is supplied, then a dummy preconditioner
      routine will be used which realizes the preconditioner matrix 
      B_r=identity.
        int     n     input  Number of vector components.
        double *z     input  Preconditioned iterate, of size n .
        double *w     output Vector which holds the matrix-vector product B_r*z,
                      i.e. the original iterate.

      void *preconl(int n, double *z, double *w);
      A pointer to the left preconditioner user routine, which computes 
      w=B_l*z, where B_l should be an approximation of the inverse of the
      matrix A. If a null pointer is supplied, then a dummy preconditioner
      routine will be used which realizes the preconditioner matrix 
      B_l=identity.
        int     n     input  Number of vector components.
        double *z     input  Residual vector, of size n .
        double *w     output Vector which holds the matrix-vector product B_l*z,
                      i.e. the preconditioned residuum.
      
      double *b :
      A pointer to an array of double values of size n.
      The pointed array must hold the right hand side of the linear system 
      to solve.
                                 
      struct ITLIN_OPT *opt:
      A pointer to an options structure. The pointed fields of the structure
      contain input options to GMRES.
      
      opt->tol is of type double and must contain the error threshold
      which the (left preconditioned) residuum norm reduction quotient must fit.
      
      opt->maxiter is of type int and must contain the maximum number of allowed
      iterations. if a zero or negative value is supplied, then the maximum
      iteration count will be set to 100.
      
      opt->i_max is of type int and must contain the maximum number of 
      iterations before a restart occurs. The main portion of used memory
      by GMRES depends on i_max, such that n*i_max elements of double
      storage will be used. If a nonpositive value is supplied, then i_max
      is set to 10.
      
      opt->termcheck is of type TERM_CHECK.
      If set to CheckOnRestart, then the residuum norm reduction quotient
      will be only checked when a restart occurs - such saving additional
      computation effort which is necessary on each intermediate iteration
      step to compute the quantity ||B_l*A*B_r*y^i||: One preconr call,
      one matvec call, and one preconl call.
      If set to CheckEachIter, then additional computations on each iteration
      will be done to compute the (left preconditioned) residuum. This roughly
      doubles the number of calls to the routines matvec, preconr and preconl.
      For additional info, refer to the description of the parameter *y above.
      
      opt->errorlevel is of type PRINT_LEVEL. If it is set to level None,
      then no error message will be printed if an error condition occurs.
      If it is set to level Minimum or any higher level, then error messages
      will be printed, if appropriate.
      
      opt->monitorlevel is of type PRINT_LEVEL. If it is set to level None,
      then no monitor output will be printed.
      If it is set to level Minimum, a few infomation will be printed.
      If set to level Verbose, then some infomation about each iteration
      step, fitting into a single line, will be printed. This only applies,
      whenopt->termcheck=CheckEachIter is set, otherwise information will
      be only printed out when a restart occurs. The higher level Debug
      is reserved for future additional information output.
      
      opt->datalevel is of type PRINT_LEVEL. If it is set to level None,
      then no data output will be printed.
      If it is set to level Minimum, then the values of the initial iteration
      vector x and the final vector x will be printed.
      If set to level Verbose, then the iteration vector x will be printed for
      each step. The higher level Debug is reserved for future additional
      information output.
      
      opt->errorfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer, stderr, stdout,
      or to another file pointer which has been initialized by a fopen call.
      If it is set to NULL, opt->errorfile will be set to stdout. The error 
      messages will be printed to opt->errorfile.
      
      opt->monitorfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer, stderr, stdout,
      or to another file pointer which has been initialized by a fopen call.
      If it is set to NULL, opt->monitorfile will be set to stdout. The monitor 
      output will be printed to opt->monitorfile.
      
      opt->datafile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer, stderr, stdout,
      or to another file pointer which has been initialized by a fopen call.
      If it is set to NULL, a file named "gmres.data" will be opened by a
      fopen call and opt->datafile will be set to the filepointer which the
      fopen returns. The data output will be printed to opt->datafile.
      
      opt->iterfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number and
      the iteration vector will be written out to the associated file, for
      each iteration step. If opt->iterfile is set to NULL, no such 
      data will be written out.
      
      opt->resfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number and
      the residuum vector will be written out to the associated file, for
      each iteration step. If opt->resfile is set to NULL, no such 
      data will be written out.
      
      opt->miscfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number, an
      identification number of the calling code (0 for GMRES), the norm of the 
      (left preconditioned) residuum, followed by three zero floating point
      values as placeholders, will be written out, for each iteration step,
      if opt->termcheck=CheckEachIter is set. Otherwise, if 
      opt->termcheck=CheckOnRestart is set, data will be only written out
      when a restart occurs. 
      If opt->miscfile is set to NULL, no such data will be written out.
     
      Note: The output to the files opt->iterfile, opt->resfile and
            opt->miscfile is written as a single for each iteration step.
            Such, the data in these files are suitable as input to the
            graphics utility GNUPLOT.
      
      struct ITLIN_INFO *info:
      A pointer to an info structure. The pointed fields of the structure
      are set output info of GMRES.
      
      info->precision is of type double and is set to the achieved residual
      reduction ||r^i||/||r^0||.
      
      info->iter is set to number of iteration steps done.
      
      info->nomatvec is set to the number of done calls to the matrix times
      vector multiplication user routine matvec. 
      
      info->noprecr is set to the number of done calls to the right 
      preconditioner user routine preconr or the dummy preconditioner routine,
      if the user didn't supply a right preconditioner routine.
      
      info->noprecl is set to the number of done calls to the left 
      preconditioner user routine preconl or the dummy preconditioner routine,
      if the user didn't supply a left preconditioner routine.
      
      info->rcode is set to the return-code of GMRES. A return-code 0
      means that GMRES has terminated sucessfully. For the meaning of a
      nonzero return-code, see the error messages list below.


      The following error conditions may occur: (returned via info->rcode)
      --------------------------------------------------------------------
      
      -999 routine zibnum_fwalloc failed to allocate double memory via malloc.
      -997 routine zibnum_pfwalloc failed to allocate double pointer memory 
           via malloc.
      -995 Internal i/o control block could not be allocated via malloc.
      -994 Internally used data structure could not be allocated via malloc.
      -989 Default data-output file could not be opened via fopen call.
       -99 NULL pointer supplied for matvec - the matrix times vector routine
           must be defined!
         2 Maximum number of iterations (as set by opt->maxiter) exceeded.
        20 Nonpositive input for dimensional parameter n.
        21 Nonpositive value for opt->tol supplied.
         
      Summary of changes:
      -------------------
      
      Version  Date        Changes
      0.1      2006/06/13  Initial Prerelease.
      
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "itlin.h"

double gmres_norm2(int n, double *v);
int gmres_qrfact(int n, int lda, double *a, double *q, int ijob);
void gmres_qrsolv(int n, int lda, double *a, double *q, double *b);

#define MAXITER_DEFAULT 100

extern struct ITLIN_IO *itlin_ioctl;

extern void gmres(int n, double *y, MATVEC *matvec,
                  PRECON *preconr, PRECON *preconl, double *b,
                  struct ITLIN_OPT *opt, struct ITLIN_INFO *info)
{ 
   int i, j, l, m, m1, m11, l1, iter=0, fail,ijob,
       nomatvec=0, nopreconl=0,
       nopreconr=0, i_max, max_iter = opt->maxiter;
   double **V, *r, *y0, *w, *vi, *vip1, *vl, *Hi, *h, *z, *q;
   double s, normvip1, beta, beta0, beta1, eta, tol=opt->tol;
   LOGICAL stop_iter, io_allocated=False;
   TERM_CHECK termcheck = opt->termcheck;
   struct ITLIN_DATA *data=malloc(sizeof(struct ITLIN_DATA));

   if (!itlin_ioctl) itlin_ioctl=malloc(sizeof(struct ITLIN_IO));
   if (!itlin_ioctl) 
     { fprintf(stderr,"\n could not allocate output controlblock\n");
       RCODE=-995; return; }
   else
     io_allocated = True;
   if (!data)
     { fprintf(stderr,"\n could not allocate struct data\n");
       RCODE=-994; return; };
   data->codeid    = GMRES;
   data->normdx    = 0.0;
   data->tau       = 0.0;
   data->t         = 0.0;
   data->mode      = Initial;
   ERRORLEVEL   = opt->errorlevel;
   MONITORLEVEL = opt->monitorlevel;
   DATALEVEL    = opt->datalevel;
   ERROR    = opt->errorfile;
   MONITOR  = opt->monitorfile;
   DATA     = opt->datafile;
   FITER    = opt->iterfile;
   FRES     = opt->resfile;
   FMISC    = opt->miscfile;
   if ( !ERROR && ERRORLEVEL>0 )     ERROR   = stdout;
   if ( !MONITOR && MONITORLEVEL>0 ) MONITOR = stdout;
   if ( !DATA && DATALEVEL>0 )
     { DATA=fopen("gmres.data","w");
       if (!DATA && ERRORLEVEL>0)
         { fprintf(ERROR,"\n fopen of file gmres.data failed\n");
           RCODE=-989; return;
         };
     };
   opt->errorfile   = ERROR;
   opt->monitorfile = MONITOR;
   opt->datafile    = DATA;
   if ( MONITORLEVEL > 0 ) fprintf(MONITOR,"\n GMRES - Version 0.1\n");
   if ( max_iter <= 0 ) max_iter = MAXITER_DEFAULT; 
   RCODE = itlin_parcheck_and_print(n,matvec,opt,0);
   if ( RCODE !=0 ) 
     { if (io_allocated) {free(itlin_ioctl); itlin_ioctl=NULL;};
       if (data) free(data);
       return;
     };
   i_max = opt->i_max;
   if (!preconr) preconr = &itlin_noprecon;
   if (!preconl) preconl = &itlin_noprecon;
   RCODE = zibnum_pfwalloc(i_max+2,&V,"V");       if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(i_max+1,&h,"h");        if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(i_max+1,&z,"z");        if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc((i_max+1)*(i_max+1),&Hi,"Hi");
   if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(n,&r,"r");              if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(n,&y0,"y0");            if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(n,&w,"w");              if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(2*(i_max+1),&q,"q");    if ( RCODE !=0 ) return;
   for (i=0;i<=i_max+1;i++)
     { RCODE = zibnum_fwalloc(n,&V[i],"V[i]");
       if ( RCODE !=0 ) return; 
     };
   /* Initialization */
   if ( MONITORLEVEL>1 )
     fprintf(MONITOR,"\n\n   iter  i  norm(res)        eta\n\n");
   data->res = r;  
   restart:
   if ( iter > 0 && MONITORLEVEL > 1 && termcheck==CheckEachIter )
      fprintf(MONITOR," > RESTART\n");
   for (j=0;j<n;j++) y0[j]=y[j];
   preconr(n,y,V[0]);  nopreconr++;
   matvec(n,V[0],w);   nomatvec++;
   preconl(n,w,r);     nopreconl++;
   for (j=0;j<n;j++) r[j] = b[j]-r[j];
   beta = gmres_norm2(n,r);
   if ( iter==0 )
     { beta0 = beta; eta = ( beta0==0.0 ? 0.0 : 1.0 ); }
   else
     eta = beta; /*gan: use absolute value to avoid error when getting beta0=0*/
   beta1 = beta;
   vi = V[0];
   for (j=0;j<n;j++) vi[j] = r[j]/beta;
   
   for (i=1; i<=i_max && eta>tol; i++)
     { 
       if ( i==1 || termcheck==CheckEachIter )
         {
          data->residuum = beta1;  itlin_dataout(iter,n,y,data);
          data->mode      = Intermediate;
          if ( MONITORLEVEL>1 )
            fprintf(MONITOR," %6i %2i %10.3e %10.3e\n",iter,i,beta1,eta);
         };

       /* I. Orthogonalization  */
       vi = V[i-1];
       preconr(n,vi,V[i]);  nopreconr++;
       matvec(n,V[i],w);    nomatvec++;
       preconl(n,w,r);      nopreconl++;
       for (l=0;l<i;l++)
         { vl = V[l];  s=0.0;
           for (j=0;j<n;j++) s += vl[j]*r[j];   h[l] = s;
         };
       vip1 = V[i];
       for (j=0;j<n;j++) 
         { s = 0.0;  for (l=0;l<i;l++)  s += V[l][j]*h[l];
           vip1[j] = r[j]-s;
         };

       /* II. Normalization  */  
       normvip1 = gmres_norm2(n,vip1);
       for (j=0;j<n;j++) vip1[j] /= normvip1;

       /* III. Update  */ 
       if ( i>1 )
         { for (l=0;l<i;l++)   Hi[(i_max+1)*(i-1)+l] = h[l];
           for (l=0;l<i-1;l++) Hi[(i_max+1)*l+i] = 0.0;
           Hi[(i_max+1)*(i-1)+i] = normvip1;
         }
       else
         {
           Hi[0] = h[0];   Hi[1] = normvip1;
         };

       /* IV. Least squares problem for z*/ 
       ijob = ( i==1 ? 1 : 2);
       fail = gmres_qrfact(i,i_max+1,Hi,q,ijob);
       if ( fail != 0 ) 
         { RCODE = 1; fprintf(ERROR,"\n gmres_qrfact failed with code=%i\n",fail);
           goto errorexit;
         };
       if ( termcheck==CheckEachIter || i==i_max )
         {
           for (l=0;l<i+1;l++) z[l]=0.0;  z[0]=beta;
           gmres_qrsolv(i,i_max+1,Hi,q,z);
    
           /* V. Approximate solution*/ 
           for (j=0;j<n;j++) 
             { s=0.0;
               for (l=0;l<i;l++) s+= V[l][j]*z[l];
               y[j] = y0[j]+s;
             };
           if ( termcheck==CheckEachIter )
             {
               preconr(n,y,V[i+1]);  nopreconr++;
               matvec(n,V[i+1],w);   nomatvec++;
               preconl(n,w,r);       nopreconl++;
               for (j=0;j<n;j++) r[j] = b[j]-r[j];
               beta1 = gmres_norm2(n,r);
               eta = beta1; /*gan: use absolute value to avoid error when getting beta0=0*/
             };
         };
       iter++;
    };
    
   if ( termcheck==CheckEachIter && MONITORLEVEL>1 ) 
     fprintf(MONITOR," %6i %2i %10.3e %10.3e\n",iter,i,beta1,eta);
   if ( eta > tol && iter < max_iter ) goto restart;
   else                                RCODE = ( eta <= tol ? 0 : 2 );
   if ( MONITORLEVEL>1 && termcheck==CheckOnRestart ) 
     fprintf(MONITOR," %6i %2i %10.3e %10.3e\n",iter,i,beta1,eta);
   data->mode = ( RCODE == 0 ? Solution : Final );
   data->residuum = beta1;  itlin_dataout(iter,n,y,data);
     
   errorexit:
   if ( ERRORLEVEL > 0 && RCODE != 0 )
     {
       switch ( RCODE )
        {
         case     2:
           fprintf(ERROR,"\n Error - Maximum allowed number of iterations exceeded\n");
           break;
         default   :
           fprintf(ERROR,"\n Error, code=%i\n",RCODE);
        };
     };
     for (i=0;i<=i_max+1;i++) free(V[i]);
    free(V);  free(h);  free(z);   free(Hi);   free(r);  free(y0);
/* next line added by Zecheng Gan, July 2, 2016 */
      free(w); free(q); free(data);
   
   info->precision = eta;
   info->iter      = iter;
   info->nomatvec  = nomatvec;
   info->noprecl   = nopreconl;
   info->noprecr   = nopreconr;
}

double gmres_norm2(int n, double *v)
{  int i;
   double rval = 0.0;
   for (i=0;i<n;i++) rval += v[i]*v[i];
   return sqrt( rval );
}

int gmres_qrfact(int n, int lda, double *a, double *q, int ijob)
/*  
    Translation of Fortran routine DHEQR from file DGMRES of
    library=slatec(slap) package.
*/
{ int i, j, j1, k,info=0;
  double t, t1, t2, c, s;
  if ( ijob <= 1 )
    { for (k=0;k<n;k++)
        { for (j=0;j<k-1;j++)
            { i=2*j;   c=q[i];    s=q[i+1];
              j1=j+k*lda;  t1=a[j1];  t2=a[j1+1];
              a[j1]=c*t1-s*t2;  a[j1+1]=s*t1+c*t2;
            };
          i=2*k;  j1=k+lda*k;  t1=a[j1];  t2=a[j1+1];
          if      ( t2==0.0 ) { c=1.0;  s=0.0; }
          else if ( fabs(t2)>=fabs(t1) ) 
            { t=t1/t2;  s=-1.0/sqrt(1.0+t*t); c=-s*t; }
          else
            { t=t2/t1;  c=1.0/sqrt(1.0+t*t);  s=-c*t; };
          q[i]=c;   q[i+1]=s;  a[j1]=c*t1-s*t2;
          if( a[j1]==0.0 ) info = k;
        };
    }
  else
    { for (k=0;k<n-1;k++)
        {
         i=2*k;  j1=k+lda*(n-1);
         t1=a[j1];    t2=a[j1+1];
         c=q[i];     s=q[i+1];
         a[j1]=c*t1-s*t2;  a[j1+1]=s*t1 + c*t2;
        };
      info=0;
      j1 = n-1+lda*(n-1);  t1 = a[j1];  t2 = a[j1+1];
      if      ( t2==0.0 ) { c=1.0;  s=0.0; }
      else if ( fabs(t2)>=fabs(t1) ) 
        { t=t1/t2;  s=-1.0/sqrt(1.0+t*t); c=-s*t; }
      else
        { t=t2/t1;  c=1.0/sqrt(1.0+t*t);  s=-c*t; };
      i=2*n-2;  q[i]=c;   q[i+1]=s;   a[j1]=c*t1-s*t2;
      if( a[j1]==0.0 )  info = n-1;
    };
  return info;
}

void gmres_qrsolv(int n, int lda, double *a, double *q, double *b)
/*  
    Translation of Fortran routine DHELS from file DGMRES of
    library=slatec(slap) package.
*/
{ int i, k, kb;
  double t, t1, t2, c, s;
  for (k=0;k<n;k++)
    { i=2*k;     c=q[i];   s=q[i+1];
      t1=b[k];   t2=b[k+1];
      b[k]=c*t1-s*t2;  b[k+1]=s*t1+c*t2;
    };
  for (kb=0;kb<n;kb++)
    {
      k=n-1-kb;   b[k] /= a[k+lda*k];
      t = -b[k];
      for (i=0;i<k;i++)  b[i] += t*a[i+k*lda];
    };
  return;
}
