/*
    Common declarations for programs of the NewtonLib package.
    (Part iterative linear solvers)
    
 *  Written by        L. Weimann 
 *  Version           1.0
 *  Revision          May 2006
 *  Latest Change     May 2006
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
 
*/
#define RCODE info->rcode
#define MIN(A,B)  ( A < B ? A : B )
#define MAX(A,B)  ( A > B ? A : B )
#define SIGN(A)   ( A > 0 ? 1 : -1 )

#define SMALL  1.0e-150
#define EPMACH 1.0e-17

typedef enum {None=0, Minimum=1, Verbose=2, Debug=3} PRINT_LEVEL ;
typedef enum {False=0, True=1} LOGICAL ;
typedef enum { CheckOnRestart=0, CheckEachIter=1 } TERM_CHECK ;
typedef enum { Absolute=0, Relative=1 } CONV_CHECK ;

typedef void MATVEC(int, double*, double*);
typedef void PRECON(int, double*, double*); 

struct ITLIN_OPT
{
   double tol, rho;
   int i_max, maxiter;
   TERM_CHECK termcheck;   /* GMRES only */
   CONV_CHECK convcheck;   /* PCG only   */
   LOGICAL rescale;        /* GBIT only  */
   PRINT_LEVEL errorlevel, monitorlevel, datalevel;
   FILE *errorfile, *monitorfile, *datafile,
        *iterfile, *resfile, *miscfile;
   double *scale;
};

struct ITLIN_INFO
{
   double precision, normdx, residuum;
   int iter, rcode, subcode, nomatvec, noprecon, noprecl, noprecr;
};

struct ITLIN_DATA
{
  double *res;
  double tau, t, normdx, residuum;
  enum { GMRES=0, GBIT=1, PCG=2 } codeid;
  enum {Initial=1,Intermediate=2,Solution=3,Final=4} mode;
};

struct ITLIN_IO
{
   FILE *errfile, *monfile, *datfile,
        *iterfile, *resfile, *miscfile;
   PRINT_LEVEL errlevel, monlevel, datlevel;
};

#define ERRORLEVEL   itlin_ioctl->errlevel
#define MONITORLEVEL itlin_ioctl->monlevel
#define DATALEVEL    itlin_ioctl->datlevel
#define ERROR        itlin_ioctl->errfile
#define MONITOR      itlin_ioctl->monfile
#define DATA         itlin_ioctl->datfile
#define FITER        itlin_ioctl->iterfile
#define FRES         itlin_ioctl->resfile
#define FMISC        itlin_ioctl->miscfile

extern void daxpy_(int *n, double *alpha, double *x, int *incx,
                   double *y, int *incy);

/* routines defined in utils.c */
int    zibnum_fwalloc(int size, double **ptr, char vname[]);
int    zibnum_iwalloc(int size, int **ptr, char vname[]);
int    zibnum_pfwalloc(int size, double ***ptr, char vname[]);
double zibnum_scaled_norm2(int n, double *v, double *scale);
double zibnum_scaled_sprod(int n, double *v1, double *v2, double *scale);
double zibnum_norm2(int n, double *v);
void   zibnum_scale(int n, double *v1, double *v2, double *scale);
void   zibnum_descale(int n, double *v1, double *v2, double *scale);
void   itlin_noprecon(int n, double *x, double *z);
void   itlin_dataout(int k, int n, double *x, struct ITLIN_DATA *data);
int itlin_parcheck_and_print(int n, MATVEC *matvec,
                             struct ITLIN_OPT *opt, int itlin_code);
