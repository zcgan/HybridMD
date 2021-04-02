#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// record radial distribution function (RDF) in file.
//bin_ii records the ion-ion rdf.
//bin_ic records the ion-colloid rdf.
//rmax is set to be Rshell/2.
//we only count particles which are rmax away from the boundary.
void record_rdf( )
{
	int i,j,k;
	double ri,rij;
	double lim2, rmax;
	int ibin, ni, nc1, nc2, nc;
	double *tmp_ii, *tmp_ic1, *tmp_ic2,\
		 *tmp_ic, *tmp_c1c1, *tmp_c1c2, *tmp_c2c2, *tmp_cc;

	rmax=Rshell/2.0;
	lim2=(Rshell-rmax)*(Rshell-rmax);
	ni=0;
	nc1=0;
	nc2=0;
	nc=0;
	
	tmp_ii=new double[bin_num];
	tmp_ic1=new double[bin_num];
	tmp_ic2=new double[bin_num];
	tmp_ic=new double[bin_num];
	tmp_c1c1=new double[bin_num];
	tmp_c1c2=new double[bin_num];
	tmp_c2c2=new double[bin_num];
	tmp_cc=new double[bin_num];

	for (i=0; i<bin_num; i++)
	  {
	    tmp_ii[i]=0.0;
	    tmp_ic1[i]=0.0;
	    tmp_ic2[i]=0.0;
	    tmp_ic[i]=0.0;
	    tmp_c1c1[i]=0.0;
	    tmp_c1c2[i]=0.0;
	    tmp_c2c2[i]=0.0;
	    tmp_cc[i]=0.0;
	  }

	numsample_rdf=numsample_rdf+1.0;

	for ( i=0; i<N_ion; i++)
	{
		ri=x[i]*x[i]+y[i]*y[i]+z[i]*z[i];
	// only count particles which are Rshell/2 away from the boundary.
		if(ri<lim2)
		{
		        ni++;
			for(j=0 ; j<N; j++)
			{
				if(j!=i)
				{
				      rij=sqrt((x[i]-x[j])*(x[i]-x[j])+\
			(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]));
				      if (rij<rmax)
				      {
					  ibin=floor(rij/bin_size);
					  if (j<N_ion)
					    tmp_ii[ibin]=tmp_ii[ibin]+1.0;
					  else if (j<N_ion+N_col1)
					  {
					    tmp_ic1[ibin]=tmp_ic1[ibin]+1.0;
					    tmp_ic[ibin]=tmp_ic[ibin]+1.0;
                                          }
					  else
					  {
					    tmp_ic2[ibin]=tmp_ic2[ibin]+1.0;
					    tmp_ic[ibin]=tmp_ic[ibin]+1.0;
					  }
				      }
				}
			}
		}
	}
	
	for ( i=N_ion; i<N_ion+N_col1; i++)
	{
		ri=x[i]*x[i]+y[i]*y[i]+z[i]*z[i];
		if(ri<lim2)
		{
		        nc1++;
			nc++;
			for(j=N_ion; j<N; j++)
			{
				if(j!=i)
				{
				  rij=sqrt((x[i]-x[j])*(x[i]-x[j])+\
			(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]));
				  if (rij<rmax)
				  {
				    ibin=floor(rij/bin_size);
				    if (j<N_ion+N_col1)
				      {
					tmp_c1c1[ibin]=tmp_c1c1[ibin]+1.0;
					tmp_cc[ibin]=tmp_cc[ibin]+1.0;
				      }
				    else
				      {
					tmp_c1c2[ibin]=tmp_c1c2[ibin]+1.0;
					tmp_cc[ibin]=tmp_cc[ibin]+1.0;
				      }
				  }
				}
			}
		}
	}

	for (i=N_ion+N_col1; i<N; i++)
	{
	        ri=x[i]*x[i]+y[i]*y[i]+z[i]*z[i];
		if (ri<lim2)
		{
		      nc2++;
		      nc++;
		      for (j=N_ion+N_col1; j<N; j++)
		      {
			if (j!=i)
			  {
			    rij=sqrt((x[i]-x[j])*(x[i]-x[j])+\
			(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]));
			    if (rij<rmax)
			      {
				ibin=floor(rij/bin_size);
				tmp_c2c2[ibin]=tmp_c2c2[ibin]+1.0;
				tmp_cc[ibin]=tmp_cc[ibin]+1.0;
			      }
			  }
		      }
		}
	}
	
	for (i=0; i<bin_num; i++)
	  {
	    if (ni>0)
	    {
	      bin_ii[i]=bin_ii[i]+tmp_ii[i]/ni/ni;
	      if (nc1>0)
	      {
		bin_ic1[i]=bin_ic1[i]+tmp_ic1[i]/ni/nc1;
	      }
	      if (nc2>0)
	      {
		bin_ic2[i]=bin_ic2[i]+tmp_ic2[i]/ni/nc2;
	      }
	      if (nc>0)
	      {
		bin_ic[i]=bin_ic[i]+tmp_ic[i]/ni/nc;
	      }
	    }
	    if (nc1>0)
	    {
	      bin_c1c1[i]=bin_c1c1[i]+tmp_c1c1[i]/nc1/nc1;
	      if (nc2>0)
	      {
		bin_c1c2[i]=bin_c1c2[i]+tmp_c1c2[i]/nc1/nc2;
	      }
	    }
	    if (nc2>0)
	      bin_c2c2[i]=bin_c2c2[i]+tmp_c2c2[i]/nc2/nc2;

	    if (nc>0)
	      bin_cc[i]=bin_cc[i]+tmp_cc[i]/nc/nc;

	  }

	delete [] tmp_ii;
	delete [] tmp_ic1;
	delete [] tmp_ic2;
	delete [] tmp_ic;
	delete [] tmp_c1c1;
	delete [] tmp_c1c2;
	delete [] tmp_c2c2;
	delete [] tmp_cc;
}
