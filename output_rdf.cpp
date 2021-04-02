#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// record radial distribution function (RDF) in file.
//bin_ii records the ion-ion rdf.
//bin_ic records the ion-colloid rdf.

void output_rdf( )
{
	int i;
	double vol;
	double pi=3.141592653589793238463;
	double fac, r2;
	double rmax;

	rmax=Rshell*0.5;
	vol=4.0*pi*(Rshell-rmax)*(Rshell-rmax)*(Rshell-rmax)/3.0;
	fac=vol/(4.0*pi*bin_size*numsample_rdf);

	for(i=0;i<bin_num;i++)
	{
	  r2=rdf[i]*rdf[i];
	  bin_ii[i]=bin_ii[i]*fac/(r2);
	  bin_ic1[i]=bin_ic1[i]*fac/(r2);
	  bin_ic2[i]=bin_ic2[i]*fac/(r2);
	  bin_ic[i]=bin_ic[i]*fac/(r2);

	  bin_c1c1[i]=bin_c1c1[i]*fac/(r2);
	  bin_c1c2[i]=bin_c1c2[i]*fac/(r2);
	  bin_c2c2[i]=bin_c2c2[i]*fac/(r2);
	  bin_cc[i]=bin_cc[i]*fac/(r2);
	}
	
	for ( i=0 ; i<bin_num ; i++ )
	{
	  fprintf(rdf_file, "%f %f %f %f %f %f %f %f %f\n",\
 rdf[i]*size_scaling,bin_ii[i],bin_ic1[i],bin_ic2[i],\
bin_ic[i],bin_c1c1[i],bin_c1c2[i],bin_c2c2[i],bin_cc[i]);

	}
	
	fclose(rdf_file);
	
}
