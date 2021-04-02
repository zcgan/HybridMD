#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
void initialize_rdf( )
{
//gan note 10/25/2019: I believe part of it is fixed by Ziwei Wang.
	int i;
	double rmax;
	//int bin_num;

	rmax=Rshell*0.5;
	bin_num=floor(rmax/bin_size)+1;
	rdf  = new double[bin_num];
	bin_ii=new double[bin_num];
	bin_ic1=new double[bin_num];
	bin_ic2=new double[bin_num];
	bin_ic=new double[bin_num];

	bin_c1c1=new double[bin_num];
	bin_c1c2=new double[bin_num];
	bin_c2c2=new double[bin_num];
	bin_cc=new double[bin_num];
	
	bin_size=rmax/bin_num;
	
	for(i=0;i<bin_num;i++)
	{   
	  rdf[i]=bin_size*(i+0.5);
	  bin_ii[i]=0.0;
	  bin_ic1[i]=0.0;
	  bin_ic2[i]=0.0;
	  bin_ic[i]=0.0;

	  bin_cc[i]=0.0;
	  bin_c1c1[i]=0.0;
	  bin_c1c2[i]=0.0;
	  bin_c2c2[i]=0.0;
	}
	
	// recording rdf in file
	rdf_file = fopen( "rdf.txt", "w" );
	fprintf(rdf_file,"# r ion-ion ion-col1 ion-col2"
" ion-col col1-col1 col1-col2 col2-col2 col-col\n");

} 
