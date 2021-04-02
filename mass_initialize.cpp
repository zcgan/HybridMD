#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void mass_initialize( int iprint)
{
	int i;
	double g1, dumb;
/*Gan note: I believe Ziwei added non-uniform mass for ions and colloids*/

	for ( i=0 ; i<N ; i++ )
        {
	  if(i<N_ion)
	    mass[i]=m_ion;
	  else if(i<N_ion+N_col1)
	    mass[i]=m_col1;
	  else
	    mass[i]=m_col2;
	}

	if ( iprint == 1 )
	{
		fprintf( stderr, "Mass matrix\n" );

		for(i=0;i<N;i++)
		{
			fprintf( stderr, " %f", mass[i] );
			if ( i>0 && 5*(i/5) == i ) printf( "\n" );	

		}
	}
	if ( iprint == 1 )
		fprintf( stderr, "\n" );
}
