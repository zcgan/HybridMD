#include"MDpara.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// record trajectory in file
void record_trajectories(int indicator)
{
	int i;
    
    if (indicator==0)
    {
	for ( i=0 ; i<N ; i++ )
	{
	  if(i<N_ion1)
	    fprintf(trajectory_annealing,\
		 "%d %d %f %f %f\n",i, 1, x[i], y[i], z[i]);

	  else if(i<N_ion)
	    fprintf(trajectory_annealing,\
		 "%d %d %f %f %f\n",i, 2, x[i], y[i], z[i]);

	  else if (i<N_ion+N_col1)
	    fprintf(trajectory_annealing,\
		 "%d %d %f %f %f\n",i, 3, x[i], y[i], z[i]);
	  else
	    fprintf(trajectory_annealing,\
		 "%d %d %f %f %f\n",i, 4, x[i], y[i], z[i]);
	}
    }
    else if (indicator==1)
    {
        for ( i=0 ; i<N ; i++ )
        {
            if(i<N_ion1)
                fprintf(trajectory,\
		 "%d %d %f %f %f\n",i, 1, x[i], y[i], z[i]);

            else if(i<N_ion)
                fprintf(trajectory,\
		 "%d %d %f %f %f\n",i, 2, x[i], y[i], z[i]);

            else if (i<N_ion+N_col1)
                fprintf(trajectory,\
		 "%d %d %f %f %f\n",i, 3, x[i], y[i], z[i]);

            else
                fprintf(trajectory,\
		 "%d %d %f %f %f\n",i, 4, x[i], y[i], z[i]);
        }

    }
    
}
