#include "MDpara.h"
void allocate_arrays()
{	
//locations, charges, and radius of the sources
	x = new double[N];
	y = new double[N];
	z = new double[N];
	q = new double[N];
	r = new double[N];


//velocities, and also half-step velocities of the sources (same array)
	vx = new double[N];
	vy = new double[N];
	vz = new double[N];

//accelerations of the sources
	acc_x = new double[N];
	acc_y = new double[N];
	acc_z = new double[N];

//forces of the sources
	Fx = new double[N];
	Fy = new double[N];
	Fz = new double[N];
	
//extra random forces for the langevin thermostat.
	Fran_x = new double[N];
	Fran_y = new double[N];
	Fran_z = new double[N];
	noise = new double[N];

	// mass of the sources
	mass = new double[N];
}
