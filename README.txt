This is the MD source code for simulation of dielectric nanoparticles with Hybrid method acceleration

1) First set the input parameters in the para.txt file. (We have prepared two parameter files in the Test folder)

2) type "make -f Makefile.icc" to compile. (there is another makefile.gnu, which also works as an alternative)

3) type “make -f Makefile.omp” to run it parallel using openmp. (but before you run the code, first specify the number of threads by using type “export OMP_NUM_THREADS=N”, where N is the number of threads you want to use.)

4) ./MD para_file_name to run the simulation. (If you do not specify the input file name, it will automatically search for para.txt)

5) output files include the trajectory & radial distribution function.** #


********************************************************************************************************
**Files in the package**

1)/FMM3dlib/: the fast multipole subrutines.

2)/GMRES/: the gmres subrutines for solving linear system

3)/header/: header file for gmres

4)/JacobiGaussQuad/: subrutines for different types of Gauss-Jacobi quadratures

5)/sht/: spherical harmonic transform functions

6)allocate.cpp: memory allocation for arrays

7)Cell_list.cpp: compute the LJ interaction using cell-list.

8)Coulomb_accelerations_Directsum.cpp: direct pair-wise coulomb summation for homogeneous systems

9)Coulomb_accelerations_FMM.cpp: FMM accelerated pair-wise coulomb summation for homogeneous systems

10)Coulomb_accelerations_Hybrid.cpp: Hybrid method for solving the force for inhomogeneous systems

11)gauss.cpp: generating random gauss distribution.

12)initial_conditions.cpp:initial conditions for the MD simulation, such as initial config, initial velocity, etc.

13)initialize_rdf.cpp: initialize the RDF recording.

14)LJ_accelerations_Directsum.cpp: LJ pair-wise direct summation.

15)LJ_Boundaryforce.cpp: restricted boundary LJ force

16)mass_initialize.cpp: initialize the masses for all the particles

17)output_force.cpp: output the forces acting on all the particles

18)output_rdf.cpp: output the RDF data.

19)output_x_v.cpp: output the location and velocities of all the particles

20)read_para.cpp: read the input parameters from the "para.txt" file

21)record_rdf.cpp: record the RDF statistics

22)record_trajectories.cpp: record the trajectories for all the particles

23)Velocity_Rescale.cpp: rescale the velocity for the system to calm down to equilibrium

24)Velocity_Verlet.cpp: the scheme for integrate the Newton's equation

25)Velocity_Verlet_Langevin.cpp: the integration scheme after adding the Langevin thermostats.

26) main.cpp: the MD main subrutine.


********************************************************************************************************
**Input parameters**
N_Particles: the total number of particles (colloids plus ions)
N_Types: type of particles
N_Colloids: total number of colloids
N_Ions: number of ions
N_Colloid1: we allow for two types of colloids with different charge, this is the number of first type colloids
N_Colloid2: the number of second type colloids
Colloid_Charge1: the charge of the first type colloid
Colloid_Charge2: the charge of the second type colloid
Ion_Charge: valence of ions
Colloid_Radius: radius of colloids
Ion_Radius: radius of ions
clj/2: the clj/2 is half the constant c defined in the LJ potential, see Eq. (45) in Journal of Computational Physics 291 (2015) 317–333.
Colloid_Dielectric: the dielectric inside the colloids
Solvent_Dielectric: the dielectric in the solvent
Colloid_Mass1: the mass of the first type colloid
Colloid_Mass2: the mass of the second type colloid
Ion_Mass: the mass of ions
Inner_Shell_Charge: please set it to be 0 (not used)
Inner_Shell_Radius: please set it to be 0 (not used)
Inner_Shell_Dielectric: please set it to be 1.0 (not used)
Outer_Shell_Radius: the radius of the outside boundary (the system is a spherical shell)
Temperature: the temperature of the system in Langevin thermostats
Initial_Timestep: the initial timestep in Velocity-verlet algorithm
Final_Timestep: the final timestep after the system get into equilibrium
Timestep_IncreaseInterval: after each such interval, we double the time step, until final timestep is reached.
Velocity_RescaleInterval: the cycle number for velocity rescale in the beginning of the MD simulation.
Langevin_gamma: the dampling parameter in langevin thermostats
Langevin_EquilibrationSteps: after reaching equilibrium, we introduce the thermostats, and first run this steps without any sampling
N_lowtemp_for_annealing: the number of steps for low(desired) temperature while annealing
N_hightemp_for_annealing: the number of steps for high temperature while annealing
N_annealing_cycle: the number of high->low temperature cycles for annealing
Temprature_annealing: the higher temperature for annealing
Production_TimeSteps: the number of steps for produce the data
Sample_Interval: sample the data for each sample interval
RDF_Binsize: the bin size for RDF sampling
precision_digits: the precision for hybrid method
read_ConfigFile: if is set to be 0, then randomly generate the initial config. If is set to 1, then automatiaclly use the file name in the next line.

Excluded-volume interactions: see Eq. (45) in Journal of Computational Physics 291 (2015) 317–333. This is identical to the interaction
in Kipton Barros and Erik Luijten. Phys. Rev. Lett. 113, 017801 (2014).
Specifically, for the case of ions and colloids, note that ion-ion
interactions are STLJ (i.e., Delta_ij=0). Epsilon_lj is set to be equal to KbT.

Straightforward extensions:
Allow colloids have different size and permittivity; add additional ion species.

********************************************************************************************************
**Output**
1) trajectory file with name "dump.lammpstrj"

2) Rdf file with name "rdf.txt"

3) energy in each step before equilibrium named "energy_equilibration"

4) energy in each step during production named "energy_production"

5) the final configuration named "finalconfig.txt"

6) trajectory file for the annealing process with name “annealing_dump.lammpstrj”
