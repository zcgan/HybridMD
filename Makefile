#	Generic Makefile - Probably requires GnuMake
#	Erik Luijten
#	$Id: Makefile,v 1.2 2000-10-27 08:36:58-04 luijten Exp luijten $
#	Adapted July 1, 2016

CC	= icc
CXX	= icpc
FC	= ifort

CFLAGS	= -Iheader -ansi -pedantic -Wno-long-long -DRDF
FFLAGS  =

#OPT	= -g
OPT	= -O3 -g

LDFLAGS = -liomp5 -lifcore

OBJDIR	= objects

OBJECTS  = main.o read_para.o MDpara.o allocate.o initialize_rdf.o mass_initialize.o initial_conditions.o output_x_v.o output_force.o energy_momentum_check.o Velocity_Verlet.o Velocity_Verlet_Langevin.o LJ_accelerations_Directsum.o record_rdf.o output_rdf.o Velocity_Rescale.o record_trajectories.o gauss.o LJ_Boundaryforce.o Central_sph_force.o Coulomb_accelerations_Hybrid.o Cell_list.o

LIBOBJECTS = jacobi_rule.o gmres.o utils.o \
l3dtrans.o laprouts3d.o rotviarecur3.o cdjseval3d.o d3mtreeplot.o d3tstrcr.o lfmm3drouts.o triasymq.o triagauc.o triquadflatlib.o l3dtrirouts.o lfmm3dtria.o second.o trilib.o  l3dterms.o lfmm3dpart.o rotproj.o triahquad.o \
dfft.o legeexps.o prini.o prinm.o sshexps.o xrecursion.o yrecursion.o

VPATH = FMM3dlib:GMRES:JacobiGaussQuad:sht

$(OBJDIR)/%.o : %.c
	$(CC) $(CFLAGS) $(OPT) -c -o $@ $<

$(OBJDIR)/%.o : %.cpp
	$(CXX) $(CFLAGS) $(OPT) -c -o $@ $<

$(OBJDIR)/%.o : %.f
	$(FC) $(FFLAGS) $(OPT) -c -o $@ $<

default: $(addprefix $(OBJDIR)/, $(OBJECTS)) $(addprefix $(OBJDIR)/, $(LIBOBJECTS))
	$(CXX) -o MD $(CFLAGS) $(OPT) $(addprefix $(OBJDIR)/, $(OBJECTS)) $(addprefix $(OBJDIR)/, $(LIBOBJECTS)) $(LDFLAGS)

### NO CHANGES BELOW THIS LINE ###

.PHONY: cleanall clean remake

clean:
	@rm -f $(addprefix $(OBJDIR)/, $(OBJECTS))

cleanall:
	@rm -f $(addprefix $(OBJDIR)/, $(OBJECTS)) $(addprefix $(OBJDIR)/, $(LIBOBJECTS)) MD

remake: cleanall default

