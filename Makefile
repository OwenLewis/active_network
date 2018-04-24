# ============================================================================
# Name        : Makefile
# Author      : Owen Lewis
# Version     :
# Copyright   : Your copyright notice
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler

#
# Makefile
#
#FC = ifort
#FFLAGS = -r8 -O2 -I/usr/local/include/ -us
#LFLAGS =  -L/usr/local/lib -lfftw

FC = /usr/local/bin/gfortran
FFLAGS = -funderscoring -I/usr/local/include -Isrc -O2 -g -fmessage-length=0 -fdefault-real-8 -ffree-form
LFLAGS = -L/usr/local/lib -lfftw3
OBJS = global_param_mod.o  useful_functions_mod.o io_mod.o stokes_2D_fft_mod.o geometry_mod.o im_bnd_intern_mod.o\
im_bnd_opers_mod.o adhesion_opers_mod.o cortex_intern_mod.o cortex_opers_mod.o fluid_sim_mod.o  spectral_cortex.o

.SUFFIXES: .f90, .f
.SILENT:


%.o : ./src/%.f90
	echo "Compiling $<"
	$(FC) $(FFLAGS) $(XTARG) -c -o $@ $<



./bin/spectral_cortex: $(OBJS)
	echo "Creating Program $@"
	$(FC) -o $@ $(OBJS) $(LFLAGS)


# dependencies
adhesion_opers_mod.o: global_param_mod.o useful_functions_mod.o 
cortex_intern_mod.o: global_param_mod.o
cortex_opers_mod.o: global_param_mod.o useful_functions_mod.o
fluid_sim_mod.o: global_param_mod.o io_mod.o im_bnd_intern_mod.o im_bnd_opers_mod.o cortex_intern_mod.o cortex_opers_mod.o adhesion_opers_mod.o stokes_2D_fft_mod.o geometry_mod.o
geometry_mod.o: global_param_mod.o useful_functions_mod.o
im_bnd_intern_mod.o: global_param_mod.o
im_bnd_opers_mod.o: global_param_mod.o useful_functions_mod.o
io_mod.o: global_param_mod.o
stokes_2D_fft_mod.o: global_param_mod.o
useful_functions_mod.o: global_param_mod.o
spectral_cortex.o: fluid_sim_mod.o global_param_mod.o


clean:
	rm -f *.o *.mod *~ *.il
