#  Makefile for supercell
#

SHELL = /bin/sh

# suffixes rules to use for .f and .o
# these lines tell how to get to obtain a .o file from a .f file

  .SUFFIXES :
  .SUFFIXES : .o .f

  .f.o:
	$(FF) $(OPT) -c $<

  .F.o:
	$(FF) $(OPT) -c $<

#  Intel compiler
#IFORT  = /opt/intel/composerxe
# MKLDIR = $(IFORT)/mkl/lib/intel64
# FF  = $(IFORT)/bin/ifort
 FF  = ifort
# OPT = -O3 -free -xhost -mkl
OPT = -O3 -free -xhost -mkl -CB -traceback
 LOPT = -FR

#  gfortran
#FF = gfortran
#OPT = -ffree-form -O3 -march=native
#LOPT =


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#    subdirectories that hold the program files:

SRCDIRS = FFT
# convert the list of SRCDIRS for VPATH to look in
# ($(EMPTY) and $(SPACE) are needed to insert a space in the substitution
EMPTY =
SPACE = $(EMPTY) $(EMPTY)
VPATH = $(subst $(SPACE),:,$(SRCDIRS))

  SRC = mod_comp.f  mod_poscar.f90  supercell.f

  OBJ = ${SRC:.f=.o}
  OBJF = ${SRC_F:.F=.o}

#------------------------------------------------------------------

supercell: $(OBJ)
	@ echo " $(FF) -o $@ $(OBJ) $(LOPT)"; \
	$(FF) -o $@ $(OBJ) $(LOPT)

clean:
	rm -f core *.o *.mod

#==================================================================

# dependences (module files also create .o files)
# Note: the *.f files are 'include'd into the main module files and
# listed here to ensure that if one is changed that the module is
# recompiled

mod_poscar.o : mod_comp.o
supercell.o : mod_comp.o mod_poscar.o
