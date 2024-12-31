#  Makefile for kproj
#

SHELL = /bin/sh

FFTW_INCLUDES = /usr/local/include/ 
#FFTW_INCLUDES = /home/jiaxinchen/bin/fftw-3.3.9/build/include

# suffixes rules to use for .f and .o
# these lines tell how to get to obtain a .o file from a .f file

  .SUFFIXES :
  .SUFFIXES : .o .f

  .f.o:
	$(FF) $(OPT) -c $< $(finclude)  

  .F.o:
	$(FF) $(OPT) -c $< $(finclude)  

#  Intel compiler
#IFORT  = /opt/intel/composerxe
# MKLDIR = $(IFORT)/mkl/lib/intel64
# FF  = $(IFORT)/bin/ifort
 FF  = ifort
# OPT = -O3 -free -xhost -mkl
OPT = -O3 -free -xhost -mkl -CB -traceback -I$(FFTW_INCLUDES)
LOPT = -FR -mkl 
#  gfortran
#FF = gfortran
#OPT = -ffree-form -O3 -march=native
#LOPT =


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#================Fox_xml==================================
finclude =   -I/home/jiaxinchen/work/code/2024_11_14/kproj/src/fox/objs/finclude 
foxlib   =   /home/jiaxinchen/work/code/2024_11_14/kproj/src/fox/objs/lib/libFoX_wkml.a /home/jiaxinchen/work/code/2024_11_14/kproj/src/fox/objs/lib/libFoX_dom.a /home/jiaxinchen/work/code/2024_11_14/kproj/src/fox/objs/lib/libFoX_sax.a /home/jiaxinchen/work/code/2024_11_14/kproj/src/fox/objs/lib/libFoX_wcml.a /home/jiaxinchen/work/code/2024_11_14/kproj/src/fox/objs/lib/libFoX_wxml.a /home/jiaxinchen/work/code/2024_11_14/kproj/src/fox/objs/lib/libFoX_common.a /home/jiaxinchen/work/code/2024_11_14/kproj/src/fox/objs/lib/libFoX_utils.a /home/jiaxinchen/work/code/2024_11_14/kproj/src/fox/objs/lib/libFoX_fsys.a 
#=========================================================

#    subdirectories that hold the program files:



SRCDIRS = FFT other_unfold  qe_related_modules abinit_related_modules


# convert the list of SRCDIRS for VPATH to look in
# ($(EMPTY) and $(SPACE) are needed to insert a space in the substitution
EMPTY =
SPACE = $(EMPTY) $(EMPTY)
VPATH = $(subst $(SPACE),:,$(SRCDIRS))

  comp = mod_comp.f dlexlib.f drdatab.f 

  SRC = mod_input.f mod_metric.f  mod_lcao.f   mod_phonon.f  mod_other.f  qes_types_module.f  mod_read_information_ab.f mod_read_information_qe.f  mod_wave.f  mod_lattice.f gfft.f mod_fft.f out_psi2.f out_pchg.f mod_step_funciton.f mod_kproj_abinit.f mod_kproj_qe.f mod_kproj.f  main.f 
 
  utils = bsplt_k_projected.f

  meig = merge_bndfiles.f

  procar = mod_zwin.f mod_struct.f mod_bs.f mod_procar.f mod_reweight.f procar2kproj.f

  OBJ = ${SRC:.f=.o}
  OBJF = ${SRC_F:.F=.o}

#------------------------------------------------------------------

kproj: $(comp) $(OBJ)
	@ echo " $(FF) -o $@ $(comp) $(OBJ) $(LOPT)"; \
	$(FF) -o $@ $(comp) $(OBJ) $(LOPT)  $(foxlib)  $(finclude)

futils: $(utils)
	@ echo " $(FF) -o $@ $(utils) $(LOPT)"; \
	$(FF) -o $@ $(utils) $(LOPT) 

mergefile: $(meig)
	@ echo " $(FF) -o $@ $(meig) $(LOPT)"; \
	$(FF) -o $@ $(meig) $(OPT) $(LOPT) 
clean:
	rm -f core *.o *.mod
	rm -r fox
	rm Makefile
	cp Makefile.in Makefile
fox  :
	sh fox.sh
# dependences (module files also create .o files)
# Note: the *.f files are 'include'd into the main module files and
# listed here to ensure that if one is changed that the module is
# recompiled

mod_fft.o: mod_comp.o mod_wave.o mod_lattice.o gfft.o   \
              fft.f  fft_allocate.f  fft_cutoffs.f \
              wf_to_fftmesh.f cp_g_to_fftmesh.f

main.o : mod_comp.o mod_input.o mod_wave.o mod_lattice.o mod_fft.o
mod_input.o : mod_comp.o 
dlexlib.o : mod_comp.o
drdatab.o : mod_comp.o
mod_wave.o : mod_comp.o  
mod_lattice.o : mod_comp.o mod_wave.o 
mod_fft.o : mod_comp.o mod_lattice.o mod_wave.o 
out_psi2.o : mod_comp.o mod_input.o mod_wave.o mod_fft.o mod_lattice.o
out_pchg.o : mod_comp.o mod_input.o mod_wave.o mod_fft.o mod_lattice.o
mod_kproj.o : mod_comp.o mod_input.o mod_wave.o mod_fft.o mod_lattice.o
mod_zwin.o : mod_comp.o 
mod_struct.o : mod_comp.o 
mod_bs.o : mod_comp.o
