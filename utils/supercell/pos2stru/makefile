FC  = ifort
SRC = mod_comp.f90 mod_poscar.f90 pos2stru.f90
OBJ = mod_comp.o  mod_poscar.o pos2stru.o

pos2stru:$(OBJ)
	$(FC) -o pos2stru  $(OBJ)
%.o : %.f90
	$(FC) -c $(SRC)


clean:
	rm -rf *.mod
	rm -rf *.o
	rm pos2stru
