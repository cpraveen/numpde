#FC = ifort
#FFLAGS = -O3 -diag-disable vec -r8 -free -Tf

FC = gfortran
FFLAGS = -O3 -fdefault-real-8

TARGETS = advect

ALL: $(TARGETS)

SRC = $(wildcard *.f95)
OBJ = $(patsubst %.f95,%.o,$(SRC))

advect: $(OBJ)
	$(FC) -o advect $(OBJ)

comvar.mod: comvar.f95
	$(FC) -c $(FFLAGS) comvar.f95

%.o: %.f95 comvar.mod
	$(FC) -c $(FFLAGS) $*.f95

clean:
	rm -f $(OBJ) $(TARGETS) *.mod *.plt
