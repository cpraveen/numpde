FC = gfortran
FFLAGS = -march=native -O3 -fdefault-real-8 -Wall

TARGETS = fdweno

all: ${TARGETS}

%: %.F90
	${FC} ${FFLAGS} -o $@ $<

clean:
	-rm -f *.o *.mod ${TARGETS}
