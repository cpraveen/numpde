FC = gfortran
FFLAGS = -march=native -O3 -fdefault-real-8 -Wall

TARGETS = euler_fo euler_ho

all: ${TARGETS}

%: %.F90
	${FC} ${FFLAGS} -o $@ $<

clean:
	-rm -f *.o *.mod ${TARGETS}
