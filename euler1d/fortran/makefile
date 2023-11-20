FC = gfortran
FFLAGS = -march=native -O3 -fdefault-real-8 -Wall

TARGETS = euler_fo euler_ho euler_sod euler_msod euler_shuosher euler_toro5

all: ${TARGETS}

euler_ho_2: num_flux.F90

%: %.F90
	${FC} ${FFLAGS} -o $@ $<

euler_ho: euler_ho_1.F90 sod.F90
	${FC} ${FFLAGS} -o euler_ho euler_ho_1.F90

euler_sod: euler_ho_2.F90 sod.F90 num_flux.F90
	${FC} ${FFLAGS} -DSOD -o euler_sod euler_ho_2.F90

euler_msod: euler_ho_2.F90 msod.F90 num_flux.F90
	${FC} ${FFLAGS} -DMSOD -o euler_msod euler_ho_2.F90

euler_shuosher: euler_ho_2.F90 shuosher.F90 num_flux.F90
	${FC} ${FFLAGS} -DSHUOSHER -o euler_shuosher euler_ho_2.F90

euler_toro5: euler_ho_2.F90 toro5.F90 num_flux.F90
	${FC} ${FFLAGS} -DTORO5 -o euler_toro5 euler_ho_2.F90

clean:
	-rm -f *.o *.mod ${TARGETS}
