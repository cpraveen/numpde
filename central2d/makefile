#
# fortran flags
#
FC      = gfortran
FFLAGS  = ${CFLAGS} -fdefault-real-8

TARGETS = central

#
# program dependencies
#

central = main.o \
			 central.o \
			 savesol.o

#
# programs
#

ALL:  central

central:  ${central}
		${FC} -o central ${central}

##############################################################################
# transformation rules
##############################################################################

.f.o:     ; ${FC} ${FFLAGS} -c $*.f
.c.o:     ; ${CC} ${CFLAGS} -c $*.c

##############################################################################
# clean things up
##############################################################################

clean:	
	${RM} -f *.o $(TARGETS)
