ifndef EIGEN_DIR
TARGETS = poisson1
else
TARGETS = poisson1 poisson2
endif

CXX = c++
CXXFLAGS = -O3 -Wall
HDR = array2d.h

ifeq ($(debug),yes)
	CXXFLAGS += -DDEBUG
endif

all: $(TARGETS)

array2d.o: array2d.h array2d.cc

.cc.o: $(HDR)
	$(CXX) $(CXXFLAGS) -o $*.o -c $*.cc

poisson1: poisson1.o array2d.o $(HDR)
	$(CXX) -o poisson1 poisson1.o array2d.o

poisson2: poisson2.cc
	$(CXX) $(CXXFLAGS) -o poisson2 poisson2.cc -I$(EIGEN_DIR)/include/eigen3

help:
	@echo "Usage:"
	@echo "   make"
	@echo "   make debug=yes"
	@echo "   make poisson1"
	@echo "   make poisson2  # Needs Eigen"
	@echo "   make clean"

clean:
	rm -rf *.o $(TARGETS)
