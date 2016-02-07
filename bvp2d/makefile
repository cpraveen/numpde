CXX = c++
CXXFLAGS = -O3
HDR = array2d.h
TARGETS = poisson1

all: $(TARGETS)

.cc.o: $(HDR)
	$(CXX) $(CXXFLAGS) -o $*.o -c $*.cc

poisson1: poisson1.o array2d.o
	$(CXX) -o poisson1 poisson1.o array2d.o

clean:
	rm -rf *.o $(TARGETS)
