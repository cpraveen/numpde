CFLAGS = -std=c++11 -O3 -Wall

TARGET = fv

all: $(TARGET)

fv: fv_euler_1d.cc
	$(CXX) $(CFLAGS) -o fv fv_euler_1d.cc

clean:
	rm -f *.o $(TARGET)
