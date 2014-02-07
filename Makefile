CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=main.cpp Input.cpp LDFE_reg.cpp Region.cpp Problem.cpp Quadrature.cpp Solver.cpp Sweep.cpp Mapping.cpp Output.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=TRANSPORT

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) -I/usr/local/hdf5/include -L/usr/local/gfortran/lib -lgfortran -L/usr/local/lib -llapack -lblas -L/usr/local/hdf5/lib -lhdf5 $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *o TRANSPORT