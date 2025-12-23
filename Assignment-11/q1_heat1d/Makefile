# heat1d Makefile
all: heat1d

#  set name of compilers
CC = gcc
CXX = g++
NVCC = nvcc

CCFLAGS = -O3

# first compile with C++
%.o : %.C
	$(CXX) $(CCFLAGS) -o $@ -c $<

# second compile the cuda file with nvcc compiler 
%.o : %.cu
	$(NVCC) $(CCFLAGS) -o $@ -c $<

AppLibraries = -L/usr/local/cuda/lib64 -lcuda -lcudart

# List of all libraries
LIBS = $(AppLibraries)

# heat1d equation 
heat1dFiles = heat1d.o heat1dSolveGPU.o
heat1d: $(heat1dFiles)
	$(CXX) $(CCFLAGS) -o $@ $(heat1dFiles) $(LIBS)

clean:
	rm *.o
	rm heat1d