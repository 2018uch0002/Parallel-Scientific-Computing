# heat1d Makefile
all: montecarlo

#  set name of compilers
CC = mpicc
CXX = mpic++
NVCC = nvcc

CCFLAGS = -O3

# first compile with C++
%.o : %.C
	$(CXX) $(CCFLAGS) -o $@ -c $<

# second compile the cuda file with nvcc compiler 
NVCCFLAGS = -O3 -arch=sm_70
%.o : %.cu
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<

AppLibraries = -L/usr/local/cuda/lib64 -lcuda -lcudart -lcurand -lcudadevrt -lstdc++ 

# List of all libraries
LIBS = $(AppLibraries)

# heat1d equation 
montecarlofiles = simulation.o simulationGPU.o
montecarlo: $(montecarlofiles)
	$(CXX) $(CCFLAGS) -o $@ $(montecarlofiles) $(LIBS)

clean:
	rm *.o
	rm montecarlo