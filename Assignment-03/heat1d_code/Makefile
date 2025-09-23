#
# THIS MAKEFILE USES THE FOLLOWING ENVIRONMENTAL VARIABLE:
#       APlusPlus
#
# The APlusPlus env variable should normally be set when you login.
#
# If not you can set it using:
#
# setenv APlusPlus /home/henshw/software/AppPpp-0.8.3/A++/install
#
# If you type:
# ls $APlusPlus/include
# you should see A++.h as well as other include files.

# First target is made by default when using "make", traditionally named "all" 
all: heat1dImp

# Set names of compilers once in case we need to change them
CC = gcc
CXX = g++

CCFLAGS= -fPIC -O3 -I$(APlusPlus)/include -I../../include

# List of libraries for A++
AppLibraries = -Wl,-rpath,$(APlusPlus)/lib -L$(APlusPlus)/lib -lApp -lApp_static

# List of all libraries
LIBS = $(AppLibraries)

# implicit (generic) rule to compile .C files
# $@ = file name of the target
# $< = name of the first prerequistite
%.o : %.C
	$(CXX) $(CCFLAGS) -o $@ -c $<

# 1D heat equation, implicit time-stepping, A++ arrays
heat1dImpFiles = heat1dImp.o tridiagonal.o
heat1dImp: $(heat1dImpFiles)
	$(CXX) $(CCFLAGS) -o $@ $(heat1dImpFiles) $(LIBS)

clean:
	rm *.o

