#include <stdio.h>
#include <mpi.h> // This is the header file for MPI functions

/*
int main(intargc,char∗argv[])

The main program is the starting point for an MPI program.
This one simply reports its own process number.

Inputs: none
Outputs: Prints "Hello World #" where # is the rank of the process.
*/

int main(int argc, char* argv[])
{
  //First thing is to initialize the MPI interface. Some arguments can
  //be passed through to the MPI interface, so the argument list is
  //sent by sending the argument list
  MPI_Init(&argc, &argv);

  //Determine the rank of the current process
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //Determine the number of processes
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  printf("HelloWorld! I am process %d out of %d\n", rank, size);

  //Must shut down the MPI system when you're done.
  //If you don’t, there can be lots of junk left behind.
  MPI_Finalize();
  return 0;
}