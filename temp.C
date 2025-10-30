#include <mpi.h>
#include <iostream>

int main(int argc, char* argv[]){
  MPI_Init(&argc, &argv);
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  int a = 1;
  int b = 0;
  int c = 2;
  int d = 0;
  std::cout << "before: " << myRank << "\t" << a << "\t" << b << "\t" << c << std::endl;


  MPI_Request request[4];
  MPI_Status status[4];

  int partnerrank;
  if (myRank == 0) partnerrank = 1;
  else partnerrank = 0;

  MPI_Irecv(&b, 1, MPI_INT, partnerrank, 0, MPI_COMM_WORLD, &request[0]);
  MPI_Irecv(&d, 1, MPI_INT, partnerrank, 1, MPI_COMM_WORLD, &request[1]);

  MPI_Isend(&a, 1, MPI_INT, partnerrank, 0, MPI_COMM_WORLD, &request[2]);
  MPI_Isend(&c, 1, MPI_INT, partnerrank, 1, MPI_COMM_WORLD, &request[3]);

  MPI_Waitall(4, request, status);

  std::cout << "After: " << myRank << "\t" << a << "\t" << b << "\t" << c << std::endl;

  MPI_Finalize();
}