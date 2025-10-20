//
// Debug file demo:
// Open separate debug files on each rank
// Examples:
// mpiexec -n 2 debugFileDemo -debug=1
// mpiexec -n 4 debugFileDemo -debug=1

// define this to indicate were are using MPI (later this will mean we ae using P++ instead of A++)
#define USE_PPP

#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <float.h>
#include <assert.h>
// sleep function is here with Linux
#include <unistd.h>

// define a new type "Real" which is equivalent to a "double"
typedef double Real;

#include <string>
using std::string;
using std::max;

// include commands tp parse command line arguments
#include "parseCommand.h"

#include <ctime>
// ---------------------------------------------
// Return the current wall-clock time in seconds
// ---------------------------------------------
inline double getCPU()
{
#ifdef USE_PPP
return MPI_Wtime(); // use MPI timer
#else
return ( 1.0*std::clock() )/CLOCKS_PER_SEC ;
#endif
}


// ======================================================================================
// Return the max value of a scalar over all processors in a communicator
// /processor: return the result to this processor (-1 equals all processors)
// ======================================================================================
Real
getMaxValue(Real value, int processor = -1 , MPI_Comm comm = MPI_COMM_WORLD )
{
Real maxValue=value;
#ifdef USE_PPP
if( processor==-1 )
MPI_Allreduce(&value, &maxValue, 1, MPI_DOUBLE, MPI_MAX, comm);
else
MPI_Reduce(&value, &maxValue, 1, MPI_DOUBLE, MPI_MAX, processor, comm);
#endif
return maxValue;
}



int main( int argc, char* argv[] )
{
  MPI_Init( &argc, &argv ); // initialize MPI
  int myRank;
  MPI_Comm_rank( MPI_COMM_WORLD, &myRank ); // my process number
  int np;
  MPI_Comm_size( MPI_COMM_WORLD, &np ); // total number of proc’s

  char procName[MPI_MAX_PROCESSOR_NAME];
  int resultlen;
  MPI_Get_processor_name( procName, &resultlen);
  int version, subversion;
  MPI_Get_version( &version, &subversion);

  int pid =getpid(); // process ID

  if( myRank==0 )
    printf("Usage: debugFileDemo -debug=i\n");

  int debug=1; // set to 1 for debug info

  string line;
  bool echo = false; // do not echo in parseCommand
  for( int i=1; i<argc; i++ )
  {
    line=argv[i];
    if( parseCommand( line,"-debug=",debug,echo) ){}
  }

  FILE *debugFile=NULL;
  if( debug>0 )
  {
    // open a debug file on each processor (in the debug folder)
    char debugFileName[80];
    sprintf(debugFileName,"debug/debugFileDemoNp%dProc%d.debug",np,myRank);
    debugFile = fopen(debugFileName,"w");
  }


  // Write header info to both stdout and debugFile
  for( int ifile=0; ifile<=1; ifile++ )
  {
    FILE *file = ifile==0 ? stdout : debugFile;
    if( ( ifile==0 && myRank==0 ) || // write to terminal if myRank==0
        ( ifile==1 && debug>0 ) ) // write to debugFile if debug>0
    {
      fprintf(file,"------------------- DebugFileDemo --------------------- \n");
      fprintf(file," np=%d, myRank=%d\n",np,myRank);
    }
  }

  if( debug>0 )
    fprintf(debugFile," procName=[%s], pid=%d, MPI version %d.%d\n",procName,pid,version,subversion);

  Real cpu0=getCPU();

  // Pretend to do some computations by sleeping for some time

  // sleep( myRank*100000 ); // sleep for some time in seconds
  usleep( myRank*100000 ); // sleep for some time in micro secs, 1e-6(s)

  Real cpu = getCPU()-cpu0;

  if( debug>0 )
  {
    fprintf(debugFile," Time to sleep = %9.2e(s)\n",cpu);
    fflush(debugFile); // flush output to file
  }

  cpu = getMaxValue( cpu );
  if( myRank==0 )
    printf(" Max time to sleep = %9.2e(s)\n",cpu);


  if( debugFile )
    fclose(debugFile);

  // close down MPI
  MPI_Finalize();

  return 0;
}
