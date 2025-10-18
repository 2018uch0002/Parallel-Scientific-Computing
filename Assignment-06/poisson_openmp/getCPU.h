#ifndef GET_CPU_H
#define GET_CPU_H

#include <ctime>
// ---------------------------------------------
// Return the current wall-clock time in seconds
// ---------------------------------------------
inline double getCPU()
{
#ifdef USE_PPP
return MPI_Wtime(); // use MPI timer

#elif defined(_OPENMP)
return omp_get_wtime(); // use OpenMP timer
#else
return ( 1.0*std::clock() )/CLOCKS_PER_SEC ;
#endif
}

#endif
