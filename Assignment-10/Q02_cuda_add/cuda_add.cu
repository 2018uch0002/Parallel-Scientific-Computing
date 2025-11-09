#include <stdio.h>
#include <math.h>
#include<cuda.h>
#include<cuda_runtime.h>

#include "getCPU.h"
#include "parseCommand.h"


static void HandleError( cudaError_t err, const char *file, int line ) {
  if (err != cudaSuccess) {
    printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
             file, line );
    exit( EXIT_FAILURE );
  }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

__global__
void add(int* a_d, int *b_d, int *c_d, int N)
{
  // get the index correxponds to block-ID and thread-ID
  int tid = threadIdx.x;
  int i = blockIdx.x * blockDim.x + tid;

  if (i < N)
  {
    c_d[i] = a_d[i] + b_d[i]; 
  }

  // syncronize the threads. 
  __syncthreads();
}


int main(int argc, char* argv[])
{
  int n = 1000000;
  for (int i = 0; i < argc; i++){
    if( parseCommand(argv[i], "-n=", n) );
  }

  int *a = new int[n];
  int *b = new int[n];
  int *c_cpu = new int [n];
  int *c = new int[n];

  for (int i =0; i < n; i++)
  {
    a[i] = -i;
    b[i] = i*i;  
  }

  const double cpu1 = getCPU(); // in sec
  for (int i =0; i < n; i++)
  {
    c_cpu[i] = a[i]+ b[i];  
  }
  const double cpu2 = getCPU(); // in sec
  const double cpu_time = (cpu2-cpu1) * 1000; // in ms

// allocate the space on device
  int *a_d, *b_d, *c_d;
  int size_arry_device = n * sizeof(int);
  cudaMalloc((void **)& a_d, size_arry_device);
  cudaMalloc((void **)& b_d, size_arry_device);
  cudaMalloc((void **)& c_d, size_arry_device);

  // copy a and b on device
  cudaMemcpy(a_d, a, size_arry_device, cudaMemcpyHostToDevice);
  cudaMemcpy(b_d, b, size_arry_device, cudaMemcpyHostToDevice);

  // printf("  Nb  \t  Nt  \t  Nb*Nt  \t  n  \t  maxErr  \t  GPU(ms)  \t  speedup\n");
  printf("%-10s %-10s %-12s %-12s %-12s %-12s %-12s\n",
           "Nb", "Nt", "Nb*Nt", "n", "maxErr", "GPU(ms)", "speedup");

  // create the cuda event for time 
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  int Nt = 1;
  for (int itr_Nt=0;  itr_Nt <= 10; itr_Nt++)
  {
    if (itr_Nt > 0) Nt *= 2;
    int Nb = ceil( 1.*n/Nt );

    cudaEventRecord(start, 0);
    add<<<Nb, Nt>>>(a_d, b_d, c_d, n);
    cudaEventRecord(stop, 0);

    cudaDeviceSynchronize();

    cudaMemcpy(c, c_d, size_arry_device, cudaMemcpyDeviceToHost);

    int maxErr = 0, err=0;
    for (int i =0; i < n; i++)
    {
      err = abs(c[i] - c_cpu[i]);
      if (err > maxErr) maxErr = err;
    }

    float gpu_time; // in ms
    cudaEventElapsedTime(&gpu_time, start, stop); 
    const double speedup = cpu_time / static_cast<double>(gpu_time);
    // printf("  %d  \t  %d  \t  %d  \t  %d  \t  %d  \t  %12.4e  \t  %12.4e\n", 
    //        Nb, Nt, Nb*Nt, n, maxErr, gpu_time, speedup);
    printf("%-10d %-10d %-12d %-12d %-12d %-12.4g %-12.4g\n",
            Nb, Nt, Nb*Nt, n, maxErr, gpu_time, speedup);
        
  }

  free(a);
  free(b);
  free(c);
  cudaFree(a_d);
  cudaFree(b_d);
  cudaFree(c_d);
  
  return 0;
}


