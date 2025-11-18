#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include "A++.h"

typedef double Real;
typedef doubleSerialArray RealArray;
typedef intSerialArray IntegerArray;

#include <float.h>
#include <limits.h>
#define REAL_EPSILON DBL_EPSILON
#define REAL_MIN DBL_MIN

// ===========================================================
// =========neccesary variable and macros ====================

constexpr int dirichlet = 1;
constexpr int numberOfDimensions = 2;

#define TRIG_DD 1
#define POLY_DD 2

#ifndef SOLUTION
// #define SOLUTION TRIG_DD
#define SOLUTION POLY_DD
#endif

constexpr Real kappa=.1;

#if SOLUTION == TRIG_DD

constexpr Real pi = M_PI; // value of PI
constexpr Real kx=2., ky=3;
constexpr Real kxp=kx*pi;
constexpr Real kyp=ky*pi;

#define UTRUE(x,y,t) sin( kxp*(x) )*sin( kyp*(y) )*exp(-kappa*(kxp*kxp+kyp*kyp)*(t) )
#define FORCE(x,y,t) ( 0. )

#elif SOLUTION == POLY_DD
// polynomial manufactured solution
constexpr Real c0=.2, c1=.1, c2=.3;
constexpr Real b0=1., b1=.5, b2=.25;
constexpr Real a0=1., a1=.3, a2=0.;
#define UTRUE(x,y,t)   (b0 + (x)*( b1 + (x)*b2  ))*(c0 + (y)*( c1 + (y)*c2  ))*( a0 + (t)*( a1 + (t)*a2 ) )
#define UTRUET(x,y,t)  (b0 + (x)*( b1 +(x)*b2   ))*(c0 + (y)*( c1 + (y)*c2  ))*(         a1 + 2.*(t)*a2 )
#define UTRUEXX(x,y,t)                   ( 2.*b2 )*(c0 + (y)*( c1 + (y)*c2  ))*( a0 + (t)*( a1 + (t)*a2 ) )
#define UTRUEYY(x,y,t) (b0 + (x)*( b1 + (x)*b2  ))*(     2.*c2               )*( a0 + (t)*( a1 + (t)*a2 ) )

#define FORCE(x,y,t) ( UTRUET(x,y,t) - kappa*( UTRUEXX(x,y,t)+ UTRUEYY(x,y,t) ) )

#endif

// ===========================================================


extern "C"
{
// for CUDA
void heat2dSolveGPU(const int& debug, const int& nx, const int& ny, const int& numSteps,
  const int& numGhost, const int& n1a, const int& n1b, const int& nd1, const int& nd1a,
  const int& nd1b, const int& n2a, const int& n2b, const int& nd2, const int& nd2a,
  const int& nd2b, const int& saveMatlab, const double& tFinal, const double& rx, 
  const double& ry, const double& dt, const double& cpuTimeStep, double& maxErr,
  IntegerArray& gridIndexRange, IntegerArray& dimension, IntegerArray& boundaryCondition, 
  Range& Rx, Range& Ry, RealArray& x, const Real *dx, RealArray *ua, 
  const std::string& matlabFileNameGPU);
}

// initalize the ua[0] for the first step
__global__
void apply_inital_condition(double *x_d, double *uc_d, double t, int* nd1_d, 
      int* nd1a_d, int* nd1b_d, int* nd2_d, int* nd2a_d, int* nd2b_d)
{
  int& nd1 = *nd1_d;
  int& nd1a = *nd1a_d;
  int& nd1b = *nd1b_d;
  int& nd2 = *nd2_d; 
  int& nd2a = *nd2a_d;
  int& nd2b = *nd2b_d;
  
  #define INDEX1(tx, bx) ((tx) + blockDim.x * (bx))
  int tid = INDEX1(threadIdx.x, blockIdx.x);
    // find i1 and i2
  int i2 = tid / nd1 + nd2a;
  int i1 = tid - nd1 * (i2-nd2a) + nd1a;

  #define X(i1, i2, axis) x_d[ (i1-nd1a) + nd1 * (i2-nd2a) + axis*nd1*nd2]
  #define UC(i1,i2) uc_d[(i1-nd1a) + nd1 * (i2-nd2a)]
  
  if ( (i1 >= nd1a && i1 <=nd1b) && (i2 >= nd2a && i2 <= nd2b) )
  {
    UC(i1, i2) = UTRUE( X(i1,i2,0), X(i1,i2,1), t);
  }
}

// Forward Euler
__global__
void heat2sSolveGPUKernel(double *x_d, double *uc_d, double *un_d, double* rx_d, 
      double* ry_d, double t, double* dt_d, int* n1a_d, int* n1b_d, int* nd1_d,  
      int* nd1a_d, int* n2a_d, int* n2b_d, int* nd2_d, int* nd2a_d)
{
  int &n1a = *n1a_d, &n1b = *n1b_d;
  int& nd1 = *nd1_d;
  int& nd1a = *nd1a_d;

  int &n2a = *n2a_d, &n2b = *n2b_d;
  int& nd2 = *nd2_d;
  int& nd2a = *nd2a_d;

  double &rx = *rx_d, &ry = *ry_d, &dt = *dt_d;
  
  #define INDEX1(tx, bx) ((tx) + blockDim.x * (bx))
  int tid = INDEX1(threadIdx.x, blockIdx.x);
    // find i1 and i2
  int i2 = tid / nd1 + nd2a;
  int i1 = tid - nd1 * (i2-nd2a) + nd1a;

  #define U(i1,i2) uc_d[(i1-nd1a) + nd1 * (i2-nd2a)]
  #define UN(i1,i2) un_d[(i1-nd1a) + nd1 * (i2-nd2a)]

  if ( (i1 >= n1a && i1 <=n1b) && (i2 >= n2a && i2 <= n2b) )
  {
    UN(i1,i2) = U(i1,i2) + rx*( U(i1+1,i2) -2.*U(i1,i2) + U(i1-1,i2) )
          + ry*( U(i1,i2+1) -2.*U(i1,i2) + U(i1,i2-1) ) 
          + dt * FORCE( X(i1,i2,0), X(i1,i2,1), t );
  }
}


// Left and right side boundary condition
__global__
void applyLeftRightBC(double *x_d, double *un_d, double t, double* dt_d, 
          int* nd1_d, int* nd1a_d, int* nd2_d, int* nd2a_d, int* nd2b_d, 
          int i1, int is)
{
  int &nd1 = *nd1_d, &nd1a = *nd1a_d;
  int &nd2 = *nd2_d, &nd2a = *nd2a_d, &nd2b = *nd2b_d;

  double &dt = *dt_d;
  
  int i2 = threadIdx.x + blockDim.x  * blockIdx.x + nd2a;

  // macros for x and UN
  #define X(i1, i2, axis) x_d[ (i1-nd1a) + nd1 * (i2-nd2a) + axis*nd1*nd2]
  #define UN(i1,i2) un_d[(i1-nd1a) + nd1 * (i2-nd2a)]

  int i1g = i1 - is; // index of ghost point
  if ((nd2a <= i2) && (i2 <= nd2b))
  {
    UN(i1,i2) = UTRUE( X(i1,i2,0), X(i1,i2,1), t+dt);
    UN(i1g,i2) = 3.*UN(i1,i2) - 3.*UN(i1+is,i2) + UN(i1+2*is,i2); // extrap ghost
  }
}

__global__
void applyTopBottomBC(double *x_d, double *un_d, double t, double* dt_d, 
        int* nd1_d, int* nd1a_d, int* nd1b_d, int* nd2_d, int* nd2a_d, 
        int i2, int is)
{
  int &nd1 = *nd1_d, &nd1a = *nd1a_d, &nd1b = *nd1b_d;
  int &nd2 = *nd2_d, &nd2a = *nd2a_d;

  double &dt = *dt_d;
  
  int i1 = threadIdx.x + blockDim.x  * blockIdx.x + nd1a;

  // macros for x and UN
  #define X(i1, i2, axis) x_d[ (i1-nd1a) + nd1 * (i2-nd2a) + axis*nd1*nd2]
  #define UN(i1,i2) un_d[(i1-nd1a) + nd1 * (i2-nd2a)]

  int i2g = i2 - is; // index of ghost point
  if ((nd1a <= i1) && (i1 <= nd1b))
  {
    UN(i1, i2) = UTRUE( X(i1, i2, 0), X(i1, i2, 1), t+dt);
    UN(i1, i2g) = 3.*UN(i1,i2) - 3.*UN(i1,i2+is) + UN(i1,i2+2*is); // extrap ghost
  }
}

void heat2dSolveGPU(const int& debug, const int& nx, const int& ny, const int& numSteps,
  const int& numGhost, const int& n1a, const int& n1b, const int& nd1, const int& nd1a,
  const int& nd1b, const int& n2a, const int& n2b, const int& nd2, const int& nd2a,
  const int& nd2b, const int& saveMatlab, const double& tFinal, const double& rx, 
  const double& ry, const double& dt, const double& cpuTimeStep, double& maxErr,
  IntegerArray& gridIndexRange, IntegerArray& dimension, IntegerArray& boundaryCondition, 
  Range& Rx, Range& Ry, RealArray& x, const Real *dx, RealArray *ua, 
  const std::string& matlabFileNameGPU)
{
  const int Nt = 128;
  const int Nb = ceil( 1. * (nd1 * nd2) / Nt);

  printf("----- GPU: Solve the Heat Equation in two dimensions ------\n");
  printf("      saveMatlab=%d, matlabFileName=%s \n",saveMatlab, matlabFileNameGPU.c_str());
  printf("   kappa=%.3g, nx=%d, ny=%d, tFinal=%6.2f, debug=%d\n",kappa,nx,ny,tFinal, debug);

  // create the variable for device
  int *n1a_d, *n1b_d, *nd1_d, *nd1a_d, *nd1b_d; 
  int *n2a_d, *n2b_d, *nd2_d, *nd2a_d, *nd2b_d;
  cudaMalloc((void**)& n1a_d, sizeof(int));
  cudaMalloc((void**)& n1b_d, sizeof(int));
  cudaMalloc((void**)& nd1_d, sizeof(int));
  cudaMalloc((void**)& nd1a_d, sizeof(int));
  cudaMalloc((void**)& nd1b_d, sizeof(int));
  cudaMalloc((void**)& n2a_d, sizeof(int));
  cudaMalloc((void**)& n2b_d, sizeof(int));
  cudaMalloc((void**)& nd2_d, sizeof(int));
  cudaMalloc((void**)& nd2a_d, sizeof(int));
  cudaMalloc((void**)& nd2b_d, sizeof(int));

  cudaMemcpy(n1a_d, &n1a, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(n1b_d, &n1b, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(nd1_d, &nd1, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(nd1a_d, &nd1a, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(nd1b_d, &nd1b, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(n2a_d, &n2a, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(n2b_d, &n2b, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(nd2_d, &nd2, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(nd2a_d, &nd2a, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(nd2b_d, &nd2b, sizeof(int), cudaMemcpyHostToDevice);

  double *rx_d, *ry_d, *dt_d;
  cudaMalloc((void**)& rx_d, sizeof(double));
  cudaMalloc((void**)& ry_d, sizeof(double));
  cudaMalloc((void**)& dt_d, sizeof(double)); 

  cudaMemcpy(rx_d, &rx, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(ry_d, &ry, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dt_d, &dt, sizeof(double), cudaMemcpyHostToDevice);

  double *x_d, *u_d[2];
  const int size_of_array = sizeof(double) * nd1 * nd2;
  cudaMalloc((void**)& x_d, size_of_array * 2);
  cudaMalloc((void**)& u_d[0], size_of_array);
  cudaMalloc((void**)& u_d[1], size_of_array);

  double *x_p = x.getDataPointer();
  cudaMemcpy(x_d, x_p, size_of_array*2, cudaMemcpyHostToDevice);

  // inital condition
  Real t=0.;
  apply_inital_condition<<<Nb, Nt>>>(x_d, u_d[0], t, nd1_d,
          nd1a_d, nd1b_d, nd2_d, nd2a_d, nd2b_d);
  cudaDeviceSynchronize();

  double* u_p0 = ua[0].getDataPointer();
  cudaMemcpy(u_p0, u_d[0], size_of_array, cudaMemcpyDeviceToHost);


  // ---------- TIME-STEPPING LOOP ---------
  // create the cuda event for time 
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  
  int cur=0, n=0, i1=0, i2=0;
  
  cudaEventRecord(start, 0);
  for( n=0; n<numSteps; n++ )
  {
    t = n*dt;
    int next = (cur+1) % 2;

    heat2sSolveGPUKernel<<<Nb, Nt>>>(x_d, u_d[cur], u_d[next], rx_d, ry_d, t, dt_d,
                      n1a_d, n1b_d, nd1_d, nd1a_d, n2a_d, n2b_d, nd2_d, nd2a_d);
    cudaDeviceSynchronize();

    // --- boundary conditions ---
    for( int axis=0; axis<numberOfDimensions; axis++ )
      for( int side=0; side<=1; side++ )
      {
        int is = 1-2*side; // is=1 on left, is=-1 on right
        if( boundaryCondition(side,axis)==dirichlet )
        {
          if( axis==0 )
          { // left or right side
            i1= gridIndexRange(side,axis);
            applyLeftRightBC<<<Nb, Nt>>>(x_d, u_d[next], t, dt_d, nd1_d, nd1a_d,
                       nd2_d, nd2a_d, nd2b_d, i1, is);
            cudaDeviceSynchronize();
          }
          else
          { // bottom or top
            i2= gridIndexRange(side,axis);
            applyTopBottomBC<<<Nb, Nt>>>(x_d, u_d[next], t, dt_d, nd1_d, nd1a_d,
                        nd1b_d, nd2_d, nd2a_d, i2, is);
            cudaDeviceSynchronize();
          }

        }
        else
        {
          printf("ERROR: unknown boundaryCondition=%d\n",boundaryCondition(side,axis));
          abort();
        }

      } // end for axis  
    cur = next;                  
  }
  cudaEventRecord(stop, 0);
  cudaDeviceSynchronize();

  float gpuTimeStep; 
  cudaEventElapsedTime(&gpuTimeStep, start, stop); 
  gpuTimeStep *= 1E-03; // in sec
  const double speedup = cpuTimeStep / static_cast<double>(gpuTimeStep);

  double* un_p = ua[cur].getDataPointer();
  cudaMemcpy(un_p, u_d[cur], size_of_array, cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();

  // --- compute errors ---
  RealArray & uc = ua[cur];
  RealArray err(Rx,Ry);

  // maxErr is passed through reference as it will be 
  // written in matlab file executed from heat2d.C
  Real maxNorm=0.; 
  for( i2=n2a; i2<=n2b; i2++ )
    for( i1=n1a; i1<=n1b; i1++ )
    {
      err(i1,i2) = fabs(uc(i1,i2) - UTRUE(x(i1,i2,0),x(i1,i2,1),tFinal));
      maxErr = max(err(i1,i2),maxErr);
      maxNorm = max(uc(i1,i2),maxNorm);
    }
  maxErr /= max(maxNorm,REAL_MIN); // relative error

  printf("GPU: numSteps=%d nx=%d maxNorm=%8.2e maxRelErr=%8.2e gpuTimeStep=%9.2e(s) \n",
    numSteps,nx,maxNorm,maxErr, gpuTimeStep);
  printf(">>> Nt= %d: cpu= %12.4e, gpu= %12.4e, speedup= %12.4e \n\n", Nt, cpuTimeStep, gpuTimeStep, speedup);

}
