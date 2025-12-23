#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <string>
#include <math.h>

#include "getCPU.h"

#ifndef Real
#define Real double
#endif

constexpr Real pi = M_PI; // value of PI

// macros for x
#define x(i) x_p[i-nd1a]

// macros for boundary conditions
#define boundaryCondition(side,axis) boundaryCondition_p[(side)+2*(axis)]

constexpr Real kappa= .1;
constexpr Real kx = 3.;
constexpr Real kxPi = kx*pi;
constexpr Real kappaPiSq = kappa*kxPi*kxPi;

#define TRIG_DD 1
#define TRIG_NN 2
#define POLY_DD 3
#define POLY_NN 4

// ===== Choose the solution here or compile with -DSOLUTION=[1|2|3|4] =====
#ifndef SOLUTION
// #define SOLUTION TRIG_DD
// #define SOLUTION TRIG_NN
#define SOLUTION POLY_DD
// #define SOLUTION POLY_NN
#endif

#if SOLUTION == TRIG_DD

// True solution for dirichlet BC’s
constexpr char solutionName[] = "trueDD";
#define UTRUE(x,t) sin(kxPi*(x))*exp( -kappaPiSq*(t) )
#define UTRUEX(x,t) kxPi*cos(kxPi*(x))*exp( -kappaPiSq*(t) )
#define FORCE(x,t) (0.)

#elif SOLUTION == TRIG_NN

// True solution for Neumann BC’s
constexpr char solutionName[] = "trueNN";
#define UTRUE(x,t) cos(kxPi*(x))*exp( -kappaPiSq*(t) )
#define UTRUEX(x,t) -kxPi*sin(kxPi*(x))*exp( -kappaPiSq*(t) )
#define FORCE(x,t) (0.)

#elif (SOLUTION == POLY_DD) || (SOLUTION == POLY_NN)

// polynomial manufactured solution
 #if SOLUTION == POLY_DD
constexpr char solutionName[] = "polyDD";
#else
constexpr char solutionName[] = "polyNN";
#endif
constexpr Real b0=1., b1=.5, b2=.25;
constexpr Real a0=1., a1=.3;
#define UTRUE(x,t) (b0 + (x)*( b1 + (x)*b2 ))*( a0 + (t)*( a1 ) )
#define UTRUEX(x,t) ( b1 + 2.*(x)*b2 )*( a0 + (t)*( a1 ) )
#define UTRUET(x,t) (b0 + (x)*( b1 + (x)*b2 ))*( a1 )
#define UTRUEXX(x,t) ( 2.*b2 )*( a0 + (t)*( a1 ) )

// force = u_t - kappa*u.xx

#define FORCE(x,t) ( UTRUET(x,t) - kappa*UTRUEXX(x,t) )

#else
printf("ERROR: unknown solution");
abort();
#endif

// Macros to define fortran like arrays
#define uc(i) u_p[cur ][i-nd1a]
#define un(i) u_p[next][i-nd1a]

// for error_p macros to define fortran like arrays
#define error(i) error_p[i-nd1a]  

extern "C"
{
  void heat1dSolveGPU(const int& debug, const int& Nx, const int& numSteps, 
    const int& numGhost, const int& nd1, const int& n1a, const int& n1b, 
    const int& nd1a, const int& nd1b, const double& rx, const double& dx, 
    const double& tFinal, const double& dt, const double& cpuTimeStep, double* x_p, 
    const int *boundaryCondition_p, double** u_p, double* error_p);
}


// Macros to define fortran like arrays in DEVICE

// initialize the solution at t=0
__global__
void initial_solution(double *x_d, double *uc_d, double t, int nd1a, int nd1b)
{
  // Macros to define fortran like arrays in DEVICE
  #define X(i) x_d[i-nd1a]
  #define uc_d(i) uc_d[i-nd1a]
  #define un_d(i) un_d[i-nd1a]

  int tid = threadIdx.x;
  int i = blockIdx.x * blockDim.x + tid + nd1a;

  // n1a <= i <= n1b
  if ((nd1a <= i) && (i <= nd1b))
  {
    uc_d(i)=UTRUE(X(i),t);  
  }
}


__global__
void heat1dSolveGPUKernel( int* n1a_d, int* n1b_d, int* nd1a_d, int n,
      double* rx_d, double* dt_d, double* x_d, double *uc_d, double *un_d)
{
  int& n1a = *n1a_d;
  int& n1b = *n1b_d;
  int& nd1a = *nd1a_d;
  double& rx = *rx_d;
  double& dt = *dt_d;
  double t = dt * n;

  // get the index correxponds to block-ID and thread-ID
  int tid = threadIdx.x;
  int i = blockIdx.x * blockDim.x + tid + nd1a;

  // n1a <= i <= n1b
  if ((n1a <= i) && (i <= n1b))
  {
    un_d(i) = uc_d(i) + rx*( uc_d(i+1) - 2.*uc_d(i) + uc_d(i-1) ) + dt*FORCE( X(i),t );
  }
}

__global__
void applyDirichletBC(double *x_d, double *un_d, int* nd1a_d, int i, int is, double t_next)
{
    int &nd1a = *nd1a_d;
    un_d(i) = UTRUE(X(i),t_next);
    un_d(i-is) = 3.*un_d(i) - 3.*un_d(i+is) + un_d(i+2*is); // extrapolate ghost
}

__global__
void applyNeumannBC(double *x_d, double *un_d, double* dx_d, int* nd1a_d, 
                  int i, int is, double t_next)
{
  int &nd1a = *nd1a_d;
  un_d(i-is) = un_d(i+is) - 2.*is*(*dx_d)*UTRUEX(X(i),t_next);
}

void heat1dSolveGPU(const int& debug, const int& Nx, const int& numSteps, 
    const int& numGhost, const int& nd1, const int& n1a, const int& n1b, 
    const int& nd1a, const int& nd1b, const double& rx, const double& dx, 
    const double& tFinal, const double& dt, const double& cpuTimeStep, double* x_p, 
    const int *boundaryCondition_p, double** u_p, double* error_p)
{
  const int dirichlet=1, neumann=2; 

  const int Nt = 128;
  const int Nb = ceil(1.*nd1 /Nt);

  printf("------------------- GPU: Solve the heat equation in 1D solution=%s --------------------- \n",
    solutionName);
    printf("  Nb=%d, Nt=%d\n", Nb, Nt);
    printf("  numGhost=%d, n1a=%d, n1b=%d, nd1a=%d, nd1b=%d\n",numGhost,n1a,n1b,nd1a,nd1b);
    printf("  numSteps=%d, Nx=%d, kappa=%g, tFinal=%g, boundaryCondition(0,0)=%d, boundaryCondition(1,0)=%d\n",
            numSteps,Nx,kappa,tFinal,boundaryCondition(0,0),boundaryCondition(1,0));

  // first copy certain variables on GPU
  int *n1a_d, *n1b_d, *nd1a_d;
  cudaMalloc((void**)& n1a_d, sizeof(int)); 
  cudaMalloc((void**)& n1b_d, sizeof(int));
  cudaMalloc((void**)& nd1a_d, sizeof(int));

  const int size_of_array = sizeof(double) * nd1;
  double* rx_d, *dt_d, *x_d, *dx_d, *u_d[2];
  cudaMalloc((void**)& rx_d, sizeof(double));
  cudaMalloc((void**)& dt_d, sizeof(double));
  cudaMalloc((void**)& dx_d, sizeof(double));
  cudaMalloc((void**)& x_d,  size_of_array );
  cudaMalloc((void**)& u_d[0], size_of_array );
  cudaMalloc((void**)& u_d[1], size_of_array );
  
  cudaMemcpy(n1a_d, &n1a, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(n1b_d, &n1b, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(nd1a_d, &nd1a, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(rx_d, &rx, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dt_d, &dt, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dx_d, &dx, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(x_d, x_p, size_of_array, cudaMemcpyHostToDevice);

  Real t = 0.;
  int cur = 0;
  int i;
  // initialize the solution
  initial_solution<<<Nb, Nt>>>(x_d, u_d[cur], t, nd1a, nd1b);
  cudaDeviceSynchronize();
  cudaMemcpy(u_p[cur], u_d[cur], size_of_array, cudaMemcpyDeviceToHost);

  // ---------- TIME-STEPPING LOOP ---------
    // create the cuda event for time 
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  
  cudaEventRecord(start, 0);
  for( int n=0; n<numSteps; n++ )
  {
    t = n*dt; // current time

    const int cur = n % 2; // current time level
    const int next = (n+1) % 2; // next time level

    // --- update the interior points ----
    heat1dSolveGPUKernel<<<Nb, Nt>>>(n1a_d, n1b_d, nd1a_d, n, 
                                   rx_d, dt_d, x_d, u_d[cur], u_d[next]);
    cudaDeviceSynchronize();
    
    // ---- boundary conditions ----
    for( int side=0; side<=1; side++ )
    {
      const int i = side==0 ? n1a : n1b; // boundary index
      const int is = 1 - 2*side; // is = 1 on left, -1 on right
      if( boundaryCondition(side,0)==dirichlet )
      {
        applyDirichletBC<<<1, 1>>>(x_d, u_d[next], nd1a_d, i, is, t+dt);
      }
      else
      {
        // Neumann BC
        applyNeumannBC<<<1, 1>>>(x_d, u_d[next], dx_d, nd1a_d, i, is, t+dt);
      }
    }

    if( debug>1 )
    {
      cudaMemcpy(u_p[next], u_d[next], size_of_array, cudaMemcpyDeviceToHost);
      printf("step %d: After update interior and real BCs\n u=[",n+1);
      for( int i=nd1a; i<=nd1b; i++ )
          printf("%12.4e, ",un(i));
      printf("]\n");
    }

    if( debug>0 )
    {
      cudaMemcpy(u_p[next], u_d[next], size_of_array, cudaMemcpyDeviceToHost);
      // compute the error
      Real maxErr=0.;
      for( int i=nd1a; i<=nd1b; i++ )
      {
          Real err = fabs( un(i) - UTRUE(x(i),t+dt) );
          maxErr = max( maxErr,err );
      }
      printf("step=%d, t=%9.3e, maxErr=%9.2e\n",n+1,t+dt,maxErr);
    }

  } // end time-stepping loop
  cudaEventRecord(stop, 0);
  cudaDeviceSynchronize();

  // copy the u_p from device to host
  cudaMemcpy(u_p[cur], u_d[cur], size_of_array, cudaMemcpyDeviceToHost);

  float gpuTimeStep; 
  cudaEventElapsedTime(&gpuTimeStep, start, stop); 
  gpuTimeStep *= 1E-03; // in sec
  const double speedup = cpuTimeStep / static_cast<double>(gpuTimeStep);

  // ---- check the error -----
  t +=dt; // tFinal;
  if( fabs(t-tFinal) > 1e-3*dt/tFinal )
  {
      printf("ERROR: AFTER TIME_STEPPING: t=%16.8e IS NOT EQUAL to tFinal=%16.8e\n",t,tFinal);
  }

  cur = numSteps % 2;
  double maxErr=0.;
  // i is already declared
  for( i=nd1a; i<=nd1b; i++ )
  {
      error(i) = uc(i) - UTRUE(x(i),t);
      maxErr = max( maxErr, abs(error(i)) );
  }

  printf("numSteps=%4d, Nx=%3d, maxErr=%9.2e, gpu=%9.2e(s)\n", numSteps, Nx, maxErr, gpuTimeStep);
  printf(">>> Nt=%d, cpuTime= %12.4e(s), gpuTime= %12.4e(s), speedup = %12.4e \n", 
        Nt, cpuTimeStep, gpuTimeStep, speedup);

  cudaFree(x_d);
  cudaFree(u_d[0]);
  cudaFree(u_d[1]);
}