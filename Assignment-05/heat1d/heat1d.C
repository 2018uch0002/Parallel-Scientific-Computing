//
// Solve the heat equation in one-dimension
//

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>

// Include file for OpenMP
#include <omp.h>


// define a new type "Real" which is equivalent to a "double"
typedef double Real;

#include <string>
using std::string;
using std::max;

// getCPU() : Return the current wall-clock time in seconds
#include "getCPU.h"

// include commands tp parse command line arguments
#include "parseCommand.h"


// ---------------------------------------------------------------------------------------
// Function to save a vector to a matlab file.
// matlabFile (input) : save vector to this file
// u_p (input) : array of vector values
// name (input) : name for array
// (nd1a:nd1b) (input) : array dimensions
// ---------------------------------------------------------------------------------------
int writeMatlabVector( FILE *matlabFile, Real *u_p, const char *name, int nd1a, int nd1b )
{
#define u(i) u_p[i-nd1a]

const int numPerLine=8; // number of entries per line
// Save the vector as:
// name = [ num num num num num ...
// num num num num num ];
fprintf(matlabFile,"%s=[",name);
for( int i=nd1a; i<=nd1b; i++ )
{
fprintf(matlabFile,"%20.15e ",u(i));
if( (i-nd1a) % numPerLine == numPerLine-1 )
fprintf(matlabFile,"...\n"); // continuation line
}
fprintf(matlabFile,"];\n");

return 0;
}

int main( int argc, char* argv[] )
{
printf("Usage: heat1d -nx=<i> -tFinal=<i> -matlabFileName=<s> -numThreads=<i>\n"
        " Nx = number of grid cells.\n"
        " matlabFileName.m : save results to this file.\n");

#define TRIG_DD 1
#define TRIG_NN 2
#define POLY_DD 3
#define POLY_NN 4
// ===== Choose the solution here or compile with -DSOLUTION=[1|2|3|4] =====
#ifndef SOLUTION
#define SOLUTION TRIG_DD
// #define SOLUTION TRIG_NN
// #define SOLUTION POLY_DD
// #define SOLUTION POLY_NN
#endif

const Real pi = M_PI;

int debug=0; // set to 1 for debug info
Real xa=0., xb=1.;
Real kappa= .1;
Real tFinal = .01;
Real cfl=.9; // time-step safety factor

int Nx=10; // default
string matlabFileName = "heat1d.m";
int numThreads=1;
int version=1;

string line;
for( int i=1; i<argc; i++ )
{
line=argv[i];
if( parseCommand( line,"-nx=",Nx) ){}
else if( parseCommand( line,"-numThreads=",numThreads) ){}
else if( parseCommand( line,"-version=",version) ){}
else if( parseCommand( line, "-tFinal=",tFinal) ){}
else if( parseCommand( line,"-matlabFileName=",matlabFileName) ){}
}

// ============= Grid and indexing==============
// xa xb
// G---X---+---+---+---+-- ... ---+---X---G
// 0 1 2 Nx
// n1a n1b
// nd1a nd1b
// C index: 0 1 2 3 ...

Real dx = (xb-xa)/Nx;
const int numGhost=1;
const int n1a = 0;
const int n1b = Nx;
const int nd1a=n1a-numGhost;
const int nd1b=n1b+numGhost;
const int nd1 = nd1b-nd1a+1; // total number of grid points;

// Create an array of grid points:
Real *x_p = new Real [nd1];
#define x(i) x_p[i-nd1a]

for( int i=nd1a; i<=nd1b; i++ )
x(i) = xa + (i-n1a)*dx;

if( debug>1 )
{
for( int i=nd1a; i<=nd1b; i++ )
printf("x(%2d)=%12.4e\n",i,x(i));
}

const int dirichlet=1, neumann=2;
const int numberOfDimensions=1;
int *boundaryCondition_p = new int [2*numberOfDimensions];
#define boundaryCondition(side,axis) boundaryCondition_p[(side)+2*(axis)]


const Real kx = 3.;
const Real kxPi = kx*pi;
const Real kappaPiSq = kappa*kxPi*kxPi;

#if SOLUTION == TRIG_DD
// True solution for dirichlet BC's
boundaryCondition(0,0) = dirichlet;
boundaryCondition(1,0) = dirichlet;

const char solutionName[] = "trueDD";

#define UTRUE(x,t) sin(kxPi*(x))*exp( -kappaPiSq*(t) )
#define UTRUEX(x,t) kxPi*cos(kxPi*(x))*exp( -kappaPiSq*(t) )
#define FORCE(x,t) (0.)

#elif SOLUTION == TRIG_NN

// True solution for Neumann BC's
boundaryCondition(0,0) = neumann;
boundaryCondition(1,0) = neumann;
const char solutionName[] = "trueNN";

#define UTRUE(x,t) cos(kxPi*(x))*exp( -kappaPiSq*(t) )
#define UTRUEX(x,t) -kxPi*sin(kxPi*(x))*exp( -kappaPiSq*(t) )
#define FORCE(x,t) (0.)

#elif (SOLUTION == POLY_DD) || (SOLUTION == POLY_NN)

// polynomial manufactured solution
#if SOLUTION == POLY_DD
const char solutionName[] = "polyDD";
boundaryCondition(0,0) = dirichlet;
boundaryCondition(1,0) = dirichlet;
#else
const char solutionName[] = "polyNN";
boundaryCondition(0,0) = neumann;
boundaryCondition(1,0) = neumann;
#endif

const Real b0=1., b1=.5, b2=.25;
const Real a0=1., a1=.3;
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


Real *u_p[2]; // two arrays will be used for current and new times
u_p[0] = new Real [nd1];
u_p[1] = new Real [nd1];

// Macros to define fortran like arrays
#define uc(i) u_p[cur ][i-nd1a]
#define un(i) u_p[next][i-nd1a]

// initial conditions
Real t=0.;
int cur = 0; // "current" solution, index into u_p[]
int i;
#pragma omp parallel for default(shared) private(i) num_threads(numThreads)
for( i=nd1a; i<=nd1b; i++ )
uc(i)=UTRUE(x(i),t);

if( debug>0 )
{
printf("After initial conditions\n u=[");
for( int i=nd1a; i<=nd1b; i++ )
printf("%10.4e, ",uc(i));
printf("]\n");
}

// Time-step restriction is kappa*dt/dx^2 < .5
const Real dx2 = dx*dx;
Real dt = cfl*.5*dx2/kappa; // dt, adjusted below
const int numSteps = ceil(tFinal/dt);
dt = tFinal/numSteps; // adjust dt to reach the final time
const Real rx = kappa*dt/dx2;

printf("------------------- Solve the heat equation in 1D solution=%s --------------------- \n",
solutionName);
printf(" **version=%d**\n",version);
printf(" numGhost=%d, n1a=%d, n1b=%d, nd1a=%d, nd1b=%d, numThreads=%d\n",numGhost,n1a,n1b,nd1a,nd1b,numThreads);
printf(" numSteps=%d, Nx=%d, kappa=%g, tFinal=%g, boundaryCondition(0,0)=%d, boundaryCondition(1,0)=%d\n",
numSteps,Nx,kappa,tFinal,boundaryCondition(0,0),boundaryCondition(1,0));

Real cpu0 = getCPU();
if( version==1 )
{
// OMP version 1
// ---------- TIME-STEPPING LOOP ---------
for( int n=0; n<numSteps; n++ )
{
t = n*dt; // current time

const int cur = n % 2; // current time level
const int next = (n+1) % 2; // next time level

// --- update the interior points ----
int i;
#pragma omp parallel for default(shared) private(i) num_threads(numThreads)
for( i=n1a; i<=n1b; i++ )
{
un(i) = uc(i) + rx*( uc(i+1) - 2.*uc(i) + uc(i-1) ) + dt*FORCE( x(i),t );
}

// ---- boundary conditions ----
for( int side=0; side<=1; side++ )
{
const int i = side==0 ? n1a : n1b; // boundary index
const int is = 1 - 2*side; // is = 1 on left, -1 on right
if( boundaryCondition(side,0)==dirichlet )
{
un(i) = UTRUE(x(i),t+dt);
un(i-is) = 3.*un(i) - 3.*un(i+is) + un(i+2*is); // extrapolate ghost
}
else
{
// Neumann BC
un(i-is) = un(i+is) - 2.*is*dx*UTRUEX(x(i),t+dt);
}
}

if( debug>1 )
{
printf("step %d: After update interior and real BCs\n u=[",n+1);
for( int i=nd1a; i<=nd1b; i++ )
printf("%12.4e, ",un(i));
printf("]\n");
}

if( debug>0 )
{
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
}
else if( version==2 )
{
// OMP version 2
// Question: Is is faster to start the parallel region outside the time-stepping loop so we
// don't incur the over-head of creating threads every time step?
// Answer: This version below seems slower

// ---------- TIME-STEPPING LOOP ---------
int i,n,cur,next;
#pragma omp parallel default(shared) private(n,i) num_threads(numThreads)
{
for( n=0; n<numSteps; n++ )
{

#pragma omp single
{
t = n*dt; // current time

cur = n % 2; // current time level
next = (n+1) % 2; // next time level
}

#pragma omp barrier

// --- update the interior points ----
#pragma omp for
for( i=n1a; i<=n1b; i++ )
{
un(i) = uc(i) + rx*( uc(i+1) - 2.*uc(i) + uc(i-1) ) + dt*FORCE( x(i),t );
}

#pragma omp barrier

// ---- boundary conditions ----
#pragma omp single
{
for( int side=0; side<=1; side++ )
{
const int i = side==0 ? n1a : n1b; // boundary index
const int is = 1 - 2*side; // is = 1 on left, -1 on right
if( boundaryCondition(side,0)==dirichlet )
{
un(i) = UTRUE(x(i),t+dt);
un(i-is) = 3.*un(i) - 3.*un(i+is) + un(i+2*is); // extrapolate ghost
}
else
{
// Neumann BC
un(i-is) = un(i+is) - 2.*is*dx*UTRUEX(x(i),t+dt);
}
}
}


} // end time-stepping loop
}


}
else
{
printf("ERROR: unknown version=%d\n",version);
abort();
}
Real cpuTimeStep = getCPU()-cpu0;
// ---- check the error -----
t +=dt; // tFinal;
if( fabs(t-tFinal) > 1e-3*dt/tFinal )
{
printf("ERROR: AFTER TIME_STEPPING: t=%16.8e IS NOT EQUAL to tFinal=%16.8e\n",t,tFinal);
}

Real *error_p = new Real [nd1];
#define error(i) error_p[i-nd1a]

cur = numSteps % 2;
Real maxErr=0.;
for( int i=nd1a; i<=nd1b; i++ )
{
error(i) = uc(i) - UTRUE(x(i),t);
maxErr = max( maxErr, abs(error(i)) );
}

printf("numThreads=%2d, numSteps=%4d, Nx=%3d, maxErr=%9.2e, cpu=%9.2e(s)\n",numThreads,numSteps,Nx,maxErr,cpuTimeStep);

// --- Write a file for plotting in matlab ---
FILE *matlabFile = fopen(matlabFileName.c_str(),"w");
fprintf(matlabFile,"%% File written by heat1d.C\n");
fprintf(matlabFile,"xa=%g; xb=%g; kappa=%g; t=%g; maxErr=%10.3e; cpuTimeStep=%10.3e;\n",xa,xb, 
        kappa,tFinal,maxErr,cpuTimeStep);
fprintf(matlabFile,"Nx=%d; dx=%14.6e; numGhost=%d; n1a=%d; n1b=%d; nd1a=%d; nd1b=%d; numThreads=%d;\n",Nx,dx,numGhost,n1a,n1b,nd1a,nd1b,numThreads);
fprintf(matlabFile,"solutionName=\'%s\';\n",solutionName);

writeMatlabVector( matlabFile, x_p, "x", nd1a, nd1b );
writeMatlabVector( matlabFile, u_p[cur], "u", nd1a, nd1b );
writeMatlabVector( matlabFile, error_p, "err", nd1a, nd1b );

fclose(matlabFile);
printf("Wrote file %s\n\n",matlabFileName.c_str());

delete [] u_p[0];
delete [] u_p[1];
delete [] x_p;
delete [] error_p;
delete [] boundaryCondition_p;

return 0;
}
