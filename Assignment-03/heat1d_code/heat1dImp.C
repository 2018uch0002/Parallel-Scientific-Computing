//===================================================================================
//
// Solve the heat equation in one-dimension with an IMPLICT METHOD
// TRAPEZOIDAL RULE IN TIME
//
//===================================================================================

#include "A++.h"

// Tridiagonal factor and solve:
#include "tridiagonal.h"

// define some types
typedef double Real;
typedef doubleSerialArray RealArray;
typedef intSerialArray IntegerArray;

#include <string>
using std::string;
using std::max;

// getCPU() : Return the current wall-clock time in seconds
#include "getCPU.h"

// include commands to parse command line arguments
#include "parseCommand.h"


// Function to save a vector to a matlab file.
#include "writeMatlabVector.h"

// True solution options:
static const int trueDD = 1;
static const int trueNN = 2;
static const int poly = 3;

static Real kappa = .1;
static Real kx = 3.;
static Real kxPi = kx*M_PI;
static Real kappaPiSq = kappa*kxPi*kxPi;

// True solution for dirichlet BC's
#define UTRUEDD(x,t) sin(kxPi*(x))*exp( -kappaPiSq*(t) )
#define UTRUEDDx(x,t) kxPi*cos(kxPi*(x))*exp( -kappaPiSq*(t) )
#define FORCEDD(x,t) (0.)

// True solution for Neumann BC's
#define UTRUENN(x,t) cos(kxPi*(x))*exp( -kappaPiSq*(t) )
#define UTRUENNx(x,t) -kxPi*sin(kxPi*(x))*exp( -kappaPiSq*(t) )
#define FORCENN(x,t) (0.)

// polynomial manufactured solution
static const Real b0=1., b1=.5, b2=.25;
static const Real a0=1., a1=.3, a2=-.1;
#define POLY(x,t) (b0 + (x)*( b1 + (x)*b2 ))*( a0 + (t)*( a1 + (t)*a2 ) )
#define POLYx(x,t) ( b1 +2.*(x)*b2 )*( a0 + (t)*( a1 + (t)*a2 ) )
#define POLYT(x,t) (b0 + (x)*( b1 +(x)*b2 ))*( a1 + 2.*(t)*a2 )
#define POLYXx(x,t) ( 2.*b2 )*( a0 + (t)*( a1 + (t)*a2 ) )
#define POLYFORCE(x,t) ( POLYT(x,t) - kappa*POLYXx(x,t) )

// ---------------------------------------------------------------------------------------
// Function to evaluate the true solution
// solutionOption (input) : true solution option
// t (input) : time
// x (input) : grid points
// I1 (input) : evaluate at these index values
// uTrue (output) : uTrue(I1) = true solution
// ---------------------------------------------------------------------------------------
int getTrue( int solutionOption, Real t, RealArray & x, Index & I1, RealArray & uTrue )
{
    if( solutionOption==trueDD )
    uTrue(I1) = UTRUEDD(x(I1),t);
    else if( solutionOption==trueNN )
    uTrue(I1) = UTRUENN(x(I1),t);
    else if( solutionOption==poly )
    uTrue(I1) = POLY(x(I1),t);
    else
    {
    printf("getTrue: unknown solutonOption=%d\n",solutionOption);
    abort();
    }

    return 0;
}

// ---------------------------------------------------------------------------------------
// Function to evaluate the x-derivative of the true solution
// ---------------------------------------------------------------------------------------
int getTruex( int solutionOption, Real t, RealArray & x, Index & I1, RealArray & uTruex )
{
    if( solutionOption==trueDD )
    uTruex(I1) = UTRUEDDx(x(I1),t);
    else if( solutionOption==trueNN )
    uTruex(I1) = UTRUENNx(x(I1),t);
    else if( solutionOption==poly )
    uTruex(I1) = POLYx(x(I1),t);
    else
    {
    printf("getTrue: unknown solutonOption=%d\n",solutionOption);
    abort();
    }
    return 0;
}

// ---------------------------------------------------------------------------------------
// Function to evaluate the PDE forcing
// ---------------------------------------------------------------------------------------
int getForce( int solutionOption, Real t, RealArray & x, Index & I1, RealArray & force )
{
    if( solutionOption==trueDD )
    force(I1) = FORCEDD(x(I1),t);
    else if( solutionOption==trueNN )
    force(I1) = FORCENN(x(I1),t);
    else if( solutionOption==poly )
    force(I1) = POLYFORCE(x(I1),t);
    else
    {
    printf("getTrue: unknown solutonOption=%d\n",solutionOption);
    abort();
    }
    return 0;
}



int main( int argc, char* argv[] )
{

    printf("Usage: heat1dImp -Nx=<i> -tFinal=<f> -sol=[true|poly] -bc1=[d|n] -bc2=[d|n] -debug=<i> - matlabFileName=<s>\n");

    const Real pi = M_PI;

    const int numberOfDimensions=1;
    int Nx=10;

    Real xa=0., xb=1.;
    Real tFinal = 1.;

    // setup boundary condition array
    const int dirichlet=1, neumann=2;
    IntegerArray boundaryCondition(2,numberOfDimensions);
    boundaryCondition(0,0)=dirichlet; // left
    boundaryCondition(1,0)=dirichlet; // right

    int debug=0;
    string matlabFileName = "heat1d.m";

    const int trueSolution=0, polynomialSolution=1;
    int sol = trueSolution;
    string solName = "true";

    string line;
    for( int i=1; i<argc; i++ )
    {
    line=argv[i];
    // printf("Input: argv[%d] = [%s]\n",i,line.c_str());

    if( parseCommand( line,"-Nx=",Nx) ){}
    else if( parseCommand( line,"-debug=",debug) ){}
    else if( parseCommand( line, "-tFinal=",tFinal) ){}
    else if( parseCommand( line,"-matlabFileName=",matlabFileName) ){}
    else if( line.substr(0,4)=="-sol" )
    {
    solName = line.substr(5); // solName = character 5 to end
    if( solName=="true" )
    sol=trueSolution;
    else if( solName=="poly" )
    sol=polynomialSolution;
    else
    {
    printf("ERROR: Unknown -sol=[%s]\n",solName.c_str());
    abort();
    }
    printf("setting solName=[%s]\n",solName.c_str());
    }
    else if( line.substr(0,5)=="-bc1=" ||
    line.substr(0,5)=="-bc2=" )
    {
    int side = line.substr(0,5)=="-bc1=" ? 0 : 1;
    string bcName = line.substr(5,1);
    if( bcName=="d" )
    {
    boundaryCondition(side,0)=dirichlet;
    printf(" SETTING boundaryCondition(%d,0)=dirichlet\n",side);
    }
    else if( bcName=="n" )
    {
    boundaryCondition(side,0)=neumann;
    printf(" SETTING boundaryCondition(%d,0)=neumann\n",side);
    }
    else
    {
    printf("Uknown bc: line=[%s]\n",line.c_str());
    }
    }
    }
    string bcName[2];
    for( int side=0; side<=1; side++ )
    {
    if( boundaryCondition(side,0)==dirichlet )
    bcName[side]="D";
    else
    bcName[side]="N";
    }

    // set the solutionOption:
    int solutionOption = trueDD;
    if( boundaryCondition(0,0)==dirichlet && boundaryCondition(1,0)==dirichlet )
    {
    solutionOption = sol==trueSolution ? trueDD : poly;
    }
    else if( boundaryCondition(0,0)==neumann && boundaryCondition(1,0)==neumann )
    {
    solutionOption = sol==trueSolution ? trueNN : poly;
    }
    else
    {
    printf("Unexpected boundary conditions with the true solution -- FINISH ME\n");
    abort();
    }
    const string solutionName = solName + bcName[0] + bcName[1];

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

    int nd1 = nd1b-nd1a+1; // total number of grid points;

    // Grid points:
    Range Rx(nd1a,nd1b);
    RealArray x(Rx);

    for( int i=nd1a; i<=nd1b; i++ )
    x(i) = xa + (i-n1a)*dx;

    if( debug>1 )
    {
    for( int i=nd1a; i<=nd1b; i++ )
    printf("x(%2d)=%12.4e\n",i,x(i));
    }

    RealArray u[2]; // two arrays will be used for current and new times
    u[0].redim(Rx);
    u[1].redim(Rx);

    // Macros to define fortran like arrays
    #define uc(i) u[cur ](i)
    #define un(i) u[next](i)

    // initial conditions
    Real t=0.;
    int cur = 0; // "current" solution, index into u_p[]
    Index I1 = Range(nd1a,nd1b);
    getTrue( solutionOption, t, x, I1, u[cur] );

    if( debug>0 )
    {
    printf("After initial conditions\n u=[");
    for( int i=nd1a; i<=nd1b; i++ )
    printf("%10.4e, ",uc(i));
    printf("]\n");
    }


    // Choose time-step
    const Real dx2 = dx*dx;
    Real dt = dx; // dt, adjusted below
    const int numSteps = ceil(tFinal/dt);
    dt = tFinal/numSteps; // adjust dt to reach the final time

    // --- Tridiagonal matrix ----
    // Ax(0:2,i1) : holds the 3 diagonals
    // [ Ax(1,0) Ax(2,0) ]
    // [ Ax(0,1) Ax(1,1) Ax(2,1) ]
    // [ Ax(0,2) Ax(1,2) Ax(2,2) ]
    //
    Range Ix(n1a,n1b); // interior and boundary points
    RealArray Ax(Range(3),Ix);
    RealArray rhsx(Ix);
    Real *rhsx_p = rhsx.getDataPointer();
    #define RHS(i) rhsx_p[i-n1a]

    // ---- Fill the tridiagonal matrix ----
    const Real rx = kappa*dt/dx2;
    for( int i1=n1a+1; i1<=n1b-1; i1++ )
    {
    Ax(0,i1) = - .5*rx; // lower diagonal
    Ax(1,i1) = 1. + rx; // diagonal
    Ax(2,i1) = - .5*rx; // upper diagonal
    }
    for( int side=0; side<=1; side++ )
    {
    int i1 = side==0 ? n1a : n1b;
    if( boundaryCondition(side,0)==dirichlet )
    { // Dirichlet BC
    Ax(0,i1)=0.; Ax(1,i1)=1.; Ax(2,i1)=0.;
    }
    else
    { // Neumann BC
    // Combine Neumann BC with interior equation on the boundary
    // to eliminate the ghost point
    int is = 1-2*side;
    Ax(1-is,i1)=0.; Ax(1,i1)=1.+rx; Ax(1+is,i1)=-rx; // Neumann BC

    }
    }

    // Factor matrix ONCE
    factorTridiagonalMatrix( Ax );

    printf("--------- Implicit Solve of the Heat Equation in 1D, solutionOption=%d, solutionName=%s ----------- \n",
    solutionOption,solutionName.c_str());
    printf("  numSteps=%d, Nx=%d, debug=%d, kappa=%g, tFinal=%g, \n"
    "  boundaryCondition(0,0)=%s, boundaryCondition(1,0)=%s\n",
    numSteps,Nx,debug,kappa,tFinal,bcName[0].c_str(),bcName[1].c_str());

    RealArray uTrue(Rx); // store uTrue here
    RealArray fn[2]; // save forcing at two time levels
    fn[0].redim(Rx);
    fn[1].redim(Rx);

    getForce( solutionOption, t, x, I1, fn[0] );

    // ---------- TIME-STEPPING LOOP ---------
    Real cpu0 = getCPU();
    for( int n=0; n<numSteps; n++ )
    {
    t = n*dt; // current time

    int cur = n % 2; // current time level
    int next = (n+1) % 2; // next time level

    // --- Assign the RHS ----

    // get forcing at t+dt
    getForce( solutionOption, t+dt, x, I1, fn[next] );

    Real *uc_p = u[cur].getDataPointer();
    #define UC(i) uc_p[i-nd1a]
    Real *un_p = u[next].getDataPointer();
    #define UN(i) un_p[i-nd1a]

    Real *fc_p = fn[cur].getDataPointer();
    #define FC(i) fc_p[i-nd1a]
    Real *fn_p = fn[next].getDataPointer();
    #define FN(i) fn_p[i-nd1a]


    for( int i=n1a; i<=n1b; i++ )
    {
    RHS(i) = UC(i) + (.5*rx)*( UC(i+1) - 2.*UC(i) + UC(i-1) ) + (.5*dt)*( FC(i) + FN(i) );
    }

    // ---- RHS for boundary conditions ----
    for( int side=0; side<=1; side++ )
    {
    const int i1 = side==0 ? n1a : n1b;
    const int is = 1-2*side;
    Index Ib = Range(i1,i1);

    if( boundaryCondition(side,0)==dirichlet )
    {
    getTrue( solutionOption, t+dt, x, Ib, uTrue );
    RHS(i1) = uTrue(i1);
    }
    else
    {
    // The Neumann BC, (U(i1+1)-U(i1-1))/(2*dx) = g,
    // is combine with the equation on th boundary to eliminate the
    // ghost point value. This results in an adjustment to the rhs. (see class notes)
    getTruex( solutionOption, t+dt, x, Ib, uTrue ); // store ux in uTrue(i1)
    RHS(i1) = RHS(i1) - (is*dx*rx)*uTrue(i1);
    }
    }

    solveTridiagonal( Ax, rhsx );

    // fill in the solution
    for( int i=n1a; i<=n1b; i++ )
    UN(i) = RHS(i);

    // Fill ghost point values
    for( int side=0; side<=1; side++ )
    {
    const int i1 = side==0 ? n1a : n1b;
    const int is = 1-2*side;
    Index Ib = Range(i1,i1);
    if( boundaryCondition(side,0)==dirichlet )
    {
    // Use extrapolation:
    UN(i1-is) = 3.*UN(i1) - 3.*UN(i1+is) + UN(i1+2*is);
    }
    else
    {
    // Neumann : (U(i1+1)-U(i1-1))/(2*dx) = g
    // getTruex( solutionOption, t+dt, x, Ib, uTrue ); // no need to re-evaluate
    UN(i1-is) = UN(i1+is) -2.*is*dx*uTrue(i1);
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
    Real maxNorm = 0.;

    getTrue( solutionOption, t+dt, x, I1, uTrue );

    for( int i=nd1a; i<=nd1b; i++ )
    {
    Real err = fabs( un(i) - uTrue(i) );
    // Real err = fabs( un(i) - uTrue(x(i),t+dt) );
    maxErr = max( maxErr,err );
    maxNorm = max( maxNorm, abs(un(i)) );
    }
    maxErr /= maxNorm;
    printf("step=%d, t=%9.3e, maxNorm=%9.2e, maxRelErr=%9.2e\n",n+1,t+dt,maxNorm,maxErr);

    }

    }

    Real cpuTimeStep = getCPU()-cpu0;

    // check the error:
    t +=dt; // tFinal;
    if( fabs(t-tFinal) > 1e-3*dt/tFinal )
    {
    printf("ERROR: AFTER TIME_STEPPING: t=%16.8e IS NOT EQUAL to tFinal=%16.8e\n",t,tFinal);
    }

    RealArray err(Rx);

    cur = numSteps % 2;
    Real maxErr=0.;
    Real maxNorm = 0.;

    getTrue( solutionOption, t, x, I1, uTrue );

    // Compute final error: do not include unused ghost point on Dirichlet boundaries
    const int n1ae = boundaryCondition(0,0)==dirichlet ? n1a : n1a-1;
    const int n1be = boundaryCondition(1,0)==dirichlet ? n1b : n1b+1;
    for( int i=n1ae; i<=n1be; i++ )
    {
    err(i) = uc(i) - uTrue(i);
    maxErr = max( maxErr, abs(err(i)) );
    maxNorm = max( maxNorm, abs(uc(i)) );
    }
    maxErr /= maxNorm;


    printf("numSteps=%4d, Nx=%3d, maxNorm=%9.2e, maxRelErr=%9.2e, cpu=%9.2e(s)\n",
    numSteps,Nx,maxNorm,maxErr,cpuTimeStep);


    // --- Write a matlab file for plotting in matlab ---
    FILE *matlabFile = fopen(matlabFileName.c_str(),"w");
    fprintf(matlabFile,"%% File written by heat1dImp.C\n");
    fprintf(matlabFile,"xa=%g; xb=%g; kappa=%g; t=%g; maxErr=%10.3e; cpuTimeStep=%10.3e;\n",xa,xb,kappa,tFinal,maxErr,cpuTimeStep);
    fprintf(matlabFile,"Nx=%d; dx=%14.6e; numGhost=%d; n1a=%d; n1b=%d; nd1a=%d; nd1b=%d;\n",Nx,dx,numGhost,n1a,n1b,nd1a,nd1b);
    fprintf(matlabFile,"solutionName=\%s\';\n",solutionName.c_str());

    writeMatlabVector( matlabFile, x, "x", nd1a, nd1b );
    writeMatlabVector( matlabFile, u[cur], "u", nd1a, nd1b );
    writeMatlabVector( matlabFile, err, "err", nd1a, nd1b );

    fclose(matlabFile);
    printf("Wrote file [%s]\n",matlabFileName.c_str());

    return 0;
}
