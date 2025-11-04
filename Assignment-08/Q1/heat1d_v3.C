//
// Solve the heat equation in one-dimension
// maxErr= 6.27e-02
//

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>

// MPI header file
#include <mpi.h>

// to parse the command ouput
#include "parseCommand.h" 

// to get the local index bounds
#include "getLocalIndexBounds.h"

// define a new type "Real" which is equivalent to a "double"
typedef double Real;

#include <string>
using std::string;
using std::max;

#include <ctime>

// to get the cpu time that uses OpenMP
#include "getCPU.h"


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
    MPI_Init(&argc, &argv);

    int myRank; 
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get the rank of the processor

    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np); // get the total number of the processor

    if (myRank == 0){
        printf("Usage: heat1d [nx] [matlabFileName.m]\n"
        "       nx = number of grid cells.\n"
        "       tFinal = final-time\n"
        "       matlabFileName.m : save results to this file.\n"
        "       debug = \n");
        
    }

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
    int debug = 0; // set to 1 for debug info

    Real xa=0., xb=1.;
    Real kappa= .1;
    Real tFinal = .2;
    Real cfl=.9; // time-step safety factor

    int Nx=10; // default
    string matlabFileName = "heat1d.m";
    string line;
    for(int n_arg=1; n_arg < argc; n_arg++){
        line = argv[n_arg];
        if(parseCommand(line, "-nx=", Nx, myRank == 0)){}
        else if(parseCommand(line, "-tFinal=", tFinal, myRank == 0)){}
        else if(parseCommand(line, "-matlabFileName=", matlabFileName, myRank == 0)){}
        else if(parseCommand(line, "-debug=", debug, myRank == 0)){}
    }

    
    FILE *debugFile=NULL;
    if( debug>0 )
    {
        // open a debug file on each processor (in the debug folder)
        char debugFileName[80];
        sprintf(debugFileName,"debug/debugFileDemoNp%dProc%d.debug",np,myRank);
        debugFile = fopen(debugFileName,"w");

        // Write header info to both stdout and debugFile  
        for( int ifile=0; ifile<=1; ifile++ )
        {
            FILE *file = ifile==0 ? stdout : debugFile;
            // if( ( ifile==0 && myRank==0 ) || // write to terminal if myRank==0
            //     ( ifile==1 && debug>0 ) ) // write to debugFile if debug>0
            if( ( ifile==1 && debug>0 ) ) // write to debugFile if debug>0
            {
            fprintf(file,"------------------- DebugFileDemo --------------------- \n");
            fprintf(file," np=%d, myRank=%d\n",np,myRank);
            }
        }
    }

    // if( argc>=2 ) // read any command line arguments
    // {
    // Nx= atoi(argv[1]);
    // printf("Setting Nx=%d\n",Nx);
    // if( argc>=3 )
    // {
    // matlabFileName = argv[2];
    // printf("Setting matlabFileName=[%s]\n",matlabFileName.c_str());
    // }
    // }

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
    
    // get the local bounds  
    int nx_l = 0, n1a_l = 0, n1b_l = 0;
    getLocalIndexBounds(myRank, np, Nx, nx_l, n1a_l, n1b_l);
    const int nd1a_l = n1a_l - numGhost;
    const int nd1b_l = n1b_l + numGhost;
    const int nd1_l = nd1b_l - nd1a_l + 1; // total number of grid points per processor

    // Create an array of grid points:
    Real *x_p = NULL;
    if (myRank == 0)
        x_p = new Real [nd1];
    else
        x_p = new Real [nd1_l];
    
    #define x(i) x_p[i-nd1a_l]
    
    for( int i=nd1a_l; i<=nd1b_l; i++ )
        x(i) = xa + (i-n1a)*dx;

    if( debug>1 )
    {
        for( int i=nd1a; i<=nd1b; i++ )
            fprintf(debugFile, "x(%2d)=%12.4e\n",i,x(i));
    }

    const int dirichlet=1, neumann=2, parallel_ghost = 3;
    const int numberOfDimensions=1;
    int *boundaryCondition_p = new int [2*numberOfDimensions];
    #define boundaryCondition(side,axis) boundaryCondition_p[(side)+2*(axis)]


    const Real kx = 3.;
    const Real kxPi = kx*pi;
    const Real kappaPiSq = kappa*kxPi*kxPi;

    #if SOLUTION == TRIG_DD
    // True solution for dirichlet BC’s
    if (myRank == 0)
    {
        boundaryCondition(0,0) = dirichlet;
        if (np != 1)
            boundaryCondition(1,0) = parallel_ghost;
        else 
            boundaryCondition(1,0) = dirichlet;
    } 
    else if (myRank == np-1 )
    {
        boundaryCondition(0,0) = parallel_ghost;
        boundaryCondition(1,0) = dirichlet;
    }  
    else 
    {
        boundaryCondition(0,0) = parallel_ghost;
        boundaryCondition(1,0) = parallel_ghost;
    }

    const char solutionName[] = "trueDD";

    #define UTRUE(x,t) sin(kxPi*(x))*exp( -kappaPiSq*(t) )
    #define UTRUEX(x,t) kxPi*cos(kxPi*(x))*exp( -kappaPiSq*(t) )
    #define FORCE(x,t) (0.)

    #elif SOLUTION == TRIG_NN

    // True solution for Neumann BC’s
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
    if (myRank == 0)
        printf("ERROR: unknown solution");
    abort();
    #endif

    Real *u_p[2]; // two arrays will be used for current and new times
    if (myRank == 0){
        u_p[0] = new Real [nd1];
        u_p[1] = new Real [nd1];
    } else {
        u_p[0] = new Real [nd1_l];
        u_p[1] = new Real [nd1_l];
    }  

    // Macros to define fortran like arrays
    #define uc(i) u_p[cur ][i-nd1a_l]
    #define un(i) u_p[next][i-nd1a_l]

    // initial conditions
    Real t=0.;
    int cur = 0; // "current" solution, index into u_p[]
    int i;
    for( i=nd1a_l; i<=nd1b_l; i++ ){
        uc(i)=UTRUE(x(i),t);
    }

    if( debug>0 )
    {
        fprintf(debugFile, "After initial conditions\n u=[");
        for( int i=nd1a_l; i<=nd1b_l; i++ )
            fprintf(debugFile, "%10.4e, ",uc(i));
        fprintf(debugFile, "]\n");
    }

    // Time-step restriction is kappa*dt/dx^2 < .5
    const Real dx2 = dx*dx;
    Real dt = cfl*.5*dx2/kappa; // dt, adjusted below
    const int numSteps = ceil(tFinal/dt);
    dt = tFinal/numSteps; // adjust dt to reach the final time
    const Real rx = kappa*dt/dx2;

    if (myRank == 0){
        printf("------------------- Solve the heat equation in 1D solution=%s --------------------- \n",
        solutionName);
        printf("  numGhost=%d, n1a=%d, n1b=%d, nd1a=%d, nd1b=%d\n",numGhost,n1a,n1b,nd1a,nd1b);
        printf("  numSteps=%d, Nx=%d, kappa=%g, tFinal=%g, boundaryCondition(0,0)=%d, boundaryCondition(1,0)=%d\n",
                numSteps,Nx,kappa,tFinal,boundaryCondition(0,0),boundaryCondition(1,0));
    }

    if (debug>0)
    {
        fprintf(debugFile, "------------------- Solve the heat equation in 1D solution=%s --------------------- \n",
        solutionName);
        fprintf(debugFile, "  numGhost=%d, n1a=%d, n1b=%d, nd1a=%d, nd1b=%d\n",numGhost,n1a,n1b,nd1a,nd1b);
        fprintf(debugFile, "  numSteps=%d, Nx=%d, kappa=%g, tFinal=%g, boundaryCondition(0,0)=%d, boundaryCondition(1,0)=%d\n",
                numSteps,Nx,kappa,tFinal,boundaryCondition(0,0),boundaryCondition(1,0));
        fprintf(debugFile, "myRank= %d: nx= %d, nx_l= %d, [n1a_l,n1b_l]=[    %d,    %d], [nd1a_l,nd1b_l]=[    %d,    %d] \n\n", 
            myRank, Nx, nx_l, n1a_l, n1b_l, nd1a_l, nd1b_l);
        
        char variablename[80];
        sprintf(variablename, "x[%d:%d]", nd1a_l, nd1b_l );
        writeMatlabVector( debugFile, x_p, variablename, nd1a_l, nd1b_l );

        sprintf(variablename, "t = %12.4e, u[%d:%d]", 0., nd1a_l, nd1b_l );
        writeMatlabVector( debugFile, u_p[cur], variablename, nd1a_l, nd1b_l );
    }
    


    // ---------- TIME-STEPPING LOOP ---------
    Real cpu0 = getCPU();
    for( int n=0; n<numSteps; n++ )
    {
        t = n*dt; // current time

        const int cur = n % 2; // current time level
        const int next = (n+1) % 2; // next time level

        // --- update the interior points ----
        int i;
        for(i=n1a_l; i<=n1b_l; i++ )
        {
            un(i) = uc(i) + rx*( uc(i+1) - 2.*uc(i) + uc(i-1) ) + dt*FORCE( x(i),t );
        }

        // ---- boundary conditions ----
        for( int side=0; side<=1; side++ )
        {
            const int i = side==0 ? n1a_l : n1b_l; // boundary index
            const int is = 1 - 2*side; // is = 1 on left, -1 on right
            if( boundaryCondition(side,0)==dirichlet )
            {
                un(i) = UTRUE(x(i),t+dt);
                un(i-is) = 3.*un(i) - 3.*un(i+is) + un(i+2*is); // extrapolate ghost
            }
            else if (boundaryCondition(side,0)==neumann)
            {
                // Neumann BC
                un(i-is) = un(i+is) - 2.*is*dx*UTRUEX(x(i),t+dt);
            }

            // handelling the parallel ghost points
            if (np != 1)
            {
                if (side == 0) // left side boundary
                {
                    // first ask all the even ranks to send the value at n1b_l to odd ranks
                    if ((myRank%2 ==0) && (myRank != (np-1)))
                    {   // even ranks are sending the right side value
                        int send_tag = 100 * myRank + 1; // send tag for the right side point 
                        MPI_Send(&uc(n1b_l), 1, MPI_DOUBLE, myRank+1, send_tag, MPI_COMM_WORLD);
                    } 
                    else 
                    {   // odd ranks are receiving value from previus rank for the left side evaluation.
                        int recv_tag = 100 * (myRank-1) + 1;
                        MPI_Status status;
                        Real uc_nd1a_l = 0.;
                        MPI_Recv(&uc_nd1a_l, 1, MPI_DOUBLE, myRank-1, recv_tag, MPI_COMM_WORLD, &status);
                        uc(i-1) = uc_nd1a_l;
                        un(i) = uc(i) + rx*( uc(i+1) - 2.*uc(i) + uc(i-1) ) + dt*FORCE( x(i),t );
                    }
                    // MPI_Barrier(MPI_COMM_WORLD);
                    
                    // now, odd rank will send the value at n1b_l and previous even ranks will receieve it.
                    if (myRank%2 == 1)
                    {   // odd ranks will send the data
                        int send_tag = 100 * myRank + 1; // send tag for the right side point 
                        if (myRank != (np-1)) // last rank cannot send the value none can receive it.
                        {
                            MPI_Send(&uc(n1b_l), 1, MPI_DOUBLE, myRank+1, send_tag, MPI_COMM_WORLD);
                        }
                    }
                    else
                    {   // even ranks will recieve.
                        int recv_tag = 100 * (myRank-1) + 1;
                        MPI_Status status;
                        Real uc_nd1a_l = 0.;
                        if (myRank != 0) // rank-0 must not evaluate the left side bc as it will not parallel-ghost-bs
                        {
                            MPI_Recv(&uc_nd1a_l, 1, MPI_DOUBLE, myRank-1, recv_tag, MPI_COMM_WORLD, &status);
                            uc(i-1) = uc_nd1a_l;
                            un(i) = uc(i) + rx*( uc(i+1) - 2.*uc(i) + uc(i-1)) + dt*FORCE( x(i),t );
                        }
                    }
                    
                    // MPI_Barrier(MPI_COMM_WORLD);
                }
                else if (side == 1) // right side boundary
                {
                    // first odd ranks will send the value at n1a_l for even rank's right side value evaluation
                    if (myRank%2 == 1)
                    {   // odd rank will send the right side value
                        int send_tag = 100 * myRank + 0; // send tag for the left side value
                        MPI_Send(&uc(n1a_l), 1, MPI_DOUBLE, myRank-1, send_tag, MPI_COMM_WORLD);
                    }
                    else 
                    {
                        int recv_tag = 100 * (myRank+1) + 0;
                        MPI_Status status;
                        Real uc_nd1b_l = 0.;
                        if (myRank != (np-1)) // np-1 th rank should not evalutate the right side bc as parallel ghost bc
                        {
                            MPI_Recv(&uc_nd1b_l, 1, MPI_DOUBLE, myRank+1, recv_tag, MPI_COMM_WORLD, &status);
                            uc(i+1) = uc_nd1b_l;
                            un(i) = uc(i) + rx*( uc(i+1) - 2.*uc(i) + uc(i-1)) + dt*FORCE( x(i),t );
                        }
                    }
                    // MPI_Barrier(MPI_COMM_WORLD);

                    // now, even ranks will send the left side value to the previous (odd) ranks
                    if (myRank%2 == 0)
                    {
                        // odd rank will send the right side value
                        int send_tag = 100 * myRank + 0; // send tag for the left side value
                        MPI_Send(&uc(n1a_l), 1, MPI_DOUBLE, myRank-1, send_tag, MPI_COMM_WORLD);
                    }
                    else
                    {
                        int recv_tag = 100 * (myRank+1) + 0;
                        MPI_Status status;
                        Real uc_nd1b_l = 0.;
                        if (myRank != (np-1)) // np-1 th rank should not evalutate the right side bc as parallel ghost bc
                        {
                            MPI_Recv(&uc_nd1b_l, 1, MPI_DOUBLE, myRank+1, recv_tag, MPI_COMM_WORLD, &status);
                            uc(i+1) = uc_nd1b_l;
                            un(i) = uc(i) + rx*( uc(i+1) - 2.*uc(i) + uc(i-1)) + dt*FORCE( x(i),t );
                        }
                    }
                    // MPI_Barrier(MPI_COMM_WORLD);
                }
            }
        }

        // if( debug>1 )
        // {
        //     fprintf(debugFile, "step %d: After update interior and real BCs\n u=[",n+1);
        //     for( int i=nd1a_l; i<=nd1b_l; i++ )
        //         fprintf(debugFile, "%12.4e, ",un(i));
        //     fprintf(debugFile, "]\n");
        // }

        // if( debug>0 )
        // {
        //     // compute the error
        //     Real maxErr=0.;
        //     for( int i=nd1a_l; i<=nd1b_l; i++ )
        //     {
        //         Real err = fabs( un(i) - UTRUE(x(i),t+dt) );
        //         maxErr = max( maxErr,err );
        //     }
        //     fprintf(debugFile, "step=%d, t=%9.3e, maxErr=%9.2e\n",n+1,t+dt,maxErr);
        // }

    } // end time-stepping loop

    Real cpuTimeStep = getCPU()-cpu0;

    // ---- check the error -----
    t +=dt; // tFinal;
    if( (fabs(t-tFinal) > 1e-3*dt/tFinal) && myRank == 0 )
    {
        printf("ERROR: AFTER TIME_STEPPING: t=%16.8e IS NOT EQUAL to tFinal=%16.8e\n",t,tFinal);
    }

    Real *error_p= NULL;
    if (myRank == 0)
        error_p = new Real [nd1];
    else 
        error_p = new Real [nd1_l];
    #define error(i) error_p[i-nd1a_l]

    cur = numSteps % 2;
    Real maxErr=0.;
    Real maxErr_local = 0.;
    // i is already declared
    int itr_start = n1a_l, itr_end = n1b_l;
    if (myRank == 0)    itr_start = nd1a_l;
    if (myRank == np-1) itr_end = nd1b_l;

    for( i=itr_start; i<=itr_end; i++ )
    {
        error(i) = uc(i) - UTRUE(x(i),t);
        maxErr_local = max( maxErr_local, abs(error(i)) );
    }
    MPI_Reduce(&maxErr_local, &maxErr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


    if (debug > 0){
            char variablename[80];
            sprintf(variablename, "t = %12.4e, u[%d:%d]", t, nd1a_l, nd1b_l );
            writeMatlabVector( debugFile, u_p[cur], variablename, nd1a_l, nd1b_l );

            sprintf(variablename, "t = %12.4e, error[%d:%d]", t, nd1a_l, nd1b_l );
            writeMatlabVector( debugFile, error_p, variablename, nd1a_l, nd1b_l );
    }

    // send data from all other ranks except 0th
    if (np != 1)
    {
        if (myRank != 0){
            int buff_size = nd1_l - 1 - 1;
            if (myRank == (np-1))
                buff_size = nd1_l-1;
            MPI_Send(&buff_size, 1, MPI_INT, 0, myRank, MPI_COMM_WORLD);
            MPI_Send(x_p + 1, buff_size, MPI_DOUBLE, 0, myRank, MPI_COMM_WORLD);
            MPI_Send(error_p + 1, buff_size, MPI_DOUBLE, 0, myRank, MPI_COMM_WORLD);
            MPI_Send(u_p[cur] + 1, buff_size, MPI_DOUBLE, 0, myRank, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Status status;
            int offset = nd1_l-1;
            int buff_size = 0;
            for (int rr = 1; rr < np; rr++){
                MPI_Recv(&buff_size, 1, MPI_INT, rr, rr, MPI_COMM_WORLD, &status);
                MPI_Recv(x_p + offset, buff_size, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
                MPI_Recv(error_p + offset, buff_size, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
                MPI_Recv(u_p[cur] + offset, buff_size, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
                offset += buff_size;
            }
        }
    }

    if (myRank == 0){
        printf("numSteps=%4d, Nx=%3d, maxErr=%9.2e, cpu=%9.2e(s) \n",numSteps,Nx,maxErr,cpuTimeStep);

        // --- Write a file for plotting in matlab ---
        FILE *matlabFile = fopen(matlabFileName.c_str(),"w");
        fprintf(matlabFile,"%% File written by heat1d.C\n");
        fprintf(matlabFile,"xa=%g; xb=%g; kappa=%g; t=%g; maxErr=%10.3e; cpuTimeStep=%10.3e;\n",xa,xb, kappa,tFinal,maxErr,cpuTimeStep);
        fprintf(matlabFile,"Nx=%d; dx=%14.6e; numGhost=%d; n1a=%d; n1b=%d; nd1a=%d; nd1b=%d;\n",Nx,dx, numGhost,n1a,n1b,nd1a,nd1b);
        fprintf(matlabFile,"solutionName=\'%s\';\n",solutionName);
        
        writeMatlabVector( matlabFile, x_p, "x", nd1a, nd1b );
        writeMatlabVector( matlabFile, u_p[cur], "u", nd1a, nd1b );
        writeMatlabVector( matlabFile, error_p, "err", nd1a, nd1b );

        fclose(matlabFile);
        printf("Wrote file %s\n\n",matlabFileName.c_str());
    }

    delete [] u_p[0];
    delete [] u_p[1];
    delete [] x_p;
    delete [] error_p;
    delete [] boundaryCondition_p;

    MPI_Finalize();
    return 0;
}