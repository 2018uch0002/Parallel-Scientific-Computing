// =====================================================================
//
// HEAT EQUATION IN TWO DIMENSIONS
// Solve with A++ arrays
//
// =====================================================================
#include "A++.h"

#include <mpi.h>

// define some types
typedef double Real;
typedef doubleSerialArray RealArray;
typedef intSerialArray IntegerArray;

#include <float.h>
#include <limits.h>
#define REAL_EPSILON DBL_EPSILON
#define REAL_MIN DBL_MIN

// getCPU() : Return the current wall-clock time in seconds
#include "getCPU.h"

// include commands tp parse command line arguments
#include "parseCommand.h"

// function to write an array to a matlab reabable file:
#include "writeMatlabArray.h"

// get the local-Index bounds based on the ranks
#include "getLocalIndexBounds.h"

// --- Declare the fortran routine as a "C" function
// Some compilers add an under-score to the name of "C" and Fortran routines
// Note also that the fortran name is all lowercase
#define heat2dUpdate heat2dupdate_

extern "C"
{
  void heat2dUpdate( const int & n1a, const int & n1b, const int & n2a, const int & n2b,
                     const int & nd1a, const int & nd1b, const int & nd2a, const int & nd2b,
                     Real & un, const Real &u, const Real & rx, const Real &ry );
}

int
main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  int myRank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get the rank of the processor

  int np;
  MPI_Comm_size(MPI_COMM_WORLD, &np); // get the total number of the processor

  const Real pi = 4.*atan2(1.,1.);

  ios::sync_with_stdio(); // Synchronize C++ and C I/O subsystems
  Index::setBoundsCheck(on); // Turn on A++ array bounds checking

  if (myRank == 0)
    printf("Usage: heat2d -nx=<i> -option=[0|1|2|3] -tFinal=<f> -debug=<i> -saveMatlab=[0|1|2] -matlabFile=<s>\n   option : 0=scalarIndexing, 1=arrayIndexing, 2=cIndexing, 3=fortranRoutine\n ");

  enum BoundaryConditionsEnum
  {
    periodic=-1,
    dirichlet=1,
    neumann=2,
    parallel_ghost = 3
  };

  enum OptionsEnum
  {
    scalarIndexing=0,
    arrayIndexing =1,
    cIndexing =2,
    fortranRoutine=3
  };

  int option=scalarIndexing;

  const int numberOfDimensions=2;

  int debug = 0;
  Real kappa=.1;
  Real xa=0., xb=1.; // domain is [xa,xb] X [ya,yb]
  Real ya=0., yb=1.;

  Real tFinal=.5;
  int nx=100, ny=nx;

  int saveMatlab=0; // 1 = save a matlab file, 2=save solution too
  string matlabFileName = "heat2d.m";

  string line;
  for( int i=1; i<argc; i++ )
  {
    line=argv[i];
    // printf("Input: argv[%d] = [%s]\n",i,line.c_str());
    if( parseCommand( line,"-nx=",nx) )
    {
      ny=nx;
    }
    else if( parseCommand( line,"-debug=",debug,myRank==0) ){}
    else if( parseCommand( line,"-option=",option,myRank==0) ){}
    else if( parseCommand( line, "-tFinal=",tFinal,myRank==0) ){}
    else if( parseCommand( line,"-saveMatlab=",saveMatlab,myRank==0) ){}
    else if( parseCommand( line,"-matlabFileName=",matlabFileName,myRank==0) ){}

  }

  if (option != 2 && np > 1)
  {
    printf("ERROR: only CIndexing option can only be run on multiple MPI-Ranks.\n");
    abort();
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
        fprintf(file," np=%d, myRank=%d, option=%d\n",np,myRank, option);
        }
    }
  }

  const int numGhost=1;
  const int n1a = 0;
  const int n1b = n1a + nx;
  const int nd1a = n1a-numGhost;
  const int nd1b = n1b+numGhost;
  const int nd1 = nd1b-nd1a+1;

  const int n2a = 0;
  const int n2b = n2a + ny;
  const int nd2a = n2a-numGhost;
  const int nd2b = n2b+numGhost;
  const int nd2 = nd2b-nd2a+1;

  // get the local bounds only for y-directions
  int ny_l = 0, n2a_l = 0, n2b_l = 0;
  getLocalIndexBounds(myRank, np, ny, ny_l, n2a_l, n2b_l);
  const int nd2a_l = n2a_l - numGhost;
  const int nd2b_l = n2b_l + numGhost;
  const int nd2_l = nd2b_l - nd2a_l + 1; // total number of grid points per processor

  IntegerArray gridIndexRange(2,numberOfDimensions);
  IntegerArray dimension(2,numberOfDimensions);
  IntegerArray boundaryCondition(2,numberOfDimensions);

  gridIndexRange(0,0)=n1a; gridIndexRange(1,0)=n1b;
  gridIndexRange(0,1)=n2a_l; gridIndexRange(1,1)=n2b_l;

  // based on the code structure it is not necessary to have n2a_l and n2b_l
  dimension(0,0)=nd1a; dimension(1,0)=nd1b;
  if (myRank == 0)
  {
    dimension(0,1)=nd2a; dimension(1,1)=nd2b;
  }
  else 
  {  
    dimension(0,1)=nd2a_l; dimension(1,1)=nd2b_l;
  }

  boundaryCondition(0,0)=dirichlet; // left
  boundaryCondition(1,0)=dirichlet; // right
  
  if (myRank == 0)
  {
    boundaryCondition(0,1)=dirichlet; // bottom
    if (np==1)
      boundaryCondition(1,1)=dirichlet; // top
    else 
      boundaryCondition(1,1)=parallel_ghost; // top
     
  }
  else
  {
    boundaryCondition(0,1)=parallel_ghost; // bottom
    if (myRank != (np-1))
      boundaryCondition(1,1)= parallel_ghost; // top
    else 
      boundaryCondition(1,1)=dirichlet; // top
  }
  

  // Grid points
  Range Rx(nd1a,nd1b), Ry;
  if (myRank == 0)
    Ry = Range(nd2a,nd2b);
  else
    Ry = Range(nd2a_l,nd2b_l);
  
  RealArray x(Rx,Ry,2);

  Real dx[2];
  dx[0] = (xb-xa)/nx;
  dx[1] = (yb-ya)/ny;

  int i1,i2;
  for( i2=nd2a_l; i2<=nd2b_l; i2++ )
    for( i1=nd1a; i1<=nd1b; i1++ )
    {
      x(i1,i2,0) = xa + (i1-n1a)*dx[0];
      x(i1,i2,1) = ya + (i2-n2a)*dx[1];
    }

  #define TRIG_DD 1
  #define POLY_DD 2
  #ifndef SOLUTION
  // #define SOLUTION TRIG_DD
  #define SOLUTION POLY_DD
  #endif

  const Real kx=2., ky=3;
  const Real kxp=kx*pi;
  const Real kyp=ky*pi;

  #if SOLUTION == TRIG_DD
  
  #define UTRUE(x,y,t) sin( kxp*(x) )*sin( kyp*(y) )*exp(-kappa*(kxp*kxp+kyp*kyp)*(t) )
  #define FORCE(x,y,t) ( 0. )
  
  #elif SOLUTION == POLY_DD
  // polynomial manufactured solution
  static const Real c0=.2, c1=.1, c2=.3;
  static const Real b0=1., b1=.5, b2=.25;
  static const Real a0=1., a1=.3, a2=0.;
  #define UTRUE(x,y,t)   (b0 + (x)*( b1 + (x)*b2  ))*(c0 + (y)*( c1 + (y)*c2  ))*( a0 + (t)*( a1 + (t)*a2 ) )
  #define UTRUET(x,y,t)  (b0 + (x)*( b1 +(x)*b2   ))*(c0 + (y)*( c1 + (y)*c2  ))*(         a1 + 2.*(t)*a2 )
  #define UTRUEXX(x,y,t)                   ( 2.*b2 )*(c0 + (y)*( c1 + (y)*c2  ))*( a0 + (t)*( a1 + (t)*a2 ) )
  #define UTRUEYY(x,y,t) (b0 + (x)*( b1 + (x)*b2  ))*(     2.*c2               )*( a0 + (t)*( a1 + (t)*a2 ) )
  
  #define FORCE(x,y,t) ( UTRUET(x,y,t) - kappa*( UTRUEXX(x,y,t)+UTRUEYY(x,y,t) ) )

  #endif

  string optionName = option==scalarIndexing ? "scalarIndexing" :
  option==arrayIndexing ? "arrayIndexing " :
  option==cIndexing ? "cIndexing     " :
  option==fortranRoutine ? "fortranRoutine" : "unknown";

  if (myRank == 0)
  {
    printf("----- Solve the Heat Equation in two dimensions ------\n");
    printf("      option=%d : %s \t",option, optionName.c_str());
    if (SOLUTION==TRIG_DD)
      printf("      True-Dirchlet-Dirchlet\n");
    else if (SOLUTION == POLY_DD)
      printf("      Polynomial-Manufactured\n");     
    printf("      saveMatlab=%d, matlabFileName=%s \n",saveMatlab,matlabFileName.c_str());
    printf("   kappa=%.3g, nx=%d, ny=%d, tFinal=%6.2f, kx=%g, ky=%g\n",kappa,nx,ny,tFinal,kx,ky);
  }
  if (debug>0){
    fprintf(debugFile, "----- Solve the Heat Equation in two dimensions ------\n");
    fprintf(debugFile, "      option=%d : %s \t",option, optionName.c_str());
    if (SOLUTION==TRIG_DD)
      fprintf(debugFile, "      True-Dirchlet-Dirchlet\n");
    else if (SOLUTION == POLY_DD)
      fprintf(debugFile, "      Polynomial-Manufactured\n");
    fprintf(debugFile, "      saveMatlab=%d, matlabFileName=%s \n",saveMatlab,matlabFileName.c_str());
    fprintf(debugFile, "   kappa=%.3g, nx=%d, ny=%d, tFinal=%6.2f, kx=%g, ky=%g\n",kappa,nx,ny,tFinal,kx,ky);
    fprintf(debugFile, "myRank= %d: nx= %d, nx_l= %d, [n1a_l,n1b_l]=[    %d,    %d], [nd1a_l,nd1b_l]=[    %d,    %d] \n\n", 
            myRank, nx, nx, n1a, n1b, nd1a, nd1b);
    fprintf(debugFile, "myRank= %d: ny= %d, ny_l= %d, [n2a_l,n2b_l]=[    %d,    %d], [nd2a_l,nd2b_l]=[    %d,    %d] \n\n", 
            myRank, ny, ny_l, n2a_l, n2b_l, nd2a_l, nd2b_l);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // we store two time levels
  RealArray ua[2];
  ua[0].redim(Rx,Ry); ua[0]=0.;
  ua[1].redim(Rx,Ry); ua[1]=0.;

  // initial conditions
  RealArray & u0 = ua[0];
  Real t=0.;
  for( i2=nd2a_l; i2<=nd2b_l; i2++ )
    for( i1=nd1a; i1<=nd1b; i1++ )
    {
      u0(i1,i2)= UTRUE(x(i1,i2,0),x(i1,i2,1),t);
    }

  // Time-step restriction:
  // Forward Euler: kappa*dt*( 1/dx^2 + 1/dy^2 ) < cfl*.5
  Real cfl=.9;
  Real dt = cfl*(.5/kappa)/( 1./(dx[0]*dx[0]) + 1./(dx[1]*dx[1]) );

  int numSteps=ceil(tFinal/dt);
  dt = tFinal/numSteps; // adjust dt to reach the final time

  Real rx = kappa*(dt/(dx[0]*dx[0]));
  Real ry = kappa*(dt/(dx[1]*dx[1]));
  int cur=0;
  int i,n;
  Index I1=Range(n1a,n1b);
  Index I2=Range(n2a,n2b);

  Real cpu1 = getCPU();
  for( n=0; n<numSteps; n++ )
  {
    t = n*dt; // cur time

    int next = (cur+1) % 2;
    RealArray & u = ua[cur];
    RealArray & un = ua[next];

    if( option==scalarIndexing )
    {
      for( i2=n2a; i2<=n2b; i2++ )
        for( i1=n1a; i1<=n1b; i1++ )
        {
          un(i1,i2) = u(i1,i2) + rx*( u(i1+1,i2) -2.*u(i1,i2) + u(i1-1,i2) )
          + ry*( u(i1,i2+1) -2.*u(i1,i2) + u(i1,i2-1) );
        }
      }
    else if( option==arrayIndexing )
    {
      un(I1,I2) = u(I1,I2) + rx*( u(I1+1,I2) -2.*u(I1,I2) + u(I1-1,I2) )
      + ry*( u(I1,I2+1) -2.*u(I1,I2) + u(I1,I2-1) );
    }
    else if( option==cIndexing )
    {
      // Index as C arrays
      const double *u_p = u.getDataPointer();
      double *un_p = un.getDataPointer();
      #define U(i1,i2) u_p[(i1-nd1a)+nd1*(i2-nd2a_l)]
      #define UN(i1,i2) un_p[(i1-nd1a)+nd1*(i2-nd2a_l)]

      for( i2=n2a_l; i2<=n2b_l; i2++ )
        for( i1=n1a; i1<=n1b; i1++ )
        {
          UN(i1,i2) = U(i1,i2) + rx*( U(i1+1,i2) -2.*U(i1,i2) + U(i1-1,i2) )
          + ry*( U(i1,i2+1) -2.*U(i1,i2) + U(i1,i2-1) ) 
          + dt * FORCE( x(i1, i2, 0), x(i1, i2, 1), t );
        }

    }
    else if( option==fortranRoutine )
    {
      // call a fortran routine
      // Note: pass first element of un and u arrays to Fortran (call by reference)
      if( true )
      {
        const double *u_p = u.getDataPointer();
        double *un_p = un.getDataPointer();
        heat2dUpdate( n1a,n1b,n2a,n2b,
        u.getBase(0),u.getBound(0), u.getBase(1),u.getBound(1), // pass array dimensions
        *un_p, *u_p, rx,ry );
      }
      else
      {
        // this will also work -- pass first element of un and u
        heat2dUpdate( n1a,n1b,n2a,n2b,
        u.getBase(0),u.getBound(0), u.getBase(1),u.getBound(1), // pass array dimensions
        un(nd1a,nd2a), u(nd1a,nd2a), rx,ry );
      }
    }
    else
    {
      printf("ERROR: unknown option=%d\n",option);
      abort();
    }

    // --- parallel ghost boundary conditions
    // axis = 1
    if (np>1)
    {
      for( int side=0; side<=1; side++ ) // only for top and bottom
      {
        //bottom or top
        int is = 1-2*side; // is=1 on bottom, is=-1 on top

        double *u_p = u.getDataPointer();
        double *un_p = un.getDataPointer();
    
        const bool can_recv = ((side==0) && (myRank!=0)) || ((side==1) && (myRank!=(np-1)));
        const bool can_send = ((side==0) && (myRank!=(np-1))) || ((side==1) && (myRank!=0));

        if ( can_recv )
        {
          i2= gridIndexRange(side,1);
          int i2g = i2 - is; // index of ghost point: 1

          const int recv_from_rank = myRank - is;
          const int recv_tag = (myRank - is) * 100 + is;
          MPI_Status status;
          
          MPI_Recv(u_p + nd1*(i2g-nd2a_l), nd1, MPI_DOUBLE, 
                    recv_from_rank, recv_tag, MPI_COMM_WORLD, &status);
        }

        if ( can_send )
        {
          const int send_to_rank = myRank + is;
          const int send_tag = myRank * 100 + is; // 0: bottom 1: top

          const int opposite_side = side*-1 + 1;
          i2= gridIndexRange(opposite_side,1) - is; // this will be index we suppose to send.
          int i2g = i2 + is; // index of ghost point: 1
          
          MPI_Send(u_p+nd1*(i2-nd2a_l), nd1, MPI_DOUBLE, 
                    send_to_rank, send_tag, MPI_COMM_WORLD);
        }

        // update the boundary
        for( i1=nd1a; i1<=nd1b; i1++ )
        {
          un(i1,i2) = u(i1,i2) + rx*( u(i1+1,i2) -2.*u(i1,i2) + u(i1-1,i2) )
                    + ry*( u(i1,i2+1) -2.*u(i1,i2) + u(i1,i2-1) ) 
                    + dt * FORCE( x(i1, i2, 0), x(i1, i2, 1), t );
        }
      }
    }

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
            int i1g = i1 - is; // index of ghost point
            for( i2=nd2a_l; i2<=nd2b_l; i2++ )
            {
              un(i1,i2) = UTRUE(x(i1,i2,0),x(i1,i2,1),t+dt);
              un(i1g,i2) = 3.*un(i1,i2) - 3.*un(i1+is,i2) + un(i1+2*is,i2); // extrap ghost
            }
          }
          else
          { // bottom or top
            i2= gridIndexRange(side,axis);
            int i2g = i2 - is; // index of ghost point
            for( i1=nd1a; i1<=nd1b; i1++ )
            {
              un(i1,i2) = UTRUE(x(i1,i2,0),x(i1,i2,1),t+dt);
              un(i1,i2g) = 3.*un(i1,i2) - 3.*un(i1,i2+is) + un(i1,i2+2*is); // extrap ghost
            }
          }
        }
        else if (boundaryCondition(side,axis) != parallel_ghost )
        {
          printf("ERROR: unknown boundaryCondition=%d\n",boundaryCondition(side,axis));
          abort();
        }
      } // end for axis

    cur = next;
  }

  Real cpuTimeStep = getCPU()-cpu1;

  // --- compute errors ---
  RealArray & uc = ua[cur];
  RealArray err(Rx,Ry);

  Real maxErr=0., maxNorm=0.;
  Real maxErr_local=0., maxNorm_local=0.;
  for( i2=n2a_l; i2<=n2b_l; i2++ )
    for( i1=n1a; i1<=n1b; i1++ )
    {
      err(i1,i2) = fabs(uc(i1,i2) - UTRUE(x(i1,i2,0),x(i1,i2,1),tFinal));
      maxErr_local = max(err(i1,i2),maxErr_local);
      maxNorm_local = max(uc(i1,i2),maxNorm_local);
    }

  if (debug > 1)
  {
    fprintf(debugFile, "t = %12.4e, maxErr =%12.4e, maxNorm= %12.4e", 
              t, maxErr_local, maxNorm_local );
  }

  MPI_Reduce(&maxErr_local, &maxErr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&maxNorm_local, &maxNorm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (myRank==0)
    maxErr /= max(maxNorm,REAL_MIN); // relative error

  if (myRank == 0)
    printf("option=%s: numSteps=%d nx=%d maxNorm=%8.2e maxRelErr=%8.2e cpuTimeStep=%9.2e(s)\n",
    optionName.c_str(),numSteps,nx,maxNorm,maxErr,cpuTimeStep);

  if( nx<=10 )
  {
    uc.display("ua[cur]");
    err.display("err");
  }

  // gather the data on rank-0
  if (np>1)
  {
    if (myRank!=0)
    {
      double* x_p = x.getDataPointer();
      // send the size of the buffer
      int buff_size = (nd2_l - 2*numGhost + (myRank == (np-1) ? numGhost : 0 )) * nd1 ;
      MPI_Send(&buff_size, 1, MPI_INT, 0, myRank, MPI_COMM_WORLD);

      // send x
      // dim-0
      MPI_Send(x_p + nd1,  buff_size, MPI_DOUBLE, 0, myRank, MPI_COMM_WORLD);
      // dim-1
      MPI_Send(x_p + nd1*nd2_l + nd1,  buff_size, MPI_DOUBLE, 0, myRank, MPI_COMM_WORLD);

      // send u
      double *u_p = ua[cur].getDataPointer();
      MPI_Send(u_p + nd1, buff_size, MPI_DOUBLE, 0, myRank, MPI_COMM_WORLD);
      
      // send err
      double* error_p = err.getDataPointer();
      MPI_Send(error_p + nd1, buff_size, MPI_DOUBLE, 0, myRank, MPI_COMM_WORLD);
    }
    else
    {
      MPI_Status status;
      int buff_size = 0;
      int offset = nd1 * (nd2_l-1);
      int xp_offset = nd1 * nd2;
      double* x_p = x.getDataPointer();
      double* u_p = ua[cur].getDataPointer();
      double* error_p = err.getDataPointer();
      for (int r = 1; r < np; r++)
      {
        
        MPI_Recv(&buff_size, 1, MPI_INT, r, r, MPI_COMM_WORLD, &status);

        // recv x
        MPI_Recv(x_p + offset , buff_size, MPI_DOUBLE, r, r, MPI_COMM_WORLD, &status);

        MPI_Recv(x_p + xp_offset + offset , buff_size, MPI_DOUBLE, r, r, MPI_COMM_WORLD, &status);

        // recv u
        MPI_Recv(u_p + offset, buff_size, MPI_DOUBLE, r, r, MPI_COMM_WORLD, &status);

        // recv err
        MPI_Recv(error_p + offset, buff_size, MPI_DOUBLE, r, r, MPI_COMM_WORLD, &status);

        offset += buff_size;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // --- OPTIONALLY write a matlab file for plotting in matlab ---
  if (myRank == 0)
  {
    if( saveMatlab )
    {
      FILE *matlabFile = fopen(matlabFileName.c_str(),"w");
      fprintf(matlabFile,"%% File written by heat2d.C\n");
      fprintf(matlabFile,"xa=%g; xb=%g; ya=%g; yb=%g; kappa=%g; t=%g; maxErr=%10.3e; cpuTimeStep=%10.3e;\n",
      xa,xb,ya,yb,kappa,tFinal,maxErr,cpuTimeStep);

      fprintf(matlabFile,"n1a=%d; n1b=%d; nd1a=%d; nd1b=%d;\n",n1a,n1b,nd1a,nd1b);
      fprintf(matlabFile,"n2a=%d; n2b=%d; nd2a=%d; nd2b=%d;\n",n2a,n2b,nd2a,nd2b);
      fprintf(matlabFile,"dx(1)=%14.6e; dx(2)=%14.6e; numGhost=%d;\n",dx[0],dx[1],numGhost);
      fprintf(matlabFile,"option=%d; optionName=\'%s\';\n",option,optionName.c_str());

      if( saveMatlab>1 )
      {
        writeMatlabArray( matlabFile, x, "x", 2, dimension );
        writeMatlabArray( matlabFile, ua[cur], "u", 1, dimension );
        writeMatlabArray( matlabFile, err, "err", 1, dimension );
      }
        fclose(matlabFile);
        printf("Wrote file [%s]\n",matlabFileName.c_str());
    }
  }

  MPI_Finalize();
  return 0;
}
