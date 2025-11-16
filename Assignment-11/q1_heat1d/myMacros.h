#ifndef MYMACROS_H
#define MYMACROS_H

#define Real double

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

#endif