//problem geometry, mesh control
// DIMS= 4 variables
// FEOrder = Order of basis function
#define DIMS 3
#define FEOrder 2

//No. of grid points = 2^(globalRefinementfactor)*XsubRf
#define globalRefinementFactor 0

//No of grid points
#define XSubRf 1000 //60

//time step controls
//dt 
#define TimeStep 1.0e-6

//Final time
#define TotalTime (1.0e+06)*TimeStep

//Write solution file at PSTEPS interval 
#define PSTEPS 1000


//output controls
#define outputFileName "solution"

//parameters
#define ALPHA (1.0/0.1938) // 1/mb
#define aa (5.13e-02)
#define bb (2.48e-05)
#define	cc (1.26e+05)
#define dd (1.45e-03)



//moving height velocity
#define Vel 0

//domain size
#define problemLength (1.0)



