//problem geometry, mesh control
// DIMS= 4 variables
// FEOrder = Order of basis function
#define DIMS 4
#define FEOrder 2

//No. of grid points = 2^(globalRefinementfactor)*XsubRf
#define globalRefinementFactor 0

//No of grid points
#define XSubRf 1000 //60

//time step controls
//dt 
#define TimeStep 1.0e-3

//Final time
#define TotalTime 5710*TimeStep

//Write solution file at PSTEPS interval 
#define PSTEPS 1


//output controls
#define outputFileName "solution"

//parameters
#define ALPHA (1.0/0.103) // 1/mb
#define betaP 7.14e-03
#define ETA 1.0


//moving height velocity
#define Vel 0.16

//domain size
#define problemLength 0.9302
