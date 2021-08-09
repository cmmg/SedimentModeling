//problem geometry, mesh control
// DIMS= 4 variables
// FEOrder = Order of basis function
#define DIMS 4
#define FEOrder 2

//No. of grid points = 2^(globalRefinementfactor)*XsubRf
#define globalRefinementFactor 0

//No of grid points
#define XSubRf 500 //60

//time step controls
//dt 
#define TimeStep 1.0e-3

//Final time
#define TotalTime 1003*TimeStep

//Write solution file at PSTEPS interval 
#define PSTEPS 1


//output controls
#define outputFileName "solution"

//parameters
#define ALPHA (1.0/0.1) // 1/mb
#define lam 4.905
#define lam_m 0.0408
#define lam_v 0.408


//moving height velocity
#define Vel 1.0

//domain size
#define problemLength 1.0
