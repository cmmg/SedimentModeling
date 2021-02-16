//problem geometry, mesh control
// DIMS= 5 variables
// FEOrder = Order of basis function
#define DIMS 5
#define FEOrder 2

//No. of grid points = 2^(globalRefinementfactor)*XsubRf*YSubRf
#define globalRefinementFactor 0
#define XSubRf 180 //60
#define YSubRf 30 //60

//time step controls
//time step controls
//dt 
#define TimeStep 1.0e-3

//Final time
#define TotalTime 1704*TimeStep

//Write solution file at PSTEPS interval 
#define PSTEPS 10


//output controls
#define outputFileName "solution"

//parameters
#define ALPHA (2.5) //1/0.975
#define betaP 2.272e-03
#define ETA 1.0


//moving height
#define Vel 1.53
#define problemLength 2.77
#define problemHeight 0.5
