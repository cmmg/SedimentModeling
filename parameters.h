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
#define TimeStep 1.0e-2

//Final time
#define TotalTime 174*TimeStep

//Write solution file at PSTEPS interval 
#define PSTEPS 1  


//output controls
#define outputFileName "solution"

//parameters
#define ALPHA (2.5) // 1/mb
#define betaP 2.272e-03
#define ETA 1.0


//moving height velocity
#define Vel 1.53

//domain size
#define problemLength 2.77
