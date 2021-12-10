//problem geometry, mesh control
// DIMS= 4 variables
// FEOrder = Order of basis function
#define DIMS 1
#define FEOrder 1

//No. of grid points = 2^(globalRefinementfactor)*XsubRf
#define globalRefinementFactor 0

//No of grid points
#define XSubRf 100 //60

//time step controls
//dt 
#define TimeStep 1.0e-00

//Final time
#define TotalTime 1802*TimeStep

//Write solution file at PSTEPS interval 
#define PSTEPS 1


//output controls
#define outputFileName "solution"

//parameters
#define Tc 53.0 //1.0

//domain size
#define problemLength 1.0
