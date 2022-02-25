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
#define TimeStep 1.0e-3

//Final time
#define TotalTime (1000*TimeStep) // (5710*TimeStep)

//Write solution file at PSTEPS interval 
#define PSTEPS 10


//output controls
#define outputFileName "solution"

//parameters
//#define ALPHA (1.0/0.2)   //(1.0/0.103) // 1/mb
//#define betaP (1.57e-03)    //(7.14e-03)
//#define ETA 1.0


//parameters
#define ALPHA (1.0/0.1938)   //(1.0/0.103) // 1/mb
#define AA (4.95e-01)    //(7.14e-03)
#define BB (2.46e-05)    //(7.14e-03)
#define CC (2.91e-02)    //(7.14e-03)
#define DD (5.02e-03)    //(7.14e-03)



//moving height velocity
//#define Vel 0.16

//domain size
#define problemLength (1.0) //0.9302
