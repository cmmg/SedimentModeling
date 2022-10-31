//problem geometry, mesh control
// DIMS= 4 variables
// FEOrder = Order of basis function
#define DIMS 3
#define FEOrder 2

//No. of grid points = 2^(globalRefinementfactor)*XsubRf
#define globalRefinementFactor 0

//No of grid points
#define XSubRf 20 //60

//time step controls
//dt 
#define TimeStep 1.0e-1

//Series of time step size
#define TimeStep_0 1.0e-5
#define TimeStep_1 1.0e-4
#define TimeStep_2 1.0e-3
#define TimeStep_3 1.0e-2

//Series of time ceiling
#define ceil_0 1.0e-2
#define ceil_1 6.0e-2
#define ceil_2 1.6e-1
#define ceil_3 1.16e-0

#define RampUp 6.94e-04 //Time at which Load is 15MPa

//Final time
#define TotalTime (ceil_3) // (Slightly more than the non-dimenisonal time 1)

//#define TotalTime (1000*TimeStep) // (5710*TimeStep)

//Write solution file at PSTEPS interval 
#define PSTEPS 10

//Number of N-R interations
#define NR_ITR 5  //increase when residual norm is not going down

//output controls
#define outputFileName "solution"

//parameters
//#define ALPHA (1.0/0.2)   //(1.0/0.103) // 1/mb
//#define betaP (1.57e-03)    //(7.14e-03)
//#define ETA 1.0


//parameters
#define ALPHA (1.0/0.1938)   //(1.0/0.103) // 1/mb
#define AA (3.52e+00)    //(7.14e-03)
#define BB (2.46e-05)    //(7.14e-03)
#define CC (2.91e-02)    //(7.14e-03)
#define DD (5.02e-03)    //(7.14e-03)



//moving height velocity
//#define Vel 0.16

//domain size
#define problemLength (1.0) //0.9302
