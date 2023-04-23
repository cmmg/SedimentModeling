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
#define TimeStep_0 1.0e-7
#define TimeStep_1 1.0e-7
#define TimeStep_2 1.0e-7
#define TimeStep_3 1.0e-7

//Series of time ceiling
#define ceil_0 1.16e-3
#define ceil_1 1.16e-2
#define ceil_2 1.16e-1
#define ceil_3 1.16e-0

#define RampUp 6.94e-04 //Time at which Load is 15MPa

//Final time
#define TotalTime (ceil_3) // (Slightly more than the non-dimenisonal time 1)

//#define TotalTime (1000*TimeStep) // (5710*TimeStep)

//Write solution file at PSTEPS interval 
#define PSTEPS 10
#define PSTEPS2 1000

//Number of N-R interations
#define NR_ITR 5  //increase when residual norm is not going down

//output controls
#define outputFileName "solution"


//parameters
#define ALPHA (1.0/0.1545)// 1/mb
#define AA (3.7034e+02)    //
#define BB (2.1073e-05)    //
#define CC0 (2.3182e+06)    //
#define DD0 (2.0029e+11)    //
#define FF0 (1.7305e+16)    //

//domain size
#define problemLength (1.0) 


//Kelvin Cells
#define kcells 4  //no. of kelvin cells + 1
#define EE {0.1545*3.4e+9, 0.1545*1.3e+9, 0.1545*9.1e+9, 0.1545*6.4e+9} //E0,E1,E2,.......
#define ttau {0, 1.0e+2, std::pow(10,3.4), std::pow(10.0,4.8)} //0,tau1,tau2,tau3, ...
#define tRatio {0, 8.64e+2, 3.44e+1, 1.37e+0}// 0, 86400/tau1, 86400/tau2, ......
