//problem geometry, mesh control
// DIMS= 4 variables
// FEOrder = Order of basis function
#define DIMS 3
#define FEOrder 2

//No. of grid points = 2^(globalRefinementfactor)*XsubRf
#define globalRefinementFactor 0

//No of grid points
#define XSubRf 10 //60

//time step controls
//dt 
#define TimeStep 1.0e-8
#define timeFactor 0.0005
#define maxdt 1.0e-4

//Series of time step size
#define TimeStep_0 5.0e-5
#define TimeStep_1 5.0e-5
#define TimeStep_2 5.0e-5
#define TimeStep_3 5.0e-5

//Series of time ceiling
#define ceil_0 1.16e-3
#define ceil_1 1.16e-2
#define ceil_2 1.16e-1
#define ceil_3 1.16e-0

#define RampUp (0.00069444444) //Time at which Load is 15MPa

//Final time
#define TotalTime (1.05) // (Slightly more than the non-dimenisonal time 1)
#define tFactor 0.1*TotalTime

//#define TotalTime (1000*TimeStep) // (5710*TimeStep)

//Write solution file at PSTEPS interval 
#define PSTEPS 10
#define PSTEPS2 10

//Number of N-R interations
#define NR_ITR 5  //increase when residual norm is not going down

//output controls
#define outputFileName "solution"


//parameters
#define ALPHA (1.0/0.1545)// 1/mb
#define AA (3.7034e+02)    //
#define BB (2.1073e-05)    //
#define CC0 (2.318e+06)
#define DD0 (2.0029e+11)
#define FF0 (1.7305e+16)

//domain size
#define problemLength (1.0) 


//Kelvin Cells
#define kcells 4  //no. of kelvin cells + 1
#define EE { 0.8*0.1545*3.38e+9, 0.8*0.1545*1.32e+9, 0.8*0.1545*9.15e+9, 0.8*0.1545*6.43e+9} //E0,E1,E2,.......
#define ttau {0, 9.86e+1, 2.51e+3, 6.22e+04} //0,tau1,tau2,tau3, ...
#define tRatio {0, 876.71, 34.37, 1.39}// 0, 86400/tau1, 86400/tau2, ......
#define method 2

