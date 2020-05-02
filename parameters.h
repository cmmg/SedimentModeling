//problem geometry, mesh control
#define DIMS 2
#define FEOrder 2
#define problemWidth 50.0
#define globalRefinementFactor 6
#define maxRefinementLevel (globalRefinementFactor+2)
#define minRefinementLevel (globalRefinementFactor-2)

//time step controls
#define TimeStep 0.1
#define TotalTime 35000*TimeStep
#define PSTEPS 50 

//output controls
#define outputFileName "solution"

//parameters 

#define M_c 1.0
#define M_eta 0.00025
#define delt 1.0
#define f_a 0.5*c[q]*c[q]/16.0
#define f_a_c 0.5*c[q]/8.0
#define f_a_c_c 0.5/8.0
#define f_b 0.5*(c[q]-1)*(c[q]-1)/16.0
#define f_b_c 0.5*(c[q]-1)/8.0
#define f_b_c_c 0.5/8.0
#define HH 3.0*eta[q]*eta[q] - 2.0*eta[q]*eta[q]*eta[q]
#define H_eta 6.0*eta[q] - 6.0*eta[q]*eta[q]
#define gamma0 1.0
#define www {50.0, 50.0, 50.0}  //w i
#define aaalpha {0.6, 0.6, 0.6}  //alpha i
#define n_orient 3
#define ooorient  {{-0.866025403784,-0.5},{0.866025403784,-0.5},{0.0,1.0}}


