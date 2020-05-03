//problem geometry, mesh control
#define DIMS 2
#define FEOrder 2
#define problemWidth 50.0
#define globalRefinementFactor 5
#define maxRefinementLevel (globalRefinementFactor+1)
#define minRefinementLevel (globalRefinementFactor-1)

//time step controls
#define TimeStep 50.0
#define TotalTime 35000*TimeStep
#define PSTEPS 50 

//output controls
#define outputFileName "solution"

//parameters 

#define M_c 1.0
#define M_eta 0.00025
#define delt 1.0
#define gamma0 1.0
#define mm 4.0
#define em 0.125
#define angle0 0.0
