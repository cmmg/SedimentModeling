//problem geometry, mesh control
#define DIMS 5
#define FEOrder 2
#define problemWidth 2.77
#define problemHeight 0.5
//#define refinementFactor 5
#define globalRefinementFactor 1
#define globalRefinementFactor2 2
#define maxRefinementLevel (globalRefinementFactor+1)
#define minRefinementLevel (globalRefinementFactor)

#define sub_x 120
#define sub_y 12


//time step controls
#define TimeStep 1.0e-03
#define TotalTime 1704*TimeStep
#define PSTEPS 10

//output controls
#define outputFileName "solution"

//parameters
#define ALPHA (2.5) //1/0.975
#define betaP 2.272e-03
#define ETA 1.0
#define Vel 1.53
