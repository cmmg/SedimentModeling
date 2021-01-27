//problem geometry, mesh control
#define DIMS 3
#define FEOrder 2

#define globalRefinementFactor 7
#define maxRefinementLevel (globalRefinementFactor+2)
#define minRefinementLevel (globalRefinementFactor-2)

//time step controls
#define TimeStep 1.0e-7
#define TotalTime 1500*TimeStep
#define PSTEPS 1 


//output controls
#define outputFileName "solution"

//parameters
#define ALPHA 0.99
#define betaP 2.272e-03
#define ETA 1.0


//moving height
#define V0 1.53
#define problemLength 1.0
