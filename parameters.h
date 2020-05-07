//problem geometry, mesh control
#define DIMS 5
#define FEOrder 2
#define problemWidth 50.0
#define globalRefinementFactor 6
#define maxRefinementLevel (globalRefinementFactor+1)
#define minRefinementLevel (globalRefinementFactor-1)

//time step controls
#define TimeStep 1
#define TotalTime 1000*TimeStep
#define PSTEPS 10

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


#define Pressure 0.001


//mechanics properties                                                                                                                      
#define elasticModulus (1.0)*(pow(10.0,0.0))
#define PoissonsRatio 0.3

