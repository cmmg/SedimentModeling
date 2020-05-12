//problem geometry, mesh control
#define DIMS 5
#define FEOrder 2
#define problemWidth 50.0
#define globalRefinementFactor 6
#define maxRefinementLevel (globalRefinementFactor+2)
#define minRefinementLevel (globalRefinementFactor-2)

//time step controls
#define TimeStep 50.0
#define TotalTime 1000*TimeStep
#define PSTEPS 50 

//output controls
#define outputFileName "solution"

//parameters 

#define M_c 1.0
#define M_eta 0.00025
#define delt 1.0
#define gamma0 1.0

//For Crystallographic Facet Anisortropy
#define www {50.0, 50.0, 50.0}  //w i
#define aaalpha {0.6, 0.6, 0.6}  //alpha i
#define n_orient 3
#define ooorient  {{-0.866025403784,-0.5},{0.866025403784,-0.5},{0.0,1.0}}


//For the four fold anisotropy
#define mm 4.0
#define em 0.3
#define angle0 0.0

//Pressure parameters
#define PressureMin 0.001
#define PressureMax 0.001  //can be changed to 0.01 and make TimeStep 1.0 to see Pressure as a function of time.

//mechanics properties                                                                                                                     
#define elasticModulus (1.0)*(pow(10.0,0.0))
#define PoissonsRatio 0.3

