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
#define www {50.0, 50.0, 50.0}  //w i
#define aaalpha {0.6, 0.6, 0.6}  //alpha i
#define n_orient 3
#define ooorient  {{-0.866025403784,-0.5},{0.866025403784,-0.5},{0.0,1.0}}

#define mm 4.0
#define em 0.125

