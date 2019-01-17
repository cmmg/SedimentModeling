//problem geometry, mesh control
#define DIMS 2
#define problemWidth 1.0
#define refinementFactor 6

//phase field properties
#define InterfaceEnergyParameter {1.0e-2, 1.0e-2, 1.0e-2} //{Kx, Ky, Kz}
#define dFdC  400*c[q]*(c[q]-1.0)*(c[q]-0.5) //derivative of the free energy
#define Mobility 1.0

//time step controls
#define TimeStep 5.0e-2
#define TotalTime 100*TimeStep

//output controls
#define outputFileName "solution"

//solidifcation parameters
#define D 1.0
#define lam 1.5957
#define newdFdC 1*(phi[q]- lam*c[q] + lam*c[q]*phi[q]*phi[q])*(1-phi[q]*phi[q])
#define W0 1.0
#define epm 0.05
#define mm 4.0
#define theta0 0.125
#define tau0 1000.0


