#ifndef HISTORYCLASS_H_
#define HISTORYCLASS_H_



template <int dim>
struct historyVariables{
public:
  historyVariables<dim>():
  integralOld(0.0), integralOldCells(kcells){}

  //using std:::map to store time history variables                                                                                                                                                       
  double integralOld;
  //double eqvstress, equvStrain;                                                                                                                                                                         
  //double elasStrain11, elasStrain22;                                                                                                                                                                    
  //dealii::Table<2, double > beta, betaIteration;                                                                                                                                                        
  dealii::Table<1, double > integralOldCells;
  };

#endif /* HISTORYCLASS_H_ */
