//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2012
//authors: rudraa (2012, 2018)
//

#ifndef CHEMO_H_
#define CHEMO_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"
#include <math.h>     

/*
template <int dim>
struct ExportVariables{
public:
  ExportVariables<dim>():
  time(0.0), Location(dim), Porosity(dim), Velocity(dim), Eff_Pressure(dim) {}
  //using std:::map to store time history variables                                                                                                                                                                                               
  double time;
  dealii::Table<1, double > Location, Porosity,Velocity,Eff_Pressure;
};
*/

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


template <class T, int dim>
  void Integration(double dt, unsigned int q, double currentTime, std::vector<historyVariables<dim>*>& History, Table<1, double>& Integral ,  double tTauRatio[], std::vector<double>  &Time, std::vector<double>  &Peff){ 

  if (currentTime=dt){
    History[q]->integralOld=0;
    for(unsigned int k=1;k<kcells;k++)
      {
	History[q]->integralOldCells[k]=0.0;
      }
    
  }
  
  //Do the numerical integral  
  //no. of steps
  int last=0,secLast=0;
  //find last and second last time and pressure if they are not negative
  if (Time.size()>1) {
    last = Time.size()-1;
    secLast=Time.size()-2; 
  }
  
  double timeSize= Time[last]-Time[secLast];
  double valIntegral=0;

  /*
  valIntegral=0.5*std::exp(Time[last]*tTauRatio)*Peff[last];
  valIntegral+=0.5*std::exp(Time[secLast]*tTauRatio)*Peff[secLast];
  valIntegral= timeSize*valIntegral;
  Integral[q]= valIntegral +   History[q]->integralOld ;
  */

  for(unsigned int k=1;k<kcells;k++){
    double valIntegral=0;  
    valIntegral=0.5*std::exp(Time[last]*tTauRatio[k])*Peff[last];
    valIntegral+=0.5*std::exp(Time[secLast]*tTauRatio[k])*Peff[secLast];
    valIntegral= timeSize*valIntegral;
    Integral[k]= valIntegral +   History[q]->integralOldCells[k] ;    
  }
  
  //update latest integral to history variables
  History[q]->integralOld=Integral[0];
  for(unsigned int k=1;k<kcells;k++)
    {
      History[q]->integralOldCells[k]=Integral[k];
    }
  
}


double trapezoidal(double dtSize, double tTauRatio, std::vector<double>  &Time, std::vector<double>  &Peff) {
  
  double valueIntegral=0;
  int rowFinal= Time.size();
  
  if (rowFinal==1) {
    valueIntegral+=0;
  }

  else if(rowFinal==2) {
    valueIntegral+= 0.5*(Time.back()-0)*std::exp(Time.back()*tTauRatio)*Peff.back();
  }
  
  else if (rowFinal>2) {    
    for (unsigned int row=1; row< rowFinal; ++row) {                                                                                                                                                                                                                            
      valueIntegral+= 0.5*(Time[row]-Time[row-1])*std::exp(Time[row]*tTauRatio)*Peff[row];
      valueIntegral+= 0.5*(Time[row]-Time[row-1])*std::exp(Time[row-1]*tTauRatio)*Peff[row-1];            
    }    
  }  
  return valueIntegral;
}

double threePointIntegral(std::vector<double>  &x, std::vector<double>  &f) {
  // x, f are 3 element vector
  // x is time
  // f is function at x0,x1,x2
  int n=x.size();
  FullMatrix<double> Mat(n,n), invMat(n,n);
  double coeff[n], func[n], value=0;

  for (unsigned int i=0; i<n; i++){
    Mat(i,0)=x[i]*x[i];
    Mat(i,1)=x[i];
    Mat(i,2)=1;
  }

  Mat.gauss_jordan();
  
  for (unsigned int i=0; i<n; i++){
    coeff[i]=Mat(i,0)*f[0]+Mat(i,1)*f[1]+Mat(i,2)*f[2];
  }

  value+=coeff[2]*(x[2]-x[0]);
  value+=(1.0/2.0)*coeff[1]*(x[2]*x[2]-x[0]*x[0]);
  value+=(1.0/3.0)*coeff[0]*(x[2]*x[2]*x[2]-x[0]*x[0]*x[0]);

  //std::cout <<"value of integral is "<< value<<"\n";
  return value ; 
}

double higherOrderIntegral(double dtSize, double tTauRatio, std::vector<double>  &Time, std::vector<double>  &Peff) {

  double valueIntegral=0;
  int rowFinal= Time.size();
  std::vector<double> xx,ff;
  
  if (rowFinal<=2) {
    valueIntegral+=trapezoidal(dtSize,tTauRatio,Time,Peff);
  }
  
  else if (rowFinal>2) {
    int row=2;
    while(row<rowFinal) {
      xx.push_back(Time[row-2]);
      xx.push_back(Time[row-1]);
      xx.push_back(Time[row]);
      ff.push_back(std::exp(Time[row-2]*tTauRatio)*Peff[row-2]);
      ff.push_back(std::exp(Time[row-1]*tTauRatio)*Peff[row-1]);
      ff.push_back(std::exp(Time[row]*tTauRatio)*Peff[row]);      
      valueIntegral+=threePointIntegral(xx,ff);
      row=row+2;
      xx.clear();
      ff.clear();
    }

    if (row==rowFinal+1) {
      valueIntegral+= 0.5*(Time[row-3]-Time[row-3])*std::exp(Time[row-3]*tTauRatio)*Peff[row-3];
      valueIntegral+= 0.5*(Time[row-2]-Time[row-2])*std::exp(Time[row-2]*tTauRatio)*Peff[row-2];            
    }
    
  }  
     
  return valueIntegral;
}




//Chemistry residual implementation
template <int dim>
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF,  const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R, std::vector<historyVariables<dim>*>& History, Table<2,std::vector<double> > &CellHist, double currentTime, double totalTime) {
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
 
  //evaluate field and gradients
  //dealii::Table<1,Sacado::Fad::DFad<double> > phi(n_q_points), vel(n_q_points) , press(n_q_points), ORDER(n_q_points);
  dealii::Table<1,Sacado::Fad::DFad<double> > phi(n_q_points), vel(n_q_points) , press(n_q_points);
  dealii::Table<2,Sacado::Fad::DFad<double> > phi_j(n_q_points, dim), vel_j(n_q_points,dim), press_j(n_q_points, dim);
  
  dealii::Table<1,double> phi_conv(n_q_points);
  dealii::Table<1,double> press_conv(n_q_points);
  //dealii::Table<1,double> ORDER_conv(n_q_points);
  
  for (unsigned int q=0; q<n_q_points; ++q) {
    phi[q]=0.0; vel[q]=0.0; press[q]=0.0; phi_conv[q]=0.0; press_conv[q]=0.0; 

    for (unsigned int j=0; j<dim; j++) {phi_j[q][j]=0.0; vel_j[q][j]=0.0; press_j[q][j]=0.0;}
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;      
      if (ck==0) { phi[q]+=fe_values.shape_value(i, q)*ULocal[i]; phi_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i]; }
      else if (ck==1){ vel[q]+=fe_values.shape_value(i, q)*ULocal[i]; }
      else if (ck==2) {press[q]+=fe_values.shape_value(i, q)*ULocal[i]; press_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i]; }

      //else if (ck==3) {ORDER[q]+=fe_values.shape_value(i, q)*ULocal[i]; ORDER_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i];}

      for (unsigned int j=0; j<dim; j++) {
	if (ck==0) {
	  phi_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];
	}
	
	else if (ck==1) {
	  vel_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];
	}
	else if (ck==2){
	  press_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];	   
	}	
      }
    }
    
  }
  
  double CC=0;
  double DD=0; 
  double E[] =EE, tR[]=tRatio, tau[]=ttau;


  CC= 1.0/E[0] ; 
  CC=CC0*CC;
  
  for (unsigned int i = 1; i < kcells; i++) {
    DD+= 1.0/E[i]/tau[i];
  }
  DD=DD0*DD;
  
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    //double ORDER;
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    dealii::Table<1,double> FF(n_q_points), Integral(kcells);
    
    for (unsigned int q=0; q<n_q_points; ++q) {
      Point<dim> qPoint=fe_values.quadrature_point(q);                
      if (ck==0) {
	//adding residual for the continuity equation
	R[i]+= (1.0)*(1.0/dt)*fe_values.shape_value(i, q)*(phi[q]-phi_conv[q])*fe_values.JxW(q);
	for (unsigned int j = 0; j < dim; j++){	
	  R[i]+=-(1.0)*fe_values.shape_value(i, q)*(ALPHA-phi[q])*(vel_j[q][j])*fe_values.JxW(q);
	}
	
      }
      
      else if(ck==1) {
	//additing residual for the velocity equation
	R[i]+= (1.0)*fe_values.shape_value(i, q)*(vel[q])*fe_values.JxW(q);  
	for (unsigned int j = 0; j < dim; j++) {
	  //constant
	  R[i]+=(AA)*(1.0)*fe_values.shape_value(i, q)*((press_j[q][j]+BB))*fe_values.JxW(q);
	  // linear
	  //R[i] +=(AA)*(pow(phi[q],1))*fe_values.shape_value(i, q)*((press_j[q][j]+BB))*fe_values.JxW(q);
	  //cubic
	  //R[i] +=(AA)*(pow(phi[q],3))*fe_values.shape_value(i, q)*((press_j[q][j]+BB))*fe_values.JxW(q);
	  //higher order
	  //R[i] +=(AA)*(pow(phi[q],5.86))*fe_values.shape_value(i, q)*((press_j[q][j]+BB))*fe_values.JxW(q);
	  
	}
      }
      
      else if(ck==2) {

	if (method==0) {
	  FF[q]=0;
	  for (unsigned int k = 1; k < kcells; k++) {
	    FF[q]+=(std::exp(-currentTime*tR[k])/E[k]/tau[k]/tau[k])*trapezoidal(dt,tR[k],CellHist[q][0],CellHist[q][1]);
	    //FF[q]+=(std::exp(-currentTime*tR[k])/E[k]/tau[k]/tau[k])*higherOrderIntegral(dt,tR[k],CellHist[q][0],CellHist[q][1]);
	  }	
	  FF[q]=FF0*FF[q];
	  
	}
	else if (method==1) {
	  FF[q]=0;
	  Integration<double, dim>(dt, q, currentTime, History, Integral ,  tR , CellHist[q][0], CellHist[q][1]);	
	  for (unsigned int k = 1; k < kcells; k++) {
	    FF[q]+=(std::exp(-currentTime*tR[k])/E[k]/tau[k]/tau[k])*Integral[k];	  
	  }
	  FF[q]=FF0*FF[q];
	}
	       	
	R[i]+=(CC)*(1.0/dt)*fe_values.shape_value(i, q)*(phi[q]*(press[q]-press_conv[q]))*fe_values.JxW(q);
	R[i]+=(DD)*fe_values.shape_value(i, q)*(phi[q]*press[q])*fe_values.JxW(q);
	R[i]+=-fe_values.shape_value(i, q)*(phi[q]*FF[q])*fe_values.JxW(q);
	
	for (unsigned int j = 0; j < dim; j++) {
	  R[i] +=(1.0)*fe_values.shape_value(i, q)*(vel_j[q][j])*fe_values.JxW(q);
	}
      }
            
    }
  } 
  
}

#endif /* CHEMO_H_ */
