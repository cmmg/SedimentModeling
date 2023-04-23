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


double Integration(double dtSize, double tTauRatio, std::vector<double>  &Time, std::vector<double>  &Peff) {
  //no. of steps
  double Integral=0;  
  for (unsigned int row=0; row<Time.size(); ++row) {    
    //func[0] is time
    //func[1] is pressure
    if(row==0) {
      Integral+=0.5*dtSize*std::exp(Time[row]*tTauRatio)*Peff[row] ;   
    }
    else if (row==Time.size()-1) {
      Integral+=0.5*dtSize*std::exp(Time[row]*tTauRatio)*Peff[row] ;   
    }

    else {
      Integral+=dtSize*std::exp(Time[row]*tTauRatio)*Peff[row];    
    }

  }
  
  return Integral;
}


//Chemistry residual implementation
template <int dim>
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF,  const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R, Table<2,std::vector<double> > &CellHist, double currentTime, double totalTime) {
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
    dealii::Table<1,double> FF(n_q_points);
    
    for (unsigned int q=0; q<n_q_points; ++q) {
      Point<dim> qPoint=fe_values.quadrature_point(q);  
      
      if (ck==0) {
	//adding residual for the continuity equation
	R[i] += (1.0)*(1.0/dt)*fe_values.shape_value(i, q)*(phi[q]-phi_conv[q])*fe_values.JxW(q);
	for (unsigned int j = 0; j < dim; j++){	
	  R[i] +=-(1.0)*fe_values.shape_value(i, q)*(ALPHA-phi[q])*(vel_j[q][j])*fe_values.JxW(q);
	}
	
      }
      
      else if(ck==1) {
	//additing residual for the velocity equation
	R[i] += (1.0)*fe_values.shape_value(i, q)*(vel[q])*fe_values.JxW(q);  
	for (unsigned int j = 0; j < dim; j++) {
	  //constant
	  R[i] +=(AA)*(1.0)*fe_values.shape_value(i, q)*((press_j[q][j]+BB))*fe_values.JxW(q);
	  // linear
	  //R[i] +=(AA)*(pow(phi[q],1))*fe_values.shape_value(i, q)*((press_j[q][j]+BB))*fe_values.JxW(q);
	  //cubic
	  //R[i] +=(AA)*(pow(phi[q],3))*fe_values.shape_value(i, q)*((press_j[q][j]+BB))*fe_values.JxW(q);
	  //higher order
	  //R[i] +=(AA)*(pow(phi[q],5.86))*fe_values.shape_value(i, q)*((press_j[q][j]+BB))*fe_values.JxW(q);
	  
	}
      }
      
      else if(ck==2) {
	FF[q]=0;
	for (unsigned int k = 1; k < kcells; k++) {
	  FF[q]+= (std::exp(-currentTime*tR[k])/E[k]/tau[k]/tau[k])*Integration(dt,tR[k],CellHist[q][0],CellHist[q][1]);
	}
	FF[q]= FF0*FF[q];

	R[i] +=(CC)*(1.0/dt)*fe_values.shape_value(i, q)*(phi[q]*(press[q]-press_conv[q]))*fe_values.JxW(q);
	R[i] +=(DD)*fe_values.shape_value(i, q)*(phi[q]*press[q])*fe_values.JxW(q);
	R[i] +=-fe_values.shape_value(i, q)*(phi[q]*FF[q])*fe_values.JxW(q);		
		
	for (unsigned int j = 0; j < dim; j++) {
	  R[i] +=(1.0)*fe_values.shape_value(i, q)*(vel_j[q][j])*fe_values.JxW(q);
	}
      }
            
    }
  } 
  
}

#endif /* CHEMO_H_ */
