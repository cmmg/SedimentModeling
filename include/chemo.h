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

//Chemistry residual implementation
template <int dim>
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF,  const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R, double currentTime, double totalTime) {
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;

 
  //evaluate field and gradients
  dealii::Table<1,Sacado::Fad::DFad<double> > phi(n_q_points), vel(n_q_points) , press(n_q_points), ORDER(n_q_points); 
  dealii::Table<2,Sacado::Fad::DFad<double> > phi_j(n_q_points, dim), vel_j(n_q_points,dim), press_j(n_q_points, dim);
   
  dealii::Table<1,double> phi_conv(n_q_points);
  dealii::Table<1,double> press_conv(n_q_points);
  // dealii::Table<1,double> ORDER_conv(n_q_points);
  
  for (unsigned int q=0; q<n_q_points; ++q) {
    phi[q]=0.0; vel[q]=0.0; press[q]=0.0; phi_conv[q]=0.0; press_conv[q]=0.0; ORDER[q]=0;

    for (unsigned int j=0; j<dim; j++) {phi_j[q][j]=0.0; vel_j[q][j]=0.0; press_j[q][j]=0.0;}
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;      
      if (ck==0) { phi[q]+=fe_values.shape_value(i, q)*ULocal[i]; phi_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i]; }
      else if (ck==1){ vel[q]+=fe_values.shape_value(i, q)*ULocal[i]; }
      else if (ck==2) {press[q]+=fe_values.shape_value(i, q)*ULocal[i]; press_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i]; }
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

  

   
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    //double ORDER;
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    Sacado::Fad::DFad<double> BB=bb;
    for (unsigned int q=0; q<n_q_points; ++q) {
      //Point<dim> qPoint=fe_values.quadrature_point(q);  
         
      if (ck==0) {
	//adding residual for the conituity equation
	R[i] += (1.0)*(1.0/dt)*fe_values.shape_value(i, q)*(phi[q]-phi_conv[q])*fe_values.JxW(q);
	for (unsigned int j = 0; j < dim; j++){	
	  R[i] +=-fe_values.shape_value(i, q)*(ALPHA-phi[q])*(vel_j[q][j])*fe_values.JxW(q);
	}
	
      }
      
      else if(ck==1) {
	//additing residual for the velocity equation
	R[i] += (1.0)*fe_values.shape_value(i, q)*(vel[q])*fe_values.JxW(q);
	
	for (unsigned int j = 0; j < dim; j++) {
	  R[i] +=(aa)*(pow(phi[q],3))*fe_values.shape_value(i, q)*(press_j[q][j]+BB)*fe_values.JxW(q);	  
	}
      }

      else if(ck==2) {
	//adding residual for the pressure equation
	R[i] +=(dd)*(1.0/dt)*fe_values.shape_value(i, q)*(phi[q]*(press[q]-press_conv[q]))*fe_values.JxW(q);
	R[i] +=(cc)*fe_values.shape_value(i, q)*(phi[q]*press[q])*fe_values.JxW(q);


	for (unsigned int j = 0; j < dim; j++) {
	  R[i] +=fe_values.shape_value(i, q)*(vel_j[q][j])*fe_values.JxW(q);
	}
      }

      else if (ck==3) {
	//adding residula for the order parameter
	//R[i]+=fe_values.shape_value(i, q)*(ORDER[q]-0.5*(1.0+std::tanh(200.0*(Vel*currentTime-qPoint[0]))))*fe_values.JxW(q);
	
      }
    }
  } 
  
}

#endif /* CHEMO_H_ */
