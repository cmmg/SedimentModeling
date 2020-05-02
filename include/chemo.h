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
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R, double currentTime, double totalTime) {
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;

  dealii::Table<1,double> eta_conv(n_q_points);
  dealii::Table<1,double> c_conv(n_q_points);

 
  //evaluate gradients 
  dealii::Table<1,Sacado::Fad::DFad<double> > phi(n_q_points), eta(n_q_points) , c(n_q_points); 
  dealii::Table<2,Sacado::Fad::DFad<double> > phi_j(n_q_points,dim), eta_j(n_q_points, dim), c_j(n_q_points, dim) ;
   
   for (unsigned int q=0; q<n_q_points; ++q) {
     phi[q]=0.0; eta[q]=0.0; c[q]=0.0; eta_conv[q]=0.0; c_conv[q]=0.0;

    for (unsigned int j=0; j<dim; j++) {phi_j[q][j]=0.0; eta_j[q][j]=0.0;  c_j[q][j]=0.0; }
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;      
      if (ck==0) { phi[q]+=fe_values.shape_value(i, q)*ULocal[i];  }
      else if (ck==1){ eta[q]+=fe_values.shape_value(i, q)*ULocal[i]; eta_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i];  }
      else if (ck==2) {c[q]+=fe_values.shape_value(i, q)*ULocal[i]; c_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i]; }

      for (unsigned int j=0; j<dim; j++) {
	if (ck==0) {
	  phi_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];
	}
	else if (ck==1)
	  { eta_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];	   
	}

	else if (ck==2)
	  { c_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];
	  }
      }
    }
        
  }

  
     
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    
    for (unsigned int q=0; q<n_q_points; ++q) {
      dealii::Table<1,Sacado::Fad::DFad<double> > bigM(dim),dgammadtheta(dim), theta(dim) ;
      Sacado::Fad::DFad<double> gamma ;
      Sacado::Fad::DFad<double> Unity=1.0 ;
      Sacado::Fad::DFad<double> f_a = 0.5*c[q]*c[q]/16.0 ;
      Sacado::Fad::DFad<double> f_a_c= 0.5*c[q]/8.0 ;
      Sacado::Fad::DFad<double> f_a_c_c = 0.5/8.0;
      Sacado::Fad::DFad<double> f_b = 0.5*(c[q]-1)*(c[q]-1)/16.0 ;
      Sacado::Fad::DFad<double> f_b_c= 0.5*(c[q]-1)/8.0;
      Sacado::Fad::DFad<double> f_b_c_c= 0.5/8.0 ;
      Sacado::Fad::DFad<double> HH = 3.0*eta[q]*eta[q] - 2.0*eta[q]*eta[q]*eta[q];
      Sacado::Fad::DFad<double> H_eta = 6.0*eta[q] - 6.0*eta[q]*eta[q];

      Sacado::Fad::DFad<double> ww[] =www;
      Sacado::Fad::DFad<double> aalpha[] =aaalpha;
      Sacado::Fad::DFad<double> oorient[][dim] = ooorient;
  
      
      Sacado::Fad::DFad<double> GradMagEta= eta_j[q][0]*eta_j[q][0]+eta_j[q][1]*eta_j[q][1] ;
      GradMagEta = std::sqrt(GradMagEta);

      
      if (GradMagEta < 1.0e-12) {
	theta[0]=0; theta[1]=0;
      }
      else {
	theta[0]=eta_j[q][0]/GradMagEta;
	theta[1]=eta_j[q][1]/GradMagEta;
      }

      dealii::Table<2,Sacado::Fad::DFad<double> > orient(n_orient,dim) ;

      Sacado::Fad::DFad<double> product=0 ;
      double HeavySide=0;
      gamma=1.0;dgammadtheta[0]=0;dgammadtheta[1]=0; 
      for (unsigned int i=0; i<n_orient; ++i) {
	Sacado::Fad::DFad<double> w= ww[i];
	Sacado::Fad::DFad<double> alpha= aalpha[i];
	orient[i][0]= oorient[i][0];
	orient[i][1]= oorient[i][1];
	
	//orient[1]= oorient[i];
	//orient[2]= oorient[i];
	
	product=orient[i][0]*theta[0] + orient[i][1]*theta[1];
	if (product > 0) HeavySide =1.0;
	gamma-=alpha*std::pow(product.val(),w)*HeavySide;

	dgammadtheta[0]-=w*alpha*std::pow(product.val(),w-1)*orient[i][0]*HeavySide;
	dgammadtheta[1]-=w*alpha*std::pow(product.val(),w-1)*orient[i][1]*HeavySide;
	HeavySide=0; //reset Heavy side
      }
      gamma*=gamma0;
      dgammadtheta[0]*=gamma0;
      dgammadtheta[1]*=gamma0;

      bigM[0]=0; bigM[1]=0;
      for (unsigned int i = 0; i < dim; i++) {
	for (unsigned int j = 0; j < dim; j++) {
	  bigM[i]+=-theta[i]*theta[j]*dgammadtheta[j]*GradMagEta;
	  if (i==j) bigM[i]+=dgammadtheta[j]*GradMagEta;
	}
	bigM[i]+=gamma*eta_j[q][i];
	
      }
      bigM[0]=gamma*bigM[0];
      bigM[1]=gamma*bigM[1];
      
      //std::cout << "theta  x  is " <<theta[0].val() <<std::endl;
      //std::cout << "theta y  is " <<theta[1].val() <<std::endl;
      
      if (ck==0) {
	R[i] += fe_values.shape_value(i, q)*(phi[q])*fe_values.JxW(q);
	for (unsigned int j = 0; j < dim; j++){	
	  R[i] += fe_values.shape_grad(i, q)[j]*eta_j[q][j]*fe_values.JxW(q);  
	}
	
      }
      
      else if(ck==1) {	
	R[i] +=  (1.0/dt)*fe_values.shape_value(i, q)*(eta[q] - eta_conv[q])*fe_values.JxW(q);  
	R[i] +=  (M_eta)*fe_values.shape_value(i, q)*(f_b-f_a)*(H_eta)*fe_values.JxW(q);

	for (unsigned int j = 0; j < dim; j++) {
	  R[i]+= (M_eta)*fe_values.shape_grad(i, q)[j]*bigM[j]*fe_values.JxW(q);
	  R[i]+=-(M_eta)*(delt*delt)*fe_values.shape_grad(i, q)[j]*phi_j[q][j]*fe_values.JxW(q);
	}
      }

      else if(ck==2) {
	R[i] +=  (1.0/dt)*fe_values.shape_value(i, q)*(c[q] - c_conv[q])*fe_values.JxW(q);
	
	for (unsigned int j = 0; j < dim; j++) {
	  R[i]+= (M_c)*fe_values.shape_grad(i, q)[j]*(f_a_c_c*(1.0-HH)+f_b_c_c*HH)*(c_j[q][j])*fe_values.JxW(q);
	  R[i]+= (M_c)*fe_values.shape_grad(i, q)[j]*(f_b_c-f_a_c)*(H_eta)*(eta_j[q][j])*fe_values.JxW(q); 
	  
	}

      }

 

    }
  } 
  
}

#endif /* CHEMO_H_ */
