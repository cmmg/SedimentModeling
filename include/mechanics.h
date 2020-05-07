//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2011
//authors: rudraa (2011, 2018)
//
#ifndef MECHANICS_H_
#define MECHANICS_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"
#include "deformationMap.h"


//Mechanics implementation
template <class T, int dim>
  void evaluateStress(const FEValues<dim>& fe_values, FEFaceValues<dim> & fe_face_values, const unsigned int DOF, const Table<1, T>& ULocal, Table<3, T>& P,Table<3, T>& PFace , const deformationMap<T, dim>& defMap, typename DoFHandler<dim>::active_cell_iterator &cell){

  //number of quadrature poits
  unsigned int n_q_points= fe_values.n_quadrature_points;
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;

  unsigned int n_q_points_face= fe_face_values.n_quadrature_points;
  const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;
  
  Table<3, Sacado::Fad::DFad<double>> gradU(n_q_points, dim, dim);
  Table<3, Sacado::Fad::DFad<double>> gradUFace(n_q_points_face, dim, dim);
  
   //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	gradU[q][i][j]=0.0;
      }
    }
    for (unsigned int k=0; k<dofs_per_cell; ++k){
      unsigned int ck = fe_values.get_fe().system_to_component_index(k).first - DOF;
      if (ck>=0 && ck<2){
	for (unsigned int i=0; i<dim; ++i){
	  gradU[q][ck][i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU
	}
      }
    }
  }

  
  //loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    //Fe
    Table<2, Sacado::Fad::DFad<double> > Fe (dim, dim);
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	Fe[i][j]=defMap.F[q][i][j];
      }
    }
    //E
    Table<2, Sacado::Fad::DFad<double> > E (dim, dim);
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){	
	//small strain: E is epsilon
	E[i][j] = 0.5*(gradU[q][i][j]+gradU[q][j][i]);
      }
    }

    //S
    //Material moduli
    double Y=elasticModulus, nu=PoissonsRatio;
    //Lame parameters
    double lambda=(nu*Y)/((1+nu)*(1-2*nu)), mu=Y/(2*(1+nu));
    //double lambda=LAM, mu=MU;
    Table<2, Sacado::Fad::DFad<double> > S (dim, dim);
  
    double C11=lambda+2*mu, C12=lambda, C44=mu;
    if(dim==2){
      //Plane strain model
      S[0][0]=C11*E[0][0]+C12*E[1][1];
      S[1][1]=C12*E[0][0]+C11*E[1][1];
      S[0][1]=S[1][0]=C44*E[0][1];
    }
    else throw "dim not equal to 2";


  
    //P
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	//small strain
	P[q][i][j]=S[i][j];

	
      }
    }
    
  }

  
}

//mechanics residual implementation
template <int dim>
void residualForMechanics(FEValues<dim>& fe_values, FEFaceValues<dim>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double> >& ULocal, Table<1, double>& ULocalConv, Table<1, Sacado::Fad::DFad<double> >& R, deformationMap<Sacado::Fad::DFad<double>, dim>& defMap, typename DoFHandler<dim>::active_cell_iterator& cell, unsigned int CURRENT){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;

  unsigned int n_q_points_face= fe_face_values.n_quadrature_points;
  const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;
  

  dealii::Table<1,Sacado::Fad::DFad<double> > phi(n_q_points), eta(n_q_points) , c(n_q_points);	
  for (unsigned int q=0; q<n_q_points; ++q) {
      eta[q]=0.0;      
     for (unsigned int i=0; i<dofs_per_cell; ++i) {
       const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
     if (ck==3){ eta[q]+=fe_values.shape_value(i, q)*ULocal[i]; }            
     }
  }


  

  double PP=0;
  if (CURRENT > 10) {PP=Pressure;}
  //temporary arrays
  Table<3,Sacado::Fad::DFad<double> > P (n_q_points, dim, dim);
  Table<3,Sacado::Fad::DFad<double> > PFace (n_q_points, dim, dim);
  //evaluate stress
  evaluateStress<Sacado::Fad::DFad<double>, dim>(fe_values, fe_face_values, DOF, ULocal, P,PFace, defMap, cell);
  
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    if (ck>=0 && ck<2){
      // R = Grad(w)*P    
      for (unsigned int q=0; q<n_q_points; ++q){
	R[i] +=(eta[q])*(PP)*fe_values.shape_grad(i, q)[ck]*fe_values.JxW(q);
	
	for (unsigned int d = 0; d < dim; d++){
	  R[i] +=(eta[q])*fe_values.shape_grad(i, q)[d]*P[q][ck][d]*fe_values.JxW(q);
	  //R[i] +=(eta[q])*(Pressure)*fe_values.shape_grad(i, q)[d]*fe_values.JxW(q);
	}
      }
    }
  }



}

#endif /* MECHANICS_H_ */
