//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Code for Bubble propogation which includes evolution of concentration,shape and mechanics 
//Created May 2020
//authors: kpbhagat(2020) rudraa (2018)
//

//deal.II headers
#include "include/headers.h"
//input parameter headers
#include "parameters.h"
//physics headers
#include "include/chemo.h"
//#include "include/historyclass.h"

//#include "include/mechanics.h"
//#include "include/writeSolutions.h"


//Namespace
namespace phaseField1
{
  using namespace dealii;

  //Initial conditions
  template <int dim>
  class InitalConditions: public Function<dim>{
  public:
    InitalConditions (): Function<dim>(DIMS){}
    void vector_value (const Point<dim>   &p, Vector<double>   &values) const{
        Assert (values.size() == DIMS, ExcDimensionMismatch (values.size(), DIMS));	

	//Here we specify the initial conditons. The interface boundary condition is same as the initial conditons
	values(0)=1.0;//porosity	      
	values(1)=0.0; // velocity
	values(2)=0.0; //effective pressure

    }
  };

  template <int dim>
  struct ExportVariables{
    //public:
    //ExportVariables<dim>():
    //Time(0.0), Location(dim), Porosity(dim), Velocity(dim), Eff_Pressure(dim) {}
    
    //using std:::map to store time history variables
    
    double Time;
    double Location, Porosity,Velocity,Eff_Pressure;
    //dealii::Table<1, double > Location, Porosity,Velocity,Eff_Pressure;
    //double  Location, Porosity ,Velocity,Eff_Pressure;
    //Vector<double>  Location(dim), Porosity (dim),Velocity(dim),Eff_Pressure(dim);
    //std::map<typename DoFHandler<dim>::active_cell_iterator, std::vector< ExportVariables<dim>* > > Export;
  };

  
  template <int dim>
  class phaseField{
  public:
    phaseField ();
    ~phaseField ();
    void run ();
    std::map<typename DoFHandler<dim>::active_cell_iterator, std::vector< historyVariables<dim>* > > History;                                                                                           
    //std::vector<double>  HistoryTable[XSubRf][3][2];     
    //dealii::Table<3,double >  HistoryTable[XSubRf][3][3];
    std::vector<std::vector<std::vector<double>>> HistoryTable;
    std::vector<std::vector<std::vector<double>>> IntegralTable;
    double Peffective;
   private:
    void applyBoundaryConditions(const unsigned int increment, const double time, const double deltaLoad);
    void setup_system ();
    void assemble_system ();
    void solveIteration (const unsigned int increment);
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int increment);
    void output_results_txt (const unsigned int increment,const double time);
    void storeHistory (const unsigned int increment,const double time, const double tR[], const double eC[], const double tauC[], const double Peffective);
    double Phi_average (const unsigned int increment, const double time, const double loadTop);
    Triangulation<dim>                        triangulation;
    FESystem<dim>                             fe;
    DoFHandler<dim>                           dof_handler;
    ConstraintMatrix                          constraints,constraints2;
    SparseMatrix<double>                      system_matrix;
    Vector<double>                            locally_relevant_solution, U, Un, UGhost, UnGhost, dU;
    Vector<double>                            system_rhs;
    ConditionalOStream                        pcout;
    TimerOutput                               computing_timer;
    SparsityPattern                           sparsity_pattern;

    //solution variables
    unsigned int currentIncrement, currentIteration;
    double totalTime, currentTime, dt;
    std::vector<std::string> nodal_solution_names; std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;    
    //std::map<typename DoFHandler<dim>::active_cell_iterator, std::vector< ExportVariables<dim>* > > Export;

    std::map<typename DoFHandler<dim>::active_cell_iterator,std::vector<double> > Export;
    //std::map<typename DoFHandler<dim>::active_cell_iterator,std::vector<double> > History;
    //std::map<typename DoFHandler<dim>::active_cell_iterator, std::vector<double> > History;
    //std::map<typename DoFHandler<dim>::active_cell_iterator, std::vector< historyVariables<dim>* > > History;

    //dealii::Table<3,double > HistoryTable(40,3,3); //cell, quadpoint, kcellno    
    std::vector<double>  stress[XSubRf][3][2]; //#cell,#quadpoint#time,press
    //const Table<3,double >  HistoryTable[XSubRf][3][3];

  };
  
  template <int dim>
  phaseField<dim>::phaseField ():
    fe(FE_Q<dim>(FEOrder),DIMS),
    dof_handler (triangulation),
    pcout (std::cout),
    computing_timer (pcout, TimerOutput::summary, TimerOutput::wall_times){
    //solution variables
    dt=TimeStep;
    totalTime=TotalTime;
    currentIncrement=0; currentTime=0;

    //nodal Solution names
    nodal_solution_names.push_back("phi"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    nodal_solution_names.push_back("vel"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    nodal_solution_names.push_back("press"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);


     
  }
  
  template <int dim>
  phaseField<dim>::~phaseField (){
    dof_handler.clear ();
    system_matrix.clear();
  }

  //Apply boundary conditions
  template <int dim>
  void phaseField<dim>::applyBoundaryConditions(const unsigned int increment, const double time, const double deltaLoad){
    constraints.clear (); constraints2.clear ();  
    //constraints.reinit (locally_relevant_dofs);
    //constraints2.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints2);

    //Setup boundary conditions
    //two dirichlet boundary condition on the top
    //two dirichlet boundary condition on the bottom
    std::vector<bool> top (DIMS, false); top[2]=true;              
    std::vector<bool> bottom (DIMS, false); bottom[1]=true; //bottom[2]=true;               


    std::vector<double> valueBottom (DIMS);    
    valueBottom[0]=0.0; //porosity nobc 
    valueBottom[1]=0.0 ; //velocity  fixed velocity
    valueBottom[2]=0.0; //effective pressure nobc  
   
    std::vector<double> valueTop (DIMS);    
    valueTop[0]=0.0; //porosity nobc
    valueTop[1]=0.0; //velocity nobc
    //valueTop[2]=1.0; //pressure p0
    valueTop[2]=deltaLoad; //pressure p0

    
    //bottom
    VectorTools::interpolate_boundary_values (dof_handler, 0, ConstantFunction<dim>(valueBottom), constraints, bottom);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(DIMS), constraints2, bottom);

    //top
    VectorTools::interpolate_boundary_values (dof_handler, 1, ConstantFunction<dim>(valueTop) , constraints, top);
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(DIMS), constraints2, top);
    
    constraints.close ();
    constraints2.close ();

  }

  
  
  //Setup
  template <int dim>
  void phaseField<dim>::setup_system (){
    TimerOutput::Scope t(computing_timer, "setup");
    dof_handler.distribute_dofs (fe);

    //Solution vectors
    system_rhs.reinit (dof_handler.n_dofs());
    U.reinit (dof_handler.n_dofs());
    Un.reinit (dof_handler.n_dofs());    
    dU.reinit (dof_handler.n_dofs());
    
    
    //call applyBoundaryConditions to setup constraints matrix needed for generating the sparsity pattern
    //applyBoundaryConditions(0);
    applyBoundaryConditions(0,0,0); 
    
    DynamicSparsityPattern dsp (dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints2, false);
    // SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit (sparsity_pattern);

    //setup history variables
    //HistoryTable.fill(0,true);
    
    //setup history variables
    const QGauss<dim>  quadrature_formula(FEOrder+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell){
      if (cell->is_locally_owned()){
	std::vector<std::vector<double>> v2d;
	
	//History[cell]=new std::vector< historyVariables<dim>* > (fe_values.n_quadrature_points );
	for (unsigned int q=0; q<fe_values.n_quadrature_points; q++){
	  History[cell].push_back(new historyVariables<dim>); //create histroy variables object at each quad point of the cell.
	  History[cell].back()->integralOld=0.0;
	  std::vector<double> v1d;
	  
	  for(unsigned int k=0;k<kcells;k++)
	    {
	      History[cell].back()->integralOldCells[k]=0.0;
	      v1d.push_back(0);
	    }
	   v2d.push_back(v1d);
	  
	}
	HistoryTable.push_back(v2d);
	IntegralTable.push_back(v2d);
      }
    }



                   
  }

  //Setup                                                                                                                                  

  //Assembly
  template <int dim>
  void phaseField<dim>::assemble_system (){
    TimerOutput::Scope t(computing_timer, "assembly");
    system_rhs=0.0; system_matrix=0.0;
    const QGauss<dim>  quadrature_formula(FEOrder+1);
    //   const QGauss<dim-1>	face_quadrature_formula (dim);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |
                             update_JxW_values);
    // FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, update_values |update_gradients |update_quadrature_points | update_JxW_values | update_normal_vectors);

    
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    unsigned int n_q_points= fe_values.n_quadrature_points;

    //std::vector<double>  CellHist[3][2]; //#cell,#quadpoint#time,press
    Table<2,std::vector<double> > CellHist(3,2); 
    Table<2,double > HistoryCellTable(3,4);
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned()) {
	fe_values.reinit (cell);
	local_matrix = 0; local_rhs = 0; 
	cell->get_dof_indices (local_dof_indices);
	 //AD variables
	Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell); Table<1, double > ULocalConv(dofs_per_cell);
	for (unsigned int i=0; i<dofs_per_cell; ++i){
	  if (std::abs(U(local_dof_indices[i]))<1.0e-16){ULocal[i]=0.0;}
	  else{ULocal[i]=U(local_dof_indices[i]);}
	  // std::cout<<" "<< ULocal[i].val(); //<<std::endl;
	  ULocal[i].diff (i, dofs_per_cell);
	  ULocalConv[i]= Un(local_dof_indices[i]);
	}
	
		
	//setup residual vector
	Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); 
	for (unsigned int i=0; i<dofs_per_cell; ++i) {
	  R[i]=0.0;
	}

	for (unsigned int q=0; q<n_q_points; ++q) {
	  CellHist[q][0]=stress[cell->active_cell_index()][q][0];
	  CellHist[q][1]=stress[cell->active_cell_index()][q][1];
	  for (unsigned int k=0; k<kcells; ++k) {	 
	    //HistoryCellTable[q][k]=HistoryTable[cell->active_cell_index()][q][k];
	    HistoryCellTable[q][k]=IntegralTable[cell->active_cell_index()][q][k];

	  }
	}

	//std::cout<<"printing dt " <<CellHist[0][0][1]<<"\n";
	//populate residual vector 
	//residualForChemo(fe_values, 0,cell, dt, ULocal, ULocalConv, R, History[cell], CellHist, currentTime, totalTime);
	//residualForChemo(fe_values, 0,cell, dt, ULocal, ULocalConv, R, History[cell], CellHist, currentIteration ,currentTime, totalTime);
	residualForChemo(fe_values, 0,cell, dt, ULocal, ULocalConv, R, History[cell], HistoryCellTable, CellHist, currentIteration ,currentTime, totalTime, Peffective);                                                              

	//evaluate Residual(R) and Jacobian(R')
	for (unsigned int i=0; i<dofs_per_cell; ++i) {
	  for (unsigned int j=0; j<dofs_per_cell; ++j){
	    // R' by AD
	    //local_matrix(i,j)= R[i].fastAccessDx(j);
	    local_matrix(i,j)= R[i].dx(j);
	  }
	  //R
	  local_rhs(i) = -R[i].val();
	}

	if ((currentIteration==0) /*&&(currentIncrement==1)*/){
          constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
        }
        else{
          constraints2.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
        }
	
      }
    system_matrix.compress (VectorOperation::add);
    system_rhs.compress (VectorOperation::add);
    //pcout << "\nK norm is: " << system_matrix.frobenius_norm() << std::endl; 

  }
  
                                                                                                                                                                                                 
  //Solve                                                                                                                                                                                                
  template <int dim>                                                                                                                                                                                     
  void phaseField<dim>::solveIteration(unsigned int currentIncrement){                                                                                                                                                                
    TimerOutput::Scope t(computing_timer, "solve");                                                                                                                                                      

    SparseDirectUMFPACK A_direct;                                                                                                                                                                       
    A_direct.initialize(system_matrix);                                                                                                                                                                 
    A_direct.vmult(dU, system_rhs);                                                                                                                                                                     
    if ((currentIteration==0) ){                                                                                                                                                                        
      constraints.distribute (dU);                                                                                                                                                                      
    }                                                                                                                                                                                                   
    else{                                                                                                                                                                                               
      constraints2.distribute (dU);                                                                                                                                                                      
    }         
  }
  
  /*  
  //Solve
  template <int dim>
 void phaseField<dim>::solveIteration(){
    // TimerOutput::Scope t(computing_timer, "solve");
    //LA::MPI::Vector completely_distributed_solution (locally_owned_dofs, mpi_communicator);
    //Direct solver MUMPS  
    SolverControl cn;
    SparseDirectMUMPS solver(cn);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, du, system_rhs);
    if ((currentIteration==0) ){                                                                                                                                                                         
      constraints.distribute (dU);                                                                                                                                                                       
    }                                                                                                                                                                                                   
    else{                                                                                                                                                                                               
      constraints2.distribute (dU);                                                                                                                                                                      
     }           
  }
*/
  
  
  //Solve
  template <int dim>
  void phaseField<dim>::solve(){
    double res=1, tol=1.0e-12, abs_tol=1.0e-14, initial_norm=0, current_norm=0;
    double machineEPS=1.0e-15;
    currentIteration=0;
    char buffer[200];
    while (true){
      if (currentIteration>=NR_ITR){sprintf(buffer, "maximum number of iterations reached without convergence. \n"); pcout<<buffer; break;}
      if (current_norm>1/std::pow(tol,2)){sprintf(buffer, "\n norm is too high. \n\n"); pcout<<buffer; break; exit (1);}
      assemble_system();
      current_norm=system_rhs.l2_norm();
      initial_norm=std::max(initial_norm, current_norm);
      res=current_norm/initial_norm;
      sprintf(buffer,"inc:%3u (time:%10.3e, dt:%10.3e), iter:%2u, abs-norm: %10.2e, rel-norm: %10.2e\n", currentIncrement, currentTime, dt,  currentIteration, current_norm, res); pcout<<buffer; 
      if ((currentIteration>1) && ((res<tol) || (current_norm<abs_tol))){sprintf(buffer,"residual converged in %u iterations.\n\n", currentIteration); pcout<<buffer; break;}
      solveIteration(currentIncrement);
      U+=dU;  
      ++currentIteration;
    }
    Un=U;    
  }

  template <int dim>
  double phaseField<dim>::Phi_average (const unsigned int ITR, const double time, const double loadTop)  {
    TimerOutput::Scope t(computing_timer, "bubbleVolume");
    const QGauss<dim>  quadrature_formula(FEOrder+1);
    //   const QGauss<dim-1> face_quadrature_formula (2);
    FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |  update_gradients |
			     update_quadrature_points|update_JxW_values);
    //unsigned int dofs_per_cell= fe_values.dofs_per_cell;
    unsigned int n_q_points= fe_values.n_quadrature_points;
    //const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;
  

    std::vector<Vector<double> > quadSolutions;
    std::vector< std::vector< Tensor< 1, dim, double >>> quadSolutionsGrad ;

    for (unsigned int q=0; q<n_q_points; ++q){
      quadSolutions.push_back(dealii::Vector<double>(DIMS)); //3 since there are three degrees of freedom per cell (phi,eta,c)
      quadSolutionsGrad.push_back(std::vector< Tensor< 1, dim, double >>(DIMS));

    }
  
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  
    double averagePhi=0,  averagePeff=0, averageGradVs=0;
   
    for (;cell!=endc; ++cell) {
      if (cell->is_locally_owned()) {
	fe_values.reinit(cell);
	fe_values.get_function_values(Un, quadSolutions);
	fe_values.get_function_gradients(Un,quadSolutionsGrad);	

	for (unsigned int q=0; q<n_q_points; ++q) {
	  averagePhi+= (quadSolutions[q][0])*fe_values.JxW(q);
	  averagePeff+= (quadSolutions[q][2])*fe_values.JxW(q);
	  averageGradVs+=(quadSolutionsGrad[q][1][0])*fe_values.JxW(q);
	}
      }
      //++t_cell;
    }

    char buffer[200];
    sprintf(buffer, "  Average of Phi:     %14.9e, Average of Peff:     %14.9e :  \n", averagePhi, averagePeff);
    pcout<<buffer;
    
    //std::ofstream myfile ("file.txt");
    std::ofstream myfile("file.txt", std::ofstream::out | std::ofstream::app);
    if (myfile.is_open()) {
      myfile<<ITR<<" "<<time<<" "<<averagePhi<<" "<<averagePeff<<" "<<averageGradVs<<" "<<loadTop<<"\n";	
      myfile.close();
    }
    else pcout << "Unable to open file";

    // totalBubbleVolume=volumeTotal;
    //surfaceAreaTotal= Utilities::MPI::sum(surfaceArea, mpi_communicator);
    //sprintf(buffer, "  Surface Area of bubble is:     %14.9e\n", surfaceAreaTotal);
    //pcout<<buffer;
    //pcout << "  Volume of bubble  is :       "
    //      <<volumeTotal <<std::endl;

    return averagePeff;
  }

  
  
  //Output
  template <int dim>
  void phaseField<dim>::output_results (const unsigned int cycle) {
    TimerOutput::Scope t(computing_timer, "output");
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (Un, nodal_solution_names, DataOut<dim>::type_dof_data, nodal_data_component_interpretation);    
    
    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector (subdomain, "subdomain");    
    data_out.build_patches ();
    
    const std::string filename = ("solution-" +
                                  Utilities::int_to_string (cycle, 2) );
    std::ofstream output ((filename + ".vtu").c_str());
    data_out.write_vtu (output);
    
  }

  //Output
  template <int dim>
  void phaseField<dim>::output_results_txt (const unsigned int cycle, const double time) {
    TimerOutput::Scope t(computing_timer, "output_text");
    //FEValues<dim> fe_values (fe);
    const QGauss<dim>  quadrature_formula(FEOrder+1);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |update_JxW_values);
    unsigned int dofs_per_cell= fe_values.dofs_per_cell;
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    const std::string filename = ("T-" +
                                Utilities::int_to_string (cycle, 2) );

    //const std::string filename = ("T-" +
    //                            Utilities::to_string (time, 2) );
   
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    //typename parallel::distributed::Triangulation<dim>::active_cell_iterator t_cell = triangulation.begin_active();
   
    for (;cell!=endc; ++cell) {
      if (cell->is_locally_owned()) {
	fe_values.reinit(cell);
	cell->get_dof_indices (local_dof_indices);
	//std::cout <<"dof per cell is " << dofs_per_cell << " vertices per cell is " << GeometryInfo<dim>::vertices_per_cell << std::endl; 

	//Export[cell][i]->Time=time;
	Export[cell].push_back(time);
	//for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex) {
	for (unsigned int ver=0; ver<1; ++ver) {
	  //Export[cell].push_back(GeometryInfo<dim>::unit_cell_vertex(vertex)[0]);
	  Export[cell].push_back(cell->vertex(ver)[0]);
	}
	
	//for (unsigned int i=0; i<dofs_per_cell; ++i) {
	for (unsigned int i=0; i<3; ++i) {
	  const unsigned int ci = fe_values.get_fe().system_to_component_index(i).first - 0;   
	  if (ci==0) {
	    //add phi
	    //Export[cell][i]->Porosity[0]=Un(local_dof_indices[i]);
	    Export[cell].push_back(Un(local_dof_indices[i]));
	  }
	  
	  if (ci==1) {
	    //add vel
	    //Export[cell][i]->Velocity[0]=Un(local_dof_indices[i]);
	    //Export[cell].Velocity=Un(local_dof_indices[i]);
	    Export[cell].push_back(Un(local_dof_indices[i]));
	  }

	  if (ci==2) {
	    //add eff press
	    //Export[cell][i]->Eff_Pressure[0]=Un(local_dof_indices[i]);
	    //Export[cell].Eff_Pressure=Un(local_dof_indices[i]);
	    Export[cell].push_back(Un(local_dof_indices[i]));	    
	  }
	  //Export[cell]=time;
	}
	//++t_cell;
      }
    }

    typename std::map<typename DoFHandler<dim>::active_cell_iterator, std::vector<double>  >::iterator it = Export.begin();	
    std::ofstream myfile((filename + ".text").c_str(), std::ofstream::out | std::ofstream::app);
    if (myfile.is_open()) {
      //char buffer[200];
      //sprintf(buffer, "  Average of Phi:     %14.9e, Average of Peff:     %14.9e :  \n", averagePhi, averagePeff);
      //pcout<<buffer;
      myfile<<"Tim "<<"Pos "<<"Phi "<<"Vel "<<"Peff "<<std::endl ;
      for(; it != Export.end(); it++) {
	//myfile<<"Ce is : "<<it->first ;
	int counter=0; 
	for(auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
	  myfile<<" "<<std::setprecision (18) <<*it2;
	  //if (counter==0) myfile<<" ti is : "<<*it2 ;
	  //if (counter==1) myfile<<" po is : "<<*it2 ;
	  //if (counter==2) myfile<<" ph is : "<<*it2 ;
	  //if (counter==3) myfile<<" ve is : "<<*it2 ;
	  //if (counter==4) myfile<<" pe is : "<<*it2 ;
	  //myfile << "This is another line.\n";
	  counter=counter+1;
	}
	myfile<<std::endl ;
      }
      myfile.close();
    }
    else pcout << "Unable to open file";
    Export.clear();

    
	/*
	FILE *fp = fopen((filename + ".text").c_str(), "w");	
	if (fp){
	  for(; it != Export.end(); it++) {
	    //fprintf(fp, "%s=%s\n", it->first, it->second);
	    //fprintf(fp, " %14.9e, %14.9e  \n", it->first, it->second);
	    for(auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
	      fprintf(fp, "Val%14.9e ", *it2);
	      if (it2 == it->second.end() ){fprintf(fp, "\n", " ");}
	    }
	  }

	}
	fclose(fp);
	*/
	


      //const std::string filename = ("solution-" +
    //Utilities::to_string (time) );

  }

  


    //Output
  template <int dim>
  void phaseField<dim>::storeHistory(const unsigned int cycle, const double time, const double tR[], const double eC[], const double tauC[], const double Peffective ) {
    TimerOutput::Scope t(computing_timer, "output_text");
    const QGauss<dim>  quadrature_formula(FEOrder+1);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |update_JxW_values);
    unsigned int dofs_per_cell= fe_values.dofs_per_cell;
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    double timeSize =0, f0 = 0, f1 = 0;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    unsigned int n_q_points= fe_values.n_quadrature_points;
    std::vector<Vector<double> > quadSolutions;
    //    std::vector<double>  stress[triangulation.n_global_active_cells()][n_q_points][2];


    for (unsigned int q=0; q<n_q_points; ++q){
      quadSolutions.push_back(dealii::Vector<double>(DIMS)); //4 since there are four degrees of freedom per cell (phi,eta,c,xi)
    }

    double checkPressure =0 ;
    for (;cell!=endc; ++cell) {
      const std::string filename = ("Cell-" +
				    Utilities::int_to_string (cell->active_cell_index(), 2) );      
      if (cell->is_locally_owned()) {
	fe_values.reinit(cell);
	cell->get_dof_indices (local_dof_indices);
	fe_values.get_function_values(Un, quadSolutions);
	
	for (unsigned int q=0; q<n_q_points; ++q) {
	  Point<dim> quadPoint=fe_values.quadrature_point(q);                    
	  stress[cell->active_cell_index()][q][0].push_back(time); //Store time
	  stress[cell->active_cell_index()][q][1].push_back(quadSolutions[q][DIMS-1]); //Store effective pressure

	  //Trapezoidal rule
	  for (unsigned int k=1; k<kcells; ++k) {
	    //timesize
	    int last = stress[cell->active_cell_index()][q][0].size();
	    if (last>=2) {
	       timeSize = stress[cell->active_cell_index()][q][0][last-1]-stress[cell->active_cell_index()][q][0][last-2];
	       f1 = (1.0/eC[k]/tauC[k]/tauC[k])*std::exp(stress[cell->active_cell_index()][q][0][last-1]*tR[k])*stress[cell->active_cell_index()][q][1][last-1];
	       f0 = (1.0/eC[k]/tauC[k]/tauC[k])*std::exp(stress[cell->active_cell_index()][q][0][last-2]*tR[k])*stress[cell->active_cell_index()][q][1][last-2];
	       checkPressure = stress[cell->active_cell_index()][q][1][last-1];
	    }
	    else if (last==1) {
	       timeSize = stress[cell->active_cell_index()][q][0][last-1];
               f1 = (1.0/eC[k]/tauC[k]/tauC[k])*std::exp(stress[cell->active_cell_index()][q][0][last-1]*tR[k])*stress[cell->active_cell_index()][q][1][last-1];
               f0 = 0 ;
	    }
	    else if (last==0) {
	      timeSize = 0;
	      f1 = 0;
	      f0 = 0;
	    }
	    
	    HistoryTable[cell->active_cell_index()][q][k]+=(0.5*timeSize*(f0+f1));
	    bool checkNan = std::isnan(HistoryTable[cell->active_cell_index()][q][k] ) ;
	    bool checkInf = std::isinf(HistoryTable[cell->active_cell_index()][q][k] ) ;
	    if (Peffective >= 1.0 ) {
	      IntegralTable[cell->active_cell_index()][q][k]=(1-std::exp(-time*tR[k]))*stress[cell->active_cell_index()][q][1][last-1]/eC[k]/tauC[k] ;			    
	    }

	    else  {
	      IntegralTable[cell->active_cell_index()][q][k]=std::exp(-stress[cell->active_cell_index()][q][0][last-1]*tR[k])*HistoryTable[cell->active_cell_index()][q][k];	      
	    }
	    
	    //std::cout << "Integral Table and Pressure is " << IntegralTable[cell->active_cell_index()][q][k]<< " and " << stress[cell->active_cell_index()][q][1][last-1] <<std::endl;
	    
	  }
	  
	  if (q==0) {
	    std::ofstream myfile((filename + ".text").c_str(), std::ofstream::out | std::ofstream::app);
	    if (myfile.is_open())
	      {
		myfile<<quadPoint[0]<<" "<<stress[cell->active_cell_index()][0][0].back()<<" "<<stress[cell->active_cell_index()][0][1].back()<<"\n";
		myfile.close();
	      }
	    else pcout << "Unable to open file";	    	    
	  }	  
	}				
      }
    }

    
  }
  
  
  //Solve problem
  template <int dim>
  void phaseField<dim>::run (){

    //setup problem geometry and mesh

    /*GridGenerator::hyper_cube (triangulation, -problemLength/2.0, problemLength/2.0,true);
    //  GridGenerator::hyper_cube (triangulation, -10/2.0, 10/2.0, true);
    triangulation.refine_global (globalRefinementFactor);
    */
        
    std::vector<unsigned int> numRepetitions;
    numRepetitions.push_back(XSubRf); // x refinement
    Point<dim> p1 (0);
    Point<dim> p2 (problemLength);
    GridGenerator::subdivided_hyper_rectangle (triangulation, numRepetitions, p1, p2, true);     

    setup_system (); //inital set up
    
    pcout << "   Number of active cells:       "
	  << triangulation.n_global_active_cells()
	  << std::endl
	  << "   Number of degrees of freedom: "
	  << dof_handler.n_dofs()
	  << std::endl;
    
    //setup initial conditions
    VectorTools::interpolate(dof_handler, InitalConditions<dim>(), U); Un=U;
    
    //sync ghost vectors to non-ghost vectors
    //  UGhost=U;  UnGhost=Un;
    //output_results (0);
    
    //Time stepping
    currentIncrement=0;
    currentTime=0;

    pcout << "   coord. of mesh:       "
    	  << triangulation.get_vertices()[0]
    	  << std::endl;

    double deltaLoad = 0, Load = 0, deltaOld = 0;
    Peffective=0;
    // while (currentTime<totalTime) {
    for (currentTime=0; currentTime<totalTime; currentTime+=dt){
      //change dt based on some logic
      //if (currentTime <=ceil_0) {dt=TimeStep_0; }
      //else if (currentTime > ceil_0 && currentTime <=ceil_1) {dt=TimeStep_1;}
      //else if (currentTime > ceil_1 && currentTime <=ceil_2) {dt=TimeStep_2;}
      //else if (currentTime > ceil_2 && currentTime <=ceil_3) {dt=TimeStep_3;}
      
      currentIncrement++;
      if (dt < maxdt) {
	dt=(1.0+timeFactor)*dt;
      }
      //currentTime+=dt;

      
      if (currentTime <= RampUp ){	
	deltaLoad = dt/RampUp;
	deltaOld=deltaLoad;
	Load +=deltaLoad;
	if (Load > 1.0) {
	  deltaLoad= Load - 1.0;
	  Load = Load-deltaOld+deltaLoad;
	}
      }
      else {
	deltaLoad=0;
	pcout << "Load is : " <<Load<<std::endl;	
      }

      applyBoundaryConditions(currentIncrement,currentTime,deltaLoad);

      double tR[]=tRatio;
      double eC[]=EE;
      double tauC[]=ttau;
      //storeHistory(currentIncrement,currentTime,tR,eC,tauC);
      solve();
      //int NSTEP=(currentTime/dt);
      int NSTEP=(currentIncrement);

      if (NSTEP%PSTEPS2==0) {
	output_results(currentIncrement);
	output_results_txt(currentIncrement,currentTime);
	//writeSolutionsToFile(Un, tag)       
      }
      
      if (NSTEP%PSTEPS==0) {
	Peffective=Phi_average (currentIncrement,currentTime,Load);
      }
      
      storeHistory(currentIncrement,currentTime,tR,eC,tauC,Peffective);
     pcout << std::endl;
     
    }


    
    //computing_timer.print_summary ();
  }
}


int main(int argc, char *argv[]){
  try
    {
      using namespace dealii;
      using namespace phaseField1;
      //  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      phaseField<1> problem;
      problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
