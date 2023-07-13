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
	values(2)=0.0; //pressure
	values(3)=0.0; //order
    }
  };

  template <int dim>
  struct ExportVariables{
    double Time;
    double Location, Porosity,Velocity,Eff_Pressure;
  };
  
  template <int dim>
  class phaseField{
  public:
    phaseField ();
    ~phaseField ();
    void run ();

   private:
    void applyBoundaryConditions(const unsigned int increment);
    void setup_system ();
    void assemble_system ();
    void solveIteration ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int increment);
    void output_results_txt (const unsigned int increment,const double time);
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

    std::map<typename DoFHandler<dim>::active_cell_iterator,std::vector<double> > Export;


  };

  template <int dim>
  phaseField<dim>::phaseField ():
    fe(FE_Q<dim>(FEOrder),DIMS),
    dof_handler (triangulation),
    pcout (std::cout),
    computing_timer (pcout, TimerOutput::summary, TimerOutput::wall_times){
    //solution variables
    dt=TimeStep; totalTime=TotalTime;
    currentIncrement=0; currentTime=0;

    //nodal Solution names
    nodal_solution_names.push_back("phi"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);

     nodal_solution_names.push_back("vel"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
     nodal_solution_names.push_back("press"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);

     nodal_solution_names.push_back("ORDER"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
     
  }
  
  template <int dim>
  phaseField<dim>::~phaseField (){
    dof_handler.clear ();
    system_matrix.clear();
  }

  //Apply boundary conditions
  template <int dim>
  void phaseField<dim>::applyBoundaryConditions(const unsigned int increment){
    constraints.clear (); constraints2.clear ();  
    //constraints.reinit (locally_relevant_dofs);
    //constraints2.reinit (locally_relevant_dofs);
    //DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints2);

    //Setup boundary conditions
    //two dirichlet boundary condition on the top
    //two dirichlet boundary condition on the bottom
    std::vector<bool> top (DIMS, false); top[0]=true; top[2]=true;              
    std::vector<bool> bottom (DIMS, false); bottom[1]=true; bottom[3]=true;              


    std::vector<double> valueBottom (DIMS);    
    valueBottom[0]=0; //porosity 
    valueBottom[1]=0.0 ; //1.53; //velocity
    valueBottom[2]=0.0; //pressure 
    valueBottom[3]=1.0; //order
    std::vector<double> valueTop (DIMS);    
    valueTop[0]=0.0; //porosity 
    valueTop[1]=0; //1.53; //velocity
    valueTop[2]=0.0; //pressure 
    valueTop[3]=0.0; //order 

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
    applyBoundaryConditions(0);
    
    DynamicSparsityPattern dsp (dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints2, false);
    // SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit (sparsity_pattern);
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
	
	//populate residual vector 
	residualForChemo(fe_values, 0,cell, dt, ULocal, ULocalConv, R, currentTime, totalTime);

	
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

	if ((currentIteration==0)&&(currentIncrement==1)){
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
  void phaseField<dim>::solveIteration(){
    TimerOutput::Scope t(computing_timer, "solve");
    //LA::MPI::Vector completely_distributed_solution (locally_owned_dofs, mpi_communicator);
    //Direct solver MUMPS

    /*
    SolverControl           cn(1000,1.0e-12);
    SolverCG<>              cg (cn);
    cg.solve (system_matrix, dU, system_rhs,
              PreconditionIdentity()); 
    if ((currentIteration==0)&&(currentIncrement==1)){
      constraints.distribute (dU);
    }
    else{
      constraints2.distribute (dU);
    }

    std::cout << "   " << cn.last_step()
              << " CG iterations needed to obtain convergence."
	              << std::endl;
    */

    
     SparseDirectUMFPACK A_direct;
     A_direct.initialize(system_matrix);
     A_direct.vmult(dU, system_rhs);
     if ((currentIteration==0)&&(currentIncrement==1)){
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
    PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, completely_distributed_solution, system_rhs);
    constraints.distribute (completely_distributed_solution);
    locally_relevant_solution = completely_distributed_solution;
    dU = completely_distributed_solution;    
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
      if (currentIteration>=15){sprintf(buffer, "maximum number of iterations reached without convergence. \n"); pcout<<buffer; break;}
      if (current_norm>1/std::pow(tol,2)){sprintf(buffer, "\n norm is too high. \n\n"); pcout<<buffer; break; exit (1);}
      assemble_system();
      current_norm=system_rhs.l2_norm();
      initial_norm=std::max(initial_norm, current_norm);
      res=current_norm/initial_norm;
      sprintf(buffer,"inc:%3u (time:%10.3e, dt:%10.3e), iter:%2u, abs-norm: %10.2e, rel-norm: %10.2e\n", currentIncrement, currentTime, dt,  currentIteration, current_norm, res); pcout<<buffer; 
      if ((currentIteration>1) && ((res<tol) || (current_norm<abs_tol))){sprintf(buffer,"residual converged in %u iterations.\n\n", currentIteration); pcout<<buffer; break;}
      solveIteration();
      U+=dU;  
      ++currentIteration;
    }
    Un=U; 
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

   
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
   
    for (;cell!=endc; ++cell) {
      if (cell->is_locally_owned()) {
	fe_values.reinit(cell);
	cell->get_dof_indices (local_dof_indices);
	//std::cout <<"dof per cell is " << dofs_per_cell << " vertices per cell is " << GeometryInfo<dim>::vertices_per_cell << std::endl; 

	Export[cell].push_back(time);
	for (unsigned int ver=0; ver<1; ++ver) {
	  Export[cell].push_back(cell->vertex(ver)[0]);
	}
	
	for (unsigned int i=0; i<4; ++i) {
	  const unsigned int ci = fe_values.get_fe().system_to_component_index(i).first - 0;   
	  if (ci==0) {
	    //add phi
	    Export[cell].push_back(Un(local_dof_indices[i]));
	  }
	  
	  if (ci==1) {
	    //add vel
	    Export[cell].push_back(Un(local_dof_indices[i]));
	  }

	  if (ci==2) {
	    //add eff press
	    Export[cell].push_back(Un(local_dof_indices[i]));	    
	  }

	  if (ci==3) {
	     //add order                                                                                                                                                                               
            Export[cell].push_back(Un(local_dof_indices[i]));
	  }	  
	}
      }
    }

    typename std::map<typename DoFHandler<dim>::active_cell_iterator, std::vector<double>  >::iterator it = Export.begin();	
    std::ofstream myfile((filename + ".text").c_str(), std::ofstream::out | std::ofstream::app);
    if (myfile.is_open()) {
      myfile<<"Tim "<<"Pos "<<"Phi "<<"Vel "<<"Peff "<<"PhaseF "<<std::endl ;
      for(; it != Export.end(); it++) {
	int counter=0; 
	for(auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
	  myfile<<" "<<std::setprecision (18) <<*it2;
	  counter=counter+1;
	}
	myfile<<std::endl ;
      }
      myfile.close();
    }
    else pcout << "Unable to open file";
    Export.clear();

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
      output_results (0);
    
    //Time stepping
    currentIncrement=0;
    for (currentTime=0; currentTime<totalTime; currentTime+=dt){
      currentIncrement++;
      solve(); 
     int NSTEP=(currentTime/dt);
     if (NSTEP%PSTEPS==0) {
       output_results(currentIncrement);
       output_results_txt(currentIncrement,currentTime);
       //writeSolutionsToFile(Un, tag);	
     }
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
