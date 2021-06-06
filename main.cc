//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Basic framework for Phase Field (Cahn-Hilliard mixed formulation)
//Created May 2018
//authors: rudraa (2018)
//

//deal.II headers
#include "include/headers.h"
//input parameter headers
#include "parameters.h"
//physics headers
#include "include/chemo.h"

//Namespace
namespace phaseField1
{
  using namespace dealii;

  //Initial conditions
  template <int dim>
  class InitalConditions: public Function<dim>{
  public:
    InitalConditions (): Function<dim>(DIMS){std::srand(5);}
    void vector_value (const Point<dim>   &p, Vector<double>   &values) const{
        Assert (values.size() == DIMS, ExcDimensionMismatch (values.size(), DIMS));
	values(0)=0.0; // velocityx
	values(1)=0.0; // velocityy 	
	values(2)=1.0;//porosity	      
	values(3)=0.0; //pressure
	values(4)=1.0; //order     
    }
  };


  /*  
  template <int dim>
  class SolutionBase
  {
  protected:
    static const std::array<Point<dim>, 2> source_centers;
    static const double                    width;
  };
  
  template <>
  const std::array<Point<1>, 2> SolutionBase<1>::source_centers = {
    {Point<1>(-1.0 / 3.0), Point<1>(0.0) }};
  template <>
  const std::array<Point<2>, 2> SolutionBase<2>::source_centers = {
    {Point<2>(-0.5, +0.5), Point<2>(-0.5, -0.5)}};

  template <int dim>
  const double SolutionBase<dim>::width = 1. / 8.;

  template <int dim>
  class Solution : public Function<dim>, protected SolutionBase<dim>
  {
  public:
    Solution()
      : Function<dim>(2)
    {}
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;
  };

  template <int dim>
  double Solution<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    double return_value = 0;
    for (const auto &center : this->source_centers)
      {
        const Tensor<1, dim> x_minus_xi = p - center;
        return_value +=
          std::exp(-x_minus_xi.norm_square() / (this->width * this->width));
      }
    return 0;
  }
  template <int dim>
  Tensor<1, dim> Solution<dim>::gradient(const Point<dim> &p,
                                         const unsigned int) const
  {
    Tensor<1, dim> return_value;
    for (const auto &center : this->source_centers)
      {
        const Tensor<1, dim> x_minus_xi = p - center;
        return_value +=
          (-2 / (this->width * this->width) *
           std::exp(-x_minus_xi.norm_square() / (this->width * this->width)) *
           x_minus_xi);
      }
    return_value[0]=0; return_value[1]=0;
    return return_value;
  }
  template <int dim>
  class RightHandSide : public Function<dim>, protected SolutionBase<dim>
  {
  public:
    RightHandSide()
      : Function<dim>()
    {}
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
  };
  template <int dim>
  double RightHandSide<dim>::value(const Point<dim> &p,
                                   const unsigned int) const
  {
    double return_value = 0;
    for (const auto &center : this->source_centers)
      {
        const Tensor<1, dim> x_minus_xi = p - center;
        return_value +=
          ((2 * dim -
            4 * x_minus_xi.norm_square() / (this->width * this->width)) /
           (this->width * this->width) *
           std::exp(-x_minus_xi.norm_square() / (this->width * this->width)));
        return_value +=
          std::exp(-x_minus_xi.norm_square() / (this->width * this->width));
      }
    return 0;
  }
 
  */

  
  
  
  
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
    void coarse_grid ();
    void output_results (const unsigned int increment);

    void applyBoundaryConditions2(const unsigned int increment);
    void setup_system2 ();
    void assemble_system2 ();
    void solveIteration2 ();
    void solve2 ();
    //void refine_grid ();
    //void coarse_grid ();
    void output_results2 (const unsigned int increment);

    
    MPI_Comm                                  mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim>                             fe;
    DoFHandler<dim>                           dof_handler;
    IndexSet                                  locally_owned_dofs;
    IndexSet                                  locally_relevant_dofs;
    ConstraintMatrix                          constraints,constraintsZero;
    LA::MPI::SparseMatrix                     system_matrix;
    LA::MPI::Vector                           locally_relevant_solution, U, Un, UGhost, UnGhost, dU;
    //LA::MPI::Vector                           U2, Un2, UGhost2, UnGhost2;
    LA::MPI::Vector                           system_rhs;
    ConditionalOStream                        pcout;
    TimerOutput                               computing_timer;

    parallel::distributed::Triangulation<dim> triangulation2;
    FESystem<dim>                             fe2;
    DoFHandler<dim>                           dof_handler2;
    IndexSet                                  locally_owned_dofs2;
    IndexSet                                  locally_relevant_dofs2;
    ConstraintMatrix                          constraints2,constraints2Zero;
    LA::MPI::SparseMatrix                     system_matrix2;
    LA::MPI::Vector                           locally_relevant_solution2, U2, Un2, UGhost2, UnGhost2, dU2;    
    LA::MPI::Vector                           system_rhs2;
    LA::MPI::Vector                           U3;
    

    
    //solution variables
    unsigned int currentIncrement, currentIteration;
    double totalTime, currentTime, dt;
    std::vector<std::string> nodal_solution_names; std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;

    std::vector<std::string> nodal_solution_names2; std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation2;
    
    
  };

  template <int dim>
  phaseField<dim>::phaseField ():
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator,
                   typename Triangulation<dim>::MeshSmoothing
                   (Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening)),
    fe(FE_Q<dim>(FEOrder),DIMS),
    dof_handler (triangulation),

    triangulation2 (mpi_communicator,
                   typename Triangulation<dim>::MeshSmoothing
                   (Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening)),
    fe2(FE_Q<dim>(FEOrder),DIMS),
    dof_handler2 (triangulation2),    
    pcout (std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator)== 0)),
    computing_timer (mpi_communicator, pcout, TimerOutput::summary, TimerOutput::wall_times){
    //solution variables
    dt=TimeStep; totalTime=TotalTime;
    currentIncrement=0; currentTime=0;
    
    //nodal Solution names
  
     
    //nodal Solution names
    for (unsigned int i=0; i<dim; ++i){
      nodal_solution_names.push_back("velocity"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
    }
    nodal_solution_names.push_back("phi"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar); 
    nodal_solution_names.push_back("press"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    nodal_solution_names.push_back("ORDER"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
 
  }
  
  template <int dim>
  phaseField<dim>::~phaseField (){
    dof_handler.clear ();
    dof_handler2.clear ();
  }

  //Apply boundary conditions
  template <int dim>
  void phaseField<dim>::applyBoundaryConditions(const unsigned int increment){
    constraints.clear (); constraintsZero.clear ();  
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraintsZero);
    //Setup boundary conditions    
    std::vector<bool> top (DIMS, false);
    top[2]=true; //porosity
    top[3]=true; //pressure
    
    std::vector<bool> bottom (DIMS, false);
    bottom[0]=true; //vx
    bottom[1]=true; //vy
    bottom[4]=true; //Order

    std::vector<bool> sides(DIMS, false);
    sides[1]=true; //vy

    std::vector<double> valueBottom (DIMS);    
    valueBottom[0]=0.0; //vx
    valueBottom[1]=0.0; //vy
    valueBottom[2]=0.0; //porosity
    valueBottom[3]=0.0; //pressure
    valueBottom[4]=0.0; //order
    
    std::vector<double> valueTop (DIMS);    
    valueTop[0]=0.0; //vx 
    valueTop[1]=0.0; //Vy
    valueTop[2]=0.0; //porosity 
    valueTop[3]=0.0; //pressure
    valueTop[4]=0.0; //order 

    std::vector<double> valueSides (DIMS);    
    valueSides[0]=0; //
    valueSides[1]=0; //Vy
    valueSides[2]=0; //
    valueSides[3]=0; //
    valueSides[4]=0; //


    //bottom
    VectorTools::interpolate_boundary_values (dof_handler, 0, ConstantFunction<dim>(valueBottom), constraints, bottom);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(DIMS), constraintsZero, bottom);
    
    //top
    VectorTools::interpolate_boundary_values (dof_handler, 1, ConstantFunction<dim>(valueTop) , constraints, top);
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(DIMS), constraintsZero, top);
    
    constraints.close ();
    constraintsZero.close ();

  }

  //Apply boundary conditions 2
  template <int dim>
  void phaseField<dim>::applyBoundaryConditions2(const unsigned int increment){
    constraints2.clear (); constraints2Zero.clear ();  
    constraints2.reinit (locally_relevant_dofs2);
    DoFTools::make_hanging_node_constraints (dof_handler2, constraints2Zero);
    //Setup boundary conditions    
    std::vector<bool> top (DIMS, false);
    top[2]=true; //porosity
    top[3]=true; //pressure
    
    std::vector<bool> bottom (DIMS, false);
    bottom[0]=true; //vx
    bottom[1]=true; //vy
    bottom[4]=true; //Order

    std::vector<bool> sides(DIMS, false);
    sides[1]=true; //vy

    std::vector<double> valueBottom (DIMS);    
    valueBottom[0]=0.0; //vx
    valueBottom[1]=0.0; //vy
    valueBottom[2]=0.0; //porosity
    valueBottom[3]=0.0; //pressure
    valueBottom[4]=0.0; //order
    
    std::vector<double> valueTop (DIMS);    
    valueTop[0]=0.0; //vx 
    valueTop[1]=0.0; //Vy
    valueTop[2]=0.0; //porosity 
    valueTop[3]=0.0; //pressure
    valueTop[4]=0.0; //order 

    std::vector<double> valueSides (DIMS);    
    valueSides[0]=0; //
    valueSides[1]=0; //Vy
    valueSides[2]=0; //
    valueSides[3]=0; //
    valueSides[4]=0; //


    //bottom
    VectorTools::interpolate_boundary_values (dof_handler2, 0, ConstantFunction<dim>(valueBottom), constraints2, bottom);
    VectorTools::interpolate_boundary_values (dof_handler2, 0, ZeroFunction<dim>(DIMS), constraints2Zero, bottom);
    
    //top
    VectorTools::interpolate_boundary_values (dof_handler2, 1, ConstantFunction<dim>(valueTop) , constraints2, top);
    VectorTools::interpolate_boundary_values (dof_handler2, 1, ZeroFunction<dim>(DIMS), constraints2Zero, top);
    
    constraints2.close ();
    constraints2Zero.close ();
    
  }

  
  
  
  //Setup
  template <int dim>
  void phaseField<dim>::setup_system (){
    TimerOutput::Scope t(computing_timer, "setup");
    dof_handler.distribute_dofs (fe);
    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler,
                                             locally_relevant_dofs);
    
    locally_relevant_solution.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    //Non-ghost vectors
    system_rhs.reinit (locally_owned_dofs, mpi_communicator);
    U.reinit (locally_owned_dofs, mpi_communicator);
    Un.reinit (locally_owned_dofs, mpi_communicator);
    dU.reinit (locally_owned_dofs, mpi_communicator);
    
    //Ghost vectors
    UGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    UnGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

    //call applyBoundaryConditions to setup constraints matrix needed for generating the sparsity pattern
    applyBoundaryConditions(0);
    
    DynamicSparsityPattern dsp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);

  }


    //Setup 2 
  template <int dim>
  void phaseField<dim>::setup_system2 (){
    TimerOutput::Scope t(computing_timer, "setup");
    dof_handler2.distribute_dofs (fe2);
    locally_owned_dofs2 = dof_handler2.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler2,
                                             locally_relevant_dofs2);
    
    locally_relevant_solution2.reinit (locally_owned_dofs2, locally_relevant_dofs2, mpi_communicator);
    //Non-ghost vectors
    system_rhs2.reinit (locally_owned_dofs2, mpi_communicator);
    U2.reinit (locally_owned_dofs2, mpi_communicator);
    Un2.reinit (locally_owned_dofs2, mpi_communicator);
    dU2.reinit (locally_owned_dofs2, mpi_communicator);
    //U3.reinit (locally_owned_dofs2, mpi_communicator);
    //Ghost vectors
    UGhost2.reinit (locally_owned_dofs2, locally_relevant_dofs2, mpi_communicator);
    UnGhost2.reinit (locally_owned_dofs2, locally_relevant_dofs2, mpi_communicator);

    //call applyBoundaryConditions to setup constraints matrix needed for generating the sparsity pattern
    applyBoundaryConditions2(0);
    
    DynamicSparsityPattern dsp2 (locally_relevant_dofs2);
    DoFTools::make_sparsity_pattern (dof_handler2, dsp2, constraints2, false);
    SparsityTools::distribute_sparsity_pattern (dsp2, dof_handler2.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs2);
    system_matrix2.reinit (locally_owned_dofs2, locally_owned_dofs2, dsp2, mpi_communicator);

  }



  //Assembly
  template <int dim>
  void phaseField<dim>::assemble_system (){
    TimerOutput::Scope t(computing_timer, "assembly");
    system_rhs=0.0; system_matrix=0.0;
    const QGauss<dim>  quadrature_formula(3);
    const QGauss<dim-1>	face_quadrature_formula (2);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |
                             update_JxW_values);
    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors);

   
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    unsigned int n_q_points= fe_values.n_quadrature_points;
  
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned()){
	fe_values.reinit (cell);
	local_matrix = 0; local_rhs = 0; 
	cell->get_dof_indices (local_dof_indices);
	 //AD variables
	Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell); Table<1, double > ULocalConv(dofs_per_cell);
	for (unsigned int i=0; i<dofs_per_cell; ++i){
	  if (std::abs(UGhost(local_dof_indices[i]))<1.0e-16){ULocal[i]=0.0;}
	  else{ULocal[i]=UGhost(local_dof_indices[i]);}
	  // std::cout<<" "<< ULocal[i].val(); //<<std::endl;
	  ULocal[i].diff (i, dofs_per_cell);
	  ULocalConv[i]= UnGhost(local_dof_indices[i]);
	}

	//setup residual vector
	Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); 
	for (unsigned int i=0; i<dofs_per_cell; ++i) {R[i]=0.0;}
	
	//populate residual vector 
	residualForChemo(fe_values, 0, fe_face_values, cell, dt, ULocal, ULocalConv, R, currentTime, totalTime);
	
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
	constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
      }
    system_matrix.compress (VectorOperation::add);
    system_rhs.compress (VectorOperation::add);
    //pcout << "\nK norm is: " << system_matrix.frobenius_norm() << std::endl;
    
  }


 //Assembly2
  template <int dim>
  void phaseField<dim>::assemble_system2 (){
    //mesh two

    system_rhs2=0.0; system_matrix2=0.0;
    const QGauss<dim>  quadrature_formula(3);
    const QGauss<dim-1>	face_quadrature_formula (2);
    FEValues<dim> fe_values (fe2, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |
                             update_JxW_values);
   FEFaceValues<dim> fe_face_values (fe2, face_quadrature_formula, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors);

    const unsigned int   dofs_per_cell2 = fe2.dofs_per_cell;
    FullMatrix<double>   local_matrix2 (dofs_per_cell2, dofs_per_cell2);
    Vector<double>       local_rhs2 (dofs_per_cell2);
    std::vector<unsigned int> local_dof_indices2 (dofs_per_cell2);
    unsigned int n_q_points= fe_values.n_quadrature_points;
  
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler2.begin_active(), endc = dof_handler2.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned()){
	fe_values.reinit (cell);
	local_matrix2 = 0; local_rhs2 = 0; 
	cell->get_dof_indices (local_dof_indices2);
	 //AD variables
	Table<1, Sacado::Fad::DFad<double> > ULocal2(dofs_per_cell2); Table<1, double > ULocalConv2(dofs_per_cell2);
	for (unsigned int i=0; i<dofs_per_cell2; ++i){
	  if (std::abs(UGhost2(local_dof_indices2[i]))<1.0e-16){ULocal2[i]=0.0;}
	  else{ULocal2[i]=UGhost2(local_dof_indices2[i]);}
	  // std::cout<<" "<< ULocal[i].val(); //<<std::endl;
	  ULocal2[i].diff (i, dofs_per_cell2);
	  ULocalConv2[i]= UnGhost2(local_dof_indices2[i]);
	}

	//setup residual vector
	Table<1, Sacado::Fad::DFad<double> > R2(dofs_per_cell2); 
	for (unsigned int i=0; i<dofs_per_cell2; ++i) {R2[i]=0.0;}
	
	//populate residual vector 
	residualForChemo(fe_values, 0, fe_face_values, cell, dt, ULocal2, ULocalConv2, R2, currentTime, totalTime);
	
	//evaluate Residual(R) and Jacobian(R')
	for (unsigned int i=0; i<dofs_per_cell2; ++i) {
	  for (unsigned int j=0; j<dofs_per_cell2; ++j){
	    // R' by AD
	    //local_matrix(i,j)= R[i].fastAccessDx(j);
	    local_matrix2(i,j)= R2[i].dx(j);
	  }
	  //R2
	  local_rhs2(i) = -R2[i].val();
	}
	constraints2.distribute_local_to_global (local_matrix2, local_rhs2, local_dof_indices2, system_matrix2, system_rhs2);
      }
    system_matrix2.compress (VectorOperation::add);
    system_rhs2.compress (VectorOperation::add);
        
  }


  
  /*
  //Solve
  template <int dim>
 void phaseField<dim>::solveIteration(){
    TimerOutput::Scope t(computing_timer, "solve");
    LA::MPI::Vector completely_distributed_solution (locally_owned_dofs, mpi_communicator);
    /*      
    //Iterative solvers from Petsc and Trilinos
    SolverControl solver_control (dof_handler.n_dofs(), 1e-12);
#ifdef USE_PETSC_LA
    LA::SolverGMRES solver(solver_control, mpi_communicator);
#else
    LA::SolverGMRES solver(solver_control);
#endif
    LA::MPI::PreconditionAMG preconditioner;
    LA::MPI::PreconditionAMG::AdditionalData data;
#ifdef USE_PETSC_LA dof_handler.get_fe().n_components()
    //data.symmetric_operator = true;
#else
    // Trilinos defaults are good 
#endif
    preconditioner.initialize(system_matrix, data);
    solver.solve (system_matrix, completely_distributed_solution, system_rhs, preconditioner);
    pcout << "   Solved in " << solver_control.last_step()
          << " iterations." << std::endl;
    */
    //Direct solver MUMPS
  /*
    SolverControl cn;
    PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, completely_distributed_solution, system_rhs);
    if ((currentIteration==0)&&(currentIncrement==1)){
      constraints.distribute (completely_distributed_solution);
    }
    else{
      constraintsZero.distribute (completely_distributed_solution);
    }
    locally_relevant_solution = completely_distributed_solution;
    dU = completely_distributed_solution; 
  }

  

  /*
  //Solve2
  template <int dim>
 void phaseField<dim>::solveIteration2(){
    TimerOutput::Scope t(computing_timer, "solve");
    LA::MPI::Vector completely_distributed_solution (locally_owned_dofs2, mpi_communicator);
    /*      
    //Iterative solvers from Petsc and Trilinos
    SolverControl solver_control (dof_handler.n_dofs(), 1e-12);
#ifdef USE_PETSC_LA
    LA::SolverGMRES solver(solver_control, mpi_communicator);
#else
    LA::SolverGMRES solver(solver_control);
#endif
    LA::MPI::PreconditionAMG preconditioner;
    LA::MPI::PreconditionAMG::AdditionalData data;
#ifdef USE_PETSC_LA dof_handler.get_fe().n_components()
    //data.symmetric_operator = true;
#else
    // Trilinos defaults are good 
#endif
    preconditioner.initialize(system_matrix, data);
    solver.solve (system_matrix, completely_distributed_solution, system_rhs, preconditioner);
    pcout << "   Solved in " << solver_control.last_step()
          << " iterations." << std::endl;
    */
    //Direct solver MUMPS

    /*
    SolverControl cn;
    PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix2, completely_distributed_solution, system_rhs2);

    if ((currentIteration==0)&&(currentIncrement==1)){
      constraints2.distribute (completely_distributed_solution);
    }
    else{
      constraints2Zero.distribute (completely_distributed_solution);
    }
    locally_relevant_solution2 = completely_distributed_solution;
    dU2 = completely_distributed_solution; 
  }
    */
  
  
  
   template <int dim>
  void phaseField<dim>::solveIteration ()
  {
    TimerOutput::Scope t(computing_timer, "solve");
    PETScWrappers::MPI::Vector
      completely_distributed_solution (locally_owned_dofs,mpi_communicator);
    //distributed_incremental_displacement = incremental_displacement;
    SolverControl           solver_control (dof_handler.n_dofs(),
                                            1e-16*system_rhs.l2_norm());
    PETScWrappers::SolverGMRES cg (solver_control,
                                mpi_communicator);
    PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);
    cg.solve (system_matrix,  completely_distributed_solution, system_rhs,
              preconditioner);

    if ((currentIteration==0)&&(currentIncrement==1)){
      constraints.distribute (completely_distributed_solution);
    }
    else{
      constraintsZero.distribute (completely_distributed_solution);
    }
    
    locally_relevant_solution = completely_distributed_solution;
    dU = completely_distributed_solution;
    pcout << "   Solved in " << solver_control.last_step()
          << " iterations." << std::endl;
  }

  
  
  
  
  template <int dim>
  void phaseField<dim>::solveIteration2 ()
  {
    TimerOutput::Scope t(computing_timer, "solve");
    PETScWrappers::MPI::Vector
      completely_distributed_solution (locally_owned_dofs2,mpi_communicator);
    //distributed_incremental_displacement = incremental_displacement;
    SolverControl           solver_control (dof_handler2.n_dofs(),
                                            1e-16*system_rhs.l2_norm());
    PETScWrappers::SolverGMRES cg (solver_control,
				   mpi_communicator);
    PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix2);
    cg.solve (system_matrix2,  completely_distributed_solution, system_rhs2,
              preconditioner);

    if ((currentIteration==0)&&(currentIncrement==1)){
      constraints2.distribute (completely_distributed_solution);
    }
    else{
      constraints2Zero.distribute (completely_distributed_solution);
    }

    locally_relevant_solution2 = completely_distributed_solution;
    dU2 = completely_distributed_solution;
    pcout << "   Solved in " << solver_control.last_step()
          << " iterations." << std::endl;
  }
  
 
  
  
  
  //Solve
  template <int dim>
  void phaseField<dim>::solve(){
    double res=1, tol=1.0e-12, abs_tol=1.0e-14, initial_norm=0, current_norm=0;
    double machineEPS=1.0e-15;
    currentIteration=0;
    char buffer[200];
    pcout << "Solving for mesh 1"<<std::endl;

    while (true){
      if (currentIteration>=25){sprintf(buffer, "maximum number of iterations reached without convergence. \n"); pcout<<buffer; break;exit (1);}
      if (current_norm>1/std::pow(tol,2)){sprintf(buffer, "\n norm is too high. \n\n"); pcout<<buffer; break; exit (1);}
      assemble_system();
      current_norm=system_rhs.l2_norm();
      initial_norm=std::max(initial_norm, current_norm);
      res=current_norm/initial_norm;
      sprintf(buffer,"inc:%3u (time:%10.3e, dt:%10.3e), iter:%2u, abs-norm: %10.2e, rel-norm: %10.2e\n", currentIncrement, currentTime, dt,  currentIteration, current_norm, res); pcout<<buffer; 
      if ((currentIteration>1) && ((res<tol) || (current_norm<abs_tol))){sprintf(buffer,"residual converged in %u iterations.\n\n", currentIteration); pcout<<buffer; break;}
      solveIteration();
      U+=dU; UGhost=U; 
      ++currentIteration;
    }
    Un=U; UnGhost=Un;
  }


   template <int dim>
  void phaseField<dim>::solve2(){
    double res=1, tol=1.0e-12, abs_tol=1.0e-14, initial_norm=0, current_norm=0;
    double machineEPS=1.0e-15;
    currentIteration=0;
    char buffer[200];
    pcout << "Solving for mesh 2"<<std::endl;
    
    while (true){
      if (currentIteration>=25){sprintf(buffer, "maximum number of iterations reached without convergence. \n"); pcout<<buffer; break;exit (1);}
      if (current_norm>1/std::pow(tol,2)){sprintf(buffer, "\n norm is too high. \n\n"); pcout<<buffer; break; exit (1);}
      assemble_system2();
      current_norm=system_rhs2.l2_norm();
      initial_norm=std::max(initial_norm, current_norm);
      res=current_norm/initial_norm;
      sprintf(buffer,"inc:%3u (time:%10.3e, dt:%10.3e), iter:%2u, abs-norm: %10.2e, rel-norm: %10.2e\n", currentIncrement, currentTime, dt,  currentIteration, current_norm, res); pcout<<buffer; 
      if ((currentIteration>1) && ((res<tol) || (current_norm<abs_tol))){sprintf(buffer,"residual converged in %u iterations.\n\n", currentIteration); pcout<<buffer; break;}
      solveIteration2();
      U2+=dU2; UGhost2=U2; 
      ++currentIteration;
    }
    Un2=U2; UnGhost2=Un2;
  }


  

  //Output
  template <int dim>
  void phaseField<dim>::output_results (const unsigned int cycle) {
    TimerOutput::Scope t(computing_timer, "output");
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (UnGhost, nodal_solution_names, DataOut<dim>::type_dof_data, nodal_data_component_interpretation);    

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector (subdomain, "subdomain");
    
    data_out.build_patches ();
    const std::string filename = ("solution-" +
                                  Utilities::int_to_string (cycle, 2) +
                                  "." +
                                  Utilities::int_to_string
                                  (triangulation.locally_owned_subdomain(), 4));
    std::ofstream output ((filename + ".vtu").c_str());
    data_out.write_vtu (output);
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
      std::vector<std::string> filenames;
      for (unsigned int i=0;
	   i<Utilities::MPI::n_mpi_processes(mpi_communicator);
	   ++i)
	filenames.push_back ("solution-" +
			     Utilities::int_to_string (cycle, 2) +
			     "." +
			     Utilities::int_to_string (i, 4) +
			     ".vtu");
      
      std::ofstream master_output (("solution-" +
				    Utilities::int_to_string (cycle, 2) +
				    ".pvtu").c_str());
      data_out.write_pvtu_record (master_output, filenames);
    }
  }


  //Output
  template <int dim>
  void phaseField<dim>::output_results2 (const unsigned int cycle) {
    TimerOutput::Scope t(computing_timer, "output");
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler2);
    data_out.add_data_vector (UnGhost2, nodal_solution_names2, DataOut<dim>::type_dof_data, nodal_data_component_interpretation2);    

    Vector<float> subdomain (triangulation2.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation2.locally_owned_subdomain();
    data_out.add_data_vector (subdomain, "subdomain");
    
    data_out.build_patches ();
    const std::string filename = ("solution2-" +
                                  Utilities::int_to_string (cycle, 2) +
                                  "." +
                                  Utilities::int_to_string
                                  (triangulation2.locally_owned_subdomain(), 4));
    std::ofstream output ((filename + ".vtu").c_str());
    data_out.write_vtu (output);
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
      std::vector<std::string> filenames;
      for (unsigned int i=0;
	   i<Utilities::MPI::n_mpi_processes(mpi_communicator);
	   ++i)
	filenames.push_back ("solution2-" +
			     Utilities::int_to_string (cycle, 2) +
			     "." +
			     Utilities::int_to_string (i, 4) +
			     ".vtu");
      
      std::ofstream master_output (("solution2-" +
				    Utilities::int_to_string (cycle, 2) +
				    ".pvtu").c_str());
      data_out.write_pvtu_record (master_output, filenames);
    }
  }





  

//Adaptive grid refinement
  template <int dim>
  void phaseField<dim>::refine_grid ()  {
    TimerOutput::Scope t(computing_timer, "adaptiveRefinement");
    const QGauss<dim>  quadrature_formula(FEOrder+2);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points);
    unsigned int dofs_per_cell= fe_values.dofs_per_cell;
    unsigned int n_q_points= fe_values.n_quadrature_points;

    std::vector<Vector<double> > quadSolutions;
    for (unsigned int q=0; q<n_q_points; ++q){
      quadSolutions.push_back(dealii::Vector<double>(DIMS)); //2 since there are two degree of freedom per cell
    }

    bool checkForFurtherRefinement=true;
    while (checkForFurtherRefinement) { 
      bool isMeshRefined=false;
      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator t_cell = triangulation.begin_active();
      for (;cell!=endc; ++cell) {
	if (cell->is_locally_owned()) {
	  fe_values.reinit (cell);
	  fe_values.get_function_values(UnGhost, quadSolutions);
	  //limit the maximal and minimal refinement depth of the mesh
	  unsigned int current_level = t_cell->level();
	  // Mark qPoins where refinement is to be done using bool.
	  bool mark_refine = false, mark_refine_solute=false;
	     for (unsigned int q=0; q<n_q_points; ++q) {
	       Point<dim> qPoint=fe_values.quadrature_point(q);	       
	       if ( 1/*quadSolutions[q][1]>-0.9 && quadSolutions[q][1]<0.9*/) {
		 mark_refine = true;
	       }
	           
	     }
	     if ( (mark_refine /*&& mark_refine_solute*/ && (current_level < (maxRefinementLevel)))){
	    cell->set_refine_flag(); isMeshRefined=true; //refine
	  }
	 
	}
	++t_cell;
      }
      
      //check for blocking in MPI
      double checkSum=0.0;
      if (isMeshRefined){checkSum=1.0;}
      checkSum= Utilities::MPI::sum(checkSum, mpi_communicator); //checkSum is greater then 0, then all processors call adative refinementsshown below
      //
      if (checkSum>0.0){
	//execute refinement
	parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> soltrans(dof_handler);
	
	// prepare the triangulation,
	triangulation.prepare_coarsening_and_refinement();
	
	// prepare the SolutionTransfer object for coarsening and refinement
	// and give the solution vector that we intend to interpolate later,
	soltrans.prepare_for_coarsening_and_refinement(UnGhost);
	
	// actually execute the refinement,
	triangulation.execute_coarsening_and_refinement ();
	
	//reset dof's, vectors, matrices, constraints, etc. all on the new mesh.
	setup_system();
	
	// and interpolate all the solutions on the new mesh from the old mesh solution
	soltrans.interpolate(Un);
	U=Un; UGhost=U; UnGhost=Un;
	UGhost.update_ghost_values();
	UnGhost.update_ghost_values();
	//set flag for another check of refinement
	checkForFurtherRefinement=true;
      }
      else{
	checkForFurtherRefinement=false;
      }

    }


    //L2 norm calculation
    FEValues<dim> fe_values2 (fe2, quadrature_formula,
			      update_values    |  update_gradients |
			      update_quadrature_points |
			      update_JxW_values);
    
    unsigned int dofs_per_cell2= fe_values2.dofs_per_cell;
    unsigned int n_q_points2= fe_values2.n_quadrature_points;
    std::vector<Vector<double> > quadSolutions2,quadSolutionsCoarse;
    //std::vector<Tensor<2, dim,double> > quadGrad2(n_q_points2),quadGradCoarse(n_q_points2) ;
    std::vector< std::vector< Tensor< 1, dim, double >>> quadGrad2,quadGradCoarse ;
    
    for (unsigned int q=0; q<n_q_points2; ++q){
      quadSolutions2.push_back(dealii::Vector<double>(DIMS)); //2 since there are two degree of freedom per cell
      quadGrad2.push_back(std::vector< Tensor< 1, dim, double >>(DIMS));
      //pcout <<quadSolutions2[q]<<std::endl;
      
    }

    for (unsigned int q=0; q<n_q_points2; ++q){
      //quadsolution for coarse mesh is already availiable
      quadSolutionsCoarse.push_back(dealii::Vector<double>(DIMS)); //2 since there are two degree of freedom per cell
      quadGradCoarse.push_back(std::vector< Tensor< 1, dim, double >>(DIMS));
      
    }

    

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler2.begin_active(), endc = dof_handler2.end();
    typename parallel::distributed::Triangulation<dim>::active_cell_iterator t_cell = triangulation2.begin_active();
    typename DoFHandler<dim>::active_cell_iterator cell_C = dof_handler.begin_active(), endc_C = dof_handler.end();
    typename parallel::distributed::Triangulation<dim>::active_cell_iterator t_cell_C = triangulation.begin_active();

    double L2ERROR_PHI=0;
    double L2ERROR_Pe=0;
    double L2ERROR_Vx=0;
    double H1SEMI_PHI=0;
    double H1SEMI_Pe=0;
    double H1SEMI_Vx=0;
    double H1_PHI=0;
    double H1_Pe=0;
    double H1_Vx=0;
    //   for (;cell_C!=endc_C; ++cell_C) {
    for (;cell!=endc; ++cell) {
      //double sumT =0 ;
      //double sumPHI=0;
   
      if (cell->is_locally_owned() /*&& cell_C->is_locally_owned()*/) {
	//pcout <<"is this loop even working " <<std::endl;
	fe_values2.reinit(cell);
	fe_values2.get_function_values(UnGhost2, quadSolutions2);
	fe_values2.get_function_gradients(UnGhost2, quadGrad2);	

	fe_values2.reinit(cell);
	fe_values2.get_function_values(UnGhost, quadSolutionsCoarse);
	fe_values2.get_function_gradients(UnGhost, quadGradCoarse);
	dealii::Table<1,double>  PHI(n_q_points2), Pe(n_q_points2),Vx(n_q_points2) ; 
	dealii::Table<1,double>  PHI_g(n_q_points2), Pe_g(n_q_points2),Vx_g(n_q_points2) ;

	/*
	for (unsigned int q=0; q<n_q_points2; ++q) {
	  for (unsigned int i=0; i<dofs_per_cell2; ++i) {
	    //for (unsigned int j=0; j<dim; j++) {T_j[q][j]=0.0;  PHI_j[q][j]=0.0; TC_j[q][j]=0.0;  PHIC_j[q][j]=0.0;}

	    const unsigned int ck = fe_values2.get_fe().system_to_component_index(i).first - 0;
	    for (unsigned int j=0; j<dim; j++) {
	    if (ck==0) {
	      T_j[q][j]+=fe_values2.shape_grad(i, q)[j]*quadSolutions2[q][0];
	      TC_j[q][j]+=fe_values2.shape_grad(i, q)[j]*quadSolutionsCoarse[q][0];
	    }
	    
	    if (ck==1) {
	      PHI_j[q][j]+=fe_values2.shape_grad(i, q)[j]*quadSolutions2[q][1];
	      PHIC_j[q][j]+=fe_values2.shape_grad(i, q)[j]*quadSolutionsCoarse[q][1];
	    }
	    
	    }
	    
	  }
	}
	*/
	
	
	for (unsigned int i=0; i<dofs_per_cell2; ++i) {
	  for (unsigned int q=0; q<n_q_points2; ++q) {
	    //   Point<dim> qPoint2=fe_values2.quadrature_point(q);	       
	    PHI[q]=0.0; Pe[q]=0.0;Vx[q]=0.0;PHI_g[q]=0;Pe_g[q]=0;Vx_g[q]=0.0;
	    
	    const unsigned int ck = fe_values2.get_fe().system_to_component_index(i).first - 0;
	    	    
	    if (ck==0) {
	      Vx[q]+=quadSolutions2[q][0];
	      Vx[q]+=-quadSolutionsCoarse[q][0];
	      Vx_g[q]=pow(quadGrad2[q][0][0]-quadGradCoarse[q][0][0],2)+pow(quadGrad2[q][0][1]-quadGradCoarse[q][0][1],2);
	      L2ERROR_Vx +=(pow(Vx[q],2))*fe_values2.JxW(q);
	      //H1SEMI_PHI+=(pow(PHI_j[q][0]-PHIC_j[q][0],2) + pow(PHI_j[q][1]-PHIC_j[q][1],2) )*fe_values2.JxW(q);
	      H1SEMI_Vx+=(Vx_g[q])*fe_values2.JxW(q);		    
	      //H1_PHI+=(pow(PHI[q],2))*fe_values2.JxW(q)+(pow(PHI_j[q][0]-PHIC_j[q][0],2) + pow(PHI_j[q][1]-PHIC_j[q][1],2) )*fe_values2.JxW(q);
	      H1_Vx+=(pow(Vx[q],2))*fe_values2.JxW(q)+(Vx_g[q])*fe_values2.JxW(q);			    	      
	    }


	    if (ck==2) {
	      PHI[q]+=quadSolutions2[q][2];
	      PHI[q]+=-quadSolutionsCoarse[q][2];
	      PHI_g[q]=pow(quadGrad2[q][2][0]-quadGradCoarse[q][2][0],2)+pow(quadGrad2[q][2][1]-quadGradCoarse[q][2][1],2);
	      
	      L2ERROR_PHI += (pow(PHI[q],2))*fe_values2.JxW(q);
	      //H1SEMI_T+=(pow(T_j[q][0]-TC_j[q][0],2) + pow(T_j[q][1]-TC_j[q][1],2) )*fe_values2.JxW(q);
	      H1SEMI_PHI+=(PHI_g[q])*fe_values2.JxW(q);
	      // H1_T+=(pow(T[q],2))*fe_values2.JxW(q)+(pow(T_j[q][0]-TC_j[q][0],2) + pow(T_j[q][1]-TC_j[q][1],2) )*fe_values2.JxW(q); 
	      H1_PHI+=(pow(PHI[q],2))*fe_values2.JxW(q)+(PHI_g[q])*fe_values2.JxW(q);			    
	    }

	    if (ck==3) {
	      Pe[q]+=quadSolutions2[q][3];
	      Pe[q]+=-quadSolutionsCoarse[q][3];
	      Pe_g[q]=pow(quadGrad2[q][3][0]-quadGradCoarse[q][3][0],2)+pow(quadGrad2[q][3][1]-quadGradCoarse[q][3][1],2);
	      L2ERROR_Pe +=(pow(Pe[q],2))*fe_values2.JxW(q);
	      //H1SEMI_PHI+=(pow(PHI_j[q][0]-PHIC_j[q][0],2) + pow(PHI_j[q][1]-PHIC_j[q][1],2) )*fe_values2.JxW(q);
	      H1SEMI_Pe+=(Pe_g[q])*fe_values2.JxW(q);		    
	      //H1_PHI+=(pow(PHI[q],2))*fe_values2.JxW(q)+(pow(PHI_j[q][0]-PHIC_j[q][0],2) + pow(PHI_j[q][1]-PHIC_j[q][1],2) )*fe_values2.JxW(q);
	      H1_Pe+=(pow(Pe[q],2))*fe_values2.JxW(q)+(Pe_g[q])*fe_values2.JxW(q);			    
	      
	    }
    
	  }
	  
	}
	
	
      }
      ++t_cell;

    }
    
    // ++t_cell_C;
    
    //  }
    
    L2ERROR_PHI= Utilities::MPI::sum(L2ERROR_PHI, mpi_communicator);
    L2ERROR_Pe= Utilities::MPI::sum(L2ERROR_Pe, mpi_communicator);
    L2ERROR_Vx= Utilities::MPI::sum(L2ERROR_Vx, mpi_communicator);
    
    H1SEMI_PHI= Utilities::MPI::sum(H1SEMI_PHI, mpi_communicator);
    H1SEMI_Pe= Utilities::MPI::sum(H1SEMI_Pe, mpi_communicator);
    H1SEMI_Vx= Utilities::MPI::sum(H1SEMI_Vx, mpi_communicator);
    
    H1_PHI= Utilities::MPI::sum(H1_PHI, mpi_communicator);
    H1_Pe= Utilities::MPI::sum(H1_Pe, mpi_communicator);
    H1_Vx= Utilities::MPI::sum(H1_Vx, mpi_communicator);
    
    pcout << "  L2_Error for Phi is :       "
	  << std::sqrt(L2ERROR_PHI)<<std::endl;
    
    
    pcout << "  L2_Error for Pe is :       "
	  << std::sqrt(L2ERROR_Pe)<<std::endl;

    pcout << "  L2_Error for Vx is :       "
	  << std::sqrt(L2ERROR_Vx)<<std::endl;

    
    pcout << "  H1SEMI_NORM for Phi is :       "
	  << std::sqrt(H1SEMI_PHI)<<std::endl;
    
    
    pcout << "  H1SEMI_NORM  for Pe is :       "
	  << std::sqrt(H1SEMI_Pe)<<std::endl;

    pcout << "  H1SEMI_NORM  for Vx is :       "
	  << std::sqrt(H1SEMI_Vx)<<std::endl;


    pcout << "  H1_NORM for Phi is :       "
          << std::sqrt(H1_PHI)<<std::endl;

    pcout << "  H1_NORM  for Pe is :       "
          << std::sqrt(H1_Pe)<<std::endl;

    pcout << "  H1_NORM  for Vx is :       "
          << std::sqrt(H1_Vx)<<std::endl;

    
  }
  
  //Adaptive grid refinement
  template <int dim>
  void phaseField<dim>::coarse_grid () {
    TimerOutput::Scope t(computing_timer, "adaptiveRefinement");
    const QGauss<dim>  quadrature_formula(FEOrder+2);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points);
    unsigned int dofs_per_cell= fe_values.dofs_per_cell;
    unsigned int n_q_points= fe_values.n_quadrature_points;


    /*
    char buffer[200];
    sprintf(buffer, "laser source at: (%7.3e, %7.3e,%7.3e)\n",laserLocationX,laserLocationY,laserLocationZ);
    pcout << buffer;
    */

    std::vector<Vector<double> > quadSolutions;
    for (unsigned int q=0; q<n_q_points; ++q){
      quadSolutions.push_back(dealii::Vector<double>(DIMS)); //2 since there are two degree of freedom per cell
    }

    bool checkForFurtherRefinement=true;
    while (checkForFurtherRefinement) { 
      bool isMeshRefined=false;
      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator t_cell = triangulation.begin_active();
      for (;cell!=endc; ++cell) {
	if (cell->is_locally_owned()) {
	  fe_values.reinit (cell);
	  fe_values.get_function_values(UnGhost, quadSolutions);
	  
	  //limit the maximal and minimal refinement depth of the mesh
	  unsigned int current_level = t_cell->level();
	  
	  // Mark qPoins where refinement is to be done using bool.
	  bool mark_refine = false, mark_refine_solute=false;
	     for (unsigned int q=0; q<n_q_points; ++q) {
	       Point<dim> qPoint=fe_values.quadrature_point(q);	       
	       if (1/*quadSolutions[q][1]>-0.9 && quadSolutions[q][1]<0.9*/) {
		 mark_refine = true;
	       }
	           
	     }
	   
	   
	     if ( (mark_refine /*&& mark_refine_solute*/ && (current_level > (minRefinementLevel)))){
	    cell->set_coarsen_flag();; isMeshRefined=true; //refine
	  }
	 
	}
	++t_cell;
      }
      
      //check for blocking in MPI
      double checkSum=0.0;
      if (isMeshRefined){checkSum=1.0;}
      checkSum= Utilities::MPI::sum(checkSum, mpi_communicator); //checkSum is greater then 0, then all processors call adative refinement shown below
      //
      if (checkSum>0.0){
	//execute refinement
	parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> soltrans(dof_handler);
	
	// prepare the triangulation,
	triangulation.prepare_coarsening_and_refinement();
	
	// prepare the SolutionTransfer object for coarsening and refinement
	// and give the solution vector that we intend to interpolate later,
	soltrans.prepare_for_coarsening_and_refinement(UnGhost);
	
	// actually execute the refinement,
	triangulation.execute_coarsening_and_refinement ();
	
	//reset dof's, vectors, matrices, constraints, etc. all on the new mesh.
	setup_system();
	
	// and interpolate all the solutions on the new mesh from the old mesh solution
	soltrans.interpolate(Un);
	U=Un; UGhost=U; UnGhost=Un;
	UGhost.update_ghost_values();
	UnGhost.update_ghost_values();
	//set flag for another check of refinement
	checkForFurtherRefinement=true;
      }
      else{
	checkForFurtherRefinement=false;
      }
    }

    
  }

  
  
  
  //Solve problem
  template <int dim>
  void phaseField<dim>::run (){
    //setup problem geometry and mesh

    std::vector<unsigned int> numRepetitions;
    numRepetitions.push_back(sub_x); // x refinement
    numRepetitions.push_back(sub_y); // y refinement
    //numRepetitions.push_back(ZSubRf); // z refinement
    
    Point<dim> p1 (0,0);
    Point<dim> p2 (problemWidth,problemHeight);
    GridGenerator::subdivided_hyper_rectangle (triangulation, numRepetitions, p1, p2, true);
    GridGenerator::subdivided_hyper_rectangle (triangulation2, numRepetitions, p1, p2, true);
    triangulation.refine_global (globalRefinementFactor);
    setup_system (); //inital set up
    
    triangulation2.refine_global (globalRefinementFactor2);
    setup_system2 (); //inital set up
     
    
    pcout << "   Number of active cells in mesh 1:       "
	  << triangulation.n_global_active_cells()
	  << std::endl
	  << "   Number of degrees of freedom: "
	  << dof_handler.n_dofs()
	  << std::endl;

    pcout << "   Number of active cells in mesh 2:       "
	  << triangulation2.n_global_active_cells()
	  << std::endl
	  << "   Number of degrees of freedom: "
	  << dof_handler2.n_dofs()
	  << std::endl;

    
    //setup initial conditions
    VectorTools::interpolate(dof_handler, InitalConditions<dim>(), U); Un=U;
    VectorTools::interpolate(dof_handler2, InitalConditions<dim>(), U2); Un2=U2;

    
    //sync ghost vectors to non-ghost vectors
    UGhost=U;  UnGhost=Un;
    output_results (0);

    //sync ghost vectors to non-ghost vectors
    UGhost2=U2;  UnGhost2=Un2;
    
   
    //Time stepping
    currentIncrement=0;
    for (currentTime=0; currentTime<totalTime; currentTime+=dt){
      currentIncrement++;
      solve();
      solve2();
      int NSTEP=(currentTime/dt);
      if (NSTEP%PSTEPS==0)
	{
	output_results(currentIncrement);
	//output_results2(currentIncrement);
	refine_grid(); //adative refinement

	/*pcout << "   Number of active cells:       "
	      << triangulation.n_global_active_cells()
	      << std::endl; */

	coarse_grid(); //adative refinement

	/*pcout << "   Number of active cells:       "
	      << triangulation.n_global_active_cells()
	      << std::endl; */

	
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
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      phaseField<2> problem;
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
