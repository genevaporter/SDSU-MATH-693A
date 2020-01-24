// CREATED USING THE DEAL.II LIBRARY

/*
 * FOR RELEASE MODE:
 * mkdir Release
 * cd Release
 * cmake -DCMAKE_BUILD_TYPE=Release ..
 * make
 */


/* 
 * dealii: 10,000 steps, 76 seconds
 * me: 10,000, 45 seconds
 * 
 * dealii euler: 582 steps, 12 seconds
 * me euler:  253 steps, 1 seconds
 *
 * currently set to solve deal.II euler
 * 
 */

#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_tools.h>

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <cmath>
#include <chrono>

namespace Step26
{
  using namespace dealii;

  double epsilon = 0.1; // don't change
  double alpha   = 0.1; // don't change
  double beta    = 0.9; // don't change

  double delta   = 9.0;
  double gamma   = 35.0;

  template <int spacedim>
  class HeatEquation
  {
  public:
    HeatEquation(const unsigned degree = 2);
    void run();

  private:
    static constexpr unsigned int dim = spacedim - 1;

    void setup_system();
    void assemble_surface();
    void solve_time_step_U();
    void solve_time_step_V();
    void my_cg_solver_u();
    void my_cg_solver_v();
    void calculate_norm();
    void output_results() const;
    void compute_error() const;
    void adapt_time_step_euler();

    Triangulation<dim, spacedim> 	triangulation;
    FE_Q<dim, spacedim>          	fe;
    DoFHandler<dim, spacedim>    	dof_handler;
    MappingQ<dim, spacedim> 		mapping;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;
    SparseMatrix<double> system_matrix_u, system_matrix_v;

    Vector<double> solution_u, solution_v;
    Vector<double> exact_u, exact_v;
    Vector<double> old_solution_u, old_solution_v;
    Vector<double> system_rhs_u, system_rhs_v;
    Vector<double> u_squared_v;
    Vector<double> normal_u, normal_v;
    Vector<double> forcing_terms;

    double       time;
    double       time_step;
    unsigned int timestep_number;
    unsigned int skip_step;
    double	 tolerance;

    std::ofstream printfile;
  };



  template <int dim>
  class InitialValues_U : public Function<dim> // @suppress("Class has a virtual method and non-virtual destructor")
  {
  public:
	  InitialValues_U () : Function<dim>() {}

	  virtual double value (const Point<dim> &p,
                        const unsigned int component = 0) const;
  };



  template <>
  double InitialValues_U<3>::value (const Point<3> &p,
                                            const unsigned int /*component*/) const
  {
    using numbers::PI;
    double sum = 0;

    for(unsigned int i=1; i<9; i++)
      {
        sum += std::cos(2*PI*i*p[0]*p[1]);
      }

    return (alpha + beta + epsilon*sum);

  }



  template <int dim>
  class InitialValues_V : public Function<dim> // @suppress("Class has a virtual method and non-virtual destructor")
  {
  public:
	  InitialValues_V () : Function<dim>() {}

	  virtual double value (const Point<dim> &p,
                        const unsigned int component = 0) const;
  };



  template <>
  double InitialValues_V<3>::value (const Point<3> &p,
                                            const unsigned int /*component*/) const
  {
    using numbers::PI;
    double sum = 0;

    for(unsigned int i=1; i<9; i++)
      {
        sum += std::cos(2*PI*i*p[2]);
      }    

    return (beta / std::pow(alpha+beta,2) + epsilon*sum);

  }



  template <int dim>
  class RightHandSide : public Function<dim> // @suppress("Class has a virtual method and non-virtual destructor")
  {
  public:
    RightHandSide()
      : Function<dim>()
    {}

    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override;
  };

  template <>
  double RightHandSide<3>::value(const Point<3> &p,
                                   const unsigned int /*component*/) const
  {
    return gamma;
  }


  template <int spacedim>
  HeatEquation<spacedim>::HeatEquation(const unsigned degree)
    : fe(degree),
      dof_handler(triangulation),
      mapping(degree),
      time(0.0),
      time_step(1e-3),
      timestep_number(0),
      skip_step(10),
      tolerance(1e-3)
  {}



  template <int spacedim>
  void HeatEquation<spacedim>::setup_system()
  {
    // Creating a surface mesh that sets boundary IDs to zero
    Triangulation<spacedim> volume_mesh;
    GridGenerator::hyper_ball(volume_mesh);

    std::set<types::boundary_id> boundary_ids;
    boundary_ids.insert(0);

    GridGenerator::extract_boundary_mesh(volume_mesh,
                                         triangulation,
                                         boundary_ids);

    triangulation.set_all_manifold_ids(0);
    triangulation.set_manifold(0, SphericalManifold<dim, spacedim>());

    triangulation.refine_global(4);

    dof_handler.distribute_dofs(fe);

    std::cout << std::endl
              << "===========================================" << std::endl
              << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl
              << "Number of degrees of freedom: " << dof_handler.n_dofs() <<std::endl
	      << "===========================================" << std::endl
              << std::endl;

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    true);
    sparsity_pattern.copy_from(dsp);

    mass_matrix.reinit(sparsity_pattern);
    laplace_matrix.reinit(sparsity_pattern);

    system_matrix_u.reinit(sparsity_pattern);
    system_matrix_v.reinit(sparsity_pattern);

    MatrixCreator::create_mass_matrix(dof_handler,
                                      QGauss<dim>(fe.degree + 1),
                                      mass_matrix);
    MatrixCreator::create_laplace_matrix(dof_handler,
                                         QGauss<dim>(fe.degree + 1),
                                         laplace_matrix);

    solution_u.reinit(dof_handler.n_dofs());
    solution_v.reinit(dof_handler.n_dofs());

    old_solution_u.reinit(dof_handler.n_dofs());
    old_solution_v.reinit(dof_handler.n_dofs());

    system_rhs_u.reinit(dof_handler.n_dofs());
    system_rhs_v.reinit(dof_handler.n_dofs());

    u_squared_v.reinit(dof_handler.n_dofs());
    normal_u.reinit(dof_handler.n_dofs());
    normal_v.reinit(dof_handler.n_dofs());

    forcing_terms.reinit(dof_handler.n_dofs());
  }



  template <int spacedim>
  void HeatEquation<spacedim>::assemble_surface()
  {

    system_matrix_u = 0;
    system_matrix_v = 0;
    system_rhs_u    = 0;
    system_rhs_v    = 0;

    const QGauss<dim>       quadrature_formula(2 * fe.degree);
    FEValues<dim, spacedim> fe_values(mapping,
				      fe,
					quadrature_formula,
					update_values | 
					update_gradients | 
					update_quadrature_points | 
					update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<double>                  rhs_values(n_q_points);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const RightHandSide<spacedim> rhs_function;

    for (const auto &cell : dof_handler.active_cell_iterators())
	  {
	    cell_matrix = 0;
	    cell_rhs    = 0;

	    fe_values.reinit(cell);

	    rhs_function.value_list(fe_values.get_quadrature_points(), rhs_values);

	    for (unsigned int i = 0; i < dofs_per_cell; ++i)
		  for (unsigned int j = 0; j < dofs_per_cell; ++j)
		    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
			  cell_matrix(i, j) += fe_values.shape_grad(i, q_point) *
					       fe_values.shape_grad(j, q_point) *
					       fe_values.JxW(q_point);

	    for (unsigned int i = 0; i < dofs_per_cell; ++i)
		  for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
		    cell_rhs(i) += fe_values.shape_value(i, q_point) *
				   rhs_values[q_point] * fe_values.JxW(q_point);

	    cell->get_dof_indices(local_dof_indices);
	    for (unsigned int i = 0; i < dofs_per_cell; ++i)
		  {
		    for (unsigned int j = 0; j < dofs_per_cell; ++j)
			  system_matrix_u.add(local_dof_indices[i],
					      local_dof_indices[j],
					      cell_matrix(i, j));

		      system_rhs_u(local_dof_indices[i]) += cell_rhs(i);
		      system_rhs_v(local_dof_indices[i]) += cell_rhs(i);
		  }
	  }

    system_matrix_v.copy_from(system_matrix_u);
    system_rhs_v=system_rhs_u;

    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(mapping, dof_handler, 0,
	                                     Functions::ZeroFunction<spacedim>(),
	                                     boundary_values);

    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix_u, solution_u, system_rhs_u, false);

    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix_v, solution_v, system_rhs_v, false);
  }



  // This solves the linear system using deal.ii's CG method
  template <int spacedim>
  void HeatEquation<spacedim>::solve_time_step_U()
  {
    SolverControl solver_control(10000, 1e-10 * system_rhs_u.l2_norm());
    SolverCG<>    cg(solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix_u, 1.0);

    cg.solve(system_matrix_u, solution_u, system_rhs_u, preconditioner);

    constraints.distribute(solution_u);
  }



  // This solves the linear system using deal.ii's CG method
  template <int spacedim>
  void HeatEquation<spacedim>::solve_time_step_V()
  {
    SolverControl solver_control(10000, 1e-10 * system_rhs_v.l2_norm());
    SolverCG<>    cg(solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix_v, 1.0);

    cg.solve(system_matrix_v, solution_v, system_rhs_v, preconditioner);

    constraints.distribute(solution_v);
  }



  template <int spacedim>
  void HeatEquation<spacedim>::my_cg_solver_u()
  {

   // For the linear system Ax=b, we use A = system_matrix and b = system_rhs

    double a=0, b=0, m=0, rnorm=0, temp_scalar=0, solution_norm = 0;
    double solution_tolerance = 1e-10, max_iterations = 10000;
    double n = solution_u.size();

    Vector<double> xvec, rvec, r2vec, pvec, temp_vec;

    // xvec = vector of zeros for initial guess
    xvec = solution_u;
    for(int i=1; i<n; i++)
    {  xvec(i) = 0;  }

    // rvec = A*xvec - b
    rvec = system_rhs_u;
    rvec *= -1;
    r2vec = rvec;

    // pvec = -A*xvec + b
    pvec = rvec;
    pvec *= -1;

    // temp_vec = vector of zeros for storing intermeiate values
    temp_vec = xvec;

    while ((solution_norm > solution_tolerance) && (m < max_iterations))
    {
      // a = (rvec*rvec)/(pvec*A*pvec)
      std::transform(rvec.begin(),
                     rvec.end(),
           	     rvec.begin(),
           	     temp_vec.begin(),
           	     std::multiplies<double>() );
      a = 0;
      for(int i=1; i<n; i++)
      { a += temp_vec(i); }

      system_matrix_u.vmult(temp_vec,pvec);
      std::transform(temp_vec.begin(),
                     temp_vec.end(),
           	     pvec.begin(),
           	     temp_vec.begin(),
           	     std::multiplies<double>() );

      temp_scalar = 0;
      for(int i=1; i<n; i++)
      { temp_scalar += temp_vec(i); }

      a *= 1/temp_scalar;

      // xvec = xvec + a * pvec
      xvec.add(a,pvec);

      // r2vec = rvec + a * A * pvec
      r2vec = rvec;
      system_matrix_u.vmult(temp_vec, pvec);
      r2vec.add(a, temp_vec);

      // calculating norm
      rnorm = 0;
      for (int i = 0; i < n; i++)
      { rnorm += std::pow(r2vec(i),2); }
      rnorm = std::pow(rnorm, 0.5);

      // b = (r2vec*r2vec) / (rvec*rvec)
      std::transform(r2vec.begin(),
                     r2vec.end(),
                     r2vec.begin(),
                     temp_vec.begin(),
                     std::multiplies<double>() );

      b = 0;
      for(int i=1; i<n; i++)
      { b += temp_vec(i); }

      std::transform(rvec.begin(),
                     rvec.end(),
                     rvec.begin(),
                     temp_vec.begin(),
                     std::multiplies<double>() );

      temp_scalar = 0;
      for(int i=1; i<n; i++)
      { temp_scalar += temp_vec(i); }
      b *= 1/temp_scalar;

      // pvec = -r2vec + b * pvec
      pvec *= b;
      pvec.add(-1, r2vec);

      // updating values
      rvec = r2vec;
      m++;
    }

    solution_u = xvec;

  }




  template <int spacedim>
  void HeatEquation<spacedim>::my_cg_solver_v()
  {

   // For the linear system Ax=b, we use A = system_matrix and b = system_rhs

    double a=0, b=0, m=0, rnorm=0, temp_scalar=0, solution_norm = 0;
    double solution_tolerance = 1e-10, max_iterations = 10000;
    double n = solution_v.size();

    Vector<double> xvec, rvec, r2vec, pvec, temp_vec;

    // xvec = vector of zeros for initial guess
    xvec = solution_v;
    for(int i=1; i<n; i++)
    {  xvec(i) = 0;  }

    // rvec = A*xvec - b
    rvec = system_rhs_v;
    rvec *= -1;
    r2vec = rvec;

    // pvec = -A*xvec + b
    pvec = rvec;
    pvec *= -1;

    // temp_vec = vector of zeros for storing intermeiate values
    temp_vec = xvec;

    while ((solution_norm > solution_tolerance) && (m < max_iterations))
    {
      // a = (rvec*rvec)/(pvec*A*pvec)
      std::transform(rvec.begin(),
                     rvec.end(),
           	     rvec.begin(),
           	     temp_vec.begin(),
           	     std::multiplies<double>() );
      a = 0;
      for(int i=1; i<n; i++)
      { a += temp_vec(i); }

      system_matrix_v.vmult(temp_vec,pvec);
      std::transform(temp_vec.begin(),
                     temp_vec.end(),
           	     pvec.begin(),
           	     temp_vec.begin(),
           	     std::multiplies<double>() );

      temp_scalar = 0;
      for(int i=1; i<n; i++)
      { temp_scalar += temp_vec(i); }

      a *= 1/temp_scalar;

      // xvec = xvec + a * pvec
      xvec.add(a,pvec);

      // r2vec = rvec + a * A * pvec
      r2vec = rvec;
      system_matrix_v.vmult(temp_vec, pvec);
      r2vec.add(a, temp_vec);

      // calculating norm
      rnorm = 0;
      for (int i = 0; i < n; i++)
      { rnorm += std::pow(r2vec(i),2); }
      rnorm = std::pow(rnorm, 0.5);

      // b = (r2vec*r2vec) / (rvec*rvec)
      std::transform(r2vec.begin(),
                     r2vec.end(),
                     r2vec.begin(),
                     temp_vec.begin(),
                     std::multiplies<double>() );

      b = 0;
      for(int i=1; i<n; i++)
      { b += temp_vec(i); }

      std::transform(rvec.begin(),
                     rvec.end(),
                     rvec.begin(),
                     temp_vec.begin(),
                     std::multiplies<double>() );

      temp_scalar = 0;
      for(int i=1; i<n; i++)
      { temp_scalar += temp_vec(i); }
      b *= 1/temp_scalar;

      // pvec = -r2vec + b * pvec
      pvec *= b;
      pvec.add(-1, r2vec);

      // updating values
      rvec = r2vec;
      m++;
    }

    solution_v = xvec;

  }



  template<int dim>
  void HeatEquation<dim>::calculate_norm()
  {
    double norm_u = 0;
    double norm_v = 0;
    normal_u = solution_u;
    normal_v = solution_v;
    normal_u.add(-1, old_solution_u);
    normal_v.add(-1, old_solution_v);
    
    for(unsigned int i = 0; i < normal_u.size(); i++)
      {
	  norm_u += std::pow(normal_u[i],2);
	  norm_v += std::pow(normal_v[i],2);
      }

    norm_u = std::pow(norm_u, 0.5);
    norm_v = std::pow(norm_v, 0.5);
    printfile << timestep_number << "," << time << "," 
              << norm_u << "," << norm_v << "," 
              << time_step << "\n";
  }



  template <int spacedim>
  void HeatEquation<spacedim>::adapt_time_step_euler()
  {

    // EULER-MIDPOINT METHOD
    Vector<double> function_f;
    Vector<double> f1u;
    Vector<double> f2u;
    Vector<double> temp_u;

    Vector<double> function_g;
    Vector<double> g1v;
    Vector<double> g2v;
    Vector<double> temp_v;

    Vector<double> temp_uuv;

    function_f.reinit(solution_u.size());
    f1u.reinit(solution_u.size());
    f2u.reinit(solution_u.size());
    temp_u.reinit(solution_u.size());

    function_g.reinit(solution_v.size());
    g1v.reinit(solution_v.size());
    g2v.reinit(solution_v.size());
    temp_v.reinit(solution_v.size());

    temp_uuv.reinit(solution_u.size());
    
    // f1u = M*U + k*gamma*(alpha*F + M*(-U + U^2*V))
    temp_u.add(gamma, u_squared_v);
    temp_u.add(-gamma, solution_u);
    mass_matrix.vmult(function_f, temp_u);
    function_f.add(alpha, forcing_terms);
    mass_matrix.vmult(f1u, solution_u);
    f1u.add(time_step, function_f);

    // g1v = M*V + k*gamma*(beta*F - M*U^2*V)
    mass_matrix.vmult(function_g, u_squared_v);
    function_g *= -gamma;
    function_g.add(beta, forcing_terms);
    mass_matrix.vmult(g1v, solution_v);
    g1v.add(time_step, function_g);    

    // temp_u = M*U + 0.5*k*f1u
    mass_matrix.vmult(temp_u, solution_u);
    temp_u.add(0.5*time_step, f1u);

    // temp_v = M*V + 0.5*k*g1v
    mass_matrix.vmult(temp_v, solution_v);
    temp_v.add(0.5*time_step, g1v);

    // temp_uuv = M*U^2V
    mass_matrix.vmult(temp_uuv, u_squared_v);
    //*/

    // f2u = M*U + k*gamma*(alpha*F - temp_u + temp_uuv)
    function_f = temp_uuv;
    function_f *= gamma;
    function_f.add(-gamma, temp_u);
    function_f.add(alpha, forcing_terms);
    mass_matrix.vmult(f2u, solution_u);
    f2u.add(time_step, function_f);

    // g2v = M*V + k*gamma*(beta*F - temp_uuv)
    function_g = temp_uuv;
    function_g *= -gamma;
    function_g.add(beta, forcing_terms);
    mass_matrix.vmult(g2v, solution_v);
    g2v.add(time_step, function_g);

    Vector<double> e_n_f, e_n_g;
    e_n_f.reinit(solution_u.size());
    e_n_g.reinit(solution_v.size());
    double chi_F, chi_G, chi, rational_f, rational_g;
    rational_f = 0;
    rational_g = 0;

    // rational_f = ||e_n_f||
    // rational_g = ||e_n_g||
    for (unsigned int i=0; i<f2u.size(); i++)
    {
      e_n_f(i) = f2u(i) - f1u(i);
      e_n_g(i) = g2v(i) - g1v(i);
      e_n_f(i) = std::pow(e_n_f(i), 2);
      e_n_g(i) = std::pow(e_n_g(i), 2);
      rational_f += e_n_f(i);
      rational_g += e_n_g(i);
    }
    rational_f = std::pow(rational_f, 0.5);
    rational_g = std::pow(rational_g, 0.5);

    // (chi_F = tolerance / rational_f)^(1/4)
    rational_f = tolerance / rational_f;
    rational_g = tolerance / rational_g;
    chi_F = std::pow(rational_f, 0.25);
    chi_G = std::pow(rational_g, 0.25);

    // chi = min(chi_F, chi_G)
    chi = chi_G;
    if (chi_F < chi_G)
      chi = chi_F;

    if (chi < 0.9)
      chi = 0.9;

    if (chi > 1.1)
      chi = 1.1;

    time_step = chi * time_step;

  }



  template <int spacedim>
  void HeatEquation<spacedim>::run()
  {
    setup_system();
    assemble_surface();

    printfile.open("norm_dealii_euler.csv");
    printfile << "0,0,0,0,0\n";

    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(mapping,
					     dof_handler,
	                                     999,
	                                     Functions::ZeroFunction<spacedim>(),
	                                     boundary_values);

    Vector<double> tmp;

    tmp.reinit(solution_u.size());

    VectorTools::interpolate(dof_handler,
                             InitialValues_U<spacedim>(),
                             old_solution_u);

    VectorTools::interpolate(dof_handler,
                             InitialValues_V<spacedim>(),
                             old_solution_v);

    solution_u = old_solution_u;
    solution_v = old_solution_v;

    std::cout << "t = " << time << std::endl;

    while (time <= 10)
      {

        time += time_step;
        ++timestep_number;

        if (timestep_number % skip_step == 0)
	{
	  std::cout << "t=" << time << std::endl;
	}

	// forcing_terms = gamma*F
	RightHandSide<spacedim> rhs_function;
        VectorTools::create_right_hand_side(mapping,
					    dof_handler,
                                            QGauss<dim>(fe.degree + 1),
                                            rhs_function,
                                            forcing_terms,
					    constraints);

        // u_squared_v = u^2*v
        std::transform(old_solution_u.begin(),
    	  	       old_solution_u.end(),
    	  	       old_solution_u.begin(),
    	  	       tmp.begin(),
    	  	       std::multiplies<double>() );
        std::transform(tmp.begin(),
    		       tmp.end(),
    	  	       old_solution_v.begin(),
    	  	       u_squared_v.begin(),
    	  	       std::multiplies<double>() );

	adapt_time_step_euler();

	///// FOR U /////

	// system_rhs_u = 
	// t*alpha*F + t*gamma*M*(-U_{n-1} + U_{n-1}^2*V_{n-1}) + M*U_{n-1}
        system_rhs_u = forcing_terms;
        system_rhs_u *= time_step*alpha;
        mass_matrix.vmult(tmp, old_solution_u);
        system_rhs_u.add(1-time_step*gamma, tmp);
        mass_matrix.vmult(tmp, u_squared_v);
        system_rhs_u.add(time_step*gamma, tmp);


	// system_matrix_u = M + t*A
	system_matrix_u.copy_from(mass_matrix);
	system_matrix_u.add(time_step, laplace_matrix);
	
	constraints.condense(system_matrix_u, system_rhs_u);

	MatrixTools::apply_boundary_values(boundary_values,
	                                   system_matrix_u,
	                                   solution_u,
	                                   system_rhs_u);


	solve_time_step_U();
	// my_cg_solver_u();

        // u_squared_v = u^2*v
        std::transform(solution_u.begin(),
    	  	       solution_u.end(),
    	  	       solution_u.begin(),
    	  	       tmp.begin(),
    	  	       std::multiplies<double>() );
        std::transform(tmp.begin(),
    		       tmp.end(),
    	  	       old_solution_v.begin(),
    	  	       u_squared_v.begin(),
    	  	       std::multiplies<double>() );


	///// FOR V /////

	// system_rhs_v = t*beta*F + M*V_{n-1} - t*gamma*M*U_{n}^2*V_{n-1}
	system_rhs_v = forcing_terms;
	system_rhs_v *= time_step*beta;
	mass_matrix.vmult(tmp, old_solution_v);
	system_rhs_v.add(1, tmp);
	mass_matrix.vmult(tmp, u_squared_v);
	system_rhs_v.add(-time_step*gamma, tmp);



	// system_matrix_v = M + d*t*A
	system_matrix_v.copy_from(mass_matrix);
	system_matrix_v.add(time_step*delta, laplace_matrix);


	constraints.condense(system_matrix_v, system_rhs_v);

	MatrixTools::apply_boundary_values(boundary_values,
					   system_matrix_v,
					   solution_v,
					   system_rhs_v);

	solve_time_step_V();
	// my_cg_solver_v();
        
	
	if (timestep_number % skip_step == 0)
	{
	  std::cout << time_step << std::endl << std::endl;
	}
	
        calculate_norm();	

        old_solution_u = solution_u;
	old_solution_v = solution_v;

      }

    std::cout << "Ended on time step number " << timestep_number << std::endl;
  }
}


int main()
{
  try
    {
      using namespace dealii;
      using namespace Step26;
      using namespace std::chrono;

      auto start = high_resolution_clock::now();

      HeatEquation<3> heat_equation_solver;
      heat_equation_solver.run();

      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<seconds>(stop - start); 
      std::cout << std::endl << "Time elapsed: " << duration.count() << " seconds." << std::endl;

    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
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
      std::cerr << std::endl
                << std::endl
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
