// Author: Jerome Droniou (jerome.droniou@monash.edu)

#include <fstream>
#include <iomanip>
#include <thread>

#include "sddr-navier-stokes.hpp"
#include <parallel_for.hpp>
#include <max_degrees_quadratures.hpp>
#include <GMpoly_cell.hpp>
#include <IntegrateTripleProduct.hpp>

#include "vtu_writer.hpp"

#include <boost/program_options.hpp>
#include <display_timer.hpp>

#define FORMAT(W)                                                       \
  std::setiosflags(std::ios_base::left) << std::setw(W) << std::setfill(' ')

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Mesh filenames
//------------------------------------------------------------------------------

const std::string mesh_dir = "../../meshes/";
std::string default_mesh = mesh_dir + "Voro-small-0/RF_fmt/voro-2";

// Max number of cells to plot a graph
constexpr size_t max_nb_cells_for_plot = 20000;

//------------------------------------------------------------------------------

int main(int argc, const char* argv[])
{
  // Program options
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Display this help message")
    ("mesh,m", boost::program_options::value<std::string>(), "Set the mesh")
    ("degree,k", boost::program_options::value<size_t>()->default_value(1), "The polynomial degree of the sequence")
    ("pthread,p", boost::program_options::value<bool>()->default_value(true), "Use thread-based parallelism")
    ("solution,s", boost::program_options::value<int>()->default_value(0), "Select the solution")
    ("boundary-conditions,b", boost::program_options::value<std::string>()->default_value("N"), "Boundary conditions (D for essential to fix u x n and p, N for natural to fix curl u x n and u.n)")
    ("reynolds,r", boost::program_options::value<double>()->default_value(1.), "Select the Reynolds number")
    ("pressure-scaling", boost::program_options::value<double>()->default_value(1.), "Select the pressure scaling")
    ("navier-scaling", boost::program_options::value<double>()->default_value(1.), "Apply the scaling to the nonlinear term (curl u) x u")
    ("plot", boost::program_options::value<std::string>(), "Save plot of the solution to the given filename")
    ("solver", boost::program_options::value<std::string>()->default_value("PardisoLU"), "Choice of solver, not case dependent. Options are: PardisoLU, UMFPACK, PaStiXLU, PaStiXLLT, EigenLU, EigenBiCGSTAB (reverts to EigenLU if the selected solver is not available)")
    ("stabilization-parameter,x", boost::program_options::value<double>(), "Set the stabilization parameter")
    ("newton-relax", boost::program_options::value<double>()->default_value(1.), "Initial relaxation parameter in Newton iterations")
    ("newton-tolerance", boost::program_options::value<double>()->default_value(1e-8), "Tolerance to exit Newton iterations")
    ("newton-verbose,v", boost::program_options::value<int>()->default_value(1), "Verbose level of output for Newton iterations");

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  // Display the help options
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  // Select the mesh
  std::string mesh_file = (vm.count("mesh") ? vm["mesh"].as<std::string>() : default_mesh);

  std::cout << "[main] Mesh file: " << mesh_file << std::endl;
  
  // Select the degree 
  size_t K = vm["degree"].as<size_t>();
  std::cout << FORMAT(25) << "[main] Degree:" << K << std::endl;

  // Select the solution
  int solution = (vm.count("solution") ? vm["solution"].as<int>() : 0);
  NavierStokes::ForcingTermType f;
  NavierStokes::VelocityType u;
  NavierStokes::VorticityType omega;
  NavierStokes::PressureType p;
  NavierStokes::PressureGradientType grad_p;
  
  pressure_scaling = vm["pressure-scaling"].as<double>();
  reynolds = vm["reynolds"].as<double>();
  navier_scaling = vm["navier-scaling"].as<double>();

  std::cout << FORMAT(25) << "[main] Reynolds number: " << reynolds << " (scalings: Navier= " << navier_scaling << ", pressure= " << pressure_scaling << ")" << std::endl;

  switch (solution) {
  case 0:
    std::cout << "[main] Test case: trigonometric solution" << std::endl;
    f = trigonometric_f;
    u = trigonometric_u;
    omega = trigonometric_curl_u;
    p = trigonometric_p;
    grad_p = trigonometric_grad_p;
    break;

  case 1:
    std::cout << "[main] Test case: linear solution" << std::endl;
    f = linear_f;
    u = linear_u;
    omega = linear_curl_u;
    p = linear_p;
    grad_p = linear_grad_p;
    break;

  case 2:
    std::cout << "[main] Test case: constant" << std::endl;
    f = constant_f;
    u = constant_u;
    omega = constant_curl_u;
    p = constant_p;
    grad_p = constant_grad_p;
    break;

  case 3:
    std::cout << "[main] Test case: linear velocity, trigonometric pressure" << std::endl;
    f = trigonometric_grad_p;
    u = linear_u;
    omega = linear_curl_u;
    p = trigonometric_p;
    grad_p = trigonometric_grad_p;
    break;

  case 4:
    std::cout << "[main] Test case: vertical flow" << std::endl;
    f = vertical_f;
    u = vertical_u;
    omega = vertical_curl_u;
    p = vertical_p;
    grad_p = vertical_grad_p;
    break;

  case 5:
    std::cout << "[main] Test case: pressure entrance, flux exit" << std::endl;
    f = pressflux_f;
    u = pressflux_u;
    omega = pressflux_curl_u;
    p = pressflux_p;
    grad_p = pressflux_grad_p;
    break;

  default:
    std::cerr << "[main] ERROR: Unknown exact solution" << std::endl;
    exit(1);    
  }

  // Build the mesh and reorder based on the BCs
  MeshBuilder meshbuilder = MeshBuilder(mesh_file);
  std::unique_ptr<Mesh> mesh_ptr = meshbuilder.build_the_mesh();
  std::string bctype = vm["boundary-conditions"].as<std::string>();
  std::cout << FORMAT(25) << "[main] Boundary conditions: " << bctype << std::endl;
  BoundaryConditions BC(bctype, *mesh_ptr.get());
  BC.reorder_vertices("start");
  BC.reorder_edges("start");
  BC.reorder_faces("start");

  // Create DDR core
  bool use_threads = (vm.count("pthread") ? vm["pthread"].as<bool>() : true);
  std::cout << "[main] " << (use_threads ? "Parallel execution" : "Sequential execution") << std:: endl;
  DDRCore ddr_core(*mesh_ptr, K, use_threads);

  // Create NS problem and set some parameters
  NavierStokes ns(ddr_core, BC, use_threads);  
  if(vm.count("stabilization-parameter")) {
    ns.stabilizationParameter() = vm["stabilization-parameter"].as<double>();
  }
  if(vm.count("iterative-solver")) {
    ns.iterativeSolver() = true;
  }

  // Select linear solver
  std::string name_solver = vm["solver"].as<std::string>();
  LinearSolver<NavierStokes::SystemMatrixType> solver(name_solver);
  std::cout << "[main] Solving the system using " << solver.name() << std::endl;

  //--------------------------------------------
  //    NEWTON ITERATIONS
  //--------------------------------------------
  boost::timer::cpu_timer timer;
  timer.start();
  
  // Newton parameters
  double tol = vm["newton-tolerance"].as<double>();      // stopping criterion on residual
  double residual = 1.;     // norm(F(X)-b)/norm(b)
  double prev_residual = residual;
  double relax = vm["newton-relax"].as<double>();   // Initial relaxation parameter
  size_t it_newton = 0;       
  size_t it_newton_MAX = 1000;
  size_t tot_newton = 0;
  size_t good_iter = 0;
  size_t bad_iter = 0;
  size_t nb_restart = 0;
  int verb = vm["newton-verbose"].as<int>();

  // Compute the linear local matrices and RHS (source and natural BCs), used to assemble the Jacobian
  ns.computeLinearComponents(f, u, omega, reynolds);

  // Need a relative value to compute the residual of F(X)=b. The norm of b would be the natural one, but sometimes the
  // RHS can be zero (in case of Dirichlet BCs and no forcing term). In that case we take a simple choice (could be improved
  // by taking the norm of F(X_b), where X_b contains only the BCs, but that requires additional calculations).
  double norm_rhs_b = ns.normGlobalRHS(ns.rhsSourceNaturalBC());
  if (norm_rhs_b < 1e-10) norm_rhs_b = 1.;
      
  // We will increase the Reynolds number progressively
  const Eigen::VectorXd list_reynolds = Eigen::VectorXd::LinSpaced(200,1,10000);
  long int idx_reynolds = 0;
  ns.Reynolds() = std::min(reynolds, list_reynolds[idx_reynolds]);
  if (navier_scaling == 0) ns.Reynolds() = reynolds;

  // Initial state of Newton (includes the Dirichlet BCs), and for restarts (either when we just increased Reynolds, or when we
  //  found a minimal residual in the current Reynolds)
  Eigen::VectorXd uph_n = ns.vectorDirBC(u, p);
   
  Eigen::VectorXd restart_state = Eigen::VectorXd::Zero(ns.nDOFs_up());
  double restart_relax = relax;
  constexpr double residual_threshold = 1e-1;   // Residual threshold for state from which we can restart
  double min_residual = residual_threshold;
  
  if (verb>=2) std::cout << "[Newton] additional outputs formatted as {good it./bad it. | ratio residuals | norm dX}" << std::endl;
      
  // Assemble system (with pseudo-temporal iterations) and analyse pattern
  double reac = (navier_scaling == 0 ? 0. : residual);
  residual = ns.assembleLinearSystem(uph_n, reac) / norm_rhs_b;
  solver.analyzePattern(ns.systemMatrix());

  // Iterations
  do {
    // 
    // Check if tolerance is reached (or if we are in the linear setting).
    //    Then either we exit the iterations, or increase Reynolds and restart
    //
    if ( residual < tol || (navier_scaling == 0 && it_newton>=1) ){
      // Exit iterations if we have the desired Reynolds
      if( std::abs(reynolds - ns.Reynolds()) < 1e-12){
        tot_newton += it_newton;
        break;
      }
      
      // Increase Reynolds to the next number
      idx_reynolds++;
      if (idx_reynolds < list_reynolds.size()){
        ns.Reynolds() = std::min(reynolds, list_reynolds[idx_reynolds]);
      } else {
        ns.Reynolds() = reynolds;
      }
      // Re-assemble the system with the new Reynolds, and restart some Newton parameters
      // The formula for the relaxation parameter comes from trial and error...
      residual = ns.assembleLinearSystem(uph_n, reac) / norm_rhs_b;
      prev_residual = residual;
      relax = (ns.Reynolds() <= 100. ? 1. : 1./(ns.Reynolds()-100));
      good_iter = 0;
      bad_iter = 0;
      tot_newton += it_newton;
      it_newton = 0;
      if (verb>=1) std::cout << "\n [Newton] Increase Reynolds to " << ns.Reynolds() << std::endl << std::endl;
      
      // Update restart, to start again from this solution before increasing Reynolds
      restart_state = uph_n;
      restart_relax = relax;
      min_residual = residual_threshold;
    } // end check residual < tol

    // 
    // Check type of iteration (good/neutral/bad) depending on what the residual history shows.
    //    The relaxation parameter is adjusted depending on the type of iteration
    //
    double norm_dX = 1.;
    double ratio_residuals = residual/prev_residual;
    if (ratio_residuals <= (1. - relax/10.) && residual < 1.){
      //---- Good iteration -----//
      good_iter++;
      uph_n += ns.newtonIncrement(solver, navier_scaling, relax, norm_dX);
      prev_residual = residual;

      // If we are below the current minimal residual, we fix that as the restart state (not sure this is worth keeping actually)
      if (residual < min_residual){
        restart_state = uph_n;
        restart_relax = relax;
        min_residual = residual;
      }

      // If we have enough good iteration, we cancel the previous bad ones, and if we really
      // have quite a bit of them we increase the relaxation parameter.
      if (good_iter >= 5) bad_iter = 0;        
      if (good_iter >= 10 || residual < 1e-3){
        // Increase relaxation depending if the residual is already small or not
        if (relax<1.){
          double factor_relax = (residual < 1e-3 ? 3. : 1.5);
          relax = std::min(1., factor_relax * relax);        
          good_iter = 0;
          if (verb>=1) std::cout << "[Newton] *********** increase relax" << std::endl;
        }
      }

    } else if (ratio_residuals < 1.05) {
      //---- Neutral iteration -----//
      uph_n += ns.newtonIncrement(solver, navier_scaling, relax, norm_dX);
      prev_residual = residual;

    } else {
      //---- Bad iteration ------//
      bad_iter++;
      good_iter = 0;
      relax /= 2;
      
      if (bad_iter >= 5){
        // Too many bad iterations, we restart
        nb_restart++;
        uph_n = restart_state;
        relax = restart_relax / 10.;
        min_residual = residual_threshold;
        bad_iter = 0;
        if (verb >= 1) std::cout << "[Newton] *********** RESTART " << std::endl;
      } else {
        // We update the state
        uph_n += ns.newtonIncrement(solver, navier_scaling, relax, norm_dX);      
      }      
      prev_residual = residual;

    } // end of test good/neutral/bad iteration
    
    it_newton++;

    if (verb>=1){
      std::cout << FORMAT(22) << "[Newton] iteration " << it_newton;
      std::cout.precision(2);
      std::cout << FORMAT(12) << " | residual " << std::scientific << residual;
      std::cout << FORMAT(8) << " | relax " << std::scientific << relax;
      std::cout.precision(0);
      std::cout << FORMAT(10) << " | reynolds " << std::fixed << ns.Reynolds();
      std::cout << std::scientific;
      if (verb >=2){
        std::cout.precision(2);
        std::cout << "     {" << good_iter << "/" << bad_iter << " | " << ratio_residuals << " | " << norm_dX << "}";
      }
      std::cout << std::endl;
      std::cout.precision(6);
    }

  // Assemble system to get ready for next iteration
  double reac = (navier_scaling == 0 ? 0. : residual);
  residual = ns.assembleLinearSystem(uph_n, reac) / norm_rhs_b;

  } while(it_newton < it_newton_MAX);
  if (it_newton >= it_newton_MAX){
    std::cout << "Newton did not converge" << std::endl;
    exit(1);
  }
  
  // For testing
//      saveMarket(ns.systemMatrix(), "A_stokes.mtx");
  
  std::cout << "[Newton] total number of iterations = " << tot_newton << ", final residual " << residual << std::endl;
  timer.stop();
  double t_wall_newton, t_proc_newton;
  std::tie(t_wall_newton, t_proc_newton) = store_times(timer, "[main] Time Newton (wall/proc) ");
  
  //----------------------------------------------
  //    OUTPUTS
  //----------------------------------------------
  
  // Interpolate of exact solution and error vector
  Eigen::VectorXd upI = Eigen::VectorXd::Zero(ns.nDOFs_up());  
  upI.head(ns.sxCurl().dimension()) = ns.sxCurl().interpolate(u, MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
  upI.tail(ns.sxGrad().dimension()) = ns.sxGrad().interpolate(p, MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
  Eigen::VectorXd eph = uph_n - upI;

  std::cout << "[main] Compute errors" << std::endl;
  timer.start();
  // Errors in discrete norms
  std::cout << "[main] Discrete norms and errors:" << std::endl;
  std::vector<StokesNorms> list_norms = ns.computeStokesNorms(std::vector<Eigen::VectorXd> {eph, upI, uph_n});
  std::cout <<"here1"<<std::endl;
  StokesNorms errors = list_norms[0];
  StokesNorms norms = list_norms[1];
  StokesNorms discrete_norms = list_norms[2];
  std::cout << FORMAT(25) << "\tNorms: Hcurl u= " << norms.hcurl_u << "; L2 grad p= " << norms.grad_p <<std::endl;
  double error_hcurl_u = (norms.hcurl_u < 1e-12 ? errors.hcurl_u : errors.hcurl_u / norms.hcurl_u);
  double error_L2grad_p = (norms.grad_p < 1e-12 ? errors.grad_p : errors.grad_p/norms.grad_p);
  std::cout << FORMAT(25) << "\tRel. errors: Hcurl u= " << error_hcurl_u << "; L2 grad p= " << error_L2grad_p << std::endl;
  std::cout << FORMAT(25) << "\tAbs. errors: Hcurl u= " << errors.hcurl_u << "; L2 grad p= " << errors.grad_p << std::endl;
  // Errors in continuous norms
  std::pair<StokesNorms, StokesNorms> c_errors_norms = ns.computeContinuousErrorsNorms(uph_n, u, omega, p, grad_p);
  StokesNorms c_errors = c_errors_norms.first;
  StokesNorms c_norms = c_errors_norms.second;
  double c_error_hcurl_u = (c_norms.hcurl_u < 1e-12 ? c_errors.hcurl_u : c_errors.hcurl_u / c_norms.hcurl_u);
  double c_error_L2grad_p = (c_norms.grad_p < 1e-12 ? c_errors.grad_p : c_errors.grad_p / c_norms.grad_p);
  std::cout << "[main] Continuous relative errors: Hcurl u= " << c_error_hcurl_u << "; L2 grad p= " << c_error_L2grad_p <<std::endl;
  timer.stop();
  double t_wall_norms, t_proc_norms;
  std::tie(t_wall_norms, t_proc_norms) = store_times(timer, "[main] Time errors (wall/proc) ");
  
  // Flux through entrance and exit for BC=M0
  double flux_entrance = 0.;
  double flux_exit = 0.;
  if (BC.id() == "M1"){
    flux_entrance = ns.computeFlux(uph_n.head(ns.sxCurl().dimension()), lower_bottom_corner_x_zero);
    flux_exit = ns.computeFlux(uph_n.head(ns.sxCurl().dimension()), lower_bottom_corner_x_one);
    std::cout << "[main] Fluxes entrance/exit: " << flux_entrance << "/" << flux_exit << std::endl;
  }
 
  std::cout << "[main] Mesh diameter " << mesh_ptr->h_max() << std::endl;

  // --------------------------------------------------------------------------
  //           Creates a .vtu file of the velocity and pressure
  // --------------------------------------------------------------------------
  // Only if we do not have too many cells
  if (vm.count("plot") && mesh_ptr->n_cells() <= max_nb_cells_for_plot) {
    std::cout << "[main] Writing solution to file" << std::endl;
    std::string filename = vm["plot"].as<std::string>();
    VtuWriter plotdata(mesh_ptr.get());
		
		// Exact pressure and velocity at the vertices
    std::vector<double> exact_p_vertex(mesh_ptr->n_vertices(), 0.);
    std::vector<VectorRd> exact_u_vertex(mesh_ptr->n_vertices(), VectorRd::Zero());
    for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
      VectorRd v_coords = mesh_ptr->vertex(iV)->coords();
      exact_p_vertex[iV] = p(v_coords);
      exact_u_vertex[iV] = u(v_coords);
    }
		
    // Approximate pressure and velocity at the vertices
    std::vector<VectorRd> u_vertex = ns.sxCurl().computeVertexValues(uph_n.head(ns.sxCurl().dimension()));
    std::vector<double> p_vertex = ns.sxGrad().computeVertexValues(uph_n.tail(ns.sxGrad().dimension()));
        
    // Plot
    plotdata.write_to_vtu(filename + "_p.vtu", 
                          std::vector<std::vector<double>> {p_vertex, exact_p_vertex},
                          std::vector<std::string>{"pressure", "exact_pressure"});
    plotdata.write_to_vtu(filename + "_u.vtu", 
                          std::vector<std::vector<VectorRd>> {u_vertex, exact_u_vertex},
                          std::vector<std::string>{"velocity", "exact_velocity"});    
  }


  // --------------------------------------------------------------------------
  //                     Creates .txt file with data and results
  // --------------------------------------------------------------------------
  std::ofstream out("results.txt");
  out << std::scientific;
  out << "Solution: " << solution << std::endl;
  out << "Boundary conditions: " << bctype << std::endl;
  out << "Reynolds: " << ns.Reynolds() << std::endl;
  out << "PressureScaling: " << pressure_scaling << std::endl;
  out << "NaviderScaling: " << navier_scaling << std::endl;
  out << "StabilisationParameter: " << ns.stabilizationParameter() << std::endl;
  out << "Mesh: " << mesh_file << std::endl;
  out << "Degree: " << K << std::endl;
  out << "MeshSize: " << mesh_ptr->h_max() << std::endl;
  out << "NbCells: " << mesh_ptr->n_cells() << std::endl;
  out << "NbFaces: " << mesh_ptr->n_faces() << std::endl;
  out << "NbEdges: " << mesh_ptr->n_edges() << std::endl;
  out << "DimSXCurl: " << ns.sxCurl().dimension() << std::endl;
  out << "DimSXGrad: " << ns.sxGrad().dimension() << std::endl;
  out << "SizeSCsystem: " << ns.systemMatrix().rows() << std::endl;
  // Discrete errors and norms
  out << "E_HcurlVel: " << error_hcurl_u << std::endl;
  out << "E_HgradPre: " << (norms.hgrad_p < 1e-12 ? errors.hgrad_p : errors.hgrad_p / norms.hgrad_p) << std::endl;
  out << "E_L2Vel: " << (norms.u < 1e-12 ? errors.u : errors.u/norms.u) << std::endl;
  out << "E_L2CurlVel: " << (norms.curl_u < 1e-12 ? errors.curl_u : errors.curl_u/norms.curl_u) << std::endl;
  out << "E_L2Pre: " << (norms.p < 1e-12 ? errors.p : errors.p/norms.p) << std::endl;
  out << "E_L2GradPre: " << error_L2grad_p << std::endl;
  out << "absE_HcurlVel: " << errors.hcurl_u << std::endl;
  out << "absE_HgradPre: " << errors.hgrad_p << std::endl;
  out << "absE_L2Vel: " << errors.u << std::endl;
  out << "absE_L2CurlVel: " << errors.curl_u << std::endl;
  out << "absE_L2Pre: " << errors.p << std::endl;
  out << "absE_L2GradPre: " << errors.grad_p << std::endl;
  out.precision(16);
  out << "N_L2Vel: " << norms.u << std::endl;
  out << "N_L2CurlVel: " << norms.curl_u << std::endl;
  out << "N_L2Pre: " << norms.p << std::endl;
  out << "N_L2GradPre: " << norms.grad_p << std::endl;
  out << "discN_L2Vel: " << discrete_norms.u << std::endl;
  out << "discN_L2CurlVel: " << discrete_norms.curl_u << std::endl;
  out << "discN_HCurlVel: " << discrete_norms.hcurl_u << std::endl;
  out << "discN_L2Pre: " << discrete_norms.p << std::endl;
  out << "discN_L2GradPre: " << discrete_norms.grad_p << std::endl;
  out << "discN_HGradPre: " << discrete_norms.hgrad_p << std::endl;
  out.precision(6);
  // Continuous errors and norms
  out << "E_cHcurlVel: " << c_error_hcurl_u << std::endl;
  out << "E_cHgradPre: " << (c_norms.hgrad_p < 1e-12 ? c_errors.hgrad_p : c_errors.hgrad_p / c_norms.hgrad_p) << std::endl;
  out << "E_cL2Vel: " << (c_norms.u < 1e-12 ? c_errors.u : c_errors.u/c_norms.u) << std::endl;
  out << "E_cL2CurlVel: " << (c_norms.curl_u < 1e-12 ? c_errors.curl_u : c_errors.curl_u/c_norms.curl_u) << std::endl;
  out << "E_cL2Pre: " << (c_norms.p < 1e-12 ? c_errors.p : c_errors.p/c_norms.p) << std::endl;
  out << "E_cL2GradPre: " << c_error_L2grad_p << std::endl;
  out << "absE_cHcurlVel: " << c_errors.hcurl_u << std::endl;
  out << "absE_cHgradPre: " << c_errors.hgrad_p << std::endl;
  out << "absE_cL2Vel: " << c_errors.u << std::endl;
  out << "absE_cL2CurlVel: " << c_errors.curl_u << std::endl;
  out << "absE_cL2Pre: " << c_errors.p << std::endl;
  out << "absE_cL2GradPre: " << c_errors.grad_p << std::endl;
  out.precision(16);
  out << "N_cL2Vel: " << c_norms.u << std::endl;
  out << "N_cL2CurlVel: " << c_norms.curl_u << std::endl;
  out << "N_cL2Pre: " << c_norms.p << std::endl;
  out << "N_cL2GradPre: " << c_norms.grad_p << std::endl;
  out.precision(6);
  // Newton parameters and times
  out << "TolNewton: " << tol << std::endl;
  out << "ResidualNewton: " << residual << std::endl;
  out << "TotNewton: " << tot_newton << std::endl;
  out << "ItNewtonMax: " << it_newton_MAX << std::endl;
  out << "TwallNewton: " << t_wall_newton << std::endl;  
  out << "TprocNewton: " << t_proc_newton << std::endl;  
  out << std::flush;
  out.close();

  std::cout << "[main] Done" << std::endl;
  return 0;
}

//------------------------------------------------------------------------------
// NavierStokes
//------------------------------------------------------------------------------

NavierStokes::NavierStokes(
               const DDRCore & ddrcore,
               const BoundaryConditions & BC,
               bool use_threads,
               std::ostream & output
               )
  : m_ddrcore(ddrcore),
    m_use_threads(use_threads),
    m_output(output),
    m_ser_pro(ddrcore, use_threads, output),
    m_sxgrad(ddrcore, m_ser_pro, use_threads),
    m_sxcurl(ddrcore, m_ser_pro, use_threads),
    m_sxdiv(ddrcore, use_threads),
    m_nloc_sc_u( m_ser_pro.nDOFs_cells_SXCurl() ),
    m_nloc_sc_p( m_ser_pro.nDOFs_cells_SXGrad() ),
    m_BC(BC),
    m_nDOFs_dir_edge_u(0),
    m_nDOFs_dir_face_u(0),
    m_nDOFs_dir_vertex_p(0),
    m_nDOFs_dir_edge_p(0),
    m_nDOFs_dir_face_p(0),
    m_nUKN_lambda( (m_BC.name() == "Neumann" ? 1 : 0) ),
    m_locSCUKN({}),
    m_locGLUKN({}),
    m_DOFtoSCUKN(-2 * Eigen::VectorXi::Ones(nDOFs_up()+nUKN_lambda())),  // These maps are initialised at -2, but create just below
    m_DOFtoGLUKN(-2 * Eigen::VectorXi::Ones(nDOFs_up()+nUKN_lambda())),
    m_stab_par(0.1),
    m_itsolver(false),
    m_CTuCTu(ddrcore.mesh().n_cells()),
    m_PTuGTp(ddrcore.mesh().n_cells()),
    m_rhs_source_naturalBC(ddrcore.mesh().n_cells()),
    m_massT(ddrcore.mesh().n_cells()),
    m_linear_computed(false),
    m_A(nGLUKN(), nGLUKN()),
    m_b(Eigen::VectorXd::Zero(nGLUKN())),
    m_sc_A(nSCUKN(), nGLUKN()),
    m_sc_b(Eigen::VectorXd::Zero(nSCUKN()))
{
  m_output << "[NavierStokes] Initializing" << std::endl;
  
  // To remove the static condensation, multiply by 0 the vectors in the initialisation of m_nloc_sc_u and m_nloc_sc_p.
  
  // Compute number of Dirichlet DOFs for velocity
  for (Edge * E : ddrcore.mesh().get_b_edges()){
    if (BC.type(*E) == "dir"){
      m_nDOFs_dir_edge_u += m_sxcurl.numLocalDofsEdge(*E);
    }
  }
  for (Face * F : ddrcore.mesh().get_b_faces()){
    if (BC.type(*F) == "dir"){
      m_nDOFs_dir_face_u += m_sxcurl.numLocalDofsFace(*F);
    }
  }
  
  // Compute number of Dirichlet DOFs for pressure
  for (Vertex * V : ddrcore.mesh().get_b_vertices()){
    if (BC.type(*V) == "dir"){
      m_nDOFs_dir_vertex_p += m_sxgrad.numLocalDofsVertex(*V);
    }
  }
  for (Edge * E : ddrcore.mesh().get_b_edges()){
    if (BC.type(*E) == "dir"){
      m_nDOFs_dir_edge_p += m_sxgrad.numLocalDofsEdge(*E);
    }
  }
  for (Face * F : ddrcore.mesh().get_b_faces()){
    if (BC.type(*F) == "dir"){
      m_nDOFs_dir_face_p += m_sxgrad.numLocalDofsFace(*F);
    }
  }

  // We compute the maps DOF->SCUKN and DOF->GLUKN by locating the positions of these unknowns through m_locSCUKN, m_locGLUKN,
  // and using replaceSectionsVector (see BChandlers). 
  // Following the construction of the m_locXXX vectors below is best done by doodling a picture of a global vector of DOFs (u,p,lambda)
  // and locating the dirichlet vertex/edge/face DOFs (at the start of each vertex/edge/face section in the vector) for each of u,p.
  const Eigen::VectorXi MinusOnes = -Eigen::VectorXi::Ones(nDOFs_up()+nUKN_lambda());
  const size_t start_p = m_sxcurl.dimension();

  // ----- DOF->SCUKN -------
  m_locSCUKN.resize(2);
  // Statically condensed cell unknowns for velocity
  m_locSCUKN[0].first = m_sxcurl.nDOFs_edges() + m_sxcurl.nDOFs_faces();
  m_locSCUKN[0].second = nSCUKN_u();
  // Statically condensed cell unknowns for pressure
  m_locSCUKN[1].first = start_p + m_sxgrad.nDOFs_vertices() + m_sxgrad.nDOFs_edges() + m_sxgrad.nDOFs_faces();
  m_locSCUKN[1].second = nSCUKN_p();
  // Create map
  size_t s = nSCUKN();
  m_DOFtoSCUKN = replaceSectionsVector(MinusOnes, Eigen::VectorXi::LinSpaced(s, 0, s-1).eval(), m_locSCUKN);
     
  // ----- DOF->GLUKN -------
  m_locGLUKN.resize(5+m_nUKN_lambda);
  // non-Dirichlet edge DOFs for velocity
  m_locGLUKN[0].first = m_nDOFs_dir_edge_u;
  m_locGLUKN[0].second = m_sxcurl.nDOFs_edges() - m_nDOFs_dir_edge_u; 
  // non-Dirichlet face DOFs for velocity and non-SC cell unknowns
  m_locGLUKN[1].first = m_sxcurl.nDOFs_edges() + m_nDOFs_dir_face_u;
  m_locGLUKN[1].second = (m_sxcurl.nDOFs_faces() - m_nDOFs_dir_face_u) + (m_sxcurl.nDOFs_cells() - nSCUKN_u());
  // non-Dirichlet vertex DOFs for pressure
  m_locGLUKN[2].first = start_p + m_nDOFs_dir_vertex_p;
  m_locGLUKN[2].second = m_sxgrad.nDOFs_vertices() - m_nDOFs_dir_vertex_p;
  // non-Dirichlet edge DOFs for pressure
  m_locGLUKN[3].first = start_p + m_sxgrad.nDOFs_vertices() + m_nDOFs_dir_edge_p;
  m_locGLUKN[3].second = m_sxgrad.nDOFs_edges() - m_nDOFs_dir_edge_p;
  // non-Dirichlet face DOFs and non-SC cell DOFs for pressure
  m_locGLUKN[4].first = start_p + m_sxgrad.nDOFs_vertices() + m_sxgrad.nDOFs_edges() + m_nDOFs_dir_face_p;
  m_locGLUKN[4].second = (m_sxgrad.nDOFs_faces() - m_nDOFs_dir_face_p) + (m_sxgrad.nDOFs_cells() - nSCUKN_p());
  // Lagrange multiplier
  if (m_nUKN_lambda == 1){
    m_locGLUKN[5].first = nDOFs_up();
    m_locGLUKN[5].second = m_nUKN_lambda;
  }

  s = nGLUKN();
  m_DOFtoGLUKN = replaceSectionsVector(MinusOnes, Eigen::VectorXi::LinSpaced(s, 0, s-1).eval(), m_locGLUKN);

}

//------------------------------------------------------------------------------
Eigen::VectorXd NavierStokes::vectorDirBC( const VelocityType & u, const PressureType & p )
{
  // We interpolate u,p, then only extract the values corresponding to the BCs
  // It's not the most efficient, but it is the easiest and this function is only called once anyway

  Eigen::VectorXd interp_u = m_sxcurl.interpolate(u, MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
  Eigen::VectorXd interp_p = m_sxgrad.interpolate(p, MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
  Eigen::VectorXd vectorBC = Eigen::VectorXd::Zero(nDOFs_up());
  
  size_t start_p = m_sxcurl.dimension();
  
  for (Vertex * V : m_ddrcore.mesh().get_b_vertices()){
    if (m_BC.type(*V) == "dir"){
      vectorBC.segment(start_p + m_sxgrad.globalOffset(*V), m_sxgrad.numLocalDofsVertex(*V))
          = interp_p.segment(m_sxgrad.globalOffset(*V), m_sxgrad.numLocalDofsVertex(*V));
    }
  }

  for (Edge * E : m_ddrcore.mesh().get_b_edges()){
    if (m_BC.type(*E) == "dir"){
      vectorBC.segment(m_sxcurl.globalOffset(*E), m_sxcurl.numLocalDofsEdge(*E))
          = interp_u.segment(m_sxcurl.globalOffset(*E), m_sxcurl.numLocalDofsEdge(*E));
      vectorBC.segment(start_p + m_sxgrad.globalOffset(*E), m_sxgrad.numLocalDofsEdge(*E))
          = interp_p.segment(m_sxgrad.globalOffset(*E), m_sxgrad.numLocalDofsEdge(*E));
    }
  }

  for (Face * F : m_ddrcore.mesh().get_b_faces()){
    if (m_BC.type(*F) == "dir"){
      vectorBC.segment(m_sxcurl.globalOffset(*F), m_sxcurl.numLocalDofsFace(*F))
          = interp_u.segment(m_sxcurl.globalOffset(*F), m_sxcurl.numLocalDofsFace(*F));
      vectorBC.segment(start_p + m_sxgrad.globalOffset(*F), m_sxgrad.numLocalDofsFace(*F))
          = interp_p.segment(m_sxgrad.globalOffset(*F), m_sxgrad.numLocalDofsFace(*F));
    }
  }
    
  return vectorBC;
}

//------------------------------------------------------------------------------

void NavierStokes::computeLinearComponents(
                                  const ForcingTermType & f,
                                  const VelocityType & u,
                                  const VorticityType & omega,
                                  const double & reynolds
                                  )
{
  // Interpolate of forcing term in SXCurl
  Eigen::VectorXd interp_f = m_sxcurl.interpolate(f, MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);

  // Compute all local matrices and RHS
  std::function<void(size_t, size_t)> construct_all_linear_components
    = [this, &interp_f, &f, &omega, &u, &reynolds](size_t start, size_t end)->void      {
	      for (size_t iT = start; iT < end; iT++) {
	        // Local contributions without boundary conditions
          std::tie(m_CTuCTu[iT], m_PTuGTp[iT], m_rhs_source_naturalBC[iT], m_massT[iT]) = this->_compute_linear_contributions(iT, interp_f, f);

          // Add natural (Neumann) boundary conditions to rhs
          const Cell & T = *m_ddrcore.mesh().cell(iT);
          if (T.is_boundary()){
            const size_t dqrbc = MAX_DOE_FACE;
            for (Face * pF : T.get_faces()){
              const Face & F = *pF;
              if (m_BC.type(F) == "neu"){
                // Unit normal vector to F pointing out of the domain
                const Cell & TF = *F.cell(0);
                Eigen::Vector3d nF = TF.face_normal(TF.index_face(&F));
                
                QuadratureRule quad_dqrbc_F = generate_quadrature_rule(F, dqrbc);
                
                // Boundary condition on the tangential component of the vorticity
                {
                  FType<Eigen::Vector3d> omega_cross_nF = [&omega, &nF](const Eigen::Vector3d & x) {
                                                                    return omega(x).cross(nF);
                                                                    };
                  // Create integral of omega x nF * PotF as a linear (row) operator on the face DOFs, to extended it afterwards to the cell DOFs.
                  Eigen::RowVectorXd bF =  (1./reynolds) * integrate(omega_cross_nF, 
                                                evaluate_quad<Function>::compute(*m_ddrcore.faceBases(F.global_index()).Polyk2, quad_dqrbc_F),
                                                quad_dqrbc_F).transpose() 
                                           * m_sxcurl.facePotential(F.global_index());
                                           
                  // Add contribution to velocity equations on the face F 
                  m_rhs_source_naturalBC[iT].head(m_sxcurl.dimensionCell(T)) += m_sxcurl.extendOperator(T, F, bF).transpose(); 
                }
                
                // Boundary condition on the normal component of the velocity    
                {
                  FType<double> u_dot_nF = [&u, &nF](const Eigen::Vector3d & x) {
                                                           return u(x).dot(nF);
                                                           };

                  // Create integral of u.nF * PotF as a linear (row) operator on the face DOFs, to extended it afterwards to the cell DOFs.
                  Eigen::RowVectorXd bF = - integrate(u_dot_nF, 
                                                  evaluate_quad<Function>::compute(*m_ddrcore.faceBases(F.global_index()).Polykpo, quad_dqrbc_F),
                                                  quad_dqrbc_F).transpose() 
                                              * m_sxgrad.facePotential(F.global_index());
                  // Add contribution to pressure equations on F
                  m_rhs_source_naturalBC[iT].segment(m_sxcurl.dimensionCell(T), m_sxgrad.dimensionCell(T)) 
                            += m_sxgrad.extendOperator(T, F, bF).transpose();
                }

              } // if m_BC.type(F)=="neu"
            } // Loop over F
 	        } // if T.is_boundary()
 	      } // Loop over T
      };

  parallel_for(m_ddrcore.mesh().n_cells(), construct_all_linear_components, m_use_threads);
  m_linear_computed = true;

}

//------------------------------------------------------------------------------ 

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd>
    NavierStokes::_compute_linear_contributions(
                                              const size_t iT,
                                              const Eigen::VectorXd & interp_f,
                                              const ForcingTermType & f
                                              )
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  size_t dim_sxcurl_T = m_sxcurl.dimensionCell(iT);
  size_t dim_sxgrad_T = m_sxgrad.dimensionCell(iT);
  size_t dim_T = dim_sxcurl_T + dim_sxgrad_T;   
  
  //-------------------
  // Local matrix
  //-------------------

  // Mass matrix for (P^k(T))^3  
  MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_ddrcore.degree());
  Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_2k);

  Eigen::MatrixXd CTuCTu = m_sxdiv.computeL2ProductCurl(iT, m_sxcurl, "both", m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd PTuGTp = m_sxcurl.computeL2ProductGradient(iT, m_sxgrad, "right", m_stab_par, mass_Pk3_T);

  Eigen::MatrixXd MassT = Eigen::MatrixXd::Zero(dim_T+nUKN_lambda(), dim_T+nUKN_lambda());
  MassT.topLeftCorner(dim_sxcurl_T, dim_sxcurl_T) = m_sxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T);
  MassT.block(dim_sxcurl_T, dim_sxcurl_T, dim_sxgrad_T, dim_sxgrad_T) = m_sxgrad.computeL2Product(iT, m_stab_par);
  if (nUKN_lambda()==1) MassT(dim_T, dim_T) = 1.;

  //-----------------------
  // Local source vector
  //-----------------------
  Eigen::VectorXd lT = Eigen::VectorXd::Zero(dim_T+nUKN_lambda());
  lT.head(dim_sxcurl_T) = m_sxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T) * m_sxcurl.restrictCell(iT, interp_f);

  // The lines below replace the line above, using a different source term, that is not pressure-robust:
  // \int f P_curl v. It should only be used for testing purposes
//  QuadratureRule quad_2kp2_T = generate_quadrature_rule(T, 2*m_ddrcore.degree()+2 );
//  lT.head(dim_sxcurl_T) = m_sxcurl.cellPotential(iT).transpose() * 
//                           integrate(f, evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polyk3, quad_2kp2_T), quad_2kp2_T);

  return std::make_tuple(CTuCTu, PTuGTp, lT, MassT);
}

//------------------------------------------------------------------------------ 

Eigen::MatrixXd NavierStokes::_compute_Stokes_matrix(const size_t iT)
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  size_t dim_sxcurl_T = m_sxcurl.dimensionCell(iT);
  size_t dim_sxgrad_T = m_sxgrad.dimensionCell(iT);
  size_t dim_T = dim_sxcurl_T + dim_sxgrad_T;   
  
  Eigen::MatrixXd AT = Eigen::MatrixXd::Zero(dim_T+nUKN_lambda(), dim_T+nUKN_lambda());

  AT.topLeftCorner(dim_sxcurl_T, dim_sxcurl_T) = (1./m_reynolds) * m_CTuCTu[iT];
  AT.block(0, dim_sxcurl_T, dim_sxcurl_T, dim_sxgrad_T ) = m_PTuGTp[iT];
  AT.block(dim_sxcurl_T, 0, dim_sxgrad_T, dim_sxcurl_T ) = -m_PTuGTp[iT].transpose();

  if (nUKN_lambda()==1){
    // Lagrange multiplier: enforce \int P^{k+1}p = 0
    FType<double> cst_fct_one = [](const Eigen::Vector3d &x) -> double { return 1.0; };
    QuadratureRule quad_kpo_T = generate_quadrature_rule(T, m_ddrcore.degree()+1 );
    // intPko is line vector representing \int_T P^{k+1}
    Eigen::RowVectorXd intPkpo = 
        integrate(cst_fct_one, evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polykpo, quad_kpo_T), quad_kpo_T).transpose()
          * m_sxgrad.cellPotential(iT);

    AT.block(dim_T, dim_sxcurl_T, 1, dim_sxgrad_T) = intPkpo;
    AT.block(dim_sxcurl_T, dim_T, dim_sxgrad_T, 1) = -intPkpo.transpose();
  }
      
  return AT;
}

//------------------------------------------------------------------------------

Eigen::VectorXd NavierStokes::newtonIncrement(LinearSolver<NavierStokes::SystemMatrixType> & solver, const double & navier_scaling, const double & relax, double & norm_dX)
{
  // Solves one linear iteration. As dX is an increment and Xn already contains the proper essential BCs, dX is solved with zero essential BCs.
  solver.factorize(systemMatrix());
  Eigen::VectorXd dX_gl = solver.solve(systemVector());

  // Check residual after resolution
    std::cout << "\tResidual linear solver " << (systemMatrix()*dX_gl - systemVector()).lpNorm<Eigen::Infinity>() / (systemVector().lpNorm<Eigen::Infinity>()+1e-8) << std::endl;
    
  // Re-create statically condensed unknowns
  Eigen::VectorXd dX_sc = scVector() + scMatrix() * dX_gl;
 
  // Re-create complete DOFs of dX starting from a zero vector (essential BCs are zero for dX), then replacing the SC and GL unknowns
  Eigen::VectorXd ZeroVec = Eigen::VectorXd::Zero(nDOFs_up()+nUKN_lambda());
  Eigen::VectorXd dX = replaceSectionsVector(replaceSectionsVector(ZeroVec, dX_gl, m_locGLUKN), 
                               dX_sc,
                               m_locSCUKN
                               ).head(nDOFs_up());
                               
  // Actual increment dX is limited based on the relaxation and on the maximum norm increase
  norm_dX = dX.lpNorm<Eigen::Infinity>();
  double fracdX = (navier_scaling == 0 ? 1. : std::min(relax, 1e2/norm_dX));
  
  return fracdX * dX;
}

//------------------------------------------------------------------------------

double NavierStokes::assembleLinearSystem(
                                          const Eigen::VectorXd & Xn,
                                          const double & reac
                                          )
{
  // Store local RHS before static condensation to compute residual
  std::vector<Eigen::VectorXd> loc_rhs(m_ddrcore.mesh().n_cells());
  
  // Assemble all local contributions
  assert(m_linear_computed);
  auto assemble_all = [this, &Xn, &reac, &loc_rhs](
                             size_t start,
                             size_t end,
                             std::list<Eigen::Triplet<double> > * triplets_gl,
                             Eigen::VectorXd * rhs_gl,
                             std::list<Eigen::Triplet<double> > * triplets_sc,
                             Eigen::VectorXd * rhs_sc
                             )->void
                      {
                        for (size_t iT = start; iT < end; iT++) {
                          std::pair<Eigen::MatrixXd, Eigen::VectorXd> locT = _compute_jacobian_rhs(iT, Xn, reac);
                          loc_rhs[iT] = locT.second;
                          this->_assemble_local_contribution(
                                                             iT,
                                                             locT,
                                                             *triplets_gl,
                                                             *rhs_gl,
                                                             *triplets_sc,
                                                             *rhs_sc
                                                             );
                        } // for iT
                      };
                      
  // Assemble the matrix and rhs
  std::tie(m_A, m_b, m_sc_A, m_sc_b) = parallel_assembly_system(m_ddrcore.mesh().n_cells(), this->nGLUKN(), std::make_pair(this->nSCUKN(), this->nGLUKN()), this->nSCUKN(), assemble_all, m_use_threads);

  // Return residual, norm of F(Xn)-b
  return normGlobalRHS(loc_rhs);

}

//------------------------------------------------------------------------------
 
std::pair<Eigen::MatrixXd, Eigen::VectorXd>
    NavierStokes::_compute_jacobian_rhs(
                                        const size_t iT,
                                        const Eigen::VectorXd & Xn,
                                        const double & reac
                                        )
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  // If the nonlinear system is F(X)=b, the Newton iteration is DF(X_n) dX = b-F(X_n) (with dX=X_{n+1}-X_n). 
  // We compute here the local contributions to DF(X_n) and b-F(X_n). 
  
  // Initialise with Stokes matrix and RHS
  Eigen::MatrixXd AT = _compute_Stokes_matrix(iT);
  Eigen::VectorXd lT = m_rhs_source_naturalBC[iT];

  // Potential and Curl of vector Un at current Newton iteration
  Eigen::VectorXd UnT = m_sxcurl.restrict(T, Xn.head(m_sxcurl.dimension()));
  Eigen::MatrixXd PotT = m_sxcurl.cellPotential(iT);
  Eigen::MatrixXd CT = m_sxcurl.cellCurl(iT);
  Eigen::VectorXd PTUnT = PotT * UnT;
  Eigen::VectorXd CTUnT = CT * UnT;

  // lT has been initialised with b, we need to subtract F(X_n). We start by subtracting the contributions of the Stokes terms,
  // by creating the vector (u,p) local, and multiplying it by AT (which only contains the Stokes matrix at this stage)
  size_t dimT = m_sxcurl.dimensionCell(iT) + m_sxgrad.dimensionCell(iT);
  Eigen::VectorXd UPnT = Eigen::VectorXd::Zero(dimT);
  UPnT.head(m_sxcurl.dimensionCell(iT)) = UnT;
  UPnT.tail(m_sxgrad.dimensionCell(iT)) = m_sxgrad.restrict(T, Xn.segment(m_sxcurl.dimension(), m_sxgrad.dimension()));   
  lT.head(dimT) -= AT.topLeftCorner(dimT, dimT) * UPnT;

  // We then adjust the local matrix with the reaction terms corresponding to the pseudo-temporal iterations
  AT += reac * m_massT[iT];

  // We compute the contributions to AT and lT of the nonlinear (Navier) term.
  // The matrices Pk3_dot_Pk3_cross_PTUnT and Pk3_dot_CTUnT_cross_Pk3 represent the trilinear 
  //  form (a,b,c) -> \int_T a . (b x c) on (P^k(T))^3, with 'c' frozen at PT UnT, or 'b' frozen at CT UnT
  Scalar3Tensor Pk3_dot_Pk3_cross_Pk3 = tripleInt(T, *m_ddrcore.cellBases(iT).Polyk3, *m_ddrcore.cellBases(iT).Polyk3);
  size_t dim_Pk3 = m_ddrcore.cellBases(iT).Polyk3->dimension();
  Eigen::MatrixXd Pk3_dot_Pk3_cross_PTUnT = Eigen::MatrixXd::Zero(dim_Pk3, dim_Pk3);
  Eigen::MatrixXd Pk3_dot_CTUnT_cross_Pk3 = Eigen::MatrixXd::Zero(dim_Pk3, dim_Pk3);
  for (size_t i=0; i < dim_Pk3; i++){
    Eigen::MatrixXd triple_sliced_i = slice(Pk3_dot_Pk3_cross_Pk3, 2, i);
    Pk3_dot_Pk3_cross_PTUnT += PTUnT(i) * triple_sliced_i;
    //  Below, minus sign because we use the same slicing w.r.t. last index, but it should be on the 2nd index (- from swapping the terms in x)
    Pk3_dot_CTUnT_cross_Pk3 -= CTUnT(i) * triple_sliced_i;  
  }

  // Adjust matrix and RHS to include the contribution of (curl u) x u to DF(X_n) and F(X_n)
  AT.topLeftCorner(m_sxcurl.dimensionCell(iT), m_sxcurl.dimensionCell(iT)) 
            += navier_scaling * PotT.transpose() * ( Pk3_dot_Pk3_cross_PTUnT * CT + Pk3_dot_CTUnT_cross_Pk3 * PotT);
  lT.head(m_sxcurl.dimensionCell(iT)) -= navier_scaling * PotT.transpose() * Pk3_dot_CTUnT_cross_Pk3 * PTUnT;

  return std::make_pair(AT, lT);
}


//------------------------------------------------------------------------------

void NavierStokes::_assemble_local_contribution(
                                          size_t iT,
                                          const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT,
                                          std::list<Eigen::Triplet<double> > & triplets_gl,
                                          Eigen::VectorXd & rhs_gl,
                                          std::list<Eigen::Triplet<double> > & triplets_sc,
                                          Eigen::VectorXd & rhs_sc
                                          )
{
  // Get information for local static condensation
  LocalStaticCondensation locSC = _compute_static_condensation(iT);

  Eigen::MatrixXd AT_gl, AT_sc;
  Eigen::VectorXd bT_gl, bT_sc;
  std::tie(AT_gl, bT_gl, AT_sc, bT_sc) = locSC.compute(lsT);
    
  // STATICALLY CONDENSED SYSTEM
  // DOFs_gl contain the globally coupled DOFs. We check if they are actual unknowns in the globally coupled system through DOFtoGLUKN
  std::vector<size_t> DOFs_gl = locSC.globalDOFs_gl();
  for (size_t i = 0; i < locSC.dim_gl(); i++){
    int ukn_i = m_DOFtoGLUKN(DOFs_gl[i]);
    if (ukn_i >= 0){
      rhs_gl(ukn_i) += bT_gl(i);
      for (size_t j = 0; j < locSC.dim_gl(); j++){ 
        int ukn_j = m_DOFtoGLUKN(DOFs_gl[j]);
        if (ukn_j >= 0){
          triplets_gl.emplace_back(ukn_i, ukn_j, AT_gl(i,j));
        }
      }
    }
  }

  // RECOVERY OPERATOR
  std::vector<size_t> DOFs_sc = locSC.globalDOFs_sc();
  for (size_t i = 0; i < locSC.dim_sc(); i++){
    int ukn_i = m_DOFtoSCUKN(DOFs_sc[i]);
    if (ukn_i >= 0){
      rhs_sc(ukn_i) += bT_sc(i);
      for (size_t j = 0; j < locSC.dim_gl(); j++){
        int ukn_j = m_DOFtoGLUKN(DOFs_gl[j]);
        if (ukn_j >= 0){
          triplets_sc.emplace_back(ukn_i, ukn_j, AT_sc(i,j));
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

LocalStaticCondensation NavierStokes::_compute_static_condensation(const size_t & iT) const
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  // Dimensions
  size_t dim_u = m_sxcurl.dimensionCell(iT) - m_nloc_sc_u(iT);     // number of local velocity unknowns after SC
  size_t dim_p = m_sxgrad.dimensionCell(iT) - m_nloc_sc_p(iT);     // number of local pressure unknowns after SC
  size_t dim_dofs = dim_u + dim_p + nUKN_lambda();      // nb of dofs remaining after SC (including Lagrange multiplier)
  size_t dim_sc = m_nloc_sc_u(iT) + m_nloc_sc_p(iT);      // nb of SC dofs in total

  // Creation of permutation matrix
  Eigen::MatrixXd Perm = Eigen::MatrixXd::Zero(dim_dofs+dim_sc, dim_dofs+dim_sc);
  Perm.topLeftCorner(dim_u, dim_u) = Eigen::MatrixXd::Identity(dim_u, dim_u);
  Perm.block(dim_u, dim_u + m_nloc_sc_u(iT), dim_p, dim_p) = Eigen::MatrixXd::Identity(dim_p, dim_p);
  if (nUKN_lambda()==1) Perm(dim_u+dim_p, size_t(dim_u + m_nloc_sc_u(iT) + dim_p+m_nloc_sc_p(iT))) = 1.;   // Lagrange multiplier
  Perm.block(dim_dofs, dim_u, m_nloc_sc_u(iT), m_nloc_sc_u(iT)) = Eigen::MatrixXd::Identity(m_nloc_sc_u(iT), m_nloc_sc_u(iT));
  Perm.block(dim_dofs + m_nloc_sc_u(iT), dim_u + m_nloc_sc_u(iT) + dim_p, m_nloc_sc_p(iT), m_nloc_sc_p(iT))
      = Eigen::MatrixXd::Identity(m_nloc_sc_p(iT), m_nloc_sc_p(iT));

  // Creation of global DOFs for system
  //   DOFs_gl and DOFs_sc respectively contains the DOFs of u,p,lambda for globally coupled system and for the recovery system
  //   These DOFs need to be read through DOFtoGLUKN and DOFtoSCUKN when assembling the systems
  std::vector<size_t> DOFs_gl(dim_dofs, 0);
  std::vector<size_t> DOFs_sc(dim_sc, 0);
  auto DOFs_T_sxcurl = m_sxcurl.globalDOFIndices(T);
  auto DOFs_T_sxgrad = m_sxgrad.globalDOFIndices(T);
   
  auto it_DOFs_gl = std::copy(DOFs_T_sxcurl.begin(), DOFs_T_sxcurl.begin()+dim_u, DOFs_gl.begin()); // put skeletal DOFs of u in DOFs_gl
  size_t offset = m_sxcurl.dimension();     // where the DOFs of p start in the global numbering of DOFs
  std::transform(DOFs_T_sxgrad.begin(), DOFs_T_sxgrad.begin()+dim_p, it_DOFs_gl, [&offset](const size_t & index) { return index + offset; });
  if (nUKN_lambda()==1) DOFs_gl[dim_dofs-1] = nDOFs_up(); // Lagrange multiplier
   
  if (dim_sc>0){
    auto it_DOFs_sc = std::copy(DOFs_T_sxcurl.begin()+dim_u, DOFs_T_sxcurl.end(), DOFs_sc.begin());
    std::transform(DOFs_T_sxgrad.begin()+dim_p, DOFs_T_sxgrad.end(), it_DOFs_sc, [&offset](const size_t & index) { return index + offset; });
  }  

  return LocalStaticCondensation(Perm, DOFs_gl, DOFs_sc);
}


//------------------------------------------------------------------------------

double NavierStokes::normGlobalRHS(
                                  const std::vector<Eigen::VectorXd> & locRHS
                                  )
{
  // Create global vector, excluding the DOFs corresponding to Dirichlet values
  Eigen::VectorXd globalRHS = Eigen::VectorXd::Zero(nDOFs_up());
  
  for (size_t iT=0; iT < m_ddrcore.mesh().n_cells(); iT++){
    const Cell & T = *m_ddrcore.mesh().cell(iT);

    std::vector<size_t> ind_vel = m_sxcurl.globalDOFIndices(T);
    for (size_t i=0; i < ind_vel.size(); i++){
      if (m_DOFtoSCUKN[ind_vel[i]]>=0 || m_DOFtoGLUKN[ind_vel[i]]>=0){
        globalRHS(ind_vel[i]) += locRHS[iT](i);
      }
    }

    std::vector<size_t> ind_pres = m_sxgrad.globalDOFIndices(T);
    for (size_t i=0; i < ind_pres.size(); i++){
      if (m_DOFtoSCUKN[m_sxcurl.dimension() + ind_pres[i]]>=0 || m_DOFtoGLUKN[m_sxcurl.dimension() + ind_pres[i]]>=0){
        globalRHS(m_sxcurl.dimension() + ind_pres[i]) += locRHS[iT](ind_vel.size() + i);
      }
    }
  }

  return globalRHS.lpNorm<Eigen::Infinity>();  
}

//------------------------------------------------------------------------------

std::vector<StokesNorms> NavierStokes::computeStokesNorms(const std::vector<Eigen::VectorXd> & list_dofs) const
{
  const size_t ncells = m_ddrcore.mesh().n_cells();
  const size_t nb_vectors = list_dofs.size();
  std::vector<Eigen::VectorXd> local_sqnorm_u(nb_vectors, Eigen::VectorXd::Zero(ncells));
  std::vector<Eigen::VectorXd> local_sqnorm_curl_u(nb_vectors, Eigen::VectorXd::Zero(ncells));
  std::vector<Eigen::VectorXd> local_sqnorm_p(nb_vectors, Eigen::VectorXd::Zero(ncells));
  std::vector<Eigen::VectorXd> local_sqnorm_grad_p(nb_vectors, Eigen::VectorXd::Zero(ncells));

  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [&, this](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *m_ddrcore.mesh().cell(iT);
        // Matrices required to compute the local contributions to the norms
        MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*m_ddrcore.degree() + 2);
        Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_2kp2);
        Eigen::MatrixXd mass_Pkpo_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polykpo, int_mono_2kp2);
        Eigen::MatrixXd L2curl_T = m_sxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T);
        Eigen::MatrixXd L2grad_T = m_sxgrad.computeL2Product(iT, m_stab_par, mass_Pkpo_T);
        Eigen::MatrixXd L2div_CC = m_sxdiv.computeL2ProductCurl(iT, m_sxcurl, "both", m_stab_par, mass_Pk3_T);
        Eigen::MatrixXd GT = m_sxgrad.cellGradient(iT);
        
        for (size_t i=0; i<nb_vectors; i++){
          // Local DOFs for the SXcurl and SXgrad spaces
          Eigen::VectorXd v_curl_T = m_sxcurl.restrict(T, list_dofs[i].head(m_sxcurl.dimension()));
          Eigen::VectorXd v_grad_T = m_sxgrad.restrict(T, list_dofs[i].tail(m_sxgrad.dimension()));

          // Contribution of L2 norms, without any weight (no viscosity)
          local_sqnorm_u[i](iT) = v_curl_T.transpose() * L2curl_T * v_curl_T;
          local_sqnorm_p[i](iT) = v_grad_T.transpose() * L2grad_T * v_grad_T;

          // Contribution of L2 norms of curl and grad
          local_sqnorm_curl_u[i](iT) = v_curl_T.transpose() * L2div_CC * v_curl_T;
          local_sqnorm_grad_p[i](iT) = v_grad_T.transpose() * GT.transpose() * mass_Pk3_T * GT * v_grad_T;
        }
      }
    };
  parallel_for(ncells, compute_local_squarednorms, m_use_threads);

  // Assemble the output
  std::vector<StokesNorms> list_norms;
  list_norms.resize(nb_vectors);
  for (size_t i=0; i<nb_vectors; i++){  
    double sqnorm_u = local_sqnorm_u[i].sum();
    double sqnorm_curl_u = local_sqnorm_curl_u[i].sum();
    double sqnorm_p = local_sqnorm_p[i].sum();
    double sqnorm_grad_p = local_sqnorm_grad_p[i].sum();
    list_norms[i] = StokesNorms(std::sqrt(std::abs(sqnorm_u)), std::sqrt(std::abs(sqnorm_curl_u)), std::sqrt(std::abs(sqnorm_p)), std::sqrt(std::abs(sqnorm_grad_p)));
  }  

  return list_norms;
}

//------------------------------------------------------------------------------

std::pair<StokesNorms, StokesNorms> NavierStokes::computeContinuousErrorsNorms(
            const Eigen::VectorXd & v,
            const VelocityType & u, 
            const VorticityType & curl_u,
            const PressureType & p,
            const PressureGradientType & grad_p
            ) const
{
  // We do not compute the L2 norm/error of p, we just use the L2 norm/error of grad p. De-commenting some lines will result
  // in computing the L2 norm/error of p and using it for the Hgrad error.
  const size_t ncells = m_ddrcore.mesh().n_cells();
  Eigen::VectorXd local_sqerrors_u = Eigen::VectorXd::Zero(ncells);
  Eigen::VectorXd local_sqerrors_curl_u = Eigen::VectorXd::Zero(ncells);
  Eigen::VectorXd local_sqerrors_p = Eigen::VectorXd::Zero(ncells);
  Eigen::VectorXd local_sqerrors_grad_p = Eigen::VectorXd::Zero(ncells);
  Eigen::VectorXd local_sqnorms_u = Eigen::VectorXd::Zero(ncells);
  Eigen::VectorXd local_sqnorms_curl_u = Eigen::VectorXd::Zero(ncells);
  Eigen::VectorXd local_sqnorms_p = Eigen::VectorXd::Zero(ncells);
  Eigen::VectorXd local_sqnorms_grad_p = Eigen::VectorXd::Zero(ncells);

  // Xcurl correspond to the first components of v, SXGrad to the last ones
  Eigen::VectorXd v_curl = v.head(m_sxcurl.dimension());
  Eigen::VectorXd v_grad = v.tail(m_sxgrad.dimension());
  
  std::function<void(size_t, size_t)> compute_local_terms
    = [&](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *m_ddrcore.mesh().cell(iT);
        // Quadrature and evaluation of basis functions
        QuadratureRule qr_high = generate_quadrature_rule(T, 10);
        auto basis_Pk3_quad = evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polyk3, qr_high);
        auto basis_Pkpo_quad = evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polykpo, qr_high);
        
        // Vectors representing (on the bases above) Pcurl u, CT u, Pgrad p and GT p
        Eigen::VectorXd v_curl_T = m_sxcurl.restrict(T, v_curl);
        Eigen::VectorXd Pu = m_sxcurl.cellPotential(iT) * v_curl_T;
        Eigen::VectorXd Cu = m_sxcurl.cellCurl(iT) * v_curl_T;
  
        Eigen::VectorXd v_grad_T = m_sxgrad.restrict(T, v_grad);
        Eigen::VectorXd Pp = m_sxgrad.cellPotential(iT) * v_grad_T;
        Eigen::VectorXd Gp = m_sxgrad.cellGradient(iT) * v_grad_T;

        // Loop over quad nodes
        for (size_t iqn = 0; iqn < qr_high.size(); iqn++){
          // norms
          local_sqnorms_u(iT) += qr_high[iqn].w * u(qr_high[iqn].vector()).squaredNorm();
          local_sqnorms_curl_u(iT) += qr_high[iqn].w * curl_u(qr_high[iqn].vector()).squaredNorm();
          local_sqnorms_p(iT) += qr_high[iqn].w * std::pow( p(qr_high[iqn].vector()), 2);
          local_sqnorms_grad_p(iT) += qr_high[iqn].w * grad_p(qr_high[iqn].vector()).squaredNorm();
          
          // values of potentials/operators at quadrature node
          VectorRd Pu_iqn = Pu(0) * basis_Pk3_quad[0][iqn];
          VectorRd Cu_iqn = Cu(0) * basis_Pk3_quad[0][iqn];
          VectorRd Gp_iqn = Gp(0) * basis_Pk3_quad[0][iqn];
          for (size_t i=1; i<m_ddrcore.cellBases(iT).Polyk3->dimension(); i++){
            Pu_iqn += Pu(i) * basis_Pk3_quad[i][iqn];
            Cu_iqn += Cu(i) * basis_Pk3_quad[i][iqn];
            Gp_iqn += Gp(i) * basis_Pk3_quad[i][iqn];
          }
          double Pp_iqn = Pp(0) * basis_Pkpo_quad[0][iqn];
          for (size_t i=1; i<m_ddrcore.cellBases(iT).Polykpo->dimension(); i++){
            Pp_iqn += Pp(i) * basis_Pkpo_quad[i][iqn];
          }

          // error
          local_sqerrors_u(iT) += qr_high[iqn].w * ( u(qr_high[iqn].vector()) - Pu_iqn ).squaredNorm();
          local_sqerrors_curl_u(iT) += qr_high[iqn].w * ( curl_u(qr_high[iqn].vector()) - Cu_iqn ).squaredNorm();
          local_sqerrors_p(iT) += qr_high[iqn].w * std::pow( p(qr_high[iqn].vector()) - Pp_iqn, 2);
          local_sqerrors_grad_p(iT) += qr_high[iqn].w * ( grad_p(qr_high[iqn].vector()) - Gp_iqn ).squaredNorm();

        }
      }
    };
  parallel_for(ncells, compute_local_terms, m_use_threads);

  double error_u = std::sqrt(std::abs(local_sqerrors_u.sum()));
  double error_curl_u = std::sqrt(std::abs(local_sqerrors_curl_u.sum()));
  double error_p = std::sqrt(std::abs(local_sqerrors_p.sum()));
  double error_grad_p = std::sqrt(std::abs(local_sqerrors_grad_p.sum()));
  StokesNorms errors(error_u, error_curl_u, error_p, error_grad_p);
  
  double norm_u = std::sqrt(std::abs(local_sqnorms_u.sum()));
  double norm_curl_u = std::sqrt(std::abs(local_sqnorms_curl_u.sum()));
  double norm_p = std::sqrt(std::abs(local_sqnorms_p.sum()));
  double norm_grad_p = std::sqrt(std::abs(local_sqnorms_grad_p.sum()));
  StokesNorms norms(norm_u, norm_curl_u, norm_p, norm_grad_p);

  return std::make_pair(errors, norms);

}


//------------------------------------------------------------------------------

double NavierStokes::computeFlux(const Eigen::VectorXd & u, const Hull & surface) const
{
  double flux = 0.;
  
  FType<VectorRd> cst_fct_normal_surf = [&surface](const VectorRd &z) -> VectorRd { return surface.normal; };
  
  // Loop over faces and integrate on those that lie inside the surface
  for (auto F : m_sxcurl.mesh().get_faces()){
    if (surface.is_in(*F)){
      // Cell on the side of F
      const Cell & T = *F->cell(0);
      
      // Integrate Pcurl_T uF.n.
      Eigen::VectorXd uT = m_sxcurl.restrict(T, u);
      QuadratureRule quad_k_F = generate_quadrature_rule(*F, m_sxcurl.degree() );
      Eigen::RowVectorXd int_nF_PkT3 = 
        integrate(cst_fct_normal_surf, evaluate_quad<Function>::compute(*m_sxcurl.cellBases(T).Polyk3, quad_k_F), quad_k_F).transpose();
      
      flux += int_nF_PkT3 * m_sxcurl.cellPotential(T) * uT;
      
    }
  }
  
  return flux;
}
