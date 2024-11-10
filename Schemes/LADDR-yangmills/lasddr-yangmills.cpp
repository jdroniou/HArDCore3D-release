// Authors: Jia Jia Qian (jia.qian@monash.edu)

#include <fstream>
#include <iomanip>
#include <thread>

#include "lasddr-yangmills.hpp"
#include <max_degrees_quadratures.hpp>
#include <GMpoly_cell.hpp>
#include <IntegrateTripleProduct.hpp>
#include <unsupported/Eigen/KroneckerProduct>

#include <boost/program_options.hpp>
#include <display_timer.hpp>

#ifdef WITH_UMFPACK
  #include <Eigen/UmfPackSupport>
#endif

#ifdef WITH_MKL
  #include <Eigen/PardisoSupport>
  #include <mkl.h>
#endif

#ifdef WITH_SPECTRA
  #include <Spectra/SymEigsSolver.h>
  #include <Spectra/MatOp/SparseSymShiftSolve.h>
  #include <Spectra/MatOp/SparseSymMatProd.h>
  #include <Spectra/SymEigsShiftSolver.h>
#endif

#define FORMAT(W)                                                       \
  std::setiosflags(std::ios_base::left) << std::setw(W) << std::setfill(' ')

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Mesh filenames
//------------------------------------------------------------------------------

const std::string mesh_dir = "../../meshes/";
std::string default_mesh = mesh_dir + "Voro-small-0/RF_fmt/voro-2";

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
    ("solution,s", boost::program_options::value<int>()->default_value(1), "Select the solution")
    ("export-matrix,e", "Export matrix to Matrix Market format")
    ("iterative-solver,i", boost::program_options::value<bool>()->default_value(false), "Use iterative linear solver")
    ("stabilization-parameter,x", boost::program_options::value<double>()->default_value(.1), "Set the stabilization parameter")
    ("iterations,n", boost::program_options::value<int>()->default_value(-5), "Set the number of iterations")
    ("stopping-value,v", boost::program_options::value<double>()->default_value(1e-6), "Set the value to stop each iteration")
    ("final-time,f", boost::program_options::value<double>()->default_value(1.), "Set the final time")
    ("nonlinear-scaling,l", boost::program_options::value<double>()->default_value(1.), "Degree of nonlinearity")
    ("theta,t", boost::program_options::value<double>()->default_value(1.), "Theta")
    ("condition-numbers,c", boost::program_options::value<bool>()->default_value(false), "Calculate condition numbers (requires the Spectra library)")
    ("exact-initial-conditions,q", boost::program_options::value<bool>()->default_value(false), "Use constrained initial conditions")
    ("bracket,b", boost::program_options::value<int>()->default_value(0), "Option 0: discrete Xdiv bracket, Option 1: Integral potential bracket");

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
  std::cout << FORMAT(25) << "[main] Degree" << K << std::endl;

  // Select the solution
  int solution = (vm.count("solution") ? vm["solution"].as<int>() : 1);
  double nonlinear_coeff = (vm.count("nonlinear-scaling") ? vm["nonlinear-scaling"].as<double>() : 1.);
  YangMills::TLAForcingTermType f;
  YangMills::TLAElectricFieldType E;
  YangMills::TLAElectricFieldType A;
  YangMills::TLAMagneticFieldType B;
  YangMills::TLAElectricFieldType dtE;

  std::cout << FORMAT(25) << "[main] Nonlinear coefficient " << nonlinear_coeff << std::endl;

  switch (std::abs(solution)) {
  case 1:
    std::cout << "[main] Trigonometric solution" << std::endl;
    f = sumLA<Eigen::Vector3d, YangMills::TForcingTermType>(trigonometric_f_linear, trigonometric_f_nonlinear, nonlinear_coeff);
    E = trigonometric_E; 
    A = trigonometric_A;
    B = sumLA<Eigen::Vector3d, YangMills::TMagneticFieldType>(trigonometric_B_linear, trigonometric_B_nonlinear, nonlinear_coeff);
    dtE = trigonometric_dtE;
    break;

  case 2:
    std::cout << "[main] Linear solution" << std::endl;
    f = sumLA<Eigen::Vector3d, YangMills::TForcingTermType>(linear_f_linear, linear_f_nonlinear, nonlinear_coeff);
    E = linear_E;
    A = linear_A;
    B = sumLA<Eigen::Vector3d, YangMills::TMagneticFieldType>(linear_B_linear, linear_B_nonlinear, nonlinear_coeff);
    dtE = linear_dtE;
    break;

  case 3:
    std::cout << "[main] Constant solution" << std::endl;
    f = sumLA<Eigen::Vector3d, YangMills::TForcingTermType>(const_f_linear, const_f_nonlinear, nonlinear_coeff);
    E = const_E;
    A = const_A;
    B = sumLA<Eigen::Vector3d, YangMills::TMagneticFieldType>(const_B_linear, const_B_nonlinear, nonlinear_coeff);
    dtE = const_dtE;
    break;

  default:
    std::cerr << "[main] ERROR: Unknown exact solution" << std::endl;
    exit(1);    
  }

  // Build the mesh
  MeshBuilder meshbuilder = MeshBuilder(mesh_file);
  std::unique_ptr<Mesh> mesh_ptr = meshbuilder.build_the_mesh();

  // Select the no. of iterations and final time
  double t_initial = 0.;
  double t_final = vm["final-time"].as<double>();
  int iter = vm["iterations"].as<int>();
  // Set the time step and iterations
  double dt;
  size_t it;
  if (iter > 0){
    it = iter;
    dt = (t_final - t_initial)/double(iter);
  } else {
    double measure3 = 0.;
    for (auto T : mesh_ptr->get_cells()){
      measure3 += T->measure();
    }
    if (iter < 0){
      it = size_t(ceil(std::abs(iter) * std::pow(double(std::cbrt(measure3))/mesh_ptr->h_max(), K+1)));
    }else{
      it = size_t(ceil(std::pow(double(std::cbrt(measure3))/mesh_ptr->h_max(), K+1)));
    }
    dt = (t_final - t_initial)/double(it);
  }
  // At least 10 iterations
  size_t min_it = 10;
  it = std::max(min_it, it);
  dt = std::min((t_final - t_initial)/double(min_it), dt);

  double stop_val = (vm.count("stopping-value") ? vm["stopping-value"].as<double>() : 1e-6);
  double theta = (vm.count("theta") ? vm["theta"].as<double>() : 1.);
  bool calculate_cond_num = (vm.count("condition-numbers") ? vm["condition-numbers"].as<bool>() : false);
  bool itersolver = (vm.count("iterative-solver") ? vm["iterative-solver"].as<bool>() : false);
  bool initial_cond = (vm.count("exact-initial-conditions") ? vm["exact-initial-conditions"].as<bool>() : false);
  int nonlinear_discretisation = (vm.count("bracket") ? vm["bracket"].as<int>() : 0);

  if (nonlinear_discretisation == 0){
    std::cout << FORMAT(25) << "[main] Using nonlinear Xdiv discretisation" << std::endl;
  }else if (nonlinear_discretisation == 1){
    std::cout << FORMAT(25) << "[main] Using nonlinear integral discretisation" << std::endl;
  }else{
    std::cerr << "[main] ERROR: Unknown bracket discretisation" << std::endl;
    exit(1);
  }

  std::cout << FORMAT(25) << "[main] Initial time" << t_initial << std::endl;
  std::cout << FORMAT(25) << "[main] Final time" << t_final << std::endl;
  std::cout << FORMAT(25) << "[main] Time step" << dt << std::endl;
  std::cout << FORMAT(25) << "[main] Stopping value" << stop_val << std::endl;
  std::cout << FORMAT(25) << "[main] Theta" << theta << std::endl;

  boost::timer::cpu_timer timer;
  // Create DDR core
  timer.start();
  bool use_threads = (vm.count("pthread") ? vm["pthread"].as<bool>() : true);
  std::cout << "[main] " << (use_threads ? "Parallel execution" : "Sequential execution") << std:: endl;
  DDRCore ddr_core(*mesh_ptr, K, use_threads);
  LieAlgebra su2;
  LADDRCore laddr_core(su2, *mesh_ptr, K);
  timer.stop();
  double t_wall_ddrcore, t_proc_ddrcore;
  std::tie(t_wall_ddrcore, t_proc_ddrcore) = store_times(timer, "[main] Time DDRCore (wall/proc) ");

  // Assemble the problem
  timer.start();

  YangMills ym(ddr_core, laddr_core, su2, nonlinear_discretisation, use_threads); 
  ym.assembleLinearSystem(dt);
  timer.stop();
  double t_wall_model, t_proc_model;
  std::tie(t_wall_model, t_proc_model) = store_times(timer, "[main] Time model (wall/proc) ");
  timer.start();
  // Initial conditions
  Eigen::VectorXd Elambdah = Eigen::VectorXd::Zero(ym.dimensionSpace());
  Eigen::VectorXd Ah = ym.laSXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ElectricFieldType>(A, t_initial), 2*K+2, 2*K+2, 2*K+2);
  Eigen::VectorXd Eh = ym.laSXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ElectricFieldType>(E, t_initial), 2*K+2, 2*K+2, 2*K+2);
  Elambdah.head(ym.laSXCurl().dimension()) = Eh;
  // Set the constrained initial conditions
  if (initial_cond){
    Elambdah.head(ym.laSXCurl().dimension()) = ym.computeInitialConditions(Eh, Ah, nonlinear_coeff, itersolver);
  }

  Eigen::VectorXd initial_constraint = ym.computeConstraint(Elambdah.head(ym.laSXCurl().dimension()), Ah, nonlinear_coeff);
  double initial_constraint_norm = ym.computeConstraintNorm(initial_constraint, itersolver);
  std::cout << FORMAT(25) << "[main] Initial constraint norm: " << initial_constraint_norm << std::endl;
  double max_res = 0;
  double max_const = 0;
  size_t nonlinear_it = 0;
  double t_wall_assemble, t_proc_assemble;
  double t_wall_solve, t_proc_solve;

  // Print solver
  if (itersolver) {
    std::cout << "[main] Solving the linear system using LeastSquaresConjugateGradient" << std::endl;
  } else { 
  #ifdef WITH_MKL
      std::cout << "[main] Solving the linear system using Pardiso" << std::endl;
      unsigned nb_threads_hint = std::thread::hardware_concurrency();
      mkl_set_dynamic(0);
      mkl_set_num_threads(nb_threads_hint);
  #elif WITH_UMFPACK
      std::cout << "[main] Solving the linear system using Umfpack" << std::endl;
  #else
      std::cout << "[main] Solving the linear system using direct solver" << std::endl;
  #endif
  }
  // File for condition numbers
  std::ofstream cn("condition-numbers.txt");
  cn << FORMAT(5) << "it" << std::setw(15) << "max_eig" << std::setw(15) << "min_eig" << std::setw(15) << "cond_num" << std::setw(15) << "residual" << std::endl;
  
  // Solve the problem
{
  size_t prog = 0;
  boost::timer::cpu_timer it_timer;
  for (size_t i = 0 ; i < it; i++){
    if (100*i/it >= prog){
      std::cout << "[main] " << 100*i/it << "%..." << std::endl;
      prog = 100*i/it + 10;
    }
    // Set the time
    double t = t_initial + dt*i;

    // Interpolate of forcing terms in LAXCurl: default at 0, but we put the ones corresponding to the exact solution if solution>0
    Eigen::VectorXd interp_f = Eigen::VectorXd::Zero(ym.laSXCurl().dimension());
    Eigen::VectorXd interp_Eipo = Eigen::VectorXd::Zero(ym.laSXCurl().dimension());
    Eigen::VectorXd interp_Ei = Eigen::VectorXd::Zero(ym.laSXCurl().dimension());
    Eigen::VectorXd interp_Aipomtheta = Eigen::VectorXd::Zero(ym.laSXCurl().dimension());
    if (solution>0){
      interp_f = ym.laSXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ForcingTermType>(f, t+dt), 2*K+2, 2*K+2, 2*K+2);
      interp_Eipo = ym.laSXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ForcingTermType>(E, t+dt), 2*K+2, 2*K+2, 2*K+2);
      interp_Ei = ym.laSXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ForcingTermType>(E, t), 2*K+2, 2*K+2, 2*K+2);
      interp_Aipomtheta = ym.laSXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ForcingTermType>(A, t+(1-theta)*dt), 2*K+2, 2*K+2, 2*K+2);
    }

    // System vector of the nonlinear problem
    Eigen::VectorXd Eh_i = Elambdah.head(ym.laSXCurl().dimension());
    ym.setSystemVector(interp_f, interp_Eipo-interp_Ei, interp_Aipomtheta, Eh_i, Ah, dt, theta, nonlinear_coeff);
    // Add boundary conditions for chosen solution, only if solution>0
    if (solution>0){
      ym.addBoundaryConditions(ym.contractParaLA<Eigen::Vector3d, YangMills::MagneticFieldType>(B, t+dt), dt);
    }
    
    // Solution vector of Newton iterations and residual
    Eigen::VectorXd dElambdak_condensed;
    ym.setNonlinearRes(Elambdah, Eh_i, Ah, dt, theta, nonlinear_coeff);

    while ((ym.systemVector().norm() < 1e-13 ? ym.nonlinearRes().norm() :  ym.nonlinearRes().norm()/ym.systemVector().norm()) > stop_val){
      nonlinear_it++;

      it_timer.start();
      ym.assembleSystemNewton(Eh_i, Ah, Elambdah, dt, theta, nonlinear_coeff);
      it_timer.stop();
      double it_wall_assemble, it_proc_assemble;
      std::tie(it_wall_assemble, it_proc_assemble) = store_times(it_timer);
      // std::tie(it_wall_assemble, it_proc_assemble) = store_times(it_timer, "[main] Time assemble (wall/proc) ");
      t_wall_assemble += it_wall_assemble;
      t_proc_assemble += it_proc_assemble;

      it_timer.start();
      if (itersolver){
        Eigen::LeastSquaresConjugateGradient<YangMills::SystemMatrixType> solver;
        // solver.preconditioner().setFillfactor(2);
        solver.compute(ym.systemMatrix());
        if (solver.info() != Eigen::Success) {
          std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
          exit(1);
        }
        dElambdak_condensed = solver.solve(ym.systemVectorNewton());
        if (solver.info() != Eigen::Success) {
          std::cerr << "[main] ERROR: Could not solve direct system" << std::endl;
          exit(1);
        }
      } else { 
        #ifdef WITH_MKL   
            Eigen::PardisoLU<YangMills::SystemMatrixType> solver;
        #elif WITH_UMFPACK   
            Eigen::UmfPackLU<YangMills::SystemMatrixType> solver;
        #else 
            Eigen::SparseLU<YangMills::SystemMatrixType> solver;
        #endif
            solver.compute(ym.systemMatrix());
            if (solver.info() != Eigen::Success) {
              std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
            }
            dElambdak_condensed = solver.solve(ym.systemVectorNewton());
            if (solver.info() != Eigen::Success) {
              std::cerr << "[main] ERROR: Could not solve linear system" << std::endl;
            }
        }
      // Re-create statically condensed unknowns
      Eigen::VectorXd dElambdak = Eigen::VectorXd::Zero(ym.dimensionSpace());
      Eigen::VectorXd sc_unknowns = ym.scVector() + ym.scMatrix() * dElambdak_condensed;
      dElambdak.head(ym.laSXCurl().dimension() - ym.nbSCDOFs_E()) = dElambdak_condensed.head(ym.laSXCurl().dimension() - ym.nbSCDOFs_E());
      dElambdak.segment(ym.laSXCurl().dimension() - ym.nbSCDOFs_E(), ym.nbSCDOFs_E()) = sc_unknowns.head(ym.nbSCDOFs_E());
      dElambdak.segment(ym.laSXCurl().dimension(), ym.laSXGrad().dimension() - ym.nbSCDOFs_lambda()) 
          = dElambdak_condensed.segment(ym.laSXCurl().dimension() - ym.nbSCDOFs_E(), ym.laSXGrad().dimension() - ym.nbSCDOFs_lambda());
      dElambdak.tail(ym.nbSCDOFs_lambda()) = sc_unknowns.tail(ym.nbSCDOFs_lambda());
      Elambdah = Elambdah + dElambdak;

      it_timer.stop();
      double it_wall_solve, it_proc_solve;
      std::tie(it_wall_solve, it_proc_solve) = store_times(it_timer);
      // std::tie(it_wall_solve, it_proc_solve) = store_times(it_timer, "[main] Time solve (wall/proc) ");
      t_wall_solve += it_wall_solve;
      t_proc_solve += it_proc_solve;

      // Update residual and condition numbers
      double res = ym.computeResidual(dElambdak_condensed);
      // std::cout << FORMAT(25) << "[main] Residual: " << res << std::endl;
      if (calculate_cond_num){
        #ifdef WITH_SPECTRA
        auto [max_eig, min_eig] = ym.computeConditionNum();
        double cond_num = std::sqrt(max_eig/double(min_eig));
        // std::cout << FORMAT(25) << "[main] Condition number: " << cond_num << std::endl;
        cn << std::setw(5) << i << std::setw(15) << max_eig << std::setw(15) << min_eig << std::setw(15) << cond_num << std::setw(15) << res << std::endl;
        #endif
      }
      max_res = std::max(max_res, res);
      ym.setNonlinearRes(Elambdah, Eh_i, Ah, dt, theta, nonlinear_coeff);
    };

    // Update A and the maximum difference from initial constraint
    Eigen::VectorXd Eh_ipo = Elambdah.head(ym.laSXCurl().dimension());
    Ah -= dt*(1-theta)*Eh_i + dt*theta*Eh_ipo;
    max_const = std::max(max_const, ym.computeConstraintNorm(ym.computeConstraint(Eh_ipo, Ah, nonlinear_coeff)-initial_constraint, itersolver));
  } // for i
}
  cn << std::flush;
  cn.close();
  timer.stop();
  double t_wall_run, t_proc_run;
  std::tie(t_wall_run, t_proc_run) = store_times(timer, "[main] Time run (wall/proc) ");
  t_wall_assemble = t_wall_assemble/double(nonlinear_it);
  t_proc_assemble = t_proc_assemble/double(nonlinear_it);
  t_wall_solve = t_wall_solve/double(nonlinear_it);
  t_proc_solve = t_proc_solve/double(nonlinear_it);

  // Export matrix if requested  
  if (vm.count("export-matrix")) {
    std::cout << "[main] Exporting matrix to Matrix Market format" << std::endl;
    saveMarket(ym.systemMatrix(), "A_maxwell.mtx");
    // saveMarket(ym.systemVector(), "B_maxwell.mtx");
  }

  // Interpolate of exact solution and error vector
  Eigen::VectorXd ElambdaAI = Eigen::VectorXd::Zero(ym.dimensionSpace() + ym.laSXCurl().dimension());
  Eigen::VectorXd ElambdaAh = Eigen::VectorXd::Zero(ym.dimensionSpace() + ym.laSXCurl().dimension());
  ElambdaAI.head(ym.laSXCurl().dimension()) = ym.laSXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ElectricFieldType>(E, t_final), 2*K+2, 2*K+2, 2*K+2);
  ElambdaAI.tail(ym.laSXCurl().dimension()) = ym.laSXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ElectricFieldType>(A, t_final), 2*K+2, 2*K+2, 2*K+2);
  ElambdaAh.head(ym.dimensionSpace()) = Elambdah;
  ElambdaAh.tail(ym.laSXCurl().dimension()) = Ah;
  Eigen::VectorXd erh = ElambdaAh - ElambdaAI;

  std::cout << "[main] Compute errors" << std::endl;
  timer.start();
  // Errors in discrete norms
  std::vector<YangMillsNorms> list_norms = ym.computeYangMillsNorms(std::vector<Eigen::VectorXd> {erh, ElambdaAI});
  YangMillsNorms errors = list_norms[0];
  YangMillsNorms norms = list_norms[1];
  std::cout << "[main] Electric field error = " << (norms.E < 1e-12 ? errors.E : errors.E/norms.E) << std::endl;
  std::cout << "[main] Potential error = " << (norms.A < 1e-12 ? errors.A : errors.A/norms.A) << std::endl;
  std::cout << "[main] Lambda = " << errors.lambda << std::endl;
  std::cout << "[main] Maximum residual = " << max_res << std::endl;
  std::cout << "[main] Maximum constraint difference = " << max_const << std::endl;
  timer.stop();
  double t_wall_norms, t_proc_norms;
  std::tie(t_wall_norms, t_proc_norms) = store_times(timer, "[main] Time norms (wall/proc) ");
  
  std::cout << "[main] Mesh diameter " << mesh_ptr->h_max() << std::endl;

  // Write results to file
  std::ofstream out("results.txt");
  out << "Solution: " << solution << std::endl;
  out << "NonlinearDiscretisation" << nonlinear_discretisation << std::endl;
  out << "StabilisationParameter: " << ym.stabilizationParameter() << std::endl;
  out << "Mesh: " << mesh_file << std::endl;
  out << "Degree: " << K << std::endl;
  out << "MeshSize: " << mesh_ptr->h_max() << std::endl;
  out << "NbCells: " << mesh_ptr->n_cells() << std::endl;
  out << "NbFaces: " << mesh_ptr->n_faces() << std::endl;
  out << "NbEdges: " << mesh_ptr->n_edges() << std::endl;
  out << "DimSXGrad: " << ym.xSGrad().dimension() << std::endl;
  out << "DimSXCurl: " << ym.xSCurl().dimension() << std::endl;
  out << "DimSXDiv: " << ym.xSDiv().dimension() << std::endl;
  out << "DimLieAlg: " << ym.lieAlg().dimension() << std::endl;
  out << "DimLASXGrad: " << ym.laSXGrad().dimension() << std::endl;
  out << "DimLASXCurl: " << ym.laSXCurl().dimension() << std::endl;
  out << "DimLASXDiv: " << ym.laSXDiv().dimension() << std::endl;
  out << "InitialTime: " << t_initial << std::endl;
  out << "FinalTime: " << t_final << std::endl;
  out << "Timestep: " << dt << std::endl;
  out << "StoppingValue: " << stop_val << std::endl;
  out << "Theta: " << theta << std::endl;
  out << "MaxResidual: " << max_res << std::endl;
  out << "ConstrainedInitialConditions: " << initial_cond << std::endl;
  out << "InitialConstraintNorm: " << initial_constraint_norm << std::endl;
  out << "MaxConstraintDiff: " << max_const << std::endl;
  out << "OuterIterations: " << it << std::endl;
  out << "AvgNonlinearIterations: " << nonlinear_it/double(it) << std::endl;
  // Discrete errors and norms
  out << "E_L2Elec: " << (norms.E < 1e-12 ? errors.E : errors.E/norms.E) << std::endl;
  out << "E_L2Pot: " << (norms.A < 1e-12 ? errors.A : errors.A/norms.A) << std::endl;
  out << "absE_L2Elec: " << errors.E << std::endl;
  out << "absE_L2Pot: " << errors.A << std::endl;
  out << "absE_L2Lambda: " << errors.lambda << std::endl;
  out << "N_L2Elec: " << norms.E << std::endl;
  out << "N_L2Pot: " << norms.A << std::endl;
  out << "TwallDDRCore: " << t_wall_ddrcore << std::endl;
  out << "TprocDDRCore: " << t_proc_ddrcore << std::endl;
  out << "TwallModel: " << t_wall_model << std::endl;
  out << "TprocModel: " << t_proc_model << std::endl;
  out << "TwallAssemble: " << t_wall_assemble << std::endl;
  out << "TprocAssemble: " << t_proc_assemble << std::endl;
  out << "TwallSolve: " << t_wall_solve << std::endl;
  out << "TprocSolve: " << t_proc_solve << std::endl;
  out << "TwallRun: " << t_wall_run << std::endl;  
  out << "TprocRun: " << t_proc_run << std::endl; 
  if (itersolver){
    out << "Solver: " << "LeastSquaresConjugateGradient" << std::endl;
  } else { 
  #ifdef WITH_MKL
    out << "Solver: " << "Pardiso" << std::endl;
  #elif WITH_UMFPACK
    out << "Solver: " << "Umfpack" << std::endl;
  #else
    out << "Solver: " << "direct solver" << std::endl;
  #endif
  }
  out << std::flush;
  out.close();

  std::cout << "[main] Done" << std::endl;
  return 0;
}

//------------------------------------------------------------------------------
// YangMills
//------------------------------------------------------------------------------

YangMills::YangMills(
               const DDRCore & ddrcore,
               const LADDRCore & laddrcore,
               const LieAlgebra & lie_algebra,
               size_t nonlinear_discretisation,
               bool use_threads,
               std::ostream & output
               )
  : m_ddrcore(ddrcore),
    m_laddrcore(laddrcore),
    m_use_threads(use_threads),
    m_output(output),
    m_ser_pro(ddrcore, use_threads, output),
    m_xdiv(ddrcore, use_threads),
    m_sxgrad(ddrcore, m_ser_pro, use_threads),
    m_sxcurl(ddrcore, m_ser_pro, use_threads),
    m_sxdiv(ddrcore, use_threads),
    m_liealg(lie_algebra),
    m_lasxgrad(lie_algebra, m_sxgrad, use_threads),
    m_lasxcurl(lie_algebra, m_sxcurl, use_threads),
    m_lasxdiv(lie_algebra, m_xdiv, m_sxdiv, use_threads),
    m_nloc_sc_E(m_liealg.dimension()*m_ser_pro.nDOFs_cells_SXCurl()),
    m_nloc_sc_lambda(m_liealg.dimension()*m_ser_pro.nDOFs_cells_SXGrad()),
    m_A(sizeSystem(), sizeSystem()),
    m_b(Eigen::VectorXd::Zero(sizeSystem())),
    m_b_i(Eigen::VectorXd::Zero(dimensionSpace())),
    m_b_k(Eigen::VectorXd::Zero(dimensionSpace())),
    m_sc_A(nbSCDOFs(), sizeSystem()),
    m_sc_b(Eigen::VectorXd::Zero(nbSCDOFs())),
    m_laxgrad_L2(m_lasxgrad.dimension(), m_lasxgrad.dimension()),
    m_laxcurl_L2(m_lasxcurl.dimension(), m_lasxcurl.dimension()),
    m_laxdiv_L2(m_lasxdiv.dimension(), m_lasxdiv.dimension()),
    m_nl_par(nonlinear_discretisation),
    m_stab_par(0.1)
{
  m_output << "[YangMills] Initializing" << std::endl;

  // To avoid performing static condensation, initialise m_nloc_sc_E and m_nloc_sc_lambda at 0
}

//------------------------------------------------------------------------------

void YangMills::assembleLinearSystem(
                                  double dt
                                  )
{
  // Assemble all local contributions
  auto assemble = [this, dt](
                            size_t start,
                            size_t end,
                            std::vector<std::list<Eigen::Triplet<double>>> * triplets,
                            std::vector<Eigen::VectorXd> * vecs
                            )->void
                      {
                        for (size_t iT = start; iT < end; iT++) {
                          this->_assemble_local_contribution(
                                                             iT,
                                                             this->_compute_local_contribution(iT, dt),
                                                             *triplets,
                                                             *vecs
                                                             );
                        } // for iT
                      };
                    
  // Assemble the product matrices
  if (m_use_threads) {
    m_output << "[YangMills] Parallel assembly" << std::endl;
  }else{
    m_output << "[YangMills] Sequential assembly" << std::endl;
  }
  
  std::vector<std::pair<size_t, size_t>> size_systems{std::make_pair(m_lasxcurl.dimension(), m_lasxcurl.dimension()), std::make_pair(m_lasxdiv.dimension(), m_lasxdiv.dimension()), std::make_pair(m_lasxgrad.dimension(), m_lasxgrad.dimension())}; 
  auto [systems, vectors] = parallel_assembly_system(m_ddrcore.mesh().n_cells(), size_systems, {}, assemble, m_use_threads);

  m_laxcurl_L2 = systems[0];
  m_laxdiv_L2 = systems[1];
  m_laxgrad_L2 = systems[2];

}

//------------------------------------------------------------------------------

void YangMills::addBoundaryConditions(
                                    const LAMagneticFieldType & B,
                                    double dt
                                    )
{
  // Assemble boundary conditions
  for (auto iF : m_ddrcore.mesh().get_b_faces()) {
    const Face & F = *iF;
    
    // Unit normal vector to F pointing out of the domain
    const Cell & TF = *F.cell(0);
    Eigen::Vector3d nF = TF.face_normal(TF.index_face(&F));

    // Degree of quadratures to compute boundary conditions
    const size_t dqrbc = 2 * m_ddrcore.degree() + 3;

    // Boundary condition on the tangential component of curl B
    {
      Eigen::MatrixXd int_Bi_cross_nF(m_liealg.dimension(), m_sxcurl.dimensionFace(F));
      QuadratureRule quad_dqrbc_F = generate_quadrature_rule(F, dqrbc);
      
      for (size_t i = 0; i < m_liealg.dimension(); i++){
        MagneticFieldType Bi = B[i];
        FType<Eigen::Vector3d> Bi_cross_nF = [&Bi, &nF](const Eigen::Vector3d & x) {
                                                        return Bi(x).cross(nF);
                                                        };
        int_Bi_cross_nF.row(i) = integrate(Bi_cross_nF, evaluate_quad<Function>::compute(*m_ddrcore.faceBases(F.global_index()).Polyk2, quad_dqrbc_F), quad_dqrbc_F).transpose() 
        * m_sxcurl.facePotential(F.global_index());
      }// for i
      Eigen::MatrixXd bF_M = m_liealg.massMatrix() * int_Bi_cross_nF;
      Eigen::Map<Eigen::VectorXd> bF(bF_M.data(), bF_M.size());
      auto I_F = m_lasxcurl.globalDOFIndices(F);
      for (size_t i = 0; i < I_F.size(); i++) {
        m_b_i(I_F[i]) -= dt*bF(i);
      } // for i
    }
  } // for iF
}

//------------------------------------------------------------------------------

std::pair<double, double> YangMills::computeConditionNum() const
{
  #ifdef WITH_SPECTRA
  Eigen::SparseMatrix<double> M = systemMatrix().transpose()*systemMatrix();
  // Calculate the largest eigenvalue
  Spectra::SparseSymMatProd<double> op(M);
  Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs_max(op, 1, 20);
  eigs_max.init();
	eigs_max.compute(Spectra::SortRule::LargestMagn);
  // Calculate the smallest eigenvalue
  Spectra::SparseSymShiftSolve<double> opp(M);
  Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> eigs_min(opp, 1, 20, 0.0);
  eigs_min.init();
  eigs_min.compute(Spectra::SortRule::LargestMagn);
  
  if(eigs_max.info() == Spectra::CompInfo::Successful && eigs_min.info() == Spectra::CompInfo::Successful){
    // std::cout << FORMAT(25) << "[main] Eigenvalues: " << eigs_max.eigenvalues()[0] << " " << eigs_min.eigenvalues()[0] << std::endl;
    return std::make_pair(eigs_max.eigenvalues()[0], eigs_min.eigenvalues()[0]);
  }
  #endif
  return std::make_pair(-1, -1);
}

//------------------------------------------------------------------------------

double YangMills::stoppingCrit(
                            const Eigen::VectorXd & v,
                            const Eigen::VectorXd & u
                            )
{
  double err = std::sqrt((v - u).transpose() * m_laxcurl_L2 * (v - u));
  double norm = std::sqrt((u.transpose() * m_laxcurl_L2 * u));
  return (norm > 1e-12) ? err/double(norm) : err;
}

//------------------------------------------------------------------------------

void YangMills::setNonlinearRes(
                                const Eigen::VectorXd & Elambda_k,
                                const Eigen::VectorXd & E_i,
                                const Eigen::VectorXd & A_i,
                                double dt,
                                double theta,
                                double nonlinear_coeff
                                )
{
  // Assemble local vectors
  std::function<void(size_t, const SystemVectors<Eigen::MatrixXd> &, std::vector<std::list<Eigen::Triplet<double> > > &, std::vector<Eigen::VectorXd> &)> assemble_local_vec
  = [this](size_t iT, const SystemVectors<Eigen::MatrixXd> & lsT, std::vector<std::list<Eigen::Triplet<double>>> & triplets_sys, std::vector<Eigen::VectorXd> & vecs)-> void
  {
    auto & [system, vectors] = lsT;
    Cell & T = *m_ddrcore.mesh().cell(iT);
    std::vector<size_t> IT_laxcurl = m_lasxcurl.globalDOFIndices(T);
    std::vector<size_t> IT_laxgrad = m_lasxgrad.globalDOFIndices(T);
    
    size_t dim_laxcurl = m_lasxcurl.dimension();
    size_t dim_laxcurl_T = m_lasxcurl.dimensionCell(iT);

    for (size_t i = 0; i < IT_laxcurl.size(); i++){
      vecs[0][IT_laxcurl[i]] += vectors[0](i);
    }
    for (size_t i = 0; i < IT_laxgrad.size(); i++){
      vecs[0][dim_laxcurl + IT_laxgrad[i]] += vectors[0](dim_laxcurl_T + i);
    }
  };
  // Assemble
  auto assemble = [&](size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs) -> void
  {
    for (size_t iT = start; iT < end; iT++){
      assemble_local_vec(iT, _compute_local_newton(iT, 1, E_i, A_i, Elambda_k, dt, theta, nonlinear_coeff), *triplets, *vecs);
    }
  };

  Eigen::VectorXd F_Ek = parallel_assembly_system(m_ddrcore.mesh().n_cells(), {}, {this->dimensionSpace()}, assemble, m_use_threads).vectors[0];
  
  m_b_k = m_b_i-F_Ek;
  // std::cout << FORMAT(25) << "[main] Stop crit: " << nonlinearRes().norm()/systemVector().norm() << std::endl;
}

//------------------------------------------------------------------------------

void YangMills::setSystemVector(
                                const Eigen::VectorXd & interp_f,
                                const Eigen::VectorXd & interp_dE,
                                const Eigen::VectorXd & interp_A,
                                const Eigen::VectorXd & E_old,
                                const Eigen::VectorXd & A_old,
                                double dt,
                                double theta,
                                double nonlinear_coeff
                                )
  { 
  // Assemble vector
  auto assemble = [&](size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs) -> void
  {
    for (size_t iT = start; iT < end; iT++){
      _assemble_local_vec(iT, _compute_local_vec(iT, interp_f, interp_dE, interp_A, E_old, A_old, dt, theta, nonlinear_coeff), *triplets, *vecs);
    }
  };
  m_b_i = parallel_assembly_system(m_ddrcore.mesh().n_cells(), {}, {this->dimensionSpace()}, assemble, m_use_threads).vectors[0];
}

//------------------------------------------------------------------------------

void YangMills::assembleSystemNewton(
                                    const Eigen::VectorXd & E_i,    
                                    const Eigen::VectorXd & A_i,    
                                    const Eigen::VectorXd & Elambda_k,
                                    double dt,                       
                                    double theta,
                                    double nonlinear_coeff
                                    )
{ 
  auto assemble = [&](size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs) -> void
  {
    for (size_t iT = start; iT < end; iT++){
      _assemble_local_newton(iT, _compute_local_newton(iT, 2, E_i, A_i, Elambda_k, dt, theta, nonlinear_coeff), *triplets, *vecs);
    }
  };
  
  // Assemble the matrix and set the rhs 
  auto [systems, vectors] = parallel_assembly_system(m_ddrcore.mesh().n_cells(), {std::make_pair(sizeSystem(), sizeSystem()), std::make_pair(this->nbSCDOFs(), this->sizeSystem())}, {sizeSystem(), this->nbSCDOFs()}, assemble, m_use_threads);

  m_A = systems[0];
  m_b = vectors[0];
  m_b.head(m_lasxcurl.dimension()-nbSCDOFs_E()) += m_b_k.head(m_lasxcurl.dimension()-nbSCDOFs_E());
  m_b.segment(m_lasxcurl.dimension()-nbSCDOFs_E(), m_lasxgrad.dimension()-nbSCDOFs_lambda()) += m_b_k.segment(m_lasxcurl.dimension(), m_lasxgrad.dimension()-nbSCDOFs_lambda());

  m_sc_A = systems[1];
  m_sc_b = vectors[1];
}

//------------------------------------------------------------------------------

Eigen::VectorXd YangMills::_compute_local_vec(size_t iT, 
                                              const Eigen::VectorXd & interp_f,
                                              const Eigen::VectorXd & interp_dE,
                                              const Eigen::VectorXd & interp_A,
                                              const Eigen::VectorXd & E_i,
                                              const Eigen::VectorXd & A_i,
                                              double dt,
                                              double theta,
                                              double nonlinear_coeff)
{

  Cell & T = *m_ddrcore.mesh().cell(iT);
  size_t dim_P2k3 = m_laddrcore.P2k3(iT).dimension();
  
  // Mass matrix for (P^k(T))^3  
  MonomialCellIntegralsType int_mono_3k = IntegrateCellMonomials(T, 3*m_ddrcore.degree());
  Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_3k);

  // Vectors
  Eigen::VectorXd X = A_i - dt * (1-theta) * E_i;
  Eigen::VectorXd Y = A_i - 0.5 * dt * (1-theta) * E_i;
  Eigen::VectorXd Z = A_i - dt * (1-theta) * (1-theta) * E_i;

  // Local vectors
  Eigen::VectorXd A_old_T = m_lasxcurl.restrict(T, A_i);
  Eigen::VectorXd E_old_T = m_lasxcurl.restrict(T, E_i);
  Eigen::VectorXd interp_f_T = m_lasxcurl.restrict(T, interp_f);
  Eigen::VectorXd interp_dE_T = m_lasxcurl.restrict(T, interp_dE);
  Eigen::VectorXd X_T = m_lasxcurl.restrict(T, X);

  // Bracket matrices: *[u,.]^div
  Scalar3Tensor crossij_Pot_T = _compute_crossij_Pot_T(iT);
  std::vector<Eigen::VectorXd> ebkt_vec_list{A_i, X, Y};
  std::vector<Eigen::MatrixXd> ebkt_list = epsBkt_v(iT, crossij_Pot_T, ebkt_vec_list);
  Eigen::MatrixXd ebkt_A = nonlinear_coeff * ebkt_list[0];
  Eigen::MatrixXd ebkt_X = nonlinear_coeff * ebkt_list[1];
  Eigen::MatrixXd ebkt_Y = nonlinear_coeff * ebkt_list[2];

  // Bracket vectors: *[u,v]^div
  Eigen::VectorXd ebkt_A_A = ebkt_A * A_old_T;
  Eigen::VectorXd ebkt_X_X = ebkt_X * X_T;

  // L2 product matrices
  Eigen::MatrixXd laxcurl_L2 = m_lasxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd laxcurl_L2Grad_left = m_lasxcurl.computeL2ProductGradient(iT, m_sxgrad, "left", m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd laxdiv_L2Curl_both = m_lasxdiv.computeL2ProductCurl(iT, m_sxcurl, "both", m_stab_par, mass_Pk3_T);
  // These depend on the nonlinear discretisation we use: (Xdiv discretisation = LaXdiv product (and curl), integral discretisation = mass matrix of La P2k3 (and cell curl))
  Eigen::MatrixXd L2_nlProduct_T;
  Eigen::MatrixXd L2Curl_nlProduct_T;
  if (m_nl_par == 0){
    L2_nlProduct_T = m_lasxdiv.computeL2Product(iT, m_stab_par, mass_Pk3_T);
    L2Curl_nlProduct_T = m_lasxdiv.computeL2ProductCurl(iT, m_sxcurl, "left", m_stab_par, mass_Pk3_T);
  }else if (m_nl_par == 1){
    Eigen::MatrixXd mass_Pk3_P2k3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, m_laddrcore.P2k3(iT), int_mono_3k);
    Eigen::MatrixXd mass_P2k3_T = Eigen::MatrixXd::Identity(dim_P2k3, dim_P2k3); // Since we assume orthonormal basis
    Eigen::MatrixXd cell_curl = m_lasxcurl.cellCurl(iT);
    L2_nlProduct_T = Eigen::KroneckerProduct(mass_P2k3_T, m_liealg.massMatrix());
    L2Curl_nlProduct_T = cell_curl.transpose() * Eigen::KroneckerProduct(mass_Pk3_P2k3_T, m_liealg.massMatrix());
  }

  // L2 product matrices with brackets
  Eigen::MatrixXd L2_ebkt_Y = ebkt_Y.transpose() * L2_nlProduct_T;
  Eigen::MatrixXd L2Curl_left_ebkt_Y = L2Curl_nlProduct_T * ebkt_Y;

  // L2 integral bracket product matrices: int <1,[2,3]>
  Scalar3Tensor int_PciPcjPgk = _integral_Pot_ijk(iT);
  std::vector<Eigen::MatrixXd> L2v_Bkt_2 = L2v_Bkt(iT, int_PciPcjPgk, {m_lasxcurl.restrictCell(iT, interp_A), 
                                                                        m_lasxcurl.restrictCell(iT, Z)}, 2);
  Eigen::MatrixXd L2bkt_interp_A = nonlinear_coeff * L2v_Bkt_2[0];
  Eigen::MatrixXd L2bkt_Z = nonlinear_coeff * L2v_Bkt_2[1];

  //------------------------------------------------------------------------------
  // System Vector 
  //------------------------------------------------------------------------------

  size_t dim_laxgrad_T = m_lasxgrad.dimensionCell(iT);
  size_t dim_laxcurl_T = m_lasxcurl.dimensionCell(iT);

  Eigen::VectorXd b(dim_laxcurl_T+dim_laxgrad_T);

  b.head(dim_laxcurl_T)
    = laxcurl_L2 * E_old_T
    + dt * (1-theta) * 
    (
      laxdiv_L2Curl_both * A_old_T
      + L2Curl_left_ebkt_Y.transpose() * A_old_T
      + 0.5 * L2Curl_nlProduct_T * ebkt_A_A
      + 0.5 * L2_ebkt_Y * ebkt_A_A
    )
    + dt * theta *
    (
      laxdiv_L2Curl_both * X_T
      + L2Curl_left_ebkt_Y.transpose() * X_T
      + 0.5 * L2Curl_nlProduct_T * ebkt_X_X
      + 0.5 * L2_ebkt_Y * ebkt_X_X
    )
    + dt * laxcurl_L2 * interp_f_T;
  
  b.segment(dim_laxcurl_T, dim_laxgrad_T)
    =  laxcurl_L2Grad_left * E_old_T
      + L2bkt_Z.transpose() * E_old_T
      + laxcurl_L2Grad_left * interp_dE_T
      + L2bkt_interp_A.transpose() * interp_dE_T;

  return b;
}

//------------------------------------------------------------------------------

void YangMills::_assemble_local_vec(size_t iT, 
                                   const Eigen::VectorXd & v, 
                                   std::vector<std::list<Eigen::Triplet<double>>> & triplets_sys, 
                                   std::vector<Eigen::VectorXd> & vecs)
{
  // Assemble local system vector
  Cell & T = *m_ddrcore.mesh().cell(iT);
  std::vector<size_t> IT_laxgrad = m_lasxgrad.globalDOFIndices(T);
  std::vector<size_t> IT_laxcurl = m_lasxcurl.globalDOFIndices(T);
  
  size_t dim_laxcurl = m_lasxcurl.dimension();
  size_t dim_laxcurl_T = m_lasxcurl.dimensionCell(iT);

  for (size_t i = 0; i < IT_laxcurl.size(); i++){
    vecs[0][IT_laxcurl[i]] += v(i);
  }
  for (size_t i = 0; i < IT_laxgrad.size(); i++){
    vecs[0][dim_laxcurl + IT_laxgrad[i]] += v(dim_laxcurl_T + i);
  }
}
//------------------------------------------------------------------------------

SystemVectors<Eigen::MatrixXd> 
YangMills::_compute_local_newton(
                                size_t iT,
                                size_t Res1_DF2,
                                const Eigen::VectorXd & E_i,    
                                const Eigen::VectorXd & A_i,    
                                const Eigen::VectorXd & Elambda_k,
                                double dt,
                                double theta,
                                double nonlinear_coeff
                                )
{
  Cell & T = *m_ddrcore.mesh().cell(iT);
  size_t dim_P2k3 = m_laddrcore.P2k3(iT).dimension();

  // Mass matrix for (P^k(T))^3  
  MonomialCellIntegralsType int_mono_3k = IntegrateCellMonomials(T, 3*m_ddrcore.degree());
  Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_3k);

  // k is the inner loop index, i is the outer loop index (constant in k)
  // Vectors
  Eigen::VectorXd E_k = Elambda_k.head(m_lasxcurl.dimension());
  Eigen::VectorXd lambda_k = Elambda_k.segment(m_lasxcurl.dimension(), m_lasxgrad.dimension());
  Eigen::VectorXd X = A_i - dt * (1-theta) * E_i;
  Eigen::VectorXd Y = A_i - 0.5 * dt * (1-theta) * E_i;
  Eigen::VectorXd Z = A_i - dt * (1-theta) * (1-theta) * E_i;

  // Local vectors
  Eigen::VectorXd E_k_T = m_lasxcurl.restrict(T, E_k);
  Eigen::VectorXd lambda_k_T = m_lasxgrad.restrict(T, lambda_k);
  Eigen::VectorXd A_i_T = m_lasxcurl.restrict(T, A_i);
  Eigen::VectorXd X_T = m_lasxcurl.restrict(T, X);
  Eigen::VectorXd Y_T = m_lasxcurl.restrict(T, Y);

  // Local bracket matrices: *[u,.]^div
  Scalar3Tensor crossij_Pot_T = _compute_crossij_Pot_T(iT);
  std::vector<Eigen::VectorXd> ebkt_vec_list{A_i, E_k, X, Y};
  std::vector<Eigen::MatrixXd> ebkt_list = epsBkt_v(iT, crossij_Pot_T, ebkt_vec_list);

  Eigen::MatrixXd ebkt_A = nonlinear_coeff * ebkt_list[0];
  Eigen::MatrixXd ebkt_E = nonlinear_coeff * ebkt_list[1];
  Eigen::MatrixXd ebkt_X = nonlinear_coeff * ebkt_list[2];
  Eigen::MatrixXd ebkt_Y = nonlinear_coeff * ebkt_list[3];

  // Local bracket vectors: *[u,v]^div
  Eigen::VectorXd ebkt_A_A = ebkt_A * A_i_T;
  Eigen::VectorXd ebkt_X_X = ebkt_X * X_T;

  // L2 product matrices
  Eigen::MatrixXd laxcurl_L2 = m_lasxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd laxcurl_L2Grad_right = m_lasxcurl.computeL2ProductGradient(iT, m_sxgrad, "right", m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd laxdiv_L2Curl_both = m_lasxdiv.computeL2ProductCurl(iT, m_sxcurl, "both", m_stab_par, mass_Pk3_T);
  // These depend on the nonlinear discretisation we use: (Xdiv discretisation = LaXdiv product (and curl), integral discretisation = mass matrix of La P2k3 (and cell curl))
  Eigen::MatrixXd L2_nlProduct_T;
  Eigen::MatrixXd L2Curl_nlProduct_T;
  if (m_nl_par == 0){
    L2_nlProduct_T = m_lasxdiv.computeL2Product(iT, m_stab_par, mass_Pk3_T);
    L2Curl_nlProduct_T = m_lasxdiv.computeL2ProductCurl(iT, m_sxcurl, "left", m_stab_par, mass_Pk3_T);
  }else if (m_nl_par == 1){
    Eigen::MatrixXd mass_Pk3_P2k3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, m_laddrcore.P2k3(iT), int_mono_3k);
    Eigen::MatrixXd mass_P2k3_T = Eigen::MatrixXd::Identity(dim_P2k3, dim_P2k3); // Since we assume orthonormal basis
    Eigen::MatrixXd cell_curl = m_lasxcurl.cellCurl(iT);
    L2_nlProduct_T = Eigen::KroneckerProduct(mass_P2k3_T, m_liealg.massMatrix());
    L2Curl_nlProduct_T = cell_curl.transpose() * Eigen::KroneckerProduct(mass_Pk3_P2k3_T, m_liealg.massMatrix());
  }

  // L2 product matrices with brackets
  std::vector<Eigen::MatrixXd> L2v_epsBkt_list = L2v_epsBkt(iT, crossij_Pot_T, {A_i_T.transpose() * L2Curl_nlProduct_T, 
                                                                                X_T.transpose() * L2Curl_nlProduct_T,
                                                                                ebkt_A_A.transpose() * L2_nlProduct_T,
                                                                                ebkt_X_X.transpose() * L2_nlProduct_T});
  Eigen::MatrixXd L2CurlA_ebkt = nonlinear_coeff * L2v_epsBkt_list[0];
  Eigen::MatrixXd L2CurlX_ebkt = nonlinear_coeff * L2v_epsBkt_list[1];
  Eigen::MatrixXd L2ebktAA_ebkt = L2v_epsBkt_list[2];
  Eigen::MatrixXd L2ebktXX_ebkt = L2v_epsBkt_list[3];
  Eigen::MatrixXd L2Curl_left_ebkt_E = L2Curl_nlProduct_T * ebkt_E;
  Eigen::MatrixXd L2Curl_left_ebkt_X = L2Curl_nlProduct_T * ebkt_X;
  Eigen::MatrixXd L2Curl_left_ebkt_Y = L2Curl_nlProduct_T * ebkt_Y;
  Eigen::MatrixXd L2_ebkt_E = ebkt_E.transpose() * L2_nlProduct_T;

  // L2 integral bracket product matrices: int <1,[2,3]>
  Scalar3Tensor int_PciPcjPgk = _integral_Pot_ijk(iT);
  std::vector<Eigen::MatrixXd> L2v_Bkt_1 = L2v_Bkt(iT, int_PciPcjPgk, {m_lasxcurl.restrictCell(iT, E_k), 
                                                                       m_lasxcurl.restrictCell(iT, E_i)}, 1);
  Eigen::MatrixXd L2E_bkt = nonlinear_coeff * L2v_Bkt_1[0];
  Eigen::MatrixXd L2Ei_bkt = nonlinear_coeff * L2v_Bkt_1[1];
  std::vector<Eigen::MatrixXd> L2v_Bkt_2 = L2v_Bkt(iT, int_PciPcjPgk, {m_lasxcurl.restrictCell(iT, Z), 
                                                                       m_lasxcurl.restrictCell(iT, E_k)}, 2);
  Eigen::MatrixXd L2bkt_Z = nonlinear_coeff * L2v_Bkt_2[0];
  Eigen::MatrixXd L2bkt_E = nonlinear_coeff * L2v_Bkt_2[1];

  //------------------------------------------------------------------------------
  // Block matrices of F
  //------------------------------------------------------------------------------

  Eigen::MatrixXd A
    = laxcurl_L2
      + dt * (1-theta) * 
      (
        0.5 * dt * theta * L2CurlA_ebkt
        + 0.25 * dt * theta * L2ebktAA_ebkt
      )
      + dt * theta * 
      (
        0.5 * dt * theta * L2CurlX_ebkt
        + dt * theta * laxdiv_L2Curl_both
        + dt * theta * L2Curl_left_ebkt_Y.transpose()
        - 0.5 * dt * dt * theta * theta * L2Curl_left_ebkt_E.transpose()
        + 0.25 * dt * theta * L2ebktXX_ebkt
        + dt * theta * L2Curl_left_ebkt_X
        + dt * theta * ebkt_Y.transpose() * L2_nlProduct_T * ebkt_X
        - 0.5 * dt * dt * theta * theta * L2_ebkt_E * ebkt_X
        - 0.5 * dt * dt * theta * theta * L2Curl_left_ebkt_E
        - 0.5 * dt * dt * theta * theta * ebkt_Y.transpose() * L2_ebkt_E.transpose()
        + 0.25 * dt * dt * dt * theta * theta * theta * L2_ebkt_E * ebkt_E
      );

  Eigen::MatrixXd B 
    = dt * laxcurl_L2Grad_right
      + dt * L2bkt_Z
      - dt * dt * (1-theta) * theta * L2bkt_E;

  Eigen::MatrixXd C 
    = laxcurl_L2Grad_right.transpose()
      + L2bkt_Z.transpose()
      + (1-theta) * dt * theta *
      (
        - L2bkt_E.transpose()
        + L2Ei_bkt.transpose()
      );

  size_t sizeSys1 = m_lasxcurl.dimensionCell(iT);
  size_t sizeSys2 = m_lasxgrad.dimensionCell(iT);
  size_t sizeSys_T = sizeSys1 + sizeSys2;

  SystemVectors<Eigen::MatrixXd> to_assemble{{}, {}};

  if (Res1_DF2 == 1){ // Find F(Elambda_k)
    to_assemble.vectors.emplace_back(Eigen::VectorXd(sizeSys_T));

    // Nonlinear problem evalulated at kth iteration
    to_assemble.vectors[0].head(sizeSys1) = A * E_k_T + B * lambda_k_T;
    to_assemble.vectors[0].segment(sizeSys1, sizeSys2) = C * E_k_T;

  } else if (Res1_DF2 == 2) {// Find DF_k
    Eigen::VectorXd ebkt_E_E = ebkt_E * E_k_T;
    Eigen::VectorXd ebkt_X_E = ebkt_X * E_k_T;
    // L2 product matrices with brackets
    std::vector<Eigen::MatrixXd> L2v_epsBkt_list_DF = L2v_epsBkt(iT, crossij_Pot_T, {E_k_T.transpose() * L2Curl_nlProduct_T,
                                                                                  ebkt_X_E.transpose() * L2_nlProduct_T,
                                                                                  ebkt_E_E.transpose() * L2_nlProduct_T});
    Eigen::MatrixXd L2CurlE_ebkt = nonlinear_coeff * L2v_epsBkt_list_DF[0];
    Eigen::MatrixXd L2ebktXE_ebkt = L2v_epsBkt_list_DF[1];
    Eigen::MatrixXd L2ebktEE_ebkt = L2v_epsBkt_list_DF[2];
    Eigen::MatrixXd L2bkt_lambda = nonlinear_coeff * L2v_Bkt(iT, int_PciPcjPgk, {m_lasxgrad.restrictCell(iT, lambda_k)}, 3)[0];

    //------------------------------------------------------------------------------
    // Block matrices of DF
    //------------------------------------------------------------------------------
    A +=  
        dt * theta * 
        (
          - 0.5 * dt * dt * theta * theta * L2CurlE_ebkt
          - 0.5 * dt * dt * theta * theta * L2ebktXE_ebkt
          - 0.5 * dt * dt * theta * theta * L2Curl_left_ebkt_E 
          - 0.5 * dt * dt * theta * theta * ebkt_Y.transpose() * L2_ebkt_E.transpose()
          + 0.25 * dt * dt * dt * theta * theta * theta * (L2_ebkt_E * ebkt_E + L2ebktEE_ebkt)
        ) 
        - dt * dt * (1-theta) * theta * L2bkt_lambda;

    C += - (1-theta) * dt * theta * L2E_bkt.transpose();

    to_assemble.systems.emplace_back(Eigen::MatrixXd::Zero(sizeSys_T, sizeSys_T));
    to_assemble.systems[0].topLeftCorner(sizeSys1, sizeSys1) = A;
    to_assemble.systems[0].block(0, sizeSys1, sizeSys1, sizeSys2) = B;
    to_assemble.systems[0].block(sizeSys1, 0, sizeSys2, sizeSys1) = C;

    //------------------------------------------------------------------------------
    // Cell unknowns of RHS b-F(X_k)
    //------------------------------------------------------------------------------

    to_assemble.vectors.emplace_back(Eigen::VectorXd::Zero(sizeSys_T));
    // Adding cell unknowns for static condensation
    to_assemble.vectors[0].segment(m_lasxcurl.dimensionCell(iT)-m_nloc_sc_E[iT], m_nloc_sc_E[iT]) += m_lasxcurl.restrictCell(iT, m_b_k.head(m_lasxcurl.dimension())).tail(m_nloc_sc_E[iT]);
    to_assemble.vectors[0].segment(m_lasxcurl.dimensionCell(iT)+m_lasxgrad.dimensionCell(iT)-m_nloc_sc_lambda[iT], m_nloc_sc_lambda[iT]) += m_lasxgrad.restrictCell(iT, m_b_k.tail(m_lasxgrad.dimension())).tail(m_nloc_sc_lambda[iT]);
  } else {
    std::cerr << "[main] ERROR: Unknown Newton calculation" << std::endl;
    exit(1);    
  }

  return to_assemble;
}

//------------------------------------------------------------------------------

LocalStaticCondensation YangMills::_compute_static_condensation(const size_t & iT) const
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  // Dimensions
  size_t dim_E = m_lasxcurl.dimensionCell(iT) - m_nloc_sc_E(iT);     // number of E unknowns after SC
  size_t dim_lambda = m_lasxgrad.dimensionCell(iT) - m_nloc_sc_lambda(iT);     // number of lambda unknowns after SC
  size_t dim_dofs = dim_E + dim_lambda;      // nb of dofs remaining after SC (including Lagrange multiplier)
  size_t dim_sc = m_nloc_sc_E(iT) + m_nloc_sc_lambda(iT);      // nb of SC dofs in total

  // Creation of permutation matrix
  Eigen::MatrixXd Perm = Eigen::MatrixXd::Zero(dim_dofs+dim_sc, dim_dofs+dim_sc);
  Perm.topLeftCorner(dim_E, dim_E) = Eigen::MatrixXd::Identity(dim_E, dim_E);
  Perm.block(dim_E, dim_E + m_nloc_sc_E(iT), dim_lambda, dim_lambda) = Eigen::MatrixXd::Identity(dim_lambda, dim_lambda);
  Perm.block(dim_E + dim_lambda, dim_E, m_nloc_sc_E(iT), m_nloc_sc_E(iT)) = Eigen::MatrixXd::Identity(m_nloc_sc_E(iT), m_nloc_sc_E(iT));
  Perm.block(dim_E + dim_lambda + m_nloc_sc_E(iT), dim_E + m_nloc_sc_E(iT) + dim_lambda, m_nloc_sc_lambda(iT), m_nloc_sc_lambda(iT))
      = Eigen::MatrixXd::Identity(m_nloc_sc_lambda(iT), m_nloc_sc_lambda(iT));

  // Creation of global DOFs for system: IT_sys contains the skeletal dofs of E and of lambda
  std::vector<size_t> IT_sys(dim_dofs, 0);
  auto IT_lasxcurl = m_lasxcurl.globalDOFIndices(T);
  auto IT_lasxgrad = m_lasxgrad.globalDOFIndices(T);
  auto it_IT_sys = std::copy(IT_lasxcurl.begin(), IT_lasxcurl.begin()+dim_E, IT_sys.begin()); // put skeletal DOFs of E in IT_sys
  size_t offset = m_lasxcurl.dimension() - nbSCDOFs_E();     // nb total of skeletal DOFs for E (where global skeletal dofs of lambda start)
  std::transform(IT_lasxgrad.begin(), IT_lasxgrad.begin()+dim_lambda, it_IT_sys, [&offset](const size_t & index) { return index + offset; });
  
  // Creation of global DOFs for SC operator: IT_sc contains global cell dofs of E (offset to start at 0) and lambda (offset to start after those of E)
  std::vector<size_t> IT_sc(dim_sc, 0);  
  if (dim_sc>0){
    auto it_IT_sc = std::transform(IT_lasxcurl.begin()+dim_E, IT_lasxcurl.end(), IT_sc.begin(), [&offset](const size_t & index) { return index - offset; });
    size_t end_E = nbSCDOFs_E();
    offset = (m_lasxgrad.dimension() - nbSCDOFs_lambda());
    std::transform(IT_lasxgrad.begin()+dim_lambda, IT_lasxgrad.end(), it_IT_sc, [&offset, &end_E](const size_t & index) { return index - offset + end_E; });
  }  

  return LocalStaticCondensation(Perm, IT_sys, IT_sc);
}


//------------------------------------------------------------------------------

void YangMills::_assemble_local_newton(
                                      size_t iT,
                                      const SystemVectors<Eigen::MatrixXd> & lsT,
                                      std::vector<std::list<Eigen::Triplet<double> > > & triplets_sys,
                                      std::vector<Eigen::VectorXd> & vecs
                                      )
{
  auto & [systems, vectors] = lsT;
  // Get information for local static condensation
  LocalStaticCondensation locSC = _compute_static_condensation(iT);

  Eigen::MatrixXd AT_sys, AT_sc;
  Eigen::VectorXd bT_sys, bT_sc;
  std::tie(AT_sys, bT_sys, AT_sc, bT_sc) = locSC.compute(std::make_pair(systems[0], vectors[0]));

  std::vector<size_t> IT_sys = locSC.globalDOFs_gl();

  for (size_t i = 0; i < locSC.dim_gl(); i++){
    for (size_t j = 0; j < locSC.dim_gl(); j++){
      triplets_sys[0].emplace_back(IT_sys[i], IT_sys[j], AT_sys(i, j));
    }
    vecs[0](IT_sys[i]) += bT_sys(i);
  }

  std::vector<size_t> IT_sc = locSC.globalDOFs_sc();

  for (size_t i = 0; i < locSC.dim_sc(); i++){
    for (size_t j = 0; j < locSC.dim_gl(); j++){
      triplets_sys[1].emplace_back(IT_sc[i], IT_sys[j], AT_sc(i, j));
    }
    vecs[1](IT_sc[i]) += bT_sc(i);
  }
}

//------------------------------------------------------------------------------

std::vector<Eigen::MatrixXd>
YangMills::_compute_local_contribution(
                                      size_t iT,
                                      double dt
                                        )
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  //------------------------------------------------------------------------------
  // Local product matrices
  //------------------------------------------------------------------------------
  // Mass matrices
  MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*m_ddrcore.degree() + 2);
  Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_2kp2);
  Eigen::MatrixXd mass_Pkpo_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polykpo, int_mono_2kp2);

  Eigen::MatrixXd L2laxcurl = m_lasxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd L2laxdiv = m_lasxdiv.computeL2Product(iT, m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd L2laxgrad = m_lasxgrad.computeL2Product(iT, m_stab_par, mass_Pkpo_T);

  return std::vector{L2laxcurl, L2laxdiv, L2laxgrad};
}

//------------------------------------------------------------------------------

void YangMills::_assemble_local_contribution(
                                          size_t iT,
                                          const std::vector<Eigen::MatrixXd> & lsT,
                                          std::vector<std::list<Eigen::Triplet<double> > > & triplets_sys,
                                          std::vector<Eigen::VectorXd> & vecs
                                          )
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);
  Eigen::MatrixXd L2laxcurl = lsT[0];
  Eigen::MatrixXd L2laxdiv = lsT[1];
  Eigen::MatrixXd L2laxgrad = lsT[2];

  std::vector<size_t> IT_sysg = m_lasxgrad.globalDOFIndices(T);
  std::vector<size_t> IT_sysc = m_lasxcurl.globalDOFIndices(T);
  std::vector<size_t> IT_sysd = m_lasxdiv.globalDOFIndices(T);

  // Global L2 product matrices
  for (size_t i = 0; i < IT_sysc.size(); i++){
    for (size_t j = 0; j < IT_sysc.size(); j++){
      triplets_sys[0].emplace_back(IT_sysc[i], IT_sysc[j], L2laxcurl(i, j));
    }
  }

  for (size_t i = 0; i < IT_sysd.size(); i++){
    for (size_t j = 0; j < IT_sysd.size(); j++){
      triplets_sys[1].emplace_back(IT_sysd[i], IT_sysd[j], L2laxdiv(i, j));
    }
  }

  for (size_t i = 0; i < IT_sysg.size(); i++){
    for (size_t j = 0; j < IT_sysg.size(); j++){
      triplets_sys[2].emplace_back(IT_sysg[i], IT_sysg[j], L2laxgrad(i, j));
    }
  }
}

//------------------------------------------------------------------------------

Scalar3Tensor YangMills::_compute_detnij_Pot_PkF(const size_t iF) const
{
  size_t dim_PolykF = PolynomialSpaceDimension<Face>::Poly(m_ddrcore.degree());
  size_t dim_xcurl_F = m_sxcurl.dimensionFace(iF);

  Eigen::MatrixXd pot_F = m_sxcurl.facePotential(iF);

  Scalar3Tensor detnij_Pot_PkF(boost::extents[dim_xcurl_F][dim_xcurl_F][dim_PolykF]);
  std::fill_n(detnij_Pot_PkF.data(), detnij_Pot_PkF.num_elements(), 0.);

  Scalar3Tensor detnij_PkF = _compute_detnij_PkF(iF);

  for (size_t k = 0; k < dim_PolykF; k++){
    // The matrix giving the kth coefficients (on PolykF) of the integrals of (phi_i x phi_j).nF
    Eigen::MatrixXd detnij = slice(detnij_PkF, 2, k);
    // The kth coefficient (on PolykF) of the integrals of potentials (P_curl e_i x P_curl e_j).nF
    slice(detnij_Pot_PkF, 2, k) = pot_F.transpose() * detnij * pot_F;
  } 
  return detnij_Pot_PkF;
}

//------------------------------------------------------------------------------

Scalar3Tensor YangMills::_compute_detnij_PkF(const size_t iF) const
{
  Face & F = *m_ddrcore.mesh().face(iF);
  size_t dim_PolykF = PolynomialSpaceDimension<Face>::Poly(m_ddrcore.degree());
  size_t dim_Polyk2F = m_ddrcore.faceBases(iF).Polyk2->dimension();

  Scalar3Tensor detnij_PkF(boost::extents[dim_Polyk2F][dim_Polyk2F][dim_PolykF]);
  std::fill_n(detnij_PkF.data(), detnij_PkF.num_elements(), 0.);

  QuadratureRule quad_3k_F = generate_quadrature_rule(F, 3*m_ddrcore.degree());
  auto basis_Pk_F_quad = evaluate_quad<Function>::compute(*m_sxcurl.faceBases(iF).Polyk, quad_3k_F);
  const Eigen::Vector3d & nF = F.normal();
  auto & Pk2 = *(m_ddrcore.faceBases(iF).Polyk2);
  Eigen::MatrixXd mass_Pk_F = compute_gram_matrix(basis_Pk_F_quad, quad_3k_F);

  // For each i, j, calculate the projection of (phi_i x phi_j).nF and set the coefficients
  for (size_t i = 0; i < dim_Polyk2F; i++){
    for (size_t j = 0; j < i; j++){
      auto cross_ij = [&Pk2, i, j, nF](const Eigen::Vector3d & x)-> double {return Pk2.function(i, x).cross(Pk2.function(j, x)).dot(nF);};
      Eigen::VectorXd proj_ij = l2_projection(cross_ij, *m_ddrcore.faceBases(iF).Polyk, quad_3k_F, basis_Pk_F_quad, mass_Pk_F);
      for (size_t k = 0; k < dim_PolykF; k++){
        detnij_PkF[i][j][k] = proj_ij[k];
        detnij_PkF[j][i][k] = -proj_ij[k];
      }
    }
  }
  return detnij_PkF;
}

//------------------------------------------------------------------------------

Scalar3Tensor YangMills::_compute_crossij_T(const size_t iT) const
{
  Cell & T = *m_ddrcore.mesh().cell(iT);

  size_t dim_Polyk3T = m_ddrcore.cellBases(iT).Polyk3->dimension();
  size_t dim_GolykmoT = PolynomialSpaceDimension<Cell>::Goly(m_ddrcore.degree()-1);
  size_t dim_GolyComplkT = PolynomialSpaceDimension<Cell>::GolyCompl(m_ddrcore.degree());

  Scalar3Tensor crossij_T(boost::extents[dim_Polyk3T][dim_Polyk3T][dim_GolykmoT+dim_GolyComplkT]);
  std::fill_n(crossij_T.data(), crossij_T.num_elements(), 0);

  MonomialCellIntegralsType int_mono_3k = IntegrateCellMonomials(T, 3*m_ddrcore.degree());
  Scalar3Tensor triple_int_Gkmo = tripleInt(T, *m_sxcurl.cellBases(iT).Golykmo, *m_ddrcore.cellBases(iT).Polyk3, int_mono_3k);
  Scalar3Tensor triple_int_GCk = tripleInt(T, *m_sxcurl.cellBases(iT).GolyComplk, *m_ddrcore.cellBases(iT).Polyk3, int_mono_3k);
  
  Eigen::MatrixXd mass_Gkmo_T = GramMatrix(T, *m_sxcurl.cellBases(iT).Golykmo, int_mono_3k);
  Eigen::MatrixXd mass_GCk_T = GramMatrix(T, *m_sxcurl.cellBases(iT).GolyComplk, int_mono_3k);
  Eigen::LDLT<Eigen::MatrixXd> cholesky_mass_Gkmo(mass_Gkmo_T);
  Eigen::LDLT<Eigen::MatrixXd> cholesky_mass_GCk(mass_GCk_T);

  // For each i, j, calculate the projection of (phi_i x phi_j) and set the coefficients
  for (size_t i = 0; i < dim_Polyk3T; i++){
    for (size_t j = 0; j < i; j++){
      Eigen::VectorXd crossij_Gkmo = slice(triple_int_Gkmo, 1, i, 2, j);
      Eigen::VectorXd crossij_GCk = slice(triple_int_GCk, 1, i, 2, j);
      Eigen::VectorXd proj_crossij_Gkmo = cholesky_mass_Gkmo.solve(crossij_Gkmo);
      Eigen::VectorXd proj_crossij_GCk = cholesky_mass_GCk.solve(crossij_GCk);
      for (size_t k = 0; k < dim_GolykmoT; k++){
        crossij_T[i][j][k] = proj_crossij_Gkmo[k];
        crossij_T[j][i][k] = -proj_crossij_Gkmo[k];
      }
      for (size_t k = 0; k < dim_GolyComplkT; k++){
        crossij_T[i][j][dim_GolykmoT+k] = proj_crossij_GCk[k];
        crossij_T[j][i][dim_GolykmoT+k] = -proj_crossij_GCk[k];
      }
    }
  }
  return crossij_T;
}

//------------------------------------------------------------------------------

Scalar3Tensor YangMills::_compute_crossij_Pot_T(const size_t iT) const
{
  size_t dim_xcurl_T = m_sxcurl.dimensionCell(iT);
  Eigen::MatrixXd pot_T = m_sxcurl.cellPotential(iT);

  if (m_nl_par == 0){
    size_t dim_xdiv_DofsCell_T = m_sxdiv.numLocalDofsCell();

    Scalar3Tensor crossij_Pot_T(boost::extents[dim_xcurl_T][dim_xcurl_T][dim_xdiv_DofsCell_T]);
    std::fill_n(crossij_Pot_T.data(), crossij_Pot_T.num_elements(), 0.);

    if (m_ddrcore.degree()>0){

      Scalar3Tensor crossij_T = _compute_crossij_T(iT);

      for (size_t k = 0; k < dim_xdiv_DofsCell_T; k++){
        // The matrix giving the kth coefficients (on Golykmo/GolyComplk) of the integrals of (phi_i x phi_j)
        Eigen::MatrixXd crossij = slice(crossij_T, 2, k);        
        slice(crossij_Pot_T, 2, k) = pot_T.transpose() * crossij * pot_T;

      }
    }
    return crossij_Pot_T;

  } else if (m_nl_par == 1){
    Cell & T = *m_ddrcore.mesh().cell(iT);

    size_t dim_Poly2k3T = m_laddrcore.P2k3(iT).dimension();

    Scalar3Tensor crossij_Pot_T(boost::extents[dim_xcurl_T][dim_xcurl_T][dim_Poly2k3T]);
    std::fill_n(crossij_Pot_T.data(), crossij_Pot_T.num_elements(), 0.);

    MonomialCellIntegralsType int_mono_3k = IntegrateCellMonomials(T, 3*m_ddrcore.degree());
    Scalar3Tensor triple_int_P2k = tripleInt(T, m_laddrcore.P2k3(iT), *m_ddrcore.cellBases(iT).Polyk3, int_mono_3k);

    for (size_t i = 0; i < dim_Poly2k3T; i++){      
      slice(crossij_Pot_T, 2, i) = pot_T.transpose() * slice(triple_int_P2k, 0, i) * pot_T;
    }

    return crossij_Pot_T;
  } else {
    std::cerr << "[main] ERROR: Unknown bracket discretisation" << std::endl;
    exit(1);
  }
}

//------------------------------------------------------------------------------

Scalar3Tensor YangMills::_integral_Pot_ijk(size_t iT) const
{
  Cell & T = *m_ddrcore.mesh().cell(iT);
  
  size_t dim_xcurl_T = m_sxcurl.dimensionCell(iT);
  size_t dim_xgrad_T = m_sxgrad.dimensionCell(iT);

  size_t dim_Pk3 = m_ddrcore.cellBases(iT).Polyk3->dimension();
  size_t dim_Pkpo = m_ddrcore.cellBases(iT).Polykpo->dimension();

  // integral of the basis functions in T: tripleInt returns ordering [k][i][j] for (int \psi_k * (\phi_i . \phi_j))
  MonomialCellIntegralsType int_mono_3kp1 = IntegrateCellMonomials(T, 3*m_ddrcore.degree() + 1);
  Scalar3Tensor integral_basis = tripleInt(T, *m_ddrcore.cellBases(iT).Polykpo, *m_ddrcore.cellBases(iT).Polyk3, int_mono_3kp1);

  // integral over the basis in T of potentials (P_curl e_i . P_curl e_j)P_grad a_k
  Scalar3Tensor integral_PcPcPg(boost::extents[dim_xcurl_T][dim_xcurl_T][dim_xgrad_T]);
  std::fill_n(integral_PcPcPg.data(), integral_PcPcPg.num_elements(), 0.);
  Eigen::MatrixXd xgrad_pot = m_sxgrad.cellPotential(iT);
  Eigen::MatrixXd xcurl_pot = m_sxcurl.cellPotential(iT);

  for (size_t i = 0; i < dim_xcurl_T; i++){
    // Calculate (P_curl e_i . \phi_i)\phi_j
    Eigen::MatrixXd integral_jk = Eigen::MatrixXd::Zero(dim_Pk3, dim_Pkpo);
    for (size_t l=0; l<dim_Pkpo; l++){
      integral_jk.col(l) += slice(integral_basis, 0, l).transpose() * xcurl_pot.col(i);
    }
    // Calculate (P_curl e_i . P_curl e_j)P_grad a_k
    for (size_t j = 0; j <= i; j++){
      for (size_t k = 0; k < dim_xgrad_T; k++){
        integral_PcPcPg[i][j][k] = xcurl_pot.col(j).transpose() * integral_jk * xgrad_pot.col(k);
        integral_PcPcPg[j][i][k] = integral_PcPcPg[i][j][k];
      }
    }
  }

  return integral_PcPcPg;
}

//------------------------------------------------------------------------------

std::vector<Eigen::MatrixXd> YangMills::L2v_Bkt(size_t iT, 
                                   Scalar3Tensor & intPciPcjPgk, 
                                   const std::vector<Eigen::VectorXd> & vec_list, 
                                   const size_t & entry) const
{
  size_t dim_laxcurl_T = m_lasxcurl.dimensionCell(iT);
  size_t dim_laxgrad_T = m_lasxgrad.dimensionCell(iT);

  size_t dim_xcurl_T = m_sxcurl.dimensionCell(iT);
  size_t dim_xgrad_T = m_sxgrad.dimensionCell(iT);

  size_t dim_la = m_liealg.dimension();

  size_t num_vec = vec_list.size();
  std::vector<Eigen::MatrixXd> L2v_Bkt_list;
  const std::vector<Eigen::MatrixXd> & C_IJ = m_liealg.structureConst();
  const Eigen::MatrixXd & massMatrix = m_liealg.massMatrix();

  if (entry == 1){
    // Compute the Lie algebra matrices to take the kronecker product with: C^K_IJ <e_K, e_L> with L fixed
    std::vector<Eigen::MatrixXd> C_JK_massMatrix_KL(dim_la, Eigen::MatrixXd::Zero(dim_la, dim_la));
    for (size_t L=0; L < dim_la; L++){
      for (size_t K=0; K < dim_la; K++){
        C_JK_massMatrix_KL[L] += massMatrix(L, K) * C_IJ[K];
      }
    }
    for (size_t iv = 0; iv < num_vec; iv++){
      // Read the laxcurl vector v as a dim_la by dim_xcurl matrix [I][i]
      const Eigen::Map<const Eigen::MatrixXd> M_v(vec_list[iv].data(), dim_la, dim_xcurl_T);
      // Multiarray v_I_intjk[I][j][k] representing (v_i)_I (P_curl e_i . P_curl e_j)P_grad a_k
      Scalar3Tensor v_I_intjk(boost::extents[dim_la][dim_xcurl_T][dim_xgrad_T]);
      std::fill_n(v_I_intjk.data(), v_I_intjk.size(), 0.);

      for (size_t j = 0; j < dim_xcurl_T; j++){
        slice(v_I_intjk, 1, j) = M_v * slice(intPciPcjPgk, 1, j);
      }
      // Calculate the Lie algebra operator: take the kronecker product and contract (v_i)_I (P_curl e_i . P_curl e_j)P_grad a_k on I-L with C^K_IJ <e_K, e_L>
      L2v_Bkt_list.emplace_back(Eigen::MatrixXd::Zero(dim_laxcurl_T, dim_laxgrad_T));
      for (size_t I = 0; I < dim_la; I++){
        L2v_Bkt_list[iv] += Eigen::KroneckerProduct(slice(v_I_intjk, 0, I), C_JK_massMatrix_KL[I]);
      }
    }
    return L2v_Bkt_list;
  }else if (entry == 2){
    // Compute the Lie algebra matrices to take the kronecker product with: C^K_IJ <e_K, e_L> with J fixed
    std::vector<Eigen::MatrixXd> C_JK_massMatrix_KL(dim_la, Eigen::MatrixXd::Zero(dim_la, dim_la));
    for (size_t L=0; L < dim_la; L++){
      Eigen::MatrixXd C_IJK_massMatrix_K(Eigen::MatrixXd::Zero(dim_la, dim_la));
      for (size_t K=0; K < dim_la; K++){
        C_IJK_massMatrix_K += massMatrix(L, K) * C_IJ[K];
      }
      for (size_t J=0; J <  dim_la; J++){
        C_JK_massMatrix_KL[J].row(L) = C_IJK_massMatrix_K.row(J);
      }
    }

    for (size_t iv = 0; iv < num_vec; iv++){
      // Read the laxcurl vector v as a dim_la by dim_xcurl matrix [J][j]
      const Eigen::Map<const Eigen::MatrixXd> M_v(vec_list[iv].data(), dim_la, dim_xcurl_T);
      // Multiarray v_J_intik[J][i][k] representing (v_j)_J (P_curl e_i . P_curl e_j)P_grad a_k
      Scalar3Tensor v_J_intik(boost::extents[dim_la][dim_xcurl_T][dim_xgrad_T]);
      std::fill_n(v_J_intik.data(), v_J_intik.size(), 0.);

      for (size_t i = 0; i < dim_xcurl_T; i++){
        slice(v_J_intik, 1, i) = M_v * slice(intPciPcjPgk, 0, i);
      }
      // Calculate the Lie algebra operator: take the kronecker product and contract (v_j)_J (P_curl e_i . P_curl e_j)P_grad a_k on J-J with C^K_IJ <e_K, e_L>
      L2v_Bkt_list.emplace_back(Eigen::MatrixXd::Zero(dim_laxcurl_T, dim_laxgrad_T));
      for (size_t J = 0; J < dim_la; J++){
        L2v_Bkt_list[iv] += Eigen::KroneckerProduct(slice(v_J_intik, 0, J), C_JK_massMatrix_KL[J]);
      }
    }
    return L2v_Bkt_list;
  }else if (entry == 3){
    // Compute the Lie algebra matrices to take the kronecker product with: C^K_IJ <e_K, e_L> with I fixed
    std::vector<Eigen::MatrixXd> C_JK_massMatrix_KL(dim_la, Eigen::MatrixXd::Zero(dim_la, dim_la));
    for (size_t L=0; L < dim_la; L++){
      Eigen::MatrixXd C_IJK_massMatrix_K(Eigen::MatrixXd::Zero(dim_la, dim_la));
      for (size_t K=0; K < dim_la; K++){
        C_IJK_massMatrix_K += massMatrix(L, K) * C_IJ[K];
      }
      for (size_t I=0; I<dim_la; I++){
        C_JK_massMatrix_KL[I].row(L) = C_IJK_massMatrix_K.row(I);
      }
    }
    for (size_t iv = 0; iv < num_vec; iv++){
      // Read the laxgrad vector q as a dim_la by dim_xgrad matrix [K][k]
      const Eigen::Map<const Eigen::MatrixXd> M_v(vec_list[iv].data(), dim_la, dim_xgrad_T);
      // Multiarray q_K_intij[K][i][j] representing (q_k)_K (P_curl e_i . P_curl e_j)P_grad a_k
      Scalar3Tensor q_K_intij(boost::extents[dim_la][dim_xcurl_T][dim_xcurl_T]);
      std::fill_n(q_K_intij.data(), q_K_intij.size(), 0.);

      for (size_t i = 0; i < dim_xcurl_T; i++){
        slice(q_K_intij, 1, i) = M_v * slice(intPciPcjPgk, 0, i).transpose();
      }
      // Calculate the Lie algebra operator: take the kronecker product and contract (q_k)_K (P_curl e_i . P_curl e_j)P_grad a_k on K-L with C^K_IJ <e_K, e_L>
      L2v_Bkt_list.emplace_back(Eigen::MatrixXd::Zero(dim_laxcurl_T, dim_laxcurl_T));
      for (size_t K = 0; K < dim_la; K++){
        L2v_Bkt_list[iv] += Eigen::KroneckerProduct(slice(q_K_intij, 0, K), C_JK_massMatrix_KL[K]);
      }
    }
    return L2v_Bkt_list;
  }else{
    std::cout << "[Main] Invalid entry in L2v_Bkt" << std::endl;
    exit(1);
  }
}

//------------------------------------------------------------------------------

std::vector<Eigen::MatrixXd> YangMills::epsBkt_v(size_t iT, Scalar3Tensor & crossij_Pot_T, const std::vector<Eigen::VectorXd> & vec_list) const
{
  Cell & T = *m_ddrcore.mesh().cell(iT);
  size_t dim_la = m_liealg.dimension();
  // Compute the Lie algebra matrices to take the kronecker product with: C^K_IJ with I fixed
  const std::vector<Eigen::MatrixXd> & C_IJ =  m_liealg.structureConst();
  std::vector<Eigen::MatrixXd> C_JK(dim_la, Eigen::MatrixXd::Zero(dim_la, dim_la));
  for (size_t I=0; I < dim_la; I++){
    for (size_t K=0; K < dim_la; K++){
      C_JK[I].col(K) = C_IJ[K].row(I).transpose();
    }
  }

  if (m_nl_par == 0){
  size_t dim_xdiv_DofsFace_F = m_sxdiv.numLocalDofsFace();
  size_t dim_laxdiv_DofsFace_F = m_lasxdiv.numLocalDofsFace();
  size_t dim_laxcurl_T = m_lasxcurl.dimensionCell(iT);

  size_t num_vec = vec_list.size();

  // List of cell operators representing *[v,.] on T
  std::vector<Eigen::MatrixXd> epsBkt_v_list(num_vec, Eigen::MatrixXd::Zero(m_lasxdiv.dimensionCell(iT), dim_laxcurl_T));

  // Add the effect of each face F operator to the cell operator
  for (auto Fp : T.get_faces()){
    Face & F = *Fp;
    size_t iF = F.global_index();
    size_t offset_F = m_lasxdiv.localOffset(T, F);
    size_t dim_xcurl_F = m_sxcurl.dimensionFace(iF);
    Scalar3Tensor detnij_Pot_PkF = _compute_detnij_Pot_PkF(iF);

    std::vector<Eigen::VectorXd> vec_list_F(num_vec);
    for (size_t iv = 0; iv < num_vec; iv++){
      vec_list_F[iv] = m_lasxcurl.restrict(F, vec_list[iv]);
    }

    for (size_t iv = 0; iv < num_vec; iv++){
      // Read the laxcurl vector v as a dim_la by dim_xcurl matrix [I][i]
      const Eigen::Map<const Eigen::MatrixXd> M_v(vec_list_F[iv].data(), dim_la, dim_xcurl_F);

      // Multiarray v_I_det_jk[I][j][k] representing (v_i)_I (P_curl e_i x P_curl e_j).nF (k coefficients of result)
      Scalar3Tensor v_I_det_jk(boost::extents[dim_la][dim_xcurl_F][dim_xdiv_DofsFace_F]);
      std::fill_n(v_I_det_jk.data(), v_I_det_jk.size(), 0.);

      for (size_t j=0; j < dim_xcurl_F; j++){
        slice(v_I_det_jk, 1, j) = M_v * slice(detnij_Pot_PkF, 1, j);
      }
      // Calculate the Lie algebra operator: contract (v_i)_I (P_curl e_i x P_curl e_j).nF on I with C^K_IJ, then extend from face to cell and sum
      for (size_t I=0; I < dim_la; I++){
        m_lasxcurl.extendOperator(T, F, epsBkt_v_list[iv].block(offset_F, 0, dim_laxdiv_DofsFace_F, dim_laxcurl_T), Eigen::KroneckerProduct(slice(v_I_det_jk, 0, I), C_JK[I]).transpose());
      }
    }//for iv
  }//for Fp

  if (m_ddrcore.degree()>0){
    size_t offset_T = m_lasxdiv.localOffset(T);
    size_t dim_laxdiv_DofsCell_T = m_lasxdiv.numLocalDofsCell();
    size_t dim_xdiv_DofsCell_T = m_sxdiv.numLocalDofsCell();
    size_t dim_xcurl_T = m_sxcurl.dimensionCell(iT);

    std::vector<Eigen::VectorXd> vec_list_T(num_vec);
    for (size_t iv = 0; iv < num_vec; iv++){
      vec_list_T[iv] = m_lasxcurl.restrict(T, vec_list[iv]);
    }

    for (size_t iv = 0; iv < num_vec; iv++){
      // Read the laxcurl vector v as a dim_la by dim_xcurl matrix [I][i]
      const Eigen::Map<const Eigen::MatrixXd> M_v(vec_list_T[iv].data(), dim_la, dim_xcurl_T);
      // Multiarray v_I_det_jk[I][j][k] representing (v_i)_I (P_curl e_i x P_curl e_j)^k (k coefficients of result)
      Scalar3Tensor v_I_cross_jk(boost::extents[dim_la][dim_xcurl_T][dim_xdiv_DofsCell_T]);
      std::fill_n(v_I_cross_jk.data(), v_I_cross_jk.size(), 0.);

      for (size_t j=0; j < dim_xcurl_T; j++){
        slice(v_I_cross_jk, 1, j) = M_v * slice(crossij_Pot_T, 1, j);
      }
      // Calculate the Lie algebra operator: contract (v_i)_I (P_curl e_i x P_curl e_j)^k on I with C^K_IJ, and set the section on the cell Dofs
      for (size_t I=0; I < dim_la; I++){
        epsBkt_v_list[iv].block(offset_T, 0, dim_laxdiv_DofsCell_T, dim_laxcurl_T) += Eigen::KroneckerProduct(slice(v_I_cross_jk, 0, I), C_JK[I]).transpose();
      }
    }
  }
  return epsBkt_v_list;
  } else if (m_nl_par == 1){

    size_t dim_xcurl_T = m_sxcurl.dimensionCell(iT);
    size_t dim_P2k3_T = m_laddrcore.P2k3(iT).dimension();
    size_t dim_laxcurl_T = m_lasxcurl.dimensionCell(iT);

    size_t num_vec = vec_list.size();
    std::vector<Eigen::VectorXd> vec_list_T(num_vec);
    for (size_t iv = 0; iv < num_vec; iv++){
      vec_list_T[iv] = m_lasxcurl.restrict(T, vec_list[iv]);
    }

    // List of cell operators representing *[v,.] on T
    std::vector<Eigen::MatrixXd> epsBkt_v_list(num_vec, Eigen::MatrixXd::Zero(dim_la*dim_P2k3_T, dim_laxcurl_T));

    for (size_t iv = 0; iv < num_vec; iv++){
      // Read the laxcurl vector v as a dim_la by dim_xcurl matrix [I][i]
      const Eigen::Map<const Eigen::MatrixXd> M_v(vec_list_T[iv].data(), dim_la, dim_xcurl_T);
      // Multiarray v_I_det_jk[I][j][k] representing (v_i)_I (P_curl e_i x P_curl e_j)^k (k coefficients of result)
      Scalar3Tensor v_I_cross_jk(boost::extents[dim_la][dim_xcurl_T][dim_P2k3_T]);
      std::fill_n(v_I_cross_jk.data(), v_I_cross_jk.size(), 0.);

      for (size_t j=0; j < dim_xcurl_T; j++){
        slice(v_I_cross_jk, 1, j) = M_v * slice(crossij_Pot_T, 1, j);
      }
      // Calculate the Lie algebra operator: contract (v_i)_I (P_curl e_i x P_curl e_j)^k on I with C^K_IJ
      for (size_t I=0; I < dim_la; I++){
        epsBkt_v_list[iv] += Eigen::KroneckerProduct(slice(v_I_cross_jk, 0, I), C_JK[I]).transpose();
      }
    }
    return epsBkt_v_list;
  } else{
    std::cerr << "[main] ERROR: Unknown bracket discretisation" << std::endl;
    exit(1);
  }

}

//------------------------------------------------------------------------------

std::vector<Eigen::MatrixXd> YangMills::L2v_epsBkt(size_t iT, Scalar3Tensor & crossij_Pot_T, const std::vector<Eigen::VectorXd> & v_L2prod_list) const
{
  if (m_nl_par == 0){
    Cell & T = *m_ddrcore.mesh().cell(iT);
    size_t num_vec = v_L2prod_list.size();

    size_t dim_la = m_liealg.dimension();
    size_t dim_xcurl_T = m_sxcurl.dimensionCell(iT);
    size_t dim_laxcurl_T = m_lasxcurl.dimensionCell(iT);

    // List of cell operators representing (v, *[.,.])_L2 on T
    std::vector<Eigen::MatrixXd> L2v_epsBkt_list(num_vec, Eigen::MatrixXd::Zero(dim_laxcurl_T, dim_laxcurl_T));
    const std::vector<Eigen::MatrixXd> & C_IJ = m_liealg.structureConst();

    for (auto Fp : T.get_faces()){
      Face & F = *Fp;
      size_t iF = F.global_index();
      size_t offset_F = m_lasxdiv.localOffset(T, F);
      size_t dim_xcurl_F = m_sxcurl.dimensionFace(iF);
      size_t dim_laxcurl_F = m_lasxcurl.dimensionFace(iF);
      size_t dim_xdiv_DofsFace_F = m_sxdiv.numLocalDofsFace();
      size_t dim_laxdiv_DofsFace_F = m_lasxdiv.numLocalDofsFace();
        
      Scalar3Tensor detnij_Pot_PkF = _compute_detnij_Pot_PkF(iF);
      
      std::vector<Eigen::VectorXd> v_L2prod_DofsCell_F(num_vec, Eigen::VectorXd::Zero(dim_laxdiv_DofsFace_F));
      for (size_t iv = 0; iv < num_vec; iv++){
        v_L2prod_DofsCell_F[iv] = v_L2prod_list[iv].segment(offset_F, dim_laxdiv_DofsFace_F);
      }// for iv

      for (size_t iv = 0; iv < num_vec; iv++){
        // Read the operator (v, .)_L2 as a dim_la by dim_xdiv matrix [L][i]
        const Eigen::Map<const Eigen::MatrixXd> M_vL2(v_L2prod_DofsCell_F[iv].data(), dim_la, dim_xdiv_DofsFace_F);
        // Multiarray vL2_L_cross_ij[L][i][j] representing ((v_i)_I, (e_i x e_j).nF)_L2 <e_I, e_L>
        Scalar3Tensor vL2_L_cross_ij(boost::extents[dim_la][dim_xcurl_F][dim_xcurl_F]);
        std::fill_n(vL2_L_cross_ij.data(), vL2_L_cross_ij.size(), 0.);

        for (size_t j=0; j<dim_xcurl_F; j++){
          slice(vL2_L_cross_ij, 2, j) = M_vL2 * slice(detnij_Pot_PkF, 1, j).transpose();
        }
        // Calculate the Lie algebra operator: take the kronecker product and contract ((v_i)_I, (e_i x e_j).nF)_L2 <e_I, e_L> on L-K with C^K_IJ, then extend from face to cell and sum
        Eigen::MatrixXd L2v_epsBkt_F(Eigen::MatrixXd::Zero(dim_laxcurl_F, dim_laxcurl_F));
        for (size_t L=0; L<dim_la; L++){
          L2v_epsBkt_F += Eigen::KroneckerProduct(slice(vL2_L_cross_ij, 0, L), C_IJ[L]);
        }
        m_lasxcurl.addInnerProductContribution(T, F, L2v_epsBkt_list[iv], L2v_epsBkt_F);
        
      }// for iv
    } // for Fp

    if (m_ddrcore.degree()>0){
      size_t offset_T = m_lasxdiv.localOffset(T);
      size_t dim_laxdiv_DofsCell_T = m_lasxdiv.numLocalDofsCell();
      size_t dim_xdiv_DofsCell_T = m_sxdiv.numLocalDofsCell();

      std::vector<Eigen::VectorXd> v_L2prod_DofsCell_T(num_vec, Eigen::VectorXd::Zero(dim_laxdiv_DofsCell_T));
      for (size_t iv = 0; iv < num_vec; iv++){
        v_L2prod_DofsCell_T[iv] = v_L2prod_list[iv].segment(offset_T, dim_laxdiv_DofsCell_T);
      }// for i
      
      for (size_t iv = 0; iv < num_vec; iv++){
        // Read the operator (v, .)_L2 as a dim_la by dim_xdiv matrix [L][i]
        const Eigen::Map<const Eigen::MatrixXd> M_vL2(v_L2prod_DofsCell_T[iv].data(), dim_la, dim_xdiv_DofsCell_T);
        // Multiarray vL2_L_cross_ij[L][i][j] representing ((v_i)_I, e_i x e_j)_L2 <e_I, e_L>
        Scalar3Tensor vL2_L_cross_ij(boost::extents[dim_la][dim_xcurl_T][dim_xcurl_T]);
        std::fill_n(vL2_L_cross_ij.data(), vL2_L_cross_ij.size(), 0.);

        for (size_t j=0; j<dim_xcurl_T; j++){
          slice(vL2_L_cross_ij, 2, j) = M_vL2 * slice(crossij_Pot_T, 1, j).transpose();
        }
        // Calculate the Lie algebra operator: take the kronecker product and contract ((v_i)_I, e_i x e_j)_L2 <e_I, e_L> on L-K with C^K_IJ
        for (size_t L=0; L<dim_la; L++){
          L2v_epsBkt_list[iv] += Eigen::KroneckerProduct(slice(vL2_L_cross_ij, 0, L), C_IJ[L]);
        }
      }// for iv
    }// if

    return L2v_epsBkt_list;
  } else if (m_nl_par == 1){

    size_t num_vec = v_L2prod_list.size();

    size_t dim_la = m_liealg.dimension();
    size_t dim_xcurl_T = m_sxcurl.dimensionCell(iT);
    size_t dim_laxcurl_T = m_lasxcurl.dimensionCell(iT);

    // List of cell operators representing (v, *[.,.])_L2 on T
    std::vector<Eigen::MatrixXd> L2v_epsBkt_list(num_vec, Eigen::MatrixXd::Zero(dim_laxcurl_T, dim_laxcurl_T));

    const std::vector<Eigen::MatrixXd> & C_IJ = m_liealg.structureConst();

    size_t dim_P2k3_T = m_laddrcore.P2k3(iT).dimension();

    for (size_t iv = 0; iv < num_vec; iv++){

      // Read the operator (v, .)_L2 as a dim_la by dim_P2k3 matrix [L][i]
      const Eigen::Map<const Eigen::MatrixXd> M_vL2(v_L2prod_list[iv].data(), dim_la, dim_P2k3_T);
      // Multiarray vL2_L_cross_ij[L][i][j] representing ((v_i)_I, e_i x e_j)_L2 <e_I, e_L>
      Scalar3Tensor vL2_L_cross_ij(boost::extents[dim_la][dim_xcurl_T][dim_xcurl_T]);
      std::fill_n(vL2_L_cross_ij.data(), vL2_L_cross_ij.size(), 0.);

      for (size_t j=0; j<dim_xcurl_T; j++){
        slice(vL2_L_cross_ij, 2, j) = M_vL2 * slice(crossij_Pot_T, 1, j).transpose();
      }

      // Calculate the Lie algebra operator: take the Kronecker product and contract ((v_i)_I, e_i x e_j)_L2 <e_I, e_L> on L-K with C^K_IJ
      for (size_t L=0; L<dim_la; L++){
        L2v_epsBkt_list[iv] += Eigen::KroneckerProduct(slice(vL2_L_cross_ij, 0, L), C_IJ[L]);
      }
    }// for iv

    return L2v_epsBkt_list;
  } else{
    std::cerr << "[main] ERROR: Unknown bracket discretisation" << std::endl;
    exit(1);
  }
}

//------------------------------------------------------------------------------

std::vector<YangMillsNorms> YangMills::computeYangMillsNorms(const std::vector<Eigen::VectorXd> & list_dofs) const
{
  const size_t ncells = m_ddrcore.mesh().n_cells();
  const size_t nb_vectors = list_dofs.size();
  std::vector<Eigen::VectorXd> local_sqnorm_E(nb_vectors, Eigen::VectorXd::Zero(ncells));
  std::vector<Eigen::VectorXd> local_sqnorm_A(nb_vectors, Eigen::VectorXd::Zero(ncells));
  std::vector<Eigen::VectorXd> local_sqnorm_lambda(nb_vectors, Eigen::VectorXd::Zero(ncells));

  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [&, this](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *m_ddrcore.mesh().cell(iT);
        // Matrices required to compute the local contributions to the norms
        MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*m_ddrcore.degree() + 2);
        Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_2kp2);
        Eigen::MatrixXd mass_Pkpo_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polykpo, int_mono_2kp2);
        Eigen::MatrixXd L2curl_T = m_lasxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T);
        Eigen::MatrixXd L2grad_T = m_lasxgrad.computeL2Product(iT, m_stab_par, mass_Pkpo_T);
       
        for (size_t i=0; i<nb_vectors; i++){
          // Local DOFs for the LAXcurl space
          Eigen::VectorXd E_curl_T = m_lasxcurl.restrict(T, list_dofs[i].head(m_lasxcurl.dimension()));
          Eigen::VectorXd lambda_grad_T = m_lasxgrad.restrict(T, list_dofs[i].segment(m_lasxcurl.dimension(), m_lasxgrad.dimension()));
          Eigen::VectorXd A_curl_T = m_lasxcurl.restrict(T, list_dofs[i].tail(m_lasxcurl.dimension()));

          // Contribution of L2 norms
          local_sqnorm_E[i](iT) = E_curl_T.transpose() * L2curl_T * E_curl_T;
          local_sqnorm_lambda[i](iT) = lambda_grad_T.transpose() * L2grad_T * lambda_grad_T;
          local_sqnorm_A[i](iT) = A_curl_T.transpose() * L2curl_T * A_curl_T;
        }
      }
    };
  parallel_for(ncells, compute_local_squarednorms, m_use_threads);

  // Assemble the output
  std::vector<YangMillsNorms> list_norms;
  list_norms.reserve(nb_vectors);
  for (size_t i=0; i<nb_vectors; i++){
    double sqnorm_E = local_sqnorm_E[i].sum();
    double sqnorm_A = local_sqnorm_A[i].sum();
    double sqnorm_lambda = local_sqnorm_lambda[i].sum();
    // Eigen::VectorXd E_curl = list_dofs[i].head(m_lasxcurl.dimension()).transpose();
    // Eigen::VectorXd lambda_grad = list_dofs[i].segment(m_lasxcurl.dimension(), m_lasxgrad.dimension());
    // Eigen::VectorXd A_curl = list_dofs[i].tail(m_lasxcurl.dimension()).transpose();
    // double sqnorm_E = E_curl.transpose() * m_laxcurl_L2 * E_curl;
    // double sqnorm_A = A_curl.transpose() * m_laxcurl_L2 * A_curl;
    // double sqnorm_lambda = lambda_grad.transpose() * m_laxgrad_L2 * lambda_grad;
    list_norms[i] = YangMillsNorms(std::sqrt(std::abs(sqnorm_E)), std::sqrt(std::abs(sqnorm_A)), std::sqrt(std::abs(sqnorm_lambda)));
  }  
  return list_norms;
}

//------------------------------------------------------------------------------

double YangMills::computeResidual(const Eigen::VectorXd & Elambda) const
{
  Eigen::VectorXd err = m_A * Elambda - m_b;
  // l2 norms
  double norm_err = err.norm();
  double norm_b = m_b.norm();

  double total_res = (norm_b < 1e-10 ? norm_err : norm_err/double(norm_b));

  return total_res;
}

// //------------------------------------------------------------------------------

double YangMills::computeConstraintNorm(const Eigen::VectorXd & Ch, const size_t itersolver) const
{
  Eigen::VectorXd v(m_lasxgrad.dimension());
  if (itersolver) {
    Eigen::LeastSquaresConjugateGradient<YangMills::SystemMatrixType> solver;
    solver.compute(m_laxgrad_L2);
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
      exit(1);
    }
    v = solver.solve(Ch);
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not solve direct system" << std::endl;
      exit(1);
    }
  } else { 
  #ifdef WITH_MKL   
      Eigen::PardisoLU<YangMills::SystemMatrixType> solver;
  #elif WITH_UMFPACK   
      Eigen::UmfPackLU<YangMills::SystemMatrixType> solver;
  #else 
      Eigen::SparseLU<YangMills::SystemMatrixType> solver;
  #endif
      solver.compute(m_laxgrad_L2);
      if (solver.info() != Eigen::Success) {
        std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
      }
      v = solver.solve(Ch);
      if (solver.info() != Eigen::Success) {
        std::cerr << "[main] ERROR: Could not solve linear system" << std::endl;
      }
  }
  return std::sqrt(v.transpose()*Ch);
}

// //------------------------------------------------------------------------------

Eigen::VectorXd YangMills::computeConstraint(const Eigen::VectorXd & Eh, const Eigen::VectorXd & Ah, const double nonlinear_coeff)
{

std::function<Eigen::VectorXd(size_t iT)> compute_local_const
  = [&, this](size_t iT)-> Eigen::VectorXd
  {
    Cell & T = *m_ddrcore.mesh().cell(iT);
    Eigen::VectorXd E_T = m_lasxcurl.restrict(T, Eh);

    MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_ddrcore.degree());
    Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_2k);
    Eigen::MatrixXd L2laxCurl_Gleft = m_lasxcurl.computeL2ProductGradient(iT, m_sxgrad, "left", m_stab_par, mass_Pk3_T);
    
    Scalar3Tensor int_PciPcjPgk = _integral_Pot_ijk(iT);
    Eigen::MatrixXd L2bkt_Ah = nonlinear_coeff * L2v_Bkt(iT, int_PciPcjPgk, {m_lasxcurl.restrictCell(iT, Ah)}, 2)[0];
    Eigen::VectorXd const_T = L2laxCurl_Gleft * E_T + L2bkt_Ah.transpose() * E_T;
  return const_T;
  };

  std::function<void(size_t, const Eigen::VectorXd &, std::vector<std::list<Eigen::Triplet<double> > > &, std::vector<Eigen::VectorXd> &)> assemble_local_vec
  = [this](size_t iT, const Eigen::VectorXd & v, std::vector<std::list<Eigen::Triplet<double>>> & triplets_sys, std::vector<Eigen::VectorXd> & vecs)-> void
  {
    Cell & T = *m_ddrcore.mesh().cell(iT);
    std::vector<size_t> IT_laxgrad = m_lasxgrad.globalDOFIndices(T);

    for (size_t i = 0; i < IT_laxgrad.size(); i++){
      vecs[0][IT_laxgrad[i]] += v(i);
    }
  };

  auto assemble = [&](size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs) -> void
  {
    for (size_t iT = start; iT < end; iT++){
      assemble_local_vec(iT, compute_local_const(iT), *triplets, *vecs);
    }
  };

  Eigen::VectorXd global_constraint = parallel_assembly_system(m_ddrcore.mesh().n_cells(), {}, {m_lasxgrad.dimension()}, assemble, m_use_threads).vectors[0];

  return global_constraint;
}

// //------------------------------------------------------------------------------

Eigen::VectorXd YangMills::computeInitialConditions(const Eigen::MatrixXd & Eh, const Eigen::MatrixXd & Ah, const double nonlinear_coeff, const size_t itersolver)
{
  std::function<SystemVectors<Eigen::MatrixXd>(size_t iT)> compute_local
  = [&, this](size_t iT)-> SystemVectors<Eigen::MatrixXd>
  {
    Cell & T = *m_ddrcore.mesh().cell(iT);

    Eigen::VectorXd Eh_T = m_lasxcurl.restrictCell(iT, Eh);
    MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*m_ddrcore.degree()+2);
    Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_2kp2);
    // Local product matrices
    Eigen::MatrixXd L2laxcurl = m_lasxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T);
    Eigen::MatrixXd L2laxcurl_Grad_right = m_lasxcurl.computeL2ProductGradient(iT, m_sxgrad, "right", m_stab_par, mass_Pk3_T);
    Scalar3Tensor int_PciPcjPgk = _integral_Pot_ijk(iT);
    Eigen::MatrixXd L2bkt_Ah = nonlinear_coeff * L2v_Bkt(iT, int_PciPcjPgk, {m_lasxcurl.restrictCell(iT, Ah)}, 2)[0];

    size_t dim_laxcurl_T = m_lasxcurl.dimensionCell(iT);
    size_t dim_laxgrad_T = m_lasxgrad.dimensionCell(iT);
    size_t sizeSys_T = dim_laxcurl_T + dim_laxgrad_T;

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(sizeSys_T, sizeSys_T);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(sizeSys_T);

    Eigen::MatrixXd C = L2laxcurl_Grad_right + L2bkt_Ah;
    M.topLeftCorner(dim_laxcurl_T, dim_laxcurl_T) = L2laxcurl;
    M.block(0, dim_laxcurl_T, dim_laxcurl_T, dim_laxgrad_T) = C;
    M.block(dim_laxcurl_T, 0, dim_laxgrad_T, dim_laxcurl_T) = C.transpose();
    b.head(dim_laxcurl_T) = L2laxcurl * Eh_T;

    return SystemVectors<Eigen::MatrixXd>{{M},{b}};
  };

  std::function<void(size_t, const SystemVectors<Eigen::MatrixXd> &, std::vector<std::list<Eigen::Triplet<double>>> &, std::vector<Eigen::VectorXd> &)> assemble_local
  = [this](size_t iT, const SystemVectors<Eigen::MatrixXd> & lsT, std::vector<std::list<Eigen::Triplet<double>>> & triplets_sys, std::vector<Eigen::VectorXd> & vecs)-> void
  {
    auto [systems, vectors] = lsT;
    Cell & T = *m_ddrcore.mesh().cell(iT);
    
    size_t dim_laxcurl_T = m_lasxcurl.dimensionCell(iT);
    size_t dim_laxcurl = m_lasxcurl.dimension();

    std::vector<size_t> IT_laxcurl = m_lasxcurl.globalDOFIndices(T);
    std::vector<size_t> IT_laxgrad = m_lasxgrad.globalDOFIndices(T);

    for (size_t i = 0; i < IT_laxcurl.size(); i++){
      for (size_t j = 0; j < IT_laxcurl.size(); j++){
        triplets_sys[0].emplace_back(IT_laxcurl[i], IT_laxcurl[j], systems[0](i, j));
      }
      vecs[0](IT_laxcurl[i]) += vectors[0](i);
    }

    for (size_t i = 0; i < IT_laxcurl.size(); i++){
      for (size_t j = 0; j < IT_laxgrad.size(); j++){
        triplets_sys[0].emplace_back(IT_laxcurl[i], dim_laxcurl + IT_laxgrad[j], systems[0](i, dim_laxcurl_T + j));
        triplets_sys[0].emplace_back(dim_laxcurl + IT_laxgrad[j], IT_laxcurl[i], systems[0](dim_laxcurl_T + j, i));
      }
    }
  };

  auto assemble = [&](size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs) -> void
  {
    for (size_t iT = start; iT < end; iT++){
      assemble_local(iT, compute_local(iT), *triplets, *vecs);
    }
  };

  size_t sizeSys = m_lasxcurl.dimension() + m_lasxgrad.dimension();
  auto [systems, vectors] = parallel_assembly_system(m_ddrcore.mesh().n_cells(), {std::make_pair(sizeSys, sizeSys)}, {sizeSys}, assemble, m_use_threads);
  YangMills::SystemMatrixType M = systems[0];
  Eigen::VectorXd b = vectors[0];

  Eigen::VectorXd Elambda(m_lasxcurl.dimension()+m_lasxdiv.dimension());
  if (itersolver) {
    Eigen::LeastSquaresConjugateGradient<YangMills::SystemMatrixType> solver;
    // solver.preconditioner().setFillfactor(2);
    solver.compute(M);
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
      exit(1);
    }
    Elambda = solver.solve(b);
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not solve direct system" << std::endl;
      exit(1);
    }
  } else { 
  #ifdef WITH_MKL   
      Eigen::PardisoLU<YangMills::SystemMatrixType> solver;
  #elif WITH_UMFPACK   
      Eigen::UmfPackLU<YangMills::SystemMatrixType> solver;
  #else 
      Eigen::SparseLU<YangMills::SystemMatrixType> solver;
  #endif
      solver.compute(M);
      if (solver.info() != Eigen::Success) {
        std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
      }
      Elambda = solver.solve(b);
      if (solver.info() != Eigen::Success) {
        std::cerr << "[main] ERROR: Could not solve linear system" << std::endl;
      }
  }
  Eigen::VectorXd E = Elambda.head(m_lasxcurl.dimension());
  return E;
}

// // //------------------------------------------------------------------------------
