// Authors: 
#include <fstream>
#include <iomanip>
#include <thread>

#include "laddr-yangmills-lm.hpp"
#include <parallel_for.hpp>
#include <max_degrees_quadratures.hpp>
#include <GMpoly_cell.hpp>
#include <unsupported/Eigen/KroneckerProduct>

#include <boost/program_options.hpp>
#include <display_timer.hpp>

#ifdef WITH_UMFPACK
  #include <Eigen/UmfPackSupport>
#endif

#ifdef WITH_MKL
  #include <Eigen/PardisoSupport>
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
    ("degree,k", boost::program_options::value<size_t>()->default_value(0), "The polynomial degree of the sequence")
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
    ("exact-initial-conditions,q", boost::program_options::value<bool>()->default_value(false), "Use constrained initial conditions");

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
  timer.stop();
  double t_wall_ddrcore, t_proc_ddrcore;
  std::tie(t_wall_ddrcore, t_proc_ddrcore) = store_times(timer, "[main] Time DDRCore (wall/proc) ");

  // Assemble the problem
  timer.start();
  LieAlgebra su2;
  YangMills ym(ddr_core, su2, use_threads);  
  ym.assembleLinearSystem(dt);
  timer.stop();
  double t_wall_model, t_proc_model;
  std::tie(t_wall_model, t_proc_model) = store_times(timer, "[main] Time model (wall/proc) ");
  timer.start();
  // Initial conditions
  Eigen::VectorXd Elambdah = Eigen::VectorXd::Zero(ym.sizeSystem());
  Eigen::VectorXd Ah = ym.laXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ElectricFieldType>(A, t_initial), MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
  Eigen::VectorXd Eh = ym.laXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ElectricFieldType>(E, t_initial), MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
  Elambdah.head(ym.laXCurl().dimension()) = Eh;
  // Set the constrained initial conditions
  if (initial_cond){
    Elambdah.head(ym.laXCurl().dimension()) = ym.computeInitialConditions(Eh, Ah, nonlinear_coeff, itersolver);
  }

  Eigen::VectorXd initial_constraint = ym.computeConstraint(Elambdah.head(ym.laXCurl().dimension()), Ah, nonlinear_coeff);
  double initial_constraint_norm = ym.computeConstraintNorm(initial_constraint, itersolver);
  std::cout << FORMAT(25) << "[main] Initial constraint norm: " << initial_constraint_norm << std::endl;
  double max_res = 0;
  double max_const = 0;
  size_t nonlinear_it = 0;
  // Print solver
  if (itersolver) {
    std::cout << "[main] Solving the linear system using LeastSquaresConjugateGradient" << std::endl;
  } else { 
  #ifdef WITH_MKL
      std::cout << "[main] Solving the linear system using Pardiso" << std::endl; 
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
  for (size_t i = 0 ; i < it; i++){
    if (100*i/it >= prog){
      std::cout << "[main] " << 100*i/it << "%..." << std::endl;
      prog = 100*i/it + 10;
    }
    // Set the time
    double t = t_initial + dt*i;

    // Interpolate of forcing terms in LAXCurl: default at 0, but we put the ones corresponding to the exact solution if solution>0
    Eigen::VectorXd interp_f = Eigen::VectorXd::Zero(ym.laXCurl().dimension());
    Eigen::VectorXd interp_Eipo = Eigen::VectorXd::Zero(ym.laXCurl().dimension());;
    Eigen::VectorXd interp_Ei = Eigen::VectorXd::Zero(ym.laXCurl().dimension());;
    Eigen::VectorXd interp_Aipomtheta = Eigen::VectorXd::Zero(ym.laXCurl().dimension());
    if (solution>0){
      interp_f = ym.laXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ForcingTermType>(f, t+dt), MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
      interp_Eipo = ym.laXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ForcingTermType>(E, t+dt), MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
      interp_Ei = ym.laXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ForcingTermType>(E, t), MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
      interp_Aipomtheta = ym.laXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ForcingTermType>(A, t+(1-theta)*dt), MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
    }


    // System vector of the nonlinear problem
    Eigen::VectorXd Eh_old = Elambdah.head(ym.laXCurl().dimension());
    ym.setSystemVector(interp_f, interp_Eipo-interp_Ei, interp_Aipomtheta, Eh_old, Ah, dt, theta, nonlinear_coeff);
    // Add boundary conditions for chosen solution, only if solution>0
    if (solution>0){
      ym.addBoundaryConditions(ym.contractParaLA<Eigen::Vector3d, YangMills::MagneticFieldType>(B, t+dt), dt);
    }
    Eigen::VectorXd sysVec = ym.systemVector();
    
    // Solution vector of Newton iterations
    Eigen::VectorXd dElambdak;

    while (ym.stoppingCrit2(Elambdah, Eh_old, Ah, sysVec, dt, theta, nonlinear_coeff) > stop_val){
      nonlinear_it++;
      ym.assembleSystemNewton(Eh_old, Ah, Elambdah, sysVec, dt, theta, nonlinear_coeff);
     
      if (itersolver){
        Eigen::LeastSquaresConjugateGradient<YangMills::SystemMatrixType> solver;
        // solver.preconditioner().setFillfactor(2);
        solver.compute(ym.systemMatrix());
        if (solver.info() != Eigen::Success) {
          std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
          exit(1);
        }
        dElambdak = solver.solve(ym.systemVector());
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
            dElambdak = solver.solve(ym.systemVector());
            if (solver.info() != Eigen::Success) {
              std::cerr << "[main] ERROR: Could not solve linear system" << std::endl;
            }
        }
      Elambdah = Elambdah + dElambdak;

      // Update residual and condition numbers
      double res = ym.computeResidual(dElambdak);
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
    };

    // Update A and the constraint violation
    Eigen::VectorXd Eh_new = Elambdah.head(ym.laXCurl().dimension());
    Ah -= dt*(1-theta)*Eh_old + dt*theta*Eh_new;
    max_const = std::max(max_const, ym.computeConstraintNorm(ym.computeConstraint(Eh_new, Ah, nonlinear_coeff)-initial_constraint, itersolver));
  } // for i
}
  cn << std::flush;
  cn.close();
  timer.stop();
  double t_wall_solve, t_proc_solve;
  std::tie(t_wall_solve, t_proc_solve) = store_times(timer, "[main] Time solve (wall/proc) ");

  // Export matrix if requested  
  if (vm.count("export-matrix")) {
    std::cout << "[main] Exporting matrix to Matrix Market format" << std::endl;
    saveMarket(ym.systemMatrix(), "A_maxwell.mtx");
    // saveMarket(ym.systemVector(), "B_maxwell.mtx");
  }

  // Interpolate of exact solution and error vector
  Eigen::VectorXd ElambdaAI = Eigen::VectorXd::Zero(ym.sizeSystem() + ym.laXCurl().dimension());
  Eigen::VectorXd ElambdaAh = Eigen::VectorXd::Zero(ym.sizeSystem() + ym.laXCurl().dimension());
  ElambdaAI.head(ym.laXCurl().dimension()) = ym.laXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ElectricFieldType>(E, t_final), MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
  ElambdaAI.tail(ym.laXCurl().dimension()) = ym.laXCurl().interpolate(ym.contractParaLA<Eigen::Vector3d, YangMills::ElectricFieldType>(A, t_final), MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
  ElambdaAh.head(ym.sizeSystem()) = Elambdah;
  ElambdaAh.tail(ym.laXCurl().dimension()) = Ah;
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
  out << "StabilisationParameter: " << ym.stabilizationParameter() << std::endl;
  out << "Mesh: " << mesh_file << std::endl;
  out << "Degree: " << K << std::endl;
  out << "MeshSize: " << mesh_ptr->h_max() << std::endl;
  out << "NbCells: " << mesh_ptr->n_cells() << std::endl;
  out << "NbFaces: " << mesh_ptr->n_faces() << std::endl;
  out << "NbEdges: " << mesh_ptr->n_edges() << std::endl;
  out << "DimXGrad: " << ym.xGrad().dimension() << std::endl;
  out << "DimXCurl: " << ym.xCurl().dimension() << std::endl;
  out << "DimXDiv: " << ym.xDiv().dimension() << std::endl;
  out << "DimLieAlg: " << ym.lieAlg().dimension() << std::endl;
  out << "DimLAXGrad: " << ym.laXGrad().dimension() << std::endl;
  out << "DimLAXCurl: " << ym.laXCurl().dimension() << std::endl;
  out << "DimLAXDiv: " << ym.laXDiv().dimension() << std::endl;
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
  out << "TwallSolve: " << t_wall_solve << std::endl;  
  out << "TprocSolve: " << t_proc_solve << std::endl; 
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
               const LieAlgebra & liealgebra,
               bool use_threads,
               std::ostream & output
               )
  : m_ddrcore(ddrcore),
    m_use_threads(use_threads),
    m_output(output),
    m_xgrad(ddrcore, use_threads),
    m_xcurl(ddrcore, use_threads),
    m_xdiv(ddrcore, use_threads),
    m_liealg(liealgebra),
    m_laxgrad(liealgebra, m_xgrad, use_threads),
    m_laxcurl(liealgebra, m_xcurl, use_threads),
    m_laxdiv(liealgebra, m_xdiv, use_threads),
    m_A(sizeSystem(), sizeSystem()),
    m_b(Eigen::VectorXd::Zero(sizeSystem())),
    m_laxcurl_L2(m_laxcurl.dimension(), m_laxcurl.dimension()),
    m_laxdiv_L2(m_laxdiv.dimension(), m_laxdiv.dimension()),
    m_laxgrad_L2(m_laxgrad.dimension(), m_laxgrad.dimension()),
    m_F_ebkt(_epsBkt()),
    m_stab_par(0.1)   
{
  // Compute the basis integrals of the bracket (LaxGrad)x(LaxCurl)->LaxCurl
  for (size_t iT = 0; iT < ddrcore.mesh().n_cells(); iT++){
    m_int_PciPcjPgk.emplace_back(boost::multi_array<double, 3>(boost::extents[m_laxcurl.dimensionCell(iT)][m_laxcurl.dimensionCell(iT)][m_laxgrad.dimensionCell(iT)]));
  }
  std::function<void(size_t, size_t)> fill_ints
  = [this](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++) {
        m_int_PciPcjPgk[iT]=_int_Pci_bktPcjPgk(iT);
      } // for iF
    };
  parallel_for(ddrcore.mesh().n_cells(), fill_ints, use_threads);
  
  m_output << "[YangMills] Initializing" << std::endl;
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
  
  std::vector<std::pair<size_t, size_t>> size_systems{std::make_pair(m_laxcurl.dimension(), m_laxcurl.dimension()), std::make_pair(m_laxdiv.dimension(), m_laxdiv.dimension()), std::make_pair(m_laxgrad.dimension(), m_laxgrad.dimension())}; 
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
      Eigen::MatrixXd int_Bi_cross_nF(m_liealg.dimension(), m_xcurl.dimensionFace(F));
      QuadratureRule quad_dqrbc_F = generate_quadrature_rule(F, dqrbc);
      
      for (size_t i = 0; i < m_liealg.dimension(); i++){
        MagneticFieldType Bi = B[i];
        FType<Eigen::Vector3d> Bi_cross_nF = [&Bi, &nF](const Eigen::Vector3d & x) {
                                                        return Bi(x).cross(nF);
                                                        };
        int_Bi_cross_nF.row(i) = integrate(Bi_cross_nF, evaluate_quad<Function>::compute(*m_ddrcore.faceBases(F.global_index()).Polyk2, quad_dqrbc_F), quad_dqrbc_F).transpose() 
        * m_xcurl.faceOperators(F.global_index()).potential;
      }// for i
      Eigen::MatrixXd bF_M = m_liealg.massMatrix() * int_Bi_cross_nF;
      Eigen::Map<Eigen::VectorXd> bF(bF_M.data(), bF_M.size());
      auto I_F = m_laxcurl.globalDOFIndices(F);
      for (size_t i = 0; i < I_F.size(); i++) {
        m_b(I_F[i]) -= dt*bF(i);
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

double YangMills::stoppingCrit2(
                                const Eigen::VectorXd & Elambda_k,
                                const Eigen::VectorXd & E_i,
                                const Eigen::VectorXd & A_i,
                                const Eigen::VectorXd & sysVec,
                                double dt,
                                double theta,
                                double nonlinear_coeff
                                )
{

std::function<Eigen::VectorXd(size_t iT, double dt)> compute_local_F
  = [&, this](size_t iT, double dt)-> Eigen::VectorXd
  {
  Cell & T = *m_ddrcore.mesh().cell(iT);
  // Mass matrix for (P^k(T))^3  
  MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_ddrcore.degree());
  Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_2k);

  // k is the inner loop index, i is the outer loop index (constant in k)
  // Vectors
  Eigen::VectorXd E_k = Elambda_k.head(m_laxcurl.dimension());
  Eigen::VectorXd lambda_k = Elambda_k.segment(m_laxcurl.dimension(), m_laxgrad.dimension());
  Eigen::VectorXd X = A_i - dt * (1-theta) * E_i;
  Eigen::VectorXd Y = A_i - 0.5 * dt * (1-theta) * E_i;
  Eigen::VectorXd Z = A_i - dt * (1-theta) * (1-theta) * E_i;
  // Local vectors
  Eigen::VectorXd E_k_T = m_laxcurl.restrict(T, E_k);
  Eigen::VectorXd lambda_k_T = m_laxgrad.restrict(T, lambda_k);
  Eigen::VectorXd A_i_T = m_laxcurl.restrict(T, A_i);
  Eigen::VectorXd X_T = m_laxcurl.restrict(T, X);
  Eigen::VectorXd Y_T = m_laxcurl.restrict(T, Y);
  // Bracket matrices: Hstar[u,.]
  Eigen::MatrixXd ebkt_A = nonlinear_coeff * epsBkt_v(iT, A_i);
  Eigen::MatrixXd ebkt_E = nonlinear_coeff * epsBkt_v(iT, E_k);
  Eigen::MatrixXd ebkt_X = nonlinear_coeff * epsBkt_v(iT, X);
  Eigen::MatrixXd ebkt_Y = nonlinear_coeff * epsBkt_v(iT, Y);
  // Bracket vectors: Hstar[u,v]
  Eigen::VectorXd ebkt_A_A = ebkt_A * A_i_T;
  Eigen::VectorXd ebkt_E_E = ebkt_E * E_k_T;
  Eigen::VectorXd ebkt_E_A = ebkt_E * A_i_T;
  Eigen::VectorXd ebkt_X_X = ebkt_X * X_T;
  Eigen::VectorXd ebkt_X_E = ebkt_X * E_k_T;
  // L2 product matrices
  Eigen::MatrixXd laxcurl_L2 = m_laxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd laxdiv_L2Curl_left = m_laxdiv.computeL2ProductCurl(iT, m_xcurl, "left", m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd laxdiv_L2Curl_both = m_laxdiv.computeL2ProductCurl(iT, m_xcurl, "both", m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd laxdiv_L2 = m_laxdiv.computeL2Product(iT, m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd laxcurl_L2Grad_right = m_laxcurl.computeL2ProductGradient(iT, m_xgrad, "right", m_stab_par, mass_Pk3_T);
  // L2 product matrices with brackets
  Eigen::MatrixXd L2CurlA_ebkt = nonlinear_coeff * L2v_epsBkt(iT, A_i_T, laxdiv_L2Curl_left);
  Eigen::MatrixXd L2CurlX_ebkt = nonlinear_coeff * L2v_epsBkt(iT, X_T, laxdiv_L2Curl_left);
  Eigen::MatrixXd L2CurlE_ebkt = nonlinear_coeff * L2v_epsBkt(iT, E_k_T, laxdiv_L2Curl_left);
  Eigen::MatrixXd L2ebktAA_ebkt = L2v_epsBkt(iT, ebkt_A_A, laxdiv_L2);
  Eigen::MatrixXd L2ebktXX_ebkt = L2v_epsBkt(iT, ebkt_X_X, laxdiv_L2);
  Eigen::MatrixXd L2ebktXE_ebkt = L2v_epsBkt(iT, ebkt_X_E, laxdiv_L2);
  Eigen::MatrixXd L2ebktEE_ebkt = L2v_epsBkt(iT, ebkt_E_E, laxdiv_L2);
  Eigen::MatrixXd L2Curl_left_ebkt_E = laxdiv_L2Curl_left * ebkt_E;
  Eigen::MatrixXd L2Curl_left_ebkt_X = laxdiv_L2Curl_left * ebkt_X;
  Eigen::MatrixXd L2Curl_left_ebkt_Y = laxdiv_L2Curl_left * ebkt_Y;
  Eigen::MatrixXd L2_ebkt_E = ebkt_E.transpose() * laxdiv_L2;
  // L2 integral bracket product matrices: int <1,[2,3]>
  Eigen::MatrixXd L2E_bkt = nonlinear_coeff * L2v_Bkt(iT, m_int_PciPcjPgk[iT], E_k, 1);
  Eigen::MatrixXd L2Ei_bkt = nonlinear_coeff * L2v_Bkt(iT, m_int_PciPcjPgk[iT], E_i, 1);
  Eigen::MatrixXd L2bkt_Z = nonlinear_coeff * L2v_Bkt(iT, m_int_PciPcjPgk[iT], Z, 2);
  Eigen::MatrixXd L2bkt_E = nonlinear_coeff * L2v_Bkt(iT, m_int_PciPcjPgk[iT], E_k, 2);

  //------------------------------------------------------------------------------
  // Vector
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
        + dt * theta * ebkt_Y.transpose() * laxdiv_L2 * ebkt_X
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
      + (1-theta) *
      (
        - dt * theta * L2bkt_E.transpose()
        + dt * theta * L2Ei_bkt.transpose()
      );

    size_t sizeSys_T = m_laxcurl.dimensionCell(iT) + m_laxgrad.dimensionCell(iT);
    size_t sizeSys1 = m_laxcurl.dimensionCell(iT);
    size_t sizeSys2 = m_laxgrad.dimensionCell(iT);
    Eigen::VectorXd b(sizeSys_T);
    b.head(sizeSys1) = A * E_k_T + B * lambda_k_T;
    b.segment(sizeSys1, sizeSys2) = C * E_k_T;

    return b;
  };
  // Assemble local vectors
  std::function<void(size_t, const Eigen::VectorXd &, std::vector<std::list<Eigen::Triplet<double> > > &, std::vector<Eigen::VectorXd> &)> assemble_local_vec
  = [this](size_t iT, const Eigen::VectorXd & v, std::vector<std::list<Eigen::Triplet<double>>> & triplets_sys, std::vector<Eigen::VectorXd> & vecs)-> void
  {
    Cell & T = *m_ddrcore.mesh().cell(iT);
    std::vector<size_t> IT_laxcurl = m_laxcurl.globalDOFIndices(T);
    std::vector<size_t> IT_laxgrad = m_laxgrad.globalDOFIndices(T);
    
    size_t dimT_laxcurl = m_laxcurl.dimensionCell(iT);
    size_t dim_laxcurl = m_laxcurl.dimension();

    for (size_t i = 0; i < IT_laxcurl.size(); i++){
      vecs[0][IT_laxcurl[i]] += v(i);
    }
    for (size_t i = 0; i < IT_laxgrad.size(); i++){
      vecs[0][dim_laxcurl + IT_laxgrad[i]] += v(dimT_laxcurl + i);
    }
  };
  // Assemble
  auto assemble = [&](size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs) -> void
  {
    for (size_t iT = start; iT < end; iT++){
      assemble_local_vec(iT, compute_local_F(iT, dt), *triplets, *vecs);
    }
  };

  Eigen::VectorXd F_Ek = parallel_assembly_system(m_ddrcore.mesh().n_cells(), {}, {sizeSystem()}, assemble, m_use_threads).vectors[0];
  // std::cout << FORMAT(25) << "[main] Stop crit: " << (F_Ek - sysVec).norm()/sysVec.norm() << std::endl;
  return sysVec.norm() < 1e-13 ? (F_Ek - sysVec).norm() :  (F_Ek - sysVec).norm()/sysVec.norm();
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
 
  std::function<Eigen::VectorXd(size_t iT, double dt)> compute_local_vec
  = [&, this](size_t iT, double dt)-> Eigen::VectorXd
  {
    Cell & T = *m_ddrcore.mesh().cell(iT);
    // Mass matrix for (P^k(T))^3  
    MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_ddrcore.degree());
    Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_2k);
    // Vectors
    Eigen::VectorXd X = A_old - dt * (1-theta) * E_old;
    Eigen::VectorXd Y = A_old - 0.5 * dt * (1-theta) * E_old;
    Eigen::VectorXd Z = A_old - dt * (1-theta) * (1-theta) * E_old;
    // Local vectors
    Eigen::VectorXd A_old_T = m_laxcurl.restrict(T, A_old);
    Eigen::VectorXd E_old_T = m_laxcurl.restrict(T, E_old);
    Eigen::VectorXd interp_f_T = m_laxcurl.restrict(T, interp_f);
    Eigen::VectorXd interp_dE_T = m_laxcurl.restrict(T, interp_dE);
    Eigen::VectorXd X_T = m_laxcurl.restrict(T, X);
    // Bracket matrices: Hstar[u,.]
    Eigen::MatrixXd ebkt_A = nonlinear_coeff * epsBkt_v(iT, A_old);
    Eigen::MatrixXd ebkt_X = nonlinear_coeff * epsBkt_v(iT, X);
    Eigen::MatrixXd ebkt_Y = nonlinear_coeff * epsBkt_v(iT, Y);
    // Bracket vectors: Hstar[u,v]
    Eigen::VectorXd ebkt_A_A = ebkt_A * A_old_T;
    Eigen::VectorXd ebkt_X_X = ebkt_X * X_T;
    // L2 product matrices
    Eigen::MatrixXd laxdiv_L2 = m_laxdiv.computeL2Product(iT, m_stab_par, mass_Pk3_T);
    Eigen::MatrixXd laxcurl_L2 = m_laxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T);
    Eigen::MatrixXd laxdiv_L2Curl_left = m_laxdiv.computeL2ProductCurl(iT, m_xcurl, "left", m_stab_par, mass_Pk3_T);
    Eigen::MatrixXd laxdiv_L2Curl_both = m_laxdiv.computeL2ProductCurl(iT, m_xcurl, "both", m_stab_par, mass_Pk3_T);
    Eigen::MatrixXd laxcurl_L2Grad_left = m_laxcurl.computeL2ProductGradient(iT, m_xgrad, "left", m_stab_par, mass_Pk3_T);
    // L2 product matrices with brackets
    Eigen::MatrixXd L2_ebkt_Y = ebkt_Y.transpose() * laxdiv_L2;
    Eigen::MatrixXd L2Curl_left_ebkt_Y = laxdiv_L2Curl_left * ebkt_Y;
    // L2 integral bracket product matrices: int <1,[2,3]> 
    Eigen::MatrixXd L2bkt_interp_A = nonlinear_coeff * L2v_Bkt(iT, m_int_PciPcjPgk[iT], interp_A, 2);
    Eigen::MatrixXd L2bkt_Z = nonlinear_coeff * L2v_Bkt(iT, m_int_PciPcjPgk[iT], Z, 2);

    //------------------------------------------------------------------------------
    // System Vector 
    //------------------------------------------------------------------------------

    size_t dimT_laxcurl = m_laxcurl.dimensionCell(iT);
    size_t dimT_laxgrad = m_laxgrad.dimensionCell(iT);

    Eigen::VectorXd b = Eigen::VectorXd::Zero(dimT_laxcurl + dimT_laxgrad + m_liealg.dimension());

    Eigen::VectorXd b1
      = laxcurl_L2 * E_old_T
      + dt * (1-theta) * 
      (
        laxdiv_L2Curl_both * A_old_T
        + L2Curl_left_ebkt_Y.transpose() * A_old_T
        + 0.5 * laxdiv_L2Curl_left * ebkt_A_A
        + 0.5 * L2_ebkt_Y * ebkt_A_A
      )
      + dt * theta *
      (
        laxdiv_L2Curl_both * X_T
        + L2Curl_left_ebkt_Y.transpose() * X_T
        + 0.5 * laxdiv_L2Curl_left * ebkt_X_X
        + 0.5 * L2_ebkt_Y * ebkt_X_X
      )
      + dt * laxcurl_L2 * interp_f_T;
    
    Eigen::VectorXd b2
      =  laxcurl_L2Grad_left * E_old_T
        + L2bkt_Z.transpose() * E_old_T
        + laxcurl_L2Grad_left * interp_dE_T
        + L2bkt_interp_A.transpose() * interp_dE_T;

    b.head(dimT_laxcurl) = b1;
    b.segment(dimT_laxcurl, dimT_laxgrad) = b2;

    return b;
  };
  // Assemble local vectors
  std::function<void(size_t, const Eigen::VectorXd &, std::vector<std::list<Eigen::Triplet<double> > > &, std::vector<Eigen::VectorXd> &)> assemble_local_vec
  = [this](size_t iT, const Eigen::VectorXd & v, std::vector<std::list<Eigen::Triplet<double>>> & triplets_sys, std::vector<Eigen::VectorXd> & vecs)-> void
  {
    Cell & T = *m_ddrcore.mesh().cell(iT);
    std::vector<size_t> IT_laxcurl = m_laxcurl.globalDOFIndices(T);
    std::vector<size_t> IT_laxgrad = m_laxgrad.globalDOFIndices(T);
    
    size_t dimT_laxcurl = m_laxcurl.dimensionCell(iT);
    size_t dim_laxcurl = m_laxcurl.dimension();

    for (size_t i = 0; i < IT_laxcurl.size(); i++){
      vecs[0][IT_laxcurl[i]] += v(i);
    }
    for (size_t i = 0; i < IT_laxgrad.size(); i++){
      vecs[0][dim_laxcurl + IT_laxgrad[i]] += v(dimT_laxcurl + i);
    }
  };
  // Assemble vector
  auto assemble = [&](size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs) -> void
  {
    for (size_t iT = start; iT < end; iT++){
      assemble_local_vec(iT, compute_local_vec(iT, dt), *triplets, *vecs);
    }
  };
  m_b = parallel_assembly_system(m_ddrcore.mesh().n_cells(), {}, {sizeSystem()}, assemble, m_use_threads).vectors[0];
}

//------------------------------------------------------------------------------

void YangMills::assembleSystemNewton(
                                    const Eigen::VectorXd & E_i,    
                                    const Eigen::VectorXd & A_i,    
                                    const Eigen::VectorXd & Elambda_k,
                                    const Eigen::VectorXd & sysVec,
                                    double dt,                       
                                    double theta,
                                    double nonlinear_coeff
                                    )
{ 
  std::function<AssembleSystems(size_t iT, double dt)> compute_local_Newton
  = [&, this](size_t iT, double dt)-> AssembleSystems
  {
  Cell & T = *m_ddrcore.mesh().cell(iT);
  // Mass matrix for (P^k(T))^3  
  MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_ddrcore.degree());
  Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_2k);

  // k is the inner loop index, i is the outer loop index (constant in k)
  // Vectors
  Eigen::VectorXd E_k = Elambda_k.head(m_laxcurl.dimension());
  Eigen::VectorXd lambda_k = Elambda_k.segment(m_laxcurl.dimension(), m_laxgrad.dimension());
  Eigen::VectorXd X = A_i - dt * (1-theta) * E_i;
  Eigen::VectorXd Y = A_i - 0.5 * dt * (1-theta) * E_i;
  Eigen::VectorXd Z = A_i - dt * (1-theta) * (1-theta) * E_i;
  // Local vectors
  Eigen::VectorXd E_k_T = m_laxcurl.restrict(T, E_k);
  Eigen::VectorXd lambda_k_T = m_laxgrad.restrict(T, lambda_k);
  Eigen::VectorXd A_i_T = m_laxcurl.restrict(T, A_i);
  Eigen::VectorXd X_T = m_laxcurl.restrict(T, X);
  Eigen::VectorXd Y_T = m_laxcurl.restrict(T, Y);
  // Bracket matrices: Hstar[u,.]
  Eigen::MatrixXd ebkt_A = nonlinear_coeff * epsBkt_v(iT, A_i);
  Eigen::MatrixXd ebkt_E = nonlinear_coeff * epsBkt_v(iT, E_k);
  Eigen::MatrixXd ebkt_X = nonlinear_coeff * epsBkt_v(iT, X);
  Eigen::MatrixXd ebkt_Y = nonlinear_coeff * epsBkt_v(iT, Y);
  // Bracket vectors: Hstar[u,v]
  Eigen::VectorXd ebkt_A_A = ebkt_A * A_i_T;
  Eigen::VectorXd ebkt_E_E = ebkt_E * E_k_T;
  Eigen::VectorXd ebkt_E_A = ebkt_E * A_i_T;
  Eigen::VectorXd ebkt_X_X = ebkt_X * X_T;
  Eigen::VectorXd ebkt_X_E = ebkt_X * E_k_T;
  // L2 product matrices
  Eigen::MatrixXd laxcurl_L2 = m_laxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd laxdiv_L2Curl_left = m_laxdiv.computeL2ProductCurl(iT, m_xcurl, "left", m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd laxdiv_L2Curl_both = m_laxdiv.computeL2ProductCurl(iT, m_xcurl, "both", m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd laxdiv_L2 = m_laxdiv.computeL2Product(iT, m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd laxcurl_L2Grad_right = m_laxcurl.computeL2ProductGradient(iT, m_xgrad, "right", m_stab_par, mass_Pk3_T);
  // L2 product matrices with brackets
  Eigen::MatrixXd L2CurlA_ebkt = nonlinear_coeff * L2v_epsBkt(iT, A_i_T, laxdiv_L2Curl_left);
  Eigen::MatrixXd L2CurlX_ebkt = nonlinear_coeff * L2v_epsBkt(iT, X_T, laxdiv_L2Curl_left);
  Eigen::MatrixXd L2CurlE_ebkt = nonlinear_coeff * L2v_epsBkt(iT, E_k_T, laxdiv_L2Curl_left);
  Eigen::MatrixXd L2ebktAA_ebkt = L2v_epsBkt(iT, ebkt_A_A, laxdiv_L2);
  Eigen::MatrixXd L2ebktXX_ebkt = L2v_epsBkt(iT, ebkt_X_X, laxdiv_L2);
  Eigen::MatrixXd L2ebktXE_ebkt = L2v_epsBkt(iT, ebkt_X_E, laxdiv_L2);
  Eigen::MatrixXd L2ebktEE_ebkt = L2v_epsBkt(iT, ebkt_E_E, laxdiv_L2);
  Eigen::MatrixXd L2Curl_left_ebkt_E = laxdiv_L2Curl_left * ebkt_E;
  Eigen::MatrixXd L2Curl_left_ebkt_X = laxdiv_L2Curl_left * ebkt_X;
  Eigen::MatrixXd L2Curl_left_ebkt_Y = laxdiv_L2Curl_left * ebkt_Y;
  Eigen::MatrixXd L2_ebkt_E = ebkt_E.transpose() * laxdiv_L2;
  // L2 integral bracket product matrices: int <1,[2,3]> 
  Eigen::MatrixXd L2E_bkt = nonlinear_coeff * L2v_Bkt(iT, m_int_PciPcjPgk[iT], E_k, 1);
  Eigen::MatrixXd L2Ei_bkt = nonlinear_coeff * L2v_Bkt(iT, m_int_PciPcjPgk[iT], E_i, 1);
  Eigen::MatrixXd L2bkt_Z = nonlinear_coeff * L2v_Bkt(iT, m_int_PciPcjPgk[iT], Z, 2);
  Eigen::MatrixXd L2bkt_E = nonlinear_coeff * L2v_Bkt(iT, m_int_PciPcjPgk[iT], E_k, 2);
  Eigen::MatrixXd L2bkt_lambda = nonlinear_coeff * L2v_Bkt(iT, m_int_PciPcjPgk[iT], lambda_k, 3);

  //------------------------------------------------------------------------------
  // System Vector 
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
        + dt * theta * ebkt_Y.transpose() * laxdiv_L2 * ebkt_X
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

  size_t sizeSys_T = m_laxcurl.dimensionCell(iT) + m_laxgrad.dimensionCell(iT);
  size_t sizeSys1 = m_laxcurl.dimensionCell(iT);
  size_t sizeSys2 = m_laxgrad.dimensionCell(iT);

  // Nonlinear problem evalulated at kth iteration
  Eigen::VectorXd b(sizeSys_T);
  b.head(sizeSys1) = A * E_k_T + B * lambda_k_T;
  b.segment(sizeSys1, sizeSys2) = C * E_k_T;

  //------------------------------------------------------------------------------
  // System Matrix
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

  Eigen::MatrixXd sysM = Eigen::MatrixXd::Zero(sizeSys_T, sizeSys_T);
  sysM.topLeftCorner(sizeSys1, sizeSys1) = A;
  sysM.block(0, sizeSys1, sizeSys1, sizeSys2) = B;
  sysM.block(sizeSys1, 0, sizeSys2, sizeSys1) = C;

  return AssembleSystems{{sysM}, {b}};
  };
  // Assemble local contributions
  std::function<void(size_t, const AssembleSystems &, std::vector<std::list<Eigen::Triplet<double>>> &, std::vector<Eigen::VectorXd> &)> assemble_local_Newton
  = [this](size_t iT, const AssembleSystems & lsT, std::vector<std::list<Eigen::Triplet<double>>> & triplets_sys, std::vector<Eigen::VectorXd> & vecs)-> void
  {
    auto [systems, vectors] = lsT;
    Cell & T = *m_ddrcore.mesh().cell(iT);
    
    size_t dimT_laxcurl = m_laxcurl.dimensionCell(iT);
    size_t dim_laxcurl = m_laxcurl.dimension();

    std::vector<size_t> IT_laxcurl = m_laxcurl.globalDOFIndices(T);
    std::vector<size_t> IT_laxgrad = m_laxgrad.globalDOFIndices(T);

    for (size_t i = 0; i < IT_laxcurl.size(); i++){
      for (size_t j = 0; j < IT_laxcurl.size(); j++){
        triplets_sys[0].emplace_back(IT_laxcurl[i], IT_laxcurl[j], systems[0](i, j));
      }
      vecs[0](IT_laxcurl[i]) += vectors[0](i);
    }

    for (size_t i = 0; i < IT_laxcurl.size(); i++){
      for (size_t j = 0; j < IT_laxgrad.size(); j++){
        triplets_sys[0].emplace_back(IT_laxcurl[i], dim_laxcurl + IT_laxgrad[j], systems[0](i, dimT_laxcurl + j));
        triplets_sys[0].emplace_back(dim_laxcurl + IT_laxgrad[j], IT_laxcurl[i], systems[0](dimT_laxcurl + j, i));
      }
    }

    for (size_t i = 0; i < IT_laxgrad.size(); i++){
      vecs[0](dim_laxcurl + IT_laxgrad[i]) += vectors[0](dimT_laxcurl + i);
    }
  };

  auto assemble = [&](size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs) -> void
  {
    for (size_t iT = start; iT < end; iT++){
      assemble_local_Newton(iT, compute_local_Newton(iT, dt), *triplets, *vecs);
    }
  };
  // Assemble the matrix and set the rhs
  auto [systems, vectors] = parallel_assembly_system(m_ddrcore.mesh().n_cells(), {std::make_pair(sizeSystem(), sizeSystem())}, {sizeSystem()}, assemble, m_use_threads);
  m_A = systems[0];
  m_b = -vectors[0] + sysVec;
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
  MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_ddrcore.degree());
  Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_2k);
  Eigen::MatrixXd mass_Pkpo_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polykpo, int_mono_2kp2);

  Eigen::MatrixXd L2laxcurl = m_laxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd L2laxdiv = m_laxdiv.computeL2Product(iT, m_stab_par, mass_Pk3_T);
  Eigen::MatrixXd L2laxgrad = m_laxgrad.computeL2Product(iT, m_stab_par, mass_Pkpo_T);

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

  std::vector<size_t> IT_sysg = m_laxgrad.globalDOFIndices(T);
  std::vector<size_t> IT_sysc = m_laxcurl.globalDOFIndices(T);
  std::vector<size_t> IT_sysd = m_laxdiv.globalDOFIndices(T);

  //Global L2 product matrices
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

std::vector<Eigen::MatrixXd> YangMills::_epsBkt()
{
  std::vector<Eigen::MatrixXd> epsBkt(m_laxdiv.dimension());
  std::vector<Eigen::MatrixXd> C_IJ = m_liealg.structureConst();

  std::function<void(size_t, size_t)> compute_ebkt
  = [&, this](size_t start, size_t end)->void
  {
    for (size_t iF = start; iF < end; iF++){
      Face & F = *m_laxcurl.mesh().face(iF);
      // Xcurl matrix to recover scalar face values
      Eigen::MatrixXd detProduct = _compute_det_bracket_face(iF);
      std::vector<size_t> IF_laXdiv = m_laxdiv.globalDOFIndices(F);
      // Kronecker product with structure constants to expand in Lie algebra indices
      for (size_t K = 0; K < m_liealg.dimension(); K++){
        epsBkt[IF_laXdiv[K]] = 1./m_xcurl.faceBases(iF).Polyk->function(0, F.center_mass()) * Eigen::KroneckerProduct(detProduct, C_IJ[K]);
      }
    }
  };
  parallel_for(m_ddrcore.mesh().n_faces(), compute_ebkt, m_use_threads);
  return epsBkt;
}

//------------------------------------------------------------------------------

Eigen::MatrixXd YangMills::_compute_det_bracket_face(const size_t iF) const
{
  Face & F = *m_xcurl.mesh().face(iF);

  Eigen::MatrixXd potM_F = _compute_face_tangent(iF);

  Eigen::MatrixXd det_nij = Eigen::MatrixXd::Zero(m_xcurl.dimensionFace(iF), m_xcurl.dimensionFace(iF));
  Eigen::MatrixXd nij = Eigen::MatrixXd::Zero(F.normal().size(), F.normal().size());
  nij.col(0) = F.normal();
  for (size_t i = 0; i < m_xcurl.dimensionFace(iF); i++){
    nij.col(1) = potM_F.col(i);
    for (size_t j = 0; j < i; j++){
      nij.col(2) = potM_F.col(j);
      det_nij(i, j) = nij.determinant();
      det_nij(j, i) = -det_nij(i, j);
    }
  }
  return det_nij;
}

//------------------------------------------------------------------------------

Eigen::MatrixXd YangMills::_compute_face_tangent(const size_t iF) const
{
  Face & F = *m_xcurl.mesh().face(iF);
  Eigen::MatrixXd pot_F = m_xcurl.faceOperators(iF).potential;

  size_t dim_F = pot_F.rows();
  size_t dim_func = 3;

  Eigen::MatrixXd basis_F = Eigen::MatrixXd::Zero(dim_func, dim_F);

  for (size_t i = 0; i < dim_F; i++){
    basis_F.col(i) = (*m_xcurl.faceBases(iF).Polyk2).function(i, F.center_mass());
  }

  return basis_F * pot_F;
}

//------------------------------------------------------------------------------

Eigen::MatrixXd YangMills::epsBkt_v(size_t iT, const Eigen::VectorXd & v) const
{
  Cell & T = *m_ddrcore.mesh().cell(iT);
  Eigen::MatrixXd epsBktv_T = Eigen::MatrixXd::Zero(m_laxdiv.dimensionCell(iT), m_laxcurl.dimensionCell(iT));
    
  std::vector<size_t> IT_laXdiv = m_laxdiv.globalDOFIndices(T);
  for (size_t iF = 0; iF < T.n_faces(); iF++){
    Face & F = *T.face(iF);
    size_t DOFIndexF = iF * m_liealg.dimension();
    // Extends the face brackets calculated in _epsBkt to T
    for (size_t K = 0; K < m_liealg.dimension(); K++){
      Eigen::RowVectorXd F_vepsBkt = m_F_ebkt[IT_laXdiv[DOFIndexF + K]] * m_laxcurl.restrict(F, v);
      epsBktv_T.row(DOFIndexF + K) = m_laxcurl.extendOperator(T, F, F_vepsBkt);
    }
  }
  return epsBktv_T;
}

//------------------------------------------------------------------------------

boost::multi_array<double, 3> YangMills::_integral_ijk(size_t iT) const
{
  Cell & T = *m_ddrcore.mesh().cell(iT);
  
  size_t dim_xcurl_T = m_xcurl.dimensionCell(iT);
  size_t dim_xgrad_T = m_xgrad.dimensionCell(iT);

  size_t dim_Pk3 = m_ddrcore.cellBases(iT).Polyk3->dimension();
  size_t dim_Pkpo = m_ddrcore.cellBases(iT).Polykpo->dimension();

  // integral of the basis functions in T (int (\phi_i . \phi_j) * \psi_k)
  boost::multi_array<double, 3> integral_basis(boost::extents[dim_Pk3][dim_Pk3][dim_Pkpo]);
  std::fill_n(integral_basis.data(), integral_basis.num_elements(), 0);
  // Calculate quadrature rule for 3k+1 degrees
  QuadratureRule quad_3kpo_T = generate_quadrature_rule(T, 3*m_ddrcore.degree()+1);
  boost::multi_array<Eigen::Vector3d, 2> B_Pk3 = evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polyk3, quad_3kpo_T);
  boost::multi_array<double, 2> B_Pkpo = evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polykpo, quad_3kpo_T);
  // Manual integration using quadrature rules
  for (size_t i = 0; i < dim_Pk3; i++){
    for (size_t j = 0; j <= i; j++){
      for (size_t k = 0; k < dim_Pk3; k++){
        for (size_t iqn=0; iqn<quad_3kpo_T.size(); iqn++){
          integral_basis[i][j][k] += quad_3kpo_T[iqn].w * B_Pk3[i][iqn].dot(B_Pk3[j][iqn]) * B_Pkpo[k][iqn];
        }
        integral_basis[j][i][k] = integral_basis[i][j][k];
      }
    }
  }

  // integral over the basis in T of (Xcurl)x(Xcurl)x(xGrad) potentials
  boost::multi_array<double, 3> integral_PcPcPg(boost::extents[dim_xcurl_T][dim_xcurl_T][dim_xgrad_T]);
  std::fill_n(integral_PcPcPg.data(), integral_PcPcPg.num_elements(), 0);
  Eigen::MatrixXd xgrad_pot = m_xgrad.cellOperators(iT).potential;
  Eigen::MatrixXd xcurl_pot = m_xcurl.cellOperators(iT).potential;

  for (size_t i = 0; i < dim_xcurl_T; i++){
    Eigen::MatrixXd integral_jk = Eigen::MatrixXd::Zero(dim_Pk3, dim_Pkpo);
    for (size_t l=0; l<dim_Pk3; l++){
      Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>  int_ljk(integral_basis.data() + l * dim_Pk3 * dim_Pkpo, dim_Pk3, dim_Pkpo, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(1, dim_Pkpo));
      integral_jk += xcurl_pot(l, i) * int_ljk;
    }
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

boost::multi_array<double, 3> YangMills::_int_Pci_bktPcjPgk(size_t iT)
{ 
  size_t dim_xcurl_T = m_xcurl.dimensionCell(iT);
  size_t dim_xgrad_T = m_xgrad.dimensionCell(iT);
  size_t dim_laxcurl_T = m_laxcurl.dimensionCell(iT);
  size_t dim_laxgrad_T = m_laxgrad.dimensionCell(iT);

  boost::multi_array<double, 3> intPciPcjPgk(boost::extents[dim_laxcurl_T][dim_laxcurl_T][dim_laxgrad_T]);
  std::fill_n(intPciPcjPgk.data(), intPciPcjPgk.num_elements(), 0);
  std::vector<Eigen::MatrixXd> C_IJ = m_liealg.structureConst();
  Eigen::MatrixXd massMatrix = m_liealg.massMatrix();

  // Combine the integral on the triple basis integrals with the Lie bracket
  boost::multi_array<double, 3> int_ijk = _integral_ijk(iT);
  for (size_t i = 0; i < dim_xcurl_T; i++){
    Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>  int_jk(int_ijk.data() + i * dim_xcurl_T * dim_xgrad_T, dim_xcurl_T, dim_xgrad_T, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(1, dim_xgrad_T));
    for (size_t K = 0; K < m_liealg.dimension(); K++){
      Eigen::MatrixXd int_ijk_cJK = Eigen::KroneckerProduct(int_jk, C_IJ[K]);
      for (size_t L = 0; L < m_liealg.dimension(); L++){
        Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(intPciPcjPgk.data() + (i*m_liealg.dimension()+L) * dim_laxcurl_T * dim_laxgrad_T, dim_laxcurl_T, dim_laxgrad_T, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(1, dim_laxgrad_T)) += massMatrix(L,K)*int_ijk_cJK;
      }
    }
  }
  return intPciPcjPgk;
}

//------------------------------------------------------------------------------

Eigen::MatrixXd YangMills::L2v_Bkt(size_t iT, boost::multi_array<double, 3> & intPciPcjPgk, const Eigen::VectorXd & v, const size_t & entry) const
{
  size_t dim_laxcurl_T = m_laxcurl.dimensionCell(iT);
  size_t dim_laxgrad_T = m_laxgrad.dimensionCell(iT);
  if (entry == 1){
    Eigen::MatrixXd L2matrix = Eigen::MatrixXd::Zero(dim_laxcurl_T, dim_laxgrad_T);
    Eigen::VectorXd v_T = m_laxcurl.restrictCell(iT, v);
    for (size_t i = 0; i < dim_laxcurl_T; i++){
      L2matrix += v_T[i]*Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(intPciPcjPgk.data() + i * dim_laxcurl_T * dim_laxgrad_T, dim_laxcurl_T, dim_laxgrad_T, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(1, dim_laxgrad_T));
    }
    return L2matrix;
  }else if (entry == 2){
     Eigen::MatrixXd L2matrix = Eigen::MatrixXd::Zero(dim_laxcurl_T, dim_laxgrad_T);
     Eigen::VectorXd v_T = m_laxcurl.restrictCell(iT, v);
    for (size_t i = 0; i < dim_laxcurl_T; i++){
      L2matrix += v_T[i]*Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(intPciPcjPgk.data() + i * dim_laxgrad_T, dim_laxcurl_T, dim_laxgrad_T, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(1, dim_laxcurl_T * dim_laxgrad_T));
    }
    return L2matrix;
  }else if (entry == 3){
     Eigen::MatrixXd L2matrix = Eigen::MatrixXd::Zero(dim_laxcurl_T, dim_laxcurl_T);
     Eigen::VectorXd v_T = m_laxgrad.restrictCell(iT, v);
    for (size_t i = 0; i < dim_laxgrad_T; i++){
      L2matrix += v_T[i]*Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(intPciPcjPgk.data() + i, dim_laxcurl_T, dim_laxcurl_T, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(dim_laxgrad_T, dim_laxcurl_T * dim_laxgrad_T));
    }
    return L2matrix;
  }else{
    std::cout << "[Main] Invalid entry in L2v_Bkt" << std::endl;
    exit(1);
  }
}

//------------------------------------------------------------------------------

Eigen::MatrixXd YangMills::L2v_epsBkt(size_t iT, const Eigen::VectorXd & v_T, const Eigen::MatrixXd & L2prod) const
{
    Cell & T = *m_ddrcore.mesh().cell(iT);

    Eigen::VectorXd T_L2prodv = v_T.transpose() * L2prod;;

    std::vector<size_t> IT_laXdiv = m_laxdiv.globalDOFIndices(T);
    Eigen::MatrixXd T_L2vepsBkt = Eigen::MatrixXd::Zero(m_laxcurl.dimensionCell(T), m_laxcurl.dimensionCell(T));

    for (size_t iF = 0; iF < T.n_faces(); iF++){
      Face & F = *T.face(iF);
      size_t DOFIndexF = iF * m_liealg.dimension();
      for (size_t K = 0; K < m_liealg.dimension(); K++){
        T_L2vepsBkt += T_L2prodv[DOFIndexF + K] * m_laxcurl.extendOperator(T, F, m_laxcurl.extendOperator(T, F, m_F_ebkt[IT_laXdiv[DOFIndexF + K]]).transpose());
      }
    }
    return T_L2vepsBkt;
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
        Eigen::MatrixXd L2curl_T = m_laxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T);
        Eigen::MatrixXd L2grad_T = m_laxgrad.computeL2Product(iT, m_stab_par, mass_Pkpo_T);
       
        for (size_t i=0; i<nb_vectors; i++){
          // Local DOFs for the LAXcurl space
          Eigen::VectorXd E_curl_T = m_laxcurl.restrict(T, list_dofs[i].head(m_laxcurl.dimension()));
          Eigen::VectorXd lambda_grad_T = m_laxgrad.restrict(T, list_dofs[i].segment(m_laxcurl.dimension(), m_laxgrad.dimension()));
          Eigen::VectorXd A_curl_T = m_laxcurl.restrict(T, list_dofs[i].tail(m_laxcurl.dimension()));

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
    list_norms[i] = YangMillsNorms(std::sqrt(std::abs(sqnorm_E)), std::sqrt(std::abs(sqnorm_A)), std::sqrt(std::abs(sqnorm_lambda)));
  }  
  return list_norms;
}

//------------------------------------------------------------------------------

double YangMills::computeResidual(const Eigen::VectorXd & Elambda) const
{
  Eigen::VectorXd err = m_A * Elambda - m_b;
  Eigen::VectorXd err_E = err.head(m_laxcurl.dimension());
  Eigen::VectorXd err_lambda = err.segment(m_laxcurl.dimension(), m_laxgrad.dimension());
  // l2 norms
  double norm_err_E = err_E.norm();
  double norm_err_lambda = err_lambda.norm();
  double norm_b_E = m_b.head(m_laxcurl.dimension()).norm();

  double total_res = (norm_b_E ? norm_err_E/double(norm_b_E) : norm_err_E) + norm_err_lambda;

  return total_res;
}

// //------------------------------------------------------------------------------

double YangMills::computeConstraintNorm(const Eigen::VectorXd & Ch, const size_t itersolver) const
{
  Eigen::VectorXd v(m_laxgrad.dimension());
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
    Eigen::VectorXd E_T = m_laxcurl.restrict(T, Eh);

    MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_ddrcore.degree());
    Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_2k);
    Eigen::MatrixXd L2laxCurl_Gleft = m_laxcurl.computeL2ProductGradient(iT, m_xgrad, "left", m_stab_par, mass_Pk3_T);

    Eigen::MatrixXd L2bkt_Ah = nonlinear_coeff * L2v_Bkt(iT, m_int_PciPcjPgk[iT], Ah, 2);
    Eigen::VectorXd const_T = L2laxCurl_Gleft * E_T + L2bkt_Ah.transpose() * E_T;
  return const_T;
  };

  std::function<void(size_t, const Eigen::VectorXd &, std::vector<std::list<Eigen::Triplet<double> > > &, std::vector<Eigen::VectorXd> &)> assemble_local_vec
  = [this](size_t iT, const Eigen::VectorXd & v, std::vector<std::list<Eigen::Triplet<double>>> & triplets_sys, std::vector<Eigen::VectorXd> & vecs)-> void
  {
    Cell & T = *m_ddrcore.mesh().cell(iT);
    std::vector<size_t> IT_laxgrad = m_laxgrad.globalDOFIndices(T);

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

  Eigen::VectorXd global_constraint = parallel_assembly_system(m_ddrcore.mesh().n_cells(), {}, {m_laxgrad.dimension()}, assemble, m_use_threads).vectors[0];

  return global_constraint;
}

// //------------------------------------------------------------------------------

Eigen::VectorXd YangMills::computeInitialConditions(const Eigen::MatrixXd & Eh, const Eigen::MatrixXd & Ah, const double nonlinear_coeff, const size_t itersolver)
{
  std::function<AssembleSystems(size_t iT)> compute_local
  = [&, this](size_t iT)-> AssembleSystems
  {
    Cell & T = *m_ddrcore.mesh().cell(iT);

    Eigen::VectorXd Eh_T = m_laxcurl.restrictCell(iT, Eh);
    MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*m_ddrcore.degree()+2);
    Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_2kp2);
    // Local product matrices
    Eigen::MatrixXd L2laxcurl = m_laxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T);
    Eigen::MatrixXd L2laxcurl_Grad_right = m_laxcurl.computeL2ProductGradient(iT, m_xgrad, "right", m_stab_par, mass_Pk3_T);
    Eigen::MatrixXd L2bkt_Ah = nonlinear_coeff * L2v_Bkt(iT, m_int_PciPcjPgk[iT], Ah, 2);

    size_t dimT_laxcurl = m_laxcurl.dimensionCell(iT);
    size_t dimT_laxgrad = m_laxgrad.dimensionCell(iT);
    size_t sizeSys_T = dimT_laxcurl + dimT_laxgrad;

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(sizeSys_T, sizeSys_T);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(sizeSys_T);

    Eigen::MatrixXd C = L2laxcurl_Grad_right + L2bkt_Ah;
    M.topLeftCorner(dimT_laxcurl, dimT_laxcurl) = L2laxcurl;
    M.block(0, dimT_laxcurl, dimT_laxcurl, dimT_laxgrad) = C;
    M.block(dimT_laxcurl, 0, dimT_laxgrad, dimT_laxcurl) = C.transpose();
    b.head(dimT_laxcurl) = L2laxcurl * Eh_T;

    return AssembleSystems{{M},{b}};
  };

  std::function<void(size_t, const AssembleSystems &, std::vector<std::list<Eigen::Triplet<double>>> &, std::vector<Eigen::VectorXd> &)> assemble_local
  = [this](size_t iT, const AssembleSystems & lsT, std::vector<std::list<Eigen::Triplet<double>>> & triplets_sys, std::vector<Eigen::VectorXd> & vecs)-> void
  {
    auto [systems, vectors] = lsT;
    Cell & T = *m_ddrcore.mesh().cell(iT);
    
    size_t dimT_laxcurl = m_laxcurl.dimensionCell(iT);
    size_t dim_laxcurl = m_laxcurl.dimension();

    std::vector<size_t> IT_laxcurl = m_laxcurl.globalDOFIndices(T);
    std::vector<size_t> IT_laxgrad = m_laxgrad.globalDOFIndices(T);

    for (size_t i = 0; i < IT_laxcurl.size(); i++){
      for (size_t j = 0; j < IT_laxcurl.size(); j++){
        triplets_sys[0].emplace_back(IT_laxcurl[i], IT_laxcurl[j], systems[0](i, j));
      }
      vecs[0](IT_laxcurl[i]) += vectors[0](i);
    }

    for (size_t i = 0; i < IT_laxcurl.size(); i++){
      for (size_t j = 0; j < IT_laxgrad.size(); j++){
        triplets_sys[0].emplace_back(IT_laxcurl[i], dim_laxcurl + IT_laxgrad[j], systems[0](i, dimT_laxcurl + j));
        triplets_sys[0].emplace_back(dim_laxcurl + IT_laxgrad[j], IT_laxcurl[i], systems[0](dimT_laxcurl + j, i));
      }
    }
  };

  auto assemble = [&](size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs) -> void
  {
    for (size_t iT = start; iT < end; iT++){
      assemble_local(iT, compute_local(iT), *triplets, *vecs);
    }
  };

  size_t sizeSys = m_laxcurl.dimension() + m_laxgrad.dimension();
  auto [systems, vectors] = parallel_assembly_system(m_ddrcore.mesh().n_cells(), {std::make_pair(sizeSys, sizeSys)}, {sizeSys}, assemble, m_use_threads);
  YangMills::SystemMatrixType M = systems[0];
  Eigen::VectorXd b = vectors[0];

  Eigen::VectorXd Elambda(m_laxcurl.dimension()+m_laxdiv.dimension());
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
  Eigen::VectorXd E = Elambda.head(m_laxcurl.dimension());
  return E;
}

// // //------------------------------------------------------------------------------
