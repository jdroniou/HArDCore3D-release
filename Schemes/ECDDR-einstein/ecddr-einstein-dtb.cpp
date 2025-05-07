// Author: Jia Jia Qian (jia.qian@monash.edu)

#include <fstream>
#include <iomanip>
#include <thread>

#include "ecddr-einstein-dtb.hpp"
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

// // Max number of cells to plot a graph
// constexpr size_t max_nb_cells_for_plot = 20000;

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
    ("solution,s", boost::program_options::value<int>()->default_value(0), "Select the solution")
    ("iterations,n", boost::program_options::value<int>()->default_value(-5), "Set the number of iterations")
    ("h-scaling,d", boost::program_options::value<size_t>()->default_value(1), "Set dt=h^{K+?}")
    ("initial-time,t", boost::program_options::value<double>()->default_value(0.5), "Set the initial time")
    ("final-time,f", boost::program_options::value<double>()->default_value(1.), "Set the final time")
    ("implicit,i", boost::program_options::value<bool>()->default_value(true), "Explicit or implicit scheme")
    ("plot", boost::program_options::value<std::string>(), "Save plot of the solution to the given filename")
    ("stabilisation,b", boost::program_options::value<double>()->default_value(0.1), "Choose stabilisation")
    ("solver", boost::program_options::value<std::string>()->default_value("PardisoLU"), "Choice of solver for implicit scheme. Options are: PardisoLU, UMFPACK, PaStiXLU, PaStiXLLT, EigenLU, EigenBiCGSTAB (reverts to EigenLU if the selected solver is not available)");

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

  // Implicit or explicit scheme
  bool is_implicit = (vm.count("implicit") ? vm["implicit"].as<bool>() : false);
  if (is_implicit){
    std::cout << "[main] Semi-implicit scheme " << mesh_file << std::endl;
  }
  else{
    std::cout << "[main] Explicit scheme" << std::endl;
  }

  Einstein::TCollection1Forms theta, E, H;
  Einstein::TCollection2Forms dH, B, D, dtD, dE;
  Einstein::TZeroFormType lapse;

  // Select the solution
  int solution = (vm.count("solution") ? vm["solution"].as<int>() : 0);

  switch (solution) {
  case 0:
    std::cout << "[main] Test case: Minkowski solution" << std::endl;
    lapse = minkowski_lapse;
    theta = minkowski_theta;
    E = minkowski_E;
    B = minkowski_B;
    D = minkowski_D;
    dtD = minkowski_dtD;
    H = minkowski_H;
    dH = minkowski_dH;
    break;
  case 1:
    std::cout << "[main] Test case: Kasner solution" << std::endl;
    std::cout << "[main] Geodesic slicing (lapse = 1)" << std::endl;
    lapse = minkowski_lapse;
    theta = kasner_theta;
    E = kasner_E;
    E[0] = minkowski_E0;
    B = kasner_B;
    D = kasner_D;
    dtD = kasner_dtD;
    H = kasner_H;
    dH = kasner_dH;
    break;
  case 2:
    std::cout << "[main] Test case: Constant solution" << std::endl;
    std::cout << "[main] Geodesic slicing (lapse = 1)" << std::endl;
    lapse = minkowski_lapse;
    theta = constant_theta;
    E = constant_E;
    E[0] = minkowski_E0;
    B = constant_B;
    D = constant_D;
    dtD = constant_dtD;
    H = constant_H;
    dH = constant_dH;
    dE = constant_dE;
    break;
  case 3:
    std::cout << "[main] Test case: Linear solution" << std::endl;
    std::cout << "[main] Geodesic slicing (lapse = 1)" << std::endl;
    lapse = minkowski_lapse;
    theta = linear_theta;
    E = linear_E;
    E[0] = minkowski_E0;
    B = linear_B;
    D = linear_D;
    dtD = linear_dtD;
    H = linear_H;
    dH = linear_dH;
    dE = linear_dE;
    break;
  case 4:
    std::cout << "[main] Test case: Test solution" << std::endl;
    std::cout << "[main] Gowdy lapse" << std::endl;
    lapse = test_lapse;
    theta = test_theta;
    E = test_E;
    B = test_B;
    D = test_D;
    dtD = test_dtD;
    H = test_H;
    dH = test_dH;
    dE = test_dE;
    break;
  case 5:
    std::cout << "[main] Test case: Gowdy solution" << std::endl;
    std::cout << "[main] Gowdy slicing" << std::endl;
    lapse = gowdy_lapse;
    theta = gowdy_theta;
    E = gowdy_E;
    B = gowdy_B;
    D = gowdy_D;
    dtD = gowdy_dtD;
    H = gowdy_H;
    dH = gowdy_dH;
    break;
  default:
    std::cerr << "[main] ERROR: Unknown exact solution" << std::endl;
    exit(1);    
  }

  // Build the mesh
  MeshBuilder meshbuilder = MeshBuilder(mesh_file);
  std::unique_ptr<Mesh> mesh_ptr = meshbuilder.build_the_mesh();
  // Stabilisation parameter
  double stab_par = vm["stabilisation"].as<double>();

  boost::timer::cpu_timer timer;
  // Create ECDDR core
  timer.start();
  bool use_threads = (vm.count("pthread") ? vm["pthread"].as<bool>() : true);
  std::cout << "[main] " << (use_threads ? "Parallel execution" : "Sequential execution") << std:: endl;
  DDR_Spaces ddr_spaces(*mesh_ptr, K, use_threads);
  timer.stop();
  double t_wall_ddrcore, t_proc_ddrcore;
  std::tie(t_wall_ddrcore, t_proc_ddrcore) = store_times(timer, "[main] Time ECDDRCore (wall/proc) ");

  // Create Einstein problem
  Einstein es(ddr_spaces, use_threads, stab_par);  

  // Select linear solver
  std::string name_solver = vm["solver"].as<std::string>();
  LinearSolver<Einstein::SystemMatrixType> solver(name_solver);
  std::cout << "[main] Solving the system using " << solver.name() << std::endl;

  timer.start();

  //--------------------------------------------
  //    Solve problem
  //--------------------------------------------

  // Select the number of iterations and final time
  double t_initial = vm["initial-time"].as<double>();
  double t_final = vm["final-time"].as<double>();
  std::cout << FORMAT(25) << "[main] Initial time: " << t_initial << std::endl;
  std::cout << FORMAT(25) << "[main] Final time: " << t_final << std::endl;

  int iter = vm["iterations"].as<int>();
  size_t h_scaling = vm["h-scaling"].as<size_t>();

  // Set the time step and iterations, positive number of iterations chosen, or if negative, generated to scale in ^K+1 with the mesh size
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
      it = size_t(ceil(std::abs(iter) * std::pow(double(std::cbrt(measure3))/mesh_ptr->h_max(), K+h_scaling)));
    }else{
      it = size_t(ceil(std::pow(double(std::cbrt(measure3))/mesh_ptr->h_max(), K+h_scaling)));
    }
    dt = (t_final - t_initial)/double(it);
  }

  std::cout << FORMAT(25) << "[main] Timestep size: " << dt << std::endl;

  // (Explicit & Implicit) Computes and sets LLT decomposition of product matrices
  es.precomputeL2();

  //--------------------------------------------
  //    Initial conditions
  //--------------------------------------------

  timer.start();

  size_t dim_X1 = es.ddrSpaces().dofspace(1).dimension();
  size_t dim_X2 = es.ddrSpaces().dofspace(2).dimension();

  Eigen::VectorXd starDtB_0(3*(dim_X1 + dim_X2 + dim_X1));
  
  for (size_t i=0; i<3; i++){
    starDtB_0.segment(dim_X1*i, dim_X1) = es.ddrSpaces().interpolate(DDR_Spaces::DDR_function_type{es.contractParaVec<Eigen::Vector3d, Einstein::TwoFormType>(starCol(D), t_initial)[i+1]}, 1);
    starDtB_0.segment(dim_X1*3 + dim_X2*i, dim_X2) = es.ddrSpaces().interpolate(DDR_Spaces::DDR_function_type{es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(starCol(theta), t_initial)[i+1]}, 2);
    starDtB_0.segment(dim_X1*3 + dim_X2*3 + dim_X1*i, dim_X1) = es.ddrSpaces().interpolate(DDR_Spaces::DDR_function_type{es.contractParaVec<Eigen::Vector3d, Einstein::TwoFormType>(starCol(B), t_initial)[i+1]}, 1);
  }
  Eigen::VectorXd starDtB_n(starDtB_0);

  double t_wall_assemble, t_proc_assemble = 0 ;
  double t_wall_solve, t_proc_solve = 0;
  boost::timer::cpu_timer it_timer;

  int n_plots = 0;
  double max_residual = 0;

  // --------------------------------------------------------------------------
  //           Creates initial .vtu file of D, theta, B at initial time
  // --------------------------------------------------------------------------
  // Only if we do not have too many cells, and only for some t
  // Max number of cells to plot a graph
  constexpr size_t max_nb_cells_for_plot = 20000;
  size_t nvertex = mesh_ptr->n_vertices();

  if (vm.count("plot") && mesh_ptr->n_cells() <= max_nb_cells_for_plot) {
    std::cout << "[main] Writing solution " << n_plots << " to file" << std::endl;
    n_plots += 1;
    std::string filename = vm["plot"].as<std::string>();
    VtuWriter plotdata(mesh_ptr.get());
    
    // Exact D, theta, B at the vertices
    std::vector<std::vector<VectorRd>> exact_D_vertex(3, {nvertex, VectorRd::Zero()});
    std::vector<std::vector<VectorRd>> exact_theta_vertex(3, {nvertex, VectorRd::Zero()});
    std::vector<std::vector<VectorRd>> exact_B_vertex(3, {nvertex, VectorRd::Zero()});

    for (size_t i=0; i<3; i++){
      for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
        VectorRd v_coords = mesh_ptr->vertex(iV)->coords();
        exact_D_vertex[i][iV] = starCol(D)[i+1](t_initial, v_coords);
        exact_theta_vertex[i][iV] = starCol(theta)[i+1](t_initial, v_coords);
        exact_B_vertex[i][iV] = starCol(B)[i+1](t_initial, v_coords);
      }
    }

    // D, theta, B at the vertices
    std::vector<std::vector<VectorRd>> D_vertex(3, {nvertex, VectorRd::Zero()});
    std::vector<std::vector<VectorRd>> theta_vertex(3, {nvertex, VectorRd::Zero()});
    std::vector<std::vector<VectorRd>> B_vertex(3, {nvertex, VectorRd::Zero()});

    for (size_t i=0; i<3; i++){
      D_vertex[i] = es.ddrSpaces().computeVertexValues(1, starDtB_n.segment(i*dim_X1, dim_X1));
      theta_vertex[i] = es.ddrSpaces().computeVertexValues(2, starDtB_n.segment(3*dim_X1+dim_X2*i, dim_X2));
      B_vertex[i] = es.ddrSpaces().computeVertexValues(1, starDtB_n.segment(3*dim_X1+3*dim_X2+i*dim_X1, dim_X1));
    }

    // Errors of D, theta, B at the vertices
    std::vector<std::vector<VectorRd>> err_D_vertex(3, {nvertex, VectorRd::Zero()});
    std::vector<std::vector<VectorRd>> err_theta_vertex(3, {nvertex, VectorRd::Zero()});
    std::vector<std::vector<VectorRd>> err_B_vertex(3, {nvertex, VectorRd::Zero()});

    for (size_t i=0; i<3; i++){
      for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
        err_D_vertex[i][iV] = exact_D_vertex[i][iV]- D_vertex[i][iV];
        err_theta_vertex[i][iV] = exact_theta_vertex[i][iV]- theta_vertex[i][iV];
        err_B_vertex[i][iV] = exact_B_vertex[i][iV]- B_vertex[i][iV];
      }
    }

    // Plot
    plotdata.write_to_vtu(filename + "_D" + std::to_string(n_plots) + ".vtu", 
                          std::vector<std::vector<VectorRd>> {D_vertex[0], D_vertex[1], D_vertex[2], exact_D_vertex[0], exact_D_vertex[1], exact_D_vertex[2], err_D_vertex[0], err_D_vertex[1], err_D_vertex[2]},
                          std::vector<std::string>{"D1", "D2", "D3", "exact_D1", "exact_D2", "exact_D3", "err_D1", "err_D2", "err_D3"});
    plotdata.write_to_vtu(filename + "_theta" + std::to_string(n_plots) + ".vtu", 
                          std::vector<std::vector<VectorRd>> {theta_vertex[0], theta_vertex[1], theta_vertex[2], exact_theta_vertex[0], exact_theta_vertex[1], exact_theta_vertex[2], err_theta_vertex[0], err_theta_vertex[1], err_theta_vertex[2]},
                          std::vector<std::string>{"theta1", "theta2", "theta3", "exact_theta1", "exact_theta2", "exact_theta3", "err_theta1", "err_theta2", "err_theta3"});
    plotdata.write_to_vtu(filename + "_B"  + std::to_string(n_plots) + ".vtu", 
                          std::vector<std::vector<VectorRd>> {B_vertex[0], B_vertex[1], B_vertex[2], exact_B_vertex[0], exact_B_vertex[1], exact_B_vertex[2], err_B_vertex[0], err_B_vertex[1], err_B_vertex[2]},
                          std::vector<std::string>{"B1", "B2", "B3", "exact_B1", "exact_B2", "exact_B3", "err_B1", "err_B2", "err_B3"}); 
  }

  //--------------------------------------------
  //    Iterations
  //--------------------------------------------

  for (size_t i = 0 ; i < it; i++){

    // Set the time
    double t = t_initial + dt*i;
    std::cout << "[main] Time: " << t << " ("<< i+1 << "/" << it << ")" << std::endl;

    // Solution at previous time
    Eigen::VectorXd starDtB_nmo(starDtB_n);

    it_timer.start();
    // Computes and sets system matrix and RHS
    es.assembleLinearSystem(
                            dt, 
                            es.contractPara<double, Einstein::TZeroFormType>(lapse, t),
                            es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(E, t)[0],
                            starDtB_n
                            );
    it_timer.stop();
    double it_wall_assemble, it_proc_assemble;
    std::tie(it_wall_assemble, it_proc_assemble) = store_times(it_timer, "[main] Time assemble (wall/proc) ");
    t_wall_assemble += it_wall_assemble;
    t_proc_assemble += it_proc_assemble;

    // Computes and adds forcing terms and boundary conditions to the RHS
    es.computeBCRHS(
                    dt,
                    es.contractPara<double, Einstein::TZeroFormType>(lapse, t),
                    es.contractParaVec<Eigen::Vector3d, Einstein::TwoFormType>(D, t),
                    es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(theta, t),
                    es.contractParaVec<Eigen::Vector3d, Einstein::TwoFormType>(B, t),
                    es.contractParaVec<Eigen::Vector3d, Einstein::TwoFormType>(dtD, t),
                    es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(E, t),
                    es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(H, t),
                    es.contractParaVec<Eigen::Vector3d, Einstein::TwoFormType>(dH, t)
                    );

    /// Solving the linear system
    if (is_implicit){ // Implicit

      it_timer.start();
      /// Solve problem (Eq. 1 & 3) using solver
      Eigen::VectorXd starDB = solver.compute_and_solve(es.systemMatrix(), es.systemVector());
      it_timer.stop();
      double it_wall_solve, it_proc_solve;
      std::tie(it_wall_solve, it_proc_solve) = store_times(it_timer, "[main] Time solve (wall/proc) ");
      t_wall_solve += it_wall_solve;
      t_proc_solve += it_proc_solve;

      double residual_it = solver.residual(es.systemMatrix(), es.systemVector(), starDB);
      if (residual_it > max_residual){
        max_residual = residual_it;
      }
      std::cout << "[main] Residual (via linearsolver) " << residual_it << std::endl;

      starDtB_n.head(3*dim_X1) = starDB.head(3*dim_X1);
      starDtB_n.tail(3*dim_X1) = starDB.tail(3*dim_X1);
      
      /// Update rhs of theta by Eq. 2 using new variable E(n+1)
      es.updateThetaRHS(starDB.head(3*dim_X1));
      
      for(size_t j=0; j<3; j++){
        LinearSolver<Einstein::SystemMatrixType> solver1(name_solver);
        Eigen::VectorXd rhs_theta_i = es.rhsTheta().segment(dim_X2*j, dim_X2);
        starDtB_n.segment(dim_X1*3+dim_X2*j, dim_X2) = solver1.compute_and_solve(es.L2_X2(), rhs_theta_i);
      }

    }
    else { // Explicit

      // Find rhs
      Eigen::VectorXd rhs1 = es.systemVector().head(dim_X1*3) - es.systemMatrix().block(0, dim_X1*3, dim_X1*3, dim_X1*3) * starDtB_n.tail(3*dim_X1);
      Eigen::VectorXd rhs2 = es.rhsTheta() + es.lhsTheta() * starDtB_n.head(3*dim_X1);
      Eigen::VectorXd rhs3 = es.systemVector().tail(dim_X1*3) - es.systemMatrix().block(dim_X1*3, 0, dim_X1*3, dim_X1*3) * starDtB_n.head(3*dim_X1);
      
      it_timer.start();
      for (size_t i=0; i<3; i++){
        Eigen::VectorXd rhs_starD_i = rhs1.segment(i*dim_X1, dim_X1);
        Eigen::VectorXd rhs_theta_i = rhs2.segment(i*dim_X2, dim_X2);
        Eigen::VectorXd rhs_B_i = rhs3.segment(i*dim_X1, dim_X1);

        starDtB_n.segment(dim_X1*i, dim_X1) = solver.compute_and_solve(es.L2_X1(), rhs_starD_i);
        LinearSolver<Einstein::SystemMatrixType> solver1(name_solver);
        starDtB_n.segment(dim_X1*3+dim_X2*i, dim_X2) = solver1.compute_and_solve(es.L2_X2(), rhs_theta_i);
        starDtB_n.segment(dim_X1*3+dim_X2*3+dim_X1*i, dim_X1) = solver.compute_and_solve(es.L2_X1(), rhs_B_i);
      }
      it_timer.stop();
      double it_wall_solve, it_proc_solve;
      std::tie(it_wall_solve, it_proc_solve) = store_times(it_timer, "[main] Time solve (wall/proc) ");
      t_wall_solve += it_wall_solve;
      t_proc_solve += it_proc_solve;
    }

    // --------------------------------------------------------------------------
    //           Computes and prints error and constraint at time t
    // --------------------------------------------------------------------------

    // Eigen::VectorXd starDtB_I = Eigen::VectorXd::Zero((dim_X1+dim_X2+dim_X1)*3);
    // for (size_t i=0; i<3; i++){
    //   starDtB_I.segment(dim_X1*i, dim_X1) = es.ddrSpaces().interpolate(DDR_Spaces::DDR_function_type{es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(starCol(D), t+dt)[i+1]}, 1);
    //   starDtB_I.segment(dim_X1*3+dim_X2*i, dim_X2) = es.ddrSpaces().interpolate(DDR_Spaces::DDR_function_type{es.contractParaVec<Eigen::Vector3d, Einstein::TwoFormType>(starCol(theta), t+dt)[i+1]}, 2);
    //   starDtB_I.segment(dim_X1*3+dim_X2*3+dim_X1*i, dim_X1)= es.ddrSpaces().interpolate(DDR_Spaces::DDR_function_type{es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(starCol(B), t+dt)[i+1]}, 1);
    // }
    // Eigen::VectorXd eDtB = starDtB_I-starDtB_n;
    // std::cout << "[main] Time: " << t << std::endl;
    // std::cout << "[main] Compute errors" << std::endl;

    // // Errors in discrete norms
    // std::cout << "[main] Discrete norms and errors:" << std::endl;
    // std::vector<EinsteinNorms> list_norms = es.computeEinsteinNorms(std::vector<Eigen::VectorXd> {eDtB, starDtB_I, starDtB_n, starDtB_nmo});

    // EinsteinNorms errors = list_norms[0];
    // EinsteinNorms norms = list_norms[1];
    // EinsteinNorms discrete_norms = list_norms[2];

    // std::cout << FORMAT(25) << "\tNorms: *D= "  << FORMAT(15) << norms.D[0]  << FORMAT(15) << norms.D[1]  << FORMAT(15) << norms.D[2] <<std::endl;
    // std::cout << FORMAT(25) << "\tNorms: *Theta= "  << FORMAT(15) << norms.theta[0]  << FORMAT(15) << norms.theta[1]  << FORMAT(15) << norms.theta[2] <<std::endl;
    // std::cout << FORMAT(25) << "\tNorms: *B= "  << FORMAT(15) << norms.B[0]  << FORMAT(15) << norms.B[1]  << FORMAT(15) << norms.B[2] <<std::endl;

    // std::cout << FORMAT(25) << "\tRel. errors: *D= "  << FORMAT(15) << (norms.D[0] < 1e-12 ? errors.D[0] : errors.D[0] / norms.D[0])  << FORMAT(15) << (norms.D[1] < 1e-12 ? errors.D[1] : errors.D[1] / norms.D[1])  << FORMAT(15) << (norms.D[2] < 1e-12 ? errors.D[2] : errors.D[2] / norms.D[2])<< std::endl;
    // std::cout << FORMAT(25) << "\tRel. errors: *Theta= " << FORMAT(15) << (norms.theta[0] < 1e-12 ? errors.theta[0] : errors.theta[0] / norms.theta[0]) << FORMAT(15) << (norms.theta[1] < 1e-12 ? errors.theta[1] : errors.theta[1] / norms.theta[1]) << FORMAT(15) << (norms.theta[2] < 1e-12 ? errors.theta[2] : errors.theta[2] / norms.theta[2])  << FORMAT(15) << std::endl;
    // std::cout << FORMAT(25) << "\tRel. errors: *B= "  << FORMAT(15)<< (norms.B[0] < 1e-12 ? errors.B[0] : errors.B[0] / norms.B[0])  << FORMAT(15) << (norms.B[1] < 1e-12 ? errors.B[1] : errors.B[1] / norms.B[1])  << FORMAT(15) << (norms.B[2] < 1e-12 ? errors.B[2] : errors.B[2] / norms.B[2])<< std::endl;  
  
    // Eigen::VectorXd constraints = es.computeConstraints(es.contractPara<double, Einstein::TZeroFormType>(lapse, t),
    //                                                     es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(theta, t_initial),
    //                                                     es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(theta, t+dt), 
    //                                                     starDtB_0, 
    //                                                     starDtB_n,
    //                                                     solver);

    // std::cout << "Constraint 1: " << FORMAT(15) << constraints[0]  << FORMAT(15) << constraints[1]  << FORMAT(15) << constraints[2] <<std::endl;
    // std::cout << "Constraint 2: " << FORMAT(15) << constraints[3]  << FORMAT(15) << constraints[4]  << FORMAT(15) << constraints[5] <<std::endl;
    // std::cout << "Minimum det P Theta: " << FORMAT(15) << es.min_detTheta() << std::endl;

    // --------------------------------------------------------------------------
    //           Creates a .vtu file of D, theta, B at time t
    // --------------------------------------------------------------------------
    // Only if we do not have too many cells, and only for some t
    // Max number of cells to plot a graph
    constexpr size_t max_nb_cells_for_plot = 20000;
    size_t nvertex = mesh_ptr->n_vertices();

    if (vm.count("plot") && mesh_ptr->n_cells() <= max_nb_cells_for_plot && double(i)/it >= double(n_plots)/100) {
      std::cout << "[main] Writing solution " << n_plots << " to file" << std::endl;
      n_plots += 1;
      std::string filename = vm["plot"].as<std::string>();
      VtuWriter plotdata(mesh_ptr.get());
      
      // Exact D, theta, B at the vertices
      std::vector<std::vector<VectorRd>> exact_D_vertex(3, {nvertex, VectorRd::Zero()});
      std::vector<std::vector<VectorRd>> exact_theta_vertex(3, {nvertex, VectorRd::Zero()});
      std::vector<std::vector<VectorRd>> exact_B_vertex(3, {nvertex, VectorRd::Zero()});

      for (size_t i=0; i<3; i++){
        for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
          VectorRd v_coords = mesh_ptr->vertex(iV)->coords();
          exact_D_vertex[i][iV] = starCol(D)[i+1](t+dt, v_coords);
          exact_theta_vertex[i][iV] = starCol(theta)[i+1](t+dt, v_coords);
          exact_B_vertex[i][iV] = starCol(B)[i+1](t+dt, v_coords);
        }
      }

      // D, theta, B at the vertices
      std::vector<std::vector<VectorRd>> D_vertex(3, {nvertex, VectorRd::Zero()});
      std::vector<std::vector<VectorRd>> theta_vertex(3, {nvertex, VectorRd::Zero()});
      std::vector<std::vector<VectorRd>> B_vertex(3, {nvertex, VectorRd::Zero()});

      for (size_t i=0; i<3; i++){
        D_vertex[i] = es.ddrSpaces().computeVertexValues(1, starDtB_n.segment(i*dim_X1, dim_X1));
        theta_vertex[i] = es.ddrSpaces().computeVertexValues(2, starDtB_n.segment(3*dim_X1+dim_X2*i, dim_X2));
        B_vertex[i] = es.ddrSpaces().computeVertexValues(1, starDtB_n.segment(3*dim_X1+3*dim_X2+i*dim_X1, dim_X1));
      }
          
      // Errors of D, theta, B at the vertices
      std::vector<std::vector<VectorRd>> err_D_vertex(3, {nvertex, VectorRd::Zero()});
      std::vector<std::vector<VectorRd>> err_theta_vertex(3, {nvertex, VectorRd::Zero()});
      std::vector<std::vector<VectorRd>> err_B_vertex(3, {nvertex, VectorRd::Zero()});

      for (size_t i=0; i<3; i++){
        for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
          err_D_vertex[i][iV] = exact_D_vertex[i][iV]- D_vertex[i][iV];
          err_theta_vertex[i][iV] = exact_theta_vertex[i][iV]- theta_vertex[i][iV];
          err_B_vertex[i][iV] = exact_B_vertex[i][iV]- B_vertex[i][iV];
        }
      }

      // Plot
      plotdata.write_to_vtu(filename + "_D" + std::to_string(n_plots) + ".vtu", 
                            std::vector<std::vector<VectorRd>> {D_vertex[0], D_vertex[1], D_vertex[2], exact_D_vertex[0], exact_D_vertex[1], exact_D_vertex[2], err_D_vertex[0], err_D_vertex[1], err_D_vertex[2]},
                            std::vector<std::string>{"D1", "D2", "D3", "exact_D1", "exact_D2", "exact_D3", "err_D1", "err_D2", "err_D3"});
      plotdata.write_to_vtu(filename + "_theta" + std::to_string(n_plots) + ".vtu", 
                            std::vector<std::vector<VectorRd>> {theta_vertex[0], theta_vertex[1], theta_vertex[2], exact_theta_vertex[0], exact_theta_vertex[1], exact_theta_vertex[2], err_theta_vertex[0], err_theta_vertex[1], err_theta_vertex[2]},
                            std::vector<std::string>{"theta1", "theta2", "theta3", "exact_theta1", "exact_theta2", "exact_theta3", "err_theta1", "err_theta2", "err_theta3"});
      plotdata.write_to_vtu(filename + "_B"  + std::to_string(n_plots) + ".vtu", 
                            std::vector<std::vector<VectorRd>> {B_vertex[0], B_vertex[1], B_vertex[2], exact_B_vertex[0], exact_B_vertex[1], exact_B_vertex[2], err_B_vertex[0], err_B_vertex[1], err_B_vertex[2]},
                            std::vector<std::string>{"B1", "B2", "B3", "exact_B1", "exact_B2", "exact_B3", "err_B1", "err_B2", "err_B3"}); 
    }
  }

  timer.stop();
  double t_wall_run, t_proc_run;
  std::tie(t_wall_run, t_proc_run) = store_times(timer, "[main] Time solution (wall/proc) ");
  t_wall_assemble = t_wall_assemble/double(it);
  t_proc_assemble = t_proc_assemble/double(it);
  t_wall_solve = t_wall_solve/double(it);
  t_proc_solve = t_proc_solve/double(it);

  // ----------------------------------------------
  //    OUTPUTS
  // ----------------------------------------------
  
  // Interpolate of exact solution and error vector
  
  Eigen::VectorXd starDtB_I((dim_X1+dim_X2+dim_X1)*3);
  for (size_t i=0; i<3; i++){
    starDtB_I.segment(dim_X1*i, dim_X1) = es.ddrSpaces().interpolate(DDR_Spaces::DDR_function_type{es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(starCol(D), t_final)[i+1]}, 1);
    starDtB_I.segment(dim_X1*3 + dim_X2*i, dim_X2) = es.ddrSpaces().interpolate(DDR_Spaces::DDR_function_type{es.contractParaVec<Eigen::Vector3d, Einstein::TwoFormType>(starCol(theta), t_final)[i+1]}, 2);
    starDtB_I.segment(dim_X1*3 + dim_X2*3 + dim_X1*i, dim_X1) = es.ddrSpaces().interpolate(DDR_Spaces::DDR_function_type{es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(starCol(B), t_final)[i+1]}, 1);
  }
  Eigen::VectorXd eDtB = starDtB_I-starDtB_n;

  std::cout << "[main] Compute errors" << std::endl;
  timer.start();
  // Errors in discrete norms
  std::cout << "[main] Discrete norms and errors:" << std::endl;
  std::vector<EinsteinNorms> list_norms = es.computeEinsteinNorms(std::vector<Eigen::VectorXd> {eDtB, starDtB_I, starDtB_n, starDtB_0});

  EinsteinNorms errors = list_norms[0];
  EinsteinNorms norms = list_norms[1];
  EinsteinNorms discrete_norms = list_norms[2];

  std::cout << FORMAT(25) << "\tNorms: *D= "  << FORMAT(15) << norms.D[0]  << FORMAT(15) << norms.D[1]  << FORMAT(15) << norms.D[2] <<std::endl;
  std::cout << FORMAT(25) << "\tNorms: *Theta= "  << FORMAT(15) << norms.theta[0]  << FORMAT(15) << norms.theta[1]  << FORMAT(15) << norms.theta[2] <<std::endl;
  std::cout << FORMAT(25) << "\tNorms: *B= "  << FORMAT(15) << norms.B[0]  << FORMAT(15) << norms.B[1]  << FORMAT(15) << norms.B[2] <<std::endl;

  std::cout << FORMAT(25) << "\tRel. errors: *D= "  << FORMAT(15) << (norms.D[0] < 1e-12 ? errors.D[0] : errors.D[0] / norms.D[0])  << FORMAT(15) << (norms.D[1] < 1e-12 ? errors.D[1] : errors.D[1] / norms.D[1])  << FORMAT(15) << (norms.D[2] < 1e-12 ? errors.D[2] : errors.D[2] / norms.D[2])<< std::endl;
  std::cout << FORMAT(25) << "\tRel. errors: *Theta= " << FORMAT(15) << (norms.theta[0] < 1e-12 ? errors.theta[0] : errors.theta[0] / norms.theta[0]) << FORMAT(15) << (norms.theta[1] < 1e-12 ? errors.theta[1] : errors.theta[1] / norms.theta[1]) << FORMAT(15) << (norms.theta[2] < 1e-12 ? errors.theta[2] : errors.theta[2] / norms.theta[2])  << FORMAT(15) << std::endl;
  std::cout << FORMAT(25) << "\tRel. errors: *B= "  << FORMAT(15)<< (norms.B[0] < 1e-12 ? errors.B[0] : errors.B[0] / norms.B[0])  << FORMAT(15) << (norms.B[1] < 1e-12 ? errors.B[1] : errors.B[1] / norms.B[1])  << FORMAT(15) << (norms.B[2] < 1e-12 ? errors.B[2] : errors.B[2] / norms.B[2])<< std::endl;

  std::pair<EinsteinNorms, EinsteinNorms> c_errors_norms = es.computeContinuousErrorsNorms(starDtB_n, 
                                  es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(starCol(D), t_final),
                                  es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(starCol(theta), t_final),
                                  es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(starCol(B), t_final));

  EinsteinNorms c_errors = c_errors_norms.first;
  EinsteinNorms c_norms = c_errors_norms.second;
  double c_error_starD1 = (c_norms.D[0] < 1e-12 ? c_errors.D[0] : c_errors.D[0] / c_norms.D[0]);
  double c_error_starD2 = (c_norms.D[1] < 1e-12 ? c_errors.D[1] : c_errors.D[1] / c_norms.D[1]);
  double c_error_starD3 = (c_norms.D[2] < 1e-12 ? c_errors.D[2] : c_errors.D[2] / c_norms.D[2]);
  double c_error_theta1 = (c_norms.theta[0] < 1e-12 ? c_errors.theta[0] : c_errors.theta[0] / c_norms.theta[0]);
  double c_error_theta2 = (c_norms.theta[1] < 1e-12 ? c_errors.theta[1] : c_errors.theta[1] / c_norms.theta[1]);
  double c_error_theta3 = (c_norms.theta[2] < 1e-12 ? c_errors.theta[2] : c_errors.theta[2] / c_norms.theta[2]);
  double c_error_starB1 = (c_norms.B[0] < 1e-12 ? c_errors.B[0] : c_errors.B[0] / c_norms.B[0]);
  double c_error_starB2 = (c_norms.B[1] < 1e-12 ? c_errors.B[1] : c_errors.B[1] / c_norms.B[1]);
  double c_error_starB3 = (c_norms.B[2] < 1e-12 ? c_errors.B[2] : c_errors.B[2] / c_norms.B[2]);
  
  std::cout << FORMAT(30) << "\tCont. norms: *D= "  << FORMAT(15) << c_norms.D[0]  << FORMAT(15) << c_norms.D[1]  << FORMAT(15) << c_norms.D[2] <<std::endl;
  std::cout << FORMAT(30) << "\tCont. norms: *Theta= "  << FORMAT(15) << c_norms.theta[0]  << FORMAT(15) << c_norms.theta[1]  << FORMAT(15) << c_norms.theta[2] <<std::endl;
  std::cout << FORMAT(30) << "\tCont. norms: *B= "  << FORMAT(15) << c_norms.B[0]  << FORMAT(15) << c_norms.B[1]  << FORMAT(15) << c_norms.B[2] <<std::endl;

  std::cout << FORMAT(30) << "\tCont. rel. errors: *D= "  << FORMAT(15) << c_error_starD1  << FORMAT(15) << c_error_starD2  << FORMAT(15) << c_error_starD3 << std::endl;
  std::cout << FORMAT(30) << "\tCont. rel. errors: *Theta= " << FORMAT(15) << c_error_theta1 << FORMAT(15) << c_error_theta2 << FORMAT(15) << c_error_theta3 << std::endl;
  std::cout << FORMAT(30) << "\tCont. rel. errors: *B= "  << FORMAT(15) << c_error_starB1  << FORMAT(15) << c_error_starB2  << FORMAT(15) << c_error_starB3 << std::endl;

  timer.stop();
  double t_wall_norms, t_proc_norms;
  std::tie(t_wall_norms, t_proc_norms) = store_times(timer, "[main] Time errors (wall/proc) ");

  Eigen::VectorXd constraints = es.computeConstraints(es.contractPara<double, Einstein::TZeroFormType>(lapse, t_final),
                                                      es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(theta, t_initial),
                                                      es.contractParaVec<Eigen::Vector3d, Einstein::OneFormType>(theta, t_final), 
                                                      starDtB_0, 
                                                      starDtB_n,
                                                      solver);

                                                      std::cout << "[main] Constraint 1: " << FORMAT(15) << constraints[0]  << FORMAT(15) << constraints[1]  << FORMAT(15) << constraints[2] <<std::endl;
                                                      std::cout << "[main] Constraint 2: " << FORMAT(15) << constraints[3]  << FORMAT(15) << constraints[4]  << FORMAT(15) << constraints[5] <<std::endl;
                                                      std::cout << "[main] Minimum det P Theta: " << FORMAT(15) << es.min_detTheta() << std::endl;
                                                      std::cout << "[main] Maximum residual: " << FORMAT(15) << max_residual << std::endl;
                                                      std::cout << "[main] Mesh diameter " << mesh_ptr->h_max() << std::endl;

  // --------------------------------------------------------------------------
  //                     Creates .txt file with data and results
  // --------------------------------------------------------------------------

  std::ofstream out("results.txt");
  out << "Solution: " << solution << std::endl;
  out << "Implicit: " << is_implicit << std::endl;
  out << "StabilisationParameter: " << es.stabilizationParameter() << std::endl;
  out << "Mesh: " << mesh_file << std::endl;
  out << "Degree: " << K << std::endl;
  out << "MeshSize: " << mesh_ptr->h_max() << std::endl;
  out << "NbCells: " << mesh_ptr->n_cells() << std::endl;
  out << "NbFaces: " << mesh_ptr->n_faces() << std::endl;
  out << "NbEdges: " << mesh_ptr->n_edges() << std::endl;
  out << "DimX1: " << es.ddrSpaces().dofspace(1).dimension() << std::endl;
  out << "DimX2: " << es.ddrSpaces().dofspace(2).dimension() << std::endl;
  out << "Sizesystem: " << es.dimensionSpace() << std::endl;
  out << "InitialTime: " << t_initial << std::endl;
  out << "FinalTime: " << t_final << std::endl;
  out << "Timestep: " << dt << std::endl;
  out << "MinDetTheta: " << es.min_detTheta() << std::endl;
  out << "MaxResidual: " << max_residual << std::endl;
  // Discrete errors and norms
  out << std::scientific;
  out << "E_L2starD1: " << (norms.D[0] < 1e-12 ? errors.D[0] : errors.D[0] / norms.D[0]) << std::endl;
  out << "E_L2starD2: " << (norms.D[1] < 1e-12 ? errors.D[1] : errors.D[1] / norms.D[1]) << std::endl;
  out << "E_L2starD3: " << (norms.D[2] < 1e-12 ? errors.D[2] : errors.D[2] / norms.D[2]) << std::endl;
  out << "E_L2starTheta1: " << (norms.theta[0] < 1e-12 ? errors.theta[0] : errors.theta[0] / norms.theta[0]) << std::endl;
  out << "E_L2starTheta2: " << (norms.theta[1] < 1e-12 ? errors.theta[1] : errors.theta[1] / norms.theta[0]) << std::endl;
  out << "E_L2starTheta3: " << (norms.theta[2] < 1e-12 ? errors.theta[2] : errors.theta[2] / norms.theta[0]) << std::endl;
  out << "E_L2starB1: " << (norms.B[0] < 1e-12 ? errors.B[0] : errors.B[0] / norms.B[0]) << std::endl;
  out << "E_L2starB2: " << (norms.B[1] < 1e-12 ? errors.B[1] : errors.B[1] / norms.B[1])  << std::endl;
  out << "E_L2starB3: " << (norms.B[2] < 1e-12 ? errors.B[2] : errors.B[2] / norms.B[2])  << std::endl;
  out << "absE_L2starD1: " << errors.D[0] << std::endl;
  out << "absE_L2starD2: " << errors.D[1]<< std::endl;
  out << "absE_L2starD3: " << errors.D[2] << std::endl;
  out << "relE_L2starD: " << ((c_norms.D[0] + c_norms.D[1] + c_norms.D[2]) < 1e-12 ? (errors.D[0] + errors.D[1] + errors.D[2]) : (errors.D[0] + errors.D[1] + errors.D[2])/(c_norms.D[0] + c_norms.D[1] + c_norms.D[2])) << std::endl;
  out << "relContE_L2starD: " << ((c_norms.D[0] + c_norms.D[1] + c_norms.D[2]) < 1e-12 ? (c_errors.D[0] + c_errors.D[1] + c_errors.D[2]) : (c_errors.D[0] + c_errors.D[1] + c_errors.D[2])/(c_norms.D[0] + c_norms.D[1] + c_norms.D[2])) << std::endl;
  out << "absE_L2starTheta1: " << errors.theta[0] << std::endl;
  out << "absE_L2starTheta2: " << errors.theta[1] << std::endl;
  out << "absE_L2starTheta3: " << errors.theta[2] << std::endl;
  out << "relE_L2starTheta: " << ((c_norms.theta[0] + c_norms.theta[1] + c_norms.theta[2]) < 1e-12 ? (errors.theta[0] + errors.theta[1] + errors.theta[2]) : (errors.theta[0] + errors.theta[1] + errors.theta[2])/(c_norms.theta[0] + c_norms.theta[1] + c_norms.theta[2])) << std::endl;
  out << "relContE_L2starTheta: " << ((c_norms.theta[0] + c_norms.theta[1] + c_norms.theta[2]) < 1e-12 ? (c_errors.theta[0] + c_errors.theta[1] + c_errors.theta[2]) : (c_errors.theta[0] + c_errors.theta[1] + c_errors.theta[2])/(c_norms.theta[0] + c_norms.theta[1] + c_norms.theta[2])) << std::endl;
  out << "absE_L2starB1: " << errors.B[0] << std::endl;
  out << "absE_L2starB2: " << errors.B[1] << std::endl;
  out << "absE_L2starB3: " << errors.B[2] << std::endl;
  out << "relE_L2starB: " << ((c_norms.B[0] + c_norms.B[1] + c_norms.B[2]) < 1e-12 ? (errors.B[0] + errors.B[1] + errors.B[2]) : (errors.B[0] + errors.B[1] + errors.B[2])/(c_norms.B[0] + c_norms.B[1] + c_norms.B[2])) << std::endl;
  out << "relContE_L2starB: " << ((c_norms.B[0] + c_norms.B[1] + c_norms.B[2]) < 1e-12 ? (c_errors.B[0] + c_errors.B[1] + c_errors.B[2]) : (c_errors.B[0] + c_errors.B[1] + c_errors.B[2])/(c_norms.B[0] + c_norms.B[1] + c_norms.B[2])) << std::endl;
  out << "N_L2starD1: " << norms.D[0] << std::endl;
  out << "N_L2starD2: " << norms.D[1] << std::endl;
  out << "N_L2starD3: " << norms.D[2] << std::endl;
  out << "N_L2starTheta1: " << norms.theta[0] << std::endl;
  out << "N_L2starTheta2: " << norms.theta[1] << std::endl;
  out << "N_L2starTheta3: " << norms.theta[2] << std::endl;
  out << "N_L2starB1: " << norms.B[0] << std::endl;
  out << "N_L2starB2: " << norms.B[1] << std::endl;
  out << "N_L2starB3: " << norms.B[2] << std::endl;
  // Continuous errors and norms
  out << "E_cstarD1: " << (c_norms.D[0] < 1e-12 ? c_errors.D[0] : c_errors.D[0] / c_norms.D[0]) << std::endl;
  out << "E_cstarD2: " << (c_norms.D[1] < 1e-12 ? c_errors.D[1] : c_errors.D[1] / c_norms.D[1]) << std::endl;
  out << "E_cstarD3: " << (c_norms.D[2] < 1e-12 ? c_errors.D[2] : c_errors.D[2] / c_norms.D[2]) << std::endl;
  out << "E_ctheta1: " << (c_norms.theta[0] < 1e-12 ? c_errors.theta[0] : c_errors.theta[0]/c_norms.theta[0]) << std::endl;
  out << "E_ctheta2: " << (c_norms.theta[1] < 1e-12 ? c_errors.theta[1] : c_errors.theta[1]/c_norms.theta[1]) << std::endl;
  out << "E_ctheta3: " << (c_norms.theta[2] < 1e-12 ? c_errors.theta[2] : c_errors.theta[2]/c_norms.theta[2]) << std::endl;
  out << "E_cstarB1: " << (c_norms.B[0] < 1e-12 ? c_errors.B[0] : c_errors.B[0] / c_norms.B[0]) << std::endl;
  out << "E_cstarB2: " << (c_norms.B[1] < 1e-12 ? c_errors.B[1] : c_errors.B[1] / c_norms.B[1]) << std::endl;
  out << "E_cstarB3: " << (c_norms.B[2] < 1e-12 ? c_errors.B[2] : c_errors.B[2] / c_norms.B[2]) << std::endl;
  out << "absE_cstarD1: " << c_errors.D[0] << std::endl;
  out << "absE_cstarD2: " << c_errors.D[1] << std::endl;
  out << "absE_cstarD3: " << c_errors.D[2] << std::endl;
  out << "absE_ctheta1: " << c_errors.theta[0] << std::endl;
  out << "absE_ctheta2: " << c_errors.theta[1] << std::endl;
  out << "absE_ctheta3: " << c_errors.theta[2] << std::endl;
  out << "absE_cstarB1: " << c_errors.B[0] << std::endl;
  out << "absE_cstarB2: " << c_errors.B[1] << std::endl;
  out << "absE_cstarB3: " << c_errors.B[2] << std::endl;
  out << "N_cL2starD1: " << c_norms.D[0] << std::endl;
  out << "N_cL2starD2: " << c_norms.D[1] << std::endl;
  out << "N_cL2starD3: " << c_norms.D[2] << std::endl;
  out << "N_cL2theta1: " << c_norms.theta[0] << std::endl;
  out << "N_cL2theta2: " << c_norms.theta[1] << std::endl;
  out << "N_cL2theta3: " << c_norms.theta[2] << std::endl;
  out << "N_cL2starB1: " << c_norms.B[0] << std::endl;
  out << "N_cL2starB2: " << c_norms.B[1] << std::endl;
  out << "N_cL2starB3: " << c_norms.B[2] << std::endl;
  out << "ConstraintTheta1: " << constraints[0] << std::endl;
  out << "ConstraintTheta2: " << constraints[1] << std::endl;
  out << "ConstraintTheta3: " << constraints[2] << std::endl;
  out << "ConstraintB1: " << constraints[3] << std::endl;
  out << "ConstraintB2: " << constraints[4] << std::endl;
  out << "ConstraintB3: " << constraints[5] << std::endl;
  out << "TwallDDRCore: " << t_wall_ddrcore << std::endl;  
  out << "TprocDDRCore: " << t_proc_ddrcore << std::endl;  
  out << "TwallAssemble: " << t_wall_assemble << std::endl;
  out << "TprocAssemble: " << t_proc_assemble << std::endl;
  out << "TwallSolve: " << t_wall_solve << std::endl;  
  out << "TprocSolve: " << t_proc_solve << std::endl;  
  out << "TwallRun: " << t_wall_run << std::endl;
  out << "TprocRun: " << t_proc_run << std::endl;
  out << "TwallNorms: " << t_wall_norms << std::endl;  
  out << "TprocNorms: " << t_proc_norms << std::endl;
  out << std::flush;
  out.close();

  std::cout << "[main] Done" << std::endl;
  return 0;
}

//------------------------------------------------------------------------------
// Einstein
//------------------------------------------------------------------------------

Einstein::Einstein(
               const DDR_Spaces & ddrspaces,
               bool use_threads,
               double stab_par,
               std::ostream & output
               )
  : m_ddrspaces(ddrspaces),
    m_use_threads(use_threads),
    m_output(output),
    m_A(6*m_ddrspaces.dofspace(1).dimension(), 6*m_ddrspaces.dofspace(1).dimension()),
    m_b(6*m_ddrspaces.dofspace(1).dimension()),
    m_L2_X1(m_ddrspaces.dofspace(1).dimension(), m_ddrspaces.dofspace(1).dimension()),
    m_L2_X2(m_ddrspaces.dofspace(2).dimension(), m_ddrspaces.dofspace(2).dimension()),
    m_A_theta(3*m_ddrspaces.dofspace(2).dimension(), 3*m_ddrspaces.dofspace(1).dimension()),
    m_b_theta(m_ddrspaces.dofspace(2).dimension()*3),
    m_stab_par(stab_par),
    m_min_detTheta(100)
{
  m_output << "[Einstein] Initializing" << std::endl;
}

//------------------------------------------------------------------------------

void Einstein::precomputeL2()
{  
  auto _assemble_local_product = [&](size_t iT, const SystemVectors<Eigen::MatrixXd> & lsT, std::vector<std::list<Eigen::Triplet<double> > > & triplets_sys, std::vector<Eigen::VectorXd> & vecs) -> void
  {
    const Cell & T = *m_ddrspaces.mesh().cell(iT);

    std::vector<size_t> IT_sysX1 = m_ddrspaces.dofspace(1).globalDOFIndices(T);
    std::vector<size_t> IT_sysX2 = m_ddrspaces.dofspace(2).globalDOFIndices(T);
      
    for (size_t i = 0; i < IT_sysX1.size(); i++){
      for (size_t j = 0; j < IT_sysX1.size(); j++){
        triplets_sys[0].emplace_back(IT_sysX1[i], IT_sysX1[j], lsT.systems[0](i, j));
      }
    }
    for (size_t i = 0; i < IT_sysX2.size(); i++){
      for (size_t j = 0; j < IT_sysX2.size(); j++){
        triplets_sys[1].emplace_back(IT_sysX2[i], IT_sysX2[j], lsT.systems[1](i, j));
      }
    }
  };

  auto _compute_local_product = [&](size_t iT) -> SystemVectors<Eigen::MatrixXd>
  {
    // Product matrices to assemble and compute rhs
    Eigen::MatrixXd L2_X1 = m_ddrspaces.computeL2Product(1, iT, m_stab_par);
    Eigen::MatrixXd L2_X2 = m_ddrspaces.computeL2Product(2, iT, m_stab_par);
    // Vector of matrices and vectors to be assembled
    SystemVectors<Eigen::MatrixXd> to_assemble{{L2_X1, L2_X2}, {}};

    return to_assemble;
  };

  // Assemble all local contributions
  auto assemble_all = [&](size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs) -> void
  {
    for (size_t iT = start; iT < end; iT++){
      _assemble_local_product(iT, _compute_local_product(iT), *triplets, *vecs);
    }
  };

  // Assemble the L2 matrices
  size_t dim_X1 = m_ddrspaces.dofspace(1).dimension();
  size_t dim_X2 = m_ddrspaces.dofspace(2).dimension();

  auto [systems, vectors] = parallel_assembly_system(m_ddrspaces.mesh().n_cells(), {std::make_pair(dim_X1, dim_X1), std::make_pair(dim_X2, dim_X2)}, {}, assemble_all, m_use_threads);

  m_L2_X1 = systems[0];
  m_L2_X2 = systems[1];
}

//------------------------------------------------------------------------------

void Einstein::assembleLinearSystem(
                                    double dt,
                                    const ZeroFormType & lapse,
                                    const OneFormType & E0,
                                    const Eigen::VectorXd & starDtB
                                    )
{  
  // Assemble all local contributions
  auto assemble_all = [&](size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vec) -> void
  {
    for (size_t iT = start; iT < end; iT++){
      _assemble_local_contribution(iT, _compute_local_contribution(dt, iT, lapse, E0, starDtB), *triplets, *vec);
    }
  };

  size_t dim_X1 = m_ddrspaces.dofspace(1).dimension();
  size_t dim_X2 = m_ddrspaces.dofspace(2).dimension();
  
  // Assemble the matrix and rhs
  auto [systems, vectors] = parallel_assembly_system(m_ddrspaces.mesh().n_cells(), {std::make_pair(6*dim_X1, 6*dim_X1), std::make_pair(3*dim_X2, 3*dim_X1)}, {6*dim_X1, 3*dim_X2, m_ddrspaces.mesh().n_cells()}, assemble_all, m_use_threads);

  m_A = systems[0];
  m_A_theta = systems[1];
  m_b = vectors[0];
  m_b_theta = vectors[1];
  if (vectors[2].minCoeff() < m_min_detTheta){
    m_min_detTheta = vectors[2].minCoeff();
  }
}

//------------------------------------------------------------------------------

void Einstein::_assemble_local_contribution(
                                          size_t iT,
                                          const SystemVectors<Eigen::MatrixXd> & lsT,
                                          std::vector<std::list<Eigen::Triplet<double>>> & triplets_sys,
                                          std::vector<Eigen::VectorXd> & vecs
                                          )
{
  const Cell & T = *m_ddrspaces.mesh().cell(iT);

  std::vector<size_t> IT_sysX1 = m_ddrspaces.dofspace(1).globalDOFIndices(T);
  std::vector<size_t> IT_sysX2 = m_ddrspaces.dofspace(2).globalDOFIndices(T);

  size_t dim_X1_T = m_ddrspaces.dofspace(1).dimensionCell(T);
  size_t dim_X2_T = m_ddrspaces.dofspace(2).dimensionCell(T);
  size_t dim_X1 = m_ddrspaces.dofspace(1).dimension();
  size_t dim_X2 = m_ddrspaces.dofspace(2).dimension();

  // Matrix for equations 1 & 3
  for (size_t i = 0; i < dim_X1_T; i++){
    for (size_t j = 0; j < dim_X1_T; j++){
      for (size_t k = 0; k < 3; k++){
        // L2 products only on the diagonals
        triplets_sys[0].emplace_back(dim_X1*k+IT_sysX1[i], dim_X1*k+IT_sysX1[j], lsT.systems[0](dim_X1_T*k+i, dim_X1_T*k+j));
        triplets_sys[0].emplace_back(3*dim_X1 + dim_X1*k+IT_sysX1[i], 3*dim_X1 + dim_X1*k+IT_sysX1[j], lsT.systems[0](3*dim_X1_T + dim_X1_T*k+i, 3*dim_X1_T + dim_X1_T*k+j));
        for (size_t l = 0; l < 3; l++){
          triplets_sys[0].emplace_back(dim_X1*k+IT_sysX1[i], 3*dim_X1 + dim_X1*l+IT_sysX1[j], lsT.systems[0](dim_X1_T*k+i, 3*dim_X1_T + dim_X1_T*l+j));
          triplets_sys[0].emplace_back(3*dim_X1 + dim_X1*k+IT_sysX1[i], dim_X1*l+IT_sysX1[j], lsT.systems[0](3*dim_X1_T + dim_X1_T*k+i, dim_X1_T*l+j));
        }
      }
    }
  }
  // Matrix for equation 2
  for (size_t i = 0; i < dim_X2_T; i++){
    for (size_t j = 0; j < dim_X1_T; j++){
      for (size_t k = 0; k < 3; k++){
        for (size_t l = 0; l < 3; l++){
          triplets_sys[1].emplace_back(dim_X2*k+IT_sysX2[i], dim_X1*l+IT_sysX1[j], lsT.systems[1](dim_X2_T*k+i, dim_X1_T*l+j));
        }
      }
    }
  }

  // Rhs for equations 1 & 3
  for (size_t i = 0; i < dim_X1_T; i++){
    for (size_t j = 0; j < 3; j++){
      vecs[0][dim_X1*j + IT_sysX1[i]] += lsT.vectors[j](i);
      vecs[0][3*dim_X1+dim_X1*j + IT_sysX1[i]] += lsT.vectors[j][dim_X1_T+dim_X2_T+i];
    }
  }
  // Rhs for equation 2
  for (size_t i = 0; i < dim_X2_T; i++){
    for (size_t j = 0; j < 3; j++){
      vecs[1][j*dim_X2+IT_sysX2[i]] += lsT.vectors[j][dim_X1_T+i];
    }
  }
  vecs[2][iT] = lsT.vectors[3][0];
}

//------------------------------------------------------------------------------

SystemVectors<Eigen::MatrixXd> Einstein::_compute_local_contribution(
                                                                    const double dt,                      
                                                                    const size_t iT,                      
                                                                    const ZeroFormType & lapse,
                                                                    const OneFormType & E0,
                                                                    const Eigen::VectorXd & starDtB
                                                                    )
{
  Cell & T = *m_ddrspaces.mesh().cell(iT);
  // Dimensions to calculate offsets
  const size_t dim_X1_T = m_ddrspaces.dofspace(1).dimensionCell(T);
  const size_t dim_X2_T = m_ddrspaces.dofspace(2).dimensionCell(T);
  const size_t dim_X1 = m_ddrspaces.dofspace(1).dimension();
  const size_t dim_X2 = m_ddrspaces.dofspace(2).dimension();

  // Hodge star for 1 and 2 forms
  Eigen::MatrixXd hodge1forms_T = m_ddrspaces.get_hodge_star(1, 3, iT);
  Eigen::MatrixXd hodge2forms_T = m_ddrspaces.get_hodge_star(2, 3, iT);

  // L2 products
  const Eigen::MatrixXd & L2_X1_T = m_ddrspaces.computeL2Product(1, iT, m_stab_par);
  const Eigen::MatrixXd & L2_X2_T = m_ddrspaces.computeL2Product(2, iT, m_stab_par);

  // Calculate the value of H, E, U at the center of mass of T
  auto [int_NU_P, int_NE0_d, int_NP, int_Nd, det_theta] = _compute_nonlinear_relations(iT, lapse, E0, starDtB);

  // Vector of matrices and vectors to be assembled
  SystemVectors<Eigen::MatrixXd> to_assemble{{Eigen::MatrixXd::Zero(6*dim_X1_T, 6*dim_X1_T), dt*int_NP}, {Eigen::VectorXd::Zero(dim_X1_T+dim_X2_T+dim_X1_T), Eigen::VectorXd::Zero(dim_X1_T+dim_X2_T+dim_X1_T), Eigen::VectorXd::Zero(dim_X1_T+dim_X2_T+dim_X1_T), Eigen::VectorXd::Zero(1)}};

  // System matrix of first equation
  to_assemble.systems[0].block(0, 0, 3*dim_X1_T, 3*dim_X1_T) = Eigen::KroneckerProduct(Eigen::Matrix3d::Identity(), L2_X1_T);
  to_assemble.systems[0].block(0, 3*dim_X1_T, 3*dim_X1_T, 3*dim_X1_T) = -dt*int_Nd;
  // System matrix of third equation
  to_assemble.systems[0].block(3*dim_X1_T, 0, 3*dim_X1_T, 3*dim_X1_T) =  dt*int_Nd;
  to_assemble.systems[0].block(3*dim_X1_T, 3*dim_X1_T, 3*dim_X1_T, 3*dim_X1_T) = Eigen::KroneckerProduct(Eigen::Matrix3d::Identity(), L2_X1_T);

  for (size_t i=0; i<3; i++){
    // RHS of first equation  (starDh, wh)_X1
    to_assemble.vectors[i].head(dim_X1_T) += L2_X1_T * m_ddrspaces.dofspace(1).restrictCell(iT, starDtB.segment(dim_X1*i, dim_X1));
    // RHS of second equation  (starThetah, vh)_X1
    to_assemble.vectors[i].segment(dim_X1_T, dim_X2_T) += L2_X2_T * m_ddrspaces.dofspace(2).restrictCell(iT, starDtB.segment(dim_X1*3 + dim_X2*i, dim_X2));
    // RHS of third equation  (starBh, uh)_X1
    to_assemble.vectors[i].tail(dim_X1_T) += L2_X1_T * m_ddrspaces.dofspace(1).restrictCell(iT, starDtB.segment(dim_X1*3 + dim_X2*3 + dim_X1*i, dim_X1));
    // RHS of first equation (starUh, wh)_X1
    to_assemble.vectors[i].head(dim_X1_T) += (dt * int_NU_P + dt * int_NE0_d).segment(dim_X1_T*i, dim_X1_T);
  }

  to_assemble.vectors[3](0) = det_theta;

  return to_assemble;
}

//------------------------------------------------------------------------------

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd, Eigen::MatrixXd, double>
    Einstein::_compute_nonlinear_relations(
                                          const size_t iT,        
                                          const ZeroFormType & lapse,
                                          const OneFormType & E0,
                                          const Eigen::VectorXd & starDtB
                                          )
{
  const Cell & T = *m_ddrspaces.mesh().cell(iT);

  // Dimensions to calculate offsets
  size_t dim_X1 = m_ddrspaces.dofspace(1).dimension();
  size_t dim_X2 = m_ddrspaces.dofspace(2).dimension();
  size_t dim_X1_T = m_ddrspaces.dofspace(1).dimensionCell(iT);
  size_t dim_X2_T = m_ddrspaces.dofspace(2).dimensionCell(iT);

  // Potentials and exterior derivative
  const Eigen::MatrixXd & Pot_X1_T = m_ddrspaces.potential(1, 3, iT);   // 2-form
  const Eigen::MatrixXd & Pot_X2_T = m_ddrspaces.potential(2, 3, iT);   // 1-form
  const Eigen::MatrixXd & d_X1_T = m_ddrspaces.full_diff(1, 3, iT);     // 1-form

  // Hodge star for 1 and 2 forms
  Eigen::MatrixXd hodge1forms_T = m_ddrspaces.get_hodge_star(1, 3, iT);
  Eigen::MatrixXd hodge2forms_T = m_ddrspaces.get_hodge_star(2, 3, iT);

  // Local D, theta
  Eigen::VectorXd D_T(3*dim_X1_T), theta_T(3*dim_X2_T), B_T(3*dim_X1_T);
  for (size_t i=0; i<3; i++){
    D_T.segment(i*dim_X1_T, dim_X1_T) = m_ddrspaces.dofspace(1).restrictCell(iT, starDtB.segment(dim_X1*i, dim_X1));
    theta_T.segment(i*dim_X2_T, dim_X2_T) = m_ddrspaces.dofspace(2).restrictCell(iT, starDtB.segment(dim_X1*3+dim_X2*i, dim_X2));
    B_T.segment(i*dim_X1_T, dim_X1_T) = m_ddrspaces.dofspace(1).restrictCell(iT, starDtB.segment(dim_X1*3+dim_X2*3+dim_X1*i, dim_X1));
  }

  Eigen::MatrixXd Id3 = Eigen::MatrixXd::Identity(3,3);
  double det_theta = 1000;

  //----------------------------------------------------------------------------
  // Quadrature rule to calculate integral matrices
  //----------------------------------------------------------------------------

  Eigen::MatrixXd int_NP(Eigen::MatrixXd::Zero(3*dim_X2_T, 3*dim_X1_T)), int_Nd(Eigen::MatrixXd::Zero(3*dim_X1_T, 3*dim_X1_T));
  Eigen::VectorXd int_NU_P(Eigen::VectorXd::Zero(3*dim_X1_T)), int_NE0_d(Eigen::VectorXd::Zero(3*dim_X1_T));

  // Degree of quadratures
  const size_t qr2k = 2 * m_ddrspaces.degree();
  QuadratureRule quad_qr2k_T = generate_quadrature_rule(T, qr2k);
  size_t dim_1PkT = Manicore::Dimension::PLDim(m_ddrspaces.degree(), 1, 3);
  {
    for (size_t iqn = 0; iqn < quad_qr2k_T.size(); iqn++){
      double lapse_iqn = lapse(quad_qr2k_T[iqn].vector());

      Eigen::Matrix3d sP2sTheta_iqn, invtheta_iqn, sP1sD_iqn, sP1sB_iqn;
      // Calculate theta and inv theta at iqn
      for (size_t i=0; i<3; i++){
        sP1sD_iqn.col(i) =  m_ddrspaces.evaluate_basis(quad_qr2k_T[iqn].vector(), 2, 3, iT, Pot_X1_T*D_T.segment(i*dim_X1_T, dim_X1_T));
        sP2sTheta_iqn.col(i) = m_ddrspaces.evaluate_basis(quad_qr2k_T[iqn].vector(), 1, 3, iT, Pot_X2_T * theta_T.segment(i*dim_X2_T, dim_X2_T));        
        sP1sB_iqn.col(i) = m_ddrspaces.evaluate_basis(quad_qr2k_T[iqn].vector(), 2, 3, iT, Pot_X1_T*B_T.segment(i*dim_X1_T, dim_X1_T));
      }
      if (det_theta > sP2sTheta_iqn.determinant()){
        det_theta = sP2sTheta_iqn.determinant();
      }
      invtheta_iqn = sP2sTheta_iqn.colPivHouseholderQr().solve(Id3);

      for (size_t i=0; i<dim_X1_T; i++){
        
        Eigen::MatrixXd evalStarPot_X1(3, dim_X1_T);
        for (size_t v_i = 0; v_i < dim_X1_T; v_i++){
          evalStarPot_X1.col(v_i) = m_ddrspaces.evaluate_basis(quad_qr2k_T[iqn].vector(), 2, 3, iT, Pot_X1_T.col(v_i));
        } 
        Eigen::MatrixXd orthPot_X1 = basisChange(2, evalStarPot_X1, invtheta_iqn);

        for (size_t j=0; j<3; j++){
          Eigen::Matrix3d basis_Di = Eigen::Matrix3d::Zero();
          Eigen::Matrix3d basis_Bi = Eigen::Matrix3d::Zero();
          basis_Bi.col(j) = orthPot_X1.col(i);
          basis_Di.col(j) = orthPot_X1.col(i);
          /// Calculating the formula using the orthonormal relations
          Eigen::Matrix3d formula_H = (hodge2forms_T * basis_Bi).transpose() - 0.5 * Eigen::MatrixXd::Identity(3,3)*(hodge2forms_T * basis_Bi).trace();
          Eigen::Matrix3d formula_E = -(hodge2forms_T * basis_Di).transpose() + 0.5 * Eigen::MatrixXd::Identity(3,3)*(hodge2forms_T * basis_Di).trace();
          // Orthonormal to canonical basis change
          Eigen::Matrix3d H_iqn = basisChange(1, formula_H, sP2sTheta_iqn);
          Eigen::Matrix3d E_iqn = basisChange(1, formula_E, sP2sTheta_iqn);

          Eigen::MatrixXd int_NH_phi(Eigen::MatrixXd::Zero(dim_1PkT, 3));
          Eigen::MatrixXd int_NE_phi(Eigen::MatrixXd::Zero(dim_1PkT, 3));

          for (size_t k = 0; k < dim_1PkT; ++k) {
            for (int l = 0; l < 3; l++) {
              int_NH_phi(k, l) += quad_qr2k_T[iqn].w*lapse_iqn*(H_iqn.col(l)).dot(m_ddrspaces.evaluate_basis(quad_qr2k_T[iqn].vector(), 1, 3, iT, k));
              int_NE_phi(k, l) += quad_qr2k_T[iqn].w*lapse_iqn*(E_iqn.col(l)).dot(m_ddrspaces.evaluate_basis(quad_qr2k_T[iqn].vector(), 1, 3, iT, k));
            }
          }
          Eigen::MatrixXd int_NH_d = d_X1_T.transpose() * int_NH_phi;
          Eigen::MatrixXd int_NE_P = Pot_X2_T.transpose() * int_NE_phi;
          for (size_t k=0; k<3; k++){
            int_Nd.col(j*dim_X1_T + i).segment(k*dim_X1_T, dim_X1_T) += int_NH_d.col(k);
            int_NP.col(j*dim_X1_T + i).segment(k*dim_X2_T, dim_X2_T) += int_NE_P.col(k);
          }
        }
      }
      
      Eigen::Vector3d E0_iqn = E0(quad_qr2k_T[iqn].vector());

      // Canonical to orthonormal basis change
      Eigen::Vector3d orthE0 = basisChange(1, E0_iqn, invtheta_iqn);
      Eigen::Matrix3d orthStarP2B = hodge2forms_T * basisChange(2, sP1sB_iqn, invtheta_iqn);
      Eigen::Matrix3d orthStarP2D = hodge2forms_T * basisChange(2, sP1sD_iqn, invtheta_iqn);

      //--------------------------------------
      // Local H, E, U
      //--------------------------------------

      Eigen::Matrix3d H, orthH, E, orthE, orthStarU, orthU, U;

      /// Calculating the H, E using the orthonormal relations
      orthH =  orthStarP2B.transpose() - 0.5 * Eigen::MatrixXd::Identity(3,3)*orthStarP2B.trace() + expand2Form(hodge2forms_T*orthE0).transpose();
      orthE = -orthStarP2D.transpose() + 0.5 * Eigen::MatrixXd::Identity(3,3)*orthStarP2D.trace();
      
      // Reconstruct D_0 using orthonormal relations
      Eigen::Vector3d orthStarD0 = {orthStarP2B(2,1)-orthStarP2B(1,2), -orthStarP2B(2,0)+orthStarP2B(0,2), orthStarP2B(1,0)-orthStarP2B(0,1)};

      // Calculate U
      orthStarU = Eigen::MatrixXd::Identity(3,3) 
                    * 0.5 * (orthE0.transpose()*orthStarD0 + (orthE.transpose()*orthStarP2D).trace() 
                                                          - (orthH.transpose()*orthStarP2B).trace())
                  - orthStarD0*orthE0.transpose() - orthStarP2D * orthE.transpose() + orthStarP2B * orthH.transpose();
      orthU = hodge1forms_T * orthStarU;

      // Orthonormal to canonical basis change
      U = basisChange(2, orthU, sP2sTheta_iqn);

      Eigen::Matrix3d Hkl_E0_T = expand2Form(hodge1forms_T * E0_iqn).transpose();

      Eigen::MatrixXd int_NU_phi(Eigen::MatrixXd::Zero(dim_1PkT, 3));
      Eigen::MatrixXd int_NE0_phi(Eigen::MatrixXd::Zero(dim_1PkT, 3));
      ///// Integrate values
      for (size_t k = 0; k < dim_1PkT; ++k) {
        for (int l = 0; l < 3; l++) {
          int_NU_phi(k, l) += quad_qr2k_T[iqn].w*lapse_iqn*U.col(l).dot(m_ddrspaces.evaluate_basis(quad_qr2k_T[iqn].vector(), 2, 3, iT, k));
          int_NE0_phi(k, l) += quad_qr2k_T[iqn].w*lapse_iqn*Hkl_E0_T.col(l).dot(m_ddrspaces.evaluate_basis(quad_qr2k_T[iqn].vector(), 1, 3, iT, k));
        }
      }

      for (int l = 0; l < 3; l++) {
        int_NU_P.segment(l*dim_X1_T, dim_X1_T) += Pot_X1_T.transpose() * int_NU_phi.col(l);
        int_NE0_d.segment(l*dim_X1_T, dim_X1_T) += d_X1_T.transpose() * int_NE0_phi.col(l);
      }
    }
  }

 return std::make_tuple(int_NU_P, int_NE0_d, int_NP, int_Nd, det_theta);
}

//------------------------------------------------------------------------------

Eigen::MatrixXd Einstein::basisChange(size_t k,
                                const Eigen::MatrixXd & forms, 
                                const Eigen::MatrixXd & basis 
                                ) const
{
  size_t num2F = forms.cols();
  Eigen::MatrixXd basisForms = Eigen::MatrixXd::Zero(3,num2F);
  if (k==1){
    basisForms = basis * forms;
  }
  else if (k==2){
    size_t pos=0;
    for (size_t I=0; I<2; I++){
      Eigen::VectorXd thetai_I = basis.row(I);
      Eigen::MatrixXd forms_I = contractIndex(0, thetai_I, forms);
      for (size_t J=I+1; J<3; J++){
        basisForms.row(pos) = basis.row(J) * forms_I;
        pos+=1;
      }
    }
  }
  return basisForms;
}                               

//------------------------------------------------------------------------------

/// Contracts a vector of 2-forms with the first or second index: (A_i)ij (indices=12); (A_i)ji (indices=13)
Eigen::VectorXd Einstein::contractIndex(size_t indices,
                                const Eigen::MatrixXd & forms
                                ) const
{
  Eigen::VectorXd contracted(3);
  contracted << -forms(0,1) - forms(1,2), forms(0,0) - forms(2,2), forms(1,0) + forms(2,1);
  if (indices==12){
    return contracted;
  }
  else if (indices==13){
    return -contracted;
  }
  else {
    std::cerr << "[main] ERROR: Function contractIndex unknown index" << std::endl;
    exit(1); 
  }
}                               

//------------------------------------------------------------------------------

Eigen::MatrixXd Einstein::contractIndex(size_t k,
                                const Eigen::MatrixXd & forms1,
                                const Eigen::MatrixXd & forms2
                                ) const
{
  size_t num1F = forms1.cols();
  size_t num2F = forms2.cols();
  if (num1F == 1 && num2F == 1){ 
    return Eigen::Vector<double,3>{-forms2(0,0)*forms1(1,0)-forms2(1,0)*forms1(2,0), forms2(0,0)*forms1(0,0)-forms2(2,0)*forms1(2,0), forms2(1,0)*forms1(0,0)+forms2(2,0)*forms1(1,0)};
  }
  else if (num1F == 1){ // collection of 2-forms with a 1-form: (A_i)_jk b^j (k=0); (A_i)_kj b^j
    Eigen::MatrixXd contracted = Eigen::MatrixXd::Zero(3, num2F);
    for (size_t i=0; i<num2F; i++){
      contracted(0,i) = -forms2(0,i)*forms1(1,0)-forms2(1,i)*forms1(2,0);
      contracted(1,i) =  forms2(0,i)*forms1(0,0)-forms2(2,i)*forms1(2,0);
      contracted(2,i) =  forms2(1,i)*forms1(0,0)+forms2(2,i)*forms1(1,0);
    }
    if (k==1){ // a 2-form with a collection of 1-forms in either 0th or 1st position: A_ij (b_k)^i (k=0); A_ji(b_k)^i (k=1)
      contracted = -contracted;
    }
    return contracted;
  }
  else if (num2F == 1){
    Eigen::MatrixXd contracted = Eigen::MatrixXd::Zero(3, num1F);
    for (size_t i=0; i<num1F; i++){
      contracted(0,i) = -forms2(0,0)*forms1(1,i)-forms2(1,0)*forms1(2,i);
      contracted(1,i) = forms2(0,0)*forms1(0,i)-forms2(2,0)*forms1(2,i);
      contracted(2,i) = forms2(1,0)*forms1(0,i)+forms2(2,0)*forms1(1,i);
    }
    return contracted;
  }
  else { // a vector of 2-forms with a vector of 2-forms: (A_i)_jk (B^i)^jl
    Eigen::MatrixXd contracted = Eigen::MatrixXd::Zero(3, 3);
    contracted << forms1.row(0)*forms2.row(0).transpose() + forms1.row(1)*forms2.row(1).transpose(), forms1.row(1)*forms2.row(2).transpose(), -forms1.row(0)*forms2.row(2).transpose(),
                  forms1.row(2)*forms2.row(1).transpose(), forms1.row(0)*forms2.row(0).transpose() + forms1.row(2)*forms2.row(2).transpose(), forms1.row(0)*forms2.row(1).transpose(),
                  -forms1.row(2)*forms2.row(0).transpose(), forms1.row(1)*forms2.row(0).transpose(), forms1.row(1)*forms2.row(1).transpose() + forms1.row(2)*forms2.row(2).transpose(); 
    return contracted;
  }
}

//------------------------------------------------------------------------------

Eigen::Vector3d Einstein::wedge11(
                                  const Eigen::Vector3d & form1,
                                  const Eigen::Vector3d & form2
                                  ) const
{
  return Eigen::Vector<double, 3>{form1(0)*form2(1)-form1(1)*form2(0), form1(0)*form2(2)-form1(2)*form2(0), form1(1)*form2(2)-form1(2)*form2(1)};
}

//------------------------------------------------------------------------------

Eigen::Matrix3d Einstein::expand2Form(
                                      const Eigen::Vector3d & form
                                      ) const
{
  return Eigen::Matrix<double, 3, 3>{{0, form(0), form(1)}, {-form(0), 0, form(2)}, {-form(1), -form(2), 0}};
}

//------------------------------------------------------------------------------

void Einstein::updateThetaRHS(
                              const Eigen::VectorXd & starD
                              )
{
  m_b_theta += m_A_theta * starD;
}

//------------------------------------------------------------------------------

void Einstein::computeBCRHS(
                            double dt,
                            const ZeroFormType & lapse,
                            const Collection2Forms & D,
                            const Collection1Forms & theta,
                            const Collection2Forms & B,
                            const Collection2Forms & dtD,
                            const Collection1Forms & E,
                            const Collection1Forms & H,
                            const Collection2Forms & dH
                            )
{
  size_t dim_X1 = m_ddrspaces.dofspace(1).dimension();
  // size_t dim_X2 = m_ddrspaces.dofspace(2).dimension();

  Eigen::MatrixXd interp_f(dim_X1, 3);
  Eigen::MatrixXd hodge12 = Eigen::Matrix<double,3,3>{{0.,0.,1.},{0.,-1.,0.},{1.,0.,0.}};

  // Construct forcing term for the first equation: star (dtD - d(NH) - N*U)
  for (size_t i=0; i<3; i++){
    auto f_cont =  [&, i, this](const Eigen::Vector3d & x) -> Eigen::Vector3d 
                    {
                      // Evaluate functions at x to calculate starU
                      Eigen::Matrix3d B_x, theta_x, invtheta_x, D_x;
                      for (size_t j=0; j<3; j++){
                        B_x.col(j) = B[j+1](x);
                        theta_x.col(j) = theta[j+1](x);
                        D_x.col(j) = D[j+1](x);
                      }
                      invtheta_x = theta_x.colPivHouseholderQr().solve(Eigen::Matrix3d::Identity());
                      // Canonical to orthonormal basis change
                      Eigen::Matrix3d orthB_x, orthD_x, orthH_x, orthE_x;
                      Eigen::Vector3d orthE0_x = basisChange(1, E[0](x), invtheta_x);
                      Eigen::Vector3d orthD0_x = basisChange(1, D[0](x), invtheta_x);
                      orthB_x = basisChange(2, B_x, invtheta_x);
                      orthD_x = basisChange(2, D_x, invtheta_x);
                      // Calculate orthonormal E, H using orthonormal relations
                      for (size_t j=0; j<3; j++){
                        orthE_x.col(j) = composite_orth_E[j+1](orthD_x);
                        orthH_x.col(j) = composite_orth_H[j+1](orthE0_x, orthB_x);
                      }
                      // Orthonormal to canonical basis change
                      Eigen::Vector3d U_x = basisChange(2, composite_orth_U[i+1](orthE0_x, orthD0_x, orthH_x, orthE_x, orthD_x, orthB_x), theta_x);
                      return hodge12 * (dtD[i+1](x) + lapse(x)*wedge11(E[0](x), H[i+1](x)) - lapse(x)*dH[i+1](x) - lapse(x)*U_x);
                    };
    interp_f.col(i) = m_ddrspaces.interpolate(DDR_Spaces::DDR_function_type{f_cont}, 1);
  }
  
  auto _construct_BC_RHS = [&](size_t iT) -> SystemVectors<Eigen::MatrixXd>
  {
    const Cell & T = *m_ddrspaces.mesh().cell(iT);
    size_t dim_X1_T = m_ddrspaces.dofspace(1).dimensionCell(iT);

    SystemVectors<Eigen::MatrixXd> to_assemble{{}, {3, Eigen::VectorXd::Zero(dim_X1_T+dim_X1_T)}};

    //-------------------
    // Local matrix
    //-------------------
    
    Eigen::MatrixXd L2_X1_T = m_ddrspaces.computeL2Product(1, iT, m_stab_par);

    //-----------------------
    // Local source vector
    //-----------------------

    for (size_t i=0; i<3; i++){
      // RHS of first equation
      to_assemble.vectors[i].head(dim_X1_T) += dt * L2_X1_T * m_ddrspaces.dofspace(1).restrictCell(iT, interp_f.col(i));
    }

    if (T.is_boundary()){
      for (size_t iF_T=0; iF_T<T.n_faces(); iF_T++){
        // Face
        const Face & F = *T.face(iF_T);
        
        if (F.is_boundary()){
          size_t iF = F.global_index();
          
          // Degree of quadratures to compute boundary conditions
          const size_t dqrbc = 2 * m_ddrspaces.degree() + 3;
          QuadratureRule quad_dqrbc_F = generate_quadrature_rule(F, dqrbc);
          size_t dim_1PkF = Manicore::Dimension::PLDim(m_ddrspaces.degree(), 1, 2);
          
          // Boundary condition on N*H
          {
            for (size_t i=0; i<3; i++){

              Eigen::VectorXd b(Eigen::VectorXd::Zero(dim_1PkF));
              for (size_t iqn = 0; iqn < quad_dqrbc_F.size(); iqn++){
                Eigen::Vector2d trace_NH_iqn = m_ddrspaces.mesh().boundary_orientation(3, iT, iF_T) * lapse(quad_dqrbc_F[iqn].vector())*m_ddrspaces.get_hodge_star(1, 2, iF)*m_ddrspaces.get_exterior_trace(1, 2, iF)*H[i+1](quad_dqrbc_F[iqn].vector());
                for (size_t j = 0; j < dim_1PkF; ++j) {
                  b(j) += quad_dqrbc_F[iqn].w*(trace_NH_iqn).dot(m_ddrspaces.evaluate_basis(quad_dqrbc_F[iqn].vector(), 1, 2, iF, j));
                }
              }
              // Boundary condition of first equation
              to_assemble.vectors[i].head(dim_X1_T) += dt*m_ddrspaces.dofspace(1).extendOperator(3, iT, 2, iF, b.transpose()*-m_ddrspaces.discrete_hodge_star(1, 2, iF)*m_ddrspaces.potential(1, 2, iF)).transpose();
            } // for i
          }

          // Boundary condition on N*E
          {
            for (size_t i=0; i<3; i++){
 
              Eigen::VectorXd b(Eigen::VectorXd::Zero(dim_1PkF));

              for (size_t iqn = 0; iqn < quad_dqrbc_F.size(); iqn++){
                Eigen::Vector2d trace_NE_iqn = m_ddrspaces.mesh().boundary_orientation(3, iT, iF_T) * lapse(quad_dqrbc_F[iqn].vector())*m_ddrspaces.get_hodge_star(1, 2, iF)*m_ddrspaces.get_exterior_trace(1,2,iF)*E[i+1](quad_dqrbc_F[iqn].vector());
                for (size_t j = 0; j < dim_1PkF; ++j) {
                  b(j) += quad_dqrbc_F[iqn].w*(trace_NE_iqn).dot(m_ddrspaces.evaluate_basis(quad_dqrbc_F[iqn].vector(), 1, 2, iF, j));
                }
              }
              // Boundary condition of third equation
              to_assemble.vectors[i].tail(dim_X1_T) += dt*m_ddrspaces.dofspace(1).extendOperator(3, iT, 2, iF, b.transpose()*-m_ddrspaces.discrete_hodge_star(1, 2, iF) *m_ddrspaces.potential(1, 2, iF)).transpose();
            }// for i
          }
        } // for iF
      } // for iT
    }
    return to_assemble;
  };

  auto assemble_contribution_BC_RHS = [&](size_t iT, const SystemVectors<Eigen::MatrixXd> & lsT, std::vector<std::list<Eigen::Triplet<double>>> & triplets_sys, std::vector<Eigen::VectorXd> & vecs) -> void
  {
    const Cell & T = *m_ddrspaces.mesh().cell(iT);

    std::vector<size_t> IT_sysX1 = m_ddrspaces.dofspace(1).globalDOFIndices(T);

    size_t dim_X1_T = m_ddrspaces.dofspace(1).dimensionCell(T);
    size_t dim_X1 = m_ddrspaces.dofspace(1).dimension();
      
    for (size_t i = 0; i < dim_X1_T; i++){
      for (size_t j = 0; j < 3; j++){
        vecs[j][IT_sysX1[i]] += lsT.vectors[j](i);
        vecs[j][dim_X1 + IT_sysX1[i]] += lsT.vectors[j][dim_X1_T+i];
      }
    }
  };

  auto assemble_BC_RHS = [&](size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs) -> void
  {
    for (size_t iT = start; iT < end; iT++){
      assemble_contribution_BC_RHS(iT, _construct_BC_RHS(iT), *triplets, *vecs);
    }
  };

  auto [_, vectors] = parallel_assembly_system(m_ddrspaces.mesh().n_cells(), {}, {2*dim_X1, 2*dim_X1, 2*dim_X1}, assemble_BC_RHS, m_use_threads);

  for (size_t i=0; i<3; i++){
    m_b.segment(i*dim_X1, dim_X1) += vectors[i].head(dim_X1);
    m_b.segment(3*dim_X1 + i*dim_X1, dim_X1) += vectors[i].tail(dim_X1);
  }

}

//------------------------------------------------------------------------------

std::vector<EinsteinNorms> Einstein::computeEinsteinNorms(const std::vector<Eigen::VectorXd> & list_dofs) const
{
  const size_t ncells = m_ddrspaces.mesh().n_cells();
  const size_t nb_vectors = list_dofs.size();
  const size_t nb_collection = 3;
  std::vector<Eigen::MatrixXd> local_sqnorm_starD(nb_vectors, Eigen::MatrixXd::Zero(ncells, nb_collection));
  std::vector<Eigen::MatrixXd> local_sqnorm_startheta(nb_vectors, Eigen::MatrixXd::Zero(ncells, nb_collection));
  std::vector<Eigen::MatrixXd> local_sqnorm_starB(nb_vectors, Eigen::MatrixXd::Zero(ncells, nb_collection));

  size_t dim_X1 = m_ddrspaces.dofspace(1).dimension();
  size_t dim_X2 = m_ddrspaces.dofspace(2).dimension();

  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [&, this](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){

        // Matrices required to compute the local contributions to the norms
        Eigen::MatrixXd L2_X1_T = m_ddrspaces.computeL2Product(1, iT, m_stab_par);
        Eigen::MatrixXd L2_X2_T = m_ddrspaces.computeL2Product(2, iT, m_stab_par);

        for (size_t i=0; i<nb_vectors; i++){
          for (size_t j=0; j<nb_collection; j++){
            Eigen::VectorXd starD_T = m_ddrspaces.dofspace(1).restrictCell(iT, list_dofs[i].segment(dim_X1*j, dim_X1));
            Eigen::VectorXd startheta_T = m_ddrspaces.dofspace(2).restrictCell(iT, list_dofs[i].segment(dim_X1*3+dim_X2*j, dim_X2));
            Eigen::VectorXd starB_T = m_ddrspaces.dofspace(1).restrictCell(iT, list_dofs[i].segment(dim_X1*3+dim_X2*3+dim_X1*j, dim_X1));

            local_sqnorm_starD[i](iT, j) = starD_T.transpose() * L2_X1_T * starD_T;
            local_sqnorm_startheta[i](iT, j) = startheta_T.transpose() * L2_X2_T * startheta_T;
            local_sqnorm_starB[i](iT, j) = starB_T.transpose() * L2_X1_T * starB_T;
          }
        }
      }
    };
  parallel_for(ncells, compute_local_squarednorms, m_use_threads);

  // Assemble the output
  std::vector<EinsteinNorms> list_norms;
  list_norms.resize(nb_vectors+1);
  for (size_t i=0; i<nb_vectors; i++){
    list_norms[i] = EinsteinNorms();
    for (size_t j=0; j<nb_collection; j++){
      list_norms[i].D(j) = std::sqrt(std::abs(local_sqnorm_starD[i].col(j).sum()));
      list_norms[i].theta(j) = std::sqrt(std::abs(local_sqnorm_startheta[i].col(j).sum()));
      list_norms[i].B(j) = std::sqrt(std::abs(local_sqnorm_starB[i].col(j).sum()));
    }
  }
  return list_norms;
}

//------------------------------------------------------------------------------

std::pair<EinsteinNorms, EinsteinNorms> Einstein::computeContinuousErrorsNorms(
            const Eigen::VectorXd & starDtB_h,
            const std::vector<OneFormType> & starD, 
            const std::vector<TwoFormType> & starTheta,
            const std::vector<OneFormType> & starB
            ) const
{
  const size_t ncells = m_ddrspaces.mesh().n_cells();

  Eigen::MatrixXd local_sqerrors_starD = Eigen::MatrixXd::Zero(ncells, 3);
  Eigen::MatrixXd local_sqerrors_starTheta = Eigen::MatrixXd::Zero(ncells, 3);
  Eigen::MatrixXd local_sqerrors_starB = Eigen::MatrixXd::Zero(ncells, 3);

  Eigen::MatrixXd local_sqnorms_starD = Eigen::MatrixXd::Zero(ncells, 3);
  Eigen::MatrixXd local_sqnorms_starTheta = Eigen::MatrixXd::Zero(ncells, 3);
  Eigen::MatrixXd local_sqnorms_starB = Eigen::MatrixXd::Zero(ncells, 3);
  
  size_t dim_X1 = m_ddrspaces.dofspace(1).dimension();
  size_t dim_X2 = m_ddrspaces.dofspace(2).dimension();
  
  std::function<void(size_t, size_t)> compute_local_terms
    = [&](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){

        // Cell
        const Cell & T = *m_ddrspaces.mesh().cell(iT);
        // Quadrature and evaluation of basis functions
        QuadratureRule qr_high = generate_quadrature_rule(T, 10);

        for (size_t i=0; i<3; i++){

          // Vectors representing P1 starD, P1 theta
          Eigen::VectorXd starD_T = m_ddrspaces.dofspace(1).restrictCell(iT, starDtB_h.segment(dim_X1*i, dim_X1));
          Eigen::VectorXd P1starD_T = m_ddrspaces.discrete_inv_hodge_star(2, 3, iT) * m_ddrspaces.potential(1, 3, iT) * starD_T;    // 1-form

          Eigen::VectorXd starTheta_T = m_ddrspaces.dofspace(2).restrictCell(iT, starDtB_h.segment(dim_X1*3 + dim_X2*i, dim_X2));
          Eigen::VectorXd P2starTheta_T = m_ddrspaces.discrete_inv_hodge_star(1, 3, iT) * m_ddrspaces.potential(2, 3, iT) * starTheta_T;    // 2-form

          Eigen::VectorXd starB_T = m_ddrspaces.dofspace(1).restrictCell(iT, starDtB_h.segment(dim_X1*3 + dim_X2*3 + dim_X1*i, dim_X1));
          Eigen::VectorXd P1starB_T = m_ddrspaces.discrete_inv_hodge_star(2, 3, iT) * m_ddrspaces.potential(1, 3, iT) * starB_T;    // 1-form

          for (size_t iqn = 0; iqn < qr_high.size(); iqn++){
            local_sqnorms_starD(iT, i) += qr_high[iqn].w * starD[i+1](qr_high[iqn].vector()).squaredNorm();
            local_sqnorms_starTheta(iT, i) += qr_high[iqn].w * starTheta[i+1](qr_high[iqn].vector()).squaredNorm();
            local_sqnorms_starB(iT, i) += qr_high[iqn].w * starB[i+1](qr_high[iqn].vector()).squaredNorm();

            // values of potentials/operators at quadrature node
            Eigen::Vector3d P1starD_iqn = m_ddrspaces.evaluate_basis(qr_high[iqn].vector(), 1, 3, iT, P1starD_T);
            Eigen::Vector3d P2starTheta_iqn = m_ddrspaces.evaluate_basis(qr_high[iqn].vector(), 2, 3, iT, P2starTheta_T);
            Eigen::Vector3d P1starB_iqn = m_ddrspaces.evaluate_basis(qr_high[iqn].vector(), 1, 3, iT, P1starB_T);

            // error
            local_sqerrors_starD(iT, i) += qr_high[iqn].w * (starD[i+1](qr_high[iqn].vector()) - P1starD_iqn).squaredNorm();
            local_sqerrors_starTheta(iT, i) += qr_high[iqn].w * (starTheta[i+1](qr_high[iqn].vector()) - P2starTheta_iqn).squaredNorm();
            local_sqerrors_starB(iT, i) += qr_high[iqn].w * (starB[i+1](qr_high[iqn].vector()) - P1starB_iqn).squaredNorm();
          }
        }
      }
    };
  parallel_for(ncells, compute_local_terms, m_use_threads);

  // Assemble the output

  Eigen::Vector3d error_starD, error_starTheta, error_starB, norm_starD, norm_starTheta, norm_starB;
  for (size_t i=0; i<3; i++){
    error_starD(i) = std::sqrt(std::abs(local_sqerrors_starD.col(i).sum()));
    error_starTheta(i) = std::sqrt(std::abs(local_sqerrors_starTheta.col(i).sum()));
    error_starB(i) = std::sqrt(std::abs(local_sqerrors_starB.col(i).sum()));
    norm_starD(i) = std::sqrt(std::abs(local_sqnorms_starD.col(i).sum()));
    norm_starTheta(i) = std::sqrt(std::abs(local_sqnorms_starTheta.col(i).sum()));
    norm_starB(i) = std::sqrt(std::abs(local_sqnorms_starB.col(i).sum()));
  }

  EinsteinNorms errors(error_starD, error_starTheta, error_starB);
  EinsteinNorms norms(norm_starD, norm_starTheta, norm_starB);

  return std::make_pair(errors, norms);

}

//------------------------------------------------------------------------------

double Einstein::computeResidual(const Eigen::VectorXd & starDB) const
{
  Eigen::VectorXd err = m_A * starDB - m_b;
  Eigen::VectorXd err_starD = err.head(3*m_ddrspaces.dofspace(1).dimension());
  Eigen::VectorXd err_starB = err.tail(3*m_ddrspaces.dofspace(1).dimension());
  // l2 norms
  double norm_err_starD = err_starD.norm();
  double norm_err_starB = err_starB.norm();
  double norm_b_starD = m_b.head(3*m_ddrspaces.dofspace(1).dimension()).norm();
  double norm_b_starB = m_b.tail(3*m_ddrspaces.dofspace(1).dimension()).norm();

  double total_res = (norm_b_starD ? norm_err_starD/double(norm_b_starD) : norm_err_starD) + (norm_b_starB ? norm_err_starB/double(norm_b_starB) : norm_err_starB);

  return total_res;
}

//------------------------------------------------------------------------------

Eigen::VectorXd Einstein::computeConstraints(const ZeroFormType & lapse,
                                             const std::vector<OneFormType> & theta0,
                                             const std::vector<OneFormType> & thetan,
                                             const Eigen::VectorXd & starDB0,
                                             const Eigen::VectorXd & starDBn,
                                             LinearSolver<Einstein::SystemMatrixType> & solver) const
{
  std::vector<Eigen::MatrixXd> constraint1(m_ddrspaces.mesh().n_cells(), Eigen::MatrixXd::Zero(m_ddrspaces.dofspace(1).dimension(), 3));
  std::vector<Eigen::MatrixXd> constraint2(m_ddrspaces.mesh().n_cells(), Eigen::MatrixXd::Zero(m_ddrspaces.dofspace(0).dimension(), 3));

  size_t dim_X1 = m_ddrspaces.dofspace(1).dimension();
  size_t dim_X2 = m_ddrspaces.dofspace(2).dimension();

  std::function<void(size_t, size_t)> compute_local_constraint
    = [&, this](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *m_ddrspaces.mesh().cell(iT);

        // Matrices required to compute the local contributions to the norms
        Eigen::MatrixXd L2_X1_T = m_ddrspaces.computeL2Product(1, iT, m_stab_par);
        Eigen::MatrixXd L2_X2_T = m_ddrspaces.computeL2Product(2, iT, m_stab_par);
        Eigen::MatrixXd d0_T = m_ddrspaces.compose_diff(0, 3, iT);
        Eigen::MatrixXd d1_T = m_ddrspaces.compose_diff(1, 3, iT);

        for (size_t i=0; i<3; i++){
          Eigen::VectorXd starDn_T = m_ddrspaces.dofspace(1).restrictCell(iT, starDBn.segment(dim_X1*i, dim_X1));
          Eigen::VectorXd starthetan_T = m_ddrspaces.dofspace(2).restrictCell(iT, starDBn.segment(dim_X1*3+dim_X2*i, dim_X2));
          Eigen::VectorXd starBn_T = m_ddrspaces.dofspace(1).restrictCell(iT, starDBn.segment(dim_X1*3+dim_X2*3+dim_X1*i, dim_X1));

          Eigen::VectorXd starD0_T = m_ddrspaces.dofspace(1).restrictCell(iT, starDB0.segment(dim_X1*i, dim_X1));
          Eigen::VectorXd startheta0_T = m_ddrspaces.dofspace(2).restrictCell(iT, starDB0.segment(dim_X1*3+dim_X2*i, dim_X2));
          Eigen::VectorXd starB0_T = m_ddrspaces.dofspace(1).restrictCell(iT, starDB0.segment(dim_X1*3+dim_X2*3+dim_X1*i, dim_X1));

          Eigen::VectorXd local_constraint1 = (starthetan_T-startheta0_T).transpose() * L2_X2_T * d1_T - (starBn_T-starB0_T).transpose() * L2_X1_T;
          Eigen::VectorXd local_constraint2 = (starBn_T-starB0_T).transpose() * L2_X1_T * d0_T;

        if (T.is_boundary()){
          for (size_t iF_T=0; iF_T<T.n_faces(); iF_T++){
            // Face
            const Face & F = *T.face(iF_T);
            
            if (F.is_boundary()){
              size_t iF = F.global_index();
              
              // Degree of quadratures to compute boundary conditions
              const size_t dqrbc = 2 * m_ddrspaces.degree() + 3;
              QuadratureRule quad_dqrbc_F = generate_quadrature_rule(F, dqrbc);
              size_t dim_1PkF = Manicore::Dimension::PLDim(m_ddrspaces.degree(), 1, 2);

              // Boundary condition on N*theta
              {
                for (size_t i=0; i<3; i++){
    
                  Eigen::VectorXd b(Eigen::VectorXd::Zero(dim_1PkF));

                  for (size_t iqn = 0; iqn < quad_dqrbc_F.size(); iqn++){
                    Eigen::Vector2d trace_NE_iqn = lapse(quad_dqrbc_F[iqn].vector())*m_ddrspaces.get_hodge_star(1, 2, iF)*m_ddrspaces.get_exterior_trace(1,2,iF)*(thetan[i+1](quad_dqrbc_F[iqn].vector())-theta0[i+1](quad_dqrbc_F[iqn].vector()));
                    if (m_ddrspaces.mesh().boundary_orientation(3, iT, iF_T) < 0){
                      trace_NE_iqn = -trace_NE_iqn;
                    }
                    for (size_t j = 0; j < dim_1PkF; ++j) {
                      b(j) += quad_dqrbc_F[iqn].w*(trace_NE_iqn).dot(m_ddrspaces.evaluate_basis(quad_dqrbc_F[iqn].vector(), 1, 2, iF, j));
                    }
                  }
                  // Boundary condition of third equation
                  local_constraint1 += m_ddrspaces.dofspace(1).extendOperator(3, iT, 2, iF, b.transpose()*-m_ddrspaces.discrete_hodge_star(1, 2, iF) *m_ddrspaces.potential(1, 2, iF)).transpose();
                }// for i
              }
            } // for iF
          } // for iT
        }

          std::vector<size_t> IT_sysX0 = m_ddrspaces.dofspace(0).globalDOFIndices(T);
          std::vector<size_t> IT_sysX1 = m_ddrspaces.dofspace(1).globalDOFIndices(T);

          for (size_t iDof=0; iDof < m_ddrspaces.dofspace(1).dimensionCell(iT); iDof++){
            constraint1[iT](IT_sysX1[iDof], i) += local_constraint1[iDof];
          }
          for (size_t iDof=0; iDof < m_ddrspaces.dofspace(0).dimensionCell(iT); iDof++){
            constraint2[iT](IT_sysX0[iDof], i) += local_constraint2[iDof];
          }
        }
      }
    };
  parallel_for(m_ddrspaces.mesh().n_cells(), compute_local_constraint, m_use_threads);

  Eigen::MatrixXd constraint1_total(Eigen::MatrixXd::Zero(m_ddrspaces.dofspace(1).dimension(), 3));
  Eigen::MatrixXd constraint2_total(Eigen::MatrixXd::Zero(m_ddrspaces.dofspace(0).dimension(), 3));
  for (size_t i=0; i<m_ddrspaces.mesh().n_cells(); i++){
    constraint1_total += constraint1[i];
    constraint2_total += constraint2[i];
  }

  return Eigen::VectorXd{{constraint1_total.col(0).norm(), constraint1_total.col(1).norm(), constraint1_total.col(2).norm(), constraint2_total.col(0).norm(), constraint2_total.col(1).norm(), constraint2_total.col(2).norm()}};
}

//------------------------------------------------------------------------------