// Author: Jerome Droniou (jerome.droniou@umontpellier.fr)

#include <fstream>
#include <iomanip>
#include <thread>
#include "ddr-elasticity_fracture.hpp"
#include <parallel_for.hpp>
#include <max_degrees_quadratures.hpp>
#include <GMpoly_cell.hpp>

#include "vtu_writer.hpp"

#include <boost/program_options.hpp>
#include <display_timer.hpp>
#include <connected_information_generator.cpp>
#include <fractureinterpolator.cpp>
#include <vtu_facewriter.cpp>
#include <typeinfo>

#define FORMAT(W)                                                       \
  std::setiosflags(std::ios_base::left) << std::setw(W) << std::setfill(' ')

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Mesh filenames
//------------------------------------------------------------------------------

const std::string mesh_dir = "../../meshes/";
// std::string default_mesh = mesh_dir + "Cubic-Cells/RF_fmt/gcube_2x2x2";
std::string default_mesh = mesh_dir + "Fracture-mesh/intersecting_fractures/unknown1";

// Max number of cells to plot a graph
constexpr size_t max_nb_cells_for_plot = 50000;

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
    ("solution,s", boost::program_options::value<int>()->default_value(11), "Select the solution")
    ("export-matrix,e", "Export matrix to Matrix Market format")
    ("boundary-conditions,b", boost::program_options::value<std::string>()->default_value("D"), "Boundary conditions (D for essential, N for natural)")
    ("mu", boost::program_options::value<double>()->default_value(1.), "Select the second Lamé coefficient")
    ("lambda", boost::program_options::value<double>()->default_value(1.), "Select the first Lamé coefficient")
    ("beta", boost::program_options::value<double>()->default_value(1.), "Select the Newton parameter")
    ("plot", boost::program_options::value<std::string>(), "Save plot of the solution to the given filename")
    // ("solver", boost::program_options::value<std::string>()->default_value("PardisoLU"), "Choice of solver, not case dependent. Options are: PardisoLU, UMFPACK, PaStiXLU, PaStiXLLT, EigenLU, EigenBiCGSTAB (reverts to EigenLU if the selected solver is not available)")
    ("solver", boost::program_options::value<std::string>()->default_value("PaStiXLU"), "Choice of solver, not case dependent. Options are: PardisoLU, UMFPACK, PaStiXLU, PaStiXLLT, EigenLU, EigenBiCGSTAB (reverts to EigenLU if the selected solver is not available)")
    ("solver-parameters,sparm",boost::program_options::value<std::vector<std::string>>()->multitoken(), SolverParamHelper)
    ("stabilization-parameter,x", boost::program_options::value<double>(), "Set the stabilization parameter");

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

  Elasticity::ForcingTermType f;
  Elasticity::DisplacementType u;
  Elasticity::FrictionType g;
  Elasticity::GradientDisplacementType grad_u;
  Elasticity::GradientDisplacementType1 sigma_u;
  
  double mu = vm["mu"].as<double>();
  double lambda = vm["lambda"].as<double>();
  double beta = vm["beta"].as<double>();

  std::cout << FORMAT(25) << "[main] Lamé coefficients: mu= " << mu << ", lambda= " << lambda << std::endl;

  switch (solution) {
  case 0:
    std::cout << "[main] Test case: trigonometric solution" << std::endl;
    f = trigonometric_f;
    u = trigonometric_u;
    grad_u = trigonometric_grad_u;
    sigma_u = trigonometric_sigma_u;
    break;

  case 1:
    std::cout << "[main] Test case: constant solution" << std::endl;
    f = constant_f;
    u = constant_u;
    grad_u = constant_grad_u;
    sigma_u = constant_sigma_u;
    g = constant_friction;
    break;

  case 2:
    std::cout << "[main] Test case: linear solution" << std::endl;
    f = linear_f;
    u = linear_u;
    grad_u = linear_grad_u;
    sigma_u = linear_sigma_u;
    break;
  
  case 3:
    std::cout << "[main] Test case: Fracture Test solution" << std::endl;
    f = fracture_test_f;
    u = fracture_test_u;
    grad_u = fracture_test_grad_u;
    sigma_u = fracture_test_sigma_u;
    break;

  case 4:
    std::cout << "[main] Test case: Continuous trignometric solution with 0 traction tangential component and 0 traction normal" << std::endl;
    f = trigonometric_tangential_f;
    u = trigonometric_tangential_u;
    grad_u = trigonometric_tangential_grad_u;
    sigma_u = trigonometric_tangential_sigma_u;
    break;
  
  case 5:
    std::cout << "[main] Test case: Continuous trignometric solution with 0 traction tangential component and negative traction normal" << std::endl;
    f = trigonometric_new_f;
    u = trigonometric_new_u;
    grad_u = trigonometric_new_grad_u;
    sigma_u = trigonometric_new_sigma_u;
    break;
  
  case 6:
    std::cout << "[main] Test case: Polynomial tangential solution" << std::endl;
    f = polynomial_new_f;
    u = polynomial_new_u;
    grad_u = polynomial_new_grad_u;
    sigma_u = polynomial_new_sigma_u;
    break;

  case 7:
    std::cout << "[main] Test case: Fracture Trigonometric solution with 0 traction" << std::endl;
    f = fracture_trigonometric_test_f;
    u = fracture_trigonometric_test_u;
    grad_u = fracture_trigonometric_test_grad_u;
    sigma_u = fracture_trigonometric_test_sigma_u;
    break;

  case 8:
    std::cout << "[main] Test case: Frictionless Contact Solution" << std::endl;
    f = frictionless_new_f;
    u = frictionless_new_u;
    grad_u = frictionless_new_grad_u;
    sigma_u = frictionless_new_sigma_u;
    g = frictionless_new_friction;
    break;

  case 9:
    std::cout << "[main] Test case: Friction Contact Solution" << std::endl;
    f = friction_f;
    u = friction_u;
    grad_u = friction_grad_u;
    sigma_u = friction_sigma_u;
    g = friction_friction;
    break;
  
  case 10:
    std::cout << "[main] Test case: Frictionless Incompressible Contact Solution" << std::endl;
    f = frictionless_incom_f;
    u = frictionless_incom_u;
    grad_u = frictionless_incom_grad_u;
    sigma_u = frictionless_incom_sigma_u;
    g = frictionless_incom_friction;
    break;

  case 11:
    std::cout << "[main] Test case: Friction Incompressible Contact Solution" << std::endl;
    f = friction_incom_f;
    u = friction_incom_u;
    grad_u = friction_incom_grad_u;
    sigma_u = friction_incom_sigma_u;
    g = friction_incom_friction;
    break;
  
  case 12:
    std::cout << "[main] Test case: Unknown continuous solution" << std::endl;
    f = unknown_f;
    u = unknown_u;
    grad_u = unknown_grad_u;
    sigma_u = unknown_sigma_u;
    g = unknown_friction;
    // need to choose E=4*10^9, \nu=0.2, \mu=1.666667e9, \lambda=1.111111e9, \beta=1e6
    break;

  default:
    std::cerr << "[main] ERROR: Unknown exact solution" << std::endl;
    exit(1);    
  }
  // Symmetric gradient
  Elasticity::GradientDisplacementType symgrad_u = [&grad_u](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d
                            { return symmetrise_matrix(grad_u(x,T,lambda)); };                          

  // Build the mesh and reorder based on the BCs
  MeshBuilder meshbuilder = MeshBuilder(mesh_file);
  std::unique_ptr<Mesh> mesh_ptr = meshbuilder.build_the_mesh();
  std::string bctype = vm["boundary-conditions"].as<std::string>();
  std::cout << FORMAT(25) << "[main] Boundary conditions: " << bctype << std::endl;

  // std::vector<size_t> fracture_making;
  // for(const auto& F:mesh_ptr->get_faces()){
  //   const double a1 = F->center_mass().x();
  //   const double a2 = F->center_mass().y();
  //   const double a3 = F->center_mass().z();
  //   if (abs(25*a1-2*a2+10*a3-11.4)<1e-6&& a1>0.2 &&a1<0.4 &&a2>0.3 && a2<0.8
  //     ) {     
  //       fracture_making.push_back(F->global_index());
  //   }else if (abs(250*a1+15*a2-75*a3-152)<1e-6 && a1>0.65 &&a1<0.8 &&a2>0.3 && a2<0.8
  //     ) {     
  //       fracture_making.push_back(F->global_index());
  //   }else if (abs((-0.1)*a1+(0.08)*a2+(0.4)*a3-0.174)<1e-6 && a1>0.1 &&a1<0.9 &&a2>0.3 && a2<0.8
  //     ) {     
  //       fracture_making.push_back(F->global_index());
  //   }else if (abs(a3-0.10)<1e-8 && a1>0.7 &&a1<0.9 &&a2>0.15 && a2<0.35
  //     ) {     
  //       fracture_making.push_back(F->global_index());
  //   } 
  // }


  // std::ofstream fcout("unknown1.frac");
  //       fcout << std::scientific << std::setprecision(17);
  //       fcout << "# *.frac file of 3D-mesh in REGN_FACE format\n";
  //       fcout << "# frac.txt\n";
  //       fcout << fracture_making.size()<<"\n";
    
  //       for (size_t i = 0; i < fracture_making.size(); ++i) {
  //           const Face * F=mesh_ptr->get_faces()[fracture_making[i]];
  //           fcout << i<<"  "<<F->n_vertices();
  //           for (size_t j = 0; j < F->n_vertices(); ++j) {
  //               fcout << "  " << F->get_vertices()[j]->global_index();
  //           }
  //           fcout << "\n";
  //       }
  //       fcout.close();

  ConnectedInformationGenerator CIG(mesh_file,*mesh_ptr.get());
  CIG.generate();

  std::vector<Face *> fracture_1_pointer;
  for(const auto & F: CIG.fracture_face_indices_pointer){
    const double a1 = F->center_mass().x();
    const double a2 = F->center_mass().y();
    const double a3 = F->center_mass().z();
    if (abs(25*a1-2*a2+10*a3-11.4)<1e-6&& a1>0.2 &&a1<0.4 &&a2>0.3 && a2<0.8
      ) {     
        fracture_1_pointer.push_back(F);
    }
  }

  BoundaryConditions BC(bctype, *mesh_ptr.get());
  BC.reorder_vertices("start");
  BC.reorder_edges("start");
  BC.reorder_faces("start");
  
  // Create DDR core
  bool use_threads = (vm.count("pthread") ? vm["pthread"].as<bool>() : true);
  std::cout << "[main] " << (use_threads ? "Parallel execution" : "Sequential execution") << std:: endl;
  DDRCore ddr_core(*mesh_ptr, K, use_threads);

  // Assemble the problem
  boost::timer::cpu_timer timer;
  timer.start();
  Elasticity elas(ddr_core, CIG, mu, lambda, beta, BC, use_threads);

  if(vm.count("stabilization-parameter")) {
    elas.stabilizationParameter() = vm["stabilization-parameter"].as<double>();
  }
  Eigen::VectorXd UDir = elas.vectorDirBC(u);
  Eigen::VectorXd U_Iguess = elas.vectorDirBC(u);
  Eigen::VectorXd normaldof(CIG.fracture_face_indices_pointer.size()), tangent1dof(CIG.fracture_face_indices_pointer.size()), tangent2dof(CIG.fracture_face_indices_pointer.size());
  for(size_t i=0; i<CIG.fracture_face_indices_pointer.size(); ++i){
        U_Iguess[elas.nDOFs()+3*i] = 0.00;
        normaldof[i] = elas.nDOFs()+3*i;
        tangent1dof[i] = elas.nDOFs()+3*i+1;
        tangent2dof[i] = elas.nDOFs()+3*i+2;
    }
  Eigen::VectorXd U_zero = Eigen::VectorXd::Zero(elas.tDOFs());
  // double error = 1;
  double tol = 1e-10; 
  size_t cI = 0;
  size_t mI = 50;
  Eigen::VectorXd uh;
  elas.assembleLinearSystem(f, g, u, sigma_u, symgrad_u, UDir, U_Iguess);
    timer.stop();
    double t_wall_model, t_proc_model;
    std::tie(t_wall_model, t_proc_model) = store_times(timer, "[main] Time model (wall/proc) ");

    // Export matrix if requested  
    if (vm.count("export-matrix")) {
      std::cout << "[main] Exporting matrix to Matrix Market format" << std::endl;
      saveMarket(elas.systemMatrix(), "A_elasticity_test.mtx");
      saveMarket(elas.systemVector(), "b_elasticity_test.mtx");
    }
  double error = (elas.systemVector()).norm();
  
  while (error > tol && cI<mI){
    double error1, error2, error3, error4;
    // if(error>tol){
      // Select solver and solve statically condensed system
      timer.start();
      std::string name_solver = vm["solver"].as<std::string>();
      LinearSolver<Elasticity::SystemMatrixType> solver(name_solver);
      std::cout << "[main] Solving the system using " << solver.name() << std::endl;
      Eigen::VectorXd u_gl = solver.compute_and_solve(elas.systemMatrix(), elas.systemVector());
      std::cout << "[main] residual linear solver: " 
                << solver.residual(elas.systemMatrix(), elas.systemVector(), u_gl) << std::endl;     
      // Re-create full vector of 
      uh = elas.reAssembleSolution(U_zero, u_gl);
      error1=uh.head(elas.nDOFs()).norm();
      error2=uh(normaldof).norm();
      error3=uh(tangent1dof).norm();
      error4=uh(tangent2dof).norm();
      uh = U_Iguess + uh;
      U_Iguess = uh;
      cI+=1;

      elas.assembleLinearSystem(f, g, u, sigma_u, symgrad_u, UDir, U_Iguess);
      timer.stop();
      double t_wall_model, t_proc_model;
      std::tie(t_wall_model, t_proc_model) = store_times(timer, "[main] Time model (wall/proc) ");

      // Export matrix if requested  
      if (vm.count("export-matrix")) {
        std::cout << "[main] Exporting matrix to Matrix Market format" << std::endl;
        saveMarket(elas.systemMatrix(), "A_elasticity_test.mtx");
        saveMarket(elas.systemVector(), "b_elasticity_test.mtx");
      }
    error = (elas.systemVector()).norm()/U_Iguess.norm();
    // }
    std::cout << "Number of Newton iterations: " << cI << ", Newton residual: " << error <<", Relative Error: " <<error1<< std::endl;
    std::cout << "Lagrange Normal: "<<error2 << ", Lagrange Tangent1: "<<error3 << ", Lagrange Tangent2: "<<error4 <<std::endl;
  }
  std::cout << "Number of Newton iterations: " << cI << ", final Newton residual: " << error << std::endl;
                               
  timer.stop();
  double t_wall_solve, t_proc_solve;
  std::tie(t_wall_solve, t_proc_solve) = store_times(timer, "[main] Time solve (wall/proc) ");  
  
  
  //----------------------------------------------
  //    Storing the jump and the z-coordinates
  //----------------------------------------------
  // std::cout<<elas.nDOFs()<<" "<<elas.tDOFs()<<"\n";
  // std::ostream & output = std::cout;
  // SerendipityProblem Ser(ddr_core, use_threads, output, false);
  // VSXGrad vsxgrad(Ser);
  // FractureInterpolator fracture_interpolator(vsxgrad,CIG,use_threads);
  // FractureDOFSpace fractureinfo(vsxgrad,CIG,fracture_interpolator,use_threads);
  // std::ofstream fcout("lambda_n_hex4.txt");
  // for(const auto F:CIG.fracture_face_indices_pointer){
  //   int face_index = CIG.findFaceIndex(*F);
  //   double lambda_n=uh(elas.nDOFs()+3*face_index);
  //   size_t iF = F->global_index();
  //   const double & xFy = F->center_mass()[1];
  //   const double & xFz = F->center_mass()[2];
  //   const Cell & T0 = *CIG.face_final_connected_pointer.at(mesh_ptr->face(iF))[0][0];
  //   int b=T0.face_orientation(T0.index_face(F));
  //   Eigen::VectorXd u_jump = b*(uh.segment(fractureinfo.globalOffset(*F),3)-uh.segment(fractureinfo.globalOffset(*F)+3,3));
  //   fcout<<xFy<<" "<<xFz<<" "<<lambda_n<<"\n";
  // }
  // fcout.close();
    
  //----------------------------------------------
  //    OUTPUTS
  //----------------------------------------------
  
  // Interpolate of exact solution and error vector
  Eigen::VectorXd uI1 = elas.fractureinfo().interpolate(u, lambda, MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
  Eigen::VectorXd uI2 = elas.fractureinfo().correction_term(u, lambda, uI1, MAX_DOE_FACE);  
  Eigen::VectorXd uI = uI1 + uI2;  
  Eigen::VectorXd eh = uh.head(elas.nDOFs()) - uI;        
   
  std::cout << "[main] Compute errors" << std::endl;
  // Errors in discrete norms
  std::cout << "[main] Discrete norms and errors:" << std::endl;
  std::vector<double> list_norms = elas.fractureinfo().computeH1norms(std::vector<Eigen::VectorXd> {eh, uI}, "symmetric");
  std::vector<double> list_l2_error_disp = elas.fractureinfo().computeL2normdisp(u,lambda,uh);
  std::vector<double> list_l2_error_dispgrad = elas.fractureinfo().computeL2normdispgrad(grad_u,lambda,uh);
  std::vector<double> list_l2_error_jumpdisp = elas.fractureinfo().computeL2normdispjump(u,lambda,uh);
  std::vector<double> list_l2_error_traction = elas.fractureinfo().computeL2normtraction(sigma_u,lambda,uh);
  std::vector<Eigen::VectorXd> list_contact = elas.fractureinfo().computecontact(g,uh,beta,lambda);
  Eigen::VectorXd discrete_normal_traction = elas.fractureinfo().computenormaltraction(fracture_1_pointer,uh);
  std::vector<double> list_reconstruct_norm = elas.fractureinfo().reconstruct_norm(uh);
  double error_symH1error = (list_norms[1] < 1e-12 ? list_norms[0] : list_norms[0]/list_norms[1]);
  double relative_l2_error_disp = (list_l2_error_disp[1] < 1e-12 ? list_l2_error_disp[0] : list_l2_error_disp[0]/list_l2_error_disp[1]);
  double relative_l2_error_dispgrad = (list_l2_error_dispgrad[1] < 1e-12 ? list_l2_error_dispgrad[0] : list_l2_error_dispgrad[0]/list_l2_error_dispgrad[1]);
  double relative_l2_error_jumpdisp = (list_l2_error_jumpdisp[1] < 1e-12 ? list_l2_error_jumpdisp[0] : list_l2_error_jumpdisp[0]/list_l2_error_jumpdisp[1]);
  double relative_l2_error_traction = (list_l2_error_traction[1] < 1e-12 ? list_l2_error_traction[0] : list_l2_error_traction[0]/list_l2_error_traction[1]);
  std::cout << FORMAT(38) << "\tSymmetric H1 relative error: " << error_symH1error << " (norm of I_h u: "<< list_norms[1] << ")" << std::endl;
  std::cout << FORMAT(38) << "\tSymmetric H1 absolute error: " << list_norms[0] << std::endl;
  std::cout << FORMAT(38) << "\tRelative L2 error Displacement: " << relative_l2_error_disp << std::endl;
  std::cout << FORMAT(38) << "\tL2 error Displacement: " << list_l2_error_disp[0] << std::endl;
  std::cout << FORMAT(38) << "\tRelative L2 error Gradient: " << relative_l2_error_dispgrad << std::endl;
  std::cout << FORMAT(38) << "\tL2 error Gradient: " << list_l2_error_dispgrad[0] << std::endl;
  std::cout << FORMAT(38) << "\tRelative L2 error Jump Displacement: " << relative_l2_error_jumpdisp << std::endl;
  std::cout << FORMAT(38) << "\tL2 error Jump Displacement: " << list_l2_error_jumpdisp[0] << std::endl;
  std::cout << FORMAT(38) << "\tRelative L2 error Traction: " << relative_l2_error_traction << std::endl;
  std::cout << FORMAT(38) << "\tL2 error Traction: " << list_l2_error_traction[0] << std::endl;

  std::cout << "[main] Mesh diameter " << mesh_ptr->h_max() << std::endl;

  // double a=0;
  // double b=0;
  // double c=0;
  // for (size_t i=0; i<fracture_1_pointer.size();i++){
  //    Face * F = fracture_1_pointer[i];
  //    a+=F->measure()*discrete_normal_traction[i];
  //    b+=discrete_normal_traction[i];
  //    c+=F->measure();
  // }
  // std::cout<<a<<" "<<b<<" "<<c<<"\n";

  // --------------------------------------------------------------------------
  //           Creates a .vtu file of the displacement
  // --------------------------------------------------------------------------
  // Only if we do not have too many cells
  if (vm.count("plot") && mesh_ptr->n_cells() <= max_nb_cells_for_plot) {
    std::cout << "[main] Writing solution to file" << std::endl;
    std::string filename = vm["plot"].as<std::string>();
    VtuWriter plotdata(mesh_ptr.get());
		
		// Exact displacement
    std::vector<VectorRd> exact_u_vertex(mesh_ptr->n_vertices(), VectorRd::Zero());
    for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
      VectorRd v_coords = mesh_ptr->vertex(iV)->coords();
      const int count = CIG.vertex_final_connected_size.at(mesh_ptr->vertex(iV));
      for(int i=0; i<count; ++i){
        const Cell * T = CIG.vertex_final_connected_pointer.at(mesh_ptr->vertex(iV))[i][0];
        exact_u_vertex[iV] += u(v_coords, T, lambda);
      }
      exact_u_vertex[iV] /= count;
    }
		
    // Approximate pressure and velocity at the vertices
    std::vector<VectorRd> u_vertex = elas.fractureinfo().computeVertexValues(uh);
        
    // Plot
    plotdata.write_to_vtu(filename + "_u.vtu", 
                          std::vector<std::vector<VectorRd>> {u_vertex, exact_u_vertex},
                          std::vector<std::string>{"displacement", "exact_displacement"});                             
  }   

  FaceVtuWriter writer(CIG.fracture_face_indices_pointer);
  writer.write_to_vtu("jump_u_normal.vtu", list_contact[0], "my_solution");
  writer.write_to_vtu("jump_u_tangent1.vtu", list_contact[1], "my_solution");
  writer.write_to_vtu("jump_u_tangent2.vtu", list_contact[2], "my_solution");
  writer.write_to_vtu("contact_condition.vtu", list_contact[5], "my_solution");
  writer.write_to_vtu("traction_tangential.vtu", list_contact[6], "my_solution");

  FaceVtuWriter writer1(fracture_1_pointer);
  writer1.write_to_vtu("traction_normal.vtu", discrete_normal_traction, "Tn");


  // --------------------------------------------------------------------------
  //                     Creates .txt file with data and results
  // --------------------------------------------------------------------------
  std::ofstream out("results.txt");
  out << std::scientific;
  out << "Solution: " << solution << std::endl;
  out << "Boundary conditions: " << bctype << std::endl;
  out << "mu: " << mu << std::endl;
  out << "lambda: " << lambda << std::endl;
  out << "Newton parameter: " << beta << std::endl;
  out << "StabilisationParameter: " << elas.stabilizationParameter() << std::endl;
  out << "Mesh: " << mesh_file << std::endl;
  out << "Degree: " << K << std::endl;
  out << "Nb_Unknown_DOF: " << elas.nUKN() << std::endl; 
  out << "MeshSize: " << mesh_ptr->h_max() << std::endl;
  out << "NbCells: " << mesh_ptr->n_cells() << std::endl;
  out << "NbFaces: " << mesh_ptr->n_faces() << std::endl;
  out << "NbEdges: " << mesh_ptr->n_edges() << std::endl;
  out << "NbNewtonIterations: " << cI << std::endl;
  out << "DimVSXGrad: " << elas.vsxGrad().dimension() << std::endl;
  out << "SizeSCsystem: " << elas.systemMatrix().rows() << std::endl;
  // Discrete errors and norms
  out << "E_symH1Dis: " << error_symH1error << std::endl;
  out << "absE_symH1Dis: " << list_norms[0] << std::endl;
  out << "relative_l2_error_displacement: " << relative_l2_error_disp <<std::endl;
  out << "l2_error_displacement: " << list_l2_error_disp[0] <<std::endl;
  out << "relative_grad_error_displacement: " << relative_l2_error_dispgrad <<std::endl;
  out << "grad_error_displacement: " << list_l2_error_dispgrad[0] <<std::endl;
  out << "relative_l2_error_jump_displacement: " << relative_l2_error_jumpdisp <<std::endl;
  out << "l2_error_jump_displacement: " << list_l2_error_jumpdisp[0] <<std::endl;
  out << "relative_l2_error_traction: " << relative_l2_error_traction <<std::endl;
  out << "l2_error_traction: " << list_l2_error_traction[0] <<std::endl;
  out << "Potential_normal_Fracture_1: " << std::setprecision(12) << list_reconstruct_norm[0] << std::endl;
  out << "Potential_normal_Fracture_2: " << std::setprecision(12) << list_reconstruct_norm[1] << std::endl;
  out << "Potential_normal_Fracture_3: " << std::setprecision(12) << list_reconstruct_norm[2] << std::endl;
  out << "Potential_normal_Fracture_4: " << std::setprecision(12) << list_reconstruct_norm[3] << std::endl;
  out << "Potential_tangential_Fracture_1: " << std::setprecision(12) << list_reconstruct_norm[4] << std::endl;
  out << "Potential_tangential_Fracture_2: " << std::setprecision(12) << list_reconstruct_norm[5] << std::endl;
  out << "Potential_tangential_Fracture_3: " << std::setprecision(12) << list_reconstruct_norm[6] << std::endl;
  out << "Potential_tangential_Fracture_4: " << std::setprecision(12) << list_reconstruct_norm[7] << std::endl;
  out << "Traction_normal_Fracture_1: " << std::setprecision(12) << list_reconstruct_norm[8] << std::endl;
  out << "Traction_normal_Fracture_2: " << std::setprecision(12) << list_reconstruct_norm[9] << std::endl;
  out << "Traction_normal_Fracture_3: " << std::setprecision(12) << list_reconstruct_norm[10] << std::endl;
  out << "Traction_normal_Fracture_4: " << std::setprecision(12) << list_reconstruct_norm[11] << std::endl;
  out << "Traction_tangential_Fracture_1: " << std::setprecision(12) << list_reconstruct_norm[12] << std::endl;
  out << "Traction_tangential_Fracture_2: " << std::setprecision(12) << list_reconstruct_norm[13] << std::endl;
  out << "Traction_tangential_Fracture_3: " << std::setprecision(12) << list_reconstruct_norm[14] << std::endl;
  out << "Traction_tangential_Fracture_4: " << std::setprecision(12) << list_reconstruct_norm[15] << std::endl;
  // out << "lambda_error: " << lambda_error << std::endl;
  // out << "TwallModel: " << t_wall_model << std::endl;  
  // out << "TprocModel: " << t_proc_model << std::endl;  
  out << "TwallSolve: " << t_wall_solve << std::endl;  
  out << "TprocSolve: " << t_proc_solve << std::endl;
  out << std::flush;
  out.close();

  std::cout << "[main] Done" << std::endl;
  return 0;
}

//------------------------------------------------------------------------------
// Elasticity
//------------------------------------------------------------------------------

Elasticity::Elasticity(
               const DDRCore & ddrcore,
               const ConnectedInformationGenerator & connectedinfo,
               const double & mu,
               const double & lambda,
               const double & beta,
               const BoundaryConditions & BC,
               bool use_threads,
               std::ostream & output
               )
  : m_ddrcore(ddrcore),
    m_connectedinfo(connectedinfo),
    m_use_threads(use_threads),
    m_output(output),
    m_ser_pro(ddrcore, use_threads, output, false),
    m_vsxgrad(m_ser_pro),
    m_fracture_interpolator(m_vsxgrad,m_connectedinfo,use_threads),
    m_fractureinfo(m_vsxgrad,m_connectedinfo,m_fracture_interpolator,use_threads),
    m_nloc_sc(0*m_ser_pro.nDOFs_cells_SXGrad() * dimspace),
    m_mu(mu),
    m_lambda(lambda),
    m_beta(beta),
    m_BC(BC),
    m_nDOFs_dir_vertices(0),
    m_nDOFs_dir_edges(0),
    m_nDOFs_dir_faces(0),
    m_locSCUKN({}),
    m_locGLUKN({}),
    m_DOFtoSCUKN(-2 * Eigen::VectorXi::Ones(tDOFs())),  // These maps are initialised at -2, but created just below
    m_DOFtoGLUKN(-2 * Eigen::VectorXi::Ones(tDOFs())),
    //m_stab_par(2*mu+lambda),
    m_stab_par(0.1),
    m_A(nGLUKN(), nGLUKN()),
    // m_fA(nGLUKN(),fDOFs()),
    m_b(Eigen::VectorXd::Zero(nGLUKN())),
    m_sc_A(nSCUKN(), nGLUKN()),
    m_sc_b(Eigen::VectorXd::Zero(nSCUKN()))
{
  m_output << "[Elasticity] Initializing" << std::endl;
  // To remove the static condensation, multiply by 0 the vectors in the initialisation of m_nloc_sc_u and m_nloc_sc_p.

  // Compute number of Dirichlet DOFs for displacement (vertices, edges, faces)
  for (Vertex * V : ddrcore.mesh().get_b_vertices()){
    if (BC.type(*V) == "dir"){
      m_nDOFs_dir_vertices += m_fractureinfo.numLocalDofsVertex(*V);
    }
  }
  for (Edge * E : ddrcore.mesh().get_b_edges()){
    if (BC.type(*E) == "dir"){
      m_nDOFs_dir_edges += m_fractureinfo.numLocalDofsEdge(*E);
    }
  }
  for (Face * F : ddrcore.mesh().get_b_faces()){
    if (BC.type(*F) == "dir"){
      m_nDOFs_dir_faces += m_fractureinfo.numLocalDofsFace(*F);
    }
  }


  // We compute the maps DOF->SCUKN and DOF->GLUKN by locating the positions of these unknowns through m_locSCUKN, m_locGLUKN,
  // and using replaceSectionsVector (see BChandlers). 
  // Following the construction of the m_locXXX vectors below is best done by doodling a picture of a global vector of DOFs
  // and locating the dirichlet vertex/edge/face DOFs in this global vector
  const Eigen::VectorXi MinusOnes = -Eigen::VectorXi::Ones(tDOFs());

  // ----- DOF->SCUKN -------
  m_locSCUKN.resize(1);
  // Statically condensed cell unknowns: start after all vertices, edges and faces
  m_locSCUKN[0].first = m_fractureinfo.nDOFs_vertices() + m_fractureinfo.nDOFs_edges() + m_fractureinfo.nDOFs_faces();
  m_locSCUKN[0].second = nSCUKN();
  // Create map
  size_t s = nSCUKN();
  m_DOFtoSCUKN = replaceSectionsVector(MinusOnes, Eigen::VectorXi::LinSpaced(s, 0, s-1).eval(), m_locSCUKN);

  // ----- DOF->GLUKN -------
  m_locGLUKN.resize(3);
  // non-Dirichlet vertex DOFs
  m_locGLUKN[0].first = m_nDOFs_dir_vertices;
  m_locGLUKN[0].second = m_fractureinfo.nDOFs_vertices() - m_nDOFs_dir_vertices;

  // non-Dirichlet edge DOFs 
  m_locGLUKN[1].first = m_fractureinfo.nDOFs_vertices() + m_nDOFs_dir_edges;
  m_locGLUKN[1].second = m_fractureinfo.nDOFs_edges() - m_nDOFs_dir_edges; 

  // non-Dirichlet face DOFs and non-SC cell unknowns
  m_locGLUKN[2].first = m_fractureinfo.nDOFs_vertices() + m_fractureinfo.nDOFs_edges() + m_nDOFs_dir_faces;
  m_locGLUKN[2].second = m_fractureinfo.nDOFs_faces() - m_nDOFs_dir_faces + (m_fractureinfo.nDOFs_cells() + m_fractureinfo.fracture_dofs() - nSCUKN()); 

  s = nGLUKN();
  m_DOFtoGLUKN = replaceSectionsVector(MinusOnes, Eigen::VectorXi::LinSpaced(s, 0, s-1).eval(), m_locGLUKN);
}

//------------------------------------------------------------------------------
Eigen::VectorXd Elasticity::vectorDirBC( const DisplacementType & u )
{
  // We interpolate u then put them in vectorDirBC (which has otherwise only zeros)
  // It's not the most efficient, but it is the easiest and this function is only called once anyway
  Eigen::VectorXd interp = m_fractureinfo.interpolate(u, m_lambda, MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
  Eigen::VectorXd vectorBC = Eigen::VectorXd::Zero(tDOFs());
    
  for (Vertex * V : m_ddrcore.mesh().get_b_vertices()){
    if (m_BC.type(*V) == "dir"){
      vectorBC.segment(m_fractureinfo.globalOffset(*V), m_fractureinfo.numLocalDofsVertex(*V))
          = interp.segment(m_fractureinfo.globalOffset(*V), m_fractureinfo.numLocalDofsVertex(*V));
    }
  }

  for (Edge * E : m_ddrcore.mesh().get_b_edges()){
    if (m_BC.type(*E) == "dir"){
      vectorBC.segment(m_fractureinfo.globalOffset(*E), m_fractureinfo.numLocalDofsEdge(*E))
          = interp.segment(m_fractureinfo.globalOffset(*E), m_fractureinfo.numLocalDofsEdge(*E));
    }
  }

  for (Face * F : m_ddrcore.mesh().get_b_faces()){
    if (m_BC.type(*F) == "dir"){
      vectorBC.segment(m_fractureinfo.globalOffset(*F), m_fractureinfo.numLocalDofsFace(*F))
          = interp.segment(m_fractureinfo.globalOffset(*F), m_fractureinfo.numLocalDofsFace(*F));
    }
  }
  return vectorBC;
}

//------------------------------------------------------------------------------

void Elasticity::assembleLinearSystem(
              const ForcingTermType & f, 
              const FrictionType &g,
              const DisplacementType & u,
              const GradientDisplacementType1 & sigma_u, 
              const GradientDisplacementType & symgrad_u, 
              const Eigen::VectorXd & UDir,
              const Eigen::VectorXd & U_Iguess )
{
  // Assemble all local contributions
  auto assemble_all = [this, &f, &g, &u, &sigma_u, &symgrad_u, &UDir, &U_Iguess](
                             size_t start,
                             size_t end,
                             std::list<Eigen::Triplet<double> > * triplets_gl,
                             Eigen::VectorXd * rhs_gl,
                             std::list<Eigen::Triplet<double> > * triplets_sc,
                             Eigen::VectorXd * rhs_sc
                             )->void
                      {
                        for (size_t iT = start; iT < end; iT++) {
                          this->_assemble_local_contribution(
                                                iT,
                                                this->_compute_local_contribution(iT, f, g, u, sigma_u, symgrad_u, UDir, U_Iguess),
                                                *triplets_gl,
                                                *rhs_gl,
                                                *triplets_sc,
                                                *rhs_sc
                                                );
                        } // for iT
                      };
                      
  // Assemble the matrix and rhs
  if (m_use_threads) {
    m_output << "[Elasticity] Parallel assembly" << std::endl;
  }else{
    m_output << "[Elasticity] Sequential assembly" << std::endl;
  }
  std::tie(m_A, m_b, m_sc_A, m_sc_b) = parallel_assembly_system(m_ddrcore.mesh().n_cells(), this->nGLUKN(), std::make_pair(this->nSCUKN(), this->nGLUKN()), this->nSCUKN(), assemble_all, m_use_threads);
}

//------------------------------------------------------------------------------
 
std::pair<Eigen::MatrixXd, Eigen::VectorXd>
    Elasticity::_compute_local_contribution (
                               const size_t iT, 
                               const ForcingTermType & f,
                               const FrictionType & g,
                               const DisplacementType & u,
                               const GradientDisplacementType1 & sigma_u,
                               const GradientDisplacementType & symgrad_u,
                               const Eigen::VectorXd & UDir,
                               const Eigen::VectorXd & U_Iguess
                               ) const
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);
  const Cell * cell_ptr=m_ddrcore.mesh().cell(iT);

  size_t dim_u = m_vsxgrad.dimensionCell(iT);
  size_t dim_fracture_T = m_fractureinfo.cell_fracture_dofs()(iT);
  size_t dim_T = dim_u + dim_fracture_T;
  
  Eigen::MatrixXd AT = Eigen::MatrixXd::Zero(dim_T, dim_T);
  Eigen::VectorXd lT = Eigen::VectorXd::Zero(dim_T);

  //------------------------------------------------------------------------------
  // Local matrix
  //------------------------------------------------------------------------------

  // Mass matrices for (P^k(T))^(3x3) and P^k
  MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_ddrcore.degree());
  Eigen::MatrixXd mass_Pk3x3_T = GramMatrix(T, *m_vsxgrad.cellBases(iT).Polyk3x3, int_mono_2k);
  Eigen::MatrixXd mass_Pk_T = GramMatrix(T, *m_vsxgrad.sXgrad().cellBases(iT).Polyk, int_mono_2k);

  // Contribution 2mu eps(u) : eps(v), including stabilisation scaled by mu
  Eigen::MatrixXd GsT = m_vsxgrad.cellSymGradient(iT);
   AT.topLeftCorner(dim_u, dim_u) += 2 * m_mu * GsT.transpose() * mass_Pk3x3_T * GsT + m_mu * m_vsxgrad.computeStabilisation(iT);
  // AT.topLeftCorner(dim_u, dim_u) += 2 * m_mu * GsT.transpose() * mass_Pk3x3_T * GsT + T.diam()*m_stab_par * m_vsxgrad.computeStabilisation(iT);

  // Contribution lambda div(u) div(v)
  Eigen::MatrixXd divT = m_vsxgrad.cellDivergence(iT);
  AT.topLeftCorner(dim_u, dim_u) += m_lambda * divT.transpose() * mass_Pk_T * divT;

  // Contribution coming from the lagrange multiplier
  for(size_t i=0;i<m_connectedinfo.Cell_fractures.at(m_ddrcore.mesh().cell(iT)).first;++i){
    const Face *pF=m_connectedinfo.Cell_fractures.at(m_ddrcore.mesh().cell(iT)).second[i];
    const Face & F = *pF;
    size_t iF = F.global_index();
    double g_face = g(pF,m_lambda);
    Eigen::MatrixXd b_fracture;
    Eigen::MatrixXd bF = Eigen::MatrixXd::Zero(3, m_vsxgrad.dimensionFace(F));
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    const Eigen::Vector3d & xF = F.center_mass();
    int a=T.face_orientation(T.index_face(&F));
    auto frame = m_fractureinfo.get_frame(F);
    std::array<Eigen::Vector3d, 3> dirs = { frame.normal, frame.tangent1, frame.tangent2 };
    for (int j = 0; j < 3; ++j) {
        for (int l = 0; l < 3; ++l) {
            bF(l, m_vsxgrad.localOffset(F) + j) += 
                a * F.measure() * scalar_product(
                    dirs[l],
                    this->m_vsxgrad.faceBases(iF).Polykmo3->function(j, xF)
                );
        }
    }                                    
    b_fracture = m_vsxgrad.extendOperator(T, F, bF).transpose();
    AT.block(0,dim_u+3*i,dim_u,3) += b_fracture;
    AT.block(dim_u+3*i,dim_u+3*i,3,3) += (0.5)*I;

    const Cell & T0 = *m_connectedinfo.face_final_connected_pointer.at(m_ddrcore.mesh().face(iF))[0][0];
    int b=T0.face_orientation(T0.index_face(&F));
    Eigen::VectorXd u_jump = b*(U_Iguess.segment(m_fractureinfo.globalOffset(F),3)-U_Iguess.segment(m_fractureinfo.globalOffset(F)+3,3));
    double u_jump_normal = m_vsxgrad.sXgrad().faceBases(iF).Polykmo->function(0, xF) * scalar_product(u_jump,dirs[0]);
    double u_jump_tangent1 = m_vsxgrad.sXgrad().faceBases(iF).Polykmo->function(0, xF) * scalar_product(u_jump,dirs[1]);
    double u_jump_tangent2 = m_vsxgrad.sXgrad().faceBases(iF).Polykmo->function(0, xF) * scalar_product(u_jump,dirs[2]);
    double lambda_prev_normal = U_Iguess(nDOFs() + 3 * m_connectedinfo.findFaceIndex(F));
    double lambda_prev_tangent1 = U_Iguess(nDOFs() + 3 * m_connectedinfo.findFaceIndex(F)+1);
    double lambda_prev_tangent2 = U_Iguess(nDOFs() + 3 * m_connectedinfo.findFaceIndex(F)+2);
    if(lambda_prev_normal + m_beta*u_jump_normal >= 0){
      AT.block(dim_u+3*i,m_vsxgrad.localOffset(T,F),1,3) -= a * m_vsxgrad.sXgrad().faceBases(iF).Polykmo->function(0, xF) * 
                          m_beta * dirs[0].transpose();
      AT(dim_u+3*i,dim_u+3*i) -= 0.5;
      lT(dim_u+3*i) += (0.5)*m_beta*(u_jump_normal);
    }else{
      lT(dim_u+3*i) -= (0.5)*lambda_prev_normal;
    }
    double tangent_part_norm = sqrt(std::pow(lambda_prev_tangent1 + m_beta*u_jump_tangent1, 2) + std::pow(lambda_prev_tangent2 + m_beta*u_jump_tangent2, 2));
    if(tangent_part_norm <= g_face){
      AT.block(dim_u+3*i+1,m_vsxgrad.localOffset(T,F),1,3) -= a*m_vsxgrad.sXgrad().faceBases(iF).Polykmo->function(0, xF) * 
                          m_beta * dirs[1].transpose();
      AT.block(dim_u+3*i+2,m_vsxgrad.localOffset(T,F),1,3) -= a*m_vsxgrad.sXgrad().faceBases(iF).Polykmo->function(0, xF) * 
                          m_beta * dirs[2].transpose();
      AT(dim_u+3*i+1,dim_u+3*i+1) -= 0.5;
      AT(dim_u+3*i+2,dim_u+3*i+2) -= 0.5;
      lT(dim_u+3*i+1) += (0.5)*m_beta*(u_jump_tangent1);
      lT(dim_u+3*i+2) += (0.5)*m_beta*(u_jump_tangent2);
    }else{
      lT(dim_u+3*i+1) -= (0.5)*(lambda_prev_tangent1-(g_face/tangent_part_norm)*(lambda_prev_tangent1 + m_beta*u_jump_tangent1));
      lT(dim_u+3*i+2) -= (0.5)*(lambda_prev_tangent2-(g_face/tangent_part_norm)*(lambda_prev_tangent2 + m_beta*u_jump_tangent2));
      double coeff1 = lambda_prev_tangent1 + m_beta*u_jump_tangent1;
      double coeff2 = lambda_prev_tangent2 + m_beta*u_jump_tangent2;
      Eigen::VectorXd grad_coeff1(5), grad_coeff2(5);
      grad_coeff1.head(3) = m_vsxgrad.sXgrad().faceBases(iF).Polykmo->function(0, xF) * 
                          m_beta * dirs[1].transpose();
      grad_coeff2.head(3) = m_vsxgrad.sXgrad().faceBases(iF).Polykmo->function(0, xF) * 
                          m_beta * dirs[2].transpose();
      grad_coeff1(3)=1;
      grad_coeff1(4)=0;
      grad_coeff2(3)=0;
      grad_coeff2(4)=1;
      auto res = m_fractureinfo.computeDf(dirs[1], dirs[2], coeff1, coeff2, grad_coeff1, grad_coeff2);
      AT.block(dim_u+3*i+1,m_vsxgrad.localOffset(T,F),1,3) -= a*g_face*res.co1.head(3).transpose();
      AT.block(dim_u+3*i+2,m_vsxgrad.localOffset(T,F),1,3) -= a*g_face*res.co2.head(3).transpose();
      AT.block(dim_u+3*i+1,dim_u+3*i+1,1,2) -= (0.5)*g_face*res.co1.tail(2).transpose();
      AT.block(dim_u+3*i+2,dim_u+3*i+1,1,2) -= (0.5)*g_face*res.co2.tail(2).transpose();
    }
  }
    
  //------------------------------------------------------------------------------
  // Local source vector
  //------------------------------------------------------------------------------
  // Forcing term
  QuadratureRule quad_2kpt_T = generate_quadrature_rule(T, 2 * m_vsxgrad.degree() + 2);
  std::function<Eigen::Vector3d(const Eigen::Vector3d&)> f_mu_lambda = [&f, this, &cell_ptr](const Eigen::Vector3d & x)->Eigen::Vector3d
                  { return f(this->m_mu,this->m_lambda,x,cell_ptr);};
  if(m_ddrcore.degree()==0){
    lT += m_vsxgrad.cellPotential(iT).transpose()
      * integrate(f_mu_lambda,
                  evaluate_quad<Function>::compute(*(this->m_vsxgrad.cellBases(iT).Polykpo3), quad_2kpt_T),
                  quad_2kpt_T);
    }
  else{
    lT.segment(m_vsxgrad.localOffset(T),m_vsxgrad.numLocalDofsCell(iT)) += integrate(f_mu_lambda,
                  evaluate_quad<Function>::compute(*(this->m_vsxgrad.cellBases(iT).Polykmo3), quad_2kpt_T),
                  quad_2kpt_T);
  }

  // //Taking into consideration the discontinuity                
  // for(size_t i=0;i<m_connectedinfo.Cell_fractures.at(m_ddrcore.mesh().cell(iT)).second.size();++i){
  //   const Face *pF=m_connectedinfo.Cell_fractures.at(m_ddrcore.mesh().cell(iT)).second[i];
  //   const Face & F = *pF;
  //   size_t iF = F.global_index();
  //   int dqr_face=2*m_ddrcore.degree()+1;
  //   const Cell * T_correct = m_connectedinfo.face_final_connected_pointer.at(m_ddrcore.mesh().face(F.global_index()))[0][0];
  //   const Cell * T_other = m_connectedinfo.face_final_connected_pointer.at(m_ddrcore.mesh().face(F.global_index()))[1][0];
  //   QuadratureRule quad_dqrbc_F = generate_quadrature_rule(F, dqr_face);
  //       Eigen::Vector3d nF = T_correct->face_normal(T_correct->index_face(&F));    
  //       std::function<Eigen::Vector3d(const Eigen::Vector3d&)> sigma_mu_lambda_jump = [&sigma_u, &nF, this, &T_correct, &T_other](const Eigen::Vector3d & x)->Eigen::Vector3d
  //                     { return (sigma_u(this->m_mu,this->m_lambda,x,T_correct)-sigma_u(this->m_mu,this->m_lambda,x,T_other))*nF;};        
  //       int a=T.face_orientation(T.index_face(&F));
  //       for(int j = 0; j < 3; ++j){
  //         Eigen::Vector3d colj = m_fractureinfo.lagrange_dof().at(m_ddrcore.mesh().face(iF)).col(j);
  //         std::function<double(const Eigen::Vector3d&)> u_local = [&u, colj, &cell_ptr](const Eigen::Vector3d & x) -> double {
  //             return scalar_product(u(x,cell_ptr), colj);
  //         };
  //         for (size_t iqn = 0; iqn < quad_dqrbc_F.size(); iqn++)
  //         {
  //           double qr_weight = quad_dqrbc_F[iqn].w;
  //           lT(dim_u+3*i+j) += a * qr_weight * u_local(quad_dqrbc_F[iqn].vector());
  //         }
  //       }        
  //       Eigen::RowVectorXd bF =  
  //                 (0.5)*integrate(sigma_mu_lambda_jump, 
  //                           evaluate_quad<Function>::compute(*m_vsxgrad.faceBases(F).Polykpo3, quad_dqrbc_F),
  //                           quad_dqrbc_F).transpose() 
  //                         * m_vsxgrad.facePotential(F);
  //       lT.head(dim_u) += m_vsxgrad.extendOperator(T, F, bF).transpose();                 
  // }
  // Dirichlet BCs
  // lT.head(dim_u) -= AT.topLeftCorner(dim_u, dim_u) * m_fractureinfo.restrict(T, UDir);
  // lT.head(dim_T) -= AT.topLeftCorner(dim_T, dim_u) * m_fractureinfo.restrict(T, UDir);
  lT.head(dim_u) -= AT.topLeftCorner(dim_u, dim_T) * m_fractureinfo.restrict(T, U_Iguess);
  
  // // Neumann BCs
  // if (T.is_boundary()){
  //   const size_t dqrbc = MAX_DOE_FACE;
  //   for (Face * pF : T.get_faces()){
  //     const Face & F = *pF;
  //     if (m_BC.type(F) == "neu"){
  //       // Unit normal vector to F pointing out of the domain
  //       const Cell & TF = *F.cell(0);
  //       Eigen::Vector3d nF = TF.face_normal(TF.index_face(&F));
        
  //       QuadratureRule quad_dqrbc_F = generate_quadrature_rule(F, dqrbc);
        
  //       // Traction boundary conditions: int_F C eps(u)n . P_F v
  //       std::function<Eigen::Vector3d(const Eigen::Vector3d &)> 
  //           Cepsu_nF = [&symgrad_u, &nF, this,&cell_ptr](const Eigen::Vector3d & x)->Eigen::Vector3d {
  //                        return 2 * this->m_mu * symgrad_u(x,cell_ptr,m_lambda) * nF + this->m_lambda * symgrad_u(x,cell_ptr,m_lambda).trace() * nF;
  //                        };
  //       // Create integral of Ceps(u)nF * PotF as a linear (row) operator on the face DOFs, to extended it afterwards to the cell DOFs.
  //       Eigen::RowVectorXd bF =  
  //             integrate(Cepsu_nF, 
  //                       evaluate_quad<Function>::compute(*m_vsxgrad.faceBases(F).Polykpo3, quad_dqrbc_F),
  //                       quad_dqrbc_F).transpose() 
  //                     * m_vsxgrad.facePotential(F);
                                 
  //       // Add contribution
  //       lT += m_vsxgrad.extendOperator(T, F, bF).transpose();
        
  //     } // if m_BC.type(F)="neu"
  //   } // Loop over F
  // } // if T.is_boundary()

  return std::make_pair(AT, lT);
}

//------------------------------------------------------------------------------

void Elasticity::_assemble_local_contribution(
                                          size_t iT,
                                          const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT,
                                          std::list<Eigen::Triplet<double> > & triplets_gl,
                                          Eigen::VectorXd & rhs_gl,
                                          // std::list<Eigen::Triplet<double> > & triplets_bdryMat,
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
  // DOFs_gl contain the globally coupled unknowns. We check if they are actual unknowns in the globally coupled system through DOFtoGLUKN
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

LocalStaticCondensation Elasticity::_compute_static_condensation(const size_t & iT) const
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  // Dimensions
  size_t dim_sc = m_nloc_sc(iT);      // nb of SC dofs
  // size_t dim_dofs = m_vsxgrad.dimensionCell(iT) - dim_sc;     // number of local displacement unknowns after SC
  size_t dim_frac_T = m_fractureinfo.cell_fracture_dofs()(iT);
  size_t dim_dofs = m_vsxgrad.dimensionCell(iT) + dim_frac_T - dim_sc; 

  // Creation of permutation matrix: trivial here since SC dofs are already at the end (cell unknowns)
  Eigen::MatrixXd Perm = Eigen::MatrixXd::Identity(dim_dofs+dim_sc, dim_dofs+dim_sc);

  // Creation of global DOFs for system
  //   DOFs_gl and DOFs_sc respectively contains the DOFs of u for globally coupled system and for the recovery system
  //   These DOFs need to be read through DOFtoGLUKN and DOFtoSCUKN when assembling the systems
  std::vector<size_t> DOFs_gl(dim_dofs, 0);
  std::vector<size_t> DOFs_sc(dim_sc, 0);
  auto DOFs_T_vsxgrad = m_fractureinfo.globalDOFIndices(T);
   
  std::copy(DOFs_T_vsxgrad.begin(), DOFs_T_vsxgrad.begin()+dim_dofs, DOFs_gl.begin()); // put skeletal DOFs in DOFs_gl
   
  if (dim_sc>0){
    std::copy(DOFs_T_vsxgrad.begin()+dim_dofs, DOFs_T_vsxgrad.end(), DOFs_sc.begin());
  }  

  return LocalStaticCondensation(Perm, DOFs_gl, DOFs_sc);
}

