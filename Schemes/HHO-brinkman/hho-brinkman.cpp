// Author: Jerome Droniou (jerome.droniou@monash.edu)
#include <fstream>
#include <iomanip>
#include <thread>

#include "hho-brinkman.hpp"
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <max_degrees_quadratures.hpp>

#include <boost/program_options.hpp>
#include <display_timer.hpp>

#include "vtu_writer.hpp"

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
    ("pthread", boost::program_options::value<bool>()->default_value(true), "Use thread-based parallelism")
    ("solution,s", boost::program_options::value<int>()->default_value(0), "Select the solution")
    ("scaling_viscosity,v", boost::program_options::value<double>()->default_value(1.), "Set the scaling of the viscosity")
    ("scaling_permeabilityinv,p", boost::program_options::value<double>()->default_value(1.), "Set the scaling of the inverse of the permeability")
    ("export-matrix,e", "Export matrix to Matrix Market format")
    ("plot", boost::program_options::value<std::string>(), "Save plot of the solution to the given filename")
    ("solver", boost::program_options::value<std::string>()->default_value("PardisoLU"), "Choice of solver, not case dependent. Options are: PardisoLU, UMFPACK, PaStiXLU, PaStiXLLT, EigenLU, EigenBiCGSTAB (reverts to EigenLU if the selected solver is not available)")
    ("stabilization-parameter-stokes,x", boost::program_options::value<double>()->default_value(3.), "Set the stabilization parameter for Stokes terms")
    ("stabilization-parameter-darcy,y", boost::program_options::value<double>()->default_value(.3), "Set the stabilization parameter for Stokes terms")
    ("condition-number,c", boost::program_options::value<int>()->default_value(0), "Calculate condition number if >0 (requires the Spectra library); if >1, compute condition number and exits.")
    ;

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
  
  // Select the degree and stabilisation parameter
  size_t K = vm["degree"].as<size_t>();
  std::pair<double,double> stabilisation_parameter 
      = std::make_pair(vm["stabilization-parameter-stokes"].as<double>(),vm["stabilization-parameter-darcy"].as<double>());
  std::cout << "[main] Degree: " << K 
            << " \t(stabilisation parameters: " << stabilisation_parameter.first << ", "
            << stabilisation_parameter.second << ")" << std::endl;

  // Select the solution
  int solution = (vm.count("solution") ? vm["solution"].as<int>() : 0);
  Brinkman::MomentumForcingTermType f;
  Brinkman::CompressibilityForcingTermType g;
  BrinkmanParameters::ViscosityType mu(1.);    // Default viscosity constant equal to 1.
  BrinkmanParameters::PermeabilityInvType kappainv(1.);    // Default permeability constant equal to 1.
  Brinkman::VelocityType u;
  Brinkman::VelocityGradientType gradu;
  Brinkman::PressureType p;
  Brinkman::PressureGradientType gradp;

  // Select the scalings of mu and kappainv (change the values of the variables in the .hpp)
  scaling_mu = vm["scaling_viscosity"].as<double>();
  scaling_kappainv = vm["scaling_permeabilityinv"].as<double>();

  switch (solution) {
  case 0:
    std::cout << "[main] Linear solution ";
    f = linear_f;
    g = linear_g;
    u = linear_u;
    p = linear_p;
    mu = scaling_mu * linear_mu;
    kappainv = scaling_kappainv * linear_kappainv;
    gradu = linear_gradu;
    gradp = linear_gradp;
    break;

  case 1:
    std::cout << "[main] Quadratic solution ";
    f = quadratic_f;
    g = quadratic_g;
    u = quadratic_u;
    p = quadratic_p;
    mu = scaling_mu * quadratic_mu;
    kappainv = scaling_kappainv * quadratic_kappainv;
    gradu = quadratic_gradu;
    gradp = quadratic_gradp;
    break;

  case 2:
    std::cout << "[main] Cubic solution ";
    f = cubic_f;
    g = cubic_g;
    u = cubic_u;
    p = cubic_p;
    mu = scaling_mu * cubic_mu;
    kappainv = scaling_kappainv * cubic_kappainv;
    gradu = cubic_gradu;
    gradp = cubic_gradp;
    break;

  case 3:
    std::cout << "[main] Trigonometric solution ";
    f = trigonometric_f;
    g = trigonometric_g;
    u = trigonometric_u;
    p = trigonometric_p;
    mu = scaling_mu * trigonometric_mu;
    kappainv = scaling_kappainv * trigonometric_kappainv;
    gradu = trigonometric_gradu;
    gradp = trigonometric_gradp;
    break;

  case 4:
    std::cout << "[main] Varying regimes ";
    f = regimes_f;
    g = regimes_g;
    u = regimes_u;
    p = regimes_p;
    mu = scaling_mu * regimes_mu;
    kappainv = scaling_kappainv * regimes_kappainv;
    gradu = regimes_gradu;
    gradp = regimes_gradp;
    break;

  case 5:
    std::cout << "[main] BCs for cracked domain ";
    f = vcracked_f;
    g = vcracked_g;
    u = vcracked_u;
    p = vcracked_p;
    mu = scaling_mu * vcracked_mu;
    kappainv = scaling_kappainv * vcracked_kappainv;
    gradu = vcracked_gradu;
    gradp = vcracked_gradp;
    break;

  case 6:
    std::cout << "[main] BCs for lid-driven cavity surrounded by porous medium ";
    f = cavity_f;
    g = cavity_g;
    u = cavity_u;
    p = cavity_p;
    mu = scaling_mu * cavity_mu;
    kappainv = scaling_kappainv * cavity_kappainv;
    gradu = cavity_gradu;
    gradp = cavity_gradp;
    break;

  case 7:
    std::cout << "[main] Two regions split at x=1/2 ";
    f = tworegions_f;
    g = tworegions_g;
    u = tworegions_u;
    p = tworegions_p;
    mu = scaling_mu * tworegions_mu;
    kappainv = scaling_kappainv * tworegions_kappainv;
    gradu = tworegions_gradu;
    gradp = tworegions_gradp;
    break;

  default:
    std::cerr << "[main] ERROR: Unknown exact solution" << std::endl;
    exit(1);
  }
  std::cout << "(scaling viscosity: " << scaling_mu << "; scaling permeability: " << scaling_kappainv << ")" << std::endl;

  // Store physical parameters and select cells for stabilisation on the boundary of the domain
  BrinkmanParameters para(mu, kappainv);
  CellSelection BoundaryStab = [&para](const Cell & T)->bool { return (para.Cf(T)<1);};

  // Build the mesh and reorder faces to handle Dirichlet BCs (Dirichlet faces at the start -- pure Dirichlet for the moment)
  MeshBuilder meshbuilder = MeshBuilder(mesh_file);
  std::unique_ptr<Mesh> mesh_ptr = meshbuilder.build_the_mesh();
  BoundaryConditions BC("D", *mesh_ptr.get());
  BC.reorder_faces("start"); 

  boost::timer::cpu_timer timer;
  // Create velocity and pressure spaces
  timer.start();
  bool use_threads = (vm.count("pthread") ? vm["pthread"].as<bool>() : true);
  std::cout << "[main] " << (use_threads ? "Parallel execution" : "Sequential execution") << std:: endl;
  VHHOSpace vhho_space(*mesh_ptr, K, BoundaryStab, use_threads);
  GlobalDOFSpace p_space(*mesh_ptr, 0, 0, 0, PolynomialSpaceDimension<Cell>::Poly(K));
  timer.stop();
  double t_wall_vhhospace, t_proc_vhhospace;
  std::tie(t_wall_vhhospace, t_proc_vhhospace) = store_times(timer, "[main] Time creation spaces (wall/proc) ");

  // Assemble the problem
  timer.start();
  Brinkman br(vhho_space, p_space, BC, use_threads);
  br.stabilizationParameter() = std::make_pair(stabilisation_parameter.first,stabilisation_parameter.second);

  Eigen::VectorXd UDir = Eigen::VectorXd::Zero(br.numDirDOFs());
  br.assembleLinearSystem(f, g, para, u, UDir);
  timer.stop();
  double t_wall_model, t_proc_model;
  std::tie(t_wall_model, t_proc_model) = store_times(timer, "[main] Time model (wall/proc) ");

  // Export matrix if requested  
  if (vm.count("export-matrix")) {
    std::cout << "[main] Exporting matrix to Matrix Market format" << std::endl;
    saveMarket(br.systemMatrix(), "A_fullgradientbr.mtx");
    saveMarket(br.systemVector(), "b_fullgradientbr.mtx");
  }

  // Compute condition number if requested
  double cond_num = -1;
  int calculate_cond_num = vm["condition-number"].as<int>();
  if (calculate_cond_num > 0){
    #ifdef WITH_SPECTRA
    auto [max_eig, min_eig] = br.computeConditionNum();
    cond_num = std::sqrt(max_eig/double(min_eig));
    std::cout << FORMAT(25) << "[main] Condition number (min/max eig): " << cond_num << " (" << min_eig << "/" << max_eig << ")" << std::endl;
    #endif
    if (calculate_cond_num > 1){
      exit(0);
    }
  }
  
  // Select linear solver
  std::string name_solver = vm["solver"].as<std::string>();
  LinearSolver<Brinkman::SystemMatrixType> solver(name_solver);
  std::cout << "[main] Solving the system using " << solver.name() << std::endl;

  // Solve the problem
  timer.start();
  Eigen::VectorXd uph_solsystem = solver.compute_and_solve(br.systemMatrix(), br.systemVector());
    
  double absolute_residual_solver = (br.systemMatrix() * uph_solsystem - br.systemVector()).norm();
  double residual_solver = absolute_residual_solver / br.systemVector().norm();
  std::cout << "...residual = "<< residual_solver << " (absolute = " << absolute_residual_solver << ")" << std::endl;
  // Reconstruct solution
  //
  // uph_solsystem only contains the non-Dirichlet and non-SC DOFs; we augment it with the Dirichlet DOFs, to then recover the SC DOFS
  Eigen::VectorXd uph_nonscdofs = Eigen::VectorXd::Zero(br.numNonSCDOFs());
  uph_nonscdofs.head(br.numDirDOFs()) = UDir;
  uph_nonscdofs.tail(br.sizeSystem()) = uph_solsystem;
  Eigen::VectorXd uph_scdofs = br.scVector() + br.scMatrix() * uph_nonscdofs;

  // Vector with all DOFs
  Eigen::VectorXd uph = Eigen::VectorXd::Zero(br.dimVelocity() + br.dimPressure());

  // Velocity DOFs: boundary and non-SC, then SC at the end.
  uph.head(br.dimVelocity() - br.numSCDOFs_u()) = uph_nonscdofs.head(br.dimVelocity() - br.numSCDOFs_u());
  uph.segment(br.dimVelocity() - br.numSCDOFs_u(), br.numSCDOFs_u()) = uph_scdofs.head(br.numSCDOFs_u());
  // Pressure DOFs: we need to work element by element, first non-SC (number nloc_nonsc_p) and then SC
  size_t nloc_nonsc_p = br.pspace().numLocalDofsCell() - br.nloc_sc_p();
  size_t skel_offset_p = br.dimVelocity() - br.numSCDOFs_u(); // where p starts in uph_nonscdofs
  for (size_t iT = 0; iT < mesh_ptr->n_cells(); iT++){
    const Cell & T = *mesh_ptr->cell(iT);
    uph.segment(br.dimVelocity() + br.pspace().globalOffset(T), nloc_nonsc_p)
          = uph_nonscdofs.segment(skel_offset_p + iT * nloc_nonsc_p, nloc_nonsc_p);
    uph.segment(br.dimVelocity() + br.pspace().globalOffset(T) + nloc_nonsc_p, br.nloc_sc_p())
          = uph_scdofs.segment(br.numSCDOFs_u() + iT * br.nloc_sc_p(), br.nloc_sc_p());
  }
  
  timer.stop();
  double t_wall_solve, t_proc_solve;
  std::tie(t_wall_solve, t_proc_solve) = store_times(timer, "[main] Time solve (wall/proc) ");
    
  // Compute the error in the L2, H1 and energy norms
  Eigen::VectorXd upI = br.interpolate(u, p, MAX_DOE_CELL, MAX_DOE_FACE);
  Eigen::VectorXd eh = uph - upI;
  std::cout << "[main] Exact solution interpolated" << std::endl;

  Eigen::VectorXd uI = upI.head(br.vhhospace().dimension());
  Eigen::VectorXd euh = eh.head(br.vhhospace().dimension());
  std::vector<std::pair<double,double>> list_disc_norms = br.vhhospace().computeNorms(std::vector<Eigen::VectorXd> {euh, uI});
  double L2_err_u = list_disc_norms[0].first / list_disc_norms[1].first;
  std::cout << "[main] Discrete L2 error on u " << L2_err_u << std::endl;
  double H1_err_u = list_disc_norms[0].second / list_disc_norms[1].second;
  std::cout << "[main] Discrete H1 error on u " << H1_err_u << std::endl;

  std::vector<BrinkmanNorms> list_energy_norms = br.computeBrinkmanNorms(std::vector<Eigen::VectorXd> {eh, upI}, para);
  double energy_err = list_energy_norms[0].energy / list_energy_norms[1].energy;
  double L2_err_p = list_energy_norms[0].p / list_energy_norms[1].p;
  std::cout << "[main] Discrete L2 error p " << L2_err_p << std::endl;
  std::cout << "[main] Discrete energy error " << energy_err << std::endl;
  std::cout << "[main] Mesh diameter " << mesh_ptr->h_max() << std::endl;
  
  // Compute flux outside cavity when we're on a mesh of the cavity
  double flux_cavity = 0.;
  double area_surface = 0.;
  if ( mesh_file.find("CubeCavityWedge") != std::string::npos ){
    // Define surface
    std::vector<VectorRd> vert(4, VectorRd::Zero());
    vert[0] = VectorRd(1., 1., 0.);
    vert[1] = VectorRd(1., 1., -.75);
    vert[2] = VectorRd(1., 0., -.75);
    vert[3] = VectorRd(1., 0., 0.);
    VectorRd normal(1., 0., 0.);
    std::tie(flux_cavity, area_surface) = br.computeFlux(uph.head(br.vhhospace().dimension()), std::make_pair(vert, normal));

    std::cout << "[main] Flux outside cavity " << flux_cavity << " (area: " << area_surface << ")" << std::endl;
  }

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
    std::vector<VectorRd> u_vertex = vhho_space.computeVertexValues(uph.head(br.dimVelocity()));
    std::vector<double> p_vertex = br.pressureVertexValues(uph.tail(br.dimPressure()));
    
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
  out << "Scheme: hho-brinkman" << std::endl;
  out << "Solution: " << solution << std::endl;
  out << "ScalingViscosity: " << scaling_mu << std::endl;
  out << "ScalingPermeabilityInv: " << scaling_kappainv << std::endl;
  out << "Mesh: " << mesh_file << std::endl;
  out << "Degree: " << K << std::endl;
  out << "StabilisationParameterStokes: " << br.stabilizationParameter().first << std::endl;
  out << "StabilisationParameterDarcy: " << br.stabilizationParameter().second << std::endl;
  out << "ConditionNumber: " << cond_num << std::endl;
  out << "ResidualSolver: " << residual_solver << std::endl;
  out << "MeshSize: " << mesh_ptr->h_max() << std::endl;
  out << "NbCells: " << mesh_ptr->n_cells() << std::endl;
  out << "NbFaces: " << mesh_ptr->n_faces() << std::endl;
  out << "DimVelocitySpace: " << vhho_space.dimension() << std::endl;
  out << "DimPressureSpace: " << p_space.dimension() << std::endl;
  out << "SizeSystem: " << br.sizeSystem() << std::endl;

  out << "L2ErrorU: " << L2_err_u << std::endl;  
  out << "abs_L2ErrorU: " << list_disc_norms[0].first << std::endl;
  out << "L2NormU: " << list_disc_norms[1].first << std::endl;
  out << "H1ErrorU: " << H1_err_u << std::endl;  
  out << "abs_H1ErrorU: " << list_disc_norms[0].second << std::endl;
  out << "H1NormU: " << list_disc_norms[1].second << std::endl;
  out << "L2ErrorP: " << L2_err_p << std::endl;  
  out << "abs_L2ErrorP: " << list_energy_norms[0].p << std::endl;
  out << "L2NormP: " << list_energy_norms[1].p << std::endl;
  out << "EnergyError: " << energy_err << std::endl;  
  out << "abs_EnergyError: " << list_energy_norms[0].energy << std::endl;
  out << "EnergyNorm: " << list_energy_norms[1].energy << std::endl;

  out << "FluxCavity: " << flux_cavity << std::endl;

  out << "TwallVHHOSpace: " << t_wall_vhhospace << std::endl;  
  out << "TprocVHHOSpace: " << t_proc_vhhospace << std::endl;  
  out << "TwallModel: " << t_wall_model << std::endl;  
  out << "TprocModel: " << t_proc_model << std::endl;  
  out << "TwallSolve: " << t_wall_solve << std::endl;  
  out << "TprocSolve: " << t_proc_solve << std::endl;  
  out << std::flush;
  out.close();

  std::cout << "[main] Done" << std::endl;
  return 0;
}

//------------------------------------------------------------------------------
// Diffusion model
//------------------------------------------------------------------------------

Brinkman::Brinkman(
               const VHHOSpace & vhho_space,
               const GlobalDOFSpace & p_space,
               const BoundaryConditions & BC,
               bool use_threads,
               std::ostream & output
               )
  : m_vhhospace(vhho_space),
    m_pspace(p_space),
    m_BC(BC),
    m_use_threads(use_threads),
    m_output(output),
    m_nloc_sc_u(m_vhhospace.numLocalDofsCell()),
    m_nloc_sc_p(m_pspace.numLocalDofsCell()-1),
//    m_nloc_sc_u(0),
//    m_nloc_sc_p(0),
    m_A(sizeSystem(), sizeSystem()),
    m_b(Eigen::VectorXd::Zero(sizeSystem())),
    m_sc_A(numSCDOFs(), sizeSystem()),
    m_sc_b(Eigen::VectorXd::Zero(numSCDOFs())),
    m_stab_par(std::make_pair(1.,1.))
    
{
  m_output << "[Brinkman] Initializing" << std::endl;
  
  // To avoid performing static condensation, initialise m_nloc_sc_u = m_nloc_sc_p = 0
  // Note that velocity can be condensed without condensing pressure, but not the converse
}

//------------------------------------------------------------------------------

void Brinkman::assembleLinearSystem(
                                  const MomentumForcingTermType & f,
                                  const CompressibilityForcingTermType & g,
                                  const BrinkmanParameters & para,
                                  const VelocityType & u,
                                  Eigen::VectorXd & UDir
                                  )
{
  // Assemble all local contributions
  auto assemble_all = [this, f, g, para, u](
                                       size_t start,
                                       size_t end,
                                       std::list<Eigen::Triplet<double> > * triplets_sys,
                                       Eigen::VectorXd * rhs_sys,
                                       std::list<Eigen::Triplet<double> > * triplets_sc,
                                       Eigen::VectorXd * rhs_sc
                                       )->void
                      {
                        for (size_t iT = start; iT < end; iT++) {
                          this->_assemble_local_contribution(
                                                             iT,
                                                             this->_compute_local_contribution(iT, f, g, para),
                                                             *triplets_sys,
                                                             *rhs_sys,
                                                             *triplets_sc,
                                                             *rhs_sc
                                                             );
                        } // for iT
                      };
                      
  // Assemble the matrix and rhs
  if (m_use_threads) {
    m_output << "[Brinkman] Parallel assembly" << std::endl;
  }else{
    m_output << "[Brinkman] Sequential assembly" << std::endl;
  }
  Brinkman::SystemMatrixType global_matrix;  // Matrix without BCs (but after static condensation)
  Eigen::VectorXd global_rhs;     // RHS without BCs (but after static condensation)
  std::tie(global_matrix, global_rhs, m_sc_A, m_sc_b) 
        = parallel_assembly_system(
                                  m_vhhospace.mesh().n_cells(), 
                                  this->numNonSCDOFs(), 
                                  std::make_pair(this->numSCDOFs(), this->sizeSystem() + this->numDirDOFs()),
                                  this->numSCDOFs(), 
                                  assemble_all, 
                                  m_use_threads
                                  ); 
    
  // System matrix without Dirichlet boundary dofs (first ones)
  m_A = global_matrix.bottomRightCorner(sizeSystem(), sizeSystem());
  
  // Adjust RHS for boundary conditions
  // Creation of Dirichlet vector: L2 projection of solution only on Dirichlet faces (reordered at the start)
  for (size_t idF = 0; idF < m_BC.n_dir_faces(); idF++){
    Face* F = mesh().face(idF);
    QuadratureRule quadF = generate_quadrature_rule(*F, 2*m_vhhospace.degree()+2);
    
    auto phiF_quadF = evaluate_quad<Function>::compute(*m_vhhospace.faceBases(idF).Polykd, quadF);
    UDir.segment(idF * m_vhhospace.numLocalDofsFace(), m_vhhospace.numLocalDofsFace()) 
        = l2_projection<VHHOSpace::PolydBasisFaceType>(u, *m_vhhospace.faceBases(idF).Polykd, quadF, phiF_quadF);
  } // for idF
  
  m_b = global_rhs.tail(sizeSystem()) - global_matrix.bottomLeftCorner(sizeSystem(), numDirDOFs()) * UDir;
   
}

//------------------------------------------------------------------------------
  
std::pair<Eigen::MatrixXd, Eigen::VectorXd>
Brinkman::_compute_local_contribution(
                                    size_t iT, 
                                    const MomentumForcingTermType & f,
                                    const CompressibilityForcingTermType & g,
                                    const BrinkmanParameters & para
                                    )
{
  const Cell & T = *mesh().cell(iT);

  size_t dim_u = m_vhhospace.dimensionCell(iT);
  size_t dim_p = m_pspace.dimensionCell(iT);
  size_t dim_T = dim_u + dim_p + 1;   // +1 for Lagrange multiplier
  Eigen::MatrixXd AT = Eigen::MatrixXd::Zero(dim_T, dim_T);
  Eigen::VectorXd lT = Eigen::VectorXd::Zero(dim_T);
  
  //------------------------------------------------------------------------------
  // Local matrix
  //------------------------------------------------------------------------------

  // Parameters constant scalar in each element for the moment
  double muT = para.mu.value(T, T.center_mass());
  double kappainvT = para.kappainv.value(T, T.center_mass());
  double CfT = para.Cf(T);
  
  // Stokes: velocity-velocity
  MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_vhhospace.degree());  
  Eigen::MatrixXd mass_Pkdxd_T = GramMatrix(T, *m_vhhospace.cellBases(iT).Polykdxd, int_mono_2k);

  AT.topLeftCorner(dim_u, dim_u) 
      += muT * m_vhhospace.operators(iT).gradient.transpose() * mass_Pkdxd_T * m_vhhospace.operators(iT).gradient
         + muT * std::min(1., std::pow(CfT, -1)) * m_stab_par.first * m_vhhospace.operators(iT).stabilisation;
  
  // Darcy: velocity-velocity
  Eigen::MatrixXd PTtilde = (CfT>=1) * m_vhhospace.operators(iT).potential_div;
  PTtilde.rightCols(m_vhhospace.numLocalDofsCell()) 
        += (CfT<1) * Eigen::MatrixXd::Identity(m_vhhospace.numLocalDofsCell(), m_vhhospace.numLocalDofsCell());
  Eigen::MatrixXd mass_Pkd_T = GramMatrix(T, *m_vhhospace.cellBases(iT).Polykd, int_mono_2k);
  AT.topLeftCorner(dim_u, dim_u) += kappainvT * PTtilde.transpose() * mass_Pkd_T * PTtilde
                           + kappainvT * std::min(1., CfT) * m_stab_par.second * m_vhhospace.operators(iT).stabilisation_div;
  
  // Velocity-pressure
  // We compute -\int D_T v q by writing D_T v q = tr(G_T v) q= qId : G_T v and decomposing qId on a MatrixFamily based on the
  // basis of q
  auto basis_qId = IsotropicMatrixFamily<VHHOSpace::PolyBasisCellType, dimspace>(*m_vhhospace.cellBases(iT).Polyk);
  Eigen::MatrixXd BT = - GramMatrix(T, basis_qId, *m_vhhospace.cellBases(iT).Polykdxd, int_mono_2k) 
                            * m_vhhospace.operators(iT).gradient;

  AT.block(0, dim_u, dim_u, dim_p) = BT.transpose();
  AT.block(dim_u, 0, dim_p, dim_u) = -BT;

  // Lagrange multiplier: enforces \int p = 0
  FType<double> cst_fct_one = [](const VectorRd &x) -> double { return 1.0; };
  QuadratureRule quad_k_T = generate_quadrature_rule(T, m_vhhospace.degree() );
  Eigen::RowVectorXd intPk = 
      integrate(cst_fct_one, evaluate_quad<Function>::compute(*m_vhhospace.cellBases(iT).Polyk, quad_k_T), quad_k_T).transpose();

  AT.block(dim_T-1, dim_u, 1, dim_p) = intPk;
  AT.block(dim_u, dim_T-1, dim_p, 1) = -intPk.transpose();
  
  //------------------------------------------------------------------------------
  // Local vector: forcing in momentum, and divergence of velocity
  //------------------------------------------------------------------------------

  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * m_vhhospace.degree());
  lT.head(m_vhhospace.dimensionCell(iT)) = 
      PTtilde.transpose() * integrate(f, evaluate_quad<Function>::compute(*m_vhhospace.cellBases(iT).Polykd, quad_2k_T), quad_2k_T);
  
  lT.segment(m_vhhospace.dimensionCell(iT), m_pspace.dimensionCell(iT)) =
      integrate(g, evaluate_quad<Function>::compute(*m_vhhospace.cellBases(iT).Polyk, quad_2k_T), quad_2k_T);
  
  return std::make_pair(AT, lT);
}

//------------------------------------------------------------------------------

LocalStaticCondensation Brinkman::_compute_static_condensation(const size_t & iT) const
{
  const Cell & T = *mesh().cell(iT);
  
  // Dimensions
  size_t dim_u = m_vhhospace.dimensionCell(iT) - m_nloc_sc_u;     // number of velocity unknowns after SC
  size_t dim_p = m_pspace.dimensionCell(iT) - m_nloc_sc_p;     // number of pressure unknowns after SC
  size_t dim_dofs = dim_u + dim_p + 1;      // nb of dofs remaining after SC (including Lagrange multiplier)
  size_t dim_sc = m_nloc_sc_u + m_nloc_sc_p;      // nb of SC dofs in total

  // Creation of permutation matrix
  Eigen::MatrixXd Perm = Eigen::MatrixXd::Zero(dim_dofs+dim_sc, dim_dofs+dim_sc);
  Perm.topLeftCorner(dim_u, dim_u) = Eigen::MatrixXd::Identity(dim_u, dim_u);
  Perm.block(dim_u, dim_u + m_nloc_sc_u, dim_p, dim_p) = Eigen::MatrixXd::Identity(dim_p, dim_p);
  Perm(dim_u+dim_p, dim_u + m_nloc_sc_u + dim_p+m_nloc_sc_p) = 1.;   // Lagrange multiplier
  Perm.block(dim_u + dim_p + 1, dim_u, m_nloc_sc_u, m_nloc_sc_u) = Eigen::MatrixXd::Identity(m_nloc_sc_u, m_nloc_sc_u);
  Perm.block(dim_u + dim_p + 1 + m_nloc_sc_u, dim_u + m_nloc_sc_u + dim_p, m_nloc_sc_p, m_nloc_sc_p)
      = Eigen::MatrixXd::Identity(m_nloc_sc_p, m_nloc_sc_p);

  // Creation of global DOFs for system: IT_sys contains the non-SC dofs of u and of p (before handling BCs)
  std::vector<size_t> IT_sys(dim_dofs, 0);
  auto IT_u = m_vhhospace.globalDOFIndices(T);
  auto IT_p = m_pspace.globalDOFIndices(T);
  auto it_IT_sys = std::copy(IT_u.begin(), IT_u.begin()+dim_u, IT_sys.begin()); // put non-SC DOFs of u in IT_sys
  
  // The offset_u is the total number of non-SC DOFs for u, where the global non-SC DOFs of p start. However, getting the global indices of these non-SC DOFs of p is trickier than for u, because they correspond to a few (one, potentially, although it can be more flexible) in each element, and their numbering is therefore not consecutive
  size_t offset_u = dimVelocity() - numSCDOFs_u();
  std::transform(IT_p.begin(), IT_p.begin()+dim_p, it_IT_sys, 
                  [&offset_u, &iT, this](const size_t & index) { return index + offset_u - iT * this->m_nloc_sc_p; });
  // Lagrange multiplier
  IT_sys[dim_dofs-1] = numNonSCDOFs()-1; 
    
  // Creation of global DOFs for SC operator: IT_sc contains global cell dofs of u (offset to start at 0) and p (offset to start after those of u, and with the same caveat as above regarding the non-consecutive numbering, which explains the offset_p)
  std::vector<size_t> IT_sc(dim_sc, 0);  
  if (dim_sc>0){
    auto it_IT_sc = std::transform(IT_u.begin()+dim_u, IT_u.end(), IT_sc.begin(), [&offset_u](const size_t & index) { return index - offset_u; });
    size_t offset_p = numSCDOFs_u() - (iT + 1) * (m_pspace.dimensionCell(iT) - m_nloc_sc_p);
    std::transform(IT_p.begin()+dim_p, IT_p.end(), it_IT_sc, 
                [&offset_p](const size_t & index) { return offset_p + index; });
  }  
  
  return LocalStaticCondensation(Perm, IT_sys, IT_sc);
}

//------------------------------------------------------------------------------

void Brinkman::_assemble_local_contribution(
                                        size_t iT,
                                        const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT,
                                        std::list<Eigen::Triplet<double> > & triplets_sys,
                                        Eigen::VectorXd & rhs_sys,
                                        std::list<Eigen::Triplet<double> > & triplets_sc,
                                        Eigen::VectorXd & rhs_sc
                                        )
{

  // Get information for local static condensation
  LocalStaticCondensation locSC = _compute_static_condensation(iT);

  Eigen::MatrixXd AT_sys, AT_sc;
  Eigen::VectorXd bT_sys, bT_sc;
  std::tie(AT_sys, bT_sys, AT_sc, bT_sc) = locSC.compute(lsT);
  
  // STATICALLY CONDENSED SYSTEM
  std::vector<size_t> IT_sys = locSC.globalDOFs_gl();
  for (size_t i = 0; i < locSC.dim_gl(); i++){
    rhs_sys(IT_sys[i]) += bT_sys(i);
    for (size_t j = 0; j < locSC.dim_gl(); j++){    
      triplets_sys.emplace_back(IT_sys[i], IT_sys[j], AT_sys(i,j));
    }
  }
  // RECOVERY OPERATOR
  std::vector<size_t> IT_sc = locSC.globalDOFs_sc();
  for (size_t i = 0; i < locSC.dim_sc(); i++){
    rhs_sc(IT_sc[i]) += bT_sc(i);
    for (size_t j = 0; j < locSC.dim_gl(); j++){
      triplets_sc.emplace_back(IT_sc[i], IT_sys[j], AT_sc(i,j));
    }
  }

}

//------------------------------------------------------------------------------
Eigen::VectorXd Brinkman::interpolate(
                    const VelocityType & u,
                    const PressureType & p,
                    const int doe_cell,
                    const int doe_face
                    ) const
{
  Eigen::VectorXd upI = Eigen::VectorXd::Zero(dimVelocity() + dimPressure());
  
  // Degrees of quadrature rules
  size_t dqr_cell = (doe_cell >= 0 ? doe_cell : 2 * m_vhhospace.degree() + 3);
  size_t dqr_face = (doe_face >= 0 ? doe_face : 2 * m_vhhospace.degree() + 3);

  // Interpolate velocity
  upI.head(dimVelocity()) = m_vhhospace.interpolate(u, dqr_cell, dqr_face);
  
  // Interpolate pressure
  size_t offset_p = dimVelocity();
  std::function<void(size_t, size_t)> interpolate_pressure
    = [this, &p, &upI, &dqr_cell, &offset_p](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          const Cell & T = *mesh().cell(iT);
          QuadratureRule quad_dqr_T = generate_quadrature_rule(T, dqr_cell);
          auto basis_Pk_T_quad = evaluate_quad<Function>::compute(*m_vhhospace.cellBases(iT).Polyk, quad_dqr_T);
          upI.segment(offset_p + m_pspace.globalOffset(T), m_pspace.numLocalDofsCell()) 
            = l2_projection(p, *m_vhhospace.cellBases(iT).Polyk, quad_dqr_T, basis_Pk_T_quad, GramMatrix(T, *m_vhhospace.cellBases(iT).Polyk));
        } // for iT
      };
  parallel_for(mesh().n_cells(), interpolate_pressure, m_use_threads);

  return upI;
}

//------------------------------------------------------------------------------

std::vector<BrinkmanNorms> Brinkman::computeBrinkmanNorms( 
                            const std::vector<Eigen::VectorXd> & list_dofs,
                            const BrinkmanParameters & para
                            ) const
{
  size_t nb_vectors = list_dofs.size();
  std::vector<Eigen::VectorXd> local_sqnorms_u(nb_vectors, Eigen::VectorXd::Zero(mesh().n_cells()));
  std::vector<Eigen::VectorXd> local_sqnorms_p(nb_vectors, Eigen::VectorXd::Zero(mesh().n_cells()));

  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &list_dofs, &para, &local_sqnorms_u, &local_sqnorms_p, &nb_vectors](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *m_vhhospace.mesh().cell(iT);
        double muT = para.mu.value(T, T.center_mass());
        double kappainvT = para.kappainv.value(T, T.center_mass());
        double CfT = para.Cf(T);
        
        // Mass matrix for P^k(T), P^k(T)^d and (P^k(T))^{dxd}
	      MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_vhhospace.degree());
        Eigen::MatrixXd mass_Pk_T = GramMatrix(T, *m_vhhospace.cellBases(iT).Polyk, int_mono_2k);
        Eigen::MatrixXd mass_Pkd_T_scaled = kappainvT * GramMatrix(T, *m_vhhospace.cellBases(iT).Polykd, int_mono_2k);
        Eigen::MatrixXd mass_Pkdxd_T_scaled = muT * GramMatrix(T, *m_vhhospace.cellBases(iT).Polykdxd, int_mono_2k);
        Eigen::MatrixXd GT = m_vhhospace.operators(iT).gradient;
        Eigen::MatrixXd ST_scaled = muT * std::min(1., std::pow(CfT, -1)) * m_stab_par.first * m_vhhospace.operators(iT).stabilisation;
        Eigen::MatrixXd PTdiv = m_vhhospace.operators(iT).potential_div;
        Eigen::MatrixXd STdiv_scaled = kappainvT * std::min(1., CfT) * m_stab_par.second *m_vhhospace.operators(iT).stabilisation_div;

        Eigen::MatrixXd PTtilde = (CfT>=1) * PTdiv;
        PTtilde.rightCols(m_vhhospace.numLocalDofsCell()) 
            += (CfT<1) * Eigen::MatrixXd::Identity(m_vhhospace.numLocalDofsCell(), m_vhhospace.numLocalDofsCell());

        for (size_t i=0; i<nb_vectors; i++){
          Eigen::VectorXd uT = m_vhhospace.restrict(T, list_dofs[i].head(dimVelocity()));
          Eigen::VectorXd pT = m_pspace.restrict(T, list_dofs[i].tail(dimPressure()));

          // Velocity
          local_sqnorms_u[i](iT) 
              += uT.transpose() * 
                  (GT.transpose() * mass_Pkdxd_T_scaled * GT + ST_scaled
                  + PTtilde.transpose() * mass_Pkd_T_scaled * PTtilde + STdiv_scaled
                  ) * uT;

          // Pressure
          local_sqnorms_p[i](iT) += pT.transpose() * mass_Pk_T * pT;
        }
      }
    };
  parallel_for(m_vhhospace.mesh().n_cells(), compute_local_squarednorms, m_use_threads);
  
  // Vector of outputs
  std::vector<BrinkmanNorms> list_norms(nb_vectors, BrinkmanNorms(0., 0.));
  for (size_t i=0; i<nb_vectors; i++){
    list_norms[i] = BrinkmanNorms( 
                          std::sqrt( std::abs(local_sqnorms_u[i].sum()) ),
                          std::sqrt( std::abs(local_sqnorms_p[i].sum()) ) 
                          );
  }
  
  return list_norms;
}

//-----------------------------------------------------------------------------------------

std::vector<double> Brinkman::pressureVertexValues(const Eigen::VectorXd & p) const
{
  std::vector<double> values(mesh().n_vertices(), 0.);

  for (Vertex * V : mesh().get_vertices()){
    size_t iV = V->global_index();
    
    // Vertex value calculated as average of all values comming from the cells around V
    for (Cell * T : V->get_cells()){
      size_t iT = T->global_index();
      Eigen::VectorXd pT = m_pspace.restrict(*T, p);
      
      for (size_t i=0; i < m_vhhospace.cellBases(iT).Polyk->dimension(); i++){
        values[iV] += pT(i) * m_vhhospace.cellBases(iT).Polyk->function(i, V->coords());
      }
    }
    
    values[iV] /= V->n_cells();
  }

  return values;
}


//------------------------------------------------------------------------------

std::pair<double, double> Brinkman::computeConditionNum() const
{
  #ifdef WITH_SPECTRA
  SystemMatrixType M = systemMatrix().transpose()*systemMatrix();
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

//-----------------------------------------------------------------------------------------

std::pair<double,double> Brinkman::computeFlux(const Eigen::VectorXd & u, const std::pair<std::vector<VectorRd>, VectorRd> & surf) const
{
  double flux = 0.;
  double area = 0.;
  
  // Vertices, normal (also as a constant function), center of the surface and point above the surface
  std::vector<VectorRd> vertices_surf = surf.first;
  VectorRd normal_surf = surf.second;
  FType<VectorRd> cst_fct_normal_surf = [&normal_surf](const VectorRd &z) -> VectorRd { return normal_surf; };
  VectorRd center_surf = VectorRd::Zero();
  for (size_t i=0; i < vertices_surf.size(); i++){
    center_surf += vertices_surf[i];
  }
  center_surf /= vertices_surf.size();
  VectorRd apex_surf = center_surf + normal_surf;
  
  // Function to determine if a point x is in the surface (check if it's in one of the triangles formed by 2 consecutive vertices and the center of the surface)
  std::function<bool(const VectorRd&)> is_in_surf = [&vertices_surf, &center_surf, &apex_surf](const VectorRd & x)->bool
      {
        double eps = 1e-10; // Tolerance, should probably be a factor of the surface diameter actually...

        // Loop over the simplices
        for (size_t i=0; i < vertices_surf.size() ; i++){
          // Two consecutive vertices
          VectorRd v1 = vertices_surf[i];
          VectorRd v2 = vertices_surf[(i+1) % vertices_surf.size()];
          // Find barycentric coordinates of x with respect to simplex formed by the apex, the center, and the two vertices
          MatrixRd A = MatrixRd::Zero();
          A.col(0) = center_surf - apex_surf;
          A.col(1) = v1 - apex_surf;
          A.col(2) = v2 - apex_surf;
          VectorRd lambda = A.partialPivLu().solve(x-apex_surf);
          
          // Check if x is in the triangle center-v1-v2
          if (lambda(0)>-eps && lambda(1)>-eps && lambda(2)>-eps && std::abs(lambda.sum() - 1.)<eps) return true;
        }
        
        return false;
      };
  
  // Loop over faces and integrate on those that have all vertices in the convex hull of the four vertices of the surface
  for (size_t iF = 0; iF < mesh().n_faces(); iF++){
    const Face & F = *mesh().face(iF);
    bool takeF = true;
    for (const Vertex * V : F.get_vertices()){
      takeF = takeF && is_in_surf(V->coords());
    }

    if (takeF){
      // Update area
      area += F.measure();
      
      // Integrate uF.n
      Eigen::VectorXd uF = m_vhhospace.restrict(F, u);
      QuadratureRule quad_k_F = generate_quadrature_rule(F, m_vhhospace.degree() );
      Eigen::RowVectorXd int_nF_PkF2 = 
        integrate(cst_fct_normal_surf, evaluate_quad<Function>::compute(*m_vhhospace.faceBases(iF).Polykd, quad_k_F), quad_k_F).transpose();
      
      flux += int_nF_PkF2 * uF;
      
    }
  }
  
  return std::make_pair(flux, area);
}



