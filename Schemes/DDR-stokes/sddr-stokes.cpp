// Author: Jerome Droniou (jerome.droniou@monash.edu)

#include <fstream>
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include <thread>

#include "sddr-stokes.hpp"
#include <parallel_for.hpp>
#include <max_degrees_quadratures.hpp>
#include <GMpoly_cell.hpp>
#include <linearsolver.hpp>

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
    ("pressure_scaling", boost::program_options::value<double>()->default_value(1), "Select the pressure scaling")
    ("export-matrix,e", "Export matrix to Matrix Market format")
    ("solver", boost::program_options::value<std::string>()->default_value("PardisoLU"), "Choice of solver, not case dependent. Options are: PardisoLU, UMFPACK, PaStiXLU, PaStiXLLT, EigenLU, EigenBiCGSTAB (reverts to EigenLU if the selected solver is not available)")
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
  std::cout << FORMAT(25) << "[main] Degree" << K << std::endl;

  // Select the solution
  int solution = (vm.count("solution") ? vm["solution"].as<int>() : 0);
  Stokes::ForcingTermType f;
  Stokes::VelocityType u;
  Stokes::VorticityType omega;
  Stokes::PressureType p;
  Stokes::PressureGradientType grad_p;
  Stokes::ViscosityType nu(1.);
  
  pressure_scaling = vm["pressure_scaling"].as<double>();

  switch (solution) {
  case 0:
    std::cout << "[main] Trigonometric solution" << std::endl;
    f = trigonometric_f;
    u = trigonometric_u;
    omega = trigonometric_curl_u;
    p = trigonometric_p;
    grad_p = trigonometric_grad_p;
    nu = trigonometric_nu;
    break;

  case 1:
    std::cout << "[main] Linear solution" << std::endl;
    f = linear_f;
    u = linear_u;
    omega = linear_curl_u;
    p = linear_p;
    grad_p = linear_grad_p;
    nu = linear_nu;
    break;

  case 2:
    std::cout << "[main] Linear velocity, trigonometric pressure" << std::endl;
    f = trigonometric_grad_p;
    u = linear_u;
    omega = linear_curl_u;
    p = trigonometric_p;
    grad_p = trigonometric_grad_p;
    nu = linear_nu;
    break;

  case 3:
    std::cout << "[main] Field solution" << std::endl;
    f = field_f;
    u = field_u;
    omega = field_curl_u;
    p = field_p;
    grad_p = field_grad_p;
    nu = field_nu;
    break;

  case 4:
    std::cout << "[main] Vertical flow" << std::endl;
    f = vertical_f;
    u = vertical_u;
    omega = vertical_curl_u;
    p = vertical_p;
    grad_p = vertical_grad_p;
    break;

  default:
    std::cerr << "[main] ERROR: Unknown exact solution" << std::endl;
    exit(1);    
  }


  // Build the mesh
  MeshBuilder meshbuilder = MeshBuilder(mesh_file);
  std::unique_ptr<Mesh> mesh_ptr = meshbuilder.build_the_mesh();

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
  Stokes st(ddr_core, use_threads);  
  if(vm.count("stabilization-parameter")) {
    st.stabilizationParameter() = vm["stabilization-parameter"].as<double>();
  }
  st.assembleLinearSystem(f, u, omega, nu);
  timer.stop();
  double t_wall_model, t_proc_model;
  std::tie(t_wall_model, t_proc_model) = store_times(timer, "[main] Time model (wall/proc) ");

  // Export matrix if requested  
  if (vm.count("export-matrix")) {
    std::cout << "[main] Exporting matrix to Matrix Market format" << std::endl;
    saveMarket(st.systemMatrix(), "A_stokes.mtx");
    saveMarket(st.systemVector(), "b_stokes.mtx");
  }

  //------ Solve the problem ------//
  timer.start();
  std::string name_solver = vm["solver"].as<std::string>();
  LinearSolver<Stokes::SystemMatrixType> solver(name_solver);
  std::cout << "[main] Solving the system using " << solver.name() << std::endl;
  Eigen::VectorXd uplambdah_condensed = solver.compute_and_solve(st.systemMatrix(), st.systemVector());
  std::cout << "[main] residual linear solver: " << solver.residual(st.systemMatrix(), st.systemVector(), uplambdah_condensed) << std::endl;

  
  // Re-create statically condensed unknowns
  Eigen::VectorXd uph = Eigen::VectorXd::Zero(st.dimensionSpace());
  Eigen::VectorXd sc_unknowns = st.scVector() + st.scMatrix() * uplambdah_condensed;
  uph.head(st.sxCurl().dimension() - st.nbSCDOFs_u()) = uplambdah_condensed.head(st.sxCurl().dimension() - st.nbSCDOFs_u()); 
  uph.segment(st.sxCurl().dimension() - st.nbSCDOFs_u(), st.nbSCDOFs_u()) = sc_unknowns.head(st.nbSCDOFs_u());
  uph.segment(st.sxCurl().dimension(), st.sxGrad().dimension() - st.nbSCDOFs_p()) 
      = uplambdah_condensed.segment(st.sxCurl().dimension() - st.nbSCDOFs_u(), st.sxGrad().dimension() - st.nbSCDOFs_p());
  uph.tail(st.nbSCDOFs_p()) = sc_unknowns.tail(st.nbSCDOFs_p());

  timer.stop();
  double t_wall_solve, t_proc_solve;
  std::tie(t_wall_solve, t_proc_solve) = store_times(timer, "[main] Time solve (wall/proc) ");
  
  // Interpolate of exact solution and error vector
  Eigen::VectorXd upI = Eigen::VectorXd::Zero(st.dimensionSpace());  
  upI.head(st.sxCurl().dimension()) = st.sxCurl().interpolate(u, MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
  upI.tail(st.sxGrad().dimension()) = st.sxGrad().interpolate(p, MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
  Eigen::VectorXd eph = uph - upI;

  std::cout << "[main] Compute errors" << std::endl;
  timer.start();
  // Errors in discrete norms
  std::vector<StokesNorms> list_norms = st.computeStokesNorms(std::vector<Eigen::VectorXd> {eph, upI});
  StokesNorms errors = list_norms[0];
  StokesNorms norms = list_norms[1];  
  std::cout << "[main] Hcurl norm u= " << norms.hcurl_u << "; Hgrad norm p= " << norms.hgrad_p <<std::endl;
  double error_hcurl_u = (norms.hcurl_u < 1e-12 ? errors.hcurl_u : errors.hcurl_u / norms.hcurl_u);
  double error_hgrad_p = (norms.hgrad_p < 1e-12 ? errors.hgrad_p : errors.hgrad_p / norms.hgrad_p);
  std::cout << "[main] Hcurl error u= " << error_hcurl_u << "; Hgrad error p= " << error_hgrad_p << std::endl;
  // Errors in continuous norms
  std::pair<StokesNorms, StokesNorms> c_errors_norms = st.computeContinuousErrorsNorms(uph, u, omega, p, grad_p);
  StokesNorms c_errors = c_errors_norms.first;
  StokesNorms c_norms = c_errors_norms.second;
  double c_error_hcurl_u = (c_norms.hcurl_u < 1e-12 ? c_errors.hcurl_u : c_errors.hcurl_u / c_norms.hcurl_u);
  double c_error_grad_p = (c_norms.grad_p < 1e-12 ? c_errors.grad_p : c_errors.grad_p / c_norms.grad_p);
  std::cout << "[main] Continuous Hcurl error= " << c_error_hcurl_u << "; Continuous Hgrad error= " << c_error_grad_p <<std::endl;
  timer.stop();
  double t_wall_norms, t_proc_norms;
  std::tie(t_wall_norms, t_proc_norms) = store_times(timer, "[main] Time norms (wall/proc) ");
  
  std::cout << "[main] Mesh diameter " << mesh_ptr->h_max() << std::endl;
  // Write results to file
  std::ofstream out("results.txt");
  out << "Solution: " << solution << std::endl;
  out << "PressureScaling: " << pressure_scaling << std::endl;
  out << "StabilisationParameter: " << st.stabilizationParameter() << std::endl;
  out << "Mesh: " << mesh_file << std::endl;
  out << "Degree: " << K << std::endl;
  out << "MeshSize: " << mesh_ptr->h_max() << std::endl;
  out << "NbCells: " << mesh_ptr->n_cells() << std::endl;
  out << "NbFaces: " << mesh_ptr->n_faces() << std::endl;
  out << "NbEdges: " << mesh_ptr->n_edges() << std::endl;
  out << "DimSXCurl: " << st.sxCurl().dimension() << std::endl;
  out << "DimSXGrad: " << st.sxGrad().dimension() << std::endl;
  out << "SizeSCsystem: " << st.systemMatrix().rows() << std::endl;
  // Discrete errors and norms
  out << "E_HcurlVel: " << error_hcurl_u << std::endl;
  out << "E_HgradPre: " << error_hgrad_p << std::endl;
  out << "E_L2Vel: " << (norms.u < 1e-12 ? errors.u : errors.u/norms.u) << std::endl;
  out << "E_L2CurlVel: " << (norms.curl_u < 1e-12 ? errors.curl_u : errors.curl_u/norms.curl_u) << std::endl;
  out << "E_L2Pre: " << (norms.p < 1e-12 ? errors.p : errors.p/norms.p) << std::endl;
  out << "E_L2GradPre: " << (norms.grad_p < 1e-12 ? errors.grad_p : errors.grad_p/norms.grad_p) << std::endl;
  out << "absE_HcurlVel: " << errors.hcurl_u << std::endl;
  out << "absE_HgradPre: " << errors.hgrad_p << std::endl;
  out << "absE_L2Vel: " << errors.u << std::endl;
  out << "absE_L2CurlVel: " << errors.curl_u << std::endl;
  out << "absE_L2Pre: " << errors.p << std::endl;
  out << "absE_L2GradPre: " << errors.grad_p << std::endl;
  out << "N_L2Vel: " << norms.u << std::endl;
  out << "N_L2CurlVel: " << norms.curl_u << std::endl;
  out << "N_L2Pre: " << norms.p << std::endl;
  out << "N_L2GradPre: " << norms.grad_p << std::endl;
  // Continuous errors and norms
  out << "E_cHcurlVel: " << c_error_hcurl_u << std::endl;
  out << "E_cHgradPre: " << c_error_grad_p << std::endl;
  out << "E_cL2Vel: " << (c_norms.u < 1e-12 ? c_errors.u : c_errors.u/c_norms.u) << std::endl;
  out << "E_cL2CurlVel: " << (c_norms.curl_u < 1e-12 ? c_errors.curl_u : c_errors.curl_u/c_norms.curl_u) << std::endl;
  out << "E_cL2Pre: " << (c_norms.p < 1e-12 ? c_errors.p : c_errors.p/c_norms.p) << std::endl;
  out << "E_cL2GradPre: " << (c_norms.grad_p < 1e-12 ? c_errors.grad_p : c_errors.grad_p/c_norms.grad_p) << std::endl;
  out << "absE_cHcurlVel: " << c_errors.hcurl_u << std::endl;
  out << "absE_cHgradPre: " << c_errors.hgrad_p << std::endl;
  out << "absE_cL2Vel: " << c_errors.u << std::endl;
  out << "absE_cL2CurlVel: " << c_errors.curl_u << std::endl;
  out << "absE_cL2Pre: " << c_errors.p << std::endl;
  out << "absE_cL2GradPre: " << c_errors.grad_p << std::endl;
  out << "N_cL2Vel: " << c_norms.u << std::endl;
  out << "N_cL2CurlVel: " << c_norms.curl_u << std::endl;
  out << "N_cL2Pre: " << c_norms.p << std::endl;
  out << "N_cL2GradPre: " << c_norms.grad_p << std::endl;
  out << "TwallDDRCore: " << t_wall_ddrcore << std::endl;  
  out << "TprocDDRCore: " << t_proc_ddrcore << std::endl;  
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
// Stokes
//------------------------------------------------------------------------------

Stokes::Stokes(
               const DDRCore & ddrcore,
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
    m_A(sizeSystem(), sizeSystem()),
    m_b(Eigen::VectorXd::Zero(sizeSystem())),
    m_sc_A(nbSCDOFs(), sizeSystem()),
    m_sc_b(Eigen::VectorXd::Zero(nbSCDOFs())),
    m_stab_par(0.1)    
{
  m_output << "[Stokes] Initializing" << std::endl;

  // To avoid performing static condensation, initialise m_nloc_sc_u and m_nloc_sc_p at 0
}

//------------------------------------------------------------------------------

void Stokes::assembleLinearSystem(
                                  const ForcingTermType & f,
                                  const VelocityType & u,
                                  const VorticityType & omega,
                                  const ViscosityType & nu
                                  )
{
  // Interpolate of forcing term in SXCurl
  Eigen::VectorXd interp_f = m_sxcurl.interpolate(f, MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);

  // Assemble all local contributions
  auto assemble_all = [this, interp_f, u, nu](
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
                                                             this->_compute_local_contribution(iT, interp_f, nu),
                                                             *triplets_sys,
                                                             *rhs_sys,
                                                             *triplets_sc,
                                                             *rhs_sc
                                                             );
                        } // for iT
                      };
                      
  // Assemble the matrix and rhs
  if (m_use_threads) {
    m_output << "[Stokes] Parallel assembly" << std::endl;
  }else{
    m_output << "[Stokes] Sequential assembly" << std::endl;
  }
  std::tie(m_A, m_b, m_sc_A, m_sc_b) = parallel_assembly_system(m_ddrcore.mesh().n_cells(), this->sizeSystem(), std::make_pair(this->nbSCDOFs(), this->sizeSystem()), this->nbSCDOFs(), assemble_all, m_use_threads);

  // Assemble boundary conditions
  for (auto iF : m_ddrcore.mesh().get_b_faces()) {
    const Face & F = *iF;
    
    // Unit normal vector to F pointing out of the domain
    const Cell & TF = *F.cell(0);
    Eigen::Vector3d nF = TF.face_normal(TF.index_face(&F));

    // Degree of quadratures to compute boundary conditions
    const size_t dqrbc = MAX_DOE_FACE;

    // Boundary condition on the tangential component of the vorticity
    {
      FType<Eigen::Vector3d> omega_cross_nF = [&omega, &nF](const Eigen::Vector3d & x) {
                                                        return omega(x).cross(nF);
                                                        };

      QuadratureRule quad_dqrbc_F = generate_quadrature_rule(F, dqrbc);
      Eigen::VectorXd bF = integrate(omega_cross_nF, evaluate_quad<Function>::compute(*m_ddrcore.faceBases(F.global_index()).Polyk2, quad_dqrbc_F), quad_dqrbc_F).transpose() * m_sxcurl.facePotential(F.global_index());
      auto I_F = m_sxcurl.globalDOFIndices(F);
      for (size_t i = 0; i < I_F.size(); i++) {
        m_b(I_F[i]) += bF(i);
      } // for i
    }

    // Boundary condition on the normal component of the velocity    
    {
      FType<double> u_dot_nF = [&u, &nF](const Eigen::Vector3d & x) {
                                               return u(x).dot(nF);
                                               };

      QuadratureRule quad_dqrbc_F = generate_quadrature_rule(F, dqrbc);
      Eigen::VectorXd bF = - integrate(u_dot_nF, evaluate_quad<Function>::compute(*m_ddrcore.faceBases(F.global_index()).Polykpo, quad_dqrbc_F), quad_dqrbc_F).transpose() * m_sxgrad.facePotential(F.global_index());
      auto I_F = m_sxgrad.globalDOFIndices(F);
      size_t dim_sxcurl = m_sxcurl.dimension() - nbSCDOFs_u();  // Gradients DOFs start here after SC
      for (size_t i = 0; i < I_F.size(); i++) {
        m_b(dim_sxcurl + I_F[i]) += bF(i);
      } // for i
    }
  } // for iF
}

//------------------------------------------------------------------------------
 
std::pair<Eigen::MatrixXd, Eigen::VectorXd>
    Stokes::_compute_local_contribution(
                                        size_t iT, 
                                        const Eigen::VectorXd & interp_f,
                                        const ViscosityType & nu
                                        )
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  size_t dim_sxcurl_T = m_sxcurl.dimensionCell(iT);
  size_t dim_sxgrad_T = m_sxgrad.dimensionCell(iT);
  size_t dim_T = dim_sxcurl_T + dim_sxgrad_T;   
  
  Eigen::MatrixXd AT = Eigen::MatrixXd::Zero(dim_T+1, dim_T+1);  // +1 for Lagrange multiplier
  Eigen::VectorXd lT = Eigen::VectorXd::Zero(dim_T+1);

  //------------------------------------------------------------------------------
  // Local matrix
  //------------------------------------------------------------------------------

  // Mass matrix for (P^k(T))^3  
  MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_ddrcore.degree());
  Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_2k);

  // aT
  AT.topLeftCorner(dim_sxcurl_T, dim_sxcurl_T) = m_sxdiv.computeL2ProductCurl(iT, m_sxcurl, "both", m_stab_par, mass_Pk3_T, nu);
  
  // bT
  Eigen::MatrixXd BT = m_sxcurl.computeL2ProductGradient(iT, m_sxgrad, "right", m_stab_par, mass_Pk3_T); 
  AT.block(0, dim_sxcurl_T, dim_sxcurl_T, dim_sxgrad_T ) = BT;
  AT.block(dim_sxcurl_T, 0, dim_sxgrad_T, dim_sxcurl_T ) = -BT.transpose();

  // Lagrange multiplier: enforce \int P^{k+1}p = 0
  FType<double> cst_fct_one = [](const Eigen::Vector3d &x) -> double { return 1.0; };
  QuadratureRule quad_kpo_T = generate_quadrature_rule(T, m_ddrcore.degree()+1 );
  // intPko is line vector representing \int_T P^{k+1}
  Eigen::RowVectorXd intPkpo = 
      integrate(cst_fct_one, evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polykpo, quad_kpo_T), quad_kpo_T).transpose()
        * m_sxgrad.cellPotential(iT);

  AT.block(dim_T, dim_sxcurl_T, 1, dim_sxgrad_T) = intPkpo;
  AT.block(dim_sxcurl_T, dim_T, dim_sxgrad_T, 1) = -intPkpo.transpose();
    
  //------------------------------------------------------------------------------
  // Local source vector
  //------------------------------------------------------------------------------
  lT.head(dim_sxcurl_T) = m_sxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T) * m_sxcurl.restrictCell(iT, interp_f);

  return std::make_pair(AT, lT);
}

//------------------------------------------------------------------------------

LocalStaticCondensation Stokes::_compute_static_condensation(const size_t & iT) const
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  // Dimensions
  size_t dim_u = m_sxcurl.dimensionCell(iT) - m_nloc_sc_u(iT);     // number of velocity unknowns after SC
  size_t dim_p = m_sxgrad.dimensionCell(iT) - m_nloc_sc_p(iT);     // number of pressure unknowns after SC
  size_t dim_dofs = dim_u + dim_p + 1;      // nb of dofs remaining after SC (including Lagrange multiplier)
  size_t dim_sc = m_nloc_sc_u(iT) + m_nloc_sc_p(iT);      // nb of SC dofs in total

  // Creation of permutation matrix
  Eigen::MatrixXd Perm = Eigen::MatrixXd::Zero(dim_dofs+dim_sc, dim_dofs+dim_sc);
  Perm.topLeftCorner(dim_u, dim_u) = Eigen::MatrixXd::Identity(dim_u, dim_u);
  Perm.block(dim_u, dim_u + m_nloc_sc_u(iT), dim_p, dim_p) = Eigen::MatrixXd::Identity(dim_p, dim_p);
  Perm(dim_u+dim_p, size_t(dim_u + m_nloc_sc_u(iT) + dim_p+m_nloc_sc_p(iT))) = 1.;   // Lagrange multiplier
  Perm.block(dim_u + dim_p + 1, dim_u, m_nloc_sc_u(iT), m_nloc_sc_u(iT)) = Eigen::MatrixXd::Identity(m_nloc_sc_u(iT), m_nloc_sc_u(iT));
  Perm.block(dim_u + dim_p + 1 + m_nloc_sc_u(iT), dim_u + m_nloc_sc_u(iT) + dim_p, m_nloc_sc_p(iT), m_nloc_sc_p(iT))
      = Eigen::MatrixXd::Identity(m_nloc_sc_p(iT), m_nloc_sc_p(iT));

  // Creation of global DOFs for system: IT_sys contains the skeletal dofs of u and of p
  std::vector<size_t> IT_sys(dim_dofs, 0);
  auto IT_sxcurl = m_sxcurl.globalDOFIndices(T);
  auto IT_sxgrad = m_sxgrad.globalDOFIndices(T);
  auto it_IT_sys = std::copy(IT_sxcurl.begin(), IT_sxcurl.begin()+dim_u, IT_sys.begin()); // put skeletal DOFs of u in IT_sys
  size_t offset = m_sxcurl.dimension() - nbSCDOFs_u();     // nb total of skeletal DOFs for u (where global skeletal dofs of p start)
  std::transform(IT_sxgrad.begin(), IT_sxgrad.begin()+dim_p, it_IT_sys, [&offset](const size_t & index) { return index + offset; });
  IT_sys[dim_dofs-1] = sizeSystem()-1; // Lagrange multiplier
    
  // Creation of global DOFs for SC operator: IT_sc contains global cell dofs of u (offset to start at 0) and p (offset to start after those of u)
  std::vector<size_t> IT_sc(dim_sc, 0);  
  if (dim_sc>0){
    auto it_IT_sc = std::transform(IT_sxcurl.begin()+dim_u, IT_sxcurl.end(), IT_sc.begin(), [&offset](const size_t & index) { return index - offset; });
    size_t end_u = nbSCDOFs_u();
    offset = (m_sxgrad.dimension() - nbSCDOFs_p());
    std::transform(IT_sxgrad.begin()+dim_p, IT_sxgrad.end(), it_IT_sc, [&offset, &end_u](const size_t & index) { return index - offset + end_u; });
  }  

  return LocalStaticCondensation(Perm, IT_sys, IT_sc);
}


//------------------------------------------------------------------------------

void Stokes::_assemble_local_contribution(
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

std::vector<StokesNorms> Stokes::computeStokesNorms(const std::vector<Eigen::VectorXd> & list_dofs) const
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

std::pair<StokesNorms, StokesNorms> Stokes::computeContinuousErrorsNorms(
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


