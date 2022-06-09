// Author: Jerome Droniou (jerome.droniou@monash.edu)
#include <fstream>
#include <iomanip>
#include <thread>

#include "sddr-magnetostatics.hpp"
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>

#include <boost/program_options.hpp>
#include <display_timer.hpp>

#ifdef WITH_UMFPACK
#include <Eigen/UmfPackSupport>
#endif

#ifdef WITH_MKL
#include <Eigen/PardisoSupport>
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
    ("solution,s", boost::program_options::value<int>()->default_value(0), "Select the solution")
    ("export-matrix,e", "Export matrix to Matrix Market format")
    ("iterative-solver,i", "Use iterative linear solver")
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
  Magnetostatics::ForcingTermType f;
  Magnetostatics::PermeabilityType mu(1.);    // Default permeability constant equal to 1.
  Magnetostatics::SolutionPotentialType u;
  Magnetostatics::SolutionCurlType sigma;

  switch (solution) {
  case 0:
    std::cout << "[main] Constant solution" << std::endl;    
    f = constant_f;
    mu = constant_mu;
    u = constant_u;
    sigma = constant_sigma;
    break;

  case 1:
    std::cout << "[main] Linear solution" << std::endl;    
    f = linear_f;
    u = linear_u;
    mu = linear_mu;
    sigma = linear_sigma;
    break;

  case 3:
    std::cout << "[main] Trigonometric solution" << std::endl;    
    f = trigonometric_f;
    u = trigonometric_u;
    mu = trigonometric_mu;
    sigma = trigonometric_sigma;
    break;

  case 4:
    std::cout << "[main] Variable permeability solution" << std::endl;
    f = variable_permeability_f;
    u = variable_permeability_u;
    mu = variable_permeability_mu;
    sigma = variable_permeability_sigma;
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
  Magnetostatics ms(ddr_core, use_threads);
  if(vm.count("stabilization-parameter")) {
    ms.stabilizationParameter() = vm["stabilization-parameter"].as<double>();
  }
  ms.assembleLinearSystem(f, mu, u);
  timer.stop();
  double t_wall_model, t_proc_model;
  std::tie(t_wall_model, t_proc_model) = store_times(timer, "[main] Time model (wall/proc) ");

  // Export matrix if requested  
  if (vm.count("export-matrix")) {
    std::cout << "[main] Exporting matrix to Matrix Market format" << std::endl;
    saveMarket(ms.systemMatrix(), "A_magnetostatics.mtx");
    saveMarket(ms.systemVector(), "b_magnetostatics.mtx");
  }

  // Solve the problem
  timer.start();
  Eigen::VectorXd uh_condensed;
  if (vm.count("iterative-solver")) {
    std::cout << "[main] Solving the linear system using BiCGSTAB" << std::endl;
    
    Eigen::BiCGSTAB<Magnetostatics::SystemMatrixType, Eigen::IncompleteLUT<double> > solver;
    // solver.preconditioner().setFillfactor(2);
    solver.compute(ms.systemMatrix());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
      exit(1);
    }
    uh_condensed = solver.solve(ms.systemVector());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not solve direct system" << std::endl;
      exit(1);
    }
  } else { 
#ifdef WITH_MKL
    std::cout << "[main] Solving the linear system using Pardiso" << std::endl;    
    Eigen::PardisoLU<Magnetostatics::SystemMatrixType> solver;
#elif WITH_UMFPACK
    std::cout << "[main] Solving the linear system using Umfpack" << std::endl;    
    Eigen::UmfPackLU<Magnetostatics::SystemMatrixType> solver;
#else
    std::cout << "[main] Solving the linear system using direct solver" << std::endl;    
    Eigen::SparseLU<Magnetostatics::SystemMatrixType> solver;
#endif
    solver.compute(ms.systemMatrix());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
    }
    uh_condensed = solver.solve(ms.systemVector());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not solve linear system" << std::endl;
    }
  }
  // Re-create statically condensed unknowns
  Eigen::VectorXd uh = Eigen::VectorXd::Zero(ms.dimensionSpace());
  uh.head(ms.sxCurl().dimension() - ms.nbSCDOFs()) = uh_condensed.head(ms.sxCurl().dimension() - ms.nbSCDOFs()); 
  uh.segment(ms.sxCurl().dimension() - ms.nbSCDOFs(), ms.nbSCDOFs()) = ms.scVector() + ms.scMatrix() * uh_condensed;
  uh.tail(ms.sxDiv().dimension()) = uh_condensed.segment(ms.sxCurl().dimension()-ms.nbSCDOFs(), ms.sxDiv().dimension());

  timer.stop();
  double t_wall_solve, t_proc_solve;
  std::tie(t_wall_solve, t_proc_solve) = store_times(timer, "[main] Time solve (wall/proc) ");
    
  // Compute the error in the energy norm
  Eigen::VectorXd uI = Eigen::VectorXd::Zero(ms.dimensionSpace());  
  uI.head(ms.sxCurl().dimension()) = ms.sxCurl().interpolate(sigma);
  uI.tail(ms.sxDiv().dimension()) = ms.sxDiv().interpolate(u);
  Eigen::VectorXd eh = uh - uI;
  // Error in Hcurl x Hdiv norm
  std::vector<double> list_norms = ms.computeNorms(std::vector<Eigen::VectorXd> {eh, uI});
  double error = list_norms[0];
  double norm = list_norms[1];
  double hcurlhdiv_err = error / norm;
  double en_err = hcurlhdiv_err;
  std::cout << "[main] Hcurl-Hdiv error " << hcurlhdiv_err << std::endl;
  std::cout << "[main] Mesh diameter " << mesh_ptr->h_max() << std::endl;
 

  // Write results to file
  std::ofstream out("results.txt");
  out << "Scheme: sddr-magnetostatic" << std::endl;
  out << "Solution: " << solution << std::endl;
  out << "Mesh: " << mesh_file << std::endl;
  out << "Degree: " << K << std::endl;
  out << "MeshSize: " << mesh_ptr->h_max() << std::endl;
  out << "NbCells: " << mesh_ptr->n_cells() << std::endl;
  out << "NbFaces: " << mesh_ptr->n_faces() << std::endl;
  out << "NbEdges: " << mesh_ptr->n_edges() << std::endl;
  out << "DimSXCurl: " << ms.sxCurl().dimension() << std::endl;
  out << "DimSXDiv: " << ms.sxDiv().dimension() << std::endl;
  out << "SizeSCsystem: " << ms.sizeSystem() << std::endl;
  out << "EnergyError: " << en_err << std::endl;  
  out << "HcurlHdivError: " << hcurlhdiv_err << std::endl;  
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
// Magnetostatics
//------------------------------------------------------------------------------

Magnetostatics::Magnetostatics(
                               const DDRCore & ddrcore,
                               bool use_threads,
                               std::ostream & output
                               )
  : m_ddrcore(ddrcore),
    m_use_threads(use_threads),
    m_output(output),
    m_ser_pro(ddrcore, use_threads, output),
    m_sxcurl(ddrcore, m_ser_pro, use_threads),
    m_sxdiv(ddrcore, use_threads),
    m_nloc_sc_H(m_ser_pro.nDOFs_cells_SXCurl()),
    m_A(sizeSystem(), sizeSystem()),
    m_b(Eigen::VectorXd::Zero(sizeSystem())),
    m_sc_A(nbSCDOFs(), sizeSystem()),
    m_sc_b(Eigen::VectorXd::Zero(nbSCDOFs())),
    m_stab_par(1.)
    
{
  m_output << "[Magnetostatics] Initializing" << std::endl;
  
  // To avoid performing static condensation, initialise m_nloc_sc_H at 0
}

//------------------------------------------------------------------------------

void Magnetostatics::assembleLinearSystem(
                                          const ForcingTermType & f,
                                          const PermeabilityType & mu,
                                          const SolutionPotentialType & u
                                          )
{
  // Assemble all local contributions
  auto assemble_all = [this, f, mu, u](
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
                                                             this->_compute_local_contribution(iT, f, mu),
                                                             *triplets_sys,
                                                             *rhs_sys,
                                                             *triplets_sc,
                                                             *rhs_sc
                                                             );
                        } // for iT
                      };
                      
  // Assemble the matrix and rhs
  if (m_use_threads) {
    m_output << "[Magnetostatics] Parallel assembly" << std::endl;
  }else{
    m_output << "[Magnetostatics] Sequential assembly" << std::endl;
  }
  std::tie(m_A, m_b, m_sc_A, m_sc_b) = parallel_assembly_system(m_ddrcore.mesh().n_cells(), this->sizeSystem(), std::make_pair(this->nbSCDOFs(), this->sizeSystem()), this->nbSCDOFs(), assemble_all, m_use_threads);
    
  // Assemble boundary conditions
  for (auto iF : m_ddrcore.mesh().get_b_faces()) {
    const Face & F = *iF;

    // Unit normal vector to F pointing out of the domain
    const Cell & TF = *F.cell(0);
    Eigen::Vector3d nF = TF.face_normal(TF.index_face(&F));
    
    std::function<Eigen::Vector3d(const Eigen::Vector3d &)> u_cross_nF
      = [&u, &nF](const Eigen::Vector3d & x) {
          return u(x).cross(nF);
        };
    
    QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * m_ddrcore.degree());    
    Eigen::VectorXd lF = -m_sxcurl.facePotential(F).transpose()
      * integrate(u_cross_nF, evaluate_quad<Function>::compute(*m_ddrcore.faceBases(F.global_index()).Polyk2, quad_2k_F), quad_2k_F);
    auto I_F = m_sxcurl.globalDOFIndices(F);
    for (size_t i = 0; i < I_F.size(); i++) {
      m_b(I_F[i]) += lF(i);
    } // for i
  } // for iF
}

//------------------------------------------------------------------------------
 
std::pair<Eigen::MatrixXd, Eigen::VectorXd>
Magnetostatics::_compute_local_contribution(
                                            size_t iT, const ForcingTermType & f,
                                            const PermeabilityType & mu
                                            )
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  size_t dim_sxcurl_T = m_sxcurl.dimensionCell(iT);
  size_t dim_sxdiv_T = m_sxdiv.dimensionCell(iT);
  size_t dim_T = dim_sxcurl_T + dim_sxdiv_T;
  
  Eigen::MatrixXd AT = Eigen::MatrixXd::Zero(dim_T, dim_T);
  Eigen::VectorXd lT = Eigen::VectorXd::Zero(dim_T);

  //------------------------------------------------------------------------------
  // Local matrix
  //------------------------------------------------------------------------------

  // Mass matrix for (P^k(T))^3
  MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_ddrcore.degree());  
  Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_2k);

  // aT
  AT.topLeftCorner(dim_sxcurl_T, dim_sxcurl_T) = m_sxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T, mu);
  
  // bT
  Eigen::MatrixXd BT = m_sxdiv.computeL2ProductCurl(iT, m_sxcurl, "right", m_stab_par, mass_Pk3_T);
  AT.topRightCorner(dim_sxcurl_T, dim_sxdiv_T) -= BT.transpose();
  AT.bottomLeftCorner(dim_sxdiv_T, dim_sxcurl_T) += BT;

  // cT
  AT.bottomRightCorner(dim_sxdiv_T, dim_sxdiv_T)
    += m_sxdiv.cellOperators(iT).divergence_rhs.transpose() * m_sxdiv.cellOperators(iT).divergence;

  //------------------------------------------------------------------------------
  // Local vector
  //------------------------------------------------------------------------------

  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * m_ddrcore.degree());
  lT.tail(dim_sxdiv_T) = m_sxdiv.cellOperators(iT).potential.transpose()
    * integrate(f, evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polyk3, quad_2k_T), quad_2k_T);
  
  return std::make_pair(AT, lT);
}

//------------------------------------------------------------------------------

LocalStaticCondensation Magnetostatics::_compute_static_condensation(const size_t & iT) const
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  // Dimensions
  size_t dim_H = m_sxcurl.dimensionCell(iT) - m_nloc_sc_H(iT);      // dimension of skeletal unknowns for H
  size_t dim_A = m_sxdiv.dimensionCell(iT);         // dimension of A
  size_t dim_dofs = dim_H+dim_A;      // nb of dofs remaining after SC

  // Creation of permutation matrix
  Eigen::MatrixXd Perm = Eigen::MatrixXd::Zero(dim_dofs+m_nloc_sc_H(iT), dim_dofs+m_nloc_sc_H(iT));
  Perm.topLeftCorner(dim_H, dim_H) = Eigen::MatrixXd::Identity(dim_H, dim_H);
  Perm.block(dim_H+dim_A, dim_H, m_nloc_sc_H(iT), m_nloc_sc_H(iT)) = Eigen::MatrixXd::Identity(m_nloc_sc_H(iT), m_nloc_sc_H(iT));
  Perm.block(dim_H, dim_H+m_nloc_sc_H(iT), dim_A, dim_A) = Eigen::MatrixXd::Identity(dim_A, dim_A);

  // Creation of global DOFs for system: IT_sys contains the skeletal dofs of H, and dofs of A
  std::vector<size_t> IT_sys(dim_dofs, 0);
  auto IT_sxcurl = m_sxcurl.globalDOFIndices(T);
  auto IT_sxdiv = m_sxdiv.globalDOFIndices(T);
  auto it_IT_sys = std::copy(IT_sxcurl.begin(), IT_sxcurl.begin()+dim_H, IT_sys.begin()); // put skeletal DOFs of H in IT_sys
  size_t offset = m_sxcurl.dimension() - nbSCDOFs();     // nb total of skeletal DOFs for H (where global dofs of A start)
  std::transform(IT_sxdiv.begin(), IT_sxdiv.end(), it_IT_sys, [&offset](const size_t & index) { return index + offset; });
  
  // Creation of global DOFs for SC operator: IT_sc contains global cell dofs of H (offset to start at 0)
  std::vector<size_t> IT_sc(m_nloc_sc_H(iT));
  std::transform(IT_sxcurl.begin()+dim_H, IT_sxcurl.end(), IT_sc.begin(), [&offset](const size_t & index) { return index - offset; });
 
  return LocalStaticCondensation(Perm, IT_sys, IT_sc);
}

//------------------------------------------------------------------------------

void Magnetostatics::_assemble_local_contribution(
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
  std::vector<size_t> IT_sys = locSC.globalDOFs_sys();
  for (size_t i = 0; i < locSC.dim_sys(); i++){
    rhs_sys(IT_sys[i]) += bT_sys(i);
    for (size_t j = 0; j < locSC.dim_sys(); j++){    
      triplets_sys.emplace_back(IT_sys[i], IT_sys[j], AT_sys(i,j));
    }
  }

  // RECOVERY OPERATOR
  std::vector<size_t> IT_sc = locSC.globalDOFs_sc();
  for (size_t i = 0; i < locSC.dim_sc(); i++){
    rhs_sc(IT_sc[i]) += bT_sc(i);
    for (size_t j = 0; j < locSC.dim_sys(); j++){
      triplets_sc.emplace_back(IT_sc[i], IT_sys[j], AT_sc(i,j));
    }
  }

}


//------------------------------------------------------------------------------

std::vector<double> Magnetostatics::computeNorms( const std::vector<Eigen::VectorXd> & list_dofs ) const
{
  size_t nb_vectors = list_dofs.size();
  std::vector<Eigen::VectorXd> local_sqnorms(nb_vectors, Eigen::VectorXd::Zero(m_ddrcore.mesh().n_cells()));

  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &list_dofs, &local_sqnorms, &nb_vectors](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *m_ddrcore.mesh().cell(iT);
        // Mass matrix for (P^k(T))^3
	      MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_ddrcore.degree());
        Eigen::MatrixXd mass_Pk3_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk3, int_mono_2k);
        Eigen::MatrixXd mass_Pk_T = GramMatrix(T, *m_ddrcore.cellBases(iT).Polyk, int_mono_2k);
        Eigen::MatrixXd L2curl_T = m_sxcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T);
        Eigen::MatrixXd L2div_T = m_sxdiv.computeL2Product(iT, m_stab_par, mass_Pk3_T);
        Eigen::MatrixXd L2div_CC_T = m_sxdiv.computeL2ProductCurl(iT, m_sxcurl, "both", m_stab_par, mass_Pk3_T);
        Eigen::MatrixXd DT = m_sxdiv.cellOperators(iT).divergence;

        for (size_t i=0; i<nb_vectors; i++){
          Eigen::VectorXd v_curl_T = m_sxcurl.restrict(T, list_dofs[i].head(m_sxcurl.dimension()));
          Eigen::VectorXd v_div_T = m_sxdiv.restrict(T, list_dofs[i].tail(m_sxdiv.dimension()));

          // Contribution of L2 norms, without any weight (no permeability)
          local_sqnorms[i](iT) += v_curl_T.transpose() * L2curl_T * v_curl_T;
          local_sqnorms[i](iT) += v_div_T.transpose() * L2div_T * v_div_T;

          // Contribution of L2 norms of curl and div
          local_sqnorms[i](iT) += v_curl_T.transpose() * L2div_CC_T * v_curl_T;
          local_sqnorms[i](iT) += v_div_T.transpose() * DT.transpose() * mass_Pk_T * DT * v_div_T;
        }
      }
    };
  parallel_for(m_ddrcore.mesh().n_cells(), compute_local_squarednorms, m_use_threads);
  
  // Vector of outputs
  std::vector<double> list_norms(nb_vectors, 0.);
  for (size_t i=0; i<nb_vectors; i++){
    list_norms[i] = std::sqrt(std::abs(local_sqnorms[i].sum()));
  }
  
  return list_norms;
}


