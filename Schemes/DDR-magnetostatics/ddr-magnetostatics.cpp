// Author: Daniele Di Pietro (daniele.di-pietro@umontpellier.fr)
#include <fstream>
#include <iomanip>
#include <thread>

#include "ddr-magnetostatics.hpp"
#include <parallel_for.hpp>

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

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
  double t_wall_ddrcore = double(timer.elapsed().wall) * pow(10, -9);
  double t_proc_ddrcore = double(timer.elapsed().user + timer.elapsed().system) * pow(10, -9);
  std::cout << "[main] Time DDRCore (wall/proc) " << t_wall_ddrcore << "/" << t_proc_ddrcore << std::endl;

  // Assemble the problem
  timer.start();
  Magnetostatics ms(ddr_core, use_threads);
  if(vm.count("stabilization-parameter")) {
    ms.stabilizationParameter() = vm["stabilization-parameter"].as<double>();
  }
  ms.assembleLinearSystem(f, mu, u);
  timer.stop();
  double t_wall_model = double(timer.elapsed().wall) * pow(10, -9);
  double t_proc_model = double(timer.elapsed().user + timer.elapsed().system) * pow(10, -9);
  std::cout << "[main] Time model (wall/proc) " << t_wall_model << "/" << t_proc_model << std::endl;

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
  uh.head(ms.xCurl().dimension() - ms.nbSCDOFs()) = uh_condensed.head(ms.xCurl().dimension() - ms.nbSCDOFs()); 
  uh.segment(ms.xCurl().dimension() - ms.nbSCDOFs(), ms.nbSCDOFs()) = ms.scVector() + ms.scMatrix() * uh_condensed;
  uh.tail(ms.xDiv().dimension()) = uh_condensed.segment(ms.xCurl().dimension()-ms.nbSCDOFs(), ms.xDiv().dimension());

  timer.stop();
  double t_wall_solve = double(timer.elapsed().wall) * pow(10, -9);
  double t_proc_solve = double(timer.elapsed().user + timer.elapsed().system) * pow(10, -9);
  std::cout << "[main] Time solve (wall/proc) " << t_wall_solve << "/" << t_proc_solve << std::endl;
  
  // Compute the error in the energy norm
  Eigen::VectorXd uI = Eigen::VectorXd::Zero(ms.dimensionSpace());  
  uI.head(ms.xCurl().dimension()) = ms.xCurl().interpolate(sigma);
  uI.tail(ms.xDiv().dimension()) = ms.xDiv().interpolate(u);
  Eigen::VectorXd eh = uh - uI;
//////  double en_err = std::sqrt( eh.transpose() * ms.systemMatrix() * eh );  // no longer valid in statically condensed system
  // Error in Hcurl x Hdiv norm
  double hcurlhdiv_err = ms.computeNorm(eh) / ms.computeNorm(uI);
  double en_err = hcurlhdiv_err;
  std::cout << "[main] Hcurl-Hdiv error " << hcurlhdiv_err << std::endl;
  std::cout << "[main] Energy error " << en_err << std::endl;
  std::cout << "[main] Mesh diameter " << mesh_ptr->h_max() << std::endl;
 

  // Write results to file
  std::ofstream out("results.txt");
  out << "Solution: " << solution << std::endl;
  out << "Mesh: " << mesh_file << std::endl;
  out << "Degree: " << K << std::endl;
  out << "MeshSize: " << mesh_ptr->h_max() << std::endl;
  out << "NbCells: " << mesh_ptr->n_cells() << std::endl;
  out << "NbFaces: " << mesh_ptr->n_faces() << std::endl;
  out << "NbEdges: " << mesh_ptr->n_edges() << std::endl;
  out << "DimXCurl: " << ms.xCurl().dimension() << std::endl;
  out << "DimXDiv: " << ms.xDiv().dimension() << std::endl;
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
    m_xcurl(ddrcore, use_threads),
    m_xdiv(ddrcore, use_threads),
    m_A(sizeSystem(), sizeSystem()),
    m_b(Eigen::VectorXd::Zero(sizeSystem())),
    m_sc_A(nbSCDOFs(), sizeSystem()),
    m_sc_b(Eigen::VectorXd::Zero(nbSCDOFs())),
    m_stab_par(1.)
    
{
  m_output << "[Magnetostatics] Initializing" << std::endl;
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
                                       std::list<Eigen::Triplet<double> > * triplets_system,
                                       Eigen::VectorXd * rhs_system,
                                       std::list<Eigen::Triplet<double> > * triplets_sc,
                                       Eigen::VectorXd * rhs_sc
                                       )->void
                      {
                        for (size_t iT = start; iT < end; iT++) {
                          this->_assemble_local_contribution(
                                                             iT,
                                                             this->_compute_local_contribution(iT, f, mu, u),
                                                             *triplets_system,
                                                             *rhs_system,
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
    Eigen::VectorXd lF = -m_xcurl.faceOperators(F.global_index()).potential.transpose()
      * integrate(u_cross_nF, evaluate_quad<Function>::compute(*m_ddrcore.faceBases(F.global_index()).Polyk2, quad_2k_F), quad_2k_F);
    auto I_F = m_xcurl.globalDOFIndices(F);
    for (size_t i = 0; i < I_F.size(); i++) {
      m_b(I_F[i]) += lF(i);
    } // for i
  } // for iF
}

//------------------------------------------------------------------------------
 
std::pair<Eigen::MatrixXd, Eigen::VectorXd>
Magnetostatics::_compute_local_contribution(
                                            size_t iT, const ForcingTermType & f,
                                            const PermeabilityType & mu,
                                            const SolutionPotentialType & u
                                            )
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  size_t dim_xcurl_T = m_xcurl.dimensionCell(iT);
  size_t dim_xdiv_T = m_xdiv.dimensionCell(iT);
  size_t dim_T = dim_xcurl_T + dim_xdiv_T;
  
  Eigen::MatrixXd AT = Eigen::MatrixXd::Zero(dim_T, dim_T);
  Eigen::VectorXd lT = Eigen::VectorXd::Zero(dim_T);

  //------------------------------------------------------------------------------
  // Local matrix
  //------------------------------------------------------------------------------

  // Mass matrix for (P^k(T))^3  
  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * m_ddrcore.degree());
  Eigen::MatrixXd mass_Pk3_T = compute_gram_matrix(evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polyk3, quad_2k_T), quad_2k_T);

  // aT
  AT.topLeftCorner(dim_xcurl_T, dim_xcurl_T) = m_xcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T, mu);
  
  // bT
  Eigen::MatrixXd BT = m_xdiv.computeL2ProductCurl(iT, m_xcurl, "right", m_stab_par, mass_Pk3_T);
  AT.topRightCorner(dim_xcurl_T, dim_xdiv_T) -= BT.transpose();
  AT.bottomLeftCorner(dim_xdiv_T, dim_xcurl_T) += BT;

  // cT
  AT.bottomRightCorner(dim_xdiv_T, dim_xdiv_T)
    += m_xdiv.cellOperators(iT).divergence_rhs.transpose() * m_xdiv.cellOperators(iT).divergence;

  //------------------------------------------------------------------------------
  // Local vector
  //------------------------------------------------------------------------------

  lT.tail(dim_xdiv_T) = m_xdiv.cellOperators(iT).potential.transpose()
    * integrate(f, evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polyk3, quad_2k_T), quad_2k_T);
  
  return std::make_pair(AT, lT);
}

//------------------------------------------------------------------------------

void Magnetostatics::_assemble_local_contribution(
                                                  size_t iT,
                                                  const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT,
                                                  std::list<Eigen::Triplet<double> > & triplets_system,
                                                  Eigen::VectorXd & rhs_system,
                                                  std::list<Eigen::Triplet<double> > & triplets_sc,
                                                  Eigen::VectorXd & rhs_sc
                                                  )
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  // AT and bT are made of 3 components: H skeleton (face+edges), H cells, and A
  // We want to condense the 2nd component H cell. We first transform AT and bT with a permutation to put
  // all these unknowns in the bottom right corner
  size_t dim1 = m_xcurl.localOffset(T);                // dimension of skeletal unknowns for H
  size_t dim_sc = m_xcurl.numLocalDofsCell();    // dimension of cell unknowns for H (sc DOFs)
  size_t dim3 = m_xdiv.dimensionCell(iT);         // dimension of A
  size_t dim_dofs = dim1+dim3;      // nb of dofs remaining after SC
  // Permutation matrix
  Eigen::MatrixXd Perm = Eigen::MatrixXd::Zero(dim_dofs+dim_sc, dim_dofs+dim_sc);
  Perm.topLeftCorner(dim1, dim1) = Eigen::MatrixXd::Identity(dim1, dim1);
  Perm.block(dim1+dim3, dim1, dim_sc, dim_sc) = Eigen::MatrixXd::Identity(dim_sc, dim_sc);
  Perm.block(dim1, dim1+dim_sc, dim3, dim3) = Eigen::MatrixXd::Identity(dim3, dim3);
  // Permuted local matrix and rhs, the velocity cell unknowns are at the end
  Eigen::MatrixXd AT = Perm * lsT.first * Perm.transpose();
  Eigen::VectorXd bT = Perm * lsT.second;

  // Extract 4 blocks of AT, bT for static condensation
  Eigen::MatrixXd A11 = AT.topLeftCorner(dim_dofs, dim_dofs);  
  Eigen::MatrixXd A12 = AT.topRightCorner(dim_dofs, dim_sc);
  Eigen::MatrixXd A21 = AT.bottomLeftCorner(dim_sc, dim_dofs);
  Eigen::MatrixXd A22 = AT.bottomRightCorner(dim_sc, dim_sc);
  Eigen::VectorXd b1 = bT.head(dim_dofs);
  Eigen::VectorXd b2 = bT.tail(dim_sc);
  
  // Create condensed system (AT_reduced, bT_reduced) and SC recovery operator (RT, cT)
  Eigen::MatrixXd A22inv_A21 = A22.ldlt().solve(A21);
  Eigen::VectorXd A22inv_b2 = A22.ldlt().solve(b2);
  Eigen::MatrixXd AT_reduced = A11 - A12 * A22inv_A21;
  Eigen::VectorXd bT_reduced = b1 - A12 * A22inv_b2;
  Eigen::MatrixXd RT = -A22inv_A21;
  Eigen::VectorXd cT = A22inv_b2;

  // STATICALLY CONDENSED SYSTEM
  //
  // Create the vector of global DOF indices for the SC system: Idofs_T contains the skeletal dofs of u, dofs of p, and lambda
  std::vector<size_t> Idofs_T(dim_dofs, 0);
  auto I_xcurl_T = m_xcurl.globalDOFIndices(T);
  auto I_xgrad_T = m_xdiv.globalDOFIndices(T);
  auto it_Idofs_T = std::copy(I_xcurl_T.begin(), I_xcurl_T.begin()+dim1, Idofs_T.begin()); // put skeletal DOFs of H in Idofs_T
  size_t offset = m_xcurl.dimension() - nbSCDOFs();     // nb total of skeletal DOFs for H (where global dofs of A start)
  std::transform(I_xgrad_T.begin(), I_xgrad_T.end(), it_Idofs_T, [&offset](const size_t & index) { return index + offset; });
  
  // Assemble statically condensed system  
  for (size_t i = 0; i < dim_dofs; i++){
    rhs_system(Idofs_T[i]) += bT_reduced(i);
    for (size_t j = 0; j < dim_dofs; j++){    
      triplets_system.emplace_back(Idofs_T[i], Idofs_T[j], AT_reduced(i,j));
    }
  }

  // RECOVERY OPERATOR
  //
  // Create the vector of DOF indices: Isc_T contains global cell dofs of H (offset to start at 0)
  std::vector<size_t> Isc_T(dim_sc);
  std::transform(I_xcurl_T.begin()+dim1, I_xcurl_T.end(), Isc_T.begin(), [&offset](const size_t & index) { return index - offset; });
 
  // Assemble recovery operator
  for (size_t i = 0; i < dim_sc; i++){
    rhs_sc(Isc_T[i]) += cT(i);
    for (size_t j = 0; j < dim_dofs; j++){
      triplets_sc.emplace_back(Isc_T[i], Idofs_T[j], RT(i,j));
    }
  }

}


//------------------------------------------------------------------------------

double Magnetostatics::computeNorm( const Eigen::VectorXd & v ) const
{
  Eigen::VectorXd local_sqnorms = Eigen::VectorXd::Zero(m_ddrcore.mesh().n_cells());

  // Xcurl correspond to the first components of v, Xdiv to the last ones
  Eigen::VectorXd v_curl = v.head(m_xcurl.dimension());
  Eigen::VectorXd v_div = v.tail(m_xdiv.dimension());
  
  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &v_curl, &v_div, &local_sqnorms](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *m_ddrcore.mesh().cell(iT);
        // Mass matrix for (P^k(T))^3
        QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2*m_ddrcore.degree() );
        Eigen::MatrixXd mass_Pk3_T = compute_gram_matrix(evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polyk3, quad_2k_T), quad_2k_T);
        Eigen::VectorXd v_curl_T = m_xcurl.restrict(T, v_curl);
        Eigen::VectorXd v_div_T = m_xdiv.restrict(T, v_div);

        // Contribution of L2 norms, without any weight (no permeability)
        local_sqnorms(iT) += v_curl_T.transpose() * m_xcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T) * v_curl_T;
        local_sqnorms(iT) += v_div_T.transpose() * m_xdiv.computeL2Product(iT, m_stab_par, mass_Pk3_T) * v_div_T;

        // Contribution of L2 norms of curl and div
        local_sqnorms(iT) += v_curl_T.transpose() * m_xdiv.computeL2ProductCurl(iT, m_xcurl, "both", m_stab_par, mass_Pk3_T) * v_curl_T;
        Eigen::MatrixXd DT = m_xdiv.cellOperators(iT).divergence;
        local_sqnorms(iT) += v_div_T.transpose() * DT.transpose() * compute_gram_matrix(evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polyk, quad_2k_T), quad_2k_T) * DT * v_div_T;
      }
    };
  parallel_for(m_ddrcore.mesh().n_cells(), compute_local_squarednorms, m_use_threads);
  
  return std::sqrt(std::abs(local_sqnorms.sum()));
}


