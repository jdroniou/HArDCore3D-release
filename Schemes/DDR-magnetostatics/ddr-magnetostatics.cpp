// Author: Daniele Di Pietro (daniele.di-pietro@umontpellier.fr)
#include <fstream>
#include <iomanip>
#include <thread>

#include "ddr-magnetostatics.hpp"

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

const std::string mesh_dir = "../HArDCore3D/meshes/";
std::string default_mesh = mesh_dir + "Voro-small-0/RF_fmt/voro-2";
std::string default_meshtype = "RF";

//------------------------------------------------------------------------------

int main(int argc, const char* argv[])
{
  // Program options
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Display this help message")
    ("mesh,m", boost::program_options::value<std::string>(), "Set the mesh")
    ("meshtype,t", boost::program_options::value<std::string>(), "Set the mesh type (TG,MSH,RF)")
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
  std::string mesh_type = (vm.count("meshtype") ? vm["meshtype"].as<std::string>() : default_meshtype);

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
  MeshBuilder meshbuilder = MeshBuilder(mesh_file, mesh_type);
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
  Eigen::VectorXd uh;
  if (vm.count("iterative-solver")) {
    std::cout << "[main] Solving the linear system using BiCGSTAB" << std::endl;
    
    Eigen::BiCGSTAB<Magnetostatics::SystemMatrixType, Eigen::IncompleteLUT<double> > solver;
    // solver.preconditioner().setFillfactor(2);
    solver.compute(ms.systemMatrix());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
      exit(1);
    }
    uh = solver.solve(ms.systemVector());
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
    uh = solver.solve(ms.systemVector());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not solve linear system" << std::endl;
    }
  }
  timer.stop();
  double t_wall_solve = double(timer.elapsed().wall) * pow(10, -9);
  double t_proc_solve = double(timer.elapsed().user + timer.elapsed().system) * pow(10, -9);
  std::cout << "[main] Time solve (wall/proc) " << t_wall_solve << "/" << t_proc_solve << std::endl;
  
  // Compute the error in the energy norm
  Eigen::VectorXd uI = Eigen::VectorXd::Zero(ms.dimension());  
  uI.head(ms.xCurl().dimension()) = ms.xCurl().interpolate(sigma);
  uI.tail(ms.xDiv().dimension()) = ms.xDiv().interpolate(u);
  Eigen::VectorXd eh = uh - uI;
  double en_err = std::sqrt( eh.transpose() * ms.systemMatrix() * eh );
  // Error in Hcurl x Hdiv norm
  double hcurlhdiv_err = ms.computeNorm(eh) / ms.computeNorm(uI);
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
    m_A(dimension(), dimension()),
    m_b(Eigen::VectorXd::Zero(dimension())),
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
                                       std::list<Eigen::Triplet<double> > * my_triplets,
                                       Eigen::VectorXd * my_rhs
                                       )->void
                      {
                        for (size_t iT = start; iT < end; iT++) {
                          this->_assemble_local_contribution(
                                                             iT,
                                                             this->_compute_local_contribution(iT, f, mu, u),
                                                             *my_triplets,
                                                             *my_rhs
                                                             );
                        } // for iT
                      };
  
  if (m_use_threads) {
    m_output << "[Magnetostatics] Parallel assembly" << std::endl;
   
    // Select the number of threads
    unsigned nb_threads_hint = std::thread::hardware_concurrency();
    unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);

    // Compute the batch size and the remainder
    unsigned nb_elements = m_ddrcore.mesh().n_cells();
    unsigned batch_size = nb_elements / nb_threads;
    unsigned batch_remainder = nb_elements % nb_threads;
 
    // Create vectors of triplets and vectors
    std::vector<std::list<Eigen::Triplet<double> > > triplets(nb_threads + 1);
    std::vector<Eigen::VectorXd> rhs(nb_threads + 1);

    for (unsigned i = 0; i < nb_threads + 1; i++) {
      rhs[i] = Eigen::VectorXd::Zero(this->dimension());
    } // for i

    // Assign a task to each thread
    std::vector<std::thread> my_threads(nb_threads);
    for (unsigned i = 0; i < nb_threads; ++i) {
      int start = i * batch_size;
      my_threads[i] = std::thread(assemble_all, start, start + batch_size, &triplets[i], &rhs[i]);
    }

    // Execute the elements left
    int start = nb_threads * batch_size;
    assemble_all(start, start + batch_remainder, &triplets[nb_threads], &rhs[nb_threads]);

    // Wait for the other threads to finish their task
    std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));

    // Create matrix from triplets
    size_t n_triplets = 0;
    for (auto triplets_thread : triplets) {
      n_triplets += triplets_thread.size();
    }
    std::vector<Eigen::Triplet<double> > all_triplets(n_triplets);
    auto triplet_index = all_triplets.begin();
    for (auto triplets_thread : triplets) {
      triplet_index = std::copy(triplets_thread.begin(), triplets_thread.end(), triplet_index);
    }
    m_A.setFromTriplets(all_triplets.begin(), all_triplets.end());
    for (auto rhs_thread : rhs) {
      m_b += rhs_thread;
    }
  } else {
    m_output << "[Magnetostatics] Sequential assembly" << std::endl;
    std::list<Eigen::Triplet<double> > triplets;
    assemble_all(0, m_ddrcore.mesh().n_cells(), &triplets, &m_b);
    m_A.setFromTriplets(triplets.begin(), triplets.end());
  }
    
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
                                                  std::list<Eigen::Triplet<double> > & my_triplets,
                                                  Eigen::VectorXd & my_rhs
                                                  )
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  size_t dim_T = m_xcurl.dimensionCell(iT) + m_xdiv.dimensionCell(iT);
  size_t dim_xcurl = m_xcurl.dimension();

  // Create the vector of DOF indices
  auto I_xcurl_T = m_xcurl.globalDOFIndices(T);
  auto I_xdiv_T = m_xdiv.globalDOFIndices(T);
  std::vector<size_t> I_T(dim_T);
  auto it_I_T = std::copy(I_xcurl_T.begin(), I_xcurl_T.end(), I_T.begin());
  std::transform(I_xdiv_T.begin(), I_xdiv_T.end(), it_I_T, [dim_xcurl](const size_t & index) { return index + dim_xcurl; });

  // Assemble
  const Eigen::MatrixXd & AT = lsT.first;
  const Eigen::VectorXd & bT = lsT.second;
  for (size_t i = 0; i < dim_T; i++) {
    my_rhs(I_T[i]) += bT(i);
    for (size_t j = 0; j < dim_T; j++) {
      my_triplets.push_back( Eigen::Triplet<double>(I_T[i], I_T[j], AT(i,j)) );
    } // for j
  } // for i
}


//------------------------------------------------------------------------------

double Magnetostatics::computeNorm( const Eigen::VectorXd & v ) const
{
  double val = 0.;

  // Xcurl correspond to the first components of v, Xdiv to the last ones
  Eigen::VectorXd v_curl = v.head(m_xcurl.dimension());
  Eigen::VectorXd v_div = v.tail(m_xdiv.dimension());
  
  // Loop over cells
  for (size_t iT = 0; iT < m_ddrcore.mesh().n_cells(); iT++){
    Cell & T = *m_ddrcore.mesh().cell(iT);
    // Mass matrix for (P^k(T))^3
    QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2*m_ddrcore.degree() );
    Eigen::MatrixXd mass_Pk3_T = compute_gram_matrix(evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polyk3, quad_2k_T), quad_2k_T);
    Eigen::VectorXd v_curl_T = m_xcurl.restrict(T, v_curl);
    Eigen::VectorXd v_div_T = m_xdiv.restrict(T, v_div);

    // Contribution of L2 norms, without any weight (no permeability)
    val += v_curl_T.transpose() * m_xcurl.computeL2Product(iT, m_stab_par, mass_Pk3_T) * v_curl_T;
    val += v_div_T.transpose() * m_xdiv.computeL2Product(iT, m_stab_par, mass_Pk3_T) * v_div_T;

    // Contribution of L2 norms of curl and div
    val += v_curl_T.transpose() * m_xdiv.computeL2ProductCurl(iT, m_xcurl, "both", m_stab_par, mass_Pk3_T) * v_curl_T;
    Eigen::MatrixXd DT = m_xdiv.cellOperators(iT).divergence;
    val += v_div_T.transpose() * DT.transpose() * compute_gram_matrix(evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polyk, quad_2k_T), quad_2k_T) * DT * v_div_T;

  }

  return std::sqrt(val);
}


