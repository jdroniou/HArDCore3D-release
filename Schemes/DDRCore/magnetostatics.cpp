// Author: Daniele Di Pietro (daniele.di-pietro@umontpellier.fr)
#include <fstream>
#include <iomanip>
#include <thread>

#include <boost/math/constants/constants.hpp>

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#ifdef WITH_UMFPACK
#include <Eigen/UmfPackSupport>
#endif

#ifdef WITH_MKL
#include <Eigen/PardisoSupport>
#endif

#include <mesh.hpp>
#include <mesh_builder.hpp>

#include <xdiv.hpp>
#include <xcurl.hpp>

#include "magnetostatics.hpp"

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

  std::cout << FORMAT(25) << "[main] Mesh file" << mesh_file << std::endl;
  
  // Select the degree 
  size_t K = vm["degree"].as<size_t>();
  std::cout << FORMAT(25) << "[main] Degree" << K << std::endl;

  // Interpolate a vector function
  int solution = (vm.count("solution") ? vm["solution"].as<int>() : 0);
  Magnetostatics::ForcingTermType f;
  Magnetostatics::SolutionPotentialType u;
  Magnetostatics::SolutionCurlType sigma;

  switch (solution) {
  case 0:
    std::cout << "[main] Constant solution" << std::endl;    
    f = constant_f;
    u = constant_u;
    sigma = constant_sigma;
    break;

  case 1:
    std::cout << "[main] Linear solution" << std::endl;    
    f = linear_f;
    u = linear_u;
    sigma = linear_sigma;
    break;

  case 3:
    std::cout << "[main] Trigonometric solution" << std::endl;    
    f = trigonometric_f;
    u = trigonometric_u;
    sigma = trigonometric_sigma;
    break;

  default:
    std::cerr << "[main] ERROR: Unknown function" << std::endl;
    exit(1);    
  }

  // Build the mesh
  MeshBuilder meshbuilder = MeshBuilder(mesh_file, mesh_type);
  std::unique_ptr<Mesh> mesh_ptr = meshbuilder.build_the_mesh();

  // Create DDR core
  bool use_threads = (vm.count("pthread") ? vm["pthread"].as<bool>() : true);
  std::cout << "[main] " << (use_threads ? "Parallel execution" : "Sequential execution") << std:: endl;
  DDRCore ddr_core(*mesh_ptr, K, use_threads);

  // Assemble the problem
  Magnetostatics ms(ddr_core, use_threads);
  if(vm.count("stabilization-parameter")) {
    ms.stabilizationParameter() = vm["stabilization-parameter"].as<double>();
  }
  ms.assembleLinearSystem(f, u);

  // Export matrix if requested  
  if (vm.count("export-matrix")) {
    std::cout << "[main] Exporting matrix to Matrix Market format" << std::endl;
    saveMarket(ms.systemMatrix(), "A_magnetostatics.mtx");
    saveMarket(ms.systemVector(), "b_magnetostatics.mtx");
  }

  // Solve the problem
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
  
  // Compute the error in the energy norm
  Eigen::VectorXd uI = Eigen::VectorXd::Zero(ms.dimension());  
  uI.head(ms.xCurl().dimension()) = ms.xCurl().interpolate(sigma);
  uI.tail(ms.xDiv().dimension()) = ms.xDiv().interpolate(u);
  Eigen::VectorXd eh = uh - uI;
  double en_err = std::sqrt( eh.transpose() * ms.systemMatrix() * eh );
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

void Magnetostatics::assembleLinearSystem(const ForcingTermType & f, const SolutionPotentialType & u)
{
  // Assemble all local contributions
  auto assemble_all = [this, f, u](
				   size_t start,
				   size_t end,
				   std::list<Eigen::Triplet<double> > * my_triplets,
				   Eigen::VectorXd * my_rhs
				   )->void
		      {
			for (size_t iT = start; iT < end; iT++) {
			  this->_assemble_local_contribution(
							     iT,
							     this->_compute_local_contribution(iT, f, u),
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

    // Wait for the other thread to finish their task
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
Magnetostatics::_compute_local_contribution(size_t iT, const ForcingTermType & f, const SolutionPotentialType & u)
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
  
  // aT
  AT.topLeftCorner(dim_xcurl_T, dim_xcurl_T) = m_xcurl.computeL2Product(iT, m_stab_par);
  
  // bT
  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * m_ddrcore.degree());
  Eigen::MatrixXd mass_Pk3_T = compute_gram_matrix(evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polyk3, quad_2k_T), quad_2k_T);
  // Consistent term
  Eigen::MatrixXd BT
    = m_xdiv.cellOperators(iT).potential.transpose() * mass_Pk3_T * m_xcurl.cellOperators(iT).curl;

  // Stabilisation terms
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    Eigen::VectorXd nF = F.normal();
    double hF = F.diam();
    
    QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * m_ddrcore.degree());
    auto basis_Pk3_T_dot_nF_quad
      = scalar_product(evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polyk3, quad_2k_F), nF);
    auto basis_Pk_F_quad
      = evaluate_quad<Function>::compute(*m_ddrcore.faceBases(F.global_index()).Polyk, quad_2k_F);

    Eigen::MatrixXd CF = m_xcurl.extendOperator(T, F, m_xcurl.faceOperators(F.global_index()).curl);

    Eigen::MatrixXd gram_Pk3_T_dot_nF_Pk_F
      = compute_gram_matrix(basis_Pk3_T_dot_nF_quad, basis_Pk_F_quad, quad_2k_F);
    BT += m_stab_par * hF
      * m_xdiv.cellOperators(iT).potential.transpose() * (
							  compute_gram_matrix(basis_Pk3_T_dot_nF_quad, quad_2k_F) * m_xcurl.cellOperators(iT).curl
							  - gram_Pk3_T_dot_nF_Pk_F * CF
							  );
    BT.block(m_xdiv.localOffset(T, F), 0, m_xdiv.faceBases(F.global_index()).Polyk->dimension(), dim_xcurl_T)
      -= m_stab_par * hF * (
			    gram_Pk3_T_dot_nF_Pk_F.transpose() * m_xcurl.cellOperators(iT).curl
			    - compute_gram_matrix(basis_Pk_F_quad, quad_2k_F) * CF
			    );
  } // for iF
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
