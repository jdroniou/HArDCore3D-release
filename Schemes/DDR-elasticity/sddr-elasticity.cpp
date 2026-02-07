// Author: Jerome Droniou (jerome.droniou@umontpellier.fr)

#include <fstream>
#include <iomanip>
#include <thread>

#include "sddr-elasticity.hpp"
#include <parallel_for.hpp>
#include <max_degrees_quadratures.hpp>
#include <GMpoly_cell.hpp>

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
    ("pthread,p", boost::program_options::value<bool>()->default_value(true), "Use thread-based parallelism")
    ("solution,s", boost::program_options::value<int>()->default_value(0), "Select the solution")
    ("export-matrix,e", "Export matrix to Matrix Market format")
    ("boundary-conditions,b", boost::program_options::value<std::string>()->default_value("D"), "Boundary conditions (D for essential, N for natural)")
    ("mu", boost::program_options::value<double>()->default_value(1.), "Select the second Lamé coefficient")
    ("lambda", boost::program_options::value<double>()->default_value(1.), "Select the first Lamé coefficient")
    ("plot", boost::program_options::value<std::string>(), "Save plot of the solution to the given filename")
    ("solver", boost::program_options::value<std::string>()->default_value("PardisoLU"), "Choice of solver, not case dependent. Options are: PardisoLU, UMFPACK, PaStiXLU, PaStiXLLT, EigenLU, EigenBiCGSTAB (reverts to EigenLU if the selected solver is not available)")
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
  Elasticity::GradientDisplacementType grad_u;
  
  double mu = vm["mu"].as<double>();
  double lambda = vm["lambda"].as<double>();

  std::cout << FORMAT(25) << "[main] Lamé coefficients: mu= " << mu << ", lambda= " << lambda << std::endl;

  switch (solution) {
  case 0:
    std::cout << "[main] Test case: trigonometric solution" << std::endl;
    f = trigonometric_f;
    u = trigonometric_u;
    grad_u = trigonometric_grad_u;
    break;

  case 1:
    std::cout << "[main] Test case: constant solution" << std::endl;
    f = constant_f;
    u = constant_u;
    grad_u = constant_grad_u;
    break;

  case 2:
    std::cout << "[main] Test case: linear solution" << std::endl;
    f = linear_f;
    u = linear_u;
    grad_u = linear_grad_u;
    break;

  case 3:
    std::cout << "[main] Test case: quadratic solution" << std::endl;
    f = quadratic_f;
    u = quadratic_u;
    grad_u = quadratic_grad_u;
    break;

  default:
    std::cerr << "[main] ERROR: Unknown exact solution" << std::endl;
    exit(1);    
  }
  // Symmetric gradient
  Elasticity::GradientDisplacementType symgrad_u = [&grad_u](const Eigen::Vector3d & x) -> Eigen::Matrix3d
                            { return symmetrise_matrix(grad_u(x)); };

  // Build the mesh and reorder based on the BCs
  MeshBuilder meshbuilder = MeshBuilder(mesh_file);
  std::unique_ptr<Mesh> mesh_ptr = meshbuilder.build_the_mesh();
  std::string bctype = vm["boundary-conditions"].as<std::string>();
  std::cout << FORMAT(25) << "[main] Boundary conditions: " << bctype << std::endl;
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
  Elasticity elas(ddr_core, mu, lambda, BC, use_threads);  
  if(vm.count("stabilization-parameter")) {
    elas.stabilizationParameter() = vm["stabilization-parameter"].as<double>();
  }
  Eigen::VectorXd UDir = elas.vectorDirBC(u);
  elas.assembleLinearSystem(f, symgrad_u, UDir);
  timer.stop();
  double t_wall_model, t_proc_model;
  std::tie(t_wall_model, t_proc_model) = store_times(timer, "[main] Time model (wall/proc) ");

  // Export matrix if requested  
  if (vm.count("export-matrix")) {
    std::cout << "[main] Exporting matrix to Matrix Market format" << std::endl;
    saveMarket(elas.systemMatrix(), "A_elasticity.mtx");
    saveMarket(elas.systemVector(), "b_elasticity.mtx");
  }

  // Select solver and solve statically condensed system
  timer.start();
  std::string name_solver = vm["solver"].as<std::string>();
  LinearSolver<Elasticity::SystemMatrixType> solver(name_solver);
  std::cout << "[main] Solving the system using " << solver.name() << std::endl;
  Eigen::VectorXd u_gl = solver.compute_and_solve(elas.systemMatrix(), elas.systemVector());
  std::cout << "[main] residual linear solver: " 
            << solver.residual(elas.systemMatrix(), elas.systemVector(), u_gl) << std::endl;
  
  // Re-create full vector of 
  Eigen::VectorXd uh = elas.reAssembleSolution(UDir, u_gl);
                               
  timer.stop();
  double t_wall_solve, t_proc_solve;
  std::tie(t_wall_solve, t_proc_solve) = store_times(timer, "[main] Time solve (wall/proc) ");
    
  //----------------------------------------------
  //    OUTPUTS
  //----------------------------------------------
  
  // Interpolate of exact solution and error vector
  Eigen::VectorXd uI = elas.vsxGrad().interpolate(u, MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
  Eigen::VectorXd eh = uh - uI;
   
  std::cout << "[main] Compute errors" << std::endl;
  // Errors in discrete norms
  std::cout << "[main] Discrete norms and errors:" << std::endl;
  std::vector<double> list_norms = elas.vsxGrad().computeH1norms(std::vector<Eigen::VectorXd> {eh, uI}, "symmetric");
  double error_symH1error = (list_norms[1] < 1e-12 ? list_norms[0] : list_norms[0]/list_norms[1]);
  std::cout << FORMAT(25) << "\tSymmetric H1 relative error: " << error_symH1error << " (norm of I_h u: "<< list_norms[1] << ")" << std::endl;
  std::cout << FORMAT(25) << "\tSymmetric H1 absolute error: " << list_norms[0] << std::endl;
   
  std::cout << "[main] Mesh diameter " << mesh_ptr->h_max() << std::endl;

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
      exact_u_vertex[iV] = u(v_coords);
    }
		
    // Approximate pressure and velocity at the vertices
    std::vector<VectorRd> u_vertex = elas.vsxGrad().computeVertexValues(uh);
        
    // Plot
    plotdata.write_to_vtu(filename + "_u.vtu", 
                          std::vector<std::vector<VectorRd>> {u_vertex, exact_u_vertex},
                          std::vector<std::string>{"displacement", "exact_displacement"});    
  }


  // --------------------------------------------------------------------------
  //                     Creates .txt file with data and results
  // --------------------------------------------------------------------------
  std::ofstream out("results.txt");
  out << std::scientific;
  out << "Solution: " << solution << std::endl;
  out << "Boundary conditions: " << bctype << std::endl;
  out << "mu: " << mu << std::endl;
  out << "lambda: " << lambda << std::endl;
  out << "StabilisationParameter: " << elas.stabilizationParameter() << std::endl;
  out << "Mesh: " << mesh_file << std::endl;
  out << "Degree: " << K << std::endl;
  out << "MeshSize: " << mesh_ptr->h_max() << std::endl;
  out << "NbCells: " << mesh_ptr->n_cells() << std::endl;
  out << "NbFaces: " << mesh_ptr->n_faces() << std::endl;
  out << "NbEdges: " << mesh_ptr->n_edges() << std::endl;
  out << "DimVSXGrad: " << elas.vsxGrad().dimension() << std::endl;
  out << "SizeSCsystem: " << elas.systemMatrix().rows() << std::endl;
  // Discrete errors and norms
  out << "E_symH1Dis: " << error_symH1error << std::endl;
  out << "absE_symH1Dis: " << list_norms[0] << std::endl;
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
// Elasticity
//------------------------------------------------------------------------------

Elasticity::Elasticity(
               const DDRCore & ddrcore,
               const double & mu,
               const double & lambda,
               const BoundaryConditions & BC,
               bool use_threads,
               std::ostream & output
               )
  : m_ddrcore(ddrcore),
    m_use_threads(use_threads),
    m_output(output),
    m_ser_pro(ddrcore, use_threads, output),
    m_vsxgrad(m_ser_pro),
    m_nloc_sc(m_ser_pro.nDOFs_cells_SXGrad() * dimspace),
    m_mu(mu),
    m_lambda(lambda),
    m_BC(BC),
    m_nDOFs_dir_vertices(0),
    m_nDOFs_dir_edges(0),
    m_nDOFs_dir_faces(0),
    m_locSCUKN({}),
    m_locGLUKN({}),
    m_DOFtoSCUKN(-2 * Eigen::VectorXi::Ones(nDOFs())),  // These maps are initialised at -2, but created just below
    m_DOFtoGLUKN(-2 * Eigen::VectorXi::Ones(nDOFs())),
    m_stab_par(0.1),
    m_A(nGLUKN(), nGLUKN()),
    m_b(Eigen::VectorXd::Zero(nGLUKN())),
    m_sc_A(nSCUKN(), nGLUKN()),
    m_sc_b(Eigen::VectorXd::Zero(nSCUKN()))
{
  m_output << "[Elasticity] Initializing" << std::endl;
  
  // To remove the static condensation, multiply by 0 the vectors in the initialisation of m_nloc_sc_u and m_nloc_sc_p.
  
  // Compute number of Dirichlet DOFs for displacement (vertices, edges, faces)
  for (Vertex * V : ddrcore.mesh().get_b_vertices()){
    if (BC.type(*V) == "dir"){
      m_nDOFs_dir_vertices += m_vsxgrad.numLocalDofsVertex(*V);
    }
  }
  for (Edge * E : ddrcore.mesh().get_b_edges()){
    if (BC.type(*E) == "dir"){
      m_nDOFs_dir_edges += m_vsxgrad.numLocalDofsEdge(*E);
    }
  }
  for (Face * F : ddrcore.mesh().get_b_faces()){
    if (BC.type(*F) == "dir"){
      m_nDOFs_dir_faces += m_vsxgrad.numLocalDofsFace(*F);
    }
  }

  // We compute the maps DOF->SCUKN and DOF->GLUKN by locating the positions of these unknowns through m_locSCUKN, m_locGLUKN,
  // and using replaceSectionsVector (see BChandlers). 
  // Following the construction of the m_locXXX vectors below is best done by doodling a picture of a global vector of DOFs
  // and locating the dirichlet vertex/edge/face DOFs in this global vector
  const Eigen::VectorXi MinusOnes = -Eigen::VectorXi::Ones(nDOFs());

  // ----- DOF->SCUKN -------
  m_locSCUKN.resize(1);
  // Statically condensed cell unknowns: start after all vertices, edges and faces
  m_locSCUKN[0].first = m_vsxgrad.nDOFs_vertices() + m_vsxgrad.nDOFs_edges() + m_vsxgrad.nDOFs_faces();
  m_locSCUKN[0].second = nSCUKN();
  // Create map
  size_t s = nSCUKN();
  m_DOFtoSCUKN = replaceSectionsVector(MinusOnes, Eigen::VectorXi::LinSpaced(s, 0, s-1).eval(), m_locSCUKN);

  // ----- DOF->GLUKN -------
  m_locGLUKN.resize(3);
  // non-Dirichlet vertex DOFs
  m_locGLUKN[0].first = m_nDOFs_dir_vertices;
  m_locGLUKN[0].second = m_vsxgrad.nDOFs_vertices() - m_nDOFs_dir_vertices; 
  // non-Dirichlet edge DOFs
  m_locGLUKN[1].first = m_vsxgrad.nDOFs_vertices() + m_nDOFs_dir_edges;
  m_locGLUKN[1].second = m_vsxgrad.nDOFs_edges() - m_nDOFs_dir_edges; 
  // non-Dirichlet face DOFs and non-SC cell unknowns
  m_locGLUKN[2].first = m_vsxgrad.nDOFs_vertices() + m_vsxgrad.nDOFs_edges() + m_nDOFs_dir_faces;
  m_locGLUKN[2].second = (m_vsxgrad.nDOFs_faces() - m_nDOFs_dir_faces) + (m_vsxgrad.nDOFs_cells() - nSCUKN());

  s = nGLUKN();
  m_DOFtoGLUKN = replaceSectionsVector(MinusOnes, Eigen::VectorXi::LinSpaced(s, 0, s-1).eval(), m_locGLUKN);

}

//------------------------------------------------------------------------------
Eigen::VectorXd Elasticity::vectorDirBC( const DisplacementType & u )
{
  // We interpolate u then put them in vectorDirBC (which has otherwise only zeros)
  // It's not the most efficient, but it is the easiest and this function is only called once anyway

  Eigen::VectorXd interp = m_vsxgrad.interpolate(u, MAX_DOE_CELL, MAX_DOE_FACE, MAX_DOE_EDGE);
  Eigen::VectorXd vectorBC = Eigen::VectorXd::Zero(nDOFs());
    
  for (Vertex * V : m_ddrcore.mesh().get_b_vertices()){
    if (m_BC.type(*V) == "dir"){
      vectorBC.segment(m_vsxgrad.globalOffset(*V), m_vsxgrad.numLocalDofsVertex(*V))
          = interp.segment(m_vsxgrad.globalOffset(*V), m_vsxgrad.numLocalDofsVertex(*V));
    }
  }

  for (Edge * E : m_ddrcore.mesh().get_b_edges()){
    if (m_BC.type(*E) == "dir"){
      vectorBC.segment(m_vsxgrad.globalOffset(*E), m_vsxgrad.numLocalDofsEdge(*E))
          = interp.segment(m_vsxgrad.globalOffset(*E), m_vsxgrad.numLocalDofsEdge(*E));
    }
  }

  for (Face * F : m_ddrcore.mesh().get_b_faces()){
    if (m_BC.type(*F) == "dir"){
      vectorBC.segment(m_vsxgrad.globalOffset(*F), m_vsxgrad.numLocalDofsFace(*F))
          = interp.segment(m_vsxgrad.globalOffset(*F), m_vsxgrad.numLocalDofsFace(*F));
    }
  }
    
  return vectorBC;
}

//------------------------------------------------------------------------------

void Elasticity::assembleLinearSystem(
              const ForcingTermType & f, 
              const GradientDisplacementType & symgrad_u, 
              const Eigen::VectorXd & UDir )
{
  // Assemble all local contributions
  auto assemble_all = [this, &f, &symgrad_u, &UDir](
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
                                                this->_compute_local_contribution(iT, f, symgrad_u, UDir),
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
                               const GradientDisplacementType & symgrad_u,
                               const Eigen::VectorXd & UDir
                               ) const
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  size_t dim_T = m_vsxgrad.dimensionCell(iT);
  
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
  AT += 2 * m_mu * GsT.transpose() * mass_Pk3x3_T * GsT + m_mu * m_vsxgrad.computeStabilisation(iT);

  // Contribution lambda div(u) div(v)
  Eigen::MatrixXd divT = m_vsxgrad.cellDivergence(iT);
  AT += m_lambda * divT.transpose() * mass_Pk_T * divT;
    
  //------------------------------------------------------------------------------
  // Local source vector
  //------------------------------------------------------------------------------
  // Forcing term
  QuadratureRule quad_2kpt_T = generate_quadrature_rule(T, 2 * m_vsxgrad.degree() + 2);
  std::function<Eigen::Vector3d(const Eigen::Vector3d&)> f_mu_lambda = [&f, this](const Eigen::Vector3d & x)->Eigen::Vector3d
                  { return f(this->m_mu,this->m_lambda,x);};
  lT += m_vsxgrad.cellPotential(iT).transpose()
    * integrate(f_mu_lambda,
                evaluate_quad<Function>::compute(*(this->m_vsxgrad.cellBases(iT).Polykpo3), quad_2kpt_T),
                quad_2kpt_T);
  
  // Dirichlet BCs
  lT -= AT * m_vsxgrad.restrict(T, UDir);
  
  // Neumann BCs
  if (T.is_boundary()){
    const size_t dqrbc = MAX_DOE_FACE;
    for (Face * pF : T.get_faces()){
      const Face & F = *pF;
      if (m_BC.type(F) == "neu"){
        // Unit normal vector to F pointing out of the domain
        const Cell & TF = *F.cell(0);
        Eigen::Vector3d nF = TF.face_normal(TF.index_face(&F));
        
        QuadratureRule quad_dqrbc_F = generate_quadrature_rule(F, dqrbc);
        
        // Traction boundary conditions: int_F C eps(u)n . P_F v
        std::function<Eigen::Vector3d(const Eigen::Vector3d &)> 
            Cepsu_nF = [&symgrad_u, &nF, this](const Eigen::Vector3d & x)->Eigen::Vector3d {
                         return 2 * this->m_mu * symgrad_u(x) * nF + this->m_lambda * symgrad_u(x).trace() * nF;
                         };
        // Create integral of Ceps(u)nF * PotF as a linear (row) operator on the face DOFs, to extended it afterwards to the cell DOFs.
        Eigen::RowVectorXd bF =  
              integrate(Cepsu_nF, 
                        evaluate_quad<Function>::compute(*m_vsxgrad.faceBases(F).Polykpo3, quad_dqrbc_F),
                        quad_dqrbc_F).transpose() 
                      * m_vsxgrad.facePotential(F);
                                 
        // Add contribution
        lT += m_vsxgrad.extendOperator(T, F, bF).transpose();
        
      } // if m_BC.type(F)="neu"
    } // Loop over F
  } // if T.is_boundary()

  return std::make_pair(AT, lT);
}

//------------------------------------------------------------------------------

void Elasticity::_assemble_local_contribution(
                                          size_t iT,
                                          const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT,
                                          std::list<Eigen::Triplet<double> > & triplets_gl,
                                          Eigen::VectorXd & rhs_gl,
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
  size_t dim_dofs = m_vsxgrad.dimensionCell(iT) - dim_sc;     // number of local displacement unknowns after SC

  // Creation of permutation matrix: trivial here since SC dofs are already at the end (cell unknowns)
  Eigen::MatrixXd Perm = Eigen::MatrixXd::Identity(dim_dofs+dim_sc, dim_dofs+dim_sc);

  // Creation of global DOFs for system
  //   DOFs_gl and DOFs_sc respectively contains the DOFs of u for globally coupled system and for the recovery system
  //   These DOFs need to be read through DOFtoGLUKN and DOFtoSCUKN when assembling the systems
  std::vector<size_t> DOFs_gl(dim_dofs, 0);
  std::vector<size_t> DOFs_sc(dim_sc, 0);
  auto DOFs_T_vsxgrad = m_vsxgrad.globalDOFIndices(T);
   
  std::copy(DOFs_T_vsxgrad.begin(), DOFs_T_vsxgrad.begin()+dim_dofs, DOFs_gl.begin()); // put skeletal DOFs in DOFs_gl
   
  if (dim_sc>0){
    std::copy(DOFs_T_vsxgrad.begin()+dim_dofs, DOFs_T_vsxgrad.end(), DOFs_sc.begin());
  }  

  return LocalStaticCondensation(Perm, DOFs_gl, DOFs_sc);
}

