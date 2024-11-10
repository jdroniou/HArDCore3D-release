HArDCore3D: Principle of implementation of a scheme
===================================================

This tutorial illustrate how the tools provided by HArDCore can be used to implement a scheme. The examples shown here come from the [HHO_fullgradientdiff](@ref HHO_fullgradientdiff) implementation of the Hybrid High-Order scheme for a diffusion equation.

* [Overview](\ref sec-overview) -- General overview and principles.
* [Spaces of degrees of freedom](\ref sec-dofs) -- Description of the classes of spaces of DOFs for the scheme.
* [Discretisation method](\ref sec-method) -- Principle for classes describing the discretisation method.
* [Assembly](\ref sec-assembly) -- How to organise the assembly of the systems describing the scheme.
* [Output](\ref sec-outputs) -- A quick word on the outputs.


<!--<a name="sec_overview">-->
\section sec-overview Overview
<!--</a>-->

The implementation of a scheme starts from a generic class (in src/) that describes the discretisation method, but is oblivious to the particular model to discretise.
For example, [HHOSpace](@ref HHOSpace) encodes the degrees of freedom and operators required to implement an HHO scheme: unknowns on the faces and cells (described as coefficients on specific polynomial bases built in the class), local potential reconstruction, gradient and stabilisation, etc.

The scheme itself is implemented in separate files (in Schemes/). These files describe in particular a class, e.g. HArDCore3D::FullGradientDiffusion, that is specific to the model under consideration and the discretisation method. That class has variables to store in particular the system to assemble, and methods to assemble this system.
The .hpp file of the scheme also contains particular test cases (choices of unknowns, to fix the source terms and boundary conditions).

A quick look at the .cpp file shows that it starts by defining a series of options that can be passed to the executable, to change some parameters of the simulation without having to re-compile the scheme. For example:

\code{.cpp}
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Display this help message")
    ("mesh,m", boost::program_options::value<std::string>(), "Set the mesh")
    ("degree,k", boost::program_options::value<size_t>()->default_value(1), "The polynomial degree of the sequence")
    ("pthread,p", boost::program_options::value<bool>()->default_value(true), "Use thread-based parallelism")
    ("solution,s", boost::program_options::value<int>()->default_value(0), "Select the solution")
    ("export-matrix,e", "Export matrix to Matrix Market format")
    ("plot", boost::program_options::value<std::string>(), "Save plot of the solution to the given filename")
    ("solver", boost::program_options::value<std::string>()->default_value("PardisoLU"), "Choice of solver, not case dependent. Options are: PardisoLU, UMFPACK, PaStiXLU, PaStiXLLT, EigenLU, EigenBiCGSTAB (reverts to EigenLU if the selected solver is not available)")    
    ("stabilization-parameter,x", boost::program_options::value<double>(), "Set the stabilization parameter");

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);
\endcode

These options are used to fix some parameters, after which the mesh is build, the boundary conditions are selected (when relevant in the scheme) -- with numbering of vertices/edges/faces can be modified, e.g. to facilitate their treatment in the implementation -- and the classes describing the method and the scheme are instantiated. Some parameters of these classes are possibly further adjusted depending on the options passed to the executable:

\code{.cpp}
  // Build the mesh and reorder faces to handle Dirichlet BCs (Dirichlet faces at the end)
  MeshBuilder meshbuilder = MeshBuilder(mesh_file);
  std::unique_ptr<Mesh> mesh_ptr = meshbuilder.build_the_mesh();
  BoundaryConditions BC("D", *mesh_ptr.get());
  BC.reorder_faces("start"); 
  
  // Create HHO space
  bool use_threads = (vm.count("pthread") ? vm["pthread"].as<bool>() : true);
  std::cout << "[main] " << (use_threads ? "Parallel execution" : "Sequential execution") << std:: endl;
  HHOSpace hho_space(*mesh_ptr, K, use_threads);

  (...)
  FullGradientDiffusion diff(hho_space, BC, use_threads);
  if(vm.count("stabilization-parameter")) {
    diff.stabilizationParameter() = vm["stabilization-parameter"].as<double>();
  }
\endcode

From that point on, the main method to assemble the system are called, the system is solved, and the complete solution is re-assembled (taking into account
the boundary conditions, and possibly the static condensation applied when the system was assembled). The solver can be chosen using the HArDCore3D::LinearSolver class, which contains various solvers (depending on the libraries found/made by cmake, see INSTALL_NOTES.txt).

\code{.cpp}
  // From hho-fullgradientdiff.cpp

  diff.assembleLinearSystem(f, kappa, u, UDir);

  (...)

  // Select linear solver (based on the parameter "solver" on the command line)
  std::string name_solver = vm["solver"].as<std::string>();
  LinearSolver<Stokes::SystemMatrixType> solver(name_solver);
  std::cout << "[main] Solving the system using " << solver.name() << std::endl;

  // Solve the problem
  Eigen::VectorXd uh_solsystem = solver.compute_and_solve(diff.systemMatrix(), diff.systemVector());

  // Re-create boundary values and statically condensed unknowns
  Eigen::VectorXd uh = Eigen::VectorXd::Zero(diff.hhospace().dimension());
  uh.head(diff.numDirDOFs()) = UDir;
  uh.segment(diff.numDirDOFs(), diff.sizeSystem()) = uh_solsystem;
  uh.tail(diff.numSCDOFs()) = diff.scVector() + diff.scMatrix() * uh.head(diff.numSkeletalDOFs());
\endcode

Finally, the output (numerical values or figures) is created.


<!--<a name="sec_dof_spaces">-->
\section sec-dofs Space of degrees of freedom
<!--</a>-->

The core class describes the degrees of freedom of the method. There are actually two such classes: HArDCore3D::GlobalDOFSpace for methods that have the same number of degrees of freedom on each mesh entity of the same dimension (all vertices have x numbers of DOFs, all edges have y numbers of DOFs, etc.), and HArDCore3D::VariableDOFSpace for methods with varying number of degree of freedom from mesh entity to mesh entity (vertex 1 has x_1 DOFs, vertex 2 has x_2 DOFs, etc.).

Each of these classes assume that the degrees of freedom are organised from the lowest-dimensional mesh entity to the highest-dimensional ones, and following the (local or global) numbering of these mesh entites. So, for example, if an element \f$T\f$ has vertices \f$V_1\f$, \f$V_2\f$..., \f$V_n\f$, edges \f$E_1\f$, \f$E_2\f$..., \f$E_r\f$, and faces \f$F_1\f$, \f$F_2\f$..., \f$F_l\f$, the degrees of freedom are organised in this order: 
  - DOFs attached to \f$V_1\f$,
  - DOFs attached to  \f$V_2\f$,
  - ...
  - DOFs attached to \f$V_n\f$,
  - DOFs attached to \f$E_1\f$,
  - DOFs attached to \f$E_2\f$,
  - ...
  - DOFs attached to \f$E_r\f$,
  - DOFs attached to \f$F_1\f$,
  - DOFs attached to  \f$F_2\f$,
  - ...
  - DOFs attached to \f$F_l\f$,
  - DOFs attached to \f$T\f$.
  
These classes provide in particular methods to:
  - easily access the location (either globally in the mesh, or locally in an element, a face, etc.) of the DOFs attached to particular mesh entities,
  - restrict a global vector of DOFs to a particular mesh entity, 
  - extend a matrix acting, say, on the DOFs of an edge into a matrix acting on the DOFs of a mesh entity and its boundary (filling with 0 for DOFs not involved in the design of this operator),
  - list the global numbering of DOFs attached to a given mesh face or element.


<!--<a name="sec_method">-->
\section sec-method Discretisation method
<!--</a>-->

A discretisation method is defined by a choice of DOFs and local operators acting on these DOFs. A class describing a discretisation method is therefore inherited from one of the two classes describing spaces of DOFs, to embed the information relevant to the choice of these degrees of freedom. The class then builds and stores (or provide functions to build on the fly) the relevant local operators.

HArDCore3D::HHOSpace is built on HArDCore3D::GlobalDOFSpace. Its instantiantion therefore fixes how many DOFs each vertex, edge, face and element should have. For a chosen polynomial degree K, the HHO method is based on unknowns that are polynomials of degree K on each face and polynomials of degree K in each element. Passing the argument K to the class as in 

\code{.cpp}
  HHOSpace hho_space(*mesh_ptr, K, use_threads);
\endcode

therefore results in a HArDCore3D::GlobalDOFSpace designed with zero DOFs on each vertex and edge, \f$dim\mathbb{P}^K(\mathbb{R}^2)\f$ DOFs on each face and \f$dim\mathbb{P}^K(\mathbb{R}^3)\f$ DOFs on each element:

\code{.cpp}
  // From hhospace.cpp

  HHOSpace::HHOSpace(const Mesh & mesh, size_t K, bool use_threads, std::ostream & output)
    : GlobalDOFSpace(mesh,
	       0,
	       0,
	       PolynomialSpaceDimension<Face>::Poly(K),
	       PolynomialSpaceDimension<Cell>::Poly(K)
	       ),
	    (...)
\endcode

Additionally, to give meaning to these DOFs on each mesh entity, the class HArDCore3D::HHOSpace builds and stores bases of the polynomial spaces of degree K on each face and each element (the construction is done here using the multi-threading capacities provided by parallel_for):

\code{.cpp}
  // From hhospace.cpp

  // Construct element bases
  std::function<void(size_t, size_t)> construct_all_cell_bases
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iT = start; iT < end; iT++) {
	        this->m_cell_bases[iT].reset( new CellBases(this->_construct_cell_bases(iT)) );
	      } // for iT
      };

  m_output << "[HHOSpace] Constructing element bases" << std::endl;
  parallel_for(mesh.n_cells(), construct_all_cell_bases, use_threads);
\endcode

Once a vector of DOFs is known, its corresponding polynomial values on a face/element can therefore be fully described through the coefficients in the vector (the restriction of these to a face/element being accessible via the methods of HArDCore3D::GlobalDOFSpace) and the polynomial basis functions embedded in the instantiation of HArDCore3D::HHOSpace.

Finally, when this class is instantiated, the local operators (potential, gradient, stabilisation) are constructed and stored:

\code{.cpp}
  // From hhospace.cpp

  // Construct gradients, potentials and stabilisation
  std::function<void(size_t, size_t)> construct_all_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_operators[iT].reset( new LocalOperators(_compute_operators(iT)) );
        } // for iT
      };

  m_output << "[HHOSpace] Constructing operators" << std::endl;
  parallel_for(mesh.n_cells(), construct_all_operators, use_threads);
\endcode

The construction of these operators is method-dependent, but most often require to set up systems of equations based on matrices obtained by integrating products of 
polynomial basis functions (these matrices are Gram-matrix like). Such integrals (on edges, faces, elements) can be done using two approaches:
  - Standard quadrature rules: weights and nodes are constructed using generate_quadrature_rule(), the basis functions are evaluated at these nodes, and the Gram matrix is built via compute_gram_matrix():
  \code{.cpp}
    // Adapted from hhospace.cpp, computes the gram matrix of the basis function of P^k(F).
    
    QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree() );
    auto basis_Pk_F_quad = evaluate_quad<Function>::compute(*faceBases(F).Polyk, quad_2k_F);
    G = compute_gram_matrix(basis_Pk_F_quad, quad_2k_F);
  \endcode
  
  - Polynomial integration rules: computationally much faster way to compute integrals of polynomial bases, and associated Gram matrix, for certain classes of polynomials (see the mainpage for the list). This uses the GramMatrix() functions.
  \code{.cpp}
    // From hhospace.cpp
    
    MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2);
    Eigen::MatrixXd MGT = GramMatrix(T, *cellBases(iT).Polykd, int_mono_2kp2);
  \endcode
  (calling IntegrateCellMonomials is optional, but avoids re-computing the integrals of the monomials on each cell if they are to be re-used in another GramMatrix() in the same scope).

<!--<a name="sec_assembly">-->
\section sec-assembly Local and global assembly
<!--</a>-->

Most numerical schemes are assembled by computing local element-wise contributions, and distributing them in the global system. The local contributions are based on the operators built in the class describing the discretisation, and are therefore very easy to compute from this class. For HArDCore3D::FullGradientDiffusion, this is done by HArDCore3D::FullGradientDiffusion::_compute_local_contribution(), e.g.:

\code{.cpp}
  // From hho-fullgradientdiff.cpp

  //------------------------------------------------------------------------------
  // Local matrix
  //------------------------------------------------------------------------------

  // Mass matrix for (P^k(T))^d
  MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_hhospace.degree());  
  Eigen::MatrixXd mass_Pkd_T = GramMatrix(T, *m_hhospace.cellBases(iT).Polykd, int_mono_2k);

  double kappaT = kappa.value(T, T.center_mass());
  AT += kappaT * m_hhospace.operators(iT).gradient.transpose() * mass_Pkd_T * m_hhospace.operators(iT).gradient
       + kappaT * m_stab_par * m_hhospace.operators(iT).stabilisation;
  
  //------------------------------------------------------------------------------
  // Local vector
  //------------------------------------------------------------------------------

  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * m_hhospace.degree());
  lT.tail(m_hhospace.numLocalDofsCell()) = 
      integrate(f, evaluate_quad<Function>::compute(*m_hhospace.cellBases(iT).Polyk, quad_2k_T), quad_2k_T);
  
\endcode

In theory, the distribution of the local contribution is then straightforward using the local-to-global map available in the class describing the method (which inherits this from the class describing the space of DOFs), see e.g. HArDCore3D::HHOSpace::globalDOFIndices(). This global assembly would look something like this:

\code{.cpp}
  // Global assembly
  for (size_t iT=0; iT<mesh.n_cells(); iT++){
    const Cell & T = *mesh.cell(iT);

    // get local contribution
    Eigen::MatrixXd AT;
    Eigen::VectorXd lT;
    std::tie(AT, lT) = _compute_local_contribution(iT);

    // Assemble in global triplets and rhs
    std::vector<size_t> IT = m_hhospace.globalDOFIndices(T);
    for (size_t i = 0; i < IT.size(); i++){
      rhs(IT[i]) += lT(i);
      for (size_t j = 0; j < IT.size(); j++){    
        triplets.emplace_back(IT[i], IT[j], AT(i,j));
      }
    }
  }
  // Create global matrix, N is the size of the global system
  Eigen::SparseMatrix<double> GlobalMatrix(N, N);
  GlobalMatrix.setFromTriplets(triplets.begin(), triplets.end());
  
\endcode

However, to make the assembly more efficient, using the multi-threaded procedure parallel_assembly_system() is recommended. This leads to a slightly more complicated local to global assembly step, based on the definition of a functor that will invoke the corresponding _compute_local_contributions() and is called in multi-threads by parallel_assembly_system(). The function assembleLinearSystem() embeds the definition of the functor and the call to parallel_assembly_system(). The situation can also be complicated when using static condensation, to reduce the size of the globally coupled system; in that case, two systems are assembled at the same time par parallel_assembly_system(): the system on the globally coupled unknowns, and the system used to recover the statically condensed unknowns.

In practice however, applying the multi-threading and the static condensation is relatively easy. For example, the skeleton of HArDCore3D::FullGradientDiffusion::assembleLinearSystem() barely has to be modifed when changing the method or the model. The process goes like this, from the innermost to the outermost procedure:
  - _compute_local_contribution() computes the local contributions, element by element,
  - These contributions are passed to _assemble_local_contribution() which:
    * calls _compute_static_condensation() which prepares the input data to instantiate a HArDCore3D::LocalStaticCondensation class (this input data identifies in particular the DOFs that are statically condensed via a permutation matrix that puts these DOFs at the end, and also contains the global DOF indices for statically condensed and non-statically condensed unknowns).
    * this HArDCore3D::LocalStaticCondensation then performs the local static condensation from the local contributions,
    * finally, these contributions are distributed in the system corresponding to globally coupled unknowns, and the system used to recover the statically condensed unknowns.




<!--<a name="sec_outputs">-->
\section sec-outputs Outputs
<!--</a>-->

Classes corresponding to discretisation methods or schemes often provide specific outputs of interest (e.g. a certain norm, etc.). Additionally, some classes provide
methods to compute values of the solution at the vertices of the mesh. In this case, a vtk plot can be generated using HArDCore3D::VtuWriter.

\code{.cpp}
  // From hho-fullgradientdiff.cpp
  
  // Only if we do not have too many cells
  if (vm.count("plot") && mesh_ptr->n_cells() <= max_nb_cells_for_plot) {
    std::cout << "[main] Writing solution to file" << std::endl;
    std::string filename = vm["plot"].as<std::string>();
    VtuWriter plotdata(mesh_ptr.get());
		
		// Exact solution at the vertices
    std::vector<double> exact_u_vertex(mesh_ptr->n_vertices(), 0.);
    for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
      exact_u_vertex[iV] = u(mesh_ptr->vertex(iV)->coords());
    }
		
    // Approximate pressure and velocity at the vertices
    std::vector<double> u_vertex = hho_space.computeVertexValues(uh);
    
    // Plot
    plotdata.write_to_vtu(filename + ".vtu", 
                          std::vector<std::vector<double>> {u_vertex, exact_u_vertex},
                          std::vector<std::string>{"approximate", "exact"});
   }
\endcode





