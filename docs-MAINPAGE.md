
HArD::Core (sources: https://github.com/jdroniou/HArDCore) provides a suite of C++ tools to implement numerical schemes whose unknowns are polynomials in the cells, on the edges, and on the faces. The focus is on dealing on generic polytopal meshes. This documentation addresses the 3D version of HArD::Core, but similar principles are valid for the 2D version. Transferring a scheme's implementation from 3D to 2D or vice-versa is very straightforward, provided that the scheme's mathematical definition does not depend on the dimension and that the generic types provided in `basis.hpp` are used; see readme file of the HArD::Core github depository https://github.com/jdroniou/HArDCore.

\tableofcontents


<a name="build">
\section build Build instructions
</a>

\subsection buildlib Building the libraries and the schemes

To build the libraries and implemented schemes, the minimal requirements are:

* CMake version 2.6 or above (https://cmake.org/)
* A C++ compiler that supports the C++14 standard, eg. GCC (https://gcc.gnu.org/) or Clang (https://clang.llvm.org/)
* Eigen C++ library, version 3.3 or above (http://eigen.tuxfamily.org/)
* The following Boost C++ libraries (http://www.boost.org/): filesystem, program options, timer, chrono.

Make sure that you have the development version of boost installed. On Linux, install `libboost-dev`, `libboost-filesystem-dev`, `libboost-program-options-dev`, `libboost-chrono-dev` and `libboost-timer-dev` from your package manager.

The linear systems resulting from the assembled scheme are solved using the BiCGStab implementation of Eigen. Alternatives are also provided, but require additional libraries (UMFPACK, SUPERLU, etc.); see the main CMakeLists.txt file.

Once you have installed all of the required dependencies, set up the build directory and generate the build files by running the following from the repository root:

```
mkdir build
cd build
cmake ..
make
```

After this, `build/Schemes` will contain the executables (e.g. `hho-diffusion`) to run the schemes. These executables need to access the typ2 meshes, which they should naturally find if you put the `typ2_meshes` directory at the root of the project's files.



\subsection doco Building the Documentation

The mesh documentation is built with Doxygen (http://www.stack.nl/~dimitri/doxygen/). If you are reading this then somebody has already built it for you. If you modify the code and wish to rebuild the documentation, simply run `doxygen` from the root directory. The HTML version of the documentation is generated inside `documentation/html` and the LaTeX version is generated inside `documentation/latex` and can be compiled using the generated Makefile.


<a name="mesh">
\section mesh Mesh module
</a>

\subsection meshpple Principles

After it is loaded, the mesh is represented by classes describing a vertex, an edge, a face, and a cell: [Vertex](@ref HArDCore3D::Vertex), [Edge](@ref HArDCore3D::Edge), [Face](@ref HArDCore3D::Face), and [Cell](@ref HArDCore3D::Cell). Each of these classes contains methods to access useful information for the corresponding element, including other geometrical quantities it is related to. The mesh itself is represented by an element of the [Mesh](@ref HArDCore3D::Mesh) class with methods to access all the vertices, edges, faces and cells (or a particular vertex, edge, face or cell). 

For example, if `mesh_ptr` is a pointer to a Mesh class, the lines
\code{.cpp}
Vertex* vertex = mesh_ptr->vertex(5);

Eigen::Vector3d vert_coord = vertex->coords()
\endcode
store the coordinates of the fifth vertex into the Eigen vector vert_coord. As a generic rule, all geometrical vectors are `Eigen::Vector3d`. We also use `Eigen::Vector{3,X}d` and `Eigen::Matrix{3,X}d` for objects on which linear algebraic operations are performed. Lists (e.g. of cells, of functions...) are usually instances of `std::vector<...>`. Finally, `Eigen::multi\_array` is used for lists of values of basis functions or their gradients at quadrature nodes.

Here is an example that loops over all cells, grabs all the faces of the cell, and loops over these faces to output their diameter. Here, `mesh_ptr` is a pointer to the mesh.

\code{.cpp}
// Loop over all cells of the mesh
for (size_t iT = 0; iT < mesh_ptr->n_cells() iT++) {

	// We grab the faces of the iT-th cell
	std::vector<Face *> faces = mesh_ptr->cell(iT)->get_faces();

	// Loop over the faces of the cell
	for (size_t ilF = 0; ilF < cell->n_faces(); ilF++) {

		// Write the face diameter on the standard output
		std::cout << "The diameter of face " << ilF+1 << " in cell " << iT+1 << " is: " << faces(ilE)->diam() << "\n";
	}

}
\endcode

The mesh classes and other auxilliary classes are located inside the namespace [HArDCore3D](@ref HArDCore3D).

There is no direct access from a high-level geometrical entity to elements purely associated with lower-level entities. For example, if `mesh_ptr` is a pointer to the mesh, there is no direct method to access the coordinates of the i-th vertex of the mesh (no `mesh_ptr->coords_vertex()` exists). Instead, this is done through `mesh_ptr->vertex(i)->coords()`. This choice is deliberate as it preserves the logical organisation of the data structure, and facilitates the memorisation of the various methods. Of course, writing a wrapper providing such a direct access is easy...


\subsection loading_mesh Loading a mesh

HArDCore3D can currently read meshes in `TG`, `RF` and `MSH`. The loader and all geometrical computations are performed using `StemMesh`, based on G. Manzini's mesh library https://github.com/gmanzini-LANL/PDE-Mesh-Manager.

A mesh file must be read using an instance of the `meshbuilder` class, and then built using `build_the_mesh`.  A working example is given below (assuming the executable will be in `build/Schemes` for example).

\code{.cpp}
#include "mesh.hpp"
#include "mesh_builder.hpp"

using namespace HArDCore3D;

int main() {

	// Mesh file to read, with type
	std::string default_mesh = "../../meshes/Voro-small-0/RF_fmt/voro-4";
	std::string default_meshtype = "RF";

	// Read the mesh file and build the mesh
	MeshBuilder meshbuilder = MeshBuilder(mesh_file, mesh_type);
	std::unique_ptr<Mesh> mesh_ptr = meshbuilder.build_the_mesh();

	std::cout << "There are " << mesh_ptr->n_cells() << " cells in the mesh.\n";
}
\endcode

Note that the builder returns a `unique_ptr` for the mesh. This ensures that, at the end of the code, the mesh destructor is called (which destroys all cells, faces, edges, vertices...). Some classes and functions use a raw pointer to the mesh, so the `.get()` method should be used when passing the mesh as argument to the class constructors or functions.

<i>Note</i>: the mesh formats allow for meshes with very generic polygonal cells, including non-convex cells.
However, the builder assumes that each cell is star-shaped with respect to the isobarycenter of its vertices -- otherwise, the calculation of the center of mass may be incorrect. Similarly, the quadrature rules (see [Quadrature rules](#quad_rules)) assume that each cell is star-shaped with respect to its center of mass.

<a name="quad_rules">
\section quad_rules Quadrature rules
</a>

HArD::Core deals with quite arbitrary cell geometries. As a consequence, no reference element can be used, and the quadrature rules have to be adjusted to each particular cell/face/edge. For the cells, for example, this is done by partitioning each cell into tetrahedra and by using classical quadrature rules on tetrahedras. The choice was also made not to pre-compute all quadrature rules for all cells, faces and edges, but rather to compute them -- with a locally chosen degree of exactness -- when needed in the code. To reduce the computational cost, quadrature rules -- and the values of basis functions at quadrature nodes -- should only be computed once when looping over each cell, before being used, e.g., to form mass matrices.

The [Quadratures](@ref HArDCore3D::Quadratures) module provides routines to do that. The key method is [generate_quadrature_rule(CFE,doe)](@ref HArDCore3D::generate_quadrature_rule) calculates quadrature nodes and weights, exact up to the polynomial degree `doe`, for an integration over cell/face/edge CFE (passed as a reference of [Cell](@ref HArDCore3D::Cell), [Face](@ref HArDCore3D::Face) or [Edge](@ref HArDCore3D::Edge) class). At present, the quadrature rules available in the code support a total degree of exactness in the cells up to 14 in the cells and 20 on the faces and the edges (the quadrature rules on the faces come from [John Burkardt's implementation of the Dunavant rules](https://people.sc.fsu.edu/~jburkardt/cpp_src/triangle_dunavant_rule/triangle_dunavant_rule.html)). The generated quadrature rule is stored in a structure [QuadratureRule](@ref HArDCore3D::QuadratureRule), which is a vector of quadrature nodes (weights and position).

<a name="basis">
\section basis_module Basis module
</a>

The [Basis](@ref Basis) module provides classes and functions that define and manipulate polynomial basis functions on a cell, face, or edge, and various derived quantities (gradient, curl...). The underlying basis functions are monomial, but families of functions defined as linear combination of monomial basis functions are also handled via the `Family` class. See examples in [HybridCore](#hybridcore).

Free functions are also available to compute basis functions at quadrature nodes (using the [Quadrature rules](#quad_rules) module), orthonormalise basis functions, and compute Gram-like matrices between various families of functions. These matrices are essential in the design of high-order methods on polytopal meshes. See examples in [HybridCore](#hybridcore).

<a name="hybridcore">
\section hybridcore Hybridcore module
</a>

The [HybridCore](@ref HArDCore3D::HybridCore) module encapsulates routines to create bases of polynomial spaces in each cell, on each face, and on each edge, and to manipulate discrete functions through the class [UVector](@ref HArDCore3D::UVector). 

\subsection basisfunc Basis functions

The instantiation of an [HybridCore](@ref HArDCore3D::HybridCore) class creates basis functions for the polynomial spaces in the cells, on the faces and on the edges, specifying the maximum degree required for each geometric entity. The basis functions are elements of the [Family](@ref HArDCore3D::Family) class and are accessed via the [CellBasis](@ref HArDCore3D::CellBasis), [FaceBasis](@ref HArDCore3D::FaceBasis) and [EdgeBasis](@ref HArDCore3D::EdgeBasis) method. For example, the following piece of code initialises an HybridCore instance with degrees \f$K+1\f$ in the cells, \f$K\f$ on the faces, and no edge basis functions, and access the value at some Eigen::Vector3d x of the i-th basis function on face iF, and the gradient of the j-th basis function in cell iT.

\code{.cpp}
  // Initialise the class
  HybridCore hho(mesh_ptr.get(), K+1, K, -1, use_threads, output);
  // Access value of face basis function
  double val_face = hho.FaceBasis(iF).function(i, x);
  // Access value of gradient of cell basis function
  Eigen::Vector3d grad_cell = hho.CellBasis(iT).gradient(j, x);
\endcode

The basis functions are hierarchical, which means that they are constructed by increasing degree. Hence, for example, if cell basis functions up to degree \f$K+1\f$ have been generated, a basis of the space of polynomials of degree \f$K\f$ in the cell is thus obtained by selecting the first \f$(K+1)(K+2)/2\f$ cell basis functions.

The [UVector](@ref HArDCore3D::UVector) class describes coefficients on cell and face basis functions. The first coefficients correspond to cell basis functions, ordered by the cells themselves, and the last coefficients correspond to face basis functions. The methods in this class provide the surrounding structure to make sense of these coefficients (degrees of considered polynomial functions in cells/on faces, restrictions to a certain cell and its faces, etc.).

\subsection usge_quad Usage of quadrature rules.

The [evaluate_quad](@ref HArDCore3D::evaluate_quad) template function evaluate basis functions, or gradients, etc. at provided quadrature nodes. The `boost::multi\_array` provided by this function can then be passed to [compute_gram_matrix](@ref HArDCore3D::compute_gram_matrix) to create a matrix of inner products of two families of basis functions. Here is an example.

\code{.cpp}
// Create basis (f_1,...,f_r) of degree k in cell T
MonomialScalarBasisCell basisT(T, k);
// Create quadrature rules of degree 2*k in cell T
QuadratureRule quadT = generate_quadrature_rule(T, 2*k);
// Compute values of gradients of basis functions at the quadrature nodes
boost::multi_array<VectorRd, 2> gradbasis_on_quadT = evaluate_quad<Gradient>::compute(basisT, quadT);
// Create Gram-like matrix (here, a stiffness matix) of (\nabla f_i,\nabla f_j)
Eigen::MatrixXd M = compute_gram_matrix(gradbasis_on_quadT, gradbasis_on_quadT, quadT, true);
\endcode

Note the usage of the type `VectorRd` defined in `basis.hpp`, which enables for a dimension-independent piece of code (easier to adapt to the 2D case). This procedure can also be applied, e.g., to cell basis functions on face quadrature nodes, etc.

However, in the [HybridCore](@ref HArDCore3D::HybridCore) class, the quadrature rules and values of basis functions (and gradients) at the quadrature nodes can be conveniently computed and stored using the [ElementQuad](@ref HArDCore3D::ElementQuad) class. Instantiating an element of this class on a cell loads these rules and values once, that can then be passed to several functions in charge of various calculations (e.g. one function computes the local cell contribution to the diffusion term, another function is in charge of computing the load term associated to the cell, etc.). This prevents recomputing these rules and values when needed by various functions. It works the following way:

\code{.cpp}
HybridCore hho(mesh_ptr.get(), K+1, K, -1, use_threads, output);    // HybridCore instantiation
size_t doeT = m_Ldeg + m_K + 1;     // degree of exactness for cell quadrature rules
size_t doeF = 2*m_K + 1;            // degree of exactness for edge quadrature rules
ElementQuad elquad(hho, iT, doeT, doeF);  // compute local quadrature rules at quadrature points in cell iT

Eigen::MatrixXd aT = diffusion_operator(hho, iT, elquad);		// compute local contribution to diffusion term
Eigen::VectorXd bT = load_operator(hho, iT, elquad);		//	compute local loading term

(...)
// Function to compute local contribution to diffusion term
Eigen::MatrixXd HHO_Diffusion::diffusion_operator(HybridCore &hho, const size_t iT, const ElementQuad &elquad) const {

(... initialise/do stuff ...)
// Cell quadrature rules and values at nodes are needed, we grab them
QuadratureRule quadT = elquad.get_quadT();
boost::multi_array<double, 2> phiT_quadT = elquad.get_phiT_quadT();
boost::multi_array<VectorRd, 2> dphiT_quadT = elquad.get_dphiT_quadT();

(... the rest as in the previous example: create mass matrices, etc. ...)

\endcode

<a name="hho3D">
\section hho_3D HHO3D general module
</a>

The [HHO3D](@ref HHO3D) module provides typedefs, a class and functions to implement Hybrid High Order (HHO) schemes. Rules (functions) to create local bilinear forms (matrices) and loading terms (vectors) are passed to the HHO3D class, that takes care of the global assembly and solving of the system.

<a name="ddr">
\section ddr DDR core module
</a>

The [DDRcore](@ref DDRcore) module provides classes and functions to implement the discrete de Rham sequence.


<a name="schemes">
\section schemes Schemes
</a>

The following schemes are currently available in HArD::Core3D. The Hybrid High-Order schemes follow the implementation principles described in Appendix B of the book available at https://hal.archives-ouvertes.fr/hal-02151813.

 - [HHO_diffusion](@ref HArDCore3D::HHO_Diffusion): Hybrid High-Order (HHO) for \f$-\mathrm{div}(K\nabla u)=f\f$, for Dirichlet, Neumann or mixed boundary conditions, with \f$K\f$ a diffusion tensor that is piecewise constant on the mesh.

 - [HHO_locvardiff](@ref HArDCore3D::HHO_LocVarDiff): HHO for \f$-\mathrm{div}(K\nabla u)=f\f$, for Dirichlet, Neumann or mixed boundary conditions, with \f$K\f$ a diffusion tensor that can vary in each cell.

 - [HHO_diffadvecreac](@ref HHO_DiffAdvecReac): HHO for \f$-\mathrm{div}(K\nabla u+\beta u)+\mu u=f\f$, for Dirichlet or mixed boundary conditions, with \f$K\f$ a diffusion tensor that can vary in each cell.

 - [DDR_magnetostatic](@ref DDR_magnetostatic): Discrete De Rham (DDR) scheme for the magnetostatic problem.

The directory `runs` contains BASH to run series of tests on families of meshes. The files `data.sh` describe the parameters of the test cases (polynomial degrees, boundary conditions, mesh families, etc.). The script produces results in the `output` directory, including a pdf file `rate.pdf` describing the rates of convergence in various energy norms.

To run the scripts as they are, you will need `pdflatex`.




