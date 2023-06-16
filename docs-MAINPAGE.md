
HArD::Core (sources: https://github.com/jdroniou/HArDCore) provides a suite of C++ tools to implement numerical schemes whose unknowns are polynomials in the cells, on the edges, and on the faces. The focus is on dealing on generic polytopal meshes. This documentation addresses the 3D version of HArD::Core, but similar principles are valid for the 2D version. Transferring a scheme's implementation from 3D to 2D or vice-versa is very straightforward, provided that the scheme's mathematical definition does not depend on the dimension and that the generic types provided in `basis.hpp` and `MeshObject.hpp` are used; see readme file of the HArD::Core github depository https://github.com/jdroniou/HArDCore.

\tableofcontents

You will find on this page the following sections:

* [Build instructions](#build) -- How to build the library.
* [Mesh module](#mesh) -- Principle to load and handle a mesh and its geometric entities, together with some functions to transform/handle meshes.
* [Common module](#common) -- Polynomial bases, and other generic simple helper functions and structures.
* [Quadratures](#quad_rules) -- Quadrature rules to integrate generic functions, and cost-effective integration methods to compute Gram matrices of polynomial functions.
* [HybridCore](#hybridcore) -- Create polynomial basis functions on mesh geometric entities, and a vector structure of degrees of freedom. Mostly useful for Hybrid High-Order (HHO) and similar methods; a little bit deprecated, [HHOSpace](@ref HHOSpace) is probably better.
* [HHO3D](#hho3D) -- Core methods to implement HHO schemes.
* [DOFSpace](#dofspace) -- Generic classes to handle methods with degrees of freedom at the vertices, edges, faces and cells of a mesh. The following modules are based on these classes:
  - [HHOSpace](#hhospace) -- Methods (basis functions, discrete spaces and operators) for the Hybrid High-Order method.
  - [DDRCore](#ddr) -- Methods (basis functions, discrete spaces and operators) for the Discrete De Rham sequence (DDR). This is extended in [LADDRCore](@ref LADDRCore) to cover Lie-algebra valued functions.
* [Schemes](#schemes) -- List of schemes currently implemented in HArDCore3D.

<a name="build">
\section build Build instructions
</a>

\subsection buildlib Building the libraries and the schemes

To build the libraries and implemented schemes, the minimal requirements are:

* CMake version 2.6 or above (https://cmake.org/)
* A C++ compiler that supports the C++14 standard, eg. GCC (https://gcc.gnu.org/) or Clang (https://clang.llvm.org/).
* Eigen C++ library, version 3.3 or above (http://eigen.tuxfamily.org/)
* The following Boost C++ libraries (http://www.boost.org/): filesystem, program options, timer, chrono.

Make sure that you have the development version of boost installed. On Linux, install `libboost-dev`, `libboost-filesystem-dev`, `libboost-program-options-dev`, `libboost-chrono-dev` and `libboost-timer-dev` from your package manager.

The linear systems resulting from the assembled scheme are solved using the BiCGStab implementation of Eigen. Alternatives are also provided, but require additional libraries (UMFPACK, SUPERLU, etc.). Some schemes have procedures to compute condition numbers of matrices; these require the Spectra library (https://spectralib.org/). See the main CMakeLists.txt file for the libraries the code attempts to find, and the 'HINTS' section of individual CMake/*.cmake files for where it tries to find them.

Once you have installed all of the required dependencies, set up the build directory and generate the build files by running the following from the repository root:

```
mkdir build
cd build
cmake ..
make
```

After this, `build/Schemes` will contain the executables (e.g. `hho-diffusion`) to run the schemes. These executables need to access the meshes, which they should naturally find if you have left the `meshes` directory at the root of the project's files.



\subsection doco Building the Documentation

The mesh documentation is built with Doxygen (http://www.stack.nl/~dimitri/doxygen/). If you are reading this then somebody has already built it for you. If you modify the code and wish to rebuild the documentation, simply run `doxygen` from the root directory. The HTML version of the documentation is generated inside `docs` and the LaTeX version is generated inside `docs/latex` and can be compiled using the generated Makefile.


<a name="mesh">
\section mesh Mesh module
</a>

\subsection meshpple Principles

After it is loaded, the mesh is represented by typedefs of [MeshObject](@ref MeshND::MeshObject) describing a Vertex, an Edge, a Face, and a Cell. Each of these classes contains methods to access useful information for the corresponding element, including other geometrical quantities it is related to. The mesh itself is represented by an element of the [Mesh](@ref MeshND::Mesh) class with methods to access all the vertices, edges, faces and cells (or a particular vertex, edge, face or cell). 


For example, if `mesh_ptr` is a pointer to a [Mesh](@ref MeshND::Mesh) instance, the lines
\code{.cpp}
Vertex* vertex = mesh_ptr->vertex(5);

Eigen::Vector3d vert_coord = vertex->coords()
\endcode
store the coordinates of the fifth vertex into the Eigen vector `vert_coord`. As a generic rule, all geometrical vectors are `Eigen::Vector3d`. We also use `Eigen::Vector{3,X}d` and `Eigen::Matrix{3,X}d` for objects on which linear algebraic operations are performed. Lists (e.g. of cells, of functions...) are usually instances of `std::vector<...>`. Finally, `Eigen::multi_array` is used for lists of values of functions at quadrature nodes.

Here is an example that loops over all cells, grabs all the faces of the cell, and loops over these faces to output their diameter. Here, `mesh_ptr` is a pointer to the mesh.

\code{.cpp}
// Loop over all cells of the mesh
for (size_t iT = 0; iT < mesh_ptr->n_cells() iT++) {

	// We grab the faces of the iT-th cell
	std::vector<Face *> faces = mesh_ptr->cell(iT)->get_faces();

	// Loop over the faces of the cell
	for (size_t ilF = 0; ilF < cell->n_faces(); ilF++) {

		// Write the face diameter on the standard output
		std::cout << "The diameter of face " << ilF+1 << " in cell " << iT+1 << " is: " << faces(ilF)->diam() << "\n";
	}

}
\endcode

There is no direct access from a high-level geometrical entity to elements purely associated with lower-level entities. For example, if `mesh_ptr` is a pointer to the mesh, there is no direct method to access the coordinates of the i-th vertex of the mesh (no `mesh_ptr->coords_vertex()` exists). Instead, this is done through `mesh_ptr->vertex(i)->coords()`. This choice is deliberate as it preserves the logical organisation of the data structure, and facilitates the memorisation of the various methods. Of course, writing a wrapper providing such a direct access is easy...


\subsection loading_mesh Loading a mesh

HArDCore3D can read meshes in `RF` format. Previous versions could read `TG` and `MSH` files and were based on G. Manzini's mesh library https://github.com/gmanzini-LANL/PDE-Mesh-Manager. From Version 4.1, HArDCore3D uses an independent mesh reader written by L. Yemm.

A mesh file must be read using an instance of the `meshbuilder` class, and then built using `build_the_mesh`.  A working example is given below (assuming the executable will be in `build/Schemes` for example).

\code{.cpp}
#include "mesh.hpp"
#include "mesh_builder.hpp"

using namespace HArDCore3D;

int main() {

  // Mesh file to read
  std::string mesh_file = "../../meshes/Voro-small-0/RF_fmt/voro-4";

  // Build the mesh
  MeshBuilder meshbuilder = MeshBuilder(mesh_file);
  std::unique_ptr<Mesh> mesh_ptr = meshbuilder.build_the_mesh();

  std::cout << "There are " << mesh_ptr->n_cells() << " cells in the mesh.\n";
}
\endcode

Note that the builder returns a `unique_ptr` for the mesh. This ensures that, at the end of the execution, the mesh destructor is called (which destroys all cells, faces, edges, vertices...). Some classes and functions use a raw pointer to the mesh, so the `.get()` method should be used when passing the mesh as argument to the class constructors or functions.

<i>Note</i>: the mesh formats allow for meshes with very generic polygonal cells, including non-convex cells.
However, the builder assumes that each cell is star-shaped with respect to the isobarycenter of its vertices -- otherwise, the calculation of the center of mass may be incorrect. Similarly, the quadrature rules (see [Quadrature rules](#quad_rules)) assume that each cell is star-shaped with respect to its center of mass.

\subsection transform_meshes Other codes associated with mesh handling.

Several other codes (in src/Mesh/Transforms) are provided to manipulate meshes:
* Create agglomerated meshes [MeshCoarsen](@ref MeshCoarsen.cpp),
* Check that an RF file produces a valid mesh [CheckMesh](@ref CheckMesh.cpp),
* Make flat faces in a mesh where some are not (e.g., generic hexahedral meshes) [MakeFlatFaces](@ref MakeFlatFaces.cpp),
* Move the vertices of a mesh [MoveVertices](@ref MoveVertices.cpp),
* Various transformations of mesh format:
  - Gmsh .msh file (version 4.11) into an RF mesh file [MSHToRF](@ref MSHToRF.cpp).
  - RF mesh into a .vtu file [RFToVTU](@ref RFToVTU.cpp),
  - RF mesh into an FVCA10 format .msh file [RFToFVCA10format](@ref RFToFVCA10format.cpp),

<a name="common">
\section common_module Common module
</a>

The main classes in the [Common](@ref Common) describe polynomial basis functions on a cell, face, or edge. These could be bases of full polynomial spaces \f$\mathbb{P}^k\f$, or other related spaces (vector-valued polynomials, subspaces, image of gradient, image of curl, complements, etc.). The underlying basis functions are monomial ([MonomialScalarBasisCell](@ref HArDCore3D::MonomialScalarBasisCell) and [MonomialScalarBasisFace](@ref HArDCore3D::MonomialScalarBasisFace)), but derived bases (or set of non-necessarily linearly independent polynomial functions) can be handled through various classes, such as [Family](@ref HArDCore3D::Family), [TensorizedVectorFamily](@ref HArDCore3D::TensorizedVectorFamily), [GradientBasis](@ref HArDCore3D::GradientBasis), etc.

Free functions are also available to compute basis functions at quadrature nodes (using the [Quadrature rules](#quad_rules) module), orthonormalise basis functions, and compute Gram-like matrices between various families of functions. These matrices are essential in the design of high-order methods on polytopal meshes. Again, see examples in the [HybridCore](#hybridcore), [DDRCore](#ddr) and the various schemes built on them.

This module also contains:
 - [PolynomialSpaceDimension](@ref HArDCore3D::PolynomialSpaceDimension): structure to compute the dimensions of various polynomial spaces on edges, faces and cell,
 - [DOFSpace](@ref HArDCore3D::DOFSpace): class to access the local degrees of freedom associated with a geometric entity (vertex, edge, face, or cell) and all its associated entities of smaller dimension. This class determines how the local degrees of freedom are ordered (in the current setting, it's by increasing dimension of the associated geometric entities: DOFs of vertices, DOFs of edges, DOFs of faces and finally DOFs of cells). 

<a name="quad_rules">
\section quad_rules Integration over mesh geometric entities
</a>

\subsection usage_quad Generic quadrature rules.

HArD::Core deals with quite arbitrary cell geometries. As a consequence, no reference element can be used, and the quadrature rules have to be adjusted to each particular cell/face/edge. For the cells, for example, this is done by partitioning each cell into tetrahedra and by using classical quadrature rules on tetrahedras. The choice was also made not to pre-compute all quadrature rules for all cells, faces and edges, but rather to compute them -- with a locally chosen degree of exactness -- when needed in the code. To reduce the computational cost, quadrature rules -- and the values of basis functions at quadrature nodes -- should only be computed once when looping over each cell, before being used, e.g., to form mass matrices.

The [Quadratures](@ref Quadratures) module provides routines to do that. The key method, [generate_quadrature_rule(CFE,doe)](@ref HArDCore3D::generate_quadrature_rule), calculates quadrature nodes and weights, exact up to the polynomial degree `doe`, for an integration over cell/face/edge CFE (passed as a reference of [Cell](@ref HArDCore3D::Cell), [Face](@ref HArDCore3D::Face) or [Edge](@ref HArDCore3D::Edge) class). At present, the quadrature rules available in the code support a total degree of exactness in the cells up to 14 in the cells and 20 on the faces and the edges (the quadrature rules on the faces come from [John Burkardt's implementation of the Dunavant rules](https://people.sc.fsu.edu/~jburkardt/cpp_src/triangle_dunavant_rule/triangle_dunavant_rule.html)). The generated quadrature rule is stored in a structure [QuadratureRule](@ref HArDCore3D::QuadratureRule), which is a vector of quadrature nodes (weights and position).


The [evaluate_quad](@ref HArDCore3D::evaluate_quad) template function evaluate basis functions (their value, gradients, etc.) at provided quadrature nodes. The `boost::multi_array` provided by this function can then be passed to [compute_gram_matrix](@ref HArDCore3D::compute_gram_matrix) to create a matrix of inner products of two families of basis functions. Here is an example.

\code{.cpp}
// Create basis (f_1,...,f_r) of degree k in cell T
MonomialScalarBasisCell basisT(T, k);
// Create quadrature rules of degree 2*k in cell T
QuadratureRule quadT = generate_quadrature_rule(T, 2*k);
// Compute values of gradients of basis functions at the quadrature nodes
boost::multi_array<VectorRd, 2> gradbasis_on_quadT = evaluate_quad<Gradient>::compute(basisT, quadT);
// Create Gram-like matrix (here, a stiffness matix) of (\nabla f_i,\nabla f_j)_ij
Eigen::MatrixXd M = compute_gram_matrix(gradbasis_on_quadT, quadT);
\endcode

Note the usage of the type `VectorRd` defined in `basis.hpp`, which enables for a dimension-independent piece of code (easier to adapt to the 2D case). This procedure can also be applied, e.g., to cell basis functions on face quadrature nodes, etc. Additionally, the values at quadrature nodes obtained via [evaluate_quad](@ref HArDCore3D::evaluate_quad) can be transformed using [transform_values_quad](@ref HArDCore3D::transform_values_quad) (see also [scalar_product](@ref HArDCore3D::scalar_product) and [vector_product](@ref HArDCore3D::vector_product)); this gives for example an easy way of constructing the values at quadrature nodes on a face of normal or tangential traces of cell polynomials.

\code{.cpp}
// Create a monomial basis in a cell T of degree k)
MonomialScalarBasisCell basis_Pk_T(T, k);
// Tensorised the basis into vector-valued polynomials
TensorizedVectorFamily<MonomialScalarBasisCell, 3> basis_Pk3_T(basis_Pk_T);
// Create quadrature nodes on a face F
QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2*k);
// Evaluate the vector-valued polynomials at the quadrature nodes, and then take the normal and tangential components of these evaluations
boost::multi_array<VectorRd, 2> basis_Pk3_T_quad = evaluate_quad<Function>::compute(basis_Pk3_T, quad_2k_F);
VectorRd nF = F.normal();
boost::multi_array<double, 2> basis_Pk3_T_quad_nF = scalar_product(basis_Pk3_T_quad, nF);
boost::multi_array<VectorRd, 2> basis_Pk3_T_quad_tangentF = 
          transform_values_quad<VectorRd>(basis_Pk3_T_quad, [&nF](const VectorRd &z)->VectorRd { return z-(z.dot(nF))*nF;});
\endcode

<a name="quad_rules">
\section monomial_integration Homogeneous Numerical Integration (HNI) of polynomials, and Gram matrices
</a>

Because of the generality of some polytopal elements and the resulting need to generate quadrature rules by splitting the elements into tetrahedra, quadrature rules can very quickly reach a very large number of nodes (more than 1000 on some Voronoi cells for higher degrees of exactness). The calculation of integrals, and in particular of Gram matrices, then represent a very large amount of the computational time. To remediate this, we have implemented in HArDCore3D the HNI approach for monomials of https://doi.org/10.1007/s00466-015-1213-7. Since all the polynomial bases in the [Common](#common) module are built from monomial polynomials, this technique actually enables us to very efficiently compute Gram matrices of many of these bases.

The function [GramMatrix](@ref HArDCore3D::GramMatrix) precisely takes care of computing gram matrices of bases and derived sets of cell- and face-polynomials. At the moment, it properly manages all the following cases.

<table>
<caption id="GM_cell">Pairs of cell polynomial bases handled by GramMatrix <i>(swapped bases and all derived bases obtained through ShiftedBasis, RestrictedBasis or Family applied to the left of these bases are also ok)</i></caption>
<tr>
  <td rowspan="3">Scalar bases:</td><td>MonomialScalarBasisCell || MonomialScalarBasisCell</td><td></td>
</tr>
<tr>
  <td>DivergenceBasis<T> || MonomialScalarBasisCell </td><td><i>where T=TensorizedVectorFamily or RolyComplBasisCell</i></td>
</tr>
<tr>
  <td>DivergenceBasis<T1> || DivergenceBasis<T2> </td><td><i>where T1,T2=TensorizedVectorFamily or RolyComplBasisCell</i></td>
</tr>
<tr>
  <td rowspan="5">Vector bases:</td><td>TensorizedVectorFamily || T </td><td><i>where T is any basis with rank=Vector (TensorizedVectorFamily, GradientBasis, GolyComplBasisCell, etc.)</i></td>
</tr>
<tr>
  <td>GradientBasis<T1> || GradientBasis<T2> </td><td><i>where T1, T2 are any scalar bases</i></td>
</tr>
<tr>
  <td>CurlBasis<T1> || CurlBasis<T2> </td><td><i>where T1,T2=TensorizedVectorFamily or GolyComplBasisCell</i></td>
</tr>
<tr>
  <td>GolyComplBasisCell || GolyComplBasisCell</td><td></td>
</tr>
<tr>
  <td>RolyComplBasisCell || RolyComplBasisCell</td><td></td>
</tr>
<tr>
  <td rowspan="4">Matrix bases:</td><td>MatrixFamily<T1, N> || MatrixFamily<T2, N></td><td><i>where T1, T2 are any scalar bases, and N is an integer</i></td>
</tr>
<tr>
  <td>DivergenceBasis<MatrixFamily<T1, dimspace>> || TensorizedVectorFamily<T2, dimspace></td><td><i>where T1, T2 are any scalar bases</i></td>
</tr>
<tr>
  <td>GradientBasis<TensorizedVectorFamily<T1, dimspace>> || MatrixFamily<T2, dimspace></td><td><i>where T1, T2 are any scalar bases</i></td>
</tr>
<tr>
  <td>GradientBasis<TensorizedVectorFamily<T1, N>> || GradientBasis<TensorizedVectorFamily<T2, N>></td><td><i>where T1, T2 are any scalar bases, N is an integer</i></td>
</tr>
</table>

<table>
<caption id="GM_face">Pairs of face polynomial bases handled by GramMatrix <i>(swapped bases all derived bases obtained through ShiftedBasis, RestrictedBasis or Family applied to the left of these bases are also ok)</i></caption>
<tr>
  <td rowspan="2">Scalar bases:</td><td>MonomialScalarBasisFace || MonomialScalarBasisFace</td><td></td>
</tr>
<tr>
  <td>DivergenceBasis<T> || MonomialScalarBasisFace</td><td>where T=TangentFamily or RolyComplBasisFace (and derived)</td>
</tr>
<tr>
  <td rowspan="7">Vector bases:</td><td>TensorizedVectorFamily || TensorizedVectorFamily </td><td></td>
</tr>
<tr>
  <td>TangentFamily || TangentFamily </td><td></td>
</tr>
<tr>
  <td>CurlBasis || CurlBasis </td><td></td>
</tr>
<tr>
  <td>CurlBasis || TangentFamily </td><td></td>
</tr>
<tr>
  <td>TangentFamily || RolyComplBasisFace </td><td></td>
</tr>
<tr>
  <td>RolyComplBasisFace || RolyComplBasisFace </td><td></td>
</tr>
<tr>
  <td>GolyComplBasisCell || GolyComplBasisCell</td><td></td>
</tr>
</table>


<a name="hybridcore">
\section hybridcore Hybridcore module
</a>

This module encapsulates routines to create bases of polynomial spaces in each cell, on each face, and on each edge, and to manipulate discrete functions through the class [UVector](@ref HArDCore3D::UVector). 

\subsection basisfunc Basis functions

The instantiation of an [HybridCore](@ref HybridCore) class creates basis functions for the full polynomial spaces \f$\mathbb{P}^k\f$ in the cells, on the faces and on the edges, specifying the maximum degree required for each geometric entity. The basis functions are elements of the [Family](@ref HArDCore3D::Family) class and are accessed via the [CellBasis](@ref HArDCore3D::CellBasis), [FaceBasis](@ref HArDCore3D::FaceBasis) and [EdgeBasis](@ref HArDCore3D::EdgeBasis) method. For example, the following piece of code initialises an HybridCore instance with degrees \f$K+1\f$ in the cells, \f$K\f$ on the faces, and no edge basis functions, and access the value at some Eigen::Vector3d x of the i-th basis function on face iF, and the gradient of the j-th basis function in cell iT.

\code{.cpp}
  // Initialise the class
  HybridCore hho(mesh_ptr.get(), K+1, K, -1, use_threads, output);
  // Access value of face basis function
  double val_face = hho.FaceBasis(iF).function(i, x);
  // Access value of gradient of cell basis function
  Eigen::Vector3d grad_cell = hho.CellBasis(iT).gradient(j, x);
\endcode

The basis functions are hierarchical, which means that they are constructed by increasing degree. Hence, for example, if cell basis functions up to degree \f$K+1\f$ have been generated, a basis of the space of polynomials of degree \f$K\f$ in the cell is thus obtained by selecting the first [PolynomialSpaceDimension<Cell>::Poly(K)](@ref HArDCore3D::PolynomialSpaceDimension) cell basis functions.

The [UVector](@ref HArDCore3D::UVector) class describes coefficients on cell and face basis functions. The first coefficients correspond to cell basis functions, ordered by the cells themselves, and the last coefficients correspond to face basis functions. The methods in this class provide the surrounding structure to make sense of these coefficients (degrees of considered polynomial functions in cells/on faces, restrictions to a certain cell and its faces, etc.).

\subsection qr_hcore Quadrature rules evaluations in HybridCore

As explained above, the [evaluate_quad](@ref HArDCore3D::evaluate_quad) template enables the evaluation of basis functions (or their gradient, curl, or divergence) at pre-computed quadrature nodes. In the HybridCore module, however, the quadrature rules and values of basis functions (and gradients) at the quadrature nodes can be conveniently computed and stored using the [ElementQuad](@ref HArDCore3D::ElementQuad) class. Instantiating an element of this class on a cell loads these rules and values once, that can then be passed to several functions in charge of various calculations (e.g. one function computes the local cell contribution to the diffusion term, another function is in charge of computing the load term associated to the cell, etc.). This prevents recomputing these rules and values when needed by various functions. It works the following way:

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

(...)

\endcode

<a name="hho3D">
\section hho_3D HHO3D general module
</a>

The [HHO3D](@ref HHO3D) module provides typedefs, a class and functions to implement Hybrid High Order (HHO) schemes. Rules (functions) to create local bilinear forms (matrices) and loading terms (vectors) are passed to the HHO3D class, that takes care of the global assembly and solving of the system.

<a name="dofspace">
\section dofspace DOFSpace
</a>

The [LocalDOFSpace](@ref LocalDOFSpace) class provides methods to manage degrees of freedom (DOFs) attached to a cell, its faces, its edges and its vertices. The global version is [GlobalDOFSpace](@ref GlobalDOFSpace), and there is also a version, [VariableDOFSpace](@ref VariableDOFSpace), to handle cases where the number of DOFs varies, e.g., from one vertex to the other, from one edge to the other, etc. These classes provide methods to access the local or global degrees of freedom (determining the correct offset), etc.

The DOFs in these classes are organised by increasing order of the mesh entities dimensions (DOFs linked to all vertices, then DOFs linked to all edges, then DOFs linked to all faces, and finally all DOFs linked to cells). These DOFs only make sense when bases of polynomial spaces have been chosen on the geometric entities (such as some of the bases created by the HHOSpace or DDRCore classes below), and correspond then to the coefficients on these bases.

The following modules are based on these classes.

<a name="hhospace">
\subsection hho HHOSpace module
</a>

The [HHOSpace](@ref HHOSpace) module provides classes and functions to construct the local basis functions and operators for HHO space, including vector-valued spaces.
 
<a name="ddr">
\subsection ddr DDRCore module
</a>

The [DDRCore](@ref DDRCore) module provides classes and functions to implement the discrete de Rham (DDR) complex. This complex is based on spaces with unknowns on all geometric entities (vertices, edges, faces and cells), and discrete differential operators acting between these spaces. It is based on the principles detailed in https://doi.org/10.1007/s10208-021-09542-8 (see also the founding work https://doi.org/10.1142/S0218202520500372, and the arxiv versions https://arxiv.org/abs/2101.04940, https://arxiv.org/abs/1911.03616).

The main elements of the DDRCore module are:

  - [DDRCore](@ref HArDCore3D::DDRCore): class to construct bases of the local polynomial spaces, on all geometric entities, that are required for DDR schemes.
  - [XGrad](@ref HArDCore3D::XGrad), [XCurl](@ref HArDCore3D::XCurl) and [XDiv](@ref HArDCore3D::XDiv): classes to compute the discrete operators, potentials, interpolators and \f$L^2\f$-inner products associated with each space in the DDR complex. Each of these classes uses some of the bases built in DDRCore, and is built on a corresponding DDRSpace (which determines how the degrees of freedom, corresponding to the bases, are organised in the space).

Note that [DDRSpace](@ref HArDCore3D::DDRSpace) could potentially be used for generic schemes (not just based on the discrete de Rham sequence), perhaps with a different ordering of the DOFs.

<b>Important note</b>: <i>a directory "DDRCore-orth" can be found in the repository. It corresponds to the DDR spaces using orthogonal complements to the images of the gradient and curl operators on polynomial spaces, as described in https://doi.org/10.1142/S0218202520500372. This directory is not commented here, and its code is not maintained any longer. The directory "DDRCore" is the one referred to in this documentation, is based on the Kozsul complements of the images of gradient and curl, as in https://arxiv.org/abs/2101.04940, and is the one that is still maintained. DDRCore-orth is only provided for comparison purposes; the complements in this directory are much more expensive to create and manipulated, as explained in https://arxiv.org/abs/2101.04946. To compile a scheme using the orthogonal complements, simply modify the main CMakeLists.txt and change all "DDRCore" into "DDRCore-orth".</i>


<a name="schemes">
\section schemes Schemes
</a>

The following schemes are currently available in HArD::Core3D. The Hybrid High-Order schemes follow the implementation principles described in Appendix B of the book available at https://hal.archives-ouvertes.fr/hal-02151813.

 - [HHO_diffusion](@ref HArDCore3D::HHO_Diffusion): Hybrid High-Order (HHO) for \f$-\mathrm{div}(K\nabla u)=f\f$, for Dirichlet, Neumann or mixed boundary conditions, with \f$K\f$ a diffusion tensor that is piecewise constant on the mesh.

 - [HHO_locvardiff](@ref HArDCore3D::HHO_LocVarDiff): HHO for \f$-\mathrm{div}(K\nabla u)=f\f$, for Dirichlet, Neumann or mixed boundary conditions, with \f$K\f$ a diffusion tensor that can vary in each cell.

 - [HHO_diffadvecreac](@ref HHO_DiffAdvecReac): HHO for \f$-\mathrm{div}(K\nabla u+\beta u)+\mu u=f\f$, for Dirichlet or mixed boundary conditions, with \f$K\f$ a diffusion tensor that can vary in each cell.

 - [DDR_magnetostatic](@ref DDR_magnetostatic): Discrete De Rham (DDR) scheme, and serendipity version, for the magnetostatic problem, as per https://doi.org/10.1016/j.jcp.2020.109991 (but using Koszul complements).

 - [DDR_stokes](@ref DDR_stokes): Discrete De Rham (DDR) scheme for the %Stokes problem in curl-curl form. See also [SDDR_stokes](@ref SDDR_stokes) for the serendipity version.

 - [VEM_stokes](@ref VEM_stokes): Virtual Element Method (VEM) scheme for the %Stokes problem in curl-curl form.

 - [HHO_MHD](@ref HHO_MHD): HHO scheme for the MHD problem.

 - [DDR_yangmills](@ref DDR_yangmills): lowest-order DDR scheme for the Yang-Mills equations (based on [LADDR](@ref LADDRCore), the extension of DDR to Lie algebra valued fields). See also [SDDR_yangmills](@ref SDDR_yangmills) for the arbitrary order serendipity version.
 
 - [HHO_fullgradientdiff](@ref HHO_fullgradientdiff): HHO scheme with full gradient; similar to [HHO_locvardiff](@ref HArDCore3D::HHO_LocVarDiff) but implemented using the [HHOSpace](#hhospace) module instead of the [HybridCore](#hybridcore) module.
 
 - [HHO_brinkman](@ref HHO_brinkman): HHO scheme for the Brinkman equation (based on [HHOSpace](@ref HHOSpace)).

The directory `runs` contains BASH scripts to run series of tests on families of meshes. The files `data.sh` describe the parameters of the test cases (polynomial degrees, boundary conditions, mesh families, etc.). The script produces results in the `output` directory, shows the convergence rate in the standard console output, and creates a pdf file `rate.pdf` describing the rates of convergence in various energy norms (you will need `pdflatex` to create this pdf file; commenting out the corresponding line is fine, the pdf will simply not be create).




