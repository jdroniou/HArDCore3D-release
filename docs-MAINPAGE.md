
This is the documentation of the 3D version (sources: https://github.com/jdroniou/HArDCore3D-release/) of the HArD::Core suite of C++ tools for schemes whose unknowns are polynomials in the cells and on the edges (in 2D) or faces (in 3D). Most of the generic details described in the main page https://jdroniou.github.io/HArDCore2D-release/ of the 2D version are applicable to this 3D version. We refer to this page for the generic principles behind HArD::Core and we focus here on the specificities of the 3D version.

Because they are built on the same principles, transferring an implementation of a scheme from 2D to 3D or vice-versa is quite straightforward (if the scheme's own principles are dimension-independent, that is, its mathematical description is the same in 2D and 3D). Details to perform such a transfer are described in the readme file of the HArD::Core github depository
https://github.com/jdroniou/HArDCore.

You will find here:

\tableofcontents

* [The mesh structure](#mesh) -- The principles of the data structure representing the mesh, and how to load a mesh.
* [Schemes](#schemes) -- The list of schemes currently implemented in HArD::Core3D, and scripts to run them.

<a name="mesh">
\section mesh The mesh
</a>

\subsection meshpple Principles

As in 2D, after it is loaded the mesh is represented by classes describing its geometric elements of dimension 0, 1, 2 and 3. There is thus a `vertex` class, an `edge` class, a `face` class and a `mesh` class: [Vertex](@ref HArDCore3D::Vertex), [Edge](@ref HArDCore3D::Edge), [Face](@ref HArDCore3D::Face), and [Cell](@ref HArDCore3D::Cell). Each of these classes contains methods to access useful information for the corresponding element, including other geometrical quantities it is related to. The mesh itself is represented by an element of the [Mesh](@ref HArDCore3D::Mesh) class with methods to access all the vertices, edges, faces and cells (or a particular vertex, edge, face or cell). In this class, each cell has a unique identifier, numbered from 0. We refer to the main page of the documentation  https://jdroniou.github.io/HArDCore2D-release/ of the 2D version for the principles on how to manipulate these classes, and to the class descriptions for specific members of the 3D version.



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


<a name="schemes">
\section schemes Schemes
</a>

The schemes currently available in HArD::Core3D are:

 - [HHO_diffusion](@ref HArDCore3D::HHO_Diffusion): Hybrid High-Order for \f$-\mathrm{div}(K\nabla u)=f\f$, for Dirichlet or Neumann boundary conditions, with \f$K\f$ a diffusion tensor that is piecewise constant on the mesh.

 - [HHO_locvardiff](@ref HArDCore3D::HHO_LocVarDiff): Hybrid High-Order for \f$-\mathrm{div}(K\nabla u)=f\f$, for Dirichlet or Neumann boundary conditions, with \f$K\f$ a diffusion tensor that can vary in each cell.

The directory `runs` contains BASH to run series of tests on families of meshes. The files `data.sh` describe the parameters of the test cases (polynomial degrees, boundary conditions, mesh families, etc.). The script produces results in the `output` directory, including a pdf file `rate.pdf` describing the rates of convergence in various energy norms.

To run the scripts as they are, you will need `pdflatex` and a FORTRAN compiler, and to adjust the `Makefile` to your compiler, to run `compute_rates.f90` and compute the rates of convergence in the various norms.




