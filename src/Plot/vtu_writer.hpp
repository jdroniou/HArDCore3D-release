// Classes and methods to ouput various file formats
//
// Currently provides:
//  - Method to plot .vtu file from a mesh and node values.
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
*
*	This library was developed around HHO methods, although some parts of it have a more
* general purpose. If you use this code or part of it in a scientific publication, 
* please mention the following book as a reference for the underlying principles
* of HHO schemes:
*
* The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications. 
*  D. A. Di Pietro and J. Droniou. Modeling, Simulation and Applications, vol. 19. 
*  Springer International Publishing, 2020, xxxi + 525p. doi: 10.1007/978-3-030-37203-3. 
*  url: https://hal.archives-ouvertes.fr/hal-02151813.
*
*/

#ifndef _VTU_WRITER_HPP
#define _VTU_WRITER_HPP

#include <memory>
#include <string>
#include <vector>
#include "mesh.hpp"

#include <Eigen/Dense>


/*!	
*	@defgroup Plot 
* @brief Classes providing tools to create vtu files to visualise solutions
*/

namespace HArDCore3D {

/*!
*	@addtogroup Plot
* @{
*/

// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/// The VtuWriter class provides methods to plot a 3D mesh and a solution on this mesh
class VtuWriter {

public:
 	/**
  * @brief Constructor for mesh writer 
  *
  * @param mesh pointer to the mesh
  */
  VtuWriter(const Mesh* mesh);

	///Writes the plot file
	bool write_to_vtu(
				const std::string file_name, ///< name of file to write to
				const Eigen::VectorXd &sol_vertex 	///< vector of values of the solution at the mesh vertices
				);

	bool write_header(FILE* pFile);
	bool write_vertices(FILE* pFile);
	bool write_solution(FILE* pFile, const Eigen::VectorXd &sol_vertex);
	bool write_cells(FILE* pFile);
	bool write_footer(FILE* pFile);

	/// overloaded writer for the mesh alone
//	bool write_to_vtu(std::string file_name); 

private:
  const Mesh* _mesh;
	size_t ncells;
	std::vector<int> vtk_type;						// "vtk types" of each cell
	std::vector<std::vector<Vertex *>> c_vertices; 		// vertices in each cell
	std::vector<size_t> c_offset;					// cell offsets
	std::vector<std::vector<std::vector<Vertex *>>> c_f_vertices;		// list of vertices in each face of each cell
	std::vector<size_t> c_f_offset;		// cell-face offsets
};

}// end of namespace HArDCore3D

#endif //_VTU_WRITER
