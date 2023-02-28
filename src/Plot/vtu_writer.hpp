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
*  @defgroup Plot 
* @brief Classes providing tools to create vtu files to visualise solutions
*/

namespace HArDCore3D {

/*!
*  @addtogroup Plot
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

	/// Writes file with a series of solutions (scalar if T=std::vector<double> (or Eigen::VectorXd), vector if T=std::vector<VectorRd>)
	template<typename T> 
	bool write_to_vtu(
				const std::string filename,                        ///< name of file to write to
				const std::vector<T> &sol_vrtx,                    ///< each instance of T represents the values of one solution at the mesh vertices
				const std::vector<std::string> &sol_names = {}     ///< each string is the name of one solution
				)
  {
    FILE* pFile=fopen(filename.c_str(),"w");
	  write_header(pFile);
	  write_vertices(pFile);
	  
	  // If sol_names is too short, we re-create it
	  std::vector<std::string> names = sol_names;
	  if (names.size() < sol_vrtx.size()){
	    names.resize(0);
	    for (size_t i=0; i < sol_vrtx.size(); i++){
	      names.push_back("solution" + std::to_string(i));
	    }
	  }
	  
    // Functions
    fprintf (pFile,"%s\n","         <PointData>");
	  for (size_t i=0; i < sol_vrtx.size(); i++){
    	write_solution(pFile, sol_vrtx[i], names[i]);
    }
    fprintf (pFile,"%s\n","         </PointData>");

	  write_cells(pFile);
	  write_footer(pFile);
	  fclose(pFile);
	  return true;
  }

  
  /// Overload to simplify the call when only one solution is involved
  template<typename T>
  inline bool write_to_vtu(
				const std::string file_name, ///< name of file to write to
				const T &sol_vertex,  ///< values of the solution at the mesh vertices
				const std::string &sol_name = "solution"  ///< name of the solution
				)
		{
		  return write_to_vtu(file_name, std::vector<T> {sol_vertex}, std::vector<std::string> {sol_name});
		}

	bool write_header(FILE* pFile);
	bool write_vertices(FILE* pFile);
	bool write_solution(FILE* pFile, const std::vector<double> &sol_vertex, const std::string &name); // scalar solution
	bool write_solution(FILE* pFile, const Eigen::VectorXd &sol_vertex, const std::string &name); // scalar solution (for legacy, values at the vertices should really be std::vector)
	bool write_solution(FILE* pFile, const std::vector<VectorRd> &sol_vertex, const std::string &name); // vector-valued solution
	bool write_cells(FILE* pFile);
	bool write_footer(FILE* pFile);


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
