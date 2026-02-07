// Classes and methods to output VTU files for selected faces
//
// Provides method to plot .vtu file from selected faces and face-centered values.
//
// Author: Adapted for face plotting

#ifndef _VTU_FACEWRITER_HPP
#define _VTU_FACEWRITER_HPP

#include <memory>
#include <string>
#include <vector>
#include <map>
#include "mesh.hpp"

#include <Eigen/Dense>

namespace HArDCore3D {

/// The FaceVtuWriter class provides methods to plot selected faces and solutions on these faces
class FaceVtuWriter {

public:
 	/**
  * @brief Constructor for face writer 
  *
  * @param faces vector of pointers to the faces to plot
  */
  FaceVtuWriter(const std::vector<Face*>& faces);

	/// Writes file with a series of solutions (scalar if T=std::vector<double> (or Eigen::VectorXd), vector if T=std::vector<VectorRd>)
	template<typename T> 
	bool write_to_vtu(
				const std::string filename,                        ///< name of file to write to
				const std::vector<T> &sol_face,                    ///< each instance of T represents the values of one solution at the faces
				const std::vector<std::string> &sol_names = {}     ///< each string is the name of one solution
				)
  {
    FILE* pFile=fopen(filename.c_str(),"w");
	  write_header(pFile);
	  write_vertices(pFile);
	  
	  // If sol_names is too short, we re-create it
	  std::vector<std::string> names = sol_names;
	  if (names.size() < sol_face.size()){
	    names.resize(0);
	    for (size_t i=0; i < sol_face.size(); i++){
	      names.push_back("solution" + std::to_string(i));
	    }
	  }
	  
	  write_cells(pFile);
	  
    // Functions - NOTE: CellData instead of PointData!
    fprintf (pFile,"%s\n","         <CellData>");
	  for (size_t i=0; i < sol_face.size(); i++){
    	write_solution(pFile, sol_face[i], names[i]);
    }
    fprintf (pFile,"%s\n","         </CellData>");

	  write_footer(pFile);
	  fclose(pFile);
	  return true;
  }

  
  /// Overload to simplify the call when only one solution is involved
  template<typename T>
  inline bool write_to_vtu(
				const std::string file_name, ///< name of file to write to
				const T &sol_face,  ///< values of the solution at the faces
				const std::string &sol_name = "solution"  ///< name of the solution
				)
		{
		  return write_to_vtu(file_name, std::vector<T> {sol_face}, std::vector<std::string> {sol_name});
		}

	bool write_header(FILE* pFile);
	bool write_vertices(FILE* pFile);
	bool write_solution(FILE* pFile, const std::vector<double> &sol_face, const std::string &name); // scalar solution
	bool write_solution(FILE* pFile, const Eigen::VectorXd &sol_face, const std::string &name); // scalar solution (for legacy)
	bool write_solution(FILE* pFile, const std::vector<VectorRd> &sol_face, const std::string &name); // vector-valued solution
	bool write_cells(FILE* pFile);
	bool write_footer(FILE* pFile);


private:
  std::vector<Face*> _faces;
	size_t nfaces;
	size_t nvertices;
	std::vector<int> vtk_type;						          // "vtk types" of each face (5=triangle, 9=quad, 7=polygon)
	std::vector<std::vector<size_t>> f_vertices; 		// vertex indices in each face (using local numbering)
	std::vector<size_t> f_offset;					          // face offsets
	std::map<Vertex*, size_t> vertex_map;           // map from vertex pointer to local index
	std::vector<Vertex*> vertices;                  // list of unique vertices
};

}// end of namespace HArDCore3D

#endif //_VTU_FACEWRITER_HPP
