// Class to build the mesh data after having read the mesh file
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#ifndef MESH_BUILDER_HPP
#define MESH_BUILDER_HPP
//#include "mesh.hpp"
//#include "cell.hpp"
//#include "edge.hpp"
//#include "vertex.hpp"
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <memory>

namespace HArDCore3D { // Forward declaration
  class Mesh;
}

namespace HArDCore3D {

  /*!
   *       @addtogroup Mesh
   * @{
   */

  // ----------------------------------------------------------------------------
  //                            Class definition
  // ----------------------------------------------------------------------------

  /// The MeshBuilder class provides build tools to create a full mesh with all connectivities
  class MeshBuilder {
  public:
    /**
     * Constructor for MeshBuilder.
     */
    MeshBuilder(
		std::string mesh_file,  ///< the mesh file to read
		std::string mesh_type           ///< type of mesh file: TG, MSH, RF
                );
    /**
     * Build a mesh from a mesh file
     *
     * @return a pointer to the mesh that was build
     */
    std::unique_ptr<Mesh> build_the_mesh();  ///< construct the connectivity in the mesh

  private:
    std::string _mesh_file;  
    std::string _mesh_type;         
  };

  /*@}*/
}
#endif /* MESH_BUILDER_HPP */

