#include "mesh_reader.hpp"
#include <algorithm> 
#include <memory> // for std::unique_pointer
#include <mesh.hpp>

#include <stdlib.h>     /* exit, EXIT_FAILURE */

#ifndef MESH_BUILDER_HPP
#define MESH_BUILDER_HPP


namespace HArDCore3D {

// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/// The MeshBuilder class provides build tools to create a full mesh with all connectivities
class MeshBuilder {
public:
    /**
    * Constructor for MeshBuilder.
    */
    MeshBuilder();
    /**
    * Overloaded constructor for MeshBuilder so import_mesh can be called from build_the_mesh().
    */
    MeshBuilder(const std::string mesh_file);
    
    /**
     *  Build mesh
     */
    std::unique_ptr<Mesh> build_the_mesh();

private:
    void build_boundary(Mesh* mesh);  
    const std::string _mesh_file;
};


} // end namespace HArDCore3D
#endif 

