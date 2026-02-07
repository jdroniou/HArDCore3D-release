// Class to read a RF mesh file

#ifndef MESH_READER_FRACTURE_HPP
#define MESH_READER_FRACTURE_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <mesh.hpp>
namespace HArDCore3D
{

    // ----------------------------------------------------------------------------
    //                            Class definition
    // ----------------------------------------------------------------------------

    /// The MeshReaderRF class provides functions to read a RF mesh file
    class MeshReaderFracture
    {
    public:
        /**
    * Constructor for mesh reader
    *
    * @param file_name name of the file name, needs to include the full path
    */
        MeshReaderFracture(std::string file_name, Mesh& mesh); ///< class to read the cells and vertices in a file
        /**
    * Reads the file into the specified containers
    *
    * @param vertices reference to a vector to hold the vertices coordinates
    * @param faces reference to a vector to hold the cell indexes
    * @param cells reference to a vector to hold the cell centers coordinates
    */

    ~MeshReaderFracture() {}
        void read_frac_file(std::vector<std::vector<Vertex *>>& fracture_pointer);

    private:
        std::string _file_name; ///< name of the file being read
        Mesh& m_mesh;
    };

} // end namespace HArDCore3D
#endif
