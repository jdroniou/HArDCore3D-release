// Class to read a RF mesh file

#ifndef MESH_READER_HPP
#define MESH_READER_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
namespace HArDCore3D
{

    // ----------------------------------------------------------------------------
    //                            Class definition
    // ----------------------------------------------------------------------------

    /// The MeshReaderRF class provides functions to read a RF mesh file
    class MeshReaderRF
    {
    public:
        /**
    * Constructor for mesh reader
    *
    * @param file_name name of the file name, needs to include the full path
    */
        MeshReaderRF(std::string file_name); ///< class to read the cells and vertices in a file
        /**
    * Reads the file into the specified containers
    *
    * @param vertices reference to a vector to hold the vertices coordinates
    * @param faces reference to a vector to hold the cell indexes
    * @param cells reference to a vector to hold the cell centers coordinates
    */

    ~MeshReaderRF() {}
        void read_node_file(std::vector<std::vector<double>> &vertices);

        void read_ele_file(std::vector<std::vector<std::vector<std::size_t>>> &cells);

        void read_partition_file(std::vector<std::vector<std::size_t>> &partition);

    private:
        std::string _file_name; ///< name of the file being read
    };

}; // end namespace HArDCore3D
#endif
