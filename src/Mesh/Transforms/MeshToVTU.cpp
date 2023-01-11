#include "vtu_writer.hpp"
#include <Eigen/Dense>
#include <mesh_builder.hpp>

/*!
 * \addtogroup TransformMeshes
 * @{
 */

/// Main executable (MeshToVTU) to transform an RF mesh into a VTU file
int main(const int argc, const char *argv[])
{
    if (argc == 1){
      std::cout << "Creates a .vtu file from a mesh file in RF format.\n Usage: rf-to-vtu <name of mesh>" << std::endl;
      exit(0);
    }
    std::cout << "Mesh: " << argv[1] << std::endl;
    const std::string mesh_file = argv[1];
    
    HArDCore3D::MeshBuilder meshbuilder = HArDCore3D::MeshBuilder(mesh_file);
    std::unique_ptr<HArDCore3D::Mesh> mesh_ptr = meshbuilder.build_the_mesh();

    HArDCore3D::VtuWriter plotdata(mesh_ptr.get());
    Eigen::VectorXd zero_vec = Eigen::VectorXd::Zero(mesh_ptr -> n_vertices());
    std::string filename = mesh_file.substr(mesh_file.find_last_of("/\\") + 1) + ".vtu";
    std::cout << "Writing vtu file to " << filename << std::endl;
    plotdata.write_to_vtu(filename, zero_vec);

    return 0;
}
