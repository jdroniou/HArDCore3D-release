#include "import_mesh.hpp"
#include "vtu_writer.hpp"
#include <Eigen/Dense>
#include <mesh_builder.hpp>

int main(const int argc, const char *argv[])
{
    // const std::string mesh_dir = "../../typ2_meshes/";
    const std::string mesh_dir = "../../typ2_meshes/aglomerated/";
    const std::string mesh_file = argv[1];
    HArDCore2D::MeshBuilder builder = HArDCore2D::MeshBuilder(mesh_dir + mesh_file + ".typ2");
    std::unique_ptr<HArDCore2D::Mesh> mesh_ptr = builder.build_the_mesh();

    HArDCore2D::VtuWriter plotdata(mesh_ptr.get());
    Eigen::VectorXd zero_vec = Eigen::VectorXd::Zero(mesh_ptr -> n_vertices());
    plotdata.write_to_vtu(mesh_file + ".vtu", zero_vec);

    return 0;
}
