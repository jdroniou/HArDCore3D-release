#include "vtu_writer.hpp"
#include <Eigen/Dense>
#include <mesh_builder.hpp>
#include <boost/program_options.hpp>

/*!
 * \addtogroup TransformMeshes
 * @{
 */

//*** Default mesh
const std::string mesh_dir = "../../../meshes/";
std::string default_mesh = mesh_dir + "Cubic-Cells/RF_fmt/gcube_2x2x2";

/// Main executable (RFToVTU) to transform an RF mesh into a VTU file
int main(const int argc, const char *argv[])
{
  // Program options
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Creates a .vtu file from a mesh file in RF format")
    ("mesh,m", boost::program_options::value<std::string>(), "Input mesh file (complete path but without extension)")
    ("outfile,o", boost::program_options::value<std::string>(), "Output mesh file (without extension)");


  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  // Display the help options
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  // Select the mesh
  std::string mesh_file = (vm.count("mesh") ? vm["mesh"].as<std::string>() : default_mesh);

  std::cout << "Mesh: " << mesh_file << std::endl;

  // Build the mesh
  HArDCore3D::MeshBuilder meshbuilder = HArDCore3D::MeshBuilder(mesh_file);
  std::unique_ptr<HArDCore3D::Mesh> mesh_ptr = meshbuilder.build_the_mesh();

  // Write the mesh
  const std::string filename = (vm.count("outfile") ? vm["outfile"].as<std::string>() : mesh_file.substr(mesh_file.find_last_of("/\\") + 1)) + ".vtu";
  HArDCore3D::VtuWriter plotdata(mesh_ptr.get());
  Eigen::VectorXd zero_vec = Eigen::VectorXd::Zero(mesh_ptr -> n_vertices());
  std::cout << "Writing vtu file to " << filename << std::endl;
  plotdata.write_to_vtu(filename, zero_vec);

  return 0;
}
