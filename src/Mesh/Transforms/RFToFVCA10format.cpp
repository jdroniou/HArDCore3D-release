#include "vtu_writer.hpp"
#include <Eigen/Dense>
#include <mesh_builder.hpp>
#include <boost/program_options.hpp>

/*!
 * \addtogroup TransformMeshes
 * @{
 */

using namespace HArDCore3D;

//*** Default mesh
const std::string mesh_dir = "../../../meshes/";
std::string default_mesh = mesh_dir + "Cubic-Cells/RF_fmt/gcube_2x2x2";

/// Main executable (RFToFVCA10format) to transform an RF mesh into a .msh FVCA10 format
int main(const int argc, const char *argv[])
{

  // Program options
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Creates a .msh (FVCA10 format) file from a mesh file in RF format")
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

  // Build mesh    
  HArDCore3D::MeshBuilder meshbuilder = HArDCore3D::MeshBuilder(mesh_file);
  std::unique_ptr<HArDCore3D::Mesh> mesh_ptr = meshbuilder.build_the_mesh();

  // Write file
  const std::string original_meshfile = mesh_file.substr(mesh_file.find_last_of("/\\") + 1);
  const std::string filename = (vm.count("outfile") ? vm["outfile"].as<std::string>() : original_meshfile) + ".msh";  
  
  std::cout << "Writing to file " << filename << std::endl;

  std::ofstream outFile(filename.c_str());
  outFile << "Mesh generated by MeshToFVCA10format.cpp, translated from RF mesh" << std::endl;
  outFile << "(from the HArDCore3D library)" << std::endl;
  outFile << "Original mesh: \n" << mesh_file << std::endl;
  outFile << "\n\n\nInformation:" << std::endl;
  outFile << "Number of vertices" << std::endl;
  outFile << "\t" << mesh_ptr->n_vertices() << std::endl;
  outFile << "Number of control volumes" << std::endl;
  outFile << "\t" << mesh_ptr->n_cells() << std::endl;
  outFile << "Number of faces" << std::endl;
  outFile << "\t" << mesh_ptr->n_faces() << std::endl;
  outFile << "Number of edges" << std::endl;
  outFile << "\t" << mesh_ptr->n_edges() << std::endl;

  outFile << "Vertices " << mesh_ptr->n_vertices() << std::endl;
  for (size_t iV=0; iV < mesh_ptr->n_vertices(); iV++){
    outFile << "\t" << mesh_ptr->vertex(iV)->coords().transpose() << std::endl;
  }

  outFile << "Volumes->Faces " << mesh_ptr->n_cells() << std::endl;
  for (size_t iT=0; iT < mesh_ptr->n_cells(); iT++){
    const Cell * T = mesh_ptr->cell(iT);
    outFile << "\t" << T->n_faces() << " ";
    for (size_t iF=0; iF < T->n_faces(); iF++){
      outFile << "\t" << T->face(iF)->global_index() + 1 << " ";
    }
    outFile << std::endl;
  }

  outFile << "Volumes->Vertices " << mesh_ptr->n_cells() << std::endl;
  for (size_t iT=0; iT < mesh_ptr->n_cells(); iT++){
    const Cell * T = mesh_ptr->cell(iT);
    outFile << "\t" << T->n_vertices() << " ";
    for (size_t iV=0; iV < T->n_vertices(); iV++){
      outFile << "\t" << T->vertex(iV)->global_index() + 1 << " ";
    }
    outFile << std::endl;
  }
   
  outFile << "Faces->Edges " << mesh_ptr->n_faces() << std::endl;
  for (size_t iF=0; iF < mesh_ptr->n_faces(); iF++){
    const Face * F = mesh_ptr->face(iF);
    outFile << "\t" << F->n_edges() << " ";
    for (size_t iE=0; iE < F->n_edges(); iE++){
      outFile << "\t" << F->edge(iE)->global_index() + 1 << " ";
    }
    outFile << std::endl;
  }

  outFile << "Faces->Vertices " << mesh_ptr->n_faces() << std::endl;
  for (size_t iF=0; iF < mesh_ptr->n_faces(); iF++){
    const Face * F = mesh_ptr->face(iF);
    outFile << "\t" << F->n_vertices() << " ";
    for (size_t iV=0; iV < F->n_vertices(); iV++){
      outFile << "\t" << F->vertex(iV)->global_index() + 1 << " ";
    }
    outFile << std::endl;
  }

  outFile << "Faces->Control Volumes " << mesh_ptr->n_faces() << std::endl;
  for (size_t iF=0; iF < mesh_ptr->n_faces(); iF++){
    const Face * F = mesh_ptr->face(iF);
    outFile << "\t" << F->cell(0)->global_index() + 1 << " ";
    if (F->is_boundary()){
      outFile << -1;
    }else{
      outFile << "\t" << F->cell(1)->global_index() + 1;
    }
    outFile << std::endl;
  }
  
  outFile << "Edges " << mesh_ptr->n_edges() << std::endl;
  for (size_t iE=0; iE < mesh_ptr->n_edges(); iE++){
    const Edge * E = mesh_ptr->edge(iE);
    outFile << "\t" << E->vertex(0)->global_index() + 1 << " " << E->vertex(1)->global_index() + 1 << std::endl;
  }

  outFile.close();
  
  return 0;
}
