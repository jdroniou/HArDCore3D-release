#include "vtu_writer.hpp"
#include <Eigen/Dense>
#include <mesh_builder.hpp>
#include <iostream>
#include <iomanip>
#include <random>
#include <boost/program_options.hpp>

/*!
 * \addtogroup TransformMeshes
 * @{
 */

using namespace HArDCore3D;

//*** Declaration of helper functions

/// Writes RF files from a mesh
void write_RF(HArDCore3D::Mesh * mesh, const std::string filename);

/// Move vertices randomly, except boundary and those with x=1/2 or z=1/2
void random_move(HArDCore3D::Mesh * mesh, double & factor);

/// Apply algebraic transformation to all vertices
void apply_transformation(HArDCore3D::Mesh * mesh, std::function<VectorRd(const VectorRd &)> trans);

//*** Default mesh
const std::string mesh_dir = "../../../meshes/";
std::string default_mesh = mesh_dir + "Cubic-Cells/RF_fmt/gcube_2x2x2";


//***
/// Main executable (MoveVertices) to modify a mesh by moving its vertices
int main(const int argc, const char *argv[])
{
  // Program options
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Apply a transformation to the vertices of the mesh, and write new RF files")
    ("mesh,m", boost::program_options::value<std::string>(), "Input mesh file (complete path but without extension)")
    ("outfile,o", boost::program_options::value<std::string>(), "Output mesh file (without extension)")
    ("typemove,t", boost::program_options::value<size_t>()->default_value(0), "Set the type of motion of the vertices")
    ("factor,f", boost::program_options::value<double>()->default_value(0.2), "By how much we move (proportion)");

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  // Display the help options
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  // Select the mesh and type of movement
  std::string mesh_file = (vm.count("mesh") ? vm["mesh"].as<std::string>() : default_mesh);
  size_t typemove = (vm.count("typemove") ? vm["typemove"].as<size_t>() : 0);
  double factor = (vm.count("factor") ? vm["factor"].as<double>() : 0.2);

  std::cout << "Mesh: " << mesh_file << std::endl;

  // Build the mesh    
  HArDCore3D::MeshBuilder meshbuilder = HArDCore3D::MeshBuilder(mesh_file);
  std::unique_ptr<HArDCore3D::Mesh> mesh_ptr = meshbuilder.build_the_mesh();
  
  switch (typemove) {
    case 0:
      // Move randomly except x=1/2, z=1/2 and brings to (-1,1)^3
      random_move(mesh_ptr.get(), factor);       
      apply_transformation(mesh_ptr.get(), [](const VectorRd & x)->VectorRd { return 2.*x - VectorRd(1.,1.,1.);}); 
      break;

  default:
    std::cerr << "[main] ERROR: Type of move unknown" << std::endl;
    exit(1);

  }
  
  // Write new mesh
  const std::string filename_core = mesh_file.substr(mesh_file.find_last_of("/\\") + 1);
  const std::string new_mesh_file = (vm.count("outfile") ? vm["outfile"].as<std::string>() : filename_core + "-moved");
  write_RF(mesh_ptr.get(), new_mesh_file);
      
  return 0;
}


//****
// Helper functions
//****

void write_RF(HArDCore3D::Mesh * mesh, const std::string filename)
{
  
  std::string filename_node = filename + ".node";
  std::string filename_ele = filename + ".ele";
  
  // write .node file
  std::cout << "Writing .node file" << std::endl;
  std::ofstream outNode(filename_node.c_str());
  outNode << "# *.node file of 3D-mesh in REGN_FACE format" << std::endl;
  outNode << mesh->n_vertices() << " " << "3  0  0" << std::endl;
  for (size_t iV = 0; iV < mesh->n_vertices(); iV++){
    outNode << std::setprecision(16) << iV;
    outNode << " ";
    outNode << std::setprecision(16) << mesh->vertex(iV)->coords().transpose();
    outNode << std::endl;
  }
  outNode.close();
   
  // write .ele file
  std::cout << "Writing .ele file" << std::endl;
  std::ofstream outEle(filename_ele.c_str());
  outEle << "# *.ele file of 3D-mesh in REGN_FACE format" << std::endl;
  outEle << mesh->n_cells() << " " << " 0" << std::endl;
  for (size_t iT = 0; iT < mesh->n_cells(); iT++){
    const Cell * T = mesh->cell(iT);

    // Header for the element
    outEle << iT << " " << T->n_faces() << std::endl;
    // Print out the faces
    for (size_t iF = 0; iF < T->n_faces(); iF++){
      const Face * F = T->face(iF);
      size_t nvertF = F->n_vertices();
      outEle << iF << " " << nvertF;
      for (size_t iV = 0; iV < nvertF; iV++){
        outEle << " " << F->vertex(iV)->global_index();
      }
      outEle << std::endl;
    }
  }
  
  outEle.close();

}

void random_move(HArDCore3D::Mesh * mesh, double & factor)
{
  // Tolerance for tests
  double tol = 1e-8 * mesh->h_max();
  
  // Random distrbution
  std::default_random_engine gen;
  std::uniform_real_distribution<double> dist(-1.0,1.0);

  // Minimum distance between V and all other vertices of the cells around
  std::vector<double> hV(mesh->n_vertices(), mesh->h_max());
  for (size_t iV=0; iV<mesh->n_vertices(); iV++){
    Vertex * V = mesh->vertex(iV);
    for (size_t iT=0; iT<V->n_cells(); iT++){
      const Cell * T = V->cell(iT);
      for (size_t ip=0; ip<T->n_vertices(); ip++){
        double distVp = (V->coords() - T->vertex(ip)->coords()).norm();
        if (distVp>tol){
          hV[iV] = std::min(hV[iV], distVp);
        }
      }
    }
  }
  
  for (size_t iV=0; iV<mesh->n_vertices(); iV++){
    Vertex * V = mesh->vertex(iV);

    // Only move internal points
    if (!V->is_boundary()){

      // Is V on one of the two fractures?
      bool is_V_fracture_x = (std::abs(V->coords()(0) - 0.5)<tol);
      bool is_V_fracture_z = (std::abs(V->coords()(2) - 0.5)<tol);

      // On a fracture, only move along that fracture
      VectorRd shift = VectorRd::Zero();
      if (is_V_fracture_x && is_V_fracture_z){
        // Both fractures, only move in y direction
        shift = factor * hV[iV] * VectorRd(0.,dist(gen),0.);
      }else if (is_V_fracture_x){
        // x fracture, only move in y,z
        shift = factor * hV[iV] * VectorRd(0.,dist(gen),dist(gen));
      }else if (is_V_fracture_z){
        // z fracture, only move in x,y
        shift = factor * hV[iV] * VectorRd(dist(gen),dist(gen),0.);
      }else{
        // Not on any fracture
        shift = factor * hV[iV] * VectorRd(dist(gen),dist(gen),dist(gen));
      }

      // Not on a fracture, move at random in any direction
      V->set_coords(V->coords() + shift);
    }

  }
  
}


void apply_transformation(HArDCore3D::Mesh * mesh, std::function<VectorRd(const VectorRd &)> trans)
{
  for (Vertex * V : mesh->get_vertices()){
    V->set_coords(trans(V->coords()));
  }
}

