#include "vtu_writer.hpp"
#include <Eigen/Dense>
#include <mesh_builder.hpp>
#include <iostream>
#include <iomanip>


/*!
 * \addtogroup TransformMeshes
 * @{
 */

constexpr size_t dimspace=3;
//typedef Eigen::Vector3d VectorRd;
using namespace HArDCore3D;

/// Given a family of points, returns the maximum index N such that the first N points are co-planar
size_t coplanar_idx(const std::vector<VectorRd> & points)
{
  assert(points.size() > 2);

  // The tolerance for the test should depend on the local scale (max distance between any two points)
  double tol = 0.;
  for (size_t i=0; i<points.size(); ++i){
    for (size_t j=i+1; j<points.size(); ++j){
      tol = std::max(tol, (points[i]-points[j]).norm());
    }
  }
  tol *= 1e-8;
  
  // Find N by looking for families of vectors (vi-v0)_i of maximal rank
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(2, dimspace);
  M.row(0) << points[1]-points[0];
  M.row(1) << points[2]-points[0];
  size_t nrows=2;  
  bool found=false;
  do {
    if (M.fullPivLu().rank() == 2){
      found=true;
    }else{
      nrows++;
      M.conservativeResize(nrows, Eigen::NoChange_t());
      M.row(nrows-1) = points[nrows]-points[0];
    }
  }
  while(!found && nrows<points.size());
  
  return nrows;
}


/// Main executable (FlattenFaces) to try and flatten the faces in a mesh
int main(const int argc, const char *argv[])
{
    if (argc == 1){
      std::cout << "Attempts to flatten the faces in the mesh, writes RF files with the flattened mesh.\n Usage: flatten-faces <name of mesh>\n\n This is not a completely full-proof method, it will fail if the faces have re-entrant corners etc. Read the source, section 'Algorithm'." << std::endl;
      exit(0);
    }
    std::cout << "Mesh: " << argv[1] << std::endl;
    const std::string mesh_file = argv[1];

    // Build the mesh    
    HArDCore3D::MeshBuilder meshbuilder = HArDCore3D::MeshBuilder(mesh_file);
    std::unique_ptr<HArDCore3D::Mesh> mesh_ptr = meshbuilder.build_the_mesh();

    // Required information to write the RF files
    std::vector<VectorRd> nodes;       // List of node coordinates
    std::vector<std::vector<std::vector<size_t>>> cell_faces; // Each cell is described by a list of face, with each face being a list of nodes
    
    // Fill in nodes with existing vertices (none will be removed)
    for (size_t iV = 0; iV < mesh_ptr->n_vertices(); iV++){
      nodes.push_back(mesh_ptr->vertex(iV)->coords());
    }
    
    // For each face, we create flat subfaces (identified by their node indices)
    std::vector<std::vector<std::vector<size_t>>> flat_face(mesh_ptr->n_faces());
    for (size_t iF = 0; iF < mesh_ptr->n_faces(); ++iF){
      const Face * F = mesh_ptr->face(iF);
      if (F->is_flat()){
        // Face is flat, only one subface
        flat_face[iF].resize(1);
        flat_face[iF][0].resize(F->n_vertices());
        for (size_t iV=0; iV < F->n_vertices(); ++iV){
          flat_face[iF][0][iV] = F->vertex(iV)->global_index();
        }
      }else{
        // ALGORITHM:
        // The face is not flat we cut it by the following process:
        //    - Starting from the first vertex, we find the maximum number of vertices that are co-planar
        //    - We use those to create a subface
        //    - Then we re-start from the last vertex taken to form that subface.

        // List of coordinates of vertices of F, the last one repeated at the end to loop
        std::vector<VectorRd> coords_vertices_F;
        for (size_t iV = 0; iV < F->n_vertices(); iV++){
          coords_vertices_F.push_back(F->vertex(iV)->coords());
        } 
        coords_vertices_F.push_back(F->vertex(0)->coords());

        // We create subfaces starting from one vertex (at position p1), grabbing the maximum number of vertices from p1 that
        // are coplanar, and creating the subface from them; then we re-start from this position
        size_t p1=0;
        do {
          std::vector<VectorRd> tail_coords_vertices_F(coords_vertices_F.begin()+p1, coords_vertices_F.end());
          size_t N = coplanar_idx(tail_coords_vertices_F);
          std::vector<size_t> subface(N);
          for (size_t i=0; i<N; ++i){
            subface[i] = F->vertex( (p1 + i) % F->n_vertices())->global_index();
          }
          flat_face[iF].push_back(subface);
          p1 += N-1;
        } while (p1 < coords_vertices_F.size()-1);
        
      }
    }    
    
    // Create the vectors of cell_faces
    cell_faces.resize(mesh_ptr->n_cells());
    for (size_t iT = 0; iT < mesh_ptr->n_cells(); ++iT){
      const Cell * T = mesh_ptr->cell(iT);
      for (size_t iF = 0; iF < T->n_faces(); ++iF){
        size_t global_iF = T->face(iF)->global_index();
        for (size_t subF = 0; subF < flat_face[global_iF].size(); ++subF){
          cell_faces[iT].push_back(flat_face[global_iF][subF]);
        }
      }
    }
  
    //------
    // File to write RF
    //------
    std::string filename_core = mesh_file.substr(mesh_file.find_last_of("/\\") + 1);
    std::string filename_node = filename_core + "-flat.node";
    std::string filename_ele = filename_core + "-flat.ele";

    // Write .node and .ele files
    std::cout << "Writing .node file" << std::endl;
    std::ofstream outNode(filename_node.c_str());
    outNode << "# *.node file of 3D-mesh in REGN_FACE format" << std::endl;
    outNode << "# " << filename_node.c_str() << " created from " << filename_core << std::endl;
    outNode << nodes.size() << " " << "3  0  0" << std::endl;
    for (size_t iV = 0; iV < nodes.size(); iV++){
      outNode << std::setprecision(16) << iV;
      outNode << " ";
      outNode << std::setprecision(16) << nodes[iV](0);
      outNode << " ";
      outNode << std::setprecision(16) << nodes[iV](1);
      outNode << " ";
      outNode << std::setprecision(16) << nodes[iV](2);
      outNode << std::endl;
    }
    outNode.close();

    // write .ele file
    std::cout << "Writing .ele file" << std::endl;
    std::ofstream outEle(filename_ele.c_str());
    outEle << "# *.ele file of 3D-mesh in REGN_FACE format" << std::endl;
    outEle << "# " << filename_ele.c_str() << " created from " << filename_core << std::endl;
    outEle << cell_faces.size() << " " << " 0" << std::endl;
    for (size_t iT = 0; iT < cell_faces.size(); ++iT){
      outEle << iT << " " << cell_faces[iT].size() << std::endl;
      for (size_t iF = 0; iF < cell_faces[iT].size(); ++iF){
        outEle << iF << " " << cell_faces[iT][iF].size() << " ";
        for (size_t iV = 0; iV < cell_faces[iT][iF].size(); ++iV){
          outEle << cell_faces[iT][iF][iV] << " ";
        }
        outEle << std::endl;
      }
    }
    outEle.close();
        
    return 0;
}
