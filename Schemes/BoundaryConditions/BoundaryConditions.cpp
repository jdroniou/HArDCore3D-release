// Class to provide description of boundary conditions
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "BoundaryConditions.hpp"

using namespace HArDCore3D;

// ----------------------------------------------------------------------------
//                          Implementation
// ----------------------------------------------------------------------------

// Class
BoundaryConditions::BoundaryConditions(const std::string bc_id, Mesh& mesh)
  : m_bc_id(bc_id),
    m_mesh(mesh),
    m_n_dir_faces(0)
  {
    // Compute number of Dirichlet faces
    for (Face* face : m_mesh.get_b_faces()){
      if (type(*face) == "dir"){
        m_n_dir_faces++;
      }
    }
  }


// Returns the type of BC of a face
const std::string BoundaryConditions::type(Face& face) const
  {
    if (face.is_boundary()){
      if (m_bc_id == "D"){
        return "dir";
      }else if (m_bc_id == "N"){
        return "neu";
      }else if (m_bc_id == "M0"){
        // Dirichlet at x=0, Neumann everywhere else.
        VectorRd v0 = face.vertex(0)->coords();
        VectorRd v1 = face.vertex(1)->coords();
        VectorRd xF = face.center_mass();
        double eps = 1e-8;
        if ( (std::abs(v0.x())<eps) && (std::abs(v1.x())<eps) && (std::abs(xF.x())<eps)){
          return "dir";
        }else{
          return "neu";
        }
      }
    }
    return "int";
  }


// Reorder faces
void BoundaryConditions::reorder_faces(const std::string pos)
  {
    // Create a vector of Dirichlet boundary faces, and a vector of other faces
    std::vector<size_t> dir_faces(m_mesh.n_faces(), 0);
    std::vector<size_t> nondir_faces(m_mesh.n_faces(), 0);
    size_t dir_idx = 0;
    size_t nondir_idx = 0;
    for (Face* face : m_mesh.get_faces()){
      if (type(*face) == "dir"){
        dir_faces[dir_idx] = face->global_index();
        dir_idx++;
      }else{
        nondir_faces[nondir_idx] = face->global_index();
        nondir_idx++;
      }
    }
    // check
    if (dir_idx + nondir_idx != m_mesh.n_faces()){
     std::cout << "Error during renumbering faces: " << dir_idx << ", " << nondir_idx << ", " << m_mesh.n_faces() << "\n";
     exit(1);
    }
    
    // Depending on "pos" we put the Dirichlet faces at the end or the start
    std::vector<size_t> new_to_old(m_mesh.n_faces(), 0);
    if (pos=="end"){
      for (size_t i=0; i < nondir_idx; i++){
        new_to_old[i] = nondir_faces[i];
      }     
      for (size_t i=nondir_idx; i < m_mesh.n_faces(); i++){
        new_to_old[i] = dir_faces[i-nondir_idx];
      }
    }else{
      for (size_t i=0; i < dir_idx; i++){
        new_to_old[i] = dir_faces[i];
      }     
      for (size_t i=dir_idx; i < m_mesh.n_faces(); i++){
        new_to_old[i] = nondir_faces[i-dir_idx];
      }
    }

    // Reordering
    m_mesh.renum('F', new_to_old);

 }


