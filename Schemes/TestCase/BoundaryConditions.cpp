// Class to provide description of boundary conditions
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "BoundaryConditions.hpp"
#include "vertex.hpp"
#include "basis.hpp" // for the VectorRd type

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
  };


// Returns the type of BC of a face
const std::string BoundaryConditions::type(const Face& face) const
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
void BoundaryConditions::reorder_faces()
  {
    // Create vector with all non-Dirichlet faces first, and all Dirichlet faces at the end
    std::vector<size_t> new_to_old(m_mesh.n_faces(), 0);
    // Index for non-Dirichlet faces start at 0 and increases, index for Dirichlet faces start at n_faces()-1
    // and decreases
    size_t idx_nondir = 0;
    size_t idx_dir = m_mesh.n_faces()-1;
    for (Face* face : m_mesh.get_faces()){
      if (idx_nondir > idx_dir){
        std::cout << "Error during creation vector to renumber faces: " << idx_nondir << ", " << idx_dir << "\n";
        exit(1);
      }
      if (type(*face) == "dir"){
        new_to_old[idx_dir] = face->global_index();
        idx_dir--;
      }else{
        new_to_old[idx_nondir] = face->global_index();
        idx_nondir++;
      }
    }
    // Check: idx_dir and idx_nondir must just have crossed
    if (idx_nondir != idx_dir + 1){
      std::cout << "Error in creating vector to renumber faces: " << idx_nondir << "/" << m_mesh.n_faces() - m_n_dir_faces << " || " << idx_dir << "/" << m_mesh.n_faces() - m_n_dir_faces << "\n";
      exit(1);
    }

    // Reordering
    m_mesh.renum('F', new_to_old);

 }


