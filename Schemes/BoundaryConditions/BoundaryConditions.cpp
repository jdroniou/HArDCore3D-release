// Class to provide description of boundary conditions
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "BoundaryConditions.hpp"

using namespace HArDCore3D;

// ----------------------------------------------------------------------------
//                          Implementation Plane and Hull
// ----------------------------------------------------------------------------
// Tolerances for tests
constexpr double tol=1e-8;

bool Plane::is_in(const VectorRd & x) const {
  return ( std::abs( (x - point).dot(normal) ) < tol );
};

bool Plane::is_on_side(const VectorRd & x, const double & side) const {
  return ( side * (x - point).dot(normal) > -tol );

};

bool Hull::is_in(const VectorRd & x) const {
  // Cut pyramid made by the hull and the apex into simplices and loop over the simplices to see if V is inside
  for (size_t i=0; i < nodes.size() ; i++){
      // Two consecutive nodes
      VectorRd v1 = nodes[i];
      VectorRd v2 = nodes[(i+1) % nodes.size()];
      // Find barycentric coordinates of x with respect to simplex formed by the apex, the center, and the two nodes
      Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
      A.col(0) = center-apex;
      A.col(1) = v1 - apex;
      A.col(2) = v2 - apex;
      VectorRd lambda = A.partialPivLu().solve(x-apex);
      
      // Check if x is in the triangle center-v1-v2
      if (lambda(0)>-tol && lambda(1)>-tol && lambda(2)>-tol && std::abs(lambda.sum() - 1.)<tol) return true;
    }
  return false;

};


// ----------------------------------------------------------------------------
//                          Implementation BoundaryConditions
// ----------------------------------------------------------------------------

// Class
BoundaryConditions::BoundaryConditions(const std::string bc_id, Mesh& mesh)
  : m_bc_id(bc_id),
    m_mesh(mesh),
    m_n_dir_faces(0),
    m_n_dir_edges(0),
    m_n_dir_vertices(0)
  {
    // Compute number of Dirichlet faces
    for (Face* face : m_mesh.get_b_faces()){
      if (type(*face) == "dir"){
        m_n_dir_faces++;
      }
    }
   // Compute number of Dirichlet edges
    for (Edge* edge : m_mesh.get_b_edges()){
      if (type(*edge) == "dir"){
        m_n_dir_edges++;
      }
    }
    // Compute number of Dirichlet vertices
    for (Vertex* vertex : m_mesh.get_b_vertices()){
      if (type(*vertex) == "dir"){
        m_n_dir_vertices++;
      }
    }
  }


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
        return ( plane_x_zero.is_in(face) ? "dir" : "neu");
      }else if (m_bc_id == "M1"){
        // Dirichlet at x=0 a small rectangle, Neumann everywhere else.
        return ( lower_bottom_corner_x_zero.is_in(face) ? "dir" : "neu");
      }else if (m_bc_id == "M2"){
        // Dirichlet between planes x=.25 and x=.5, Neumann everywhere else.
        return ( plane_x_quarter.is_on_side(face, 1.) && plane_x_half.is_on_side(face, -1.) ? "dir" : "neu");
      }

    }
    return "int";
  }

// Returns the type of BC of an edge
const std::string BoundaryConditions::type(const Edge& edge) const
  {
    if (edge.is_boundary()){
      if (m_bc_id == "D"){
        return "dir";
      }else if (m_bc_id == "N"){
        return "neu";
      }else if (m_bc_id == "M0"){
        // Dirichlet at x=0, Neumann everywhere else.
        return ( plane_x_zero.is_in(edge) ? "dir" : "neu");
      }else if (m_bc_id == "M1"){
        // Dirichlet at x=0 a small rectangle, Neumann everywhere else.
        return ( lower_bottom_corner_x_zero.is_in(edge) ? "dir" : "neu");
      }else if (m_bc_id == "M2"){
        // Dirichlet between planes x=.25 and x=.5, Neumann everywhere else.
        return ( plane_x_quarter.is_on_side(edge, 1.) && plane_x_half.is_on_side(edge, -1.) ? "dir" : "neu");
      }
    }
    return "int";
  }

// Returns the type of BC of a vertex
const std::string BoundaryConditions::type(const Vertex& vertex) const
  {
    if (vertex.is_boundary()){
      if (m_bc_id == "D"){
        return "dir";
      }else if (m_bc_id == "N"){
        return "neu";
      }else if (m_bc_id == "M0"){
        // Dirichlet at x=0, Neumann everywhere else.
        return ( plane_x_zero.is_in(vertex) ? "dir" : "neu");
      }else if (m_bc_id == "M1"){
        // Dirichlet at x=0 a small rectangle, Neumann everywhere else.
        return ( lower_bottom_corner_x_zero.is_in(vertex) ? "dir" : "neu");
      }else if (m_bc_id == "M2"){
        // Dirichlet between planes x=.25 and x=.5, Neumann everywhere else.
        return ( plane_x_quarter.is_on_side(vertex, 1.) && plane_x_half.is_on_side(vertex, -1.) ? "dir" : "neu");
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

// Reorder edges
void BoundaryConditions::reorder_edges(const std::string pos)
  {
    // Create a vector of Dirichlet boundary edges, and a vector of other edges
    std::vector<size_t> dir_edges(m_mesh.n_edges(), 0);
    std::vector<size_t> nondir_edges(m_mesh.n_edges(), 0);
    size_t dir_idx = 0;
    size_t nondir_idx = 0;
    for (Edge* edge : m_mesh.get_edges()){
      if (type(*edge) == "dir"){
        dir_edges[dir_idx] = edge->global_index();
        dir_idx++;
      }else{
        nondir_edges[nondir_idx] = edge->global_index();
        nondir_idx++;
      }
    }
    // check
    if (dir_idx + nondir_idx != m_mesh.n_edges()){
     std::cout << "Error during renumbering edges: " << dir_idx << ", " << nondir_idx << ", " << m_mesh.n_edges() << "\n";
     exit(1);
    }
    
    // Depending on "pos" we put the Dirichlet edges at the end or the start
    std::vector<size_t> new_to_old(m_mesh.n_edges(), 0);
    if (pos=="end"){
      for (size_t i=0; i < nondir_idx; i++){
        new_to_old[i] = nondir_edges[i];
      }     
      for (size_t i=nondir_idx; i < m_mesh.n_edges(); i++){
        new_to_old[i] = dir_edges[i-nondir_idx];
      }
    }else{
      for (size_t i=0; i < dir_idx; i++){
        new_to_old[i] = dir_edges[i];
      }     
      for (size_t i=dir_idx; i < m_mesh.n_edges(); i++){
        new_to_old[i] = nondir_edges[i-dir_idx];
      }
    }

    // Reordering
    m_mesh.renum('E', new_to_old);

 }

// Reorder vertices
void BoundaryConditions::reorder_vertices(const std::string pos)
  {
    // Create a vector of Dirichlet boundary vertices, and a vector of other vertices
    std::vector<size_t> dir_vertices(m_mesh.n_vertices(), 0);
    std::vector<size_t> nondir_vertices(m_mesh.n_vertices(), 0);
    size_t dir_idx = 0;
    size_t nondir_idx = 0;
    for (Vertex* vertex : m_mesh.get_vertices()){
      if (type(*vertex) == "dir"){
        dir_vertices[dir_idx] = vertex->global_index();
        dir_idx++;
      }else{
        nondir_vertices[nondir_idx] = vertex->global_index();
        nondir_idx++;
      }
    }
    // check
    if (dir_idx + nondir_idx != m_mesh.n_vertices()){
     std::cout << "Error during renumbering vertices: " << dir_idx << ", " << nondir_idx << ", " << m_mesh.n_vertices() << "\n";
     exit(1);
    }
    
    // Depending on "pos" we put the Dirichlet vertices at the end or the start
    std::vector<size_t> new_to_old(m_mesh.n_vertices(), 0);
    if (pos=="end"){
      for (size_t i=0; i < nondir_idx; i++){
        new_to_old[i] = nondir_vertices[i];
      }     
      for (size_t i=nondir_idx; i < m_mesh.n_vertices(); i++){
        new_to_old[i] = dir_vertices[i-nondir_idx];
      }
    }else{
      for (size_t i=0; i < dir_idx; i++){
        new_to_old[i] = dir_vertices[i];
      }     
      for (size_t i=dir_idx; i < m_mesh.n_vertices(); i++){
        new_to_old[i] = nondir_vertices[i-dir_idx];
      }
    }

    // Reordering
    m_mesh.renum('V', new_to_old);

 }





