// Class to build the mesh data after having read the mesh file
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "mesh_builder.hpp"
#include "mesh.hpp"
#include "cell.hpp"
#include "face.hpp"
#include "edge.hpp"
#include "vertex.hpp"

#include "StemMesh/mesh3D_greader.hh"
#include "StemMesh/mesh3D.hh"
#include <iostream>
#include <deque>

using namespace HArDCore3D;
MeshBuilder::MeshBuilder(std::string mesh_file, std::string mesh_type)
  : _mesh_file(mesh_file),
    _mesh_type(mesh_type) 
{}

std::unique_ptr<Mesh> MeshBuilder::build_the_mesh() {

  // Create the stem of the mesh
  StemMesh3D::mesh_3Dv stem_mesh;  
  StemMesh3D::ExtFileInput input_mesh;
  int offset = 1;
  if (_mesh_type=="TG") {
    input_mesh.TETGEN_format(stem_mesh, _mesh_file, offset);
  }
  else if (_mesh_type=="MSH") {
    input_mesh.MSH_format(stem_mesh, _mesh_file, offset);
  }
  else if (_mesh_type=="RF") {
    offset = 0;
    input_mesh.REGN_FACE_format(stem_mesh, _mesh_file, offset);
  }
  else {
    std::cout << "Unknown mesh type: " << _mesh_type << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Create the mesh, empty first
  std::unique_ptr<Mesh> mesh = std::make_unique<Mesh>();

  // Create vertices
  for (size_t iV = 0; iV < stem_mesh.n_vertex(); iV++){
    Vector3d coords = Vector3d(stem_mesh.coords_V(iV, 0), stem_mesh.coords_V(iV, 1), stem_mesh.coords_V(iV, 2));
    bool boundary = stem_mesh.is_boundary_vrtx(iV);
    Vertex* vertex = new Vertex(iV, mesh.get(), coords, boundary);
    mesh->add_vertex(vertex);
    if (boundary) {
      mesh->add_b_vertex(vertex);
    } else {
      mesh->add_i_vertex(vertex);
    }
  }
  // list connected vertices
  for (auto& v : mesh->get_vertices()){
    size_t iV = v->global_index();
    std::vector<size_t> idx_v(0);
    stem_mesh.get_vrtx_vrtx(iV, idx_v);
    for (size_t ilV = 0; ilV < idx_v.size(); ilV++){
      Vertex* vert = mesh->vertex(idx_v[ilV]);
      v->add_vertex(vert);
    }
  }

  // Create edges
  for (size_t iE = 0; iE < stem_mesh.n_edge(); iE++){
    // parameters for constructor
    std::vector<Vertex *> vertices(0);
    std::vector<size_t> idx_v(0);
    stem_mesh.get_edge_vrtx(iE, idx_v);
    for (size_t ilV = 0; ilV < idx_v.size(); ilV++){
      Vertex* vertex = mesh->vertex(idx_v[ilV]);
      vertices.push_back(vertex);
    }
    bool boundary = stem_mesh.is_boundary_edge(iE);
    double measure = stem_mesh.get_edge_measure(iE);
    Vector3d center_mass = Vector3d(stem_mesh.coords_E(iE, 0), stem_mesh.coords_E(iE, 1), stem_mesh.coords_E(iE, 2));

    // Construct new edge
    Edge* edge = new Edge(iE, mesh.get(), vertices, boundary, measure, center_mass);
    mesh->add_edge(edge);
    if (boundary) {
      mesh->add_b_edge(edge);
    } else {
      mesh->add_i_edge(edge);
    }

    // Add the edge to its vertices
    for (auto& v : vertices){
      v->add_edge(edge);
    }

  }

  // Create faces
  for (size_t iF = 0; iF < stem_mesh.n_face(); iF++){
    // parameters for constructor
    std::vector<Edge *> edges(0);
    std::vector<size_t> idx_e(0);
    stem_mesh.get_face_edge(iF, idx_e);
    for (size_t ilE = 0; ilE < idx_e.size(); ilE++){
      Edge* edge = mesh->edge(idx_e[ilE]);
      edges.push_back(edge);
    }
    std::vector<Vertex *> vertices(0);
    std::vector<size_t> idx_v(0);
    stem_mesh.get_face_vrtx(iF, idx_v);
    for (size_t ilV = 0; ilV < idx_v.size(); ilV++){
      Vertex* vertex = mesh->vertex(idx_v[ilV]);
      vertices.push_back(vertex);
    }
    bool boundary = stem_mesh.is_boundary_face(iF);
    double measure = stem_mesh.get_face_measure(iF);
    double diam = stem_mesh.get_face_diam(iF);
    Vector3d center_mass = Vector3d(stem_mesh.coords_F(iF, 0), stem_mesh.coords_F(iF, 1), stem_mesh.coords_F(iF, 2));
    Vector3d normal = Vector3d(stem_mesh.get_nor(iF, 0), stem_mesh.get_nor(iF, 1), stem_mesh.get_nor(iF, 2));


    // Construct new face
    Face* face = new Face(iF, mesh.get(), edges, vertices, boundary, measure, diam, center_mass, normal);
    mesh->add_face(face);
    if (boundary) {
      mesh->add_b_face(face);
    } else {
      mesh->add_i_face(face);
    }

    // Add the face to its vertices, and edges
    for (auto& v : vertices){
      v->add_face(face);
    }
    for (auto& e : edges){
      e->add_face(face);
    }
  }

  // Create cells
  for (size_t iC = 0; iC < stem_mesh.n_region(); iC++){
    // parameters for constructor
    std::vector<Face *> faces(0);
    std::vector<Vector3d> face_normals(0);
    std::vector<size_t> idx_f(0);
    stem_mesh.get_regn_face(iC, idx_f);
    for (size_t ilF = 0; ilF < idx_f.size(); ilF++){
      Face* face = mesh->face(idx_f[ilF]);
      faces.push_back(face);

      // normal to the face, and then we check orientation
      size_t iF = face->global_index();
      Vector3d nor_face = Vector3d(stem_mesh.get_nor(iF, 0), stem_mesh.get_nor(iF, 1), stem_mesh.get_nor(iF, 2));
      if (!stem_mesh.ok_regn_face(iC, ilF)){
        nor_face = -nor_face;
      }
      face_normals.push_back(nor_face);
    }
    std::vector<Edge *> edges(0);
    std::vector<size_t> idx_e(0);
    stem_mesh.get_regn_edge(iC, idx_e);
    for (size_t ilE = 0; ilE < idx_e.size(); ilE++){
      Edge* edge = mesh->edge(idx_e[ilE]);
      edges.push_back(edge);
    }
    std::vector<Vertex *> vertices(0);
    std::vector<size_t> idx_v(0);
    stem_mesh.get_regn_vrtx(iC, idx_v);
    for (size_t ilV = 0; ilV < idx_v.size(); ilV++){
      Vertex* vertex = mesh->vertex(idx_v[ilV]);
      vertices.push_back(vertex);
    }

    bool boundary = stem_mesh.is_boundary_regn(iC);
    double measure = stem_mesh.get_regn_measure(iC);
    double diam = stem_mesh.get_regn_diam(iC);
    Vector3d center_mass = Vector3d(stem_mesh.coords_R(iC, 0), stem_mesh.coords_R(iC, 1), stem_mesh.coords_R(iC, 2));

    // Construct new cell
    Cell* cell = new Cell(iC, mesh.get(), faces, edges, vertices, face_normals, boundary, measure, diam, center_mass);
    mesh->add_cell(cell);
    if (boundary) {
      mesh->add_b_cell(cell);
    } else {
      mesh->add_i_cell(cell);
    }

    // Add the cell to its vertices, edges and faces
    for (auto& v : vertices){
      v->add_cell(cell);
    }
    for (auto& e : edges){
      e->add_cell(cell);
    }
    for (auto& f : faces){
      f->add_cell(cell);
    }

  }

  // Create cell neighbours
  for (auto& c : mesh->get_cells()){
    size_t iC = c->global_index();
    std::vector<size_t> idx_n(0);
    stem_mesh.get_regn_regn(iC, idx_n);
    for (size_t ilN = 0; ilN < idx_n.size(); ilN++){
      Cell* neigh = mesh->cell(idx_n[ilN]);
      c->add_neighbour(neigh);
    }
  }

  // Mesh size
  mesh->set_h_max(stem_mesh.h_max()); 

  return mesh;
}



