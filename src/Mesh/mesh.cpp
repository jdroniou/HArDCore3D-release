// Class to describe a mesh.
//              Members: cells, faces, edges, vertices...
//              Methods: h_max...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include "mesh.hpp"
#include "cell.hpp"
#include "face.hpp"
#include "edge.hpp"
#include "vertex.hpp"
#include <iostream>
#include <set>
#include <Eigen/Dense>  //Vector3d

using namespace HArDCore3D;
using Eigen::Vector3d;
Mesh::Mesh() : _mesh_name("mesh-3d"),
	       _cells(0),
	       _faces(0),
	       _edges(0),
	       _vertices(0),
	       _b_cells(0),
	       _b_faces(0),
	       _b_edges(0),
	       _b_vertices(0),
	       _i_cells(0),
	       _i_faces(0),
	       _i_edges(0),
	       _i_vertices(0),
	       _h_max(0)
{
  // do nothing
}

Mesh::~Mesh(){
  for (auto& cell : _cells){
    delete cell;
  }
  for (auto& face : _faces){
    delete face;
  }
  for (auto& edge : _edges){
    delete edge;
  }
  for (auto& vertex : _vertices){
    delete vertex;
  }
}

Cell* Mesh::cell(size_t iC) const {
  if (iC < n_cells()) {
    return _cells[iC];
  } else {
    throw "Trying to access cell at global index which does not exist";
  }
}

Face* Mesh::face(size_t iF) const {
  if (iF < n_faces()) {
    return _faces[iF];
  } else {
    throw "Trying to access cell at global index which does not exist";
  }
}

Edge* Mesh::edge(size_t iE) const {
  if (iE < n_edges()) {
    return _edges[iE];
  } else {
    throw "Trying to access edge at global index which does not exist";
  }
}

Vertex* Mesh::vertex(size_t iV) const {
  if (iV < n_vertices()) {
    return _vertices[iV];
  } else {
    throw "Trying to access vertex at global index which does not exist";
  }
}

Cell* Mesh::b_cell(size_t iC) const {
  if (iC < n_b_cells()) {
    return _b_cells[iC];
  } else {
    throw "Trying to access boundary cell at global index which does not exist";
  }
}

Face* Mesh::b_face(size_t iF) const {
  if (iF < n_b_faces()) {
    return _b_faces[iF];
  } else {
    throw "Trying to access boundary cell at global index which does not exist";
  }
}

Edge* Mesh::b_edge(size_t iE) const {
  if (iE < n_b_edges()) {
    return _b_edges[iE];
  } else {
    throw "Trying to access boundary edge at global index which does not exist";
  }
}

Vertex* Mesh::b_vertex(size_t iV) const {
  if (iV < n_b_vertices()) {
    return _b_vertices[iV];
  } else {
    throw "Trying to access boundary vertex at global index which does not exist";
  }
}

Cell* Mesh::i_cell(size_t iC) const {
  if (iC < n_i_cells()) {
    return _i_cells[iC];
  } else {
    throw "Trying to access interior cell at global index which does not exist";
  }
}

Face* Mesh::i_face(size_t iF) const {
  if (iF < n_i_faces()) {
    return _i_faces[iF];
  } else {
    throw "Trying to access interior cell at global index which does not exist";
  }
}

Edge* Mesh::i_edge(size_t iE) const {
  if (iE < n_i_edges()) {
    return _i_edges[iE];
  } else {
    throw "Trying to access interior edge at global index which does not exist";
  }
}

Vertex* Mesh::i_vertex(size_t iV) const {
  if (iV < n_i_vertices()) {
    return _i_vertices[iV];
  } else {
    throw "Trying to access interior vertex at global index which does not exist";
  }
}


size_t Mesh::n_b_cells() const {
  return _b_cells.size();
}
size_t Mesh::n_b_faces() const {
  return _b_faces.size();
}
size_t Mesh::n_b_edges() const {
  return _b_edges.size();
}
size_t Mesh::n_b_vertices() const {
  return _b_vertices.size();
}

size_t Mesh::n_i_cells() const {
  return _i_cells.size();
}
size_t Mesh::n_i_faces() const {
  return _i_faces.size();
}
size_t Mesh::n_i_edges() const {
  return _i_edges.size();
}
size_t Mesh::n_i_vertices() const {
  return _i_vertices.size();
}

double Mesh::regularity(){
  /// Regularity factor = maximum of
  ///                     * diameter of cell / (measure of cell)^{1/3}
  ///                     * diameter of cell / diameter of face  [for each face of the cell]
  ///                     * diameter of face / (measure of face)^{1/2}
  //[NO]                  * diameter of cell / diameter of edge  [for each edge of the cell]
  //[NO]                  * diameter of face / diameter of edge  [for each edge of the face]


  double value = 0.0;
  for (auto& icell : get_cells()){
    double hC = icell->diam();

    value = std::max(value, hC / pow(icell->measure(), 1/this->dim()));

    for (auto& iface : icell->get_faces()){
      double hF = iface->diam();

      value = std::max(value, hC / hF);
      value = std::max(value, hF / pow(iface->measure(), 1/(this->dim() - 1)));

      //                              for (auto& iedge : iface->get_edges()){
      //                                      double hE = iedge->diam();

      //                                      value = std::max(value, hC / hE);
      //                                      value = std::max(value, hF / hE);
      //                              }
    }
  }

  return value;

}

void Mesh::renum(const char B, const std::vector<size_t> new_to_old){
        
  switch (B) {
  case 'C': {
    std::vector<Cell*> old_index = _cells;
    for (size_t i=0; i < n_cells(); i++){
      old_index[new_to_old[i]]->set_global_index(i);
      _cells[i] = old_index[new_to_old[i]];
    }
    break;
  }

  case 'F': {
    std::vector<Face*> old_index = _faces;
    for (size_t i=0; i < n_faces(); i++){
      old_index[new_to_old[i]]->set_global_index(i);
      _faces[i] = old_index[new_to_old[i]];
    }
    break;
  }

  case 'E': {
    std::vector<Edge*> old_index = _edges;
    for (size_t i=0; i < n_edges(); i++){
      old_index[new_to_old[i]]->set_global_index(i);
      _edges[i] = old_index[new_to_old[i]];
    }
    break;
  }

  case 'V': {
    std::vector<Vertex*> old_index = _vertices;
    for (size_t i=0; i < n_vertices(); i++){
      old_index[new_to_old[i]]->set_global_index(i);
      _vertices[i] = old_index[new_to_old[i]];
    }
    break;
  }

  }

}

