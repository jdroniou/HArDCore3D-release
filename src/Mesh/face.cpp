// Class to define a face
//              Members: cells, vertices...
//              Methods: index, diameter, center of mass...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "face.hpp"
#include "edge.hpp"
//#include "cell.hpp"
//#include "vertex.hpp"
//#include "mesh.hpp"
#include <iostream>

using namespace HArDCore3D;
Face::Face(size_t iF, 
	   Mesh *mesh,
	   std::vector<Edge *> edges,
	   std::vector<Vertex *> vertices,
	   bool boundary,
	   double measure,
	   double diam,
	   Vector3d center_mass,
     Vector3d normal)
  : _iF(iF), 
    _mesh(mesh), 
    _cells(0),
    _edges(edges),
    _vertices(vertices),
    _boundary(boundary),
    _measure(measure),
    _diam(diam),
    _center_mass(center_mass),
    _normal(normal) {
  // Do nothing
}

Face::~Face() {}

Cell *Face::cell(size_t i) const {
  if (i < _cells.size()) {
    return _cells[i];
  } else {
    throw "No cell at edge local index";
  }
}

Edge *Face::edge(size_t i) const {
  if (i < _edges.size()) {
    return _edges[i];
  } else {
    throw "No vertex at edge local index";
  }
}

Vertex *Face::vertex(size_t i) const {
  if (i < _vertices.size()) {
    return _vertices[i];
  } else {
    throw "No vertex at edge local index";
  }
}

size_t Face::index_edge(const Edge* E) const {
  size_t i = 0;
  size_t nedg = n_edges();
  while(i < nedg && edge(i) != E){
    i++;
  }
  if (i >= nedg || edge(i) != E){
    throw "Edge does not belong to face";
  }
  return i;
}

size_t Face::index_vertex(const Vertex* V) const {
  size_t i = 0;
  size_t nvert = n_vertices();
  while(i < nvert && vertex(i) != V){
    i++;
  }
  if (i >= nvert || vertex(i) != V){
    throw "Vertex does not belong to face";
  }
  return i;
}

Vector3d Face::edge_normal(size_t i) const {
  if (i < _edges.size()) {
    return _normal.cross(_edges[i]->tangent());
  } else {
    throw "No vertex at edge local index";
  }
}

int Face::edge_orientation(size_t i) const {
  if (i < _edges.size()) {
    // Here, we assume that the face is star-shaped with respect to its center of mass, 
    // so that the vector from this center of mass to the edge midpoint points outside the face. 
    double dotprod = (_edges[i]->center_mass() - _center_mass).dot(edge_normal(i));
    if (dotprod == 0){
      throw "Problem in face orientation: dotprod is zero";
    }
    return (dotprod > 0) ? 1 : -1;
  } else {
    throw "No vertex at edge local index";
  }
}

std::vector<Cell *> Face::get_cells() const { return _cells; }
std::vector<Edge *> Face::get_edges() const { return _edges; }
std::vector<Vertex *> Face::get_vertices() const { return _vertices; }

Vector3d Face::center_mass() const { return _center_mass; }
Vector3d Face::normal() const { return _normal; }

void Face::add_cell(Cell *cell) {
  _cells.push_back(cell);
}
void Face::add_edge(Edge *edge) {
  _edges.push_back(edge);
}
void Face::add_vertex(Vertex *vertex) {
  _vertices.push_back(vertex);
}


void Face::set_boundary(bool val) {
  _boundary = val;
}
void Face::set_global_index(size_t idx) {
  _iF = idx;
}


