// Class to define an edge
//              Members: cells, vertices...
//              Methods: index, measure, center of mass...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "edge.hpp"
//#include "cell.hpp"
#include "vertex.hpp"
//#include "mesh.hpp"
#include <iostream>

using namespace HArDCore3D;
Edge::Edge(size_t iE, 
	   Mesh *mesh, 
	   std::vector<Vertex *> vertices,
	   bool boundary,
	   double measure, 
	   Vector3d center_mass)
  : _iE(iE), 
    _boundary(boundary),
    _cells(0),
    _faces(0),
    _vertices(vertices),
    _measure(measure),
    _center_mass(center_mass)
{
  _line = vertices[1]->coords() - vertices[0]->coords();
}

Edge::~Edge() {}

Cell *Edge::cell(size_t i) const {
  if (i < _cells.size()) {
    return _cells[i];
  } else {
    throw "No cell at edge local index";
  }
}

Face *Edge::face(size_t i) const {
  if (i < _faces.size()) {
    return _faces[i];
  } else {
    throw "No face at edge local index";
  }
}

Vertex *Edge::vertex(size_t i) const {
  if (i < _vertices.size()) {
    return _vertices[i];
  } else {
    throw "No vertex at edge local index";
  }
}

size_t Edge::index_vertex(const Vertex* V) const {
  size_t i = 0;
  while(i < 2 && vertex(i) != V){
    i++;
  }
  if (i >= 2 || vertex(i) != V){
    throw "Vertex does not belong to edge";
  }
  return i;
}

std::vector<Cell *> Edge::get_cells() const { return _cells; }
std::vector<Face *> Edge::get_faces() const { return _faces; }
std::vector<Vertex *> Edge::get_vertices() const { return _vertices; }

void Edge::add_cell(Cell *cell) {
  _cells.push_back(cell);
}
void Edge::add_face(Face *face) {
  _faces.push_back(face);
}
void Edge::add_vertex(Vertex *vertex) {
  _vertices.push_back(vertex);
}


void Edge::set_boundary(bool val) {
  _boundary = val;
}
void Edge::set_global_index(size_t idx) {
  _iE = idx;
}


