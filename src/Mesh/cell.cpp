// Class to define a cell
//              Members: vertices, edges, neighbouring cells...
//              Methods: index, diameter, area, center of mass...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include "cell.hpp"
#include "face.hpp"
//#include "edge.hpp"
//#include "vertex.hpp"
//#include <cmath>
//#include <math.h>
#include <iostream>

using namespace HArDCore3D;

Cell::Cell(size_t iC, 
	   Mesh *mesh,
	   std::vector<Face *> faces,
	   std::vector<Edge *> edges,
	   std::vector<Vertex *> vertices,
	   std::vector<Vector3d> face_normals,
	   bool boundary,
	   double measure, 
	   double diam, 
	   Vector3d center_mass)
  : _iC(iC),
    _mesh(mesh),
    _faces(faces),
    _edges(edges),
    _vertices(vertices),
    _neighbours(0),
    _face_normals(std::move(face_normals)),
    _boundary(boundary),
    _measure(measure),
    _diam(diam),
    _center_mass(center_mass) {
  // Do nothing
}

Cell::~Cell() {}

std::vector<Face *> Cell::get_faces() const { return _faces; }
std::vector<Edge *> Cell::get_edges() const { return _edges; }
std::vector<Vertex *> Cell::get_vertices() const { return _vertices; }
std::vector<Cell *> Cell::get_neighbours() const { return _neighbours; }

Face *Cell::face(size_t i) const {
  if (i < _faces.size()) {
    return _faces[i];
  } else {
    throw "No face at local index";
  }
}

Edge *Cell::edge(size_t i) const {
  if (i < _edges.size()) {
    return _edges[i];
  } else {
    throw "No edge at local index";
  }
}

Vertex *Cell::vertex(size_t i) const {
  if (i < _vertices.size()) {
    return _vertices[i];
  } else {
    throw "No vertex at local index";
  }
}

Cell *Cell::neighbour(size_t i) const {
  if (i < _neighbours.size()) {
    return _neighbours[i];
  } else {
    throw "No neighbour at local index";
  }
}

size_t Cell::index_face(const Face* F) const {
  size_t i = 0;
  size_t nfac = n_faces();
  while(i < nfac && face(i) != F){
    i++;
  }
  if (i >= nfac || face(i) != F){
    throw "Face does not belong to cell";
  }
  return i;
}

size_t Cell::index_edge(const Edge* E) const {
  size_t i = 0;
  size_t nedg = n_edges();
  while(i < nedg && edge(i) != E){
    i++;
  }
  if (i >= nedg || edge(i) != E){
    throw "Edge does not belong to cell";
  }
  return i;
}

size_t Cell::index_vertex(const Vertex* V) const {
  size_t i = 0;
  size_t nvert = n_vertices();
  while(i < nvert && vertex(i) != V){
    i++;
  }
  if (i >= nvert || vertex(i) != V){
    throw "Vertex does not belong to cell";
  }
  return i;
}

Vector3d Cell::face_normal(size_t i) const {
  if (i < _face_normals.size()) {
    return _face_normals[i];
  } else {
    throw "No face at local index when searching for normal";
  }
}

int Cell::face_orientation(size_t i) const {
  if (i < _face_normals.size()) {
    double dotprod = (_faces[i]->normal()).dot(face_normal(i));
    if (dotprod == 0){
      throw "Problem in face orientation: dotprod is zero";
    }
    return (dotprod > 0) ? 1 : -1;
  } else {
    throw "No face at local index when searching for normal";
  }
}

bool Cell::add_face(Face *face) {
  _faces.push_back(face);
  return true;
}
bool Cell::add_edge(Edge *edge) {
  _edges.push_back(edge);
  return true;
}
bool Cell::add_vertex(Vertex *vertex) {
  _vertices.push_back(vertex);
  return true;
}
bool Cell::add_neighbour(Cell *neigh) {
  _neighbours.push_back(neigh);
  return true;
}
bool Cell::add_normal(Vector3d normal) {
  _face_normals.push_back(normal);
  return true;
}


void Cell::set_boundary(bool val) {
  _boundary = val;
}
void Cell::set_global_index(size_t idx) {
  _iC = idx;
}


