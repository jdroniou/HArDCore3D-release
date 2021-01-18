// Class to define a face
//              Members: cells, vertices...
//              Methods: index, diameter, center of mass...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#ifndef FACE_HPP
#define FACE_HPP
#include <vector>
#include <Eigen/Dense>
namespace HArDCore3D {  // forward declaration
  class Mesh;
  class Vertex;
  class Edge;
  class Cell;
}

namespace HArDCore3D {
  using Eigen::Vector3d;
  /*!
   * @addtogroup Mesh
   * @{
   */

  // ----------------------------------------------------------------------------
  //                            Class definition
  // ----------------------------------------------------------------------------

  /// The Face class provides description of an edge
  class Face {
    /**
     * A class representing a face of a 2D mesh. Contains a pointer to the mesh, as well as to 
     * the vertices, edges and cells connected to that face
     */
  public:
    /**
     * Constructor. Faces are created after vertices and edges, but before cells. So all their connectivites
     *               except cells can be filled in
     *
     * @param iF global face number
     * @param mesh pointer to the mesh
     * @param edges edges of the face
     * @param vertices vertices of the face
     * @param measure area of the face
     * @param diam diam of the face
     * @param center_mass center of mass of the face
     */
    Face(size_t iF, 
   Mesh *mesh,
   std::vector<Edge *> edges,
   std::vector<Vertex *> vertices,
   bool boundary,
   double measure,
   double diam,
   Vector3d center_mass,
   Vector3d normal);
    ~Face(); //  destructor, nothing special

    inline size_t global_index() const;  ///< returns the face global index
    inline size_t n_cells() const;                      ///< returns the number of cells neighbouring the face
    inline size_t n_edges() const;                      ///< returns the number of edges neighbouring the face
    inline size_t n_vertices() const;           ///< returns the number of vertices neighbouring the face

    std::vector<Cell *> get_cells() const; ///< list of cells that are neighbours of the face
    std::vector<Edge *> get_edges() const; ///< list of edges of the face
    std::vector<Vertex *> get_vertices() const; ///< list of vertices of the face
    Cell *cell(size_t i) const;  ///< returns pointer to the i-th cell neighbour of the face
    Edge *edge(size_t i) const;  ///< returns a pointer to the i-th edge of the face
    Vertex *vertex(size_t i) const;  ///< returns a pointer to the i-th vertex of the face
    size_t index_edge(const Edge* E) const; ///< reciprocal of edge(i): returns the local index of edge E in the face
    size_t index_vertex(const Vertex* V) const; ///< reciprocal of vertex(i): returns the local index of vertex V in the face
    Vector3d edge_normal(size_t i) const; ///< returns the normal to the i-th edge in the plane spanned by the face. This is not necessarily the outer normal to the face. To obtain this outer normal, multiply by edge_orientation(i).
    int edge_orientation(size_t i) const; ///< returns the orientation of the i-th edge. It is the number such that edge_orientation * edge_normal is the outer normal to the face.

    inline double measure() const;   ///< measure of the face
    inline double diam() const;  ///< returns diameter of face
    Vector3d center_mass() const;  ///< get the center of mass of the face
    Vector3d normal() const;  ///< get the normal to the face
    inline bool is_boundary() const;  ///< getter to see if the face is a boundary face

    void add_cell(Cell *cell);      ///< Add a cell to the face
    void add_edge(Edge *edge);      ///< Add a edge to the face
    void add_vertex(Vertex *vertex);      ///< Add a vertex to the face

    void set_boundary(bool val); ///< Set the _boundary value of the face to val
    void set_global_index(size_t idx);  ///< Set the global index of the face to idx. Used to re-index the faces, should essentially only be used inside Mesh::renum

  private:
    size_t _iF;
    std::vector<Cell *> _cells;
    std::vector<Edge *> _edges;
    std::vector<Vertex *> _vertices;
    bool _boundary;
    double _measure;
    double _diam;
    Vector3d _center_mass;
    Vector3d _normal;
        
  };


  // ----------------------------------------------------------------------------
  //                            Implementations
  // ----------------------------------------------------------------------------

  size_t Face::global_index() const { return _iF; }
  size_t Face::n_cells() const { return _cells.size();}
  size_t Face::n_edges() const { return _edges.size();}
  size_t Face::n_vertices() const { return _vertices.size();}
  bool Face::is_boundary() const { return _boundary; }
  inline double Face::measure() const { return _measure; }
  inline double Face::diam() const { return _diam; }


  /*@}*/
}
#endif /* FACE_HPP */
