// Class to define an edge
//              Members: cells, vertices...
//              Methods: index, measure, center of mass...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#ifndef EDGE_HPP
#define EDGE_HPP
#include <vector>
#include <Eigen/Dense>
namespace HArDCore3D {  // forward declaration
  class Mesh;
  class Vertex;
  class Face;
  class Cell;
}

namespace HArDCore3D {
  using Eigen::Vector3d;
  /*!
   *       @addtogroup Mesh
   * @{
   */

  // ----------------------------------------------------------------------------
  //                            Class definition
  // ----------------------------------------------------------------------------

  /// The Edge class provides description of an edge
  class Edge {
    /**
     * A class representing an edge of a 3D mesh. Contains a
     * pointer to the mesh, as well as to the cells, faces and vertices connected to that edge
     */
  public:
    /**
     * Constructor. Edges are constructed after vertices, so these can be feed into the constructor.
     *
     * @param iE global edge number
     * @param mesh pointer to the mesh
     * @param vertices vertices of the edge
     * @param measure length of the edge
     * @param center_mass midpoint of the edge
     */
    Edge(size_t iE, 
	 Mesh *mesh, 
	 std::vector<Vertex *> vertices,
	 bool boundary,
	 double measure, 
	 Vector3d center_mass);
    ~Edge(); //  destructor, nothing special

    inline size_t global_index() const;  ///< returns the edge global index
    inline size_t n_cells() const;                      ///< returns the number of cells neighbouring the edge
    inline size_t n_faces() const;                      ///< returns the number of faces neighbouring the edge

    std::vector<Cell *> get_cells() const; ///< list of cells that are neighbours of the edge
    std::vector<Face *> get_faces() const; ///< list of faces of the edge
    std::vector<Vertex *> get_vertices() const; ///< list of vertices of the edge
    Cell *cell(size_t i) const;  ///< returns pointer to the i-th cell neighbour of the edge
    Face *face(size_t i) const;  ///< returns pointer to the i-th face neighbour of the edge
    Vertex *vertex(size_t i) const;  ///< returns a pointer to the i-th vertex of the edge
    size_t index_vertex(Vertex* V) const; ///< reciprocal of vertex(i): returns the local index of vertex V in the edge

    inline double measure() const;   ///< length of the edge
    inline double diam() const;   ///< length of the edge
    inline Vector3d center_mass() const;  ///< get the midpoint of the edge
    inline Vector3d tangent() const    ///< get the normalised tangent to the edge, oriented from the first vertex to the second vertex.
    {
      return _line.normalized();
    }
    inline bool is_boundary() const;  ///< getter to see if edge is boundary edge

    void add_cell(Cell *cell);      ///< Add a new cell to the edge
    void add_face(Face *face);      ///< Add a new face to the edge
    void add_vertex(Vertex *vertex);      ///< Add a new vertex to the edge

    void set_boundary(bool val); ///< Set the _boundary value of the face to val
    void set_global_index(size_t idx);  ///< Set the global index of the edge to idx. Used to re-index the edges, should essentially only be used inside Mesh::renum

  private:
    size_t _iE;
    Mesh *_mesh;
    bool _boundary;
    std::vector<Cell *> _cells;
    std::vector<Face *> _faces;
    std::vector<Vertex *> _vertices;
    double _measure;
    Vector3d _center_mass;
    Vector3d _line;
  };


  // ----------------------------------------------------------------------------
  //                            Implementations
  // ----------------------------------------------------------------------------

  size_t Edge::global_index() const { return _iE; }
  size_t Edge::n_cells() const { return _cells.size();}
  size_t Edge::n_faces() const { return _faces.size();}
  bool Edge::is_boundary() const { return _boundary; }

  inline double Edge::measure() const { return _measure; }
  inline double Edge::diam() const { return _measure; }
  inline Vector3d Edge::center_mass() const { return _center_mass; }


  /*@}*/
}
#endif /* EDGE_HPP */
