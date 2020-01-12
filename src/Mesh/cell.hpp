// Class to define a cell
//              Members: vertices, edges, neighbouring cells...
//              Methods: index, diameter, area, center of mass...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#ifndef CELL_HPP
#define CELL_HPP
#include <vector>
#include <Eigen/Dense>
namespace HArDCore3D {  // forward declaration
  class Mesh;
  class Vertex;
  class Edge;
  class Face;
}
namespace HArDCore3D {
  using Eigen::Vector3d;

  /*!
   *    @addtogroup Mesh
   * @{
   */

  // ----------------------------------------------------------------------------
  //                            Class definition
  // ----------------------------------------------------------------------------

  /// The Cell class provides description of a cell
  class Cell {
  public:
    /**
     * Constructor. Cells are created last, so all their connectivities are available when constructing them.
     *
     * @param iG global index of the Cell
     * @param mesh pointer to the mesh that the cell is contained in
     * @param faces faces of the cell
     * @param edges edges of the cell
     * @param vertices vertices of the cell
     * @param neighbours neighbouring cells of the cell
     * @param face_normals unit outer normal to the faces, in the same order as the faces
     * @param measure volume of the cell
     * @param diam diameter of the cell
     * @param center_mass center of mass of the cell
     */
    Cell(size_t iC, 
         Mesh *mesh,
         std::vector<Face *> faces,
         std::vector<Edge *> edges,
         std::vector<Vertex *> vertices,
         std::vector<Vector3d> face_normals,
         bool boundary,
         double measure, 
         double diam, 
         Vector3d center_mass);
    ~Cell();

    inline size_t global_index() const;         ///< cell index
    inline size_t n_faces() const;  ///<  returns number of faces of the cell
    inline size_t n_edges() const;  ///<  returns number of edges of the cell
    inline size_t n_vertices() const;  ///< returns number of vertices of the cell

    std::vector<Face *> get_faces() const;      ///< returns the list of faces of the cell
    std::vector<Edge *> get_edges() const;      ///< returns the list of edges of the cell
    std::vector<Vertex *> get_vertices() const;         ///< returns the list of vertices of the cell
    std::vector<Cell *> get_neighbours() const;          ///< returns the list of neighbours of the cell
    Face *face(size_t iL) const;        ///< returns the iL-th face of the cell
    Edge *edge(size_t iL) const;        ///< returns the iL-th edge of the cell
    Vertex *vertex(size_t iL) const;   ///< returns the iL-th edge of the cell
    Cell *neighbour(size_t iL) const; ///< returns the iL-th neighbour of the cell
    size_t index_face(const Face* F) const; ///< reciprocal of face(i): returns the local index of face F in the cell    
    size_t index_edge(const Edge* E) const; ///< reciprocal of edge(i): returns the local index of edge E in the cell
    size_t index_vertex(const Vertex* V) const; ///< reciprocal of vertex(i): returns the local index of vertex V in the cell

    inline bool is_boundary() const; ///< returns true if cell touches the boundary

    inline double measure() const; ///< returns area of cell
    inline double diam() const;  ///< returns diameter of cell
    Vector3d face_normal(size_t i) const;  ///< returns the outer normal to the i-th face
    inline Vector3d center_mass() const;  ///< returns the center of mass of the cell
    int face_orientation(size_t i) const; ///< returns the relative orientation of the i-th face with respect to the cell (that is, +1 if the normal to the face is the outer normal to the cell, -1 otherwise).

    bool add_face(Face *face);  ///< add a face to the cell
    bool add_edge(Edge *edge);  ///< add an edge to the cell
    bool add_vertex(Vertex *vertex);  ///< add a vertex to the cell
    bool add_neighbour(Cell *neigh);  ///< add a cell to the neighbour
    bool add_normal(Vector3d normal);  ///< add an outer normal to the cell

    void set_boundary(bool val); ///< Set the _boundary value of the cell to val
    void set_global_index(size_t idx);  ///< Set the global index of the cell to idx. Used to re-index the cells, should essentially only be used inside Mesh::renum


  private:
    size_t _iC;         ///< cell global index
    Mesh *_mesh;                     ///< pointer to the owner mesh
    std::vector<Face *> _faces;      ///< list of cell faces
    std::vector<Edge *> _edges;      ///< list of cell edges
    std::vector<Vertex *> _vertices;      ///< a list of cell vertices
    std::vector<Cell *> _neighbours;  ///< list of cell neighbours

    std::vector<Vector3d> _face_normals;        ///< list the outer normals to the face (in the order of _faces)

    bool _boundary;                    ///< flag is cell boundary?

    double _measure;              ///< area of the cell
    double _diam;                 ///< diameter of the cell
    Vector3d _center_mass;        ///< center of mass of the cell

  };


  // ----------------------------------------------------------------------------
  //                            Implementations
  // ----------------------------------------------------------------------------

  Vector3d Cell::center_mass() const { return _center_mass; }
  bool Cell::is_boundary() const { return _boundary; }
  double Cell::measure() const { return _measure; }
  double Cell::diam() const { return _diam; }
  size_t Cell::global_index() const { return _iC; }
  size_t Cell::n_faces() const { return _faces.size(); }
  size_t Cell::n_edges() const { return _edges.size(); }
  size_t Cell::n_vertices() const { return _vertices.size(); }

  /*@}*/
}
#endif /* CELL_HPP */
