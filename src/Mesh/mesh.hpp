// Class to describe a mesh.
//              Members: cells, faces, edges, vertices...
//              Methods: h_max, get cells, faces, and edges...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
 *
 *       This library was developed around HHO methods, although some parts of it have a more
 * general purpose. If you use this code or part of it in a scientific publication, 
 * please mention the following book as a reference for the underlying principles
 * of HHO schemes:
 *
 * The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications. 
 *  D. A. Di Pietro and J. Droniou. Modeling, Simulation and Applications, vol. 19. 
 *  Springer International Publishing, 2020, xxxi + 525p. doi: 10.1007/978-3-030-37203-3. 
 *  url: https://hal.archives-ouvertes.fr/hal-02151813.
 *
 */

#ifndef MESH_HPP
#define MESH_HPP

#include <cstddef>
#include <string>
#include <vector>
#include <array>
#include <Eigen/Dense>
namespace HArDCore3D { // Forward declaration
  class Vertex;
  class Edge;
  class Face;
  class Cell;
}

/*!     
 * @defgroup Mesh 
 * @brief Classes to describe and access the elements of a 3D mesh
 */

namespace HArDCore3D {

  /*!
   * @addtogroup Mesh
   * @{
   */

  /**
   * Class which represents a 3D mesh. Contains cells, faces, edges, vertices.
   */

  using Eigen::Vector3d;

  // ----------------------------------------------------------------------------
  //                            Class definition
  // ----------------------------------------------------------------------------

  /// The Mesh class provides description of a mesh

  class Mesh {
  public:
    /**
     * default constructor for an empty mesh
     */
    Mesh();
    ~Mesh();

    inline void set_name(std::string name);  ///< set the name of the mesh
    inline std::string get_name();  ///< getter for the edge name
    inline size_t n_cells() const;     ///< number of cells in the mesh
    inline size_t n_faces() const;     ///< number of faces in the mesh
    inline size_t n_edges() const;     ///< number of edges in the mesh
    inline size_t n_vertices() const;  ///< number of vertices in the mesh
    inline double h_max() const;                     ///< max of diameter of cells
    inline size_t dim() const;         ///< dimension of the mesh (3)
    size_t n_b_cells() const;           ///< number of boundary cells
    size_t n_b_faces() const;           ///< number of boundary faces
    size_t n_b_edges() const;           ///< number of boundary edges
    size_t n_b_vertices() const;        ///< number of boundary vertices
    size_t n_i_cells() const;           ///< number of interior cells
    size_t n_i_faces() const;           ///< number of interior faces
    size_t n_i_edges() const;           ///< number of interior edges
    size_t n_i_vertices() const;        ///< number of interior vertices

    inline std::vector<Cell*> get_cells() const;  ///< lists the cells in the mesh.
    inline std::vector<Face*> get_faces() const;  ///< lists the faces in the mesh.
    inline std::vector<Edge*> get_edges() const;  ///< lists the edges in the mesh.
    inline std::vector<Vertex*> get_vertices() const;  ///< lists the vertices in the mesh.
    Cell* cell(size_t iC) const;  ///< get a constant pointer to a cell using its global index
    Face* face(size_t iF) const;   ///< get a constant pointer to a face using its global index
    Edge* edge(size_t iE) const;   ///< get a constant pointer to an edge using its global index
    Vertex* vertex(size_t iV) const;   ///< get a constant pointer to a vertex using its global index

    inline std::vector<Cell*> get_b_cells() const;  ///< lists the boundary cells in the mesh.
    inline std::vector<Face*> get_b_faces() const;  ///< lists the boundary faces in the mesh.
    inline std::vector<Edge*> get_b_edges() const;  ///< lists the boundary edges in the mesh.
    inline std::vector<Vertex*> get_b_vertices() const;  ///< lists the boundary vertices in the mesh.
    Cell* b_cell(size_t iC) const;  ///< get a constant pointer to the iC-th boundary cell
    Face* b_face(size_t iF) const;   ///< get a constant pointer to the iF-th boundary face
    Edge* b_edge(size_t iE) const;   ///< get a constant pointer to the iE-th boundary edge
    Vertex* b_vertex(size_t iV) const;   ///< get a constant pointer to the iV-th boundary vertex

    inline std::vector<Cell*> get_i_cells() const;  ///< lists the interior cells in the mesh.
    inline std::vector<Face*> get_i_faces() const;  ///< lists the interior faces in the mesh.
    inline std::vector<Edge*> get_i_edges() const;  ///< lists the interior edges in the mesh.
    inline std::vector<Vertex*> get_i_vertices() const;  ///< lists the interior vertices in the mesh.
    Cell* i_cell(size_t iC) const;  ///< get a constant pointer to the iC-th interior cell
    Face* i_face(size_t iF) const;   ///< get a constant pointer to the iF-th interior face
    Edge* i_edge(size_t iE) const;   ///< get a constant pointer to the iE-th interior edge
    Vertex* i_vertex(size_t iV) const;   ///< get a constant pointer to the iV-th interior vertex

    inline bool add_cell(Cell* cell);  ///<  adds a cell to the mesh
    inline bool add_face(Face* face);  ///<  adds a face to the mesh
    inline bool add_edge(Edge* edge);  ///<  adds an edge to the mesh
    inline bool add_vertex(Vertex* vertex);  ///<  adds a vertex to the mesh

    inline bool add_b_cell(Cell* cell);  ///<  adds a boundary cell to the mesh
    inline bool add_b_face(Face* face);  ///<  adds a boundary face to the mesh
    inline bool add_b_edge(Edge* edge);  ///<  adds a boundary edge to the mesh
    inline bool add_b_vertex(Vertex* vertex);  ///<  adds a boundary vertex to the mesh

    inline bool add_i_cell(Cell* cell);  ///<  adds an interior cell to the mesh
    inline bool add_i_face(Face* face);  ///<  adds an interior face to the mesh
    inline bool add_i_edge(Edge* edge);  ///<  adds an interior edge to the mesh
    inline bool add_i_vertex(Vertex* vertex);  ///<  adds an interior vertex to the mesh

    inline bool set_h_max(double h_max);    ///< set the mesh size to h_max
                
    double regularity(); ///< returns a regularity factor
                
    /// Re-index the cells, edges or vertices
    void renum(                                                                                     
	       const char B,                                                                                           ///< T for cells, F for faces, E for edges, V for vertices
	       const std::vector<size_t> new_to_old    ///< Vector of new indices to old ones (new_to_old[i]=j: the index formerly j will become i)
												    );

                
  private:
    std::string _mesh_name;

    // primary data: list of cells, faces, edges, vertices...
    std::vector<Cell*> _cells;
    std::vector<Face*> _faces;
    std::vector<Edge*> _edges;
    std::vector<Vertex*> _vertices;
    std::vector<Cell*> _b_cells;
    std::vector<Face*> _b_faces;
    std::vector<Edge*> _b_edges;
    std::vector<Vertex*> _b_vertices;
    std::vector<Cell*> _i_cells;
    std::vector<Face*> _i_faces;
    std::vector<Edge*> _i_edges;
    std::vector<Vertex*> _i_vertices;
    double _h_max;
    const size_t _dim = 3;
        
  };



  // ----------------------------------------------------------------------------
  //                            Implementations
  // ----------------------------------------------------------------------------

  size_t Mesh::n_cells() const { return _cells.size(); }
  size_t Mesh::n_faces() const { return _faces.size(); }
  size_t Mesh::n_edges() const { return _edges.size(); }
  size_t Mesh::n_vertices() const { return _vertices.size(); }
  double Mesh::h_max() const { return _h_max; }
  size_t Mesh::dim() const { return _dim; }
  void Mesh::set_name(std::string name) {
    _mesh_name = name;
    return;
  }
  std::string Mesh::get_name() { return _mesh_name; }

  bool Mesh::add_cell(Cell* cell) {
    _cells.push_back(cell);

    return true;
  }

  bool Mesh::add_b_cell(Cell* cell) {
    _b_cells.push_back(cell);

    return true;
  }

  bool Mesh::add_i_cell(Cell* cell) {
    _i_cells.push_back(cell);

    return true;
  }

  bool Mesh::add_face(Face* face) {
    _faces.push_back(face);

    return true;
  }

  bool Mesh::add_b_face(Face* face) {
    _b_faces.push_back(face);

    return true;
  }

  bool Mesh::add_i_face(Face* face) {
    _i_faces.push_back(face);

    return true;
  }


  bool Mesh::add_edge(Edge* edge) {
    _edges.push_back(edge);

    return true;
  }

  bool Mesh::add_b_edge(Edge* edge) {
    _b_edges.push_back(edge);

    return true;
  }

  bool Mesh::add_i_edge(Edge* edge) {
    _i_edges.push_back(edge);

    return true;
  }

  bool Mesh::add_vertex(Vertex* vertex) {
    _vertices.push_back(vertex);

    return true;
  }

  bool Mesh::add_b_vertex(Vertex* vertex) {
    _b_vertices.push_back(vertex);

    return true;
  }

  bool Mesh::add_i_vertex(Vertex* vertex) {
    _i_vertices.push_back(vertex);

    return true;
  }

  bool Mesh::set_h_max(double h_max) {
    _h_max = h_max;

    return true;
  }

  std::vector<Cell*> Mesh::get_cells() const { return _cells; }
  std::vector<Face*> Mesh::get_faces() const { return _faces; }
  std::vector<Edge*> Mesh::get_edges() const { return _edges; }
  std::vector<Vertex*> Mesh::get_vertices() const { return _vertices; }
  std::vector<Cell*> Mesh::get_b_cells() const { return _b_cells; }
  std::vector<Face*> Mesh::get_b_faces() const { return _b_faces; }
  std::vector<Edge*> Mesh::get_b_edges() const { return _b_edges; }
  std::vector<Vertex*> Mesh::get_b_vertices() const { return _b_vertices; }
  std::vector<Cell*> Mesh::get_i_cells() const { return _i_cells; }
  std::vector<Face*> Mesh::get_i_faces() const { return _i_faces; }
  std::vector<Edge*> Mesh::get_i_edges() const { return _i_edges; }
  std::vector<Vertex*> Mesh::get_i_vertices() const { return _i_vertices; }
  /*@}*/
}

#endif /* MESH_HPP */

