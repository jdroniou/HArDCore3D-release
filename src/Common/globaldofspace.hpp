#ifndef GLOBALDOFSPACE_HPP
#define GLOBALDOFSPACE_HPP

#include <localdofspace.hpp>

namespace HArDCore3D
{
  /*!
   * \addtogroup Common
   * @{
   */

  /// Base class for global DOF spaces. Provides functions to manipulate global DOFs (the local version being provided by DOFSpace).
  /** The DOFs are organised by increasing geometric entities dimensions: DOFs of vertices, DOFs of edges, DOFs of faces, DOFs of cells. */
  
  class GlobalDOFSpace : public LocalDOFSpace {
  public:
    /// Constructor
    GlobalDOFSpace(
             const Mesh & mesh,
             size_t n_local_vertex_dofs,
             size_t n_local_edge_dofs,
             size_t n_local_face_dofs,
             size_t n_local_cell_dofs
             );

    //------------------------------------------------------------------------------
    // Global offsets
    //------------------------------------------------------------------------------
    
    /// Return the global offset for the unknowns on the vertex V
    inline size_t globalOffset(const Vertex & V) const
    {
      return V.global_index() * m_n_local_vertex_dofs;

    }

    /// Return the global offset for the unknowns on the edge E
    inline size_t globalOffset(const Edge & E) const
    {
      return m_mesh.n_vertices() * m_n_local_vertex_dofs
        + E.global_index() * m_n_local_edge_dofs;
    }

    /// Return the global offset for the unknowns on the face F
    inline size_t globalOffset(const Face & F) const
    {
      return m_mesh.n_vertices() * m_n_local_vertex_dofs
        + m_mesh.n_edges() * m_n_local_edge_dofs
        + F.global_index() * m_n_local_face_dofs;
    }

    /// Return the global offset for the unknowns on the cell T
    inline size_t globalOffset(const Cell & T) const
    {
      return m_mesh.n_vertices() * m_n_local_vertex_dofs
        + m_mesh.n_edges() * m_n_local_edge_dofs
        + m_mesh.n_faces() * m_n_local_face_dofs
        + T.global_index() * m_n_local_cell_dofs;
    }

    /// Return the global offset for the unknows on the i-th element of dimension d
    inline size_t globalOffset(size_t d,size_t i) const 
    {
      size_t rv = i*numLocalDofs(d);
      switch(d) {
        case(3):
          rv += m_mesh.n_faces() * m_n_local_face_dofs;
        case(2):
          rv += m_mesh.n_edges() * m_n_local_edge_dofs;
        case(1):
          rv += m_mesh.n_vertices() * m_n_local_vertex_dofs;
        default:
          ;
      }
      return rv;
    }
    
    //------------------------------------------------------------------------------
    // Restrictions
    //------------------------------------------------------------------------------

    /// Restrict to the edge (including its vertices) of index iE
    Eigen::VectorXd restrictEdge(size_t iE, const Eigen::VectorXd & vh) const;

    /// Restrict to the face (including vertices and edges) of index iF
    Eigen::VectorXd restrictFace(size_t iF, const Eigen::VectorXd & vh) const;

    /// Restrict to the cell (including vertices, edges and faces) of index iT
    Eigen::VectorXd restrictCell(size_t iT, const Eigen::VectorXd & vh) const;

    /// Restrict to an edge
    inline Eigen::VectorXd restrict(const Edge & E, const Eigen::VectorXd vh) const
    {
      return restrictEdge(E.global_index(), vh);
    }

    /// Restrict to a face
    inline Eigen::VectorXd restrict(const Face & F, const Eigen::VectorXd vh) const
    {
      return restrictFace(F.global_index(), vh);
    }

    /// Restrict to a cell
    inline Eigen::VectorXd restrict(const Cell & T, const Eigen::VectorXd vh) const
    {
      return restrictCell(T.global_index(), vh);
    }
    
    //------------------------------------------------------------------------------
    // Extensions
    //------------------------------------------------------------------------------

    /// Extend a face operator to a cell: starting from a matrix acting on the DOFs of F (and its edges and vertices),
    /// redistribute the coefficients in a matrix acting on the DOFs viewed from T.
    Eigen::MatrixXd extendOperator(const Cell & T, const Face & F, const Eigen::MatrixXd & opF) const;

    /// Extend an edge operator to a cell
    Eigen::MatrixXd extendOperator(const Cell & T, const Edge & E, const Eigen::MatrixXd & opE) const;

    /// Extend an edge operator to a face
    Eigen::MatrixXd extendOperator(const Face & F, const Edge & E, const Eigen::MatrixXd & opE) const;

    /// Generic extension operator from the i2-th d2-cell to the i1-th d1-cell
    Eigen::MatrixXd extendOperator(size_t d1, size_t i1, size_t d2, size_t i2, const Eigen::MatrixXd & op) const;

    /// Takes an inner product prodF on a face F, and adds its contributions to the inner product prodT on the element T (distributes the contributions according to the DOFs as seen from T)
    void addInnerProductContribution(const Cell & T, const Face & F, Eigen::MatrixXd & prodT, const Eigen::MatrixXd & prodF) const;

    //------------------------------------------------------------------------------
    // Global DOF indices for an element T
    //------------------------------------------------------------------------------

    /// Returns a vector listing the global DOFs attached to the element T: vertex DOFs, edge DOFs, face DOFs and element DOFs
    std::vector<size_t> globalDOFIndices(const Cell & T) const;

    /// Returns a vector listing the global DOFs attached to the face F: vertex DOFs, edge DOFs, face DOFs
    std::vector<size_t> globalDOFIndices(const Face & F) const;    
  };

} // namespace HArDCore3D

#endif
