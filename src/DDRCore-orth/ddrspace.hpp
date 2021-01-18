#ifndef DDRSPACE_HPP
#define DDRSPACE_HPP

#include <dofspace.hpp>

namespace HArDCore3D
{
  /*!
   *  \addtogroup DDRcore
   * @{
   */

  /// Base class for DDR spaces
  class DDRSpace : public DOFSpace {
  public:
    /// Constructor
    DDRSpace(
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
    
    //------------------------------------------------------------------------------
    // Restrictions
    //------------------------------------------------------------------------------

    /// Restrict to the edge of index iE
    Eigen::VectorXd restrictEdge(size_t iE, const Eigen::VectorXd & vh) const;

    /// Restrict to the face of index iF
    Eigen::VectorXd restrictFace(size_t iF, const Eigen::VectorXd & vh) const;

    /// Restrict to the cell of index iT
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

    /// Extend a face operator to a cell
    Eigen::MatrixXd extendOperator(const Cell & T, const Face & F, const Eigen::MatrixXd & opF) const;

    /// Extend an edge operator to a cell
    Eigen::MatrixXd extendOperator(const Cell & T, const Edge & E, const Eigen::MatrixXd & opE) const;

    //------------------------------------------------------------------------------
    // Global DOF indices for an element T
    //------------------------------------------------------------------------------

    std::vector<size_t> globalDOFIndices(const Cell & T) const;
    std::vector<size_t> globalDOFIndices(const Face & F) const;    
  };

} // namespace HArDCore3D

#endif
