#ifndef DDRSPACE_HPP
#define DDRSPACE_HPP

#include <mesh.hpp>
#include <vertex.hpp>
#include <edge.hpp>
#include <face.hpp>
#include <cell.hpp>

namespace HArDCore3D
{
  /*!
   *	\addtogroup DDRcore
   * @{
   */

  /// Base class for DDR spaces
  class DDRSpace {
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
    // Number of local DOFs
    //------------------------------------------------------------------------------

    /// Return the number of local vertex DOFs
    inline size_t numLocalDofsVertex() const
    {
      return m_n_local_vertex_dofs;
    }
    
    /// Return the number of local vertex DOFs
    inline size_t numLocalDofsEdge() const
    {
      return m_n_local_edge_dofs;
    }
    
    /// Return the number of local vertex DOFs
    inline size_t numLocalDofsFace() const
    {
      return m_n_local_face_dofs;
    }
    
    /// Return the number of local vertex DOFs
    inline size_t numLocalDofsCell() const
    {
      return m_n_local_cell_dofs;
    }

    //------------------------------------------------------------------------------
    // Dimensions
    //------------------------------------------------------------------------------
    
    /// Return the dimension of the global space
    inline size_t dimension() const
    {
      return m_mesh.n_vertices() * m_n_local_vertex_dofs
	+ m_mesh.n_edges() * m_n_local_edge_dofs
	+ m_mesh.n_faces() * m_n_local_face_dofs
	+ m_mesh.n_cells() * m_n_local_cell_dofs;
    }

    /// Return the dimension of the local space on the vertex V
    inline size_t dimensionVertex(const Vertex & V) const
    {
      return m_n_local_vertex_dofs;
    }

    /// Return the dimension of the local space on the vertex of index iV
    inline size_t dimensionVertex(size_t iV) const
    {
      return dimensionVertex(*m_mesh.vertex(iV));
    }

    /// Return the dimension of the local space on the edge E
    inline size_t dimensionEdge(const Edge & E) const
    {
      return 2 * m_n_local_vertex_dofs
	+ m_n_local_edge_dofs;
    }

    /// Return the dimension of the local space on the edge of index iE
    inline size_t dimensionEdge(size_t iE) const
    {
      return dimensionEdge(*m_mesh.edge(iE));
    }

    /// Return the dimension of the local space on the face F
    inline size_t dimensionFace(const Face & F) const
    {
      return F.n_vertices() * m_n_local_vertex_dofs
	+ F.n_edges() * m_n_local_edge_dofs
	+ m_n_local_face_dofs;
    }
    
    /// Return the dimension of the local space on the face of index iF
    inline size_t dimensionFace(size_t iF) const
    {
      return dimensionFace(*m_mesh.face(iF));
    }

    /// Return the dimension of the local space on the face F
    inline size_t dimensionCell(const Cell & T) const
    {
      return T.n_vertices() * m_n_local_vertex_dofs
	+ T.n_edges() * m_n_local_edge_dofs
	+ T.n_faces() * m_n_local_face_dofs
	+ m_n_local_cell_dofs;
    }
    
    /// Return the dimension of the local space on the face of index iF
    inline size_t dimensionCell(size_t iT) const
    {
      return dimensionCell(*m_mesh.cell(iT));
    }

    //------------------------------------------------------------------------------
    // Offsets
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

    /// Return the local offset of the vertex V with respect to the face F
    inline size_t localOffset(const Face & F, const Vertex & V) const
    {
      return F.index_vertex(&V) * m_n_local_vertex_dofs;
    }

    /// Return the local offset of the edge E with respect to the face F
    inline size_t localOffset(const Face & F, const Edge & E) const
    {
      return F.n_vertices() * m_n_local_vertex_dofs
    	+ F.index_edge(&E) * m_n_local_edge_dofs; 
    }

    /// Return the local offset of the unknowns attached to the face F
    inline size_t localOffset(const Face & F) const
    {
      return F.n_vertices() * m_n_local_vertex_dofs
	+ F.n_edges() * m_n_local_edge_dofs;
    }

    /// Return the local offset of the vertex V with respect to the cell T
    inline size_t localOffset(const Cell & T, const Vertex & V) const
    {
      return T.index_vertex(&V) * m_n_local_vertex_dofs;
    }

    /// Return the local offset of the edge V with respect to the cell T
    inline size_t localOffset(const Cell & T, const Edge & E) const
    {
      return T.n_vertices() * m_n_local_vertex_dofs
	+ T.index_edge(&E) * m_n_local_edge_dofs;
    }

    /// Return the local offset of the face F with respect to the cell T
    inline size_t localOffset(const Cell & T, const Face & F) const
    {
      return T.n_vertices() * m_n_local_vertex_dofs
	+ T.n_edges() * m_n_local_edge_dofs
	+ T.index_face(&F) * m_n_local_face_dofs;
    }

    /// Return the local offset of the unknowns attached to the element T
    inline size_t localOffset(const Cell & T) const
    {
      return T.n_vertices() * m_n_local_vertex_dofs
	+ T.n_edges() * m_n_local_edge_dofs
	+ T.n_faces() * m_n_local_face_dofs;
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
    Eigen::MatrixXd extendOperator(const Cell & T, const Face & F, const Eigen::MatrixXd & opF);

    //------------------------------------------------------------------------------
    // Global DOF indices for an element T
    //------------------------------------------------------------------------------

    std::vector<size_t> globalDOFIndices(const Cell & T) const;
    std::vector<size_t> globalDOFIndices(const Face & F) const;
    
  protected:
    const Mesh & m_mesh;
    size_t m_n_local_vertex_dofs;
    size_t m_n_local_edge_dofs;
    size_t m_n_local_face_dofs;
    size_t m_n_local_cell_dofs;
  };

} // namespace HArDCore3D

#endif
