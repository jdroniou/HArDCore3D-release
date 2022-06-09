#ifndef VARIABLEDOFSPACE_HPP
#define VARIABLEDOFSPACE_HPP

#include<mesh.hpp>

namespace HArDCore3D
{
  /*!
   * \addtogroup Common
   * @{
   */

  /// Base class for global DOF spaces. 
  /** The DOFs are organised by increasing geometric entities dimensions: DOFs of vertices, DOFs of edges, DOFs of faces, DOFs of cells.
  This version allows for cases where the number of DOF on each mesh entity varies from one mesh entity to the next */
  
  class VariableDOFSpace {
  public:
    /// Constructor
    VariableDOFSpace(
             const Mesh & mesh,
             const Eigen::VectorXd n_local_vertex_dofs,
             const Eigen::VectorXd n_local_edge_dofs,
             const Eigen::VectorXd n_local_face_dofs,
             const Eigen::VectorXd n_local_cell_dofs
             );

    /// Simpler constructor if all vertices have the same number of DOFs
    VariableDOFSpace(
             const Mesh & mesh,
             size_t n_local_vertex_dofs,
             const Eigen::VectorXd n_local_edge_dofs,
             const Eigen::VectorXd n_local_face_dofs,
             const Eigen::VectorXd n_local_cell_dofs
             );

    /// Simpler constructor if all vertices/edges have the same number of DOFs
    VariableDOFSpace(
             const Mesh & mesh,
             size_t n_local_vertex_dofs,
             size_t n_local_edge_dofs,
             const Eigen::VectorXd n_local_face_dofs,
             const Eigen::VectorXd n_local_cell_dofs
             );

    /// Simpler constructor if all vertices/edges/faces have the same number of DOFs
    VariableDOFSpace(
             const Mesh & mesh,
             size_t n_local_vertex_dofs,
             size_t n_local_edge_dofs,
             size_t n_local_face_dofs,
             const Eigen::VectorXd n_local_cell_dofs
             );

    /// Simpler constructor if all vertices/edges/faces/cells have the same number of DOFs
    VariableDOFSpace(
             const Mesh & mesh,
             size_t n_local_vertex_dofs,
             size_t n_local_edge_dofs,
             size_t n_local_face_dofs,
             size_t n_local_cell_dofs
             );


    //------------------------------------------------------------------------------
    // Accessors
    //------------------------------------------------------------------------------
    
    /// Returns the mesh
    const Mesh & mesh() const
    {
      return m_mesh;
    }

    /// Returns the number of local DOFs on vertex of index iV
    inline size_t numLocalDofsVertex(const size_t iV) const
    {
      assert(iV<m_mesh.n_vertices());
      return m_n_local_vertex_dofs[iV];
    }
    
    /// Returns the number of local DOFs on vertex V
    inline size_t numLocalDofsVertex(const Vertex & V) const
    {
      return numLocalDofsVertex(V.global_index());
    }

    /// Returns the number of local DOFs on edge of index iE
    inline size_t numLocalDofsEdge(const size_t iE) const
    {
      assert(iE<m_mesh.n_edges());
      return m_n_local_edge_dofs[iE];
    }
    
    /// Returns the number of local DOFs on edge E
    inline size_t numLocalDofsEdge(const Edge & E) const
    {
      return numLocalDofsEdge(E.global_index());
    }

    /// Returns the number of local DOFs on face of index iF
    inline size_t numLocalDofsFace(const size_t iF) const
    {
      assert(iF<m_mesh.n_faces());
      return m_n_local_face_dofs[iF];
    }
    
    /// Returns the number of local DOFs on face F
    inline size_t numLocalDofsFace(const Face & F) const
    {
      return numLocalDofsFace(F.global_index());
    }

    /// Returns the number of local DOFs on cell of index iT
    inline size_t numLocalDofsCell(const size_t iT) const
    {
      assert(iT<m_mesh.n_cells());
      return m_n_local_cell_dofs[iT];
    }

    /// Returns the number of local DOFs on cell T
    inline size_t numLocalDofsCell(const Cell & T) const
    {
      return numLocalDofsCell(T.global_index());
    }

    //------------------------------------------------------------------------------
    // Dimensions
    //------------------------------------------------------------------------------
    
///// We could consider pre-calculating and storing a number of these elements below... will see later if necessary in terms of cost 
   
    /// Total number of vertices DOFs
    inline size_t nDOFs_vertices() const
    {
      return m_n_local_vertex_dofs.sum();
    }
   
    /// Total number of edges DOFs
    inline size_t nDOFs_edges() const
    {
      return m_n_local_edge_dofs.sum();
    }

    /// Total number of faces DOFs
    inline size_t nDOFs_faces() const
    {
      return m_n_local_face_dofs.sum();
    }

    /// Total number of cells DOFs
    inline size_t nDOFs_cells() const
    {
      return m_n_local_cell_dofs.sum();
    }

    /// Returns the dimension of the global space (all DOFs for all geometric entities)
    inline size_t dimension() const
    {
      return nDOFs_vertices() + nDOFs_edges() + nDOFs_faces() + nDOFs_cells();
    }

    //--- Vertices -----//
    /// Returns the dimension of the local space on the vertex V
    inline size_t dimensionVertex(const Vertex & V) const
    {
      return m_n_local_vertex_dofs[V.global_index()];
    }

    /// Returns the dimension of the local space on the vertex of index iV
    inline size_t dimensionVertex(size_t iV) const
    {
      assert(iV<m_mesh.n_vertices());
      return dimensionVertex(*m_mesh.vertex(iV));
    }

    //--- Edges ----//
    /// Returns the dimension of the local space on the edge E (including vertices)
    inline size_t dimensionEdge(const Edge & E) const
    {
      return m_n_local_vertex_dofs[E.vertex(0)->global_index()] + m_n_local_vertex_dofs[E.vertex(1)->global_index()]
        + m_n_local_edge_dofs[E.global_index()];
    }

    /// Returns the dimension of the local space on the edge of index iE (including vertices)
    inline size_t dimensionEdge(size_t iE) const
    {
      assert(iE<m_mesh.n_edges());
      return dimensionEdge(*m_mesh.edge(iE));
    }

    //--- Faces -----//
    /// Returns the dimension of the local space on the face F (including edges and vertices)
    inline size_t dimensionFace(const Face & F) const
    {
      size_t nb_dofs = 0;
      for (Vertex * V : F.get_vertices()){
        nb_dofs += m_n_local_vertex_dofs[V->global_index()];
      }
      for (Edge * E : F.get_edges()){
        nb_dofs += m_n_local_edge_dofs[E->global_index()];
      }
      return nb_dofs + m_n_local_face_dofs[F.global_index()];
    }
    
    /// Returns the dimension of the local space on the face of index iF (including edges and vertices)
    inline size_t dimensionFace(size_t iF) const
    {
      return dimensionFace(*m_mesh.face(iF));
    }

    //--- Cell -----//
    /// Returns the dimension of the local space on the cell T (including faces, edges and vertices)
    inline size_t dimensionCell(const Cell & T) const
    {
      size_t nb_dofs = 0;
      for (Vertex * V : T.get_vertices()){
        nb_dofs += m_n_local_vertex_dofs[V->global_index()];
      }
      for (Edge * E : T.get_edges()){
        nb_dofs += m_n_local_edge_dofs[E->global_index()];
      }
      for (Face * F : T.get_faces()){
        nb_dofs += m_n_local_face_dofs[F->global_index()];
      }
      return nb_dofs + m_n_local_cell_dofs[T.global_index()];
    }
    
    /// Returns the dimension of the local space on the cell of index iT (including faces, edges and vertices)
    inline size_t dimensionCell(size_t iT) const
    {
      return dimensionCell(*m_mesh.cell(iT));
    }

    //------------------------------------------------------------------------------
    // Local offsets
    //------------------------------------------------------------------------------

    /// Returns the local offset of the vertex V with respect to the edge E
    inline size_t localOffset(const Edge & E, const Vertex & V) const
    {
      size_t nb_dofs = 0;
      for (int iV=0; iV<E.index_vertex(&V); iV++){
        nb_dofs += m_n_local_vertex_dofs[E.vertex(iV)->global_index()];
      }
      return nb_dofs;
    }

    /// Returns the local offset of the unknowns attached to the edge E
    inline size_t localOffset(const Edge & E) const
    {
      return m_n_local_vertex_dofs[E.vertex(0)->global_index()] + m_n_local_vertex_dofs[E.vertex(1)->global_index()];
    }
    
    /// Returns the local offset of the vertex V with respect to the face F
    inline size_t localOffset(const Face & F, const Vertex & V) const
    {
      size_t nb_dofs = 0;
      for (int iV=0; iV<F.index_vertex(&V); iV++){
        nb_dofs += m_n_local_vertex_dofs[F.vertex(iV)->global_index()];
      }
      return nb_dofs;
    }

    /// Returns the local offset of the edge E with respect to the face F
    inline size_t localOffset(const Face & F, const Edge & E) const
    {
      size_t nb_dofs = 0;
      for (size_t iV=0; iV<F.n_vertices(); iV++){
        nb_dofs += m_n_local_vertex_dofs[F.vertex(iV)->global_index()];
      }
      for (int iE=0; iE<F.index_edge(&E); iE++){
        nb_dofs += m_n_local_edge_dofs[F.edge(iE)->global_index()];
      }
      return nb_dofs;
    }

    /// Returns the local offset of the unknowns attached to the face F
    inline size_t localOffset(const Face & F) const
    {
      size_t nb_dofs = 0;
      for (size_t iV=0; iV<F.n_vertices(); iV++){
        nb_dofs += m_n_local_vertex_dofs[F.vertex(iV)->global_index()];
      }
      for (size_t iE=0; iE<F.n_edges(); iE++){
        nb_dofs += m_n_local_edge_dofs[F.edge(iE)->global_index()];
      }
      return nb_dofs;
    }

    /// Returns the local offset of the vertex V with respect to the cell T
    inline size_t localOffset(const Cell & T, const Vertex & V) const
    {
      size_t nb_dofs = 0;
      for (int iV=0; iV<T.index_vertex(&V); iV++){
        nb_dofs += m_n_local_vertex_dofs[T.vertex(iV)->global_index()];
      }
      return nb_dofs;
    }

    /// Returns the local offset of the edge E with respect to the cell T
    inline size_t localOffset(const Cell & T, const Edge & E) const
    {
      size_t nb_dofs = 0;
      for (size_t iV=0; iV<T.n_vertices(); iV++){
        nb_dofs += m_n_local_vertex_dofs[T.vertex(iV)->global_index()];
      }
      for (int iE=0; iE<T.index_edge(&E); iE++){
        nb_dofs += m_n_local_edge_dofs[T.edge(iE)->global_index()];
      }
      return nb_dofs;
    }

    /// Returns the local offset of the face F with respect to the cell T
    inline size_t localOffset(const Cell & T, const Face & F) const
    {
      size_t nb_dofs = 0;
      for (size_t iV=0; iV<T.n_vertices(); iV++){
        nb_dofs += m_n_local_vertex_dofs[T.vertex(iV)->global_index()];
      }
      for (size_t iE=0; iE<T.n_edges(); iE++){
        nb_dofs += m_n_local_edge_dofs[T.edge(iE)->global_index()];
      }
      for (int iF=0; iF<T.index_face(&F); iF++){
        nb_dofs += m_n_local_face_dofs[T.face(iF)->global_index()];
      }
      return nb_dofs;
    }

    /// Returns the local offset of the unknowns attached to the element T
    inline size_t localOffset(const Cell & T) const
    {
      size_t nb_dofs = 0;
      for (size_t iV=0; iV<T.n_vertices(); iV++){
        nb_dofs += m_n_local_vertex_dofs[T.vertex(iV)->global_index()];
      }
      for (size_t iE=0; iE<T.n_edges(); iE++){
        nb_dofs += m_n_local_edge_dofs[T.edge(iE)->global_index()];
      }
      for (size_t iF=0; iF<T.n_faces(); iF++){
        nb_dofs += m_n_local_face_dofs[T.face(iF)->global_index()];
      }
      return nb_dofs;
    }

    //------------------------------------------------------------------------------
    // Global offsets
    //------------------------------------------------------------------------------
    
    /// Return the global offset for the unknowns on the vertex V
    inline size_t globalOffset(const Vertex & V) const
    {
      return ( m_n_local_vertex_dofs.head(V.global_index()) ).sum();
    }

    /// Return the global offset for the unknowns on the edge E
    inline size_t globalOffset(const Edge & E) const
    {
      return m_n_local_vertex_dofs.sum()
          + ( m_n_local_edge_dofs.head(E.global_index()) ).sum();
    }

    /// Return the global offset for the unknowns on the face F
    inline size_t globalOffset(const Face & F) const
    {
      return m_n_local_vertex_dofs.sum()
          + m_n_local_edge_dofs.sum()
          + ( m_n_local_face_dofs.head(F.global_index()) ).sum();
    }

    /// Return the global offset for the unknowns on the cell T
    inline size_t globalOffset(const Cell & T) const
    {
      return m_n_local_vertex_dofs.sum()
          + m_n_local_edge_dofs.sum()
          + m_n_local_face_dofs.sum()
          + ( m_n_local_cell_dofs.head(T.global_index()) ).sum();
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

    //------------------------------------------------------------------------------
    // Global DOF indices for an element T
    //------------------------------------------------------------------------------

    /// Returns a vector listing the global DOFs attached to the element T: vertex DOFs, edge DOFs, face DOFs and element DOFs
    std::vector<size_t> globalDOFIndices(const Cell & T) const;

    /// Returns a vector listing the global DOFs attached to the face F: vertex DOFs, edge DOFs, face DOFs
    std::vector<size_t> globalDOFIndices(const Face & F) const;  
    
    private:
    const Mesh & m_mesh;
    Eigen::VectorXd m_n_local_vertex_dofs;
    Eigen::VectorXd m_n_local_edge_dofs;
    Eigen::VectorXd m_n_local_face_dofs;
    Eigen::VectorXd m_n_local_cell_dofs;

      
  };

} // namespace HArDCore3D

#endif
