#include "variabledofspace.hpp"

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructors
//------------------------------------------------------------------------------

VariableDOFSpace::VariableDOFSpace(
                   const Mesh & mesh,
                   const Eigen::VectorXd n_local_vertex_dofs,
                   const Eigen::VectorXd n_local_edge_dofs,
                   const Eigen::VectorXd n_local_face_dofs,
                   const Eigen::VectorXd n_local_cell_dofs		   
                   )
  : m_mesh(mesh),
    m_n_local_vertex_dofs(n_local_vertex_dofs),
    m_n_local_edge_dofs(n_local_edge_dofs),
    m_n_local_face_dofs(n_local_face_dofs),
    m_n_local_cell_dofs(n_local_cell_dofs)
{
  // Do nothing
}

VariableDOFSpace::VariableDOFSpace(
                   const Mesh & mesh,
                   size_t n_local_vertex_dofs,
                   const Eigen::VectorXd n_local_edge_dofs,
                   const Eigen::VectorXd n_local_face_dofs,
                   const Eigen::VectorXd n_local_cell_dofs		   
                   )
  : m_mesh(mesh),
    m_n_local_vertex_dofs(n_local_vertex_dofs * Eigen::VectorXd::Ones(mesh.n_vertices())),
    m_n_local_edge_dofs(n_local_edge_dofs),
    m_n_local_face_dofs(n_local_face_dofs),
    m_n_local_cell_dofs(n_local_cell_dofs)
    
{
  // Do nothing
}

VariableDOFSpace::VariableDOFSpace(
                   const Mesh & mesh,
                   size_t n_local_vertex_dofs,
                   size_t n_local_edge_dofs,
                   const Eigen::VectorXd n_local_face_dofs,
                   const Eigen::VectorXd n_local_cell_dofs		   
                   )
  : m_mesh(mesh),
    m_n_local_vertex_dofs(n_local_vertex_dofs * Eigen::VectorXd::Ones(mesh.n_vertices())),
    m_n_local_edge_dofs(n_local_edge_dofs * Eigen::VectorXd::Ones(mesh.n_edges())),
    m_n_local_face_dofs(n_local_face_dofs),
    m_n_local_cell_dofs(n_local_cell_dofs)
    
{
  // Do nothing
}

VariableDOFSpace::VariableDOFSpace(
                   const Mesh & mesh,
                   size_t n_local_vertex_dofs,
                   size_t n_local_edge_dofs,
                   size_t n_local_face_dofs,
                   const Eigen::VectorXd n_local_cell_dofs		   
                   )
  : m_mesh(mesh),
    m_n_local_vertex_dofs(n_local_vertex_dofs * Eigen::VectorXd::Ones(mesh.n_vertices())),
    m_n_local_edge_dofs(n_local_edge_dofs * Eigen::VectorXd::Ones(mesh.n_edges())),
    m_n_local_face_dofs(n_local_face_dofs * Eigen::VectorXd::Ones(mesh.n_faces())),
    m_n_local_cell_dofs(n_local_cell_dofs)
    
{
  // Do nothing
}

VariableDOFSpace::VariableDOFSpace(
                   const Mesh & mesh,
                   size_t n_local_vertex_dofs,
                   size_t n_local_edge_dofs,
                   size_t n_local_face_dofs,
                   size_t n_local_cell_dofs		   
                   )
  : m_mesh(mesh),
    m_n_local_vertex_dofs(n_local_vertex_dofs * Eigen::VectorXd::Ones(mesh.n_vertices())),
    m_n_local_edge_dofs(n_local_edge_dofs * Eigen::VectorXd::Ones(mesh.n_edges())),
    m_n_local_face_dofs(n_local_face_dofs * Eigen::VectorXd::Ones(mesh.n_faces())),
    m_n_local_cell_dofs(n_local_cell_dofs * Eigen::VectorXd::Ones(mesh.n_cells()))
    
{
  // Do nothing
}

//------------------------------------------------------------------------------
// Restrictions
//------------------------------------------------------------------------------

Eigen::VectorXd VariableDOFSpace::restrictEdge(size_t iE, const Eigen::VectorXd & vh) const
{
  Eigen::VectorXd vE = Eigen::VectorXd::Zero(dimensionEdge(iE));
  const Edge & E = *m_mesh.edge(iE);

  vE.head(numLocalDofsVertex(*E.vertex(0)))
    = vh.segment(globalOffset(*E.vertex(0)), numLocalDofsVertex(*E.vertex(0)));
  vE.segment(numLocalDofsVertex(*E.vertex(0)), numLocalDofsVertex(*E.vertex(1)))
    = vh.segment(globalOffset(*E.vertex(1)), numLocalDofsVertex(*E.vertex(1)));

  vE.tail(numLocalDofsEdge(E)) = vh.segment(globalOffset(E), numLocalDofsEdge(E));
  
  return vE;
}

//------------------------------------------------------------------------------

Eigen::VectorXd VariableDOFSpace::restrictFace(size_t iF, const Eigen::VectorXd & vh) const
{
  Eigen::VectorXd vF = Eigen::VectorXd::Zero(dimensionFace(iF));
  const Face & F = *m_mesh.face(iF);

  for (size_t iV = 0; iV < F.n_vertices(); iV++) {
    const Vertex & V = *F.vertex(iV);
    vF.segment(localOffset(F, V), numLocalDofsVertex(V))
      = vh.segment(globalOffset(V), numLocalDofsVertex(V));
  } // for iV

  for (size_t iE = 0; iE < F.n_edges(); iE++) {
    const Edge & E = *F.edge(iE);      
    vF.segment(localOffset(F, E), numLocalDofsEdge(E))
      = vh.segment(globalOffset(E), numLocalDofsEdge(E));
  } // for iE

   vF.tail(numLocalDofsFace(F))
      = vh.segment(globalOffset(F), numLocalDofsFace(F));
  
  return vF;
}


//------------------------------------------------------------------------------

Eigen::VectorXd VariableDOFSpace::restrictCell(size_t iT, const Eigen::VectorXd & vh) const
{
  Eigen::VectorXd vT = Eigen::VectorXd::Zero(dimensionCell(iT));
  const Cell & T = *m_mesh.cell(iT);

  for (size_t iV = 0; iV < T.n_vertices(); iV++) {
    const Vertex & V = *T.vertex(iV);
    vT.segment(localOffset(T, V), numLocalDofsVertex(V))
      = vh.segment(globalOffset(V), numLocalDofsVertex(V));
  } // for iV

  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);      
    vT.segment(localOffset(T, E), numLocalDofsEdge(E))
      = vh.segment(globalOffset(E), numLocalDofsEdge(E));
  } // for iE

  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    vT.segment(localOffset(T, F), numLocalDofsFace(F))
      = vh.segment(globalOffset(F), numLocalDofsFace(F));
  } // for iF

  vT.tail(numLocalDofsCell(T))
    = vh.segment(globalOffset(T), numLocalDofsCell(T));
  
  return vT;
}


//------------------------------------------------------------------------------

Eigen::MatrixXd VariableDOFSpace::extendOperator(const Cell & T, const Face & F, const Eigen::MatrixXd & opF) const
{
  Eigen::MatrixXd opT = Eigen::MatrixXd::Zero(opF.rows(), dimensionCell(T));

  // Vertex DOFs
  for (size_t iV = 0; iV < F.n_vertices(); iV++) { 
    const Vertex & V = *F.vertex(iV);     
    opT.block(0, localOffset(T, V), opF.rows(), numLocalDofsVertex(V))
      = opF.block(0, localOffset(F, V), opF.rows(), numLocalDofsVertex(V));
  } // for iV
  
  // Edge DOFs
  for (size_t iE = 0; iE < F.n_edges(); iE++) {
    const Edge & E = *F.edge(iE);
    opT.block(0, localOffset(T, E), opF.rows(), numLocalDofsEdge(E))
      = opF.block(0, localOffset(F, E), opF.rows(), numLocalDofsEdge(E));
  } // for iE

  // Face DOFs
  opT.block(0, localOffset(T, F), opF.rows(), numLocalDofsFace(F))
    = opF.block(0, localOffset(F), opF.rows(), numLocalDofsFace(F));
  
  return opT;
}

//------------------------------------------------------------------------------

Eigen::MatrixXd VariableDOFSpace::extendOperator(const Cell & T, const Edge & E, const Eigen::MatrixXd & opE) const
{
  Eigen::MatrixXd opT = Eigen::MatrixXd::Zero(opE.rows(), dimensionCell(T));

  // Vertex DOFs
  for (size_t iV = 0; iV < 2; iV++) { 
    const Vertex & V = *E.vertex(iV);
    opT.block(0, localOffset(T, V), opE.rows(), numLocalDofsVertex(V))
      = opE.block(0, localOffset(E, V), opE.rows(), numLocalDofsVertex(V));
  } // for iV
  
  // Edge DOFs
    opT.block(0, localOffset(T, E), opE.rows(), numLocalDofsEdge(E))
      = opE.block(0, localOffset(E), opE.rows(), numLocalDofsEdge(E));
  
  return opT;
}
//------------------------------------------------------------------------------

Eigen::MatrixXd VariableDOFSpace::extendOperator(const Face & F, const Edge & E, const Eigen::MatrixXd & opE) const
{
  Eigen::MatrixXd opF = Eigen::MatrixXd::Zero(opE.rows(), dimensionFace(F));

  // Vertex DOFs
  for (size_t iV = 0; iV < 2; iV++) { 
    const Vertex & V = *E.vertex(iV);
    opF.block(0, localOffset(F, V), opE.rows(), numLocalDofsVertex(V))
      = opE.block(0, localOffset(E, V), opE.rows(), numLocalDofsVertex(V));
  } // for iV
  
  // Edge DOFs
    opF.block(0, localOffset(F, E), opE.rows(), numLocalDofsEdge(E))
      = opE.block(0, localOffset(E), opE.rows(), numLocalDofsEdge(E));
  
  return opF;
}

//------------------------------------------------------------------------------

std::vector<size_t> VariableDOFSpace::globalDOFIndices(const Cell & T) const
{
  std::vector<size_t> I(dimensionCell(T));

  size_t dof_index = 0;
  
  // Vertex DOFs
  for (size_t iV = 0; iV < T.n_vertices(); iV++) {
    const Vertex & V = *T.vertex(iV);
    size_t offset_V = globalOffset(V);
    for (size_t i = 0; i < numLocalDofsVertex(V); i++, dof_index++) {
      I[dof_index] = offset_V + i;
    } // for i
  } // for iV

  // Edge DOFs
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    size_t offset_E = globalOffset(E);
    for (size_t i = 0; i < numLocalDofsEdge(E); i++, dof_index++) {
      I[dof_index] = offset_E + i;
    } // for i
  } // for iE

  // Face DOFs
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    size_t offset_F = globalOffset(F);
    for (size_t i = 0; i < numLocalDofsFace(F); i++, dof_index++) {
      I[dof_index] = offset_F + i;
    } // for i
  } // for iF

  // Cell DOFs
  size_t offset_T = globalOffset(T);
  for (size_t i = 0; i < numLocalDofsCell(T); i++, dof_index++) {
    I[dof_index] = offset_T + i;
  } // for i

  assert( dimensionCell(T) == dof_index );
  
  return I;
}

//------------------------------------------------------------------------------

std::vector<size_t> VariableDOFSpace::globalDOFIndices(const Face & F) const
{
  std::vector<size_t> I(dimensionFace(F));

  size_t dof_index = 0;
  
  // Vertex DOFs
  for (size_t iV = 0; iV < F.n_vertices(); iV++) {
    const Vertex & V = *F.vertex(iV);
    size_t offset_V = globalOffset(V);
    for (size_t i = 0; i < numLocalDofsVertex(V); i++, dof_index++) {
      I[dof_index] = offset_V + i;
    } // for i
  } // for iV

  // Edge DOFs
  for (size_t iE = 0; iE < F.n_edges(); iE++) {
    const Edge & E = *F.edge(iE);
    size_t offset_E = globalOffset(E);
    for (size_t i = 0; i < numLocalDofsEdge(E); i++, dof_index++) {
      I[dof_index] = offset_E + i;
    } // for i
  } // for iE

  // Face DOFs
  size_t offset_F = globalOffset(F);
  for (size_t i = 0; i < numLocalDofsFace(F); i++, dof_index++) {
    I[dof_index] = offset_F + i;
  } // for i

  assert( dimensionFace(F) == dof_index );
  
  return I;
}
