#include "globaldofspace.hpp"

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

GlobalDOFSpace::GlobalDOFSpace(
                   const Mesh & mesh,
                   size_t n_local_vertex_dofs,
                   size_t n_local_edge_dofs,
                   size_t n_local_face_dofs,
                   size_t n_local_cell_dofs		   
                   )
  : LocalDOFSpace(mesh, n_local_vertex_dofs, n_local_edge_dofs, n_local_face_dofs, n_local_cell_dofs)
{
  // Do nothing
}

//------------------------------------------------------------------------------
// Restrictions
//------------------------------------------------------------------------------

Eigen::VectorXd GlobalDOFSpace::restrictEdge(size_t iE, const Eigen::VectorXd & vh) const
{
  Eigen::VectorXd vE = Eigen::VectorXd::Zero(dimensionEdge(iE));
  const Edge & E = *m_mesh.edge(iE);

  if (m_n_local_vertex_dofs > 0) {
    vE.head(m_n_local_vertex_dofs)
      = vh.segment(globalOffset(*E.vertex(0)), m_n_local_vertex_dofs);
    vE.segment(m_n_local_vertex_dofs, m_n_local_vertex_dofs)
      = vh.segment(globalOffset(*E.vertex(1)), m_n_local_vertex_dofs);
  }

  if (m_n_local_edge_dofs > 0) {
    vE.tail(m_n_local_edge_dofs) = vh.segment(globalOffset(E), m_n_local_edge_dofs);
  }
  
  return vE;
}

//------------------------------------------------------------------------------

Eigen::VectorXd GlobalDOFSpace::restrictFace(size_t iF, const Eigen::VectorXd & vh) const
{
  Eigen::VectorXd vF = Eigen::VectorXd::Zero(dimensionFace(iF));
  const Face & F = *m_mesh.face(iF);

  if (m_n_local_vertex_dofs > 0) {
    for (size_t iV = 0; iV < F.n_vertices(); iV++) {
      const Vertex & V = *F.vertex(iV);
      vF.segment(localOffset(F, V), m_n_local_vertex_dofs)
        = vh.segment(globalOffset(V), m_n_local_vertex_dofs);
    } // for iV
  }

  if (m_n_local_edge_dofs > 0) {
    for (size_t iE = 0; iE < F.n_edges(); iE++) {
      const Edge & E = *F.edge(iE);      
      vF.segment(localOffset(F, E), m_n_local_edge_dofs)
        = vh.segment(globalOffset(E), m_n_local_edge_dofs);
    } // for iE
  }

  if (m_n_local_face_dofs > 0) {
    vF.tail(m_n_local_face_dofs)
      = vh.segment(globalOffset(F), m_n_local_face_dofs);
  }
  
  return vF;
}

//------------------------------------------------------------------------------

Eigen::VectorXd GlobalDOFSpace::restrictCell(size_t iT, const Eigen::VectorXd & vh) const
{
  Eigen::VectorXd vT = Eigen::VectorXd::Zero(dimensionCell(iT));
  const Cell & T = *m_mesh.cell(iT);

  if (m_n_local_vertex_dofs > 0) {
    for (size_t iV = 0; iV < T.n_vertices(); iV++) {
      const Vertex & V = *T.vertex(iV);
      vT.segment(localOffset(T, V), m_n_local_vertex_dofs)
        = vh.segment(globalOffset(V), m_n_local_vertex_dofs);
    } // for iV
  }

  if (m_n_local_edge_dofs > 0) {
    for (size_t iE = 0; iE < T.n_edges(); iE++) {
      const Edge & E = *T.edge(iE);      
      vT.segment(localOffset(T, E), m_n_local_edge_dofs)
        = vh.segment(globalOffset(E), m_n_local_edge_dofs);
    } // for iE
  }

  if (m_n_local_face_dofs > 0) {
    for (size_t iF = 0; iF < T.n_faces(); iF++) {
      const Face & F = *T.face(iF);
      vT.segment(localOffset(T, F), m_n_local_face_dofs)
        = vh.segment(globalOffset(F), m_n_local_face_dofs);
    } // for iF
  }

  vT.tail(m_n_local_cell_dofs)
    = vh.segment(globalOffset(T), m_n_local_cell_dofs);
  
  return vT;
}

//------------------------------------------------------------------------------

Eigen::MatrixXd GlobalDOFSpace::extendOperator(const Cell & T, const Face & F, const Eigen::MatrixXd & opF) const
{
  Eigen::MatrixXd opT = Eigen::MatrixXd::Zero(opF.rows(), dimensionCell(T));

  // Vertex DOFs
  if (m_n_local_vertex_dofs > 0) {
    for (size_t iV = 0; iV < F.n_vertices(); iV++) { 
      const Vertex & V = *F.vertex(iV);     
      opT.block(0, localOffset(T, V), opF.rows(), m_n_local_vertex_dofs)
        = opF.block(0, localOffset(F, V), opF.rows(), m_n_local_vertex_dofs);
    } // for iV
  }
  
  // Edge DOFs
  if (m_n_local_edge_dofs > 0) {
    for (size_t iE = 0; iE < F.n_edges(); iE++) {
      const Edge & E = *F.edge(iE);
      opT.block(0, localOffset(T, E), opF.rows(), m_n_local_edge_dofs)
        = opF.block(0, localOffset(F, E), opF.rows(), m_n_local_edge_dofs);
    } // for iE
  }

  // Face DOFs
  if (m_n_local_face_dofs > 0) {
    opT.block(0, localOffset(T, F), opF.rows(), m_n_local_face_dofs)
      = opF.block(0, localOffset(F), opF.rows(), m_n_local_face_dofs);
  }
  
  return opT;
}

//------------------------------------------------------------------------------

Eigen::MatrixXd GlobalDOFSpace::extendOperator(const Cell & T, const Edge & E, const Eigen::MatrixXd & opE) const
{
  Eigen::MatrixXd opT = Eigen::MatrixXd::Zero(opE.rows(), dimensionCell(T));

  // Vertex DOFs
  if (m_n_local_vertex_dofs > 0) {
    for (size_t iV = 0; iV < 2; iV++) { 
      const Vertex & V = *E.vertex(iV);
      opT.block(0, localOffset(T, V), opE.rows(), m_n_local_vertex_dofs)
        = opE.block(0, localOffset(E, V), opE.rows(), m_n_local_vertex_dofs);
    } // for iV
  }
  
  // Edge DOFs
  if (m_n_local_edge_dofs > 0) {
      opT.block(0, localOffset(T, E), opE.rows(), m_n_local_edge_dofs)
        = opE.block(0, localOffset(E), opE.rows(), m_n_local_edge_dofs);
  }
  
  return opT;
}
//------------------------------------------------------------------------------

Eigen::MatrixXd GlobalDOFSpace::extendOperator(const Face & F, const Edge & E, const Eigen::MatrixXd & opE) const
{
  Eigen::MatrixXd opF = Eigen::MatrixXd::Zero(opE.rows(), dimensionFace(F));

  // Vertex DOFs
  if (m_n_local_vertex_dofs > 0) {
    for (size_t iV = 0; iV < 2; iV++) { 
      const Vertex & V = *E.vertex(iV);
      opF.block(0, localOffset(F, V), opE.rows(), m_n_local_vertex_dofs)
        = opE.block(0, localOffset(E, V), opE.rows(), m_n_local_vertex_dofs);
    } // for iV
  }
  
  // Edge DOFs
  if (m_n_local_edge_dofs > 0) {
      opF.block(0, localOffset(F, E), opE.rows(), m_n_local_edge_dofs)
        = opE.block(0, localOffset(E), opE.rows(), m_n_local_edge_dofs);
  }
  
  return opF;
}

Eigen::MatrixXd GlobalDOFSpace::extendOperator(size_t d1, size_t i1, size_t d2, size_t i2, const Eigen::MatrixXd & op) const 
{
  assert((d1==3 && (d2==1||d2==2)) || (d1==2 && d2==1) || (d1==1 && d2==0));
  if (d1==1) {
    Eigen::MatrixXd opExtended = Eigen::MatrixXd::Zero(op.rows(),dimensionEdge(i1));
    assert(i1 < m_mesh.n_edges());
    const Edge &E = *m_mesh.edge(i1);
    if (m_n_local_vertex_dofs > 0) {
      const Vertex &V = *m_mesh.vertex(i2);
      opExtended.middleCols(localOffset(E,V), m_n_local_vertex_dofs) = op.middleCols(0,m_n_local_vertex_dofs);
    }
    return opExtended;
  } else if (d1==2) {
    assert(i1 < m_mesh.n_faces());
    assert(i2 < m_mesh.n_edges());
    const Face &F = *m_mesh.face(i1);
    const Edge &E = *m_mesh.edge(i2);
    return extendOperator(F,E,op);
  } else {
    assert(d1==3);
    assert(i1 < m_mesh.n_cells());
    const Cell &T = *m_mesh.cell(i1);
    if (d2==2) {
      assert(i2 < m_mesh.n_faces());
      const Face &F = *m_mesh.face(i2);
      return extendOperator(T,F,op);
    } else {
      assert(d2==1);
      assert(i2 < m_mesh.n_edges());
      const Edge &E = *m_mesh.edge(i2);
      return extendOperator(T,E,op);
    }
  }
}

//------------------------------------------------------------------------------

void GlobalDOFSpace::addInnerProductContribution(const Cell & T, const Face & F, Eigen::MatrixXd & prodT, const Eigen::MatrixXd & prodF) const
{
  // vertex contributions
  for (size_t iV=0; iV<F.n_vertices(); iV++){
    const Vertex & V = *F.vertex(iV);
    
    if (numLocalDofsVertex()>0){
      // vertex-vertex contributions (same vertex, and two different vertices)
      prodT.block(localOffset(T, V), localOffset(T, V), numLocalDofsVertex(), numLocalDofsVertex()) 
          += prodF.block(localOffset(F, V), localOffset(F, V), numLocalDofsVertex(), numLocalDofsVertex());

      for (size_t iVp=0; iVp<iV; iVp++){
        const Vertex & Vp = *F.vertex(iVp);
        prodT.block(localOffset(T, V), localOffset(T, Vp), numLocalDofsVertex(), numLocalDofsVertex()) 
          += prodF.block(localOffset(F, V), localOffset(F, Vp), numLocalDofsVertex(), numLocalDofsVertex());
        prodT.block(localOffset(T, Vp), localOffset(T, V), numLocalDofsVertex(), numLocalDofsVertex()) 
          += prodF.block(localOffset(F, Vp), localOffset(F, V), numLocalDofsVertex(), numLocalDofsVertex());
      }
    
      // vertex-edge contributions
      for (size_t iE=0; iE<F.n_edges(); iE++){
        const Edge & E = *F.edge(iE);
        prodT.block(localOffset(T, E), localOffset(T, V), numLocalDofsEdge(), numLocalDofsVertex()) 
          += prodF.block(localOffset(F, E), localOffset(F, V), numLocalDofsEdge(), numLocalDofsVertex());
        prodT.block(localOffset(T, V), localOffset(T, E), numLocalDofsVertex(), numLocalDofsEdge()) 
          += prodF.block(localOffset(F, V), localOffset(F, E), numLocalDofsVertex(), numLocalDofsEdge());
      }

      // vertex-face contributions
      prodT.block(localOffset(T, F), localOffset(T, V), numLocalDofsFace(), numLocalDofsVertex()) 
        += prodF.block(localOffset(F), localOffset(F, V), numLocalDofsFace(), numLocalDofsVertex());
      prodT.block(localOffset(T, V), localOffset(T, F), numLocalDofsVertex(), numLocalDofsFace()) 
        += prodF.block(localOffset(F, V), localOffset(F), numLocalDofsVertex(), numLocalDofsFace());
    }
  }
  
  // edge contributions
  for (size_t iE=0; iE<F.n_edges(); iE++){
    const Edge & E = *F.edge(iE);

    if (numLocalDofsEdge()>0){
      // edge-edge contributions (same edge, and two different edges)
      prodT.block(localOffset(T, E), localOffset(T, E), numLocalDofsEdge(), numLocalDofsEdge()) 
          += prodF.block(localOffset(F, E), localOffset(F, E), numLocalDofsEdge(), numLocalDofsEdge());

      for (size_t iEp=0; iEp<iE; iEp++){
        const Edge & Ep = *F.edge(iEp);
        prodT.block(localOffset(T, E), localOffset(T, Ep), numLocalDofsEdge(), numLocalDofsEdge()) 
          += prodF.block(localOffset(F, E), localOffset(F, Ep), numLocalDofsEdge(), numLocalDofsEdge());
        prodT.block(localOffset(T, Ep), localOffset(T, E), numLocalDofsEdge(), numLocalDofsEdge()) 
          += prodF.block(localOffset(F, Ep), localOffset(F, E), numLocalDofsEdge(), numLocalDofsEdge());
      }
    
      // edge-face contributions
      prodT.block(localOffset(T, F), localOffset(T, E), numLocalDofsFace(), numLocalDofsEdge()) 
        += prodF.block(localOffset(F), localOffset(F, E), numLocalDofsFace(), numLocalDofsEdge());
      prodT.block(localOffset(T, E), localOffset(T, F), numLocalDofsEdge(), numLocalDofsFace()) 
        += prodF.block(localOffset(F, E), localOffset(F), numLocalDofsEdge(), numLocalDofsFace());
    }
  }
  
  // face-face contibutions
  prodT.block(localOffset(T, F), localOffset(T, F), numLocalDofsFace(), numLocalDofsFace()) 
          += prodF.block(localOffset(F), localOffset(F), numLocalDofsFace(), numLocalDofsFace());
          
}

//------------------------------------------------------------------------------

std::vector<size_t> GlobalDOFSpace::globalDOFIndices(const Cell & T) const
{
  std::vector<size_t> I(dimensionCell(T));

  size_t dof_index = 0;
  
  // Vertex DOFs
  if (m_n_local_vertex_dofs) {
    for (size_t iV = 0; iV < T.n_vertices(); iV++) {
      size_t offset_V = globalOffset(*T.vertex(iV));
      for (size_t i = 0; i < m_n_local_vertex_dofs; i++, dof_index++) {
        I[dof_index] = offset_V + i;
      } // for i
    } // for iV
  }

  // Edge DOFs
  if (m_n_local_edge_dofs) {
    for (size_t iE = 0; iE < T.n_edges(); iE++) {
      size_t offset_E = globalOffset(*T.edge(iE));
      for (size_t i = 0; i < m_n_local_edge_dofs; i++, dof_index++) {
        I[dof_index] = offset_E + i;
      } // for i
    } // for iE
  }

  // Face DOFs
  if (m_n_local_face_dofs) {
    for (size_t iF = 0; iF < T.n_faces(); iF++) {
      size_t offset_F = globalOffset(*T.face(iF));
      for (size_t i = 0; i < m_n_local_face_dofs; i++, dof_index++) {
        I[dof_index] = offset_F + i;
      } // for i
    } // for iF
  }

  // Cell DOFs
  if (m_n_local_cell_dofs) {
    size_t offset_T = globalOffset(T);
    for (size_t i = 0; i < m_n_local_cell_dofs; i++, dof_index++) {
      I[dof_index] = offset_T + i;
    } // for i
  }

  assert( dimensionCell(T) == dof_index );
  
  return I;
}

//------------------------------------------------------------------------------

std::vector<size_t> GlobalDOFSpace::globalDOFIndices(const Face & F) const
{
  std::vector<size_t> I(dimensionFace(F));

  size_t dof_index = 0;
  
  // Vertex DOFs
  if (m_n_local_vertex_dofs) {
    for (size_t iV = 0; iV < F.n_vertices(); iV++) {
      size_t offset_V = globalOffset(*F.vertex(iV));
      for (size_t i = 0; i < m_n_local_vertex_dofs; i++, dof_index++) {
        I[dof_index] = offset_V + i;
      } // for i
    } // for iV
  }

  // Edge DOFs
  if (m_n_local_edge_dofs) {
    for (size_t iE = 0; iE < F.n_edges(); iE++) {
      size_t offset_E = globalOffset(*F.edge(iE));
      for (size_t i = 0; i < m_n_local_edge_dofs; i++, dof_index++) {
        I[dof_index] = offset_E + i;
      } // for i
    } // for iE
  }

  // Face DOFs
  if (m_n_local_face_dofs) {
    size_t offset_F = globalOffset(F);
    for (size_t i = 0; i < m_n_local_face_dofs; i++, dof_index++) {
      I[dof_index] = offset_F + i;
    } // for i
  }

  assert( dimensionFace(F) == dof_index );
  
  return I;
}
