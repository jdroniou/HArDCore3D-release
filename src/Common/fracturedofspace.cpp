#include "fracturedofspace.hpp"
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
// #include <fractureinterpolator.cpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructors
//------------------------------------------------------------------------------

FractureDOFSpace::FractureDOFSpace(
                   const VSXGrad & vsxgrad,
                   const ConnectedInformationGenerator & connectedinfo,		
                   const FractureInterpolator & fracture_interpolator,
                   bool use_threads   
                   )
  : VariableDOFSpace(
	     vsxgrad.mesh(),
	     [&] {
          Eigen::VectorXi v(vsxgrad.mesh().n_vertices());
          for (size_t i = 0; i < vsxgrad.mesh().n_vertices(); ++i){
              v[i] = connectedinfo.vertex_final_connected_size.at(vsxgrad.mesh().vertex(i));
          }
          return Eigen::VectorXi((dimspace * v).matrix());
      }(),
        [&] {
        Eigen::VectorXi e(vsxgrad.mesh().n_edges());
        for (size_t i = 0; i < vsxgrad.mesh().n_edges(); ++i){
            e[i] = connectedinfo.edge_final_connected_size.at(vsxgrad.mesh().edge(i));
        }
        return Eigen::VectorXi(dimspace * PolynomialSpaceDimension<Edge>::Poly(vsxgrad.degree()-1) * e.matrix());
        }(),
        [&] {
        Eigen::VectorXi f(vsxgrad.mesh().n_faces());
        for (size_t i = 0; i < vsxgrad.mesh().n_faces(); ++i){
            f[i] = connectedinfo.face_final_connected_size.at(vsxgrad.mesh().face(i));
        }
        return Eigen::VectorXi(dimspace * PolynomialSpaceDimension<Face>::Poly(vsxgrad.degree()-1) * f.matrix());
        }(),
	     dimspace * PolynomialSpaceDimension<Cell>::Poly(vsxgrad.degree()-1) * Eigen::VectorXi::Ones(vsxgrad.mesh().n_cells())),
    m_vsxgrad(vsxgrad),
    // m_sxgrad(sxgrad),
    m_connectedinfo(connectedinfo),
    m_fracture_interpolator(fracture_interpolator),
    m_use_threads(use_threads)
    
{
  // Do nothing
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd FractureDOFSpace::interpolate(const FunctionType & q, const double & lambda, const int doe_cell, const int doe_face, const int doe_edge) const
{
  Eigen::MatrixXd q_M = Eigen::MatrixXd(dimspace, m_fracture_interpolator.dimension());

  for (size_t i = 0; i < dimspace; i++){
    std::function<double(const VectorRd &, const Cell *)> q_i = 
              [&q, i, &lambda](const Eigen::Vector3d &x, const Cell *T) -> double { return q(x,T,lambda).dot(basisRd[i]); };
    q_M.row(i) = m_fracture_interpolator.interpolate(q_i, doe_cell, doe_face, doe_edge);
  }

  Eigen::VectorXd qh(Eigen::Map<Eigen::VectorXd>(q_M.data(), q_M.size()));
  
  return qh;
}

//------------------------------------------------------------------------------
// Correction_Term
//------------------------------------------------------------------------------
Eigen::VectorXd FractureDOFSpace::correction_term(const FunctionType & q,const double & lambda, const Eigen::VectorXd & v, const int doe_face) const
{
  Eigen::VectorXd corr = Eigen::VectorXd::Zero(dimension());
  size_t dqr_face = (doe_face >= 0 ? doe_face : 2 * degree() + 3);

  for(size_t iT=0;iT<mesh().n_cells();++iT){
    const Cell & T = *mesh().cell(iT);
    const Cell * cell_ptr=mesh().cell(iT);  
    //------------------------------------------------------------------------------
    // Left-hand side matrix
    MonomialCellIntegralsType int_mono_2km2 = IntegrateCellMonomials(T, 2*degree()-2);
    ShiftedBasis<DDRCore::PolyBasisCellType> phi_basis(*m_vsxgrad.sXgrad().cellBases(iT).Polyk, 1);
    GradientBasis<ShiftedBasis<DDRCore::PolyBasisCellType>> grad_Pk_T(phi_basis);
    Eigen::MatrixXd M=GramMatrix(T,grad_Pk_T,*m_vsxgrad.cellBases(iT).Polykmo3,int_mono_2km2);
    //------------------------------------------------------------------------------
    // Right-hand side matrix
    Eigen::VectorXd b=Eigen::VectorXd::Zero(phi_basis.dimension());
    for(size_t iF=0; iF<T.n_faces(); ++iF){
      const Face & F=*T.face(iF);
      size_t global_iF=F.global_index();
      // int num=m_connectedinfo.face_component(global_iF,iT);
      QuadratureRule quad_F= generate_quadrature_rule(F,dqr_face);
      auto phi_quad_F = evaluate_quad<Function>::compute(phi_basis, quad_F);
      Eigen::Vector3d nF = T.face_normal(T.index_face(&F));
      std::function<double(const Eigen::Vector3d &)> 
            q_nF = [q, nF, cell_ptr,lambda](const Eigen::Vector3d & x)->double {
                         return q(x,cell_ptr,lambda).dot(nF);
                         };
      b+=-integrate(q_nF,phi_quad_F,quad_F);
      auto quad_P2F3_F=evaluate_quad<Function>::compute(*m_vsxgrad.faceBases(global_iF).Polykpo3,quad_F);
      auto quad_P2F3_nF_F=scalar_product(quad_P2F3_F,nF);
      Eigen::VectorXd selected(globalDOFIndices(F,T).size());
      for (size_t i = 0; i < globalDOFIndices(F,T).size(); ++i) {
          selected(i) = v(globalDOFIndices(F,T)[i]);
      }
      b+=compute_gram_matrix(phi_quad_F,quad_P2F3_nF_F,quad_F)*m_vsxgrad.facePotential(F)*selected;
    }
    corr.segment(globalOffset(T),m_vsxgrad.numLocalDofsCell(T))=M.partialPivLu().solve(b);
  }
  return corr;
}

//------------------------------------------------------------------------------
// Number of Fracture DOFs
//------------------------------------------------------------------------------
size_t FractureDOFSpace::fracture_dofs() const{
  return dimspace * PolynomialSpaceDimension<Face>::Poly(m_vsxgrad.degree()-1) * m_connectedinfo.fracture_face_indices_pointer.size();
}

//------------------------------------------------------------------------------
// Number of Fracture DOFs on each cell
//------------------------------------------------------------------------------
Eigen::VectorXd FractureDOFSpace::cell_fracture_dofs() const {
    Eigen::VectorXd cell_fracture(mesh().n_cells());
    for (size_t iT = 0; iT < mesh().n_cells(); ++iT) {
        cell_fracture(iT) = dimspace 
                  * PolynomialSpaceDimension<Face>::Poly(m_vsxgrad.degree() - 1) 
                  * m_connectedinfo.Cell_fractures.at(mesh().cell(iT)).first;
    }
    return cell_fracture;
}


/// Computes the values of the potential reconstruction at the mesh vertices
std::vector<VectorRd> FractureDOFSpace::computeVertexValues(const Eigen::VectorXd & u) const
{
  std::vector<VectorRd> values(mesh().n_vertices(), VectorRd::Zero());
  
  // Value at each vertex obtained averaging the values from all the cells around
  for (Vertex * V : mesh().get_vertices()){
    size_t iV = V->global_index();
 
    for (Cell * T : V->get_cells()){
      size_t iT = T->global_index();
      Eigen::VectorXd pTuT = m_vsxgrad.cellPotential(iT) * restrictnew(*T, u);
      
      for (size_t i=0; i < m_vsxgrad.cellBases(iT).Polykpo3->dimension(); i++){
        values[iV] += pTuT(i) * m_vsxgrad.cellBases(iT).Polykpo3->function(i, V->coords());
      }
    }
    
    values[iV] /= V->n_cells();

  }
  
  return values;
} 


//------------------------------------------------------------------------------
// Restrictions
//------------------------------------------------------------------------------

Eigen::VectorXd FractureDOFSpace::restrictEdge(size_t iE, const Eigen::VectorXd & vh, size_t iT) const
{
  Eigen::VectorXd vE = Eigen::VectorXd::Zero(m_vsxgrad.dimensionEdge(iE));
  const Edge & E = *mesh().edge(iE);
  for (size_t iV = 0; iV < E.n_vertices(); iV++) {
    const Vertex & V = *E.vertex(iV);
    int VertID=V.global_index();
    int num=m_connectedinfo.vertex_component(VertID,iT);
    vE.segment(m_vsxgrad.localOffset(E, V), m_vsxgrad.numLocalDofsVertex(VertID))
      = vh.segment(globalOffset(V)+num*m_vsxgrad.numLocalDofsVertex(VertID), m_vsxgrad.numLocalDofsVertex(VertID));
  }

    int num=m_connectedinfo.edge_component(iE,iT);
    vE.segment(m_vsxgrad.localOffset(E), m_vsxgrad.numLocalDofsEdge(iE))
      = vh.segment(globalOffset(E)+num*m_vsxgrad.numLocalDofsEdge(iE), m_vsxgrad.numLocalDofsEdge(iE));
  
  return vE;
}

//------------------------------------------------------------------------------

Eigen::VectorXd FractureDOFSpace::restrictFace(size_t iF, const Eigen::VectorXd & vh, size_t iT) const
{
  Eigen::VectorXd vF = Eigen::VectorXd::Zero(m_vsxgrad.dimensionFace(iF));
  const Face & F = *(mesh().face(iF));

  for (size_t iV = 0; iV < F.n_vertices(); iV++) {
    const Vertex & V = *F.vertex(iV);
    int VertID=V.global_index();
    int num=m_connectedinfo.vertex_component(VertID,iT);
    vF.segment(m_vsxgrad.localOffset(F, V), m_vsxgrad.numLocalDofsVertex(VertID))
      = vh.segment(globalOffset(V)+num*m_vsxgrad.numLocalDofsVertex(VertID), m_vsxgrad.numLocalDofsVertex(VertID));
  } // for iV

  for (size_t iE = 0; iE < F.n_edges(); iE++) {
    const Edge & E = *F.edge(iE);      
    int edgeID=E.global_index();
    int num=m_connectedinfo.edge_component(edgeID,iT);
    vF.segment(m_vsxgrad.localOffset(F, E), m_vsxgrad.numLocalDofsEdge(edgeID))
      = vh.segment(globalOffset(E)+num*m_vsxgrad.numLocalDofsEdge(edgeID), m_vsxgrad.numLocalDofsEdge(edgeID));
  } // for iE

  int num=m_connectedinfo.face_component(iF,iT);
  vF.segment(m_vsxgrad.localOffset(F), m_vsxgrad.numLocalDofsFace(iF))
      = vh.segment(globalOffset(F)+num*m_vsxgrad.numLocalDofsFace(iF), m_vsxgrad.numLocalDofsFace(iF));
  
  return vF;
}


//------------------------------------------------------------------------------

Eigen::VectorXd FractureDOFSpace::restrictCell(size_t iT, const Eigen::VectorXd & vh) const
{
  // Eigen::VectorXd vT = Eigen::VectorXd::Zero(m_vsxgrad.dimensionCell(iT));
  Eigen::VectorXd vT = Eigen::VectorXd::Zero(m_vsxgrad.dimensionCell(iT)+3*m_connectedinfo.Cell_fractures.at(m_vsxgrad.mesh().cell(iT)).first);
  const Cell & T = *(mesh().cell(iT));

  for (size_t iV = 0; iV < T.n_vertices(); iV++) {
    const Vertex & V = *T.vertex(iV);
    int VertID=V.global_index();
    int num=m_connectedinfo.vertex_component(VertID,iT);
    vT.segment(m_vsxgrad.localOffset(T, V), m_vsxgrad.numLocalDofsVertex(VertID))
      = vh.segment(globalOffset(V)+num*m_vsxgrad.numLocalDofsVertex(VertID), m_vsxgrad.numLocalDofsVertex(VertID));
  } // for iV

  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);      
    int edgeID=E.global_index();
    int num=m_connectedinfo.edge_component(edgeID,iT);
    vT.segment(m_vsxgrad.localOffset(T, E), m_vsxgrad.numLocalDofsEdge(edgeID))
      = vh.segment(globalOffset(E)+num*m_vsxgrad.numLocalDofsEdge(edgeID), m_vsxgrad.numLocalDofsEdge(edgeID));
  } // for iE

  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);      
    int faceID=F.global_index();
    int num=m_connectedinfo.face_component(faceID,iT);
    vT.segment(m_vsxgrad.localOffset(T, F), m_vsxgrad.numLocalDofsFace(faceID))
      = vh.segment(globalOffset(F)+num*m_vsxgrad.numLocalDofsFace(faceID), m_vsxgrad.numLocalDofsFace(faceID));
  } // for iF

  vT.segment(m_vsxgrad.localOffset(T), m_vsxgrad.numLocalDofsCell(T))
    = vh.segment(globalOffset(T), m_vsxgrad.numLocalDofsCell(T));


  // Lagrange Fracture DOFs
  size_t offset_frac = dimension();
  for (size_t i = 0; i < m_connectedinfo.Cell_fractures.at(m_vsxgrad.mesh().cell(iT)).first; ++i){
    const Face *pF=m_connectedinfo.Cell_fractures.at(m_vsxgrad.mesh().cell(iT)).second[i];
    const Face & F = *pF;
    size_t global_offset_frac = 3*m_connectedinfo.findFaceIndex(F);
    for (size_t j = 0; j < 3; ++j){
      vT(m_vsxgrad.dimensionCell(iT)+3*i+j)=vh(offset_frac + global_offset_frac + j);
    }
  }


  return vT;
}

//------------------------------------------------------------------------------

Eigen::VectorXd FractureDOFSpace::restrictCellnew(size_t iT, const Eigen::VectorXd & vh) const
{
  Eigen::VectorXd vT = Eigen::VectorXd::Zero(m_vsxgrad.dimensionCell(iT));
  // Eigen::VectorXd vT = Eigen::VectorXd::Zero(m_vsxgrad.dimensionCell(iT)+3*m_connectedinfo.Cell_fractures.at(m_vsxgrad.mesh().cell(iT)).first);
  const Cell & T = *(mesh().cell(iT));

  for (size_t iV = 0; iV < T.n_vertices(); iV++) {
    const Vertex & V = *T.vertex(iV);
    int VertID=V.global_index();
    int num=m_connectedinfo.vertex_component(VertID,iT);
    vT.segment(m_vsxgrad.localOffset(T, V), m_vsxgrad.numLocalDofsVertex(VertID))
      = vh.segment(globalOffset(V)+num*m_vsxgrad.numLocalDofsVertex(VertID), m_vsxgrad.numLocalDofsVertex(VertID));
  } // for iV

  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);      
    int edgeID=E.global_index();
    int num=m_connectedinfo.edge_component(edgeID,iT);
    vT.segment(m_vsxgrad.localOffset(T, E), m_vsxgrad.numLocalDofsEdge(edgeID))
      = vh.segment(globalOffset(E)+num*m_vsxgrad.numLocalDofsEdge(edgeID), m_vsxgrad.numLocalDofsEdge(edgeID));
  } // for iE

  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);      
    int faceID=F.global_index();
    int num=m_connectedinfo.face_component(faceID,iT);
    vT.segment(m_vsxgrad.localOffset(T, F), m_vsxgrad.numLocalDofsFace(faceID))
      = vh.segment(globalOffset(F)+num*m_vsxgrad.numLocalDofsFace(faceID), m_vsxgrad.numLocalDofsFace(faceID));
  } // for iF

  vT.segment(m_vsxgrad.localOffset(T), m_vsxgrad.numLocalDofsCell(T))
    = vh.segment(globalOffset(T), m_vsxgrad.numLocalDofsCell(T));


  return vT;
}

std::vector<double> FractureDOFSpace::computeH1norms(const std::vector<Eigen::VectorXd> & list_dofs, const std::string typegrad) const
{
  size_t nb_vectors = list_dofs.size();
  std::vector<Eigen::VectorXd> local_sqnorms(nb_vectors, Eigen::VectorXd::Zero(mesh().n_cells()));

  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &list_dofs, &local_sqnorms, &nb_vectors, &typegrad](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *(mesh().cell(iT));

        // Full or symmetric gradient
        Eigen::MatrixXd GT = Eigen::MatrixXd::Zero(m_vsxgrad.cellBases(iT).Polyk3x3->dimension(), m_vsxgrad.dimensionCell(iT));
        if (typegrad == "symmetric"){
          GT = m_vsxgrad.cellSymGradient(iT);
        }else{
          GT = m_vsxgrad.cellGradient(iT);
        }
        // Contribution of cell T
	      MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_vsxgrad.degree());
        Eigen::MatrixXd mass_Pk3x3_T = GramMatrix(T, *m_vsxgrad.cellBases(iT).Polyk3x3, int_mono_2k);
        Eigen::MatrixXd ST = m_vsxgrad.computeStabilisation(iT);
        for (size_t i=0; i<nb_vectors; i++){
          Eigen::VectorXd uT = restrictnew(T, list_dofs[i]);
          local_sqnorms[i](iT) += uT.transpose() * (GT.transpose() * mass_Pk3x3_T * GT + ST) * uT;
        }
        
      }
    };
  parallel_for(mesh().n_cells(), compute_local_squarednorms, m_use_threads);
  
  // Vector of outputs
  std::vector<double> list_norms(nb_vectors, 0.);
  for (size_t i=0; i<nb_vectors; i++){
    list_norms[i] = std::sqrt(std::abs(local_sqnorms[i].sum()));
  }
  
  return list_norms;

}

std::vector<double> FractureDOFSpace::computeL2normdisp(const FunctionType & exact_solution, const double & lambda, const Eigen::VectorXd & discrete_solution) const
{
  std::vector<Eigen::VectorXd> local_sqnorms(2,Eigen::VectorXd::Zero(mesh().n_cells()));
  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &exact_solution, &local_sqnorms, &discrete_solution, &lambda](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *(mesh().cell(iT));
        int dqr_cell = 8;
        QuadratureRule quad_dqrbc_T = generate_quadrature_rule(T, dqr_cell);
        // Storing exact solution for the refernece element as u_local
        std::function<Eigen::Vector3d(const Eigen::Vector3d &)> u_local = 
          [&exact_solution,&T,&lambda](const Eigen::Vector3d & x) -> Eigen::Vector3d {
            return exact_solution(x,&T,lambda);
          };
        // Calling the reference potential reconstruction
        const Eigen::VectorXd dofs = restrictnew(T,discrete_solution);
        std::function<Eigen::Vector3d(const Eigen::Vector3d&)> potential_function = 
          [this, iT, &dofs](const Eigen::Vector3d& x) {
              return evaluateCellPotential(iT, x, dofs);
          };       
        // Creating error function
        std::function<Eigen::Vector3d(const Eigen::Vector3d&)> error_function =
            [&u_local, &potential_function](const Eigen::Vector3d& x) -> Eigen::Vector3d {
                return u_local(x) - potential_function(x);
            };
          
        for (size_t iqn = 0; iqn < quad_dqrbc_T.size(); iqn++){
          double qr_weight = quad_dqrbc_T[iqn].w;
          local_sqnorms[0](iT) += qr_weight * error_function(quad_dqrbc_T[iqn].vector()).squaredNorm();
          local_sqnorms[1](iT) += qr_weight * u_local(quad_dqrbc_T[iqn].vector()).squaredNorm();
        }
      }
    };
    parallel_for(mesh().n_cells(), compute_local_squarednorms, m_use_threads); 

  // Vector of outputs
  std::vector<double> list_norms(2, 0.);
  for (size_t i=0; i<2; i++){
    list_norms[i] = std::sqrt(std::abs(local_sqnorms[i].sum()));
  }
  
  return list_norms;
}

std::vector<double> FractureDOFSpace::computeL2normdispgrad(const GradientDisplacementType & exact_gradient, const double & lambda, const Eigen::VectorXd & discrete_solution) const
{
  std::vector<Eigen::VectorXd> local_sqnorms(2,Eigen::VectorXd::Zero(mesh().n_cells()));
  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &exact_gradient, &local_sqnorms, &discrete_solution, &lambda](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *(mesh().cell(iT));
        int dqr_cell = 8;
        QuadratureRule quad_dqrbc_T = generate_quadrature_rule(T, dqr_cell);
        // Storing exact solution for the refernece element as u_local
        std::function<Eigen::Matrix3d(const Eigen::Vector3d &)> u_grad = 
          [&exact_gradient,&T, &lambda](const Eigen::Vector3d & x) -> Eigen::Matrix3d {
            return exact_gradient(x,&T,lambda);
          };
        // Calling the reference potential reconstruction
        const Eigen::VectorXd dofs = restrictnew(T,discrete_solution);
        std::function<Eigen::Matrix3d(const Eigen::Vector3d&)> potential_gradient = 
          [this, iT, &dofs](const Eigen::Vector3d& x) {
              return evaluateCellGradient(iT, x, dofs);
          };       
        // Creating error function
        std::function<Eigen::Matrix3d(const Eigen::Vector3d&)> error_function =
            [&u_grad, &potential_gradient](const Eigen::Vector3d& x) -> Eigen::Matrix3d {
                return u_grad(x) - potential_gradient(x);
            };
          
        for (size_t iqn = 0; iqn < quad_dqrbc_T.size(); iqn++){
          double qr_weight = quad_dqrbc_T[iqn].w;
          local_sqnorms[0](iT) += qr_weight * error_function(quad_dqrbc_T[iqn].vector()).squaredNorm();
          local_sqnorms[1](iT) += qr_weight * u_grad(quad_dqrbc_T[iqn].vector()).squaredNorm();
        }
      }
    };
    parallel_for(mesh().n_cells(), compute_local_squarednorms, m_use_threads); 

  // Vector of outputs
  std::vector<double> list_norms(2, 0.);
  for (size_t i=0; i<2; i++){
    list_norms[i] = std::sqrt(std::abs(local_sqnorms[i].sum()));
  }
  
  return list_norms;
}

std::vector<double> FractureDOFSpace::computeL2normdispjump(const FunctionType & exact_solution, const double & lambda, const Eigen::VectorXd & discrete_solution) const
{
  std::vector<Eigen::VectorXd> local_sqnorms(2,Eigen::VectorXd::Zero(m_connectedinfo.fracture_face_indices_pointer.size()));
  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &exact_solution, &local_sqnorms, &discrete_solution, &lambda](size_t start, size_t end)->void
    {
      for (size_t i = start; i < end; i++){
        Face * F = m_connectedinfo.fracture_face_indices_pointer[i];
        size_t iF = F->global_index();
        int dqr_face=8;
        QuadratureRule quad_dqrbc_F = generate_quadrature_rule(*F, dqr_face);
        const Cell & cell_ptr_0 = *m_connectedinfo.face_final_connected_pointer.at(mesh().face(iF))[0][0];
        const Cell & cell_ptr_1 = *m_connectedinfo.face_final_connected_pointer.at(mesh().face(iF))[1][0];
        // Storing exact solution for the refernece element as u_local
        std::function<Eigen::Vector3d(const Eigen::Vector3d &)> u_jump = 
          [&exact_solution, &cell_ptr_0, &cell_ptr_1,&lambda](const Eigen::Vector3d & x) -> Eigen::Vector3d {
            return exact_solution(x,&cell_ptr_0,lambda)-exact_solution(x,&cell_ptr_1,lambda);
          };
        // Calling the reference potential reconstruction
        const Eigen::VectorXd dofs_0 = restrict(*F,discrete_solution,cell_ptr_0);
        const Eigen::VectorXd dofs_1 = restrict(*F,discrete_solution,cell_ptr_1);
        std::function<Eigen::Vector3d(const Eigen::Vector3d&)> potential_jump = 
          [this, iF, &dofs_0, &dofs_1](const Eigen::Vector3d& x) {
              return evaluateFacePotential(iF, x, dofs_0)-evaluateFacePotential(iF, x, dofs_1);
          };       
        // Creating error function
        std::function<Eigen::Vector3d(const Eigen::Vector3d&)> error_function =
            [&u_jump, &potential_jump](const Eigen::Vector3d& x) -> Eigen::Vector3d {
                return u_jump(x) - potential_jump(x);
            };          
        for (size_t iqn = 0; iqn < quad_dqrbc_F.size(); iqn++){
          double qr_weight = quad_dqrbc_F[iqn].w;
          local_sqnorms[0](i) += qr_weight * error_function(quad_dqrbc_F[iqn].vector()).squaredNorm();
          local_sqnorms[1](i) += qr_weight * u_jump(quad_dqrbc_F[iqn].vector()).squaredNorm();
        }
      }
    };
    parallel_for(m_connectedinfo.fracture_face_indices_pointer.size(), compute_local_squarednorms, m_use_threads); 
  // Vector of outputs
  std::vector<double> list_norms(2, 0.);
  for (size_t i=0; i<2; i++){
    list_norms[i] = std::sqrt(std::abs(local_sqnorms[i].sum()));
  }
  
  return list_norms;
}

std::vector<double> FractureDOFSpace::computeL2normtraction(const GradientDisplacementType1 & exact_traction, const double & lambda, const Eigen::VectorXd & discrete_traction) const
{
  std::vector<Eigen::VectorXd> local_sqnorms(2,Eigen::VectorXd::Zero(m_connectedinfo.fracture_face_indices_pointer.size()));
  std::function<void(size_t, size_t)> compute_local_squarednorms
  = [this, &exact_traction, &local_sqnorms, &discrete_traction, &lambda](size_t start, size_t end)->void
    {
      for (size_t i = start; i < end; i++){
        Face * F = m_connectedinfo.fracture_face_indices_pointer[i];
        double lambda_n=discrete_traction(dimension()+3*i);
        size_t iF = F->global_index();
        int dqr_face=8;
        QuadratureRule quad_dqrbc_F = generate_quadrature_rule(*F, dqr_face);
        const Cell & cell_ptr = *m_connectedinfo.face_final_connected_pointer.at(mesh().face(iF))[0][0];
        Eigen::Vector3d nF = F->normal();
        std::function<double(const Eigen::Vector3d&)> traction = 
              [&exact_traction, &nF, &cell_ptr,&lambda](const Eigen::Vector3d & x) -> double {
                double a = 1.0; // or whatever values are appropriate
                return scalar_product(exact_traction(a, lambda, x, &cell_ptr) * nF, nF);
              };
        for (size_t iqn = 0; iqn < quad_dqrbc_F.size(); iqn++){
          double qr_weight = quad_dqrbc_F[iqn].w;
          local_sqnorms[0](i) += qr_weight * pow(traction(quad_dqrbc_F[iqn].vector())+lambda_n,2);
          local_sqnorms[1](i) += qr_weight * pow(traction(quad_dqrbc_F[iqn].vector()),2);
        }
      }
    };
    parallel_for(m_connectedinfo.fracture_face_indices_pointer.size(), compute_local_squarednorms, m_use_threads); 
  // Vector of outputs
  std::vector<double> list_norms(2, 0.);
  for (size_t i=0; i<2; i++){
    list_norms[i] = std::sqrt(std::abs(local_sqnorms[i].sum()));
  }
  
  return list_norms;
}

std::vector<Eigen::VectorXd> FractureDOFSpace::computecontact(const FrictionType & friction, const Eigen::VectorXd & discrete_solution, const double & beta, const double & lambda) const
{
  std::vector<Eigen::VectorXd> contact(8,Eigen::VectorXd::Zero(m_connectedinfo.fracture_face_indices_pointer.size()));
  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &contact, &discrete_solution, &friction, &beta,&lambda](size_t start, size_t end)->void
    {
      for (size_t i = start; i < end; i++){
        Face * F = m_connectedinfo.fracture_face_indices_pointer[i];
        size_t iF = F->global_index();
        auto frame = get_frame(*F);
        Eigen::Vector3d xF= F->center_mass();
        const Cell & cell_ptr_0 = *m_connectedinfo.face_final_connected_pointer.at(mesh().face(iF))[0][0];
        int direction=cell_ptr_0.face_orientation(cell_ptr_0.index_face(F));
        // Calling the reference 
        Eigen::VectorXd potential_jump = direction*(discrete_solution.segment(globalOffset(*F),3)-discrete_solution.segment(globalOffset(*F)+3,3));  
        contact[0](i) = m_vsxgrad.sXgrad().faceBases(iF).Polykmo->function(0, xF) * potential_jump.dot(frame.normal);
        contact[1](i) = m_vsxgrad.sXgrad().faceBases(iF).Polykmo->function(0, xF) * potential_jump.dot(frame.tangent1);
        contact[2](i) = m_vsxgrad.sXgrad().faceBases(iF).Polykmo->function(0, xF) * potential_jump.dot(frame.tangent2);
        Eigen::VectorXd discrete_traction = discrete_solution.segment(dimension()+3*i,3);
        const double traction_tangential_norm = std::sqrt(std::pow(discrete_traction(1), 2) + std::pow(discrete_traction(2), 2));
        Eigen::Vector3d bnfrac = Eigen::Vector3d::Zero();
        for(size_t j = 0; j < 3; ++j){
          bnfrac(j) = discrete_traction(j) + beta * contact[j](i);
        }
        const double bnfrac_tangential_norm = std::sqrt(std::pow(bnfrac(1), 2) + std::pow(bnfrac(2), 2));
        if (bnfrac(0) > 0.0){
          contact[3](i) = 1.0;
        }else{
          contact[3](i) = 0.0;
        }

        if (bnfrac_tangential_norm < friction(F,lambda)){
          contact[4](i) = 0.0;
        }else{
          contact[4](i) = 1.0;
        }

        contact[5](i) = contact[3](i) + 2 * contact[4](i);

        contact[6](i) = traction_tangential_norm;

        contact[7](i) = discrete_traction(0);

      }
    };
    parallel_for(m_connectedinfo.fracture_face_indices_pointer.size(), compute_local_squarednorms, m_use_threads); 
  return contact;
}


Eigen::VectorXd FractureDOFSpace::computenormaltraction(const std::vector<Face *> & fracture_pointer, const Eigen::VectorXd & discrete_solution) const
{
  Eigen::VectorXd contact(fracture_pointer.size());
  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &contact, &discrete_solution, &fracture_pointer](size_t start, size_t end)->void
    {
      for (size_t i = start; i < end; i++){
        Face * F = fracture_pointer[i];
        int j = m_connectedinfo.findFaceIndex(*F);
        contact(i) = discrete_solution(dimension()+3*j);
      }
    };
    parallel_for(fracture_pointer.size(), compute_local_squarednorms, m_use_threads); 
  return contact;
}




std::vector<double>  FractureDOFSpace::reconstruct_norm(const Eigen::VectorXd & discrete_solution) const
{
  std::vector<Eigen::VectorXd> contact(16,Eigen::VectorXd::Zero(m_connectedinfo.fracture_face_indices_pointer.size()));
  std::vector<double> counts(4, 0.0);
  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &contact, &discrete_solution, &counts](size_t start, size_t end)->void
    {
      std::vector<double> local_counts(4, 0.0);
      for (size_t i = start; i < end; i++){
        Face * F = m_connectedinfo.fracture_face_indices_pointer[i];
        size_t iF = F->global_index();
        auto frame = get_frame(*F);
        Eigen::Vector3d xF= F->center_mass();
        const Cell & T0 = *m_connectedinfo.face_final_connected_pointer.at(mesh().face(iF))[0][0];
        int direction=T0.face_orientation(T0.index_face(F));
        Eigen::VectorXd u_jump = direction * (discrete_solution.segment(globalOffset(*F),3)-discrete_solution.segment(globalOffset(*F)+3,3));
        double u_jump_normal = F->measure()*m_vsxgrad.sXgrad().faceBases(iF).Polykmo->function(0, xF) * u_jump.dot(frame.normal);;
        double u_jump_tangent1 = m_vsxgrad.sXgrad().faceBases(iF).Polykmo->function(0, xF) * u_jump.dot(frame.tangent1);;
        double u_jump_tangent2 = m_vsxgrad.sXgrad().faceBases(iF).Polykmo->function(0, xF) * u_jump.dot(frame.tangent2);;

        double lambda_prev_normal = F->measure() * discrete_solution(dimension() + 3 * i);
        double lambda_prev_tangent1 = discrete_solution(dimension() + 3 * i + 1);
        double lambda_prev_tangent2 = discrete_solution(dimension() + 3 * i + 2);
        const double a1 = F->center_mass().x();
        const double a2 = F->center_mass().y();
        const double a3 = F->center_mass().z();
        if (abs(25*a1-2*a2+10*a3-11.4)<1e-6 && a1>0.2 &&a1<0.4 &&a2>0.3 && a2<0.8
          ) {     
            contact[0](i) = u_jump_normal;
            contact[4](i) = F->measure() * std::sqrt(std::pow(u_jump_tangent1,2)+std::pow(u_jump_tangent2,2));
            contact[8](i) = lambda_prev_normal;
            contact[12](i) = F->measure() * std::sqrt(std::pow(lambda_prev_tangent1,2)+std::pow(lambda_prev_tangent2,2));
            local_counts[0] += F->measure();
        }else if (abs(250*a1+15*a2-75*a3-152)<1e-6 && a1>0.65 &&a1<0.8 &&a2>0.3 && a2<0.8
          ) {     
            contact[1](i) = u_jump_normal;
            contact[5](i) = F->measure() * std::sqrt(std::pow(u_jump_tangent1,2)+std::pow(u_jump_tangent2,2));
            contact[9](i) = lambda_prev_normal;
            contact[13](i) = F->measure() * std::sqrt(std::pow(lambda_prev_tangent1,2)+std::pow(lambda_prev_tangent2,2));
            local_counts[1] += F->measure();
        }else if (abs((-0.1)*a1+(0.08)*a2+(0.4)*a3-0.174)<1e-6 && a1>0.1 &&a1<0.9 &&a2>0.3 && a2<0.8
          ) {     
            contact[2](i) = u_jump_normal;
            contact[6](i) = F->measure() * std::sqrt(std::pow(u_jump_tangent1,2)+std::pow(u_jump_tangent2,2));
            contact[10](i) = lambda_prev_normal;
            contact[14](i) = F->measure() * std::sqrt(std::pow(lambda_prev_tangent1,2)+std::pow(lambda_prev_tangent2,2));
            local_counts[2] += F->measure();
        }else if (abs(a3-0.10)<1e-8 && a1>0.7 &&a1<0.9 &&a2>0.15 && a2<0.35
          ) {     
            contact[3](i) = u_jump_normal;
            contact[7](i) = F->measure() * std::sqrt(std::pow(u_jump_tangent1,2)+std::pow(u_jump_tangent2,2));
            contact[11](i) = lambda_prev_normal;
            contact[15](i) = F->measure() * std::sqrt(std::pow(lambda_prev_tangent1,2)+std::pow(lambda_prev_tangent2,2));
            local_counts[3] += F->measure();
        } 

      }

      {
        for (int j = 0; j < 4; j++) {
            counts[j] += local_counts[j];
        }
      }
    };
    parallel_for(m_connectedinfo.fracture_face_indices_pointer.size(), compute_local_squarednorms, m_use_threads); 
  // Vector of outputs
  std::vector<double> list_norms(16, 0.0);
  for (size_t i = 0; i < 4; i++){
    list_norms[4*i] = (contact[4*i].sum()) / counts[0];
    list_norms[4*i+1] = (contact[4*i+1].sum()) / counts[1];
    list_norms[4*i+2] = (contact[4*i+2].sum()) / counts[2];
    list_norms[4*i+3] = (contact[4*i+3].sum()) / counts[3];
  }
  return list_norms;
}



std::vector<size_t> FractureDOFSpace::globalDOFIndices(const Cell & T) const
{
  int iT=T.global_index();
  std::vector<size_t> I(m_vsxgrad.dimensionCell(T)+3*m_connectedinfo.Cell_fractures.at(m_vsxgrad.mesh().cell(iT)).first);

  size_t dof_index = 0;
  
  // Vertex DOFs
  for (size_t iV = 0; iV < T.n_vertices(); iV++) {
    const Vertex & V = *T.vertex(iV);
    int VertID=V.global_index();
    int num=m_connectedinfo.vertex_component(VertID,iT);
    size_t offset_V = globalOffset(V)+num*m_vsxgrad.numLocalDofsVertex(VertID);
    for (size_t i = 0; i < m_vsxgrad.numLocalDofsVertex(V); i++, dof_index++) {
      I[dof_index] = offset_V + i;
    } // for i
  } // for iV

  // Edge DOFs
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    int EdgeID=E.global_index();
    int num=m_connectedinfo.edge_component(EdgeID,iT);
    size_t offset_E = globalOffset(E)+num*m_vsxgrad.numLocalDofsEdge(EdgeID);
    for (size_t i = 0; i < m_vsxgrad.numLocalDofsEdge(E); i++, dof_index++) {
      I[dof_index] = offset_E + i;
    } // for i
  } // for iE

  // Face DOFs
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    int FaceID=F.global_index();
    int num=m_connectedinfo.face_component(FaceID,iT);
    size_t offset_F = globalOffset(F)+num*m_vsxgrad.numLocalDofsFace(FaceID);
    for (size_t i = 0; i < m_vsxgrad.numLocalDofsFace(F); i++, dof_index++) {
      I[dof_index] = offset_F + i;
    } // for i
  } // for iF

  // Cell DOFs
  size_t offset_T = globalOffset(T);
  for (size_t i = 0; i < m_vsxgrad.numLocalDofsCell(T); i++, dof_index++) {
    I[dof_index] = offset_T + i;
  } // for i

  // Lagrange Fracture DOFs
  size_t offset_frac = dimension();
  for (size_t i = 0; i < m_connectedinfo.Cell_fractures.at(m_vsxgrad.mesh().cell(iT)).first; ++i){
    const Face *pF=m_connectedinfo.Cell_fractures.at(m_vsxgrad.mesh().cell(iT)).second[i];
    const Face & F = *pF;
    size_t global_offset_frac = 3*m_connectedinfo.findFaceIndex(F);
    for (size_t j = 0; j < 3; j++, dof_index++){
      I[dof_index]=offset_frac + global_offset_frac + j;
    }
  }
  assert( m_vsxgrad.dimensionCell(T) + 3*m_connectedinfo.Cell_fractures.at(m_vsxgrad.mesh().cell(iT)).first == dof_index );
  
  return I;
}

//------------------------------------------------------------------------------

std::vector<size_t> FractureDOFSpace::globalDOFIndices(const Face & F, const Cell & T) const
{
  std::vector<size_t> I(m_vsxgrad.dimensionFace(F));
  int iF=F.global_index();
  int iT=T.global_index();

  size_t dof_index = 0;
  
  // Vertex DOFs
  for (size_t iV = 0; iV < F.n_vertices(); iV++) {
    const Vertex & V = *F.vertex(iV);
    int VertID=V.global_index();
    int num=m_connectedinfo.vertex_component(VertID,iT);
    size_t offset_V = globalOffset(V)+num*m_vsxgrad.numLocalDofsVertex(VertID);
    for (size_t i = 0; i < m_vsxgrad.numLocalDofsVertex(V); i++, dof_index++) {
      I[dof_index] = offset_V + i;
    } // for i
  } // for iV

  // Edge DOFs
  for (size_t iE = 0; iE < F.n_edges(); iE++) {
    const Edge & E = *F.edge(iE);
    int EdgeID=E.global_index();
    int num=m_connectedinfo.edge_component(EdgeID,iT);
    size_t offset_E = globalOffset(E)+num*m_vsxgrad.numLocalDofsEdge(EdgeID);
    for (size_t i = 0; i < m_vsxgrad.numLocalDofsEdge(E); i++, dof_index++) {
      I[dof_index] = offset_E + i;
    } // for i
  } // for iE

  // Face DOFs
  int num=m_connectedinfo.face_component(iF,iT);
  size_t offset_F = globalOffset(F)+num*m_vsxgrad.numLocalDofsFace(iF);
  for (size_t i = 0; i < m_vsxgrad.numLocalDofsFace(F); i++, dof_index++) {
    I[dof_index] = offset_F + i;
  } // for iF

  assert( m_vsxgrad.dimensionFace(F) == dof_index );
  
  return I;
}

std::vector<size_t> FractureDOFSpace::globalDOFIndices(const Edge & E, const Cell & T) const
{
  std::vector<size_t> I(m_vsxgrad.dimensionEdge(E));
  int iE=E.global_index();
  int iT=T.global_index();

  size_t dof_index = 0;
  
  // Vertex DOFs
  for (size_t iV = 0; iV < E.n_vertices(); iV++) {
    const Vertex & V = *E.vertex(iV);
    int VertID=V.global_index();
    int num=m_connectedinfo.vertex_component(VertID,iT);
    size_t offset_V = globalOffset(V)+num*m_vsxgrad.numLocalDofsVertex(VertID);
    for (size_t i = 0; i < m_vsxgrad.numLocalDofsVertex(V); i++, dof_index++) {
      I[dof_index] = offset_V + i;
    } // for i
  } // for iV

  // Edge DOFs
  int num=m_connectedinfo.edge_component(iE,iT);
  size_t offset_E = globalOffset(E)+num*m_vsxgrad.numLocalDofsEdge(iE);
  for (size_t i = 0; i < m_vsxgrad.numLocalDofsEdge(E); i++, dof_index++) {
    I[dof_index] = offset_E + i;
  } // for iE

  assert( m_vsxgrad.dimensionEdge(E) == dof_index );
  
  return I;
}

std::vector<size_t> FractureDOFSpace::globalDOFIndices(const Vertex & V, const Cell & T) const
{
  std::vector<size_t> I(m_vsxgrad.dimensionVertex(V));
  int iV=V.global_index();
  int iT=T.global_index();

  size_t dof_index = 0;
  
  // Vertex DOFs
  int num=m_connectedinfo.vertex_component(iV,iT);
  size_t offset_V = globalOffset(V)+num*m_vsxgrad.numLocalDofsVertex(iV);
  for (size_t i = 0; i < m_vsxgrad.numLocalDofsVertex(V); i++, dof_index++) {
    I[dof_index] = offset_V + i;
  } // for iE

  assert( m_vsxgrad.dimensionVertex(V) == dof_index );
  
  return I;
}
