
#include <vgrad.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_face.hpp>
#include <GMpoly_edge.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

VGrad::VGrad(const VEMCore & vem_core, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(vem_core.mesh(),
	     1,
	     PolynomialSpaceDimension<Edge>::Poly(vem_core.degree() - 1),
	     PolynomialSpaceDimension<Face>::RolyCompl(vem_core.degree() + 1),
	     PolynomialSpaceDimension<Cell>::RolyCompl(vem_core.degree())
	     ),
    m_vem_core(vem_core),
    m_use_threads(use_threads),
    m_output(output),
    m_edge_operators(vem_core.mesh().n_edges()),
    m_face_operators(vem_core.mesh().n_faces()),
    m_cell_operators(vem_core.mesh().n_cells())
{
  m_output << "[VGrad] Initializing" << std::endl;
  if (use_threads) {
    m_output << "[VGrad] Parallel execution" << std::endl;
  } else {
    m_output << "[VGrad] Sequential execution" << std::endl;
  }
  
  // Construct edge operators
  std::function<void(size_t, size_t)> construct_all_edge_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iE = start; iE < end; iE++) {
          m_edge_operators[iE].reset( new LocalOperators(_compute_edge_operators(iE)) );
        } // for iE
      };

  m_output << "[VGrad] Constructing edge operators" << std::endl;
  parallel_for(mesh().n_edges(), construct_all_edge_operators, use_threads);

  // Construct face operators
  std::function<void(size_t, size_t)> construct_all_face_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iF = start; iF < end; iF++) {
          m_face_operators[iF].reset( new LocalOperators(_compute_face_operators(iF)) );
        } // for iF
      };

  m_output << "[VGrad] Constructing face operators" << std::endl;
  parallel_for(mesh().n_faces(), construct_all_face_operators, use_threads);

  // Construct cell gradients and potentials
  std::function<void(size_t, size_t)> construct_all_cell_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_cell_operators[iT].reset( new LocalOperators(_compute_cell_operators(iT)) );
        } // for iT
      };

  m_output << "[VGrad] Constructing cell operators" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_operators, use_threads);
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd VGrad::interpolate(const FunctionType & q, const GradientType & grad_q, const int doe_cell, const int doe_face, const int doe_edge) const
{
  Eigen::VectorXd qh = Eigen::VectorXd::Zero(dimension());

  // Degrees of quadrature rules
  size_t dqr_cell = (doe_cell >= 0 ? doe_cell : 2 * degree() + 3);
  size_t dqr_face = (doe_face >= 0 ? doe_face : 2 * degree() + 3);
  size_t dqr_edge = (doe_edge >= 0 ? doe_edge : 2 * degree() + 3);
  
  // Interpolate at vertices
  std::function<void(size_t, size_t)> interpolate_vertices
    = [this, &qh, q](size_t start, size_t end)->void
      {
        for (size_t iV = start; iV < end; iV++) {
          qh(iV) = q(mesh().vertex(iV)->coords());
        } // for iV
      };
  parallel_for(mesh().n_vertices(), interpolate_vertices, m_use_threads);

  if (degree() > 0) {
    
    // Interpolate at edges
    std::function<void(size_t, size_t)> interpolate_edges
      = [this, &qh, q, &dqr_edge](size_t start, size_t end)->void
        {
          for (size_t iE = start; iE < end; iE++) {
            const Edge & E = *mesh().edge(iE);
            QuadratureRule quad_dqr_E = generate_quadrature_rule(E, dqr_edge);
            auto basis_Pkmo_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polykmo, quad_dqr_E);
            qh.segment(globalOffset(E), PolynomialSpaceDimension<Edge>::Poly(degree() - 1)) 
              = l2_projection(q, *edgeBases(iE).Polykmo, quad_dqr_E, basis_Pkmo_E_quad);
          } // for iE
        };
    parallel_for(mesh().n_edges(), interpolate_edges, m_use_threads);
  } // if degree() > 0 
    
  // Interpolate at faces
  std::function<void(size_t, size_t)> interpolate_faces
    = [this, &qh, grad_q, &dqr_face](size_t start, size_t end)->void
      {
        for (size_t iF = start; iF < end; iF++) {
          const Face & F = *mesh().face(iF);
          QuadratureRule quad_dqr_F = generate_quadrature_rule(F, dqr_face);
          auto basis_Rckpo_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).RolyComplkpo, quad_dqr_F);
          qh.segment(globalOffset(F), PolynomialSpaceDimension<Face>::RolyCompl(degree() + 1)) 
            = l2_projection(grad_q, *faceBases(iF).RolyComplkpo, quad_dqr_F, basis_Rckpo_F_quad);
        } // for iF
      };
  parallel_for(mesh().n_faces(), interpolate_faces, m_use_threads);

  if (degree() > 0) {
    // Interpolate at cells
    std::function<void(size_t, size_t)> interpolate_cells
      = [this, &qh, grad_q, &dqr_cell](size_t start, size_t end)->void
        {
          for (size_t iT = start; iT < end; iT++) {
            const Cell & T = *mesh().cell(iT);
            QuadratureRule quad_dqr_T = generate_quadrature_rule(T, dqr_cell);
            auto basis_Rck_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).RolyComplk, quad_dqr_T);
            MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2);
            qh.segment(globalOffset(T), PolynomialSpaceDimension<Cell>::RolyCompl(degree())) 
              = l2_projection(grad_q, *cellBases(iT).RolyComplk, quad_dqr_T, basis_Rck_T_quad, GramMatrix(T, *cellBases(iT).RolyComplk, int_mono_2kp2));
          } // for iT
        };
    parallel_for(mesh().n_cells(), interpolate_cells, m_use_threads);

  } // if degree() > 0 

  return qh;
}

//------------------------------------------------------------------------------
// Gradient and function projection
//------------------------------------------------------------------------------

VGrad::LocalOperators VGrad::_compute_edge_operators(size_t iE)
{
  const Edge & E = *mesh().edge(iE);
  
  //------------------------------------------------------------------------------
  // Gradient: projection and DOFs in Vcurl (they're the same)
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  
  MonomialEdgeIntegralsType int_mono_2kp2_E = IntegrateEdgeMonomials(E, 2*degree()+2);
  Eigen::MatrixXd MGE = GramMatrix(E, *edgeBases(iE).Polyk, int_mono_2kp2_E);

  //------------------------------------------------------------------------------
  // Right-hand side matrix
  
  Eigen::MatrixXd BGE
    = Eigen::MatrixXd::Zero(edgeBases(iE).Polyk->dimension(), dimensionEdge(iE));
  for (size_t i = 0; i < edgeBases(iE).Polyk->dimension(); i++) {
    BGE(i, 0) = -edgeBases(iE).Polyk->function(i, mesh().edge(iE)->vertex(0)->coords());
    BGE(i, 1) = edgeBases(iE).Polyk->function(i, mesh().edge(iE)->vertex(1)->coords());
  } // for i
  
  if (degree() > 0) {
    GradientBasis<VEMCore::PolyBasisEdgeType> grad_Pk_E(*edgeBases(iE).Polyk);
    BGE.rightCols(PolynomialSpaceDimension<Edge>::Poly(degree() - 1)) 
          = -GramMatrix(E, grad_Pk_E, *edgeBases(iE).Polykmo, int_mono_2kp2_E);
  }
  
  //------------------------------------------------------------------------------
  // Projection of function in \f$P^{k+1}\f$
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BPE
    = Eigen::MatrixXd::Zero(PolynomialSpaceDimension<Edge>::Poly(degree()) + 1, dimensionEdge(iE));

  // Enforce the gradient of the potential reconstruction
  BPE.topRows(PolynomialSpaceDimension<Edge>::Poly(degree())) = BGE;

  // Enforce the average value of the potential reconstruction
  if (degree() == 0) {
    // We set the average equal to the mean of vertex values
    BPE.bottomRows(1)(0, 0) = 0.5 * E.measure();
    BPE.bottomRows(1)(0, 1) = 0.5 * E.measure();
  } else {
    QuadratureRule quad_kmo_E = generate_quadrature_rule(E, degree() - 1);
    auto basis_Pkmo_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polykmo, quad_kmo_E);
    
    // We set the average value of the potential equal to the average of the edge unknown
    for (size_t i = 0; i < PolynomialSpaceDimension<Edge>::Poly(degree() - 1); i++) {
      for (size_t iqn = 0; iqn < quad_kmo_E.size(); iqn++) {
        BPE.bottomRows(1)(0, 2 + i) += quad_kmo_E[iqn].w * basis_Pkmo_E_quad[i][iqn];
      } // for iqn
    } // for i
  }
  
  //------------------------------------------------------------------------------
  // Left-hand side matrix
  
  Eigen::MatrixXd MPE
    = Eigen::MatrixXd::Zero(PolynomialSpaceDimension<Edge>::Poly(degree()) + 1, PolynomialSpaceDimension<Edge>::Poly(degree() + 1));

  GradientBasis<VEMCore::PolyBasisEdgeType> grad_Pkpo_E(*edgeBases(iE).Polykpo);				    
  MPE.topRows(PolynomialSpaceDimension<Edge>::Poly(degree()))
        = GramMatrix(E, *edgeBases(iE).Polyk, grad_Pkpo_E, int_mono_2kp2_E);

  MonomialScalarBasisEdge basis_P0_E(E, 0);
  MPE.bottomRows(1) = GramMatrix(E, basis_P0_E, *edgeBases(iE).Polykpo, int_mono_2kp2_E);

  Eigen::MatrixXd proj_gradient = MGE.ldlt().solve(BGE);
  
  return LocalOperators(proj_gradient, proj_gradient, MPE.partialPivLu().solve(BPE));
}

//------------------------------------------------------------------------------

VGrad::LocalOperators VGrad::_compute_face_operators(size_t iF)
{
  const Face & F = *mesh().face(iF);

  MonomialFaceIntegralsType int_mono_2kp2_F = IntegrateFaceMonomials(F, 2*degree()+2);
  
  //------------------------------------------------------------------------------
  // Projection of function on P^k
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  
  DivergenceBasis<VEMCore::RolyComplBasisFaceType> div_Rckpo_F(*faceBases(iF).RolyComplkpo);
  Eigen::MatrixXd MPF = GramMatrix(F, div_Rckpo_F, *faceBases(iF).Polyk, int_mono_2kp2_F);

  //------------------------------------------------------------------------------
  // Right-hand side matrix
  Eigen::MatrixXd BPF = Eigen::MatrixXd::Zero(faceBases(iF).Polyk->dimension(), dimensionFace(iF));

  // Face contribution: we need integrals up to 2k+3 here because Polyk2 is a restriction of a basis of degree k+1
  BPF.rightCols(faceBases(iF).RolyComplkpo->dimension()) = -GramMatrix(F, *faceBases(iF).RolyComplkpo, int_mono_2kp2_F);
  
  // Boundary contribution
  for (size_t iE = 0; iE < F.n_edges(); iE++) {
    const Edge & E = *F.edge(iE);
    
    QuadratureRule quad_2kp2_E = generate_quadrature_rule(E, 2 * (degree() + 1));
    auto basis_Rckpo_nFE_E_quad
      = scalar_product(evaluate_quad<Function>::compute(*faceBases(iF).RolyComplkpo, quad_2kp2_E), F.edge_normal(iE));

    auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polykpo, quad_2kp2_E);
    Eigen::MatrixXd PE = extendOperator(F, E, edgeOperators(E).proj_function);

    BPF += F.edge_orientation(iE) * compute_gram_matrix(basis_Rckpo_nFE_E_quad, basis_Pkpo_E_quad, quad_2kp2_E) * PE;
  } // for iE
  
  Eigen::MatrixXd proj_function = MPF.partialPivLu().solve(BPF);
  
  //------------------------------------------------------------------------------
  // Projection of gradient on P^k
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix

  Eigen::MatrixXd MGF = GramMatrix(F, *faceBases(iF).Polyk2, int_mono_2kp2_F);

  //------------------------------------------------------------------------------
  // Right-hand side matrix
  
  Eigen::MatrixXd BGF
    = Eigen::MatrixXd::Zero(faceBases(iF).Polyk2->dimension(), dimensionFace(iF));

  // Boundary contribution
  for (size_t iE = 0; iE < F.n_edges(); iE++) {
    const Edge & E = *F.edge(iE);
    
    QuadratureRule quad_2kp2_E = generate_quadrature_rule(E, 2 * (degree() + 1));
    auto basis_Pk2_nFE_E_quad
      = scalar_product(evaluate_quad<Function>::compute(*faceBases(iF).Polyk2, quad_2kp2_E), F.edge_normal(iE));
    auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polykpo, quad_2kp2_E);
    Eigen::MatrixXd PE = extendOperator(F, E, edgeOperators(E).proj_function);
    BGF += F.edge_orientation(iE) * compute_gram_matrix(basis_Pk2_nFE_E_quad, basis_Pkpo_E_quad, quad_2kp2_E) * PE;
  } // for iE

  // Face contribution
  if (degree() > 0) {
    DivergenceBasis<VEMCore::Poly2BasisFaceType> div_Pk2_F(*faceBases(iF).Polyk2);
    BGF -= GramMatrix(F, div_Pk2_F, *faceBases(iF).Polyk, int_mono_2kp2_F) * proj_function;
  } // if degree() > 0

  Eigen::MatrixXd proj_gradient = MGF.ldlt().solve(BGF);
  
  
  //------------------------------------------------------------------------------
  // DOFs of gradient in Vcurl(F)
  //------------------------------------------------------------------------------
  size_t n_dofs_edge = PolynomialSpaceDimension<Edge>::Poly(degree());
  size_t dim_Rckpo_F = PolynomialSpaceDimension<Face>::RolyCompl(degree()+1);
  size_t n_dofs_face = F.n_edges() * n_dofs_edge
                  + dim_Rckpo_F + PolynomialSpaceDimension<Face>::Poly(degree())-1;
  Eigen::MatrixXd dofs_gradient = Eigen::MatrixXd::Zero(n_dofs_face, dimensionFace(iF));
  for (size_t iE = 0; iE < F.n_edges(); iE++) {
    const Edge & E = *F.edge(iE);
    Eigen::MatrixXd GE = extendOperator(F, E, edgeOperators(E).dofs_gradient);
    dofs_gradient.middleRows(iE * n_dofs_edge, n_dofs_edge) = GE;
  }
  size_t dofs_offsetF = F.n_edges() * n_dofs_edge;
  dofs_gradient.block(dofs_offsetF, localOffset(F), dim_Rckpo_F, dim_Rckpo_F) = Eigen::MatrixXd::Identity(dim_Rckpo_F, dim_Rckpo_F);
  
  return LocalOperators(proj_gradient, dofs_gradient, proj_function);
}

//------------------------------------------------------------------------------

VGrad::LocalOperators VGrad::_compute_cell_operators(size_t iT)
{
  const Cell & T = *mesh().cell(iT);

  MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2);
  
  // If k=0, the "projection" will just return an average of vertex unkowns
  size_t dim_range_proj_function = PolynomialSpaceDimension<Cell>::Poly(degree()-1);
  if (degree()==0){
    dim_range_proj_function = 1;
  }
  Eigen::MatrixXd proj_function = Eigen::MatrixXd::Zero(dim_range_proj_function, dimensionCell(iT));
  if (degree() > 0){
    //------------------------------------------------------------------------------
    // Projection of function on P^{k-1}
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Left-hand side matrix
    DivergenceBasis<VEMCore::RolyComplBasisCellType> div_Rck_basis(*cellBases(iT).RolyComplk);
    Eigen::MatrixXd MPT = GramMatrix(T, div_Rck_basis, *cellBases(iT).Polykmo, int_mono_2kp2);

    //------------------------------------------------------------------------------
    // Right-hand side matrix
    Eigen::MatrixXd BPT = Eigen::MatrixXd::Zero(cellBases(iT).RolyComplk->dimension(), dimensionCell(iT));
    
    // Cell contribution
    BPT.rightCols(cellBases(iT).RolyComplk->dimension()) = -GramMatrix(T, *cellBases(iT).RolyComplk, int_mono_2kp2);

    // Boundary contribution
    for (size_t iF = 0; iF < T.n_faces(); iF++) {
     const Face & F = *T.face(iF);

     MonomialScalarBasisFace basis_Pk_F(F, degree());
     DecomposePoly dec(F, basis_Pk_F);
     auto RckT_dot_nF_nodes = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).RolyComplk, dec.get_nodes()), F.normal());
     Family<MonomialScalarBasisFace> RckT_dot_nF_family_PkF = dec.family(RckT_dot_nF_nodes);
     auto PF = extendOperator(T, F, faceOperators(F).proj_function);
     MonomialFaceIntegralsType int_mono_2k_F = IntegrateFaceMonomials(F, 2*degree());
     BPT += T.face_orientation(iF) * GramMatrix(F, RckT_dot_nF_family_PkF, *faceBases(F).Polyk, int_mono_2k_F) * PF;
        
    } // for iF                     
    
    proj_function = MPT.partialPivLu().solve(BPT);
  }else{
    // Brute average of vertex values, can later be made better (e.g. exact on interpolate of affine functions) 
    proj_function.leftCols(T.n_vertices()) = Eigen::VectorXd::LinSpaced(T.n_vertices(), 1, 1).transpose() / double(T.n_vertices());
  } // if degree() > 0
    
  //------------------------------------------------------------------------------
  // Projection of gradient on P^k
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  Eigen::MatrixXd MGT = GramMatrix(T, *cellBases(iT).Polyk3, int_mono_2kp2);

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BGT = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(), dimensionCell(iT));

  // Boundary contribution
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    
    DecomposePoly dec(F, MonomialScalarBasisFace(F, degree()));
    auto Pk3T_dot_nTF_nodes = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, dec.get_nodes()), T.face_normal(iF));
    auto Pk3T_dot_nTF_family_PkF = dec.family(Pk3T_dot_nTF_nodes);
    Eigen::MatrixXd PF = extendOperator(T, F, faceOperators(F).proj_function);
    MonomialFaceIntegralsType int_mono_2kpo_F = IntegrateFaceMonomials(F, 2*degree()+1);
    BGT += GramMatrix(F, Pk3T_dot_nTF_family_PkF, *faceBases(F).Polyk, int_mono_2kpo_F) * PF;
    
  } // for iF

  // Cell contribution
  if (degree() > 0) {
    DivergenceBasis<VEMCore::Poly3BasisCellType> div_Pk3_basis(*cellBases(iT).Polyk3);
    BGT -= GramMatrix(T, div_Pk3_basis, *cellBases(iT).Polykmo, int_mono_2kp2) * proj_function;
  } // if degree() > 0

  Eigen::MatrixXd proj_gradient = MGT.ldlt().solve(BGT);

  //------------------------------------------------------------------------------
  // DOFs of gradient in Vcurl(T)
  //------------------------------------------------------------------------------
  size_t n_dofs_edge = PolynomialSpaceDimension<Edge>::Poly(degree());
  size_t n_dofs_face = PolynomialSpaceDimension<Face>::RolyCompl(degree()+1)
                        + PolynomialSpaceDimension<Face>::Poly(degree())-1;
  size_t n_dofs_element = T.n_edges() * n_dofs_edge + T.n_faces() * n_dofs_face
                  + PolynomialSpaceDimension<Cell>::RolyCompl(degree())
                  + PolynomialSpaceDimension<Cell>::GolyCompl(degree()+1);
  Eigen::MatrixXd dofs_gradient = Eigen::MatrixXd::Zero(n_dofs_element, dimensionCell(iT));
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    Eigen::MatrixXd GE = extendOperator(T, E, edgeOperators(E).dofs_gradient);
    dofs_gradient.middleRows(iE*n_dofs_edge, n_dofs_edge) = GE;
  }
  size_t dofs_offsetF =  T.n_edges() * n_dofs_edge;
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    Eigen::MatrixXd GF = extendOperator(T, F, faceOperators(F).dofs_gradient);
    dofs_gradient.middleRows(dofs_offsetF + iF * n_dofs_face, n_dofs_face) = GF.bottomRows(n_dofs_face);
  }
  size_t dofs_offsetT = T.n_edges() * n_dofs_edge + T.n_faces() * n_dofs_face;
  size_t dim_Rck_T = PolynomialSpaceDimension<Cell>::RolyCompl(degree());
  dofs_gradient.block(dofs_offsetT, localOffset(T), dim_Rck_T, dim_Rck_T) = Eigen::MatrixXd::Identity(dim_Rck_T, dim_Rck_T);


  return LocalOperators(proj_gradient, dofs_gradient, proj_function);
}

//-----------------------------------------------------------------------------
// local L2 inner product
//-----------------------------------------------------------------------------
Eigen::MatrixXd VGrad::computeL2Product(
                                        const size_t iT,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pkmo_T,
                                        const IntegralWeight & weight
                                        ) const
{
  const Cell & T = *mesh().cell(iT); 
  
  // create the weighted mass matrix, with simple product if weight is constant
  Eigen::MatrixXd w_mass_Pkmo_T;
  if (degree()>0){
    if (weight.deg(T)==0){
      // constant weight
      if (mass_Pkmo_T.rows()==1){
        // We have to compute the mass matrix
        MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*degree());
        w_mass_Pkmo_T = weight.value(T, T.center_mass()) * GramMatrix(T, *cellBases(iT).Polykmo, int_mono_2k);
      }else{
        w_mass_Pkmo_T = weight.value(T, T.center_mass()) * mass_Pkmo_T;
      }
    }else{
      // weight is not constant, we create a weighted mass matrix
      QuadratureRule quad_2kpo_pw_T = generate_quadrature_rule(T, 2 * (degree() + 1) + weight.deg(T));
      auto basis_Pkmo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, quad_2kpo_pw_T);
      std::function<double(const Eigen::Vector3d &)> weight_T 
                = [&T, &weight](const Eigen::Vector3d &x)->double {
                    return weight.value(T, x);
                  };
      w_mass_Pkmo_T = compute_weighted_gram_matrix(weight_T, basis_Pkmo_T_quad, basis_Pkmo_T_quad, quad_2kpo_pw_T, "sym");
    }
  }else{
    // Probably not very good...
    w_mass_Pkmo_T = T.measure() * weight.value(T, T.center_mass()) * Eigen::MatrixXd::Identity(1,1);
  }

  // Compute matrix of L2 product, just consistent part...
  Eigen::MatrixXd L2P = Eigen::MatrixXd::Zero(dimensionCell(iT), dimensionCell(iT));

  L2P += cellOperators(iT).proj_function.transpose() * w_mass_Pkmo_T * cellOperators(iT).proj_function;

  return L2P;

}



//////////-----------------------------------------------------------------------------
////////// local L2 inner product
//////////-----------------------------------------------------------------------------
////////Eigen::MatrixXd VGrad::computeL2Product(
////////                                        const size_t iT,
////////                                        const double & penalty_factor,
////////                                        const Eigen::MatrixXd & mass_Pkpo_T,
////////                                        const IntegralWeight & weight
////////                                        ) const
////////{
////////  const Cell & T = *mesh().cell(iT); 
////////  
////////  // create the weighted mass matrix, with simple product if weight is constant
////////  Eigen::MatrixXd w_mass_Pkpo_T;
////////  if (weight.deg(T)==0){
////////    // constant weight
////////    if (mass_Pkpo_T.rows()==1){
////////      // We have to compute the mass matrix
////////      MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2);
////////      w_mass_Pkpo_T = weight.value(T, T.center_mass()) * GramMatrix(T, *cellBases(iT).Polykpo, int_mono_2kp2);
////////    }else{
////////      w_mass_Pkpo_T = weight.value(T, T.center_mass()) * mass_Pkpo_T;
////////    }
////////  }else{
////////    // weight is not constant, we create a weighted mass matrix
////////    QuadratureRule quad_2kpo_pw_T = generate_quadrature_rule(T, 2 * (degree() + 1) + weight.deg(T));
////////    auto basis_Pkpo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, quad_2kpo_pw_T);
////////    std::function<double(const Eigen::Vector3d &)> weight_T 
////////              = [&T, &weight](const Eigen::Vector3d &x)->double {
////////                  return weight.value(T, x);
////////                };
////////    w_mass_Pkpo_T = compute_weighted_gram_matrix(weight_T, basis_Pkpo_T_quad, basis_Pkpo_T_quad, quad_2kpo_pw_T, "sym");
////////  }


////////  // Compute matrix of L2 product  
////////  Eigen::MatrixXd L2P = Eigen::MatrixXd::Zero(dimensionCell(iT), dimensionCell(iT));

////////  // We need the potential in the cell
////////  Eigen::MatrixXd Potential_T = cellOperators(iT).potential;

////////  // Edge penalty terms
////////  for (size_t iE = 0; iE < T.n_edges(); iE++) {
////////    const Edge & E = *T.edge(iE);
////////        
////////    QuadratureRule quad_2kpo_E = generate_quadrature_rule(E, 2 * (degree()+1) );
////////    
////////    // weight and scaling hE^2
////////    double max_weight_quad_E = weight.value(T, quad_2kpo_E[0].vector());
////////    // If the weight is not constant, we want to take the largest along the edge
////////    if (weight.deg(T)>0){
////////      for (size_t iqn = 1; iqn < quad_2kpo_E.size(); iqn++) {
////////        max_weight_quad_E = std::max(max_weight_quad_E, weight.value(T, quad_2kpo_E[iqn].vector()));
////////      } // for
////////    }
////////    double w_hE2 = max_weight_quad_E * std::pow(E.measure(), 2);

////////    // The penalty term int_E (PT q - q_E) * (PT r - r_E) is computed by developping.
////////    auto basis_Pkpo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, quad_2kpo_E);
////////    auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(*edgeBases(E.global_index()).Polykpo, quad_2kpo_E);
////////    Eigen::MatrixXd gram_PkpoT_PkpoE = compute_gram_matrix(basis_Pkpo_T_quad, basis_Pkpo_E_quad, quad_2kpo_E);
////////    
////////    Eigen::MatrixXd Potential_E = extendOperator(T, E, edgeOperators(E).potential);

////////    // Contribution of edge E
////////    Eigen::MatrixXd PTtrans_mass_PE = Potential_T.transpose() * gram_PkpoT_PkpoE * Potential_E;
////////    L2P += w_hE2 * ( Potential_T.transpose() * compute_gram_matrix(basis_Pkpo_T_quad, quad_2kpo_E) * Potential_T
////////                   - PTtrans_mass_PE - PTtrans_mass_PE.transpose()
////////                   + Potential_E.transpose() * compute_gram_matrix(basis_Pkpo_E_quad, quad_2kpo_E) * Potential_E );
////////  } // for iE

////////  // Face penalty terms
////////  for (size_t iF = 0; iF < T.n_faces(); iF++) {
////////    const Face & F = *T.face(iF);
////////    QuadratureRule quad_2kpo_F = generate_quadrature_rule(F, 2 * (degree()+1) );

////////    // weight and scaling hF (we use quadrature nodes to evaluate the maximum of the weight)
////////    double max_weight_quad_F = weight.value(T, quad_2kpo_F[0].vector());
////////    // If the weight is not constant, we want to take the largest along the face
////////    if (weight.deg(T)>0){
////////      for (size_t iqn = 1; iqn < quad_2kpo_F.size(); iqn++) {
////////        max_weight_quad_F = std::max(max_weight_quad_F, weight.value(T, quad_2kpo_F[iqn].vector()));
////////      } // for
////////    }
////////    double w_hF = max_weight_quad_F * F.diam();

////////    // The penalty term int_F (PT q - gammaF q) * (PT r - gammaF r) is computed by developping.
////////    MonomialFaceIntegralsType int_monoF_2kp2 = IntegrateFaceMonomials(F, 2*degree()+2);
////////    DecomposePoly dec(F, MonomialScalarBasisFace(F, degree()+1));
////////    auto PkpoT_nodes = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, dec.get_nodes());
////////    auto PkpoT_family_PkpoF = dec.family(PkpoT_nodes);
////////    Eigen::MatrixXd gram_PkpoT_PkpoF = GramMatrix(F, PkpoT_family_PkpoF, *faceBases(F).Polykpo, int_monoF_2kp2);
////////  
////////    // Contribution of face F
////////    Eigen::MatrixXd Potential_F = extendOperator(T, F, faceOperators(F).potential);
////////    Eigen::MatrixXd PTtrans_mass_PF = Potential_T.transpose() * gram_PkpoT_PkpoF * Potential_F;
////////    L2P += w_hF * ( Potential_T.transpose() * GramMatrix(F, PkpoT_family_PkpoF, int_monoF_2kp2) * Potential_T
////////                 - PTtrans_mass_PF - PTtrans_mass_PF.transpose()
////////                 + Potential_F.transpose() * GramMatrix(F, *faceBases(F).Polykpo, int_monoF_2kp2) * Potential_F );
////////                   
////////    // Following commented block does the same as above, but without DecomposePoly (which sometimes increases errors)
////////    /*
////////      auto basis_Pkpo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, quad_2kpo_F);
////////      auto basis_Pkpo_F_quad = evaluate_quad<Function>::compute(*faceBases(F.global_index()).Polykpo, quad_2kpo_F);
////////      Eigen::MatrixXd gram_PkpoT_PkpoF = compute_gram_matrix(basis_Pkpo_T_quad, basis_Pkpo_F_quad, quad_2kpo_F);
////////      
////////      Eigen::MatrixXd Potential_F = extendOperator(T, F, faceOperators(F).potential);
////////      Eigen::MatrixXd PTtrans_mass_PF = Potential_T.transpose() * gram_PkpoT_PkpoF * Potential_F;

////////      // Contribution of face F
////////      L2P += w_hF * ( Potential_T.transpose() * compute_gram_matrix(basis_Pkpo_T_quad, quad_2kpo_F) * Potential_T
////////                     - PTtrans_mass_PF - PTtrans_mass_PF.transpose()
////////                     + Potential_F.transpose() * GramMatrix(F, *faceBases(F).Polykpo) * Potential_F );
////////    */       
////////    
////////    
////////  } // for iF

////////  L2P *= penalty_factor;

////////  // Cell term
////////  L2P += Potential_T.transpose() * w_mass_Pkpo_T * Potential_T;

////////  return L2P;

////////}



