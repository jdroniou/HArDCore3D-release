
#include <xgrad.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_face.hpp>
#include <GMpoly_edge.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

XGrad::XGrad(const DDRCore & ddr_core, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(ddr_core.mesh(),
	     1,
	     PolynomialSpaceDimension<Edge>::Poly(ddr_core.degree() - 1),
	     PolynomialSpaceDimension<Face>::Poly(ddr_core.degree() - 1),
	     PolynomialSpaceDimension<Cell>::Poly(ddr_core.degree() - 1)
	     ),
    m_ddr_core(ddr_core),
    m_use_threads(use_threads),
    m_output(output),
    m_edge_operators(ddr_core.mesh().n_edges()),
    m_face_operators(ddr_core.mesh().n_faces()),
    m_cell_operators(ddr_core.mesh().n_cells())
{
  m_output << "[XGrad] Initializing" << std::endl;
  if (use_threads) {
    m_output << "[XGrad] Parallel execution" << std::endl;
  } else {
    m_output << "[XGrad] Sequential execution" << std::endl;
  }
  
  // Construct edge gradients and potentials
  std::function<void(size_t, size_t)> construct_all_edge_gradients_potentials
    = [this](size_t start, size_t end)->void
      {
        for (size_t iE = start; iE < end; iE++) {
          m_edge_operators[iE].reset( new LocalOperators(_compute_edge_gradient_potential(iE)) );
        } // for iE
      };

  m_output << "[XGrad] Constructing edge gradients and potentials" << std::endl;
  parallel_for(mesh().n_edges(), construct_all_edge_gradients_potentials, use_threads);

  // Construct face gradients and potentials
  std::function<void(size_t, size_t)> construct_all_face_gradients_potentials
    = [this](size_t start, size_t end)->void
      {
        for (size_t iF = start; iF < end; iF++) {
          m_face_operators[iF].reset( new LocalOperators(_compute_face_gradient_potential(iF)) );
        } // for iF
      };

  m_output << "[XGrad] Constructing face gradients and potentials" << std::endl;
  parallel_for(mesh().n_faces(), construct_all_face_gradients_potentials, use_threads);

  // Construct cell gradients and potentials
  std::function<void(size_t, size_t)> construct_all_cell_gradients_potentials
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_cell_operators[iT].reset( new LocalOperators(_compute_cell_gradient_potential(iT)) );
        } // for iT
      };

  m_output << "[XGrad] Constructing cell gradients and potentials" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_gradients_potentials, use_threads);
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd XGrad::interpolate(const FunctionType & q, const int doe_cell, const int doe_face, const int doe_edge) const
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
    
    // Interpolate at faces
    std::function<void(size_t, size_t)> interpolate_faces
      = [this, &qh, q, &dqr_face](size_t start, size_t end)->void
        {
          for (size_t iF = start; iF < end; iF++) {
            const Face & F = *mesh().face(iF);
            QuadratureRule quad_dqr_F = generate_quadrature_rule(F, dqr_face);
            auto basis_Pkmo_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).Polykmo, quad_dqr_F);
            qh.segment(globalOffset(F), PolynomialSpaceDimension<Face>::Poly(degree() - 1)) 
              = l2_projection(q, *faceBases(iF).Polykmo, quad_dqr_F, basis_Pkmo_F_quad);
          } // for iF
        };
    parallel_for(mesh().n_faces(), interpolate_faces, m_use_threads);

    // Interpolate at cells
    std::function<void(size_t, size_t)> interpolate_cells
      = [this, &qh, q, &dqr_cell](size_t start, size_t end)->void
        {
          for (size_t iT = start; iT < end; iT++) {
            const Cell & T = *mesh().cell(iT);
            QuadratureRule quad_dqr_T = generate_quadrature_rule(T, dqr_cell);
            auto basis_Pkmo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykmo, quad_dqr_T);
            MonomialCellIntegralsType int_mono_2km2 = IntegrateCellMonomials(T, 2*degree()-2);
            qh.segment(globalOffset(T), PolynomialSpaceDimension<Cell>::Poly(degree() - 1)) 
              = l2_projection(q, *cellBases(iT).Polykmo, quad_dqr_T, basis_Pkmo_T_quad, GramMatrix(T, *cellBases(iT).Polykmo, int_mono_2km2));
          } // for iT
        };
    parallel_for(mesh().n_cells(), interpolate_cells, m_use_threads);
  } // if degree() > 0 

  return qh;
}

//------------------------------------------------------------------------------
// Gradient and potential reconstructions
//------------------------------------------------------------------------------

XGrad::LocalOperators XGrad::_compute_edge_gradient_potential(size_t iE)
{
  const Edge & E = *mesh().edge(iE);
  
  //------------------------------------------------------------------------------
  // Gradient
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  
  MonomialEdgeIntegralsType int_mono_2kp1_E = IntegrateEdgeMonomials(E, 2*degree()+1);
  Eigen::MatrixXd MGE = GramMatrix(E, *edgeBases(iE).Polyk, int_mono_2kp1_E);

  //------------------------------------------------------------------------------
  // Right-hand side matrix
  
  Eigen::MatrixXd BGE
    = Eigen::MatrixXd::Zero(edgeBases(iE).Polyk->dimension(), dimensionEdge(iE));
  for (size_t i = 0; i < edgeBases(iE).Polyk->dimension(); i++) {
    BGE(i, 0) = -edgeBases(iE).Polyk->function(i, mesh().edge(iE)->vertex(0)->coords());
    BGE(i, 1) = edgeBases(iE).Polyk->function(i, mesh().edge(iE)->vertex(1)->coords());
  } // for i
  
  if (degree() > 0) {
    GradientBasis<DDRCore::PolyBasisEdgeType> grad_Pk_E(*edgeBases(iE).Polyk);
    BGE.rightCols(PolynomialSpaceDimension<Edge>::Poly(degree() - 1)) 
          = -GramMatrix(E, grad_Pk_E, *edgeBases(iE).Polykmo, int_mono_2kp1_E);
  }
  
  //------------------------------------------------------------------------------
  // Potential
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

  GradientBasis<DDRCore::PolyBasisEdgeType> grad_Pkpo_E(*edgeBases(iE).Polykpo);				    
  MPE.topRows(PolynomialSpaceDimension<Edge>::Poly(degree()))
        = GramMatrix(E, *edgeBases(iE).Polyk, grad_Pkpo_E, int_mono_2kp1_E);

  MonomialScalarBasisEdge basis_P0_E(E, 0);
  MPE.bottomRows(1) = GramMatrix(E, basis_P0_E, *edgeBases(iE).Polykpo, int_mono_2kp1_E);

  return LocalOperators(MGE.ldlt().solve(BGE), MPE.partialPivLu().solve(BPE));
}

//------------------------------------------------------------------------------

XGrad::LocalOperators XGrad::_compute_face_gradient_potential(size_t iF)
{
  const Face & F = *mesh().face(iF);

  //------------------------------------------------------------------------------
  // Gradient
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  MonomialFaceIntegralsType int_mono_2kp3_F = IntegrateFaceMonomials(F, 2*degree()+3);
  Eigen::MatrixXd MGF = GramMatrix(F, *faceBases(iF).Polyk2, int_mono_2kp3_F);

  //------------------------------------------------------------------------------
  // Right-hand side matrix
  
  Eigen::MatrixXd BGF
    = Eigen::MatrixXd::Zero(faceBases(iF).Polyk2->dimension(), dimensionFace(iF));

  // Boundary contribution
  for (size_t iE = 0; iE < F.n_edges(); iE++) {
    const Edge & E = *F.edge(iE);
    
    QuadratureRule quad_2kpo_E = generate_quadrature_rule(E, 2 * (degree() + 1));
    auto basis_Pk2_nFE_E_quad
      = scalar_product(evaluate_quad<Function>::compute(*faceBases(iF).Polyk2, quad_2kpo_E), F.edge_normal(iE));
    auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polykpo, quad_2kpo_E);
    Eigen::MatrixXd PE = extendOperator(F, E, edgeOperators(E).potential);
    BGF += F.edge_orientation(iE) * compute_gram_matrix(basis_Pk2_nFE_E_quad, basis_Pkpo_E_quad, quad_2kpo_E) * PE;
  } // for iE

  // Face contribution
  if (degree() > 0) {
    DivergenceBasis<DDRCore::Poly2BasisFaceType> div_Pk2_F(*faceBases(iF).Polyk2);
    BGF.rightCols(PolynomialSpaceDimension<Face>::Poly(degree() - 1))
      -= GramMatrix(F, div_Pk2_F, *faceBases(iF).Polykmo, int_mono_2kp3_F);
  } // if degree() > 0

  Eigen::MatrixXd GF = MGF.ldlt().solve(BGF);
  
  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  
  DivergenceBasis<DDRCore::RolyComplBasisFaceType> div_Rckp2_F(*faceBases(iF).RolyComplkp2);
  Eigen::MatrixXd MPF = GramMatrix(F, div_Rckp2_F, *faceBases(iF).Polykpo, int_mono_2kp3_F);

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  // Face contribution: we need integrals up to 2k+3 here because Polyk2 is a restriction of a basis of degree k+1
  Eigen::MatrixXd BPF = -GramMatrix(F, *faceBases(iF).RolyComplkp2, *faceBases(iF).Polyk2, int_mono_2kp3_F) * GF;
  
  // Boundary contribution
  for (size_t iE = 0; iE < F.n_edges(); iE++) {
    const Edge & E = *F.edge(iE);
    
    QuadratureRule quad_2kp4_E = generate_quadrature_rule(E, 2 * (degree() + 2));
    auto basis_Rckp2_nFE_E_quad
      = scalar_product(evaluate_quad<Function>::compute(*faceBases(iF).RolyComplkp2, quad_2kp4_E), F.edge_normal(iE));
    auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polykpo, quad_2kp4_E);
    Eigen::MatrixXd PE = extendOperator(F, E, edgeOperators(E).potential);
    BPF += F.edge_orientation(iE) * compute_gram_matrix(basis_Rckp2_nFE_E_quad, basis_Pkpo_E_quad, quad_2kp4_E) * PE;
  } // for iE
  
  return LocalOperators(GF, MPF.partialPivLu().solve(BPF));
}

//------------------------------------------------------------------------------

XGrad::LocalOperators XGrad::_compute_cell_gradient_potential(size_t iT)
{
  const Cell & T = *mesh().cell(iT);

  //------------------------------------------------------------------------------
  // Gradient
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix

  // Compute all integrals of monomial powers to degree 2k+2 and the mass matrix
  MonomialCellIntegralsType int_mono_2kp3 = IntegrateCellMonomials(T, 2*degree()+3);
  Eigen::MatrixXd MGT = GramMatrix(T, *cellBases(iT).Polyk3, int_mono_2kp3);

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BGT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(), dimensionCell(iT));

  // Boundary contribution
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    
    DecomposePoly dec(F, MonomialScalarBasisFace(F, degree()));
    auto Pk3T_dot_nTF_nodes = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, dec.get_nodes()), T.face_normal(iF));
    auto Pk3T_dot_nTF_family_PkF = dec.family(Pk3T_dot_nTF_nodes);
    Eigen::MatrixXd PF = extendOperator(T, F, faceOperators(F).potential);
    MonomialFaceIntegralsType int_mono_2kp1_F = IntegrateFaceMonomials(F, 2*degree()+1);
    BGT += GramMatrix(F, Pk3T_dot_nTF_family_PkF, *faceBases(F).Polykpo, int_mono_2kp1_F) * PF;

    // Following commented block could replace the block above, without using DecomposePoly (more expensive, but sometimes better rounding errors)
    /*
      QuadratureRule quad_2kpo_F = generate_quadrature_rule(F, 2 * degree() + 1);
      auto basis_Pk3_nTF_F_quad = scalar_product(
					         evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2kpo_F),
					         T.face_normal(iF)
					         );
      auto basis_Pkpo_F_quad = evaluate_quad<Function>::compute(*faceBases(F).Polykpo, quad_2kpo_F);
      Eigen::MatrixXd PF = extendOperator(T, F, faceOperators(F).potential);
      BGT += compute_gram_matrix(basis_Pk3_nTF_F_quad, basis_Pkpo_F_quad, quad_2kpo_F) * PF;
    */
    
  } // for iF

  // Cell contribution
  if (degree() > 0) {
    DivergenceBasis<DDRCore::Poly3BasisCellType> div_Pk3_basis(*cellBases(iT).Polyk3);
    BGT.rightCols(PolynomialSpaceDimension<Cell>::Poly(degree() - 1)) -= GramMatrix(T, div_Pk3_basis, *cellBases(iT).Polykmo, int_mono_2kp3);
  } // if degree() > 0

  Eigen::MatrixXd GT = MGT.ldlt().solve(BGT);
  
  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix

  DivergenceBasis<DDRCore::RolyComplBasisCellType> div_Rckp2_basis(*cellBases(iT).RolyComplkp2);
  Eigen::MatrixXd MPT = GramMatrix(T, div_Rckp2_basis, *cellBases(iT).Polykpo, int_mono_2kp3);

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  // Cell contribution: we need integrals up to degree 2k+3 because Polyk3 comes from the restriction of a basis of degree k+1
  Eigen::MatrixXd BPT
   = -GramMatrix(T, *cellBases(iT).RolyComplkp2, *cellBases(iT).Polyk3, int_mono_2kp3) * GT;

  // Boundary contribution
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
   const Face & F = *T.face(iF);

   MonomialScalarBasisFace basis_Pkp2_F(F, degree()+2);
   DecomposePoly dec(F, basis_Pkp2_F);
   auto Rckp2T_dot_nF_nodes = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).RolyComplkp2, dec.get_nodes()), F.normal());
   Family<MonomialScalarBasisFace> Rckp2T_dot_nF_family_Pkp2F = dec.family(Rckp2T_dot_nF_nodes);
   auto PF = extendOperator(T, F, faceOperators(F).potential);
   MonomialFaceIntegralsType int_mono_2kp3_F = IntegrateFaceMonomials(F, 2*degree()+3);
   BPT += T.face_orientation(iF) * GramMatrix(F, Rckp2T_dot_nF_family_Pkp2F, *faceBases(F).Polykpo, int_mono_2kp3_F) * PF;

   // Following commented block does the same as above, but with DecomposePoly and seems to lead to increased errors
   /*
     QuadratureRule quad_2kp2_F = generate_quadrature_rule(F, 2 * (degree() + 2));
     auto basis_Rckp2_nF_F_quad
       = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).RolyComplkp2, quad_2kp2_F), F.normal());
     auto basis_Pkpo_F_quad = evaluate_quad<Function>::compute(*faceBases(F).Polykpo, quad_2kp2_F);
     auto PF = extendOperator(T, F, faceOperators(F).potential);
     BPT += T.face_orientation(iF) * compute_gram_matrix(basis_Rckp2_nF_F_quad, basis_Pkpo_F_quad, quad_2kp2_F) * PF;
   */
   
  } // for iF                                              

  return LocalOperators(GT, MPT.partialPivLu().solve(BPT));
}

//-----------------------------------------------------------------------------
// local L2 inner product
//-----------------------------------------------------------------------------
Eigen::MatrixXd XGrad::computeL2Product(
                                        const size_t iT,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pkpo_T,
                                        const IntegralWeight & weight
                                        ) const
{
  const Cell & T = *mesh().cell(iT); 
  
  // create the weighted mass matrix, with simple product if weight is constant
  Eigen::MatrixXd w_mass_Pkpo_T;
  if (weight.deg(T)==0){
    // constant weight
    if (mass_Pkpo_T.rows()==1){
      // We have to compute the mass matrix
      MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2);
      w_mass_Pkpo_T = weight.value(T, T.center_mass()) * GramMatrix(T, *cellBases(iT).Polykpo, int_mono_2kp2);
    }else{
      w_mass_Pkpo_T = weight.value(T, T.center_mass()) * mass_Pkpo_T;
    }
  }else{
    // weight is not constant, we create a weighted mass matrix
    QuadratureRule quad_2kpo_pw_T = generate_quadrature_rule(T, 2 * (degree() + 1) + weight.deg(T));
    auto basis_Pkpo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, quad_2kpo_pw_T);
    std::function<double(const Eigen::Vector3d &)> weight_T 
              = [&T, &weight](const Eigen::Vector3d &x)->double {
                  return weight.value(T, x);
                };
    w_mass_Pkpo_T = compute_weighted_gram_matrix(weight_T, basis_Pkpo_T_quad, basis_Pkpo_T_quad, quad_2kpo_pw_T, "sym");
  }

  // Potential in the cell
  Eigen::MatrixXd Potential_T = cellOperators(iT).potential;

  // Consistent term and stabilisation
  return Potential_T.transpose() * w_mass_Pkpo_T * Potential_T 
          + penalty_factor * computeStabilisation(iT, weight);

}

Eigen::MatrixXd XGrad::computeStabilisation(
                                        const size_t iT,
                                        const IntegralWeight & weight
                                        ) const
{
  const Cell & T = *mesh().cell(iT); 
  
  Eigen::MatrixXd ST = Eigen::MatrixXd::Zero(dimensionCell(iT), dimensionCell(iT));

  // We need the potential in the cell
  Eigen::MatrixXd Potential_T = cellOperators(iT).potential;

  // Edge penalty terms
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
        
    QuadratureRule quad_2kpo_E = generate_quadrature_rule(E, 2 * (degree()+1) );
    
    // weight and scaling hE^2
    double max_weight_quad_E = weight.value(T, quad_2kpo_E[0].vector());
    // If the weight is not constant, we want to take the largest along the edge
    if (weight.deg(T)>0){
      for (size_t iqn = 1; iqn < quad_2kpo_E.size(); iqn++) {
        max_weight_quad_E = std::max(max_weight_quad_E, weight.value(T, quad_2kpo_E[iqn].vector()));
      } // for
    }
    double w_hE2 = max_weight_quad_E * std::pow(E.diam(), 2);

    // The penalty term int_E (PT q - q_E) * (PT r - r_E) is computed by developping.
    auto basis_Pkpo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, quad_2kpo_E);
    auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(*edgeBases(E.global_index()).Polykpo, quad_2kpo_E);
    Eigen::MatrixXd gram_PkpoT_PkpoE = compute_gram_matrix(basis_Pkpo_T_quad, basis_Pkpo_E_quad, quad_2kpo_E);
    
    Eigen::MatrixXd Potential_E = extendOperator(T, E, edgeOperators(E).potential);

    // Contribution of edge E
    Eigen::MatrixXd PTtrans_mass_PE = Potential_T.transpose() * gram_PkpoT_PkpoE * Potential_E;
    ST += 
        w_hE2 * ( Potential_T.transpose() * compute_gram_matrix(basis_Pkpo_T_quad, quad_2kpo_E) * Potential_T
                  - PTtrans_mass_PE - PTtrans_mass_PE.transpose()
                  + Potential_E.transpose() * compute_gram_matrix(basis_Pkpo_E_quad, quad_2kpo_E) * Potential_E );
  } // for iE

  // Face penalty terms
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    QuadratureRule quad_2kpo_F = generate_quadrature_rule(F, 2 * (degree()+1) );

    // weight and scaling hF (we use quadrature nodes to evaluate the maximum of the weight)
    double max_weight_quad_F = weight.value(T, quad_2kpo_F[0].vector());
    // If the weight is not constant, we want to take the largest along the face
    if (weight.deg(T)>0){
      for (size_t iqn = 1; iqn < quad_2kpo_F.size(); iqn++) {
        max_weight_quad_F = std::max(max_weight_quad_F, weight.value(T, quad_2kpo_F[iqn].vector()));
      } // for
    }
    double w_hF = max_weight_quad_F * F.diam();

    // The penalty term int_F (PT q - gammaF q) * (PT r - gammaF r) is computed by developping.
    MonomialFaceIntegralsType int_monoF_2kp2 = IntegrateFaceMonomials(F, 2*degree()+2);
    DecomposePoly dec(F, MonomialScalarBasisFace(F, degree()+1));
    auto PkpoT_nodes = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, dec.get_nodes());
    auto PkpoT_family_PkpoF = dec.family(PkpoT_nodes);
    Eigen::MatrixXd gram_PkpoT_PkpoF = GramMatrix(F, PkpoT_family_PkpoF, *faceBases(F).Polykpo, int_monoF_2kp2);
  
    // Contribution of face F
    Eigen::MatrixXd Potential_F = extendOperator(T, F, faceOperators(F).potential);
    Eigen::MatrixXd PTtrans_mass_PF = Potential_T.transpose() * gram_PkpoT_PkpoF * Potential_F;
    ST += w_hF * ( Potential_T.transpose() * GramMatrix(F, PkpoT_family_PkpoF, int_monoF_2kp2) * Potential_T
                 - PTtrans_mass_PF - PTtrans_mass_PF.transpose()
                 + Potential_F.transpose() * GramMatrix(F, *faceBases(F).Polykpo, int_monoF_2kp2) * Potential_F );
                   
    // Following commented block does the same as above, but without DecomposePoly (which sometimes increases errors)
    /*
      auto basis_Pkpo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, quad_2kpo_F);
      auto basis_Pkpo_F_quad = evaluate_quad<Function>::compute(*faceBases(F.global_index()).Polykpo, quad_2kpo_F);
      Eigen::MatrixXd gram_PkpoT_PkpoF = compute_gram_matrix(basis_Pkpo_T_quad, basis_Pkpo_F_quad, quad_2kpo_F);
      
      Eigen::MatrixXd Potential_F = extendOperator(T, F, faceOperators(F).potential);
      Eigen::MatrixXd PTtrans_mass_PF = Potential_T.transpose() * gram_PkpoT_PkpoF * Potential_F;

      // Contribution of face F
      ST += w_hF * ( Potential_T.transpose() * compute_gram_matrix(basis_Pkpo_T_quad, quad_2kpo_F) * Potential_T
                     - PTtrans_mass_PF - PTtrans_mass_PF.transpose()
                     + Potential_F.transpose() * GramMatrix(F, *faceBases(F).Polykpo) * Potential_F );
    */       
    
    
  } // for iF

  return ST;

}


//------------------------------------------------------------------------------------
// Build the components of the gradient operator (probably never useful, actually...)
//------------------------------------------------------------------------------------

Eigen::MatrixXd XGrad::buildGradientComponentsCell(size_t iT) const
{
  const Cell & T = *m_mesh.cell(iT);
  
  size_t dim_xcurl_T
    = T.n_edges() * PolynomialSpaceDimension<Edge>::Poly(degree())
    + T.n_faces() * (PolynomialSpaceDimension<Face>::Roly(degree()-1) + PolynomialSpaceDimension<Face>::RolyCompl(degree()))
    + PolynomialSpaceDimension<Cell>::Roly(degree()-1) + PolynomialSpaceDimension<Cell>::RolyCompl(degree());

  size_t dim_xgrad_T
    = T.n_vertices()
    + T.n_edges() * PolynomialSpaceDimension<Edge>::Poly(degree()-1)
    + T.n_faces() * PolynomialSpaceDimension<Face>::Poly(degree()-1)
    + PolynomialSpaceDimension<Cell>::Poly(degree()-1);

  Eigen::MatrixXd uGT = Eigen::MatrixXd::Zero(dim_xcurl_T, dim_xgrad_T);

  // Edge components
  size_t offset = 0;
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    uGT.block(offset, 0, PolynomialSpaceDimension<Edge>::Poly(degree()), dim_xgrad_T)
      = extendOperator(T, E, edgeOperators(E).gradient);
    offset += PolynomialSpaceDimension<Edge>::Poly(degree());
  } // for iE

  if (m_ddr_core.degree() > 0) {
    // Face components
    for (size_t iF = 0; iF < T.n_faces(); iF++) {
      const Face & F = *T.face(iF);

      MonomialFaceIntegralsType int_monoF_2k = IntegrateFaceMonomials(F, 2*degree());
      Eigen::MatrixXd mass_Rck_F = GramMatrix(F, *faceBases(F).RolyComplk, int_monoF_2k);
      Eigen::MatrixXd mass_Rkmo_F = GramMatrix(F, *faceBases(F).Rolykmo, int_monoF_2k);

      auto GF = extendOperator(T, F, faceOperators(F).gradient);
      Eigen::MatrixXd pi_Rkmo_GF_F = mass_Rkmo_F.ldlt().solve(GramMatrix(F, *faceBases(F).Rolykmo, *faceBases(F).Polyk2, int_monoF_2k) * GF);

      Eigen::MatrixXd pi_Rck_GF_F = mass_Rck_F.ldlt().solve(GramMatrix(F, *faceBases(F).RolyComplk, *faceBases(F).Polyk2, int_monoF_2k) * GF);

      uGT.block(offset, 0, PolynomialSpaceDimension<Face>::Roly(degree()-1), dim_xgrad_T) = pi_Rkmo_GF_F;
      offset += PolynomialSpaceDimension<Face>::Roly(degree()-1);
      uGT.block(offset, 0, PolynomialSpaceDimension<Face>::RolyCompl(degree()), dim_xgrad_T) = pi_Rck_GF_F;
      offset += PolynomialSpaceDimension<Face>::RolyCompl(degree());                                                     
    } // for iF

    // Cell component
    MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*degree());
    Eigen::MatrixXd mass_Rkmo_T = GramMatrix(T, *m_ddr_core.cellBases(iT).Rolykmo, int_mono_2k);
    Eigen::MatrixXd mass_Rck_T = GramMatrix(T, *m_ddr_core.cellBases(iT).RolyComplk, int_mono_2k);

    Eigen::MatrixXd pi_Rkmo_GT_T = mass_Rkmo_T.ldlt().solve(GramMatrix(T, *m_ddr_core.cellBases(iT).Rolykmo, *m_ddr_core.cellBases(iT).Polyk3, int_mono_2k) * cellOperators(iT).gradient);
    Eigen::MatrixXd pi_Rck_GT_T = mass_Rck_T.ldlt().solve(
                                      GramMatrix(T, *m_ddr_core.cellBases(iT).RolyComplk, *m_ddr_core.cellBases(iT).Polyk3, int_mono_2k)
                                      * cellOperators(iT).gradient
                                  );

    uGT.block(offset, 0, PolynomialSpaceDimension<Cell>::Roly(degree()-1), dim_xgrad_T) = pi_Rkmo_GT_T;
    offset += PolynomialSpaceDimension<Cell>::Roly(degree()-1);
    uGT.block(offset, 0, PolynomialSpaceDimension<Cell>::RolyCompl(degree()), dim_xgrad_T) = pi_Rck_GT_T;
    offset += PolynomialSpaceDimension<Cell>::RolyCompl(degree());
  }
  
  return uGT;
}
