
#include <xgrad.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

XGrad::XGrad(const DDRCore & ddr_core, bool use_threads, std::ostream & output)
  : DDRSpace(ddr_core.mesh(),
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

Eigen::VectorXd XGrad::interpolate(const FunctionType & q) const
{
  Eigen::VectorXd qh = Eigen::VectorXd::Zero(dimension());

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
      = [this, &qh, q](size_t start, size_t end)->void
        {
          for (size_t iE = start; iE < end; iE++) {
            const Edge & E = *mesh().edge(iE);
            QuadratureRule quad = generate_quadrature_rule(E, 2 * degree());
            auto basis_Pkmo_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polykmo, quad);
            qh.segment(globalOffset(E), PolynomialSpaceDimension<Edge>::Poly(degree() - 1)) 
              = l2_projection(q, *edgeBases(iE).Polykmo, quad, basis_Pkmo_E_quad);
          } // for iE
        };
    parallel_for(mesh().n_edges(), interpolate_edges, m_use_threads);
    
    // Interpolate at faces
    std::function<void(size_t, size_t)> interpolate_faces
      = [this, &qh, q](size_t start, size_t end)->void
        {
          for (size_t iF = start; iF < end; iF++) {
            const Face & F = *mesh().face(iF);
            QuadratureRule quad = generate_quadrature_rule(F, 2 * degree());
            auto basis_Pkmo_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).Polykmo, quad);
            qh.segment(globalOffset(F), PolynomialSpaceDimension<Face>::Poly(degree() - 1)) 
              = l2_projection(q, *faceBases(iF).Polykmo, quad, basis_Pkmo_F_quad);
          } // for iF
        };
    parallel_for(mesh().n_faces(), interpolate_faces, m_use_threads);

    // Interpolate at cells
    std::function<void(size_t, size_t)> interpolate_cells
      = [this, &qh, q](size_t start, size_t end)->void
        {
          for (size_t iT = start; iT < end; iT++) {
            const Cell & T = *mesh().cell(iT);
            QuadratureRule quad = generate_quadrature_rule(T, 2 * degree());
            auto basis_Pkmo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykmo, quad);
            qh.segment(globalOffset(T), PolynomialSpaceDimension<Cell>::Poly(degree() - 1)) 
              = l2_projection(q, *cellBases(iT).Polykmo, quad, basis_Pkmo_T_quad);
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
  
  QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2 * degree());
  auto basis_Pk_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polyk, quad_2k_E);
  auto MGE = compute_gram_matrix(basis_Pk_E_quad, basis_Pk_E_quad, quad_2k_E, "sym");

  //------------------------------------------------------------------------------
  // Right-hand side matrix
  
  Eigen::MatrixXd BGE
    = Eigen::MatrixXd::Zero(edgeBases(iE).Polyk->dimension(), dimensionEdge(iE));
  for (size_t i = 0; i < edgeBases(iE).Polyk->dimension(); i++) {
    BGE(i, 0) = -edgeBases(iE).Polyk->function(i, mesh().edge(iE)->vertex(0)->coords());
    BGE(i, 1) = edgeBases(iE).Polyk->function(i, mesh().edge(iE)->vertex(1)->coords());
  } // for i

  QuadratureRule quad_2kmo_E = generate_quadrature_rule(E, 2 * (degree() - 1));
  
  if (degree() > 0) {    
    auto grad_Pk_tE_E_quad = scalar_product(
                                            evaluate_quad<Gradient>::compute(*edgeBases(iE).Polyk, quad_2kmo_E),
                                            E.tangent()
                                            );
    auto basis_Pkmo_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polykmo, quad_2kmo_E);
    BGE.rightCols(PolynomialSpaceDimension<Edge>::Poly(degree() - 1))
      = -compute_gram_matrix(grad_Pk_tE_E_quad, basis_Pkmo_E_quad, quad_2kmo_E);
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
  auto grad_Pkpo_tE_E_quad = scalar_product(
					    evaluate_quad<Gradient>::compute(*edgeBases(iE).Polykpo, quad_2k_E),
					    E.tangent()
					    );
  MPE.topRows(PolynomialSpaceDimension<Edge>::Poly(degree()))
    = compute_gram_matrix(basis_Pk_E_quad, grad_Pkpo_tE_E_quad, quad_2k_E);
  
  QuadratureRule quad_kpo_E = generate_quadrature_rule(E, degree() + 1);  
  auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polykpo, quad_kpo_E);  
  for (size_t i = 0; i < PolynomialSpaceDimension<Edge>::Poly(degree() + 1); i++) {
    for (size_t iqn = 0; iqn < quad_kpo_E.size(); iqn++) {
      MPE.bottomRows(1)(0, i) += quad_kpo_E[iqn].w * basis_Pkpo_E_quad[i][iqn];
    } // for iqn
  } // for i
  
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
  
  QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree());
  auto basis_Pk2_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).Polyk2, quad_2k_F);
  auto MGF = compute_gram_matrix(basis_Pk2_F_quad, basis_Pk2_F_quad, quad_2k_F, "sym");

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
    Eigen::MatrixXd BGF_E
      = F.edge_orientation(iE) * compute_gram_matrix(basis_Pk2_nFE_E_quad, basis_Pkpo_E_quad, quad_2kpo_E) * edgeOperators(E).potential;

    // Assemble local contribution
    BGF.col(localOffset(F, *E.vertex(0))) += BGF_E.col(0);
    BGF.col(localOffset(F, *E.vertex(1))) += BGF_E.col(1);
    if (degree() > 0) {
      BGF.block(0, localOffset(F, E), faceBases(iF).Polyk2->dimension(), PolynomialSpaceDimension<Edge>::Poly(degree() - 1))
        = BGF_E.rightCols(PolynomialSpaceDimension<Edge>::Poly(degree() - 1));
    } // if degree() > 0
  } // for iE

  // Face contribution
  if (degree() > 0) {
    QuadratureRule quad_2kmo_F = generate_quadrature_rule(F, 2 * (degree() - 1));
    auto div_Pk2_F_quad = evaluate_quad<Divergence>::compute(*faceBases(iF).Polyk2, quad_2kmo_F);
    auto basis_Pkmo_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).Polykmo, quad_2kmo_F);

    BGF.rightCols(PolynomialSpaceDimension<Face>::Poly(degree() - 1))
      -= compute_gram_matrix(div_Pk2_F_quad, basis_Pkmo_F_quad, quad_2kmo_F);
  } // if degree() > 0

  Eigen::MatrixXd GF = MGF.ldlt().solve(BGF);
  
  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  
  QuadratureRule quad_2kpo_F = generate_quadrature_rule(F, 2 * (degree() + 1));
  auto basis_Pkpo_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).Polykpo, quad_2kpo_F);
  auto div_Rckp2_F_quad = evaluate_quad<Divergence>::compute(*faceBases(iF).RolyComplkp2, quad_2kpo_F);
  Eigen::MatrixXd MPF = compute_gram_matrix(div_Rckp2_F_quad, basis_Pkpo_F_quad, quad_2kpo_F);

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  // Face contribution
  Eigen::MatrixXd BPF
    = -compute_gram_matrix(
			   evaluate_quad<Function>::compute(*faceBases(iF).RolyComplkp2, quad_2kpo_F),
			   evaluate_quad<Function>::compute(*faceBases(iF).Polyk2, quad_2kpo_F),
			   quad_2kpo_F
			   ) * GF;
  
  // Boundary contribution
  for (size_t iE = 0; iE < F.n_edges(); iE++) {
    const Edge & E = *F.edge(iE);
    
    QuadratureRule quad_2kp2_E = generate_quadrature_rule(E, 2 * (degree() + 2));
    auto basis_Rckp2_nFE_E_quad
      = scalar_product(evaluate_quad<Function>::compute(*faceBases(iF).RolyComplkp2, quad_2kp2_E), F.edge_normal(iE));
    auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polykpo, quad_2kp2_E);
    Eigen::MatrixXd BPF_E
      = F.edge_orientation(iE) * compute_gram_matrix(basis_Rckp2_nFE_E_quad, basis_Pkpo_E_quad, quad_2kp2_E) * edgeOperators(E).potential;

    // Assemble local contribution
    BPF.col(localOffset(F, *E.vertex(0))) += BPF_E.col(0);
    BPF.col(localOffset(F, *E.vertex(1))) += BPF_E.col(1);
    if (degree() > 0) {
      BPF.block(0, localOffset(F, E), faceBases(iF).RolyComplkp2->dimension(), PolynomialSpaceDimension<Edge>::Poly(degree() - 1))
        += BPF_E.rightCols(PolynomialSpaceDimension<Edge>::Poly(degree() - 1));
    } // if degree() > 0
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

  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * degree());
  auto basis_Pk3_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_T);
  auto MGT = compute_gram_matrix(basis_Pk3_T_quad, basis_Pk3_T_quad, quad_2k_T, "sym");

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BGT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(), dimensionCell(iT));

  // Boundary contribution
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    
    QuadratureRule quad_2kpo_F = generate_quadrature_rule(F, 2 * degree() + 1);
    auto basis_Pk3_nTF_F_quad = scalar_product(
					       evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2kpo_F),
					       T.face_normal(iF)
					       );
    auto basis_Pkpo_F_quad = evaluate_quad<Function>::compute(*faceBases(F).Polykpo, quad_2kpo_F);
    Eigen::MatrixXd BGT_F
      = compute_gram_matrix(basis_Pk3_nTF_F_quad, basis_Pkpo_F_quad, quad_2kpo_F) * faceOperators(F).potential;

    // Check that we got the dimension right
    assert( (size_t)BGT_F.cols() == dimensionFace(F.global_index()) );
    
    // Assemble local contribution
    for (size_t iV = 0; iV < F.n_vertices(); iV++) {
      BGT.col(localOffset(T, *F.vertex(iV))) += BGT_F.col(iV);
    } // for iV

    if (degree() > 0) {
      for (size_t iE = 0; iE < F.n_edges(); iE++) {
	const Edge & E = *F.edge(iE);
        BGT.block(0, localOffset(T, E), cellBases(iF).Polyk3->dimension(), PolynomialSpaceDimension<Edge>::Poly(degree() - 1))
          += BGT_F.block(0, localOffset(F, E), cellBases(iF).Polyk3->dimension(), PolynomialSpaceDimension<Edge>::Poly(degree() - 1));
      } // for iE     

      BGT.block(0, localOffset(T, F), cellBases(iF).Polyk3->dimension(), PolynomialSpaceDimension<Face>::Poly(degree() - 1))
        += BGT_F.rightCols(PolynomialSpaceDimension<Face>::Poly(degree() - 1));
    } // if degree() > 0
  } // for iF

  // Cell contribution
  if (degree() > 0) {
    QuadratureRule quad_2kmo_T = generate_quadrature_rule(T, 2 * (degree() - 1)); 
    auto div_Pk3_T_quad = evaluate_quad<Divergence>::compute(*cellBases(iT).Polyk3, quad_2kmo_T);
    auto basis_Pkmo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykmo, quad_2kmo_T);

    BGT.rightCols(PolynomialSpaceDimension<Cell>::Poly(degree() - 1))
      -= compute_gram_matrix(div_Pk3_T_quad, basis_Pkmo_T_quad, quad_2kmo_T);
  } // if degree() > 0

  Eigen::MatrixXd GT = MGT.ldlt().solve(BGT);
  
  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix

  QuadratureRule quad_2kpo_T = generate_quadrature_rule(T, 2 * (degree() + 1));
  auto basis_Pkpo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, quad_2kpo_T);
  auto div_Rckp2_T_quad = evaluate_quad<Divergence>::compute(*cellBases(iT).RolyComplkp2, quad_2kpo_T);
  Eigen::MatrixXd MPT = compute_gram_matrix(div_Rckp2_T_quad, basis_Pkpo_T_quad, quad_2kpo_T);

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  // Cell contribution
  Eigen::MatrixXd BPT
   = -compute_gram_matrix(
                          evaluate_quad<Function>::compute(*cellBases(iT).RolyComplkp2, quad_2kpo_T),
                          evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2kpo_T),
                          quad_2kpo_T
                          ) * GT;

  // Boundary contribution
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
   const Face & F = *T.face(iF);

   QuadratureRule quad_2kp2_F = generate_quadrature_rule(F, 2 * (degree() + 2));
   auto basis_Rckp2_nF_F_quad
     = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).RolyComplkp2, quad_2kp2_F), F.normal());
   auto basis_Pkpo_F_quad = evaluate_quad<Function>::compute(*faceBases(F).Polykpo, quad_2kp2_F);
   auto PF = extendOperator(T, F, faceOperators(F.global_index()).potential);
   BPT += T.face_orientation(iF) * compute_gram_matrix(basis_Rckp2_nF_F_quad, basis_Pkpo_F_quad, quad_2kp2_F) * PF;
  } // for iF                                              

  return LocalOperators(GT, MPT.partialPivLu().solve(BPT));
}

//------------------------------------------------------------------------------
// Build the components of the gradient operator
//------------------------------------------------------------------------------

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
      = extendOperator(T, E, edgeOperators(E.global_index()).gradient);
    offset += PolynomialSpaceDimension<Edge>::Poly(degree());
  } // for iE

  if (m_ddr_core.degree() > 0) {
    // Face components
    for (size_t iF = 0; iF < T.n_faces(); iF++) {
      const Face & F = *T.face(iF);

      QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree());
      auto basis_Pk2_F_quad = evaluate_quad<Function>::compute(*m_ddr_core.faceBases(F.global_index()).Polyk2, quad_2k_F);
      auto basis_Rkmo_F_quad = evaluate_quad<Function>::compute(*m_ddr_core.faceBases(F.global_index()).Rolykmo, quad_2k_F);
      auto basis_Rck_F_quad = evaluate_quad<Function>::compute(*m_ddr_core.faceBases(F.global_index()).RolyComplk, quad_2k_F);
    
      Eigen::MatrixXd mass_Rkmo_F = compute_gram_matrix(basis_Rkmo_F_quad, quad_2k_F);
      Eigen::MatrixXd mass_Rck_F = compute_gram_matrix(basis_Rck_F_quad, quad_2k_F);

      auto GF = extendOperator(T, F, faceOperators(F.global_index()).gradient);

      Eigen::MatrixXd pi_Rkmo_GF_F = mass_Rkmo_F.ldlt().solve(compute_gram_matrix(basis_Rkmo_F_quad, basis_Pk2_F_quad, quad_2k_F) * GF);
      Eigen::MatrixXd pi_Rck_GF_F = mass_Rck_F.ldlt().solve(compute_gram_matrix(basis_Rck_F_quad, basis_Pk2_F_quad, quad_2k_F) * GF);

      uGT.block(offset, 0, PolynomialSpaceDimension<Face>::Roly(degree()-1), dim_xgrad_T) = pi_Rkmo_GF_F;
      offset += PolynomialSpaceDimension<Face>::Roly(degree()-1);
      uGT.block(offset, 0, PolynomialSpaceDimension<Face>::RolyCompl(degree()), dim_xgrad_T) = pi_Rck_GF_F;
      offset += PolynomialSpaceDimension<Face>::RolyCompl(degree());                                                     
    } // for iF

    // Cell component
    QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * degree());
    auto basis_Pk3_T_quad = evaluate_quad<Function>::compute(*m_ddr_core.cellBases(iT).Polyk3, quad_2k_T);
    auto basis_Rkmo_T_quad = evaluate_quad<Function>::compute(*m_ddr_core.cellBases(iT).Rolykmo, quad_2k_T);
    auto basis_Rck_T_quad = evaluate_quad<Function>::compute(*m_ddr_core.cellBases(iT).RolyComplk, quad_2k_T);
    
    Eigen::MatrixXd mass_Rkmo_T = compute_gram_matrix(basis_Rkmo_T_quad, quad_2k_T);
    Eigen::MatrixXd mass_Rck_T = compute_gram_matrix(basis_Rck_T_quad, quad_2k_T);

    Eigen::MatrixXd pi_Rkmo_GT_T = mass_Rkmo_T.ldlt().solve(compute_gram_matrix(basis_Rkmo_T_quad, basis_Pk3_T_quad, quad_2k_T) * cellOperators(iT).gradient);
    Eigen::MatrixXd pi_Rck_GT_T = mass_Rck_T.ldlt().solve(compute_gram_matrix(basis_Rck_T_quad, basis_Pk3_T_quad, quad_2k_T) * cellOperators(iT).gradient);

    uGT.block(offset, 0, PolynomialSpaceDimension<Cell>::Roly(degree()-1), dim_xgrad_T) = pi_Rkmo_GT_T;
    offset += PolynomialSpaceDimension<Cell>::Roly(degree()-1);
    uGT.block(offset, 0, PolynomialSpaceDimension<Cell>::RolyCompl(degree()), dim_xgrad_T) = pi_Rck_GT_T;
    offset += PolynomialSpaceDimension<Cell>::RolyCompl(degree());
  }
  
  return uGT;
}
