#include <xcurl.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

XCurl::XCurl(const DDRCore & ddr_core, bool use_threads, std::ostream & output)
  : DDRSpace(
	     ddr_core.mesh(),
	     0,
	     PolynomialSpaceDimension<Edge>::Poly(ddr_core.degree()),
	     PolynomialSpaceDimension<Face>::Roly(ddr_core.degree() - 1) + PolynomialSpaceDimension<Face>::RolyOrth(ddr_core.degree()),
	     PolynomialSpaceDimension<Cell>::Roly(ddr_core.degree() - 1) + PolynomialSpaceDimension<Cell>::RolyOrth(ddr_core.degree())	     
	     ),
    m_ddr_core(ddr_core),
    m_use_threads(use_threads),
    m_output(output),
    m_cell_operators(ddr_core.mesh().n_cells()),
    m_face_operators(ddr_core.mesh().n_faces())
{
  output << "[XCurl] Initializing" << std::endl;
  if (use_threads) {
    m_output << "[XCurl] Parallel execution" << std::endl;
  } else {
    m_output << "[XCurl] Sequential execution" << std::endl;
  }

  // Construct face curls and potentials
  std::function<void(size_t, size_t)> construct_all_face_curls_potentials
    = [this](size_t start, size_t end)->void
      {
        for (size_t iF = start; iF < end; iF++) {
          m_face_operators[iF].reset( new LocalOperators(_compute_face_curl_potential(iF)) );
        } // for iF
      };

  m_output << "[XCurl] Constructing face curls and potentials" << std::endl;
  parallel_for(mesh().n_faces(), construct_all_face_curls_potentials, use_threads);

  // Construct cell curls and potentials
  std::function<void(size_t, size_t)> construct_all_cell_curls_potentials
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_cell_operators[iT].reset( new LocalOperators(_compute_cell_curl_potential(iT)) );
        } // for iT
      };

  m_output << "[XCurl] Constructing cell curls and potentials" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_curls_potentials, use_threads);  
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd XCurl::interpolate(const FunctionType & v) const
{
  Eigen::VectorXd vh = Eigen::VectorXd::Zero(dimension());

  // Interpolate at edges
  std::function<void(size_t, size_t)> interpolate_edges
    = [this, &vh, v](size_t start, size_t end)->void
      {
	for (size_t iE = start; iE < end; iE++) {
	  const Edge & E = *mesh().edge(iE);

	  Eigen::Vector3d tE = E.tangent();
	  auto v_dot_tE = [&tE, v](const Eigen::Vector3d & x)->double {
			    return v(x).dot(tE);
			  };

	  QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2 * degree());
	  auto basis_Pk_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polyk, quad_2k_E);
	  vh.segment(globalOffset(E), edgeBases(iE).Polyk->dimension())
	    = l2_projection(v_dot_tE, *edgeBases(iE).Polyk, quad_2k_E, basis_Pk_E_quad);
	} // for iE
      };
  parallel_for(mesh().n_edges(), interpolate_edges, m_use_threads);

  if (degree() > 0 ) {
    // Interpolate at faces
    std::function<void(size_t, size_t)> interpolate_faces
      = [this, &vh, v](size_t start, size_t end)->void
	{
	  for (size_t iF = start; iF < end; iF++) {
	    const Face & F = *mesh().face(iF);

	    Eigen::Vector3d nF = F.normal();
	    auto nF_cross_v_cross_nF = [&nF, v](const Eigen::Vector3d & x)->Eigen::Vector3d {
					 return nF.cross(v(x)).cross(nF);
				       };

	    QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree());

	    size_t offset_F = globalOffset(F);
	    auto basis_Rkmo_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).Rolykmo, quad_2k_F);
	    vh.segment(offset_F, PolynomialSpaceDimension<Face>::Roly(degree() - 1))
	      = l2_projection(nF_cross_v_cross_nF, *faceBases(iF).Rolykmo, quad_2k_F, basis_Rkmo_F_quad);

	    offset_F += PolynomialSpaceDimension<Face>::Roly(degree() - 1);
	    auto basis_ROk_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).RolyOrthk, quad_2k_F);
	    vh.segment(offset_F, PolynomialSpaceDimension<Face>::RolyOrth(degree()))
	      = l2_projection(nF_cross_v_cross_nF, *faceBases(iF).RolyOrthk, quad_2k_F, basis_ROk_F_quad);
	  } // for iF
	};
    parallel_for(mesh().n_faces(), interpolate_faces, m_use_threads);

    // Interpolate at cells
    std::function<void(size_t, size_t)> interpolate_cells
      = [this, &vh, v](size_t start, size_t end)->void
	{
	  for (size_t iT = start; iT < end; iT++) {
	    const Cell & T = *mesh().cell(iT);

	    QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * degree());

	    Eigen::Index offset_T = globalOffset(T);
	    auto basis_Rkmo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Rolykmo, quad_2k_T);
	    vh.segment(offset_T, PolynomialSpaceDimension<Cell>::Roly(degree() - 1))
	      = l2_projection(v, *cellBases(iT).Rolykmo, quad_2k_T, basis_Rkmo_T_quad);

	    offset_T += PolynomialSpaceDimension<Cell>::Roly(degree() - 1);
	    auto basis_ROk_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).RolyOrthk, quad_2k_T);
	    vh.segment(offset_T, PolynomialSpaceDimension<Cell>::RolyOrth(degree()))
	      = l2_projection(v, *cellBases(iT).RolyOrthk, quad_2k_T, basis_ROk_T_quad);
	  } // for iT
	};
    parallel_for(mesh().n_cells(), interpolate_cells, m_use_threads);
  } // if degree() > 0
  
  return vh;
}

//------------------------------------------------------------------------------
// Curl and potential reconstruction
//------------------------------------------------------------------------------

XCurl::LocalOperators XCurl::_compute_face_curl_potential(size_t iF)
{
  const Face & F = *mesh().face(iF);
  
  //------------------------------------------------------------------------------
  // Curl
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix

  QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree());

  auto basis_Pk_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).Polyk, quad_2k_F);
  Eigen::MatrixXd MCF = compute_gram_matrix(basis_Pk_F_quad, quad_2k_F);
  
  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BCF
    = Eigen::MatrixXd::Zero(faceBases(iF).Polyk->dimension(), dimensionFace(iF));

  for (size_t iE = 0; iE < F.n_edges(); iE++) {
    const Edge & E = *F.edge(iE);
    QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2 * degree());
    BCF.block(0, localOffset(F, E), faceBases(iF).Polyk->dimension(), edgeBases(E.global_index()).Polyk->dimension())
      -= F.edge_orientation(iE) * compute_gram_matrix(
						      evaluate_quad<Function>::compute(*faceBases(iF).Polyk, quad_2k_E),
						      evaluate_quad<Function>::compute(*edgeBases(E.global_index()).Polyk, quad_2k_E),
						      quad_2k_E
						      );
  } // for iE

  if (degree() > 0) {
    QuadratureRule quad_2kmo_F = generate_quadrature_rule(F, 2 * (degree() - 1));

    BCF.block(0, localOffset(F), faceBases(iF).Polyk->dimension(), faceBases(iF).Rolykmo->dimension())
      += compute_gram_matrix(
			     evaluate_quad<Curl>::compute(*faceBases(iF).Polyk, quad_2kmo_F),
			     evaluate_quad<Function>::compute(*faceBases(iF).Rolykmo, quad_2kmo_F),
			     quad_2kmo_F
			     );
  } // if degree() > 0
 
  Eigen::MatrixXd CF = MCF.ldlt().solve(BCF);
  
  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------

  auto basis_Pkpo0_F = ShiftedBasis<typename DDRCore::PolyFaceBasisType>(*faceBases(iF).Polykpo, 1);

  Eigen::MatrixXd MPF
    = Eigen::MatrixXd::Zero(faceBases(iF).Polyk2->dimension(),
			    faceBases(iF).Polyk2->dimension());
  Eigen::MatrixXd BPF
    = Eigen::MatrixXd::Zero(faceBases(iF).Polyk2->dimension(), dimensionFace(iF));

  auto basis_Pk2_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).Polyk2, quad_2k_F);
  MPF.topLeftCorner(basis_Pkpo0_F.dimension(), faceBases(iF).Polyk2->dimension())
    = compute_gram_matrix(
			  evaluate_quad<Curl>::compute(basis_Pkpo0_F, quad_2k_F),
			  basis_Pk2_F_quad, quad_2k_F
			  );

  if (degree() > 0) {
    auto basis_ROk_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).RolyOrthk, quad_2k_F);
    MPF.bottomLeftCorner(faceBases(iF).RolyOrthk->dimension(), faceBases(iF).Polyk2->dimension())
      = compute_gram_matrix(basis_ROk_F_quad, basis_Pk2_F_quad, quad_2k_F);
    BPF.bottomRightCorner(faceBases(iF).RolyOrthk->dimension(), faceBases(iF).RolyOrthk->dimension())
      += compute_gram_matrix(basis_ROk_F_quad, quad_2k_F);    
  } // if degree() > 0
 
  auto quad_2kpo_F = generate_quadrature_rule(F, 2 * degree() + 1);
  BPF.topLeftCorner(basis_Pkpo0_F.dimension(), dimensionFace(iF))
    += compute_gram_matrix(
			   evaluate_quad<Function>::compute(basis_Pkpo0_F, quad_2kpo_F),
			   evaluate_quad<Function>::compute(*faceBases(iF).Polyk, quad_2kpo_F),
			   quad_2kpo_F
			   ) * CF;
 
  for (size_t iE = 0; iE < F.n_edges(); iE++) {
    const Edge & E = *F.edge(iE);
    QuadratureRule quad_2kpo_E = generate_quadrature_rule(E, 2 * degree() + 1);
    BPF.block(0, localOffset(F, E), basis_Pkpo0_F.dimension(), edgeBases(E.global_index()).Polyk->dimension())
      += F.edge_orientation(iE) * compute_gram_matrix(
						      evaluate_quad<Function>::compute(basis_Pkpo0_F, quad_2kpo_E),
						      evaluate_quad<Function>::compute(*edgeBases(E.global_index()).Polyk, quad_2kpo_E),
						      quad_2kpo_E
						      );    
  } // for iE

  return LocalOperators(CF, MPF.partialPivLu().solve(BPF));
}

//------------------------------------------------------------------------------

XCurl::LocalOperators XCurl::_compute_cell_curl_potential(size_t iT)
{
  const Cell & T = *mesh().cell(iT);

  //------------------------------------------------------------------------------
  // Curl
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix

  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * degree());

  auto basis_Pk3_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_T);
  Eigen::MatrixXd gram_Pk3_T = compute_gram_matrix(basis_Pk3_T_quad, quad_2k_T);
  Eigen::LDLT<Eigen::MatrixXd> ldlt_gram_Pk3_T(gram_Pk3_T);

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BCT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(), dimensionCell(iT));

  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    Eigen::Vector3d nF = F.normal();
    QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree());

    Eigen::MatrixXd BCT_F
      = T.face_orientation(iF) * compute_gram_matrix(
						     vector_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_F), nF),
						     evaluate_quad<Function>::compute(*faceBases(F.global_index()).Polyk2, quad_2k_F),
						     quad_2k_F
						     ) * faceOperators(F).potential;
    // Assemble local contribution
    for (size_t iE = 0; iE < F.n_edges(); iE++) {
      const Edge & E = *F.edge(iE);
      BCT.block(0, localOffset(T, E), cellBases(iT).Polyk3->dimension(), edgeBases(E.global_index()).Polyk->dimension())
	+= BCT_F.block(0, localOffset(F, E), cellBases(iT).Polyk3->dimension(), edgeBases(E.global_index()).Polyk->dimension());
    } // for iE
    BCT.block(0, localOffset(T, F), cellBases(iT).Polyk3->dimension(), numLocalDofsFace())
      += BCT_F.rightCols(numLocalDofsFace());
  } // for iF

  if (degree() > 0) {
    QuadratureRule quad_2kmo_T = generate_quadrature_rule(T, 2 * (degree() - 1));
    BCT.block(0, localOffset(T), cellBases(iT).Polyk3->dimension(), cellBases(iT).Rolykmo->dimension())
      += compute_gram_matrix(
			     evaluate_quad<Curl>::compute(*cellBases(iT).Polyk3, quad_2kmo_T),
			     evaluate_quad<Function>::compute(*cellBases(iT).Rolykmo, quad_2kmo_T),
			     quad_2kmo_T
			     );
  } // if degree() > 0 

  Eigen::MatrixXd CT = ldlt_gram_Pk3_T.solve(BCT);
  
  //------------------------------------------------------------------------------
  // Discrete curl-curl product
  //------------------------------------------------------------------------------

  // Consistent term
  Eigen::MatrixXd AT = CT.transpose() * BCT;

  // Stabilisation
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    Eigen::Vector3d nF = F.normal();
    QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree());
    auto basis_Pk3_T_dot_nF_quad
      = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_F), nF);
    auto basis_Pk_F_quad
      = evaluate_quad<Function>::compute(*faceBases(F.global_index()).Polyk, quad_2k_F);
    Eigen::MatrixXd M_Pk3_T_dot_nF_Pk_F
      = compute_gram_matrix(basis_Pk3_T_dot_nF_quad, basis_Pk_F_quad, quad_2k_F);
    Eigen::MatrixXd CF
      = extendOperator(T, F, m_face_operators[F.global_index()]->curl);
    AT +=
      CT.transpose() * compute_gram_matrix(basis_Pk3_T_dot_nF_quad, quad_2k_F) * CT
      - CT.transpose() * M_Pk3_T_dot_nF_Pk_F * CF
      - CF.transpose() * M_Pk3_T_dot_nF_Pk_F.transpose() * CT
      + CF.transpose() * compute_gram_matrix(basis_Pk_F_quad, quad_2k_F) * CF;
  } // for iF
  
  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------
  
  Eigen::MatrixXd MPT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(),
			    cellBases(iT).Polyk3->dimension());  
  
  Eigen::MatrixXd BPT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(), dimensionCell(iT));

  MPT.topRows(cellBases(iT).GolyOrthkpo->dimension())
    = compute_gram_matrix(
			  evaluate_quad<Curl>::compute(*cellBases(iT).GolyOrthkpo, quad_2k_T),
			  basis_Pk3_T_quad,
			  quad_2k_T
			  );
  
  if (degree() > 0) {
    auto basis_ROk_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).RolyOrthk, quad_2k_T);
    MPT.bottomRows(cellBases(iT).RolyOrthk->dimension())
      = compute_gram_matrix(basis_ROk_T_quad, basis_Pk3_T_quad, quad_2k_T);
    BPT.bottomRightCorner(cellBases(iT).RolyOrthk->dimension(), cellBases(iT).RolyOrthk->dimension())
      = compute_gram_matrix(basis_ROk_T_quad, quad_2k_T);
  } // if degree() > 0

  
  auto quad_2kpo_T = generate_quadrature_rule(T, 2 * degree() + 1);
  BPT.topRows(cellBases(iT).GolyOrthkpo->dimension())
    += compute_gram_matrix(
			   evaluate_quad<Function>::compute(*cellBases(iT).GolyOrthkpo, quad_2kpo_T),
			   evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2kpo_T),
			   quad_2kpo_T
			   ) * CT;

  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    Eigen::Vector3d nF = F.normal();
    QuadratureRule quad_2kpo_F = generate_quadrature_rule(F, 2 * degree() + 1);
    Eigen::MatrixXd PF
      = extendOperator(T, F, m_face_operators[F.global_index()]->potential);
    BPT.topRows(cellBases(iT).GolyOrthkpo->dimension())
      -= T.face_orientation(iF) * compute_gram_matrix(
						      vector_product(evaluate_quad<Function>::compute(*cellBases(iT).GolyOrthkpo, quad_2kpo_F), nF),
						      evaluate_quad<Function>::compute(*faceBases(F.global_index()).Polyk2, quad_2kpo_F),
						      quad_2kpo_F
						      ) * PF;
  } // for iF

  Eigen::MatrixXd PT = MPT.partialPivLu().solve(BPT);
  
  // Correction to enforce that the L2-orthogonal projection of PT
  // on Rk-1(T) is equal to the cell unknown
  if (degree() > 0) {
    QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * degree());
    auto basis_Rkmo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Rolykmo, quad_2k_T);

    Eigen::MatrixXd gram_Rkmo_T = compute_gram_matrix(basis_Rkmo_T_quad, quad_2k_T);
    Eigen::MatrixXd gram_Rkmo_T_Pk3_T = compute_gram_matrix(basis_Rkmo_T_quad, basis_Pk3_T_quad, quad_2k_T);
    Eigen::MatrixXd proj_Rkmo_T_Pk3_T = ldlt_gram_Pk3_T.solve(gram_Rkmo_T_Pk3_T.transpose());
							      
    // Remove the L2-orthogonal projection of PT on Rk-1(T) and replace
    // it with the cell unknown
    PT -= proj_Rkmo_T_Pk3_T * gram_Rkmo_T.ldlt().solve(gram_Rkmo_T_Pk3_T * PT);
    PT.middleCols(localOffset(T), PolynomialSpaceDimension<Cell>::Roly(degree() - 1)) += proj_Rkmo_T_Pk3_T;
  }

  return LocalOperators(CT, PT);
}

//------------------------------------------------------------------------------

Eigen::MatrixXd XCurl::computeL2Product(size_t iT, const double & penalty_factor)
{
  const Cell & T = *mesh().cell(iT); 
  const Eigen::MatrixXd & PT = cellOperators(iT).potential;

  Eigen::MatrixXd L2P = Eigen::MatrixXd::Zero(dimensionCell(iT), dimensionCell(iT));
  
  // Edge penalty terms
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    Eigen::VectorXd tE = E.tangent();
    double hE2 = std::pow(E.measure(), 2);

    // Shortcuts
    size_t dim_Pk_E = edgeBases(E.global_index()).Polyk->dimension();
    size_t offset_E = localOffset(T, E);
    
    QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2 * degree());
    auto basis_Pk3_T_dot_tE_quad = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_E), tE);
    auto basis_Pk_E_quad = evaluate_quad<Function>::compute(*edgeBases(E.global_index()).Polyk, quad_2k_E);
    Eigen::MatrixXd gram_Pk3_T_dot_tE_Pk_E = compute_gram_matrix(basis_Pk3_T_dot_tE_quad, basis_Pk_E_quad, quad_2k_E);
    
    L2P += hE2 * PT.transpose() * compute_gram_matrix(basis_Pk3_T_dot_tE_quad, quad_2k_E) * PT;
    L2P.middleCols(offset_E, dim_Pk_E) -= hE2 * PT.transpose() * gram_Pk3_T_dot_tE_Pk_E;
    L2P.middleRows(offset_E, dim_Pk_E) -= hE2 * gram_Pk3_T_dot_tE_Pk_E.transpose() * PT;     
    L2P.block(offset_E, offset_E, dim_Pk_E, dim_Pk_E) += hE2 * compute_gram_matrix(basis_Pk_E_quad, quad_2k_E);
  } // for iE

  if (degree() > 0) {
    // Face penalty terms
    for (size_t iF = 0; iF < T.n_faces(); iF++) {
      const Face & F = *T.face(iF);
      double hF = F.diam();

      size_t offset_F = localOffset(T, F);

      QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree());
      auto basis_Rkmo_F_quad = evaluate_quad<Function>::compute(*faceBases(F.global_index()).Rolykmo, quad_2k_F);
      auto basis_ROk_F_quad = evaluate_quad<Function>::compute(*faceBases(F.global_index()).RolyOrthk, quad_2k_F);
      auto basis_Pk3_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_F);

      size_t dim_Rkmo_F = PolynomialSpaceDimension<Face>::Roly(degree() - 1);
      Eigen::MatrixXd gram_Rkmo_F = compute_gram_matrix(basis_Rkmo_F_quad, quad_2k_F);
      Eigen::MatrixXd gram_Rkmo_F_PT = compute_gram_matrix(basis_Rkmo_F_quad, basis_Pk3_T_quad, quad_2k_F) * PT;
      L2P += hF * gram_Rkmo_F_PT.transpose() * gram_Rkmo_F.ldlt().solve(gram_Rkmo_F_PT);
      L2P.middleCols(offset_F, dim_Rkmo_F) -= hF * gram_Rkmo_F_PT.transpose();
      L2P.middleRows(offset_F, dim_Rkmo_F) -= hF * gram_Rkmo_F_PT;
      L2P.block(offset_F, offset_F, dim_Rkmo_F, dim_Rkmo_F) += hF * gram_Rkmo_F;

      size_t dim_ROk_F = PolynomialSpaceDimension<Face>::RolyOrth(degree());
      offset_F += dim_Rkmo_F;
      Eigen::MatrixXd gram_ROk_F = compute_gram_matrix(basis_ROk_F_quad, quad_2k_F);
      Eigen::MatrixXd gram_ROk_F_PT = compute_gram_matrix(basis_ROk_F_quad, basis_Pk3_T_quad, quad_2k_F) * PT;
      L2P += hF * gram_ROk_F_PT.transpose() * gram_ROk_F.ldlt().solve(gram_ROk_F_PT);
      L2P.middleCols(offset_F, dim_ROk_F) -= hF * gram_ROk_F_PT.transpose();
      L2P.middleRows(offset_F, dim_ROk_F) -= hF * gram_ROk_F_PT;
      L2P.block(offset_F, offset_F, dim_ROk_F, dim_ROk_F) += hF * gram_ROk_F;
    } // for iF
  } // if degree() > 0

  L2P *= penalty_factor;
  
    // Consistent term
  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * degree());
  L2P += PT.transpose() * compute_gram_matrix(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_T), quad_2k_T) * PT;

  return L2P;
}
