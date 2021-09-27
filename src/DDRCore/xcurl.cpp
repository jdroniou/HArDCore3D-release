#include <xcurl.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_face.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

XCurl::XCurl(const DDRCore & ddr_core, bool use_threads, std::ostream & output)
  : DDRSpace(
	     ddr_core.mesh(),
	     0,
	     PolynomialSpaceDimension<Edge>::Poly(ddr_core.degree()),
	     PolynomialSpaceDimension<Face>::Roly(ddr_core.degree() - 1) + PolynomialSpaceDimension<Face>::RolyCompl(ddr_core.degree()),
	     PolynomialSpaceDimension<Cell>::Roly(ddr_core.degree() - 1) + PolynomialSpaceDimension<Cell>::RolyCompl(ddr_core.degree())	     
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

Eigen::VectorXd XCurl::interpolate(const FunctionType & v, const int doe_cell, const int doe_face, const int doe_edge) const
{
  Eigen::VectorXd vh = Eigen::VectorXd::Zero(dimension());
  
  // Degrees of quadrature rules
  size_t dqr_cell = (doe_cell >= 0 ? doe_cell : 2 * degree() + 3);
  size_t dqr_face = (doe_face >= 0 ? doe_face : 2 * degree() + 3);
  size_t dqr_edge = (doe_edge >= 0 ? doe_edge : 2 * degree() + 3);
  
  // Interpolate at edges
  std::function<void(size_t, size_t)> interpolate_edges
    = [this, &vh, v, &dqr_edge](size_t start, size_t end)->void
      {
	      for (size_t iE = start; iE < end; iE++) {
	        const Edge & E = *mesh().edge(iE);

	        Eigen::Vector3d tE = E.tangent();
	        auto v_dot_tE = [&tE, v](const Eigen::Vector3d & x)->double {
			          return v(x).dot(tE);
			        };

	        QuadratureRule quad_dqr_E = generate_quadrature_rule(E, dqr_edge);
	        auto basis_Pk_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polyk, quad_dqr_E);
	        vh.segment(globalOffset(E), edgeBases(iE).Polyk->dimension())
	          = l2_projection(v_dot_tE, *edgeBases(iE).Polyk, quad_dqr_E, basis_Pk_E_quad);
	      } // for iE
      };
  parallel_for(mesh().n_edges(), interpolate_edges, m_use_threads);

  if (degree() > 0 ) {
    // Interpolate at faces
    std::function<void(size_t, size_t)> interpolate_faces
      = [this, &vh, v, &dqr_face](size_t start, size_t end)->void
	      {
	        for (size_t iF = start; iF < end; iF++) {
	          const Face & F = *mesh().face(iF);

	          Eigen::Vector3d nF = F.normal();
	          auto nF_cross_v_cross_nF = [&nF, v](const Eigen::Vector3d & x)->Eigen::Vector3d {
                                         return nF.cross(v(x).cross(nF));
                                       };

	          QuadratureRule quad_dqr_F = generate_quadrature_rule(F, dqr_face);

	          size_t offset_F = globalOffset(F);
	          auto basis_Rkmo_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).Rolykmo, quad_dqr_F);
	          vh.segment(offset_F, PolynomialSpaceDimension<Face>::Roly(degree() - 1))
	            = l2_projection(nF_cross_v_cross_nF, *faceBases(iF).Rolykmo, quad_dqr_F, basis_Rkmo_F_quad);

	          offset_F += PolynomialSpaceDimension<Face>::Roly(degree() - 1);
	          auto basis_Rck_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).RolyComplk, quad_dqr_F);
	          vh.segment(offset_F, PolynomialSpaceDimension<Face>::RolyCompl(degree()))
	            = l2_projection(nF_cross_v_cross_nF, *faceBases(iF).RolyComplk, quad_dqr_F, basis_Rck_F_quad);
	        } // for iF
	      };
    parallel_for(mesh().n_faces(), interpolate_faces, m_use_threads);

    // Interpolate at cells
    std::function<void(size_t, size_t)> interpolate_cells
      = [this, &vh, v, &dqr_cell](size_t start, size_t end)->void
	      {
	        for (size_t iT = start; iT < end; iT++) {
	          const Cell & T = *mesh().cell(iT);

	          QuadratureRule quad_dqr_T = generate_quadrature_rule(T, dqr_cell);

	          Eigen::Index offset_T = globalOffset(T);
	          auto basis_Rkmo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Rolykmo, quad_dqr_T);
	          vh.segment(offset_T, PolynomialSpaceDimension<Cell>::Roly(degree() - 1))
	            = l2_projection(v, *cellBases(iT).Rolykmo, quad_dqr_T, basis_Rkmo_T_quad, GramMatrix(T, *cellBases(iT).Rolykmo));

	          offset_T += PolynomialSpaceDimension<Cell>::Roly(degree() - 1);
	          auto basis_Rck_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).RolyComplk, quad_dqr_T);
	          vh.segment(offset_T, PolynomialSpaceDimension<Cell>::RolyCompl(degree()))
	            = l2_projection(v, *cellBases(iT).RolyComplk, quad_dqr_T, basis_Rck_T_quad, GramMatrix(T, *cellBases(iT).RolyComplk));
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

  MonomialFaceIntegralsType int_monoF_2kpo = IntegrateFaceMonomials(F, 2*degree()+1);
  Eigen::MatrixXd MCF = GramMatrix(F, *faceBases(iF).Polyk, int_monoF_2kpo);
  
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
    CurlBasis<RestrictedBasis<DDRCore::PolyBasisFaceType>> rot_PkF(*faceBases(iF).Polyk);
    BCF.block(0, localOffset(F), faceBases(iF).Polyk->dimension(), faceBases(iF).Rolykmo->dimension())
      += GramMatrix(F, rot_PkF, *faceBases(iF).Rolykmo);
  } // if degree() > 0
 
  Eigen::MatrixXd CF = MCF.ldlt().solve(BCF);
  
  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------

  auto basis_Pkpo0_F = ShiftedBasis<typename DDRCore::PolyBasisFaceType>(*faceBases(iF).Polykpo, 1);

  Eigen::MatrixXd MPF
    = Eigen::MatrixXd::Zero(faceBases(iF).Polyk2->dimension(),
			    faceBases(iF).Polyk2->dimension());
  Eigen::MatrixXd BPF
    = Eigen::MatrixXd::Zero(faceBases(iF).Polyk2->dimension(), dimensionFace(iF));

  CurlBasis<decltype(basis_Pkpo0_F)> rot_Pkp0_F(basis_Pkpo0_F);
  MPF.topLeftCorner(basis_Pkpo0_F.dimension(), faceBases(iF).Polyk2->dimension())
    = GramMatrix(F, rot_Pkp0_F, *faceBases(iF).Polyk2, int_monoF_2kpo);

  if (degree() > 0) {
    MPF.bottomLeftCorner(faceBases(iF).RolyComplk->dimension(), faceBases(iF).Polyk2->dimension())
      = GramMatrix(F, *faceBases(iF).RolyComplk, *faceBases(iF).Polyk2, int_monoF_2kpo);
    BPF.bottomRightCorner(faceBases(iF).RolyComplk->dimension(), faceBases(iF).RolyComplk->dimension())
      += GramMatrix(F, *faceBases(iF).RolyComplk, int_monoF_2kpo);    
  } // if degree() > 0
 
  BPF.topLeftCorner(basis_Pkpo0_F.dimension(), dimensionFace(iF))
    += GramMatrix(F, basis_Pkpo0_F, *faceBases(iF).Polyk, int_monoF_2kpo) * CF;
 
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

  // Compute all integrals of monomial powers to degree 2k and the mass matrix
  MonomialIntegralsType int_mono_2kpo = IntegrateCellMonomials(T, 2*degree()+1);
  Eigen::MatrixXd gram_Pk3_T = GramMatrix(T, *cellBases(iT).Polyk3, int_mono_2kpo);
  Eigen::LDLT<Eigen::MatrixXd> ldlt_gram_Pk3_T(gram_Pk3_T);

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BCT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(), dimensionCell(iT));

  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    Eigen::Vector3d nF = F.normal();

    MonomialScalarBasisFace basis_Pk_F(F, degree());
    DecomposePoly dec(F, TangentFamily<MonomialScalarBasisFace>(basis_Pk_F, basis_Pk_F.coordinates_system()));
    boost::multi_array<VectorRd, 2> Pk3_T_cross_nF_nodes 
          = vector_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, dec.get_nodes()), nF);
    auto Pk3_T_cross_nF_family = dec.family(Pk3_T_cross_nF_nodes);
    Eigen::MatrixXd PF = extendOperator(T, F, m_face_operators[F.global_index()]->potential);
    BCT += T.face_orientation(iF) * GramMatrix(F, Pk3_T_cross_nF_family, *faceBases(F).Polyk2) * PF;

    // The following commented block replaces the previous one, without using DecomposePoly
    /*
      QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree());
      Eigen::MatrixXd PF = extendOperator(T, F, m_face_operators[F.global_index()]->potential);
      BCT += T.face_orientation(iF) * compute_gram_matrix(
					    vector_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_F), nF),
			        evaluate_quad<Function>::compute(*faceBases(F).Polyk2, quad_2k_F),
					    quad_2k_F
					    ) * PF;
    */
    
  } // for iF

  if (degree() > 0) {
    CurlBasis<DDRCore::Poly3BasisCellType> curl_Pk3_basis(*cellBases(iT).Polyk3);
    BCT.block(0, localOffset(T), cellBases(iT).Polyk3->dimension(), cellBases(iT).Rolykmo->dimension())
      += GramMatrix(T, curl_Pk3_basis, *cellBases(iT).Rolykmo, int_mono_2kpo);
  } // if degree() > 0 

  Eigen::MatrixXd CT = ldlt_gram_Pk3_T.solve(BCT);
  
  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------
  
  Eigen::MatrixXd MPT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(),
			    cellBases(iT).Polyk3->dimension());  
  
  Eigen::MatrixXd BPT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(), dimensionCell(iT));

  CurlBasis<DDRCore::GolyComplBasisCellType> Rolyk_basis(*cellBases(iT).GolyComplkpo);
  MPT.topRows(cellBases(iT).GolyComplkpo->dimension()) = GramMatrix(T, Rolyk_basis, *cellBases(iT).Polyk3, int_mono_2kpo);
  
  if (degree() > 0) {
    MPT.bottomRows(cellBases(iT).RolyComplk->dimension()) = GramMatrix(T, *cellBases(iT).RolyComplk, *cellBases(iT).Polyk3, int_mono_2kpo);
    BPT.bottomRightCorner(cellBases(iT).RolyComplk->dimension(), cellBases(iT).RolyComplk->dimension())
      = GramMatrix(T, *cellBases(iT).RolyComplk, int_mono_2kpo);
  } // if degree() > 0

  
  BPT.topRows(cellBases(iT).GolyComplkpo->dimension())
    += GramMatrix(T, *cellBases(iT).GolyComplkpo, *cellBases(iT).Polyk3, int_mono_2kpo) * CT;
    
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    Eigen::Vector3d nF = F.normal();
    Eigen::MatrixXd PF = extendOperator(T, F, m_face_operators[F.global_index()]->potential);
  
    MonomialScalarBasisFace basis_Pkpo_F(F, degree()+1);
    DecomposePoly dec(F, TangentFamily<MonomialScalarBasisFace>(basis_Pkpo_F, basis_Pkpo_F.coordinates_system()));
    boost::multi_array<VectorRd, 2> Gkpo_T_cross_nF_nodes 
            = vector_product(evaluate_quad<Function>::compute(*cellBases(iT).GolyComplkpo, dec.get_nodes()), nF);
    auto Gkpo_T_cross_nF_family = dec.family(Gkpo_T_cross_nF_nodes);
    BPT.topRows(cellBases(iT).GolyComplkpo->dimension())
      -= T.face_orientation(iF) * GramMatrix(F, Gkpo_T_cross_nF_family, *faceBases(F).Polyk2) * PF;

    // The following commented block replaces the previous one, without DecomposePoly
    /*
      QuadratureRule quad_2kpo_F = generate_quadrature_rule(F, 2 * degree() + 1);
      BPT.topRows(cellBases(iT).GolyComplkpo->dimension())
        -= T.face_orientation(iF) * compute_gram_matrix(
						      vector_product(evaluate_quad<Function>::compute(*cellBases(iT).GolyComplkpo, quad_2kpo_F), nF),
						      evaluate_quad<Function>::compute(*faceBases(F).Polyk2, quad_2kpo_F),
						      quad_2kpo_F
						      ) * PF;
    */

  } // for iF
  
  return LocalOperators(CT, MPT.partialPivLu().solve(BPT));
}

//------------------------------------------------------------------------------
//        Functions to compute matrices for local L2 products on Xcurl
//------------------------------------------------------------------------------

Eigen::MatrixXd XCurl::computeL2Product(
                                        const size_t iT,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pk3_T,
                                        const IntegralWeight & weight
                                        ) const
{
  const Cell & T = *mesh().cell(iT); 
  
  // create the weighted mass matrix, with simple product if weight is constant
  Eigen::MatrixXd w_mass_Pk3_T;
  if (weight.deg(T)==0){
    // constant weight
    if (mass_Pk3_T.rows()==1){
      // We have to compute the mass matrix
      w_mass_Pk3_T = weight.value(T, T.center_mass()) * GramMatrix(T, *cellBases(iT).Polyk3);
    }else{
      w_mass_Pk3_T = weight.value(T, T.center_mass()) * mass_Pk3_T;
    }
  }else{
    // weight is not constant, we create a weighted mass matrix
    QuadratureRule quad_2kpw_T = generate_quadrature_rule(T, 2 * degree() + weight.deg(T));
    auto basis_Pk3_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2kpw_T);
    std::function<double(const Eigen::Vector3d &)> weight_T 
              = [&T, &weight](const Eigen::Vector3d &x)->double {
                  return weight.value(T, x);
                };
    w_mass_Pk3_T = compute_weighted_gram_matrix(weight_T, basis_Pk3_T_quad, basis_Pk3_T_quad, quad_2kpw_T, "sym");
  }

  
  // The leftOp and rightOp will come from the potentials
  std::vector<Eigen::MatrixXd> potentialOp(T.n_edges()+T.n_faces()+1);
  for (size_t iE = 0; iE < T.n_edges(); iE++){
    const Edge & E = *T.edge(iE);
    potentialOp[iE] = extendOperator(T, E, Eigen::MatrixXd::Identity(dimensionEdge(E),dimensionEdge(E)));
  }
  for (size_t iF = 0; iF < T.n_faces(); iF++){
    const Face & F = *T.face(iF);
    potentialOp[T.n_edges()+iF] = extendOperator(T, F, m_face_operators[F.global_index()]->potential);
  }
  potentialOp[T.n_edges()+T.n_faces()] = m_cell_operators[iT]->potential;


  return computeL2Product_with_Ops(iT, potentialOp, potentialOp, penalty_factor, w_mass_Pk3_T, weight);

}

Eigen::MatrixXd XCurl::computeL2ProductGradient(
                                        const size_t iT,
                                        const XGrad & x_grad,
                                        const std::string & side,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pk3_T,
                                        const IntegralWeight & weight
                                        ) const
{
  const Cell & T = *mesh().cell(iT);

  // create the weighted mass matrix, with simple product if weight is constant
  Eigen::MatrixXd w_mass_Pk3_T;
  if (weight.deg(T)==0){
    // constant weight
    if (mass_Pk3_T.rows()==1){
      // We have to compute the mass matrix
      w_mass_Pk3_T = weight.value(T, T.center_mass()) * GramMatrix(T, *cellBases(iT).Polyk3);
    }else{
      w_mass_Pk3_T = weight.value(T, T.center_mass()) * mass_Pk3_T;
    }
  }else{
    // weight is not constant, we create a weighted mass matrix
    QuadratureRule quad_2kpw_T = generate_quadrature_rule(T, 2 * degree() + weight.deg(T));
    auto basis_Pk3_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2kpw_T);
    std::function<double(const Eigen::Vector3d &)> weight_T 
              = [&T, &weight](const Eigen::Vector3d &x)->double {
                  return weight.value(T, x);
                };
    w_mass_Pk3_T = compute_weighted_gram_matrix(weight_T, basis_Pk3_T_quad, basis_Pk3_T_quad, quad_2kpw_T, "sym");
  }


  // list of full gradients
  std::vector<Eigen::MatrixXd> gradientOp(T.n_edges()+T.n_faces()+1);
  for (size_t iE = 0; iE < T.n_edges(); iE++){
    const Edge & E = *T.edge(iE);
    gradientOp[iE] = x_grad.extendOperator(T, E, x_grad.edgeOperators(E).gradient);
  }
  for (size_t iF = 0; iF < T.n_faces(); iF++){
    const Face & F = *T.face(iF);
    gradientOp[T.n_edges()+iF] = x_grad.extendOperator(T, F, x_grad.faceOperators(F).gradient);
  }
  gradientOp[T.n_edges()+T.n_faces()] = x_grad.cellOperators(iT).gradient;
  
  // If we apply the gradient on one side only we'll need the potentials
  if (side != "both"){
    // list of potentials
    std::vector<Eigen::MatrixXd> potentialOp(T.n_edges()+T.n_faces()+1);
    for (size_t iE = 0; iE < T.n_edges(); iE++){
      const Edge & E = *T.edge(iE);
      potentialOp[iE] = extendOperator(T, E, Eigen::MatrixXd::Identity(dimensionEdge(E),dimensionEdge(E)));
    }
    for (size_t iF = 0; iF < T.n_faces(); iF++){
      const Face & F = *T.face(iF);
      potentialOp[T.n_edges()+iF] = extendOperator(T, F, m_face_operators[F.global_index()]->potential);
    }
    potentialOp[T.n_edges()+T.n_faces()] = m_cell_operators[iT]->potential;
  
    // Depending on side of gradient
    if (side == "left"){
      return computeL2Product_with_Ops(iT, gradientOp, potentialOp, penalty_factor, w_mass_Pk3_T, weight);
    }else{
      return computeL2Product_with_Ops(iT, potentialOp, gradientOp, penalty_factor, w_mass_Pk3_T, weight);
    }
    
  }

  // Default: gradient on both sides
  return computeL2Product_with_Ops(iT, gradientOp, gradientOp, penalty_factor, w_mass_Pk3_T, weight);

}


Eigen::MatrixXd XCurl::computeL2Product_with_Ops(
                                        const size_t iT,
                                        const std::vector<Eigen::MatrixXd> & leftOp,
                                        const std::vector<Eigen::MatrixXd> & rightOp,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & w_mass_Pk3_T,
                                        const IntegralWeight & weight
                                        ) const
{
  const Cell & T = *mesh().cell(iT); 

  // leftOp and rightOp must list the operators acting on the DOFs, and which we want to
  // use for the L2 product. Specifically, each one lists operators (matrices) returning
  // values in edges space P^k(E), faces space P^k(F)^2 (tangent) and element space P^k(T)^3.
  // For the standard Xcurl L2 product, these will respectively be identity (for each edge),
  // gamma_tF (for each F) and PT. 
  // To compute the Xcurl L2 product applied (left or right) to the discrete gradient,
  // leftOp or rightOp must list the edge, face and element (full) gradient operators.
  // All these operators must have the same domain, so possibly being extended appropriately
  // using extendOperator from ddrspace.

  Eigen::MatrixXd L2P = Eigen::MatrixXd::Zero(leftOp[0].cols(), rightOp[0].cols());
  
  size_t offset_F = T.n_edges();
  size_t offset_T = T.n_edges() + T.n_faces();

  // Edge penalty terms
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    VectorRd tE = E.tangent();
        
    QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2 * degree());
    
    // weight and scaling hE^2
    double max_weight_quad_E = weight.value(T, quad_2k_E[0].vector());
    // If the weight is not constant, we want to take the largest along the edge
    if (weight.deg(T)>0){
      for (size_t iqn = 1; iqn < quad_2k_E.size(); iqn++) {
        max_weight_quad_E = std::max(max_weight_quad_E, weight.value(T, quad_2k_E[iqn].vector()));
      } // for
    }
    double w_hE2 = max_weight_quad_E * std::pow(E.measure(), 2);

    // The penalty term int_E (PT w . tE - w_E) * (PT v . tE - v_E) is computed by developping.
    auto basis_Pk3_T_dot_tE_quad = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_E), tE);
    auto basis_Pk_E_quad = evaluate_quad<Function>::compute(*edgeBases(E.global_index()).Polyk, quad_2k_E);
    Eigen::MatrixXd gram_Pk3T_dot_tE_PkE = compute_gram_matrix(basis_Pk3_T_dot_tE_quad, basis_Pk_E_quad, quad_2k_E);
    
    // Contribution of edge E
    L2P += w_hE2 * ( leftOp[offset_T].transpose() * compute_gram_matrix(basis_Pk3_T_dot_tE_quad, quad_2k_E) * rightOp[offset_T]
                   - leftOp[offset_T].transpose() * gram_Pk3T_dot_tE_PkE * rightOp[iE]
                   - leftOp[iE].transpose() * gram_Pk3T_dot_tE_PkE.transpose() * rightOp[offset_T]
                   + leftOp[iE].transpose() * compute_gram_matrix(basis_Pk_E_quad, quad_2k_E) * rightOp[iE]);

  } // for iE

  // Face penalty terms
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree());

    // Compute mass-matrices: Polyk2-Polyk2, Polyk2-Polyk3 (also serves to take the tangential component of PT)
    MonomialScalarBasisFace basis_Pk_F(F, degree());
    DecomposePoly dec(F, TangentFamily<MonomialScalarBasisFace>(basis_Pk_F, basis_Pk_F.coordinates_system()));
    const VectorRd nF = F.normal();
    // Values of Pk3T at the nodes and projected on the tangent space to F
    auto Pk3T_tangent_nodes = transform_values_quad<VectorRd>(
                                evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, dec.get_nodes()),
                                [&nF](const VectorRd &z)->VectorRd { return z-(z.dot(nF))*nF;});
    auto Pk3T_tangent_family_Pk2F = dec.family(Pk3T_tangent_nodes);
    MonomialFaceIntegralsType int_monoF_2k = IntegrateFaceMonomials(F, 2*degree());
    Eigen::MatrixXd mass_Pk2F_Pk2F = GramMatrix(F, *faceBases(F).Polyk2, int_monoF_2k);
    Eigen::MatrixXd gram_Pk2F_Pk3T = GramMatrix(F, *faceBases(F).Polyk2, Pk3T_tangent_family_Pk2F, int_monoF_2k);

    // This commented block does the same as above, without DecomposePoly
    /*
      MonomialFaceIntegralsType int_monoF_2k = IntegrateFaceMonomials(F, 2*degree());
      Eigen::MatrixXd mass_Pk2F_Pk2F = GramMatrix(F, *faceBases(F).Polyk2, int_monoF_2k);
      auto basis_Pk2_F_quad = evaluate_quad<Function>::compute(*faceBases(F.global_index()).Polyk2, quad_2k_F);
      auto basis_Pk3_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_F);
      Eigen::MatrixXd gram_Pk2F_Pk3T = compute_gram_matrix(basis_Pk2_F_quad, basis_Pk3_T_quad, quad_2k_F, "nonsym");
    */
    
    // Weight coming from permeability and scaling hF
    double max_weight_quad_F = weight.value(T, quad_2k_F[0].vector());
    // If the weight is not constant, we want to take the largest along the edge
    if (weight.deg(T)>0){
      for (size_t iqn = 1; iqn < quad_2k_F.size(); iqn++) {
        max_weight_quad_F = std::max(max_weight_quad_F, weight.value(T, quad_2k_F[iqn].vector()));
      } // for
    }
    double w_hF = max_weight_quad_F * F.diam();

    // The penalty term int_T ((leftOp)_{tF}-leftOp_{tF}) * ((rightOp)_{tF}-rightOp_{tF}) is computed by developping,
    // and using the fact that (mass_Pk2F_Pk2F)^{-1} * gram_Pk2F_Pk3T is the matrix that represents the projection from
    // P^k(T)^3 to the tangent space P^k(F)^2
 
    // Contribution of face F
    L2P += w_hF * ( leftOp[offset_T].transpose() * gram_Pk2F_Pk3T.transpose() * mass_Pk2F_Pk2F.ldlt().solve(gram_Pk2F_Pk3T) * rightOp[offset_T]
                - leftOp[offset_T].transpose() * gram_Pk2F_Pk3T.transpose() * rightOp[offset_F+iF] 
                - leftOp[offset_F+iF].transpose() * gram_Pk2F_Pk3T * rightOp[offset_T]
                + leftOp[offset_F+iF].transpose() * mass_Pk2F_Pk2F * rightOp[offset_F+iF]                  
                  );

  } // for iF

  L2P *= penalty_factor;
  
  // Consistent (cell) term
  L2P += leftOp[offset_T].transpose() * w_mass_Pk3_T * rightOp[offset_T];
  
   
  return L2P;
}



//------------------------------------------------------------------------------

Eigen::MatrixXd XCurl::computeL2ProductDOFs(
                                        size_t iT,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pk3_T,
                                        const IntegralWeight & weight
                                        )
{
  const Cell & T = *mesh().cell(iT); 
  const Eigen::MatrixXd & PT = cellOperators(iT).potential;

  Eigen::MatrixXd L2P = Eigen::MatrixXd::Zero(dimensionCell(iT), dimensionCell(iT));
  
  // Edge penalty terms
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    VectorRd tE = E.tangent();

    // Shortcuts
    size_t dim_Pk_E = edgeBases(E.global_index()).Polyk->dimension();
    size_t offset_E = localOffset(T, E);
    
    QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2 * degree());
    
    double max_weight_quad_E = weight.value(T, quad_2k_E[0].vector());
    // If the weight is not constant, we want to take the largest along the edge
    if (weight.deg(T)>0){
      for (size_t iqn = 1; iqn < quad_2k_E.size(); iqn++) {
        max_weight_quad_E = std::max(max_weight_quad_E, weight.value(T, quad_2k_E[iqn].vector()));
      } // for
    }

    double w_hE2 = max_weight_quad_E * std::pow(E.measure(), 2);
                                      
    auto basis_Pk3_T_dot_tE_quad = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_E), tE);
    auto basis_Pk_E_quad = evaluate_quad<Function>::compute(*edgeBases(E.global_index()).Polyk, quad_2k_E);
    Eigen::MatrixXd gram_Pk3_T_dot_tE_Pk_E = compute_gram_matrix(basis_Pk3_T_dot_tE_quad, basis_Pk_E_quad, quad_2k_E);
    
    L2P += w_hE2 * PT.transpose() * compute_gram_matrix(basis_Pk3_T_dot_tE_quad, quad_2k_E) * PT;
    L2P.middleCols(offset_E, dim_Pk_E) -= w_hE2 * PT.transpose() * gram_Pk3_T_dot_tE_Pk_E;
    L2P.middleRows(offset_E, dim_Pk_E) -= w_hE2 * gram_Pk3_T_dot_tE_Pk_E.transpose() * PT;     
    L2P.block(offset_E, offset_E, dim_Pk_E, dim_Pk_E) += w_hE2 * compute_gram_matrix(basis_Pk_E_quad, quad_2k_E);
  } // for iE

  if (degree() > 0) {
    // Face penalty terms
    for (size_t iF = 0; iF < T.n_faces(); iF++) {
      const Face & F = *T.face(iF);

      size_t offset_F = localOffset(T, F);

      QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree());

      double max_weight_quad_F = weight.value(T, quad_2k_F[0].vector());
      // If the weight is not constant, we want to take the largest along the edge
      if (weight.deg(T)>0){
        for (size_t iqn = 1; iqn < quad_2k_F.size(); iqn++) {
          max_weight_quad_F = std::max(max_weight_quad_F, weight.value(T, quad_2k_F[iqn].vector()));
        } // for
      }
      
      MonomialScalarBasisFace basis_Pk_F(F, degree());
      DecomposePoly dec(F, TangentFamily<MonomialScalarBasisFace>(basis_Pk_F, basis_Pk_F.coordinates_system()));
      const VectorRd nF = F.normal();
      // Values of Pk3T at the nodes and projected on the tangent space to F
      auto Pk3T_tangent_nodes = transform_values_quad<VectorRd>(
                                  evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, dec.get_nodes()),
                                  [&nF](const VectorRd &z)->VectorRd { return z-(z.dot(nF))*nF;});
      auto Pk3T_tangent_family_Pk2F = dec.family(Pk3T_tangent_nodes);
      MonomialFaceIntegralsType int_monoF_2k = IntegrateFaceMonomials(F, 2*degree());
      Eigen::MatrixXd gram_Rkmo_F = GramMatrix(F, *faceBases(F).Rolykmo, int_monoF_2k);
      Eigen::MatrixXd gram_Rkmo_F_PT = GramMatrix(F, *faceBases(F).Rolykmo, Pk3T_tangent_family_Pk2F, int_monoF_2k) * PT;
      Eigen::MatrixXd gram_Rck_F_PT = GramMatrix(F, *faceBases(F).RolyComplk, Pk3T_tangent_family_Pk2F, int_monoF_2k) * PT;

      // This commented block does the same as above, without using DecomposePoly
      /*
        auto basis_Rkmo_F_quad = evaluate_quad<Function>::compute(*faceBases(F.global_index()).Rolykmo, quad_2k_F);
        auto basis_Rck_F_quad = evaluate_quad<Function>::compute(*faceBases(F.global_index()).RolyComplk, quad_2k_F);
        auto basis_Pk3_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_F);
        Eigen::MatrixXd gram_Rkmo_F = compute_gram_matrix(basis_Rkmo_F_quad, quad_2k_F);
        Eigen::MatrixXd gram_Rkmo_F_PT = compute_gram_matrix(basis_Rkmo_F_quad, basis_Pk3_T_quad, quad_2k_F) * PT;
        Eigen::MatrixXd gram_Rck_F_PT = compute_gram_matrix(basis_Rck_F_quad, basis_Pk3_T_quad, quad_2k_F) * PT;
      */
      
      double w_hF = max_weight_quad_F * F.diam();

      size_t dim_Rkmo_F = PolynomialSpaceDimension<Face>::Roly(degree() - 1);
      L2P += w_hF * gram_Rkmo_F_PT.transpose() * gram_Rkmo_F.ldlt().solve(gram_Rkmo_F_PT);
      L2P.middleCols(offset_F, dim_Rkmo_F) -= w_hF * gram_Rkmo_F_PT.transpose();
      L2P.middleRows(offset_F, dim_Rkmo_F) -= w_hF * gram_Rkmo_F_PT;
      L2P.block(offset_F, offset_F, dim_Rkmo_F, dim_Rkmo_F) += w_hF * gram_Rkmo_F;

      size_t dim_Rck_F = PolynomialSpaceDimension<Face>::RolyCompl(degree());
      offset_F += dim_Rkmo_F;
      Eigen::MatrixXd gram_Rck_F = GramMatrix(F, *faceBases(F.global_index()).RolyComplk);
      L2P += w_hF * gram_Rck_F_PT.transpose() * gram_Rck_F.ldlt().solve(gram_Rck_F_PT);
      L2P.middleCols(offset_F, dim_Rck_F) -= w_hF * gram_Rck_F_PT.transpose();
      L2P.middleRows(offset_F, dim_Rck_F) -= w_hF * gram_Rck_F_PT;
      L2P.block(offset_F, offset_F, dim_Rck_F, dim_Rck_F) += w_hF * gram_Rck_F;
    } // for iF
  } // if degree() > 0

  L2P *= penalty_factor;
  
  // Consistent term, two calculations depending if weight is constant or not.
  if (weight.deg(T)==0){
    // weight is constant so the calculation is based on the standard mass matrix multiplied by the weight. 
    Eigen::MatrixXd M_Pk3_T = mass_Pk3_T;
    if (M_Pk3_T.rows()==1){
      // Mass matrix not passed as parameter, we compute it
      M_Pk3_T.resize(cellBases(iT).Polyk3->dimension(), cellBases(iT).Polyk3->dimension());
      M_Pk3_T = GramMatrix(T, *cellBases(iT).Polyk3);
    }
    L2P += weight.value(T, T.center_mass()) * PT.transpose() * M_Pk3_T * PT;
  }else{
    // Weight is not constant, we use a weighted gram matrix
    QuadratureRule quad_2kpw_T = generate_quadrature_rule(T, 2 * degree() + weight.deg(T));
    auto basis_Pk3_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2kpw_T);
    std::function<double(const Eigen::Vector3d &)> weight_T 
              = [&T, &weight](const Eigen::Vector3d &x)->double {
                  return weight.value(T, x);
                };
    L2P += PT.transpose() * compute_weighted_gram_matrix(weight_T, basis_Pk3_T_quad, basis_Pk3_T_quad, quad_2kpw_T, "sym") * PT;
  }
 
  return L2P;
}
