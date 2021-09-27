#include <xdiv.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_face.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

XDiv::XDiv(const DDRCore & ddr_core, bool use_threads, std::ostream & output)
  : DDRSpace(
	     ddr_core.mesh(), 0, 0,
	     PolynomialSpaceDimension<Face>::Poly(ddr_core.degree()),
	     PolynomialSpaceDimension<Cell>::Goly(ddr_core.degree() - 1) + PolynomialSpaceDimension<Cell>::GolyCompl(ddr_core.degree())
	      ),
    m_ddr_core(ddr_core),
    m_use_threads(use_threads),
    m_output(output),
    m_cell_operators(ddr_core.mesh().n_cells())
{
  m_output << "[XDiv] Initializing" << std::endl;
  if (use_threads) {
    m_output << "[XDiv] Parallel execution" << std::endl;
  } else {
    m_output << "[XDiv] Sequential execution" << std::endl;
  }

  // Construct cell divergences and potentials
  std::function<void(size_t, size_t)> construct_all_cell_divergences_potentials
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_cell_operators[iT].reset( new LocalOperators(_compute_cell_divergence_potential(iT)) );
        } // for iT
      };

  m_output << "[XDiv] Constructing cell divergences and potentials" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_divergences_potentials, use_threads);
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd XDiv::interpolate(const FunctionType & v) const
{
  Eigen::VectorXd vh = Eigen::VectorXd::Zero(dimension());

  // Interpolate at faces
  std::function<void(size_t, size_t)> interpolate_faces
    = [this, &vh, v](size_t start, size_t end)->void
      {
	for (size_t iF = start; iF < end; iF++) {
	  const Face & F = *mesh().face(iF);

	  Eigen::Vector3d nF = F.normal();
	  auto v_dot_nF = [&nF, v](const Eigen::Vector3d & x)->double {
			    return v(x).dot(nF);
			  };
	  
	  QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree());
	  auto basis_Pk_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).Polyk, quad_2k_F);
	  vh.segment(globalOffset(F), PolynomialSpaceDimension<Face>::Poly(degree()))
	    = l2_projection(v_dot_nF, *faceBases(iF).Polyk, quad_2k_F, basis_Pk_F_quad);
	} // for iF
      };
  parallel_for(mesh().n_faces(), interpolate_faces, m_use_threads);
  
  // Interpolate at cells
  if (degree() > 0) {  
    std::function<void(size_t, size_t)> interpolate_cells
      = [this, &vh, v](size_t start, size_t end)->void
	{
	  for (size_t iT = start; iT < end; iT++) {		 
	    const Cell & T = *mesh().cell(iT);

	    size_t offset_T = globalOffset(T);

	    QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * degree());	    
	    auto basis_Gkmo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Golykmo, quad_2k_T);
	    vh.segment(offset_T, PolynomialSpaceDimension<Cell>::Goly(degree() - 1))
	      = l2_projection(v, *cellBases(iT).Golykmo, quad_2k_T, basis_Gkmo_T_quad);

	    offset_T += PolynomialSpaceDimension<Cell>::Goly(degree() - 1);	    
	    auto basis_GOk_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).GolyComplk, quad_2k_T);
	    vh.segment(offset_T, PolynomialSpaceDimension<Cell>::GolyCompl(degree()))
	      = l2_projection(v, *cellBases(iT).GolyComplk, quad_2k_T, basis_GOk_T_quad);
	  } // for iT
	};
    parallel_for(mesh().n_cells(), interpolate_cells, m_use_threads);
  } // if degree() > 0
	
  return vh;
}

//------------------------------------------------------------------------------
// Divergence and potential reconstruction
//------------------------------------------------------------------------------

XDiv::LocalOperators XDiv::_compute_cell_divergence_potential(size_t iT)
{
  const Cell & T = *mesh().cell(iT);

  //------------------------------------------------------------------------------
  // Divergence
  //------------------ ------------------------------------------------------------

  //-----------------------------------------------------------------------------------------
  // Left-hand side matrix

  MonomialIntegralsType int_mono_2kpo = IntegrateCellMonomials(T, 2*degree()+1);
  Eigen::MatrixXd MDT = GramMatrix(T, *cellBases(iT).Polyk, int_mono_2kpo);

  //------------------------------------------------------------------------------
  // Right-hand side matrix
  
  Eigen::MatrixXd BDT
    = Eigen::MatrixXd::Zero(PolynomialSpaceDimension<Cell>::Poly(degree()), dimensionCell(iT));

  // Boundary contribution
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);

    QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree());
    BDT.block(0, iF * PolynomialSpaceDimension<Face>::Poly(degree()), PolynomialSpaceDimension<Cell>::Poly(degree()), PolynomialSpaceDimension<Face>::Poly(degree()))
      += T.face_orientation(iF) * compute_gram_matrix(
						      evaluate_quad<Function>::compute(*cellBases(iT).Polyk, quad_2k_F),
						      evaluate_quad<Function>::compute(*faceBases(F).Polyk, quad_2k_F),
						      quad_2k_F
						      );
  } // for iF

  // Element contribution
  if (degree() > 0) {
    GradientBasis<RestrictedBasis<DDRCore::PolyBasisCellType>> grad_Polyk(*cellBases(iT).Polyk);
    BDT.block(0, localOffset(T), PolynomialSpaceDimension<Cell>::Poly(degree()), PolynomialSpaceDimension<Cell>::Goly(degree() - 1))
      			           -= GramMatrix(T, grad_Polyk, *cellBases(iT).Golykmo);
  } // if degree() > 0

  Eigen::MatrixXd DT = MDT.ldlt().solve(BDT);
  
  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------

  auto Poly0kpo = ShiftedBasis<DDRCore::PolyBasisCellType>(*cellBases(iT).Polykpo, 1);

  Eigen::MatrixXd MPT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(),
			    cellBases(iT).Polyk3->dimension());

  Eigen::MatrixXd BPT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(), dimensionCell(iT));

////  auto basis_Pk3_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_T);
  
  GradientBasis<ShiftedBasis<DDRCore::PolyBasisCellType>> grad_Poly0kpo(Poly0kpo);
  MPT.topRows(Poly0kpo.dimension()) = GramMatrix(T, grad_Poly0kpo, *cellBases(iT).Polyk3, int_mono_2kpo);

  if (degree() > 0) {
////    auto basis_GOk_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).GolyComplk, quad_2k_T);
////    MPT.bottomRows(cellBases(iT).GolyComplk->dimension())
////      = compute_gram_matrix(basis_GOk_T_quad, basis_Pk3_T_quad, quad_2k_T);
////    BPT.bottomRightCorner(cellBases(iT).GolyComplk->dimension(), cellBases(iT).GolyComplk->dimension())
////      = compute_gram_matrix(basis_GOk_T_quad, quad_2k_T);		   
    MPT.bottomRows(cellBases(iT).GolyComplk->dimension()) = GramMatrix(T, *cellBases(iT).GolyComplk, *cellBases(iT).Polyk3, int_mono_2kpo);
    BPT.bottomRightCorner(cellBases(iT).GolyComplk->dimension(), cellBases(iT).GolyComplk->dimension())
          = GramMatrix(T, *cellBases(iT).GolyComplk, int_mono_2kpo);
  } // if degree() > 0

  BPT.topRows(Poly0kpo.dimension()) -= GramMatrix(T, Poly0kpo, *cellBases(iT).Polyk, int_mono_2kpo) * DT;

  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    QuadratureRule quad_2kpo_F = generate_quadrature_rule(F, 2 * degree() + 1);
    BPT.block(0, localOffset(T, F), Poly0kpo.dimension(), dimensionFace(F))
      += T.face_orientation(iF) * compute_gram_matrix(
						      evaluate_quad<Function>::compute(Poly0kpo, quad_2kpo_F),
						      evaluate_quad<Function>::compute(*faceBases(F.global_index()).Polyk, quad_2kpo_F),
						      quad_2kpo_F
						      );
  } // for iF

  Eigen::MatrixXd PT = MPT.partialPivLu().solve(BPT);

  // // Correction to enforce that the L2-orthogonal projection of PT
  // // on Gk-1(T) is equal to the cell unknown
  // if (degree() > 0) {
  //   auto basis_Gkmo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Golykmo, quad_2k_T);

  //   Eigen::MatrixXd gram_Gkmo_T = compute_gram_matrix(basis_Gkmo_T_quad, quad_2k_T);
  //   Eigen::MatrixXd gram_Gkmo_T_Pk3_T = compute_gram_matrix(basis_Gkmo_T_quad, basis_Pk3_T_quad, quad_2k_T);
  //   Eigen::MatrixXd gram_Pk3_T = compute_gram_matrix(basis_Pk3_T_quad, quad_2k_T);
  //   Eigen::MatrixXd proj_Gkmo_T_Pk3_T = gram_Pk3_T.ldlt().solve(gram_Gkmo_T_Pk3_T.transpose());

  //   // Remove the L2-orthogonal projection of PT on Gk-1(T) and replace
  //   // it with the cell unknown
  //   PT -= proj_Gkmo_T_Pk3_T * gram_Gkmo_T.ldlt().solve(gram_Gkmo_T_Pk3_T * PT);
  //   PT.middleCols(localOffset(T), PolynomialSpaceDimension<Cell>::Goly(degree() - 1)) += proj_Gkmo_T_Pk3_T;
  // }
  
  return LocalOperators(DT, BDT, PT);
}

//------------------------------------------------------------------------------
//        Functions to compute matrices for local L2 products on Xdiv
//------------------------------------------------------------------------------

Eigen::MatrixXd XDiv::computeL2Product(
                                        const size_t iT,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pk3_T,
                                        const IntegralWeight & weight
                                        ) const
{
  const Cell & T = *mesh().cell(iT); 
  
  // leftOp and rightOp come from the potentials
  std::vector<Eigen::MatrixXd> leftOp(T.n_faces()+1);
  for (size_t iF = 0; iF < T.n_faces(); iF++){
    const Face & F = *T.face(iF);
    leftOp[iF] = extendOperator(T, F, Eigen::MatrixXd::Identity(dimensionFace(F),dimensionFace(F)));
  }
  leftOp[T.n_faces()] = m_cell_operators[iT]->potential;
  std::vector<Eigen::MatrixXd> rightOp = leftOp;

  return computeL2Product_with_Ops(iT, leftOp, rightOp, penalty_factor, mass_Pk3_T, weight);

}

Eigen::MatrixXd XDiv::computeL2ProductCurl(
                                        const size_t iT,
                                        const XCurl & x_curl,
                                        const std::string side,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pk3_T,
                                        const IntegralWeight & weight
                                        ) const
{
  const Cell & T = *mesh().cell(iT);

  // list of full curls
  std::vector<Eigen::MatrixXd> curlOp(T.n_faces()+1);
  for (size_t iF = 0; iF < T.n_faces(); iF++){
    const Face & F = *T.face(iF);
    curlOp[iF] = x_curl.extendOperator(T, F, x_curl.faceOperators(F).curl);
  }
  curlOp[T.n_faces()] = x_curl.cellOperators(iT).curl;
  
  // If we apply the curl on one side only we'll need the potentials
  if (side != "both"){
    // list of potentials
    std::vector<Eigen::MatrixXd> potentialOp(T.n_faces()+1);
    for (size_t iF = 0; iF < T.n_faces(); iF++){
      const Face & F = *T.face(iF);
      potentialOp[iF] = extendOperator(T, F, Eigen::MatrixXd::Identity(dimensionFace(F),dimensionFace(F)));
    }
    potentialOp[T.n_faces()] = m_cell_operators[iT]->potential;

  
    // Depending on side of curl
    if (side == "left"){
      return computeL2Product_with_Ops(iT, curlOp, potentialOp, penalty_factor, mass_Pk3_T, weight);
    }else{
      return computeL2Product_with_Ops(iT, potentialOp, curlOp, penalty_factor, mass_Pk3_T, weight);
    }
    
  }

  // Default: gradient on both sides
  return computeL2Product_with_Ops(iT, curlOp, curlOp, penalty_factor, mass_Pk3_T, weight);

}


Eigen::MatrixXd XDiv::computeL2Product_with_Ops(
                                        const size_t iT,
                                        std::vector<Eigen::MatrixXd> & leftOp,
                                        std::vector<Eigen::MatrixXd> & rightOp,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & w_mass_Pk3_T,
                                        const IntegralWeight & weight
                                        ) const
{
  const Cell & T = *mesh().cell(iT); 

  // leftOp and rightOp must list the operators acting on the DOFs, and which we want to
  // use for the L2 product. Specifically, each one lists operators (matrices) returning
  // values in faces space P^k(F) and element space P^k(T)^3.
  // For the standard Xdiv L2 product, these will respectively be identity (for each face), and PT. 
  // To compute the Xdiv L2 product applied (left or right) to the discrete curl,
  // leftOp or rightOp must list the face and element (full) curl operators.
  // All these operators must have the same domain, so possibly being extended appropriately
  // using extendOperator from ddrspace.

  Eigen::MatrixXd L2P = Eigen::MatrixXd::Zero(leftOp[0].cols(), rightOp[0].cols());
  
  size_t offset_T = T.n_faces();

  // Face penalty terms
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);

    // Compute gram matrices
////    QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree());
////    auto basis_Pk3_T_dot_nF_quad
////      = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_F), F.normal());
////    auto basis_Pk_F_quad
////      = evaluate_quad<Function>::compute(*faceBases(F.global_index()).Polyk, quad_2k_F);
////    Eigen::MatrixXd mass_PkF_PkF = compute_gram_matrix(basis_Pk_F_quad, quad_2k_F);
////    Eigen::MatrixXd gram_PkF_Pk3T_dot_nF = compute_gram_matrix(basis_Pk_F_quad, basis_Pk3_T_dot_nF_quad, quad_2k_F, "nonsym");
    MonomialFaceIntegralsType int_monoF_2k = IntegrateFaceMonomials(F, 2*degree());
    DecomposePoly dec(F, MonomialScalarBasisFace(F, degree()));
    auto Pk3_T_dot_nF_nodes = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, dec.get_nodes()), F.normal());
    auto Pk3_T_dot_nF_family_PkF = dec.family(Pk3_T_dot_nF_nodes);
    Eigen::MatrixXd mass_PkF_PkF = GramMatrix(F, *faceBases(F).Polyk, int_monoF_2k);
    Eigen::MatrixXd gram_PkF_Pk3T_dot_nF = GramMatrix(F, *faceBases(F).Polyk, Pk3_T_dot_nF_family_PkF, int_monoF_2k);

    // Weight including scaling hF
    QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree());
    double max_weight_quad_F = weight.value(T, quad_2k_F[0].vector());
    // If the weight is not constant, we want to take the largest along the edge
    if (weight.deg(T)>0){
      for (size_t iqn = 1; iqn < quad_2k_F.size(); iqn++) {
        max_weight_quad_F = std::max(max_weight_quad_F, weight.value(T, quad_2k_F[iqn].vector()));
      } // for
    }
    double w_hF = max_weight_quad_F * F.diam();

    // The penalty term int_T (leftOp.nF - (leftOp)_F) * (rightOp.nF - (rightOp)_F) is computed by developping
    // Contribution of face F
    L2P += w_hF * ( leftOp[offset_T].transpose() * gram_PkF_Pk3T_dot_nF.transpose() * mass_PkF_PkF.ldlt().solve(gram_PkF_Pk3T_dot_nF) * rightOp[offset_T]
                - leftOp[offset_T].transpose() * gram_PkF_Pk3T_dot_nF.transpose() * rightOp[iF] 
                - leftOp[iF].transpose() * gram_PkF_Pk3T_dot_nF * rightOp[offset_T]
                + leftOp[iF].transpose() * mass_PkF_PkF * rightOp[iF]                  
                  );

  } // for iF

  L2P *= penalty_factor;
  
  // Consistent term, two calculations depending if weight is constant or not.
////  if (weight.deg(T)==0){
////    // weight is constant so the calculation is based on the standard mass matrix multiplied by the weight. 
////    Eigen::MatrixXd M_Pk3_T = mass_Pk3_T;
////    if (M_Pk3_T.rows()==1){
////      // Mass matrix not passed as parameter, we compute it
////      QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * degree());
////      M_Pk3_T.resize(cellBases(iT).Polyk3->dimension(), cellBases(iT).Polyk3->dimension());
////      M_Pk3_T = compute_gram_matrix(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_T), quad_2k_T);
////    }
////    L2P += weight.value(T, T.center_mass()) * leftOp[offset_T].transpose() * M_Pk3_T * rightOp[offset_T];
////  }else{
////    // Weight is not constant, we use a weighted gram matrix
////    QuadratureRule quad_2kpw_T = generate_quadrature_rule(T, 2 * degree() + weight.deg(T));
////    auto basis_Pk3_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2kpw_T);
////    std::function<double(const Eigen::Vector3d &)> weight_T 
////              = [&T, &weight](const Eigen::Vector3d &x)->double {
////                  return weight.value(T, x);
////                };
////    L2P += leftOp[offset_T].transpose() * compute_weighted_gram_matrix(weight_T, basis_Pk3_T_quad, basis_Pk3_T_quad, quad_2kpw_T, "sym") * rightOp[offset_T];

////  }
   L2P += leftOp[offset_T].transpose() * w_mass_Pk3_T * rightOp[offset_T];

 
  return L2P;
}


