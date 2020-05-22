#include <xdiv.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

XDiv::XDiv(const DDRCore & ddr_core, bool use_threads, std::ostream & output)
  : DDRSpace(
	     ddr_core.mesh(), 0, 0,
	     PolynomialSpaceDimension<Face>::Poly(ddr_core.degree()),
	     PolynomialSpaceDimension<Cell>::Goly(ddr_core.degree() - 1) + PolynomialSpaceDimension<Cell>::GolyOrth(ddr_core.degree())
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
	    auto basis_GOk_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).GolyOrthk, quad_2k_T);
	    vh.segment(offset_T, PolynomialSpaceDimension<Cell>::GolyOrth(degree()))
	      = l2_projection(v, *cellBases(iT).GolyOrthk, quad_2k_T, basis_GOk_T_quad);
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

  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * degree());						       
  auto basis_Pk_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk, quad_2k_T);
  Eigen::MatrixXd MDT = compute_gram_matrix(basis_Pk_T_quad, basis_Pk_T_quad, quad_2k_T, "sym");

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
    QuadratureRule quad_2kmo_T = generate_quadrature_rule(T, 2 * (degree() - 1));
    BDT.block(0, localOffset(T), PolynomialSpaceDimension<Cell>::Poly(degree()), PolynomialSpaceDimension<Cell>::Goly(degree() - 1))
      -= compute_gram_matrix(
			     evaluate_quad<Gradient>::compute(*cellBases(iT).Polyk, quad_2kmo_T),
			     evaluate_quad<Function>::compute(*cellBases(iT).Golykmo, quad_2kmo_T),
			     quad_2kmo_T
			     );
  } // if degree() > 0

  Eigen::MatrixXd DT = MDT.ldlt().solve(BDT);
  
  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------

  auto Poly0kpo = ShiftedBasis<DDRCore::PolyCellBasisType>(*cellBases(iT).Polykpo, 1);

  Eigen::MatrixXd MPT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(),
			    cellBases(iT).Polyk3->dimension());

  Eigen::MatrixXd BPT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(), dimensionCell(iT));

  auto basis_Pk3_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_T);
    
  MPT.topRows(Poly0kpo.dimension())
    =  compute_gram_matrix(
			   evaluate_quad<Gradient>::compute(Poly0kpo, quad_2k_T),
			   basis_Pk3_T_quad,
			   quad_2k_T
			   );

  if (degree() > 0) {
    auto basis_GOk_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).GolyOrthk, quad_2k_T);
    MPT.bottomRows(cellBases(iT).GolyOrthk->dimension())
      = compute_gram_matrix(basis_GOk_T_quad, basis_Pk3_T_quad, quad_2k_T);
    BPT.bottomRightCorner(cellBases(iT).GolyOrthk->dimension(), cellBases(iT).GolyOrthk->dimension())
      = compute_gram_matrix(basis_GOk_T_quad, quad_2k_T);		   
  } // if degree() > 0

  auto quad_2kpo_T = generate_quadrature_rule(T, 2 * degree() + 1);
  BPT.topRows(Poly0kpo.dimension())
    -= compute_gram_matrix(
			   evaluate_quad<Function>::compute(Poly0kpo, quad_2kpo_T),
			   evaluate_quad<Function>::compute(*cellBases(iT).Polyk, quad_2kpo_T),
			   quad_2kpo_T
			   ) * DT;

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
  // // on Gk-1(T) is equal to the cell uknown
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
