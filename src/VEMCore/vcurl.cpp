#include <vcurl.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_face.hpp>
#include <GMpoly_edge.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

VCurl::VCurl(const VEMCore & vem_core, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(
	     vem_core.mesh(),
	     0,
	     PolynomialSpaceDimension<Edge>::Poly(vem_core.degree()),
	     PolynomialSpaceDimension<Face>::RolyCompl(vem_core.degree() + 1) + PolynomialSpaceDimension<Face>::Poly(vem_core.degree())-1,
	     PolynomialSpaceDimension<Cell>::RolyCompl(vem_core.degree()) + PolynomialSpaceDimension<Cell>::GolyCompl(vem_core.degree()+1)	     
	     ),
    m_vem_core(vem_core),
    m_use_threads(use_threads),
    m_output(output),
    m_cell_operators(vem_core.mesh().n_cells()),
    m_face_operators(vem_core.mesh().n_faces())
{
  output << "[VCurl] Initializing" << std::endl;
  if (use_threads) {
    m_output << "[VCurl] Parallel execution" << std::endl;
  } else {
    m_output << "[VCurl] Sequential execution" << std::endl;
  }

  // Construct face operators
  std::function<void(size_t, size_t)> construct_all_face_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iF = start; iF < end; iF++) {
          m_face_operators[iF].reset( new LocalOperators(_compute_face_operators(iF)) );
        } // for iF
      };

  m_output << "[VCurl] Constructing face operators" << std::endl;
  parallel_for(mesh().n_faces(), construct_all_face_operators, use_threads);

  // Construct cell operators
  std::function<void(size_t, size_t)> construct_all_cell_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_cell_operators[iT].reset( new LocalOperators(_compute_cell_operators(iT)) );
        } // for iT
      };

  m_output << "[VCurl] Constructing cell operators" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_operators, use_threads);  
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd VCurl::interpolate(const FunctionType & v, const FunctionType & curl_v, const int doe_cell, const int doe_face, const int doe_edge) const
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

  // Interpolate at faces
  std::function<void(size_t, size_t)> interpolate_faces
    = [this, &vh, v, curl_v, &dqr_face](size_t start, size_t end)->void
      {
        for (size_t iF = start; iF < end; iF++) {
          const Face & F = *mesh().face(iF);

          Eigen::Vector3d nF = F.normal();
          auto nF_cross_v_cross_nF = [&nF, v](const Eigen::Vector3d & x)->Eigen::Vector3d {
                                       return nF.cross(v(x).cross(nF));
                                     };

          QuadratureRule quad_dqr_F = generate_quadrature_rule(F, dqr_face);

          size_t offset_F = globalOffset(F);
          auto basis_Rckpo_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).RolyComplkpo, quad_dqr_F);
          vh.segment(offset_F, PolynomialSpaceDimension<Face>::RolyCompl(degree() + 1))
            = l2_projection(nF_cross_v_cross_nF, *faceBases(iF).RolyComplkpo, quad_dqr_F, basis_Rckpo_F_quad);

          if (degree() > 0 ) {
            offset_F += PolynomialSpaceDimension<Face>::RolyCompl(degree() + 1);
            auto basis_Pk0_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).Polyk0, quad_dqr_F);
            auto rot_F_v = [&nF, curl_v](const Eigen::Vector3d & x)->double {
                                       return nF.dot(curl_v(x));
                                     };
            vh.segment(offset_F, faceBases(iF).Polyk0->dimension() )
              = l2_projection(rot_F_v, *faceBases(iF).Polyk0, quad_dqr_F, basis_Pk0_F_quad);
          } // if degree()>0
        } // for iF
      };
  parallel_for(mesh().n_faces(), interpolate_faces, m_use_threads);

  // Interpolate at cells
  std::function<void(size_t, size_t)> interpolate_cells
    = [this, &vh, v, curl_v, &dqr_cell](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          const Cell & T = *mesh().cell(iT);

          QuadratureRule quad_dqr_T = generate_quadrature_rule(T, dqr_cell);
          MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2);

          Eigen::Index offset_T = globalOffset(T);
          if (degree() > 0 ) {
            auto basis_Rck_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).RolyComplk, quad_dqr_T);
            vh.segment(offset_T, PolynomialSpaceDimension<Cell>::RolyCompl(degree()))
              = l2_projection(v, *cellBases(iT).RolyComplk, quad_dqr_T, basis_Rck_T_quad, GramMatrix(T, *cellBases(iT).RolyComplk, int_mono_2kp2));
          }
          
          offset_T += PolynomialSpaceDimension<Cell>::RolyCompl(degree());
          auto basis_Gckpo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).GolyComplkpo, quad_dqr_T);
          vh.segment(offset_T, PolynomialSpaceDimension<Cell>::GolyCompl(degree()+1))
            = l2_projection(curl_v, *cellBases(iT).GolyComplkpo, quad_dqr_T, basis_Gckpo_T_quad, GramMatrix(T, *cellBases(iT).GolyComplkpo, int_mono_2kp2));
        } // for iT
      };
  parallel_for(mesh().n_cells(), interpolate_cells, m_use_threads);
  
  return vh;
}

//------------------------------------------------------------------------------
// Operators
//------------------------------------------------------------------------------

VCurl::LocalOperators VCurl::_compute_face_operators(size_t iF)
{
  const Face & F = *mesh().face(iF);
  
  MonomialFaceIntegralsType int_monoF_2kp3 = IntegrateFaceMonomials(F, 2*degree()+3);

  //------------------------------------------------------------------------------
  // Computation of rot in P^k
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  Eigen::MatrixXd MCF = GramMatrix(F, *faceBases(iF).Polyk, int_monoF_2kp3);
  // The first row should consist in the basis functions integrated vs. 1. The first
  // basis function is constant but may not be 1, so we have to re-scale this row accordingly
  double scaling = 1./faceBases(iF).Polyk->function(0,F.center_mass());
  MCF.row(0) *= scaling;
  
  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BCF
    = Eigen::MatrixXd::Zero(faceBases(iF).Polyk->dimension(), dimensionFace(iF));

  // First row
  for (size_t iE = 0; iE < F.n_edges(); iE++) {
    const Edge & E = *F.edge(iE);
    // Basis made of the constant function 1
    MonomialScalarBasisEdge basis_P0_E(E, 0);
    BCF.block(0, localOffset(F, E), 1, edgeBases(E).Polyk->dimension())
      -= F.edge_orientation(iE) * GramMatrix(E, basis_P0_E, *edgeBases(E).Polyk);
  } // for iE

  // Other rows
  if (degree()>0){
    BCF.bottomRightCorner(faceBases(iF).Polyk0->dimension(),faceBases(iF).Polyk0->dimension())
        = GramMatrix(F, *faceBases(iF).Polyk0, int_monoF_2kp3);
  }
 
  Eigen::MatrixXd proj_curl = MCF.partialPivLu().solve(BCF);
  
  //------------------------------------------------------------------------------
  // Projection of function on P^{k+1}
  //------------------------------------------------------------------------------

  auto basis_Pkp20_F = ShiftedBasis<typename VEMCore::PolyBasisFaceType>(*faceBases(iF).Polykp2, 1);

  // LHS & RHS matrix
  Eigen::MatrixXd MPF
    = Eigen::MatrixXd::Zero(faceBases(iF).Polykpo2->dimension(), faceBases(iF).Polykpo2->dimension());
  Eigen::MatrixXd BPF
    = Eigen::MatrixXd::Zero(faceBases(iF).Polykpo2->dimension(), dimensionFace(iF));

  // Contribution rot_F r
  CurlBasis<decltype(basis_Pkp20_F)> rot_Pkp20_F(basis_Pkp20_F);
  MPF.topRows(basis_Pkp20_F.dimension()) = GramMatrix(F, rot_Pkp20_F, *faceBases(iF).Polykpo2, int_monoF_2kp3);

  if (degree()>0){
    BPF.topRightCorner(faceBases(iF).Polyk0->dimension(), faceBases(iF).Polyk0->dimension())
      = GramMatrix(F, *faceBases(iF).Polyk0, int_monoF_2kp3);
  }
  
  // Edges
  for (size_t iE = 0; iE < F.n_edges(); iE++) {
    const Edge & E = *F.edge(iE);
    QuadratureRule quad_2kp2_E = generate_quadrature_rule(E, 2 * degree() + 2);
    BPF.block(0, localOffset(F, E), basis_Pkp20_F.dimension(), edgeBases(E).Polyk->dimension())
      += F.edge_orientation(iE) * compute_gram_matrix(
						      evaluate_quad<Function>::compute(basis_Pkp20_F, quad_2kp2_E),
						      evaluate_quad<Function>::compute(*edgeBases(E).Polyk, quad_2kp2_E),
						      quad_2kp2_E
						      );    
  } // for iE

  // Contribution x_F p
  MPF.bottomRows(faceBases(iF).RolyComplkpo->dimension())
      = GramMatrix(F, *faceBases(iF).RolyComplkpo, *faceBases(iF).Polykpo2, int_monoF_2kp3);
 
  BPF.block(basis_Pkp20_F.dimension(), localOffset(F), faceBases(iF).RolyComplkpo->dimension(), faceBases(iF).RolyComplkpo->dimension())
    += GramMatrix(F, *faceBases(iF).RolyComplkpo, int_monoF_2kp3);

  Eigen::MatrixXd proj_function = MPF.partialPivLu().solve(BPF);

  // DOFs of curl.n=rot on the face is just its projection
  return LocalOperators(proj_curl, proj_curl, proj_function);
}

//------------------------------------------------------------------------------

VCurl::LocalOperators VCurl::_compute_cell_operators(size_t iT)
{
  const Cell & T = *mesh().cell(iT);

  MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2);

  //------------------------------------------------------------------------------
  // Projection of curl on P^k
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  Eigen::MatrixXd MCT = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(), cellBases(iT).Polyk3->dimension());

  GradientBasis<std::remove_reference<decltype(*cellBases(iT).Polykpo0)>::type>
      grad_Pkpo0_T(*cellBases(iT).Polykpo0);
  MCT.topRows(cellBases(iT).Polykpo0->dimension()) 
      = GramMatrix(T, grad_Pkpo0_T, *cellBases(iT).Polyk3, int_mono_2kp2);
  if (degree()>0){
    MCT.bottomRows(cellBases(iT).GolyComplk->dimension()) 
      = GramMatrix(T, *cellBases(iT).GolyComplk, *cellBases(iT).Polyk3, int_mono_2kp2);
  }
  
  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BCT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(), dimensionCell(iT));

  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);

    MonomialScalarBasisFace basis_Pkpo_F(F, degree()+1);
    DecomposePoly dec(F, basis_Pkpo_F);
    auto basis_Pkpo0_T_nodes = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo0, dec.get_nodes());
    auto Pkpo0T_on_PkpoF = dec.family(basis_Pkpo0_T_nodes);
    Eigen::MatrixXd CF = extendOperator(T, F, faceOperators(F).proj_curl);
    MonomialFaceIntegralsType int_mono_2kp1_F = IntegrateFaceMonomials(F, 2*degree()+1);
    BCT.topRows(cellBases(iT).Polykpo0->dimension()) += T.face_orientation(iF) * GramMatrix(F, Pkpo0T_on_PkpoF, *faceBases(F).Polyk, int_mono_2kp1_F) * CF;
    
  } // for iF

  if (degree()>0){
    BCT.bottomRightCorner(cellBases(iT).GolyComplk->dimension(), cellBases(iT).GolyComplkpo->dimension())
      += GramMatrix(T, *cellBases(iT).GolyComplk, *cellBases(iT).GolyComplkpo, int_mono_2kp2);
  }
  
  Eigen::MatrixXd proj_curl = MCT.partialPivLu().solve(BCT);
  
  //------------------------------------------------------------------------------
  // Projection of function on P^k
  //------------------------------------------------------------------------------
  
  Eigen::MatrixXd MPT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(), cellBases(iT).Polyk3->dimension());  
  
  Eigen::MatrixXd BPT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(), dimensionCell(iT));

  CurlBasis<VEMCore::GolyComplBasisCellType> Rolyk_basis(*cellBases(iT).GolyComplkpo);
  MPT.topRows(cellBases(iT).GolyComplkpo->dimension()) = GramMatrix(T, Rolyk_basis, *cellBases(iT).Polyk3, int_mono_2kp2);
  
  if (degree() > 0) {
    MPT.bottomRows(cellBases(iT).RolyComplk->dimension()) = GramMatrix(T, *cellBases(iT).RolyComplk, *cellBases(iT).Polyk3, int_mono_2kp2);
    BPT.block(cellBases(iT).GolyComplkpo->dimension(), localOffset(T), cellBases(iT).RolyComplk->dimension(), cellBases(iT).RolyComplk->dimension())
      = GramMatrix(T, *cellBases(iT).RolyComplk, int_mono_2kp2);
  } // if degree() > 0

  BPT.topRightCorner(cellBases(iT).GolyComplkpo->dimension(), cellBases(iT).GolyComplkpo->dimension())
    += GramMatrix(T, *cellBases(iT).GolyComplkpo, int_mono_2kp2);
    
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    Eigen::Vector3d nF = F.normal();
    Eigen::MatrixXd PF = extendOperator(T, F, faceOperators(F).proj_function);
  
    MonomialScalarBasisFace basis_Pkpo_F(F, degree()+1);
    DecomposePoly dec(F, TangentFamily<MonomialScalarBasisFace>(basis_Pkpo_F, basis_Pkpo_F.coordinates_system()));
    boost::multi_array<VectorRd, 2> Gkpo_T_cross_nF_nodes 
            = vector_product(evaluate_quad<Function>::compute(*cellBases(iT).GolyComplkpo, dec.get_nodes()), nF);
    auto Gkpo_T_cross_nF_family = dec.family(Gkpo_T_cross_nF_nodes);
    MonomialFaceIntegralsType int_mono_2kp4_F = IntegrateFaceMonomials(F, 2*degree()+4);
    BPT.topRows(cellBases(iT).GolyComplkpo->dimension())
      -= T.face_orientation(iF) * GramMatrix(F, Gkpo_T_cross_nF_family, *faceBases(F).Polykpo2, int_mono_2kp4_F) * PF;

  } // for iF
  
  Eigen::MatrixXd proj_function = MPT.partialPivLu().solve(BPT);
  
  //------------------------------------------------------------------------------
  // DOFs of curl in Vdiv(T)
  //------------------------------------------------------------------------------
  size_t dim_Gckpo_T = PolynomialSpaceDimension<Cell>::GolyCompl(degree()+1);
  size_t n_dofs_face = PolynomialSpaceDimension<Face>::Poly(degree());
  size_t n_dofs_element = T.n_faces() * n_dofs_face 
                      + PolynomialSpaceDimension<Cell>::Poly(degree()) - 1
                      + dim_Gckpo_T;

  Eigen::MatrixXd dofs_curl = Eigen::MatrixXd::Zero(n_dofs_element, dimensionCell(iT));
  for (size_t iF=0; iF<T.n_faces(); iF++){
    const Face & F = *T.face(iF);
    Eigen::MatrixXd CF = extendOperator(T, F, faceOperators(F).dofs_curl);
    dofs_curl.middleRows(iF * n_dofs_face, n_dofs_face) = CF;
  }
  dofs_curl.bottomRightCorner(dim_Gckpo_T, dim_Gckpo_T) = Eigen::MatrixXd::Identity(dim_Gckpo_T, dim_Gckpo_T);
  
  return LocalOperators(proj_curl, dofs_curl, proj_function);
}

//------------------------------------------------------------------------------
//        Compute matrices for local L2 product on VCurl
//------------------------------------------------------------------------------

Eigen::MatrixXd VCurl::computeL2Product(
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
      MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2);
      w_mass_Pk3_T = weight.value(T, T.center_mass()) * GramMatrix(T, *cellBases(iT).Polyk3, int_mono_2kp2);
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

  
  Eigen::MatrixXd L2T = Eigen::MatrixXd::Zero(dimensionCell(iT), dimensionCell(iT));

  // STABILISATION: see implementation notes

  // Edges
  for (size_t iE=0; iE < T.n_edges(); iE++){
    const Edge & E = *T.edge(iE);
    
    // Tangential projections of Pk(T)^3 on E
    const VectorRd tE = E.tangent();
    
    QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2*degree());
    boost::multi_array<double, 2> Pk3T_dot_tE_quad
        = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_E), tE);
    boost::multi_array<double, 2> PkE_quad
        = evaluate_quad<Function>::compute(*edgeBases(E).Polyk, quad_2k_E);

    // Gram matrices on E:
    // Polyk(E)-Polyk(E)
    Eigen::MatrixXd MEE = GramMatrix(E, *edgeBases(E).Polyk);
    // Polyk(E)-Pk3(T).tE
    Eigen::MatrixXd MET = compute_gram_matrix(PkE_quad, Pk3T_dot_tE_quad, quad_2k_E);
    // Pk3(T).tE-Pk3(T).tE
    Eigen::MatrixXd MTT = compute_gram_matrix(Pk3T_dot_tE_quad, quad_2k_E);
    
    // The stabilisation on E is obtained developing \int_E (dofE(v) - (pi^k v).tE).(dofE(w) - (pi^k w).tE)
    // (we don't include scaling coming from weight here)
    L2T.block(localOffset(T, E), localOffset(T, E), edgeBases(E).Polyk->dimension(), edgeBases(E).Polyk->dimension())
        += std::pow(T.diam(), 2) * MEE;
    L2T.middleRows(localOffset(T, E), edgeBases(E).Polyk->dimension())
        -= std::pow(T.diam(), 2) * MET * cellOperators(iT).proj_function;
    L2T.middleCols(localOffset(T, E), edgeBases(E).Polyk->dimension())
        -= (std::pow(T.diam(), 2) * MET * cellOperators(iT).proj_function).transpose();
    L2T += std::pow(T.diam(), 2) *  cellOperators(iT).proj_function.transpose() * MTT *  cellOperators(iT).proj_function;
    
  }
  
  // Faces
  for (size_t iF=0; iF < T.n_faces(); iF++){
    const Face & F = *T.face(iF);
    
    // Tangential projections of Pk(T)^3 on F
    const VectorRd nF = F.normal();
    
    QuadratureRule quad_2kp2_F = generate_quadrature_rule(F, 2*degree()+2);
    boost::multi_array<VectorRd, 2> Pk3T_tF_quad
        = vector_product(
            vector_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2kp2_F), nF), 
            -nF);
    boost::multi_array<VectorRd, 2> Pkpo2F_quad
        = evaluate_quad<Function>::compute(*faceBases(F).Polykpo2, quad_2kp2_F);

    // Gram matrices on F:
    // Polykpo2(F)-Polykpo2(F)
    Eigen::MatrixXd MFF = GramMatrix(F, *faceBases(F).Polykpo2);
    // Polykpo2(F)-Pk3(T)_tF
    Eigen::MatrixXd MFT = compute_gram_matrix(Pkpo2F_quad, Pk3T_tF_quad, quad_2kp2_F);
    // Pk3(T)_tF-Pk3(T)_tF
    Eigen::MatrixXd MTT = compute_gram_matrix(Pk3T_tF_quad, quad_2kp2_F);
    
    // The stabilisation on F is obtained developing \int_F (dofF(v) - (pi^k v)_tF).(dofF(w) - (pi^k w)_tF)
    // (we don't include scaling coming from weight here)
    Eigen::MatrixXd proj_kpo_F = extendOperator(T, F, faceOperators(F).proj_function);
    Eigen::MatrixXd mixed_product = T.diam() * (proj_kpo_F.transpose() * MFT * cellOperators(iT).proj_function);
    L2T += T.diam() * proj_kpo_F.transpose() * MFF * proj_kpo_F;
    L2T -= mixed_product;
    L2T -= mixed_product.transpose();
    L2T += T.diam() *  cellOperators(iT).proj_function.transpose() * MTT *  cellOperators(iT).proj_function;
  }

  L2T *= penalty_factor;
  
  // CONSISTENT TERM
  L2T += cellOperators(iT).proj_function.transpose() * w_mass_Pk3_T * cellOperators(iT).proj_function;

  return L2T;  
}



