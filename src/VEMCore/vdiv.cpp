#include <vdiv.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_face.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

VDiv::VDiv(const VEMCore & vem_core, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(
	     vem_core.mesh(),
	     0,
	     0,
	     PolynomialSpaceDimension<Face>::Poly(vem_core.degree()),
	     PolynomialSpaceDimension<Cell>::Poly(vem_core.degree()) -1
	      + PolynomialSpaceDimension<Cell>::GolyCompl(vem_core.degree()+1)
	     ),
    m_vem_core(vem_core),
    m_use_threads(use_threads),
    m_output(output),
    m_cell_operators(vem_core.mesh().n_cells()),
    m_face_operators(vem_core.mesh().n_faces())
{
  output << "[VDiv] Initializing" << std::endl;
  if (use_threads) {
    m_output << "[VDiv] Parallel execution" << std::endl;
  } else {
    m_output << "[VDiv] Sequential execution" << std::endl;
  }

  // Construct cell operators
  std::function<void(size_t, size_t)> construct_all_cell_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_cell_operators[iT].reset( new LocalOperators(_compute_cell_operators(iT)) );
        } // for iT
      };

  m_output << "[VDiv] Constructing cell operators" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_operators, use_threads);  
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd VDiv::interpolate(const FunctionType & v, const DivergenceType & div_v, const int doe_cell, const int doe_face) const
{
  Eigen::VectorXd vh = Eigen::VectorXd::Zero(dimension());
  
  // Degrees of quadrature rules
  size_t dqr_cell = (doe_cell >= 0 ? doe_cell : 2 * degree() + 3);
  size_t dqr_face = (doe_face >= 0 ? doe_face : 2 * degree() + 3);
  
  // Interpolate at faces
  std::function<void(size_t, size_t)> interpolate_faces
    = [this, &vh, v, &dqr_face](size_t start, size_t end)->void
      {
        for (size_t iF = start; iF < end; iF++) {
          const Face & F = *mesh().face(iF);

          Eigen::Vector3d nF = F.normal();
          auto nF_dot_v = [&nF, v](const Eigen::Vector3d & x)->double {
                                     return nF.dot(v(x));
                                     };

          QuadratureRule quad_dqr_F = generate_quadrature_rule(F, dqr_face);

          auto basis_Pk_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).Polyk, quad_dqr_F);
          vh.segment(globalOffset(F), PolynomialSpaceDimension<Face>::Poly(degree()))
            = l2_projection(nF_dot_v, *faceBases(iF).Polyk, quad_dqr_F, basis_Pk_F_quad);
        } // for iF
      };
  parallel_for(mesh().n_faces(), interpolate_faces, m_use_threads);

  // Interpolate at cells
  std::function<void(size_t, size_t)> interpolate_cells
    = [this, &vh, v, div_v, &dqr_cell](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          const Cell & T = *mesh().cell(iT);

          QuadratureRule quad_dqr_T = generate_quadrature_rule(T, dqr_cell);
          MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2);

          Eigen::Index offset_T = globalOffset(T);
          if (degree() > 0 ) {
            auto basis_Pk0_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk0, quad_dqr_T);
            vh.segment(offset_T, cellBases(iT).Polyk0->dimension())
              = l2_projection(div_v, *cellBases(iT).Polyk0, quad_dqr_T, basis_Pk0_quad, GramMatrix(T, *cellBases(iT).Polyk0, int_mono_2kp2));
          }
          
          offset_T += PolynomialSpaceDimension<Cell>::Poly(degree())-1;
          auto basis_Gckpo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).GolyComplkpo, quad_dqr_T);
          vh.segment(offset_T, PolynomialSpaceDimension<Cell>::GolyCompl(degree()+1))
            = l2_projection(v, *cellBases(iT).GolyComplkpo, quad_dqr_T, basis_Gckpo_T_quad, GramMatrix(T, *cellBases(iT).GolyComplkpo, int_mono_2kp2));
        } // for iT
      };
  parallel_for(mesh().n_cells(), interpolate_cells, m_use_threads);
  
  return vh;
}

//------------------------------------------------------------------------------
// Operators
//------------------------------------------------------------------------------

VDiv::LocalOperators VDiv::_compute_cell_operators(size_t iT)
{
  const Cell & T = *mesh().cell(iT);
  
  MonomialCellIntegralsType int_monoT_2kp3 = IntegrateCellMonomials(T, 2*degree()+3);

  //------------------------------------------------------------------------------
  // Computation of div in P^k
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  Eigen::MatrixXd MDT = GramMatrix(T, *cellBases(iT).Polyk, int_monoT_2kp3);
  // The first row should consist in the basis functions integrated vs. 1. The first
  // basis function is constant but may not be 1, so we have to re-scale this row accordingly
  double scaling = 1./cellBases(iT).Polyk->function(0,T.center_mass());
  MDT.row(0) *= scaling;
  
  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BDT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk->dimension(), dimensionCell(iT));

  // First row
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    // Basis made of the constant function 1
    MonomialScalarBasisFace basis_P0_F(F, 0);
    BDT.block(0, localOffset(T, F), 1, faceBases(F).Polyk->dimension())
      = T.face_orientation(iF) * GramMatrix(F, basis_P0_F, *faceBases(F).Polyk);
  } // for iF

  // Other rows
  if (degree()>0){
    BDT.block(1, localOffset(T), cellBases(iT).Polyk0->dimension(), cellBases(iT).Polyk0->dimension())
      = GramMatrix(T, *cellBases(iT).Polyk0, int_monoT_2kp3);
  }
 
  Eigen::MatrixXd proj_div = MDT.partialPivLu().solve(BDT);
  
  //------------------------------------------------------------------------------
  // Projection of function on P^k
  //------------------------------------------------------------------------------
  
  //------------------------------------------------------------------------------
  // Left-hand side matrix
  Eigen::MatrixXd MPT = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(), cellBases(iT).Polyk3->dimension());

  GradientBasis<std::remove_reference<decltype(*cellBases(iT).Polykpo0)>::type>
          grad_Pkpo0_T(*cellBases(iT).Polykpo0);
  MPT.topRows(cellBases(iT).Polykpo0->dimension()) 
      = GramMatrix(T, grad_Pkpo0_T, *cellBases(iT).Polyk3, int_monoT_2kp3);
  if (degree()>0){
    MPT.bottomRows(cellBases(iT).GolyComplk->dimension()) 
      = GramMatrix(T, *cellBases(iT).GolyComplk, *cellBases(iT).Polyk3, int_monoT_2kp3);
  }

  //------------------------------------------------------------------------------
  // Right-hand side matrix  
  Eigen::MatrixXd BPT = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension(), dimensionCell(iT));

  if (degree() > 0) {
    BPT.block(0, localOffset(T), cellBases(iT).Polykpo0->dimension(), cellBases(iT).Polyk0->dimension())
      -= GramMatrix(T, *cellBases(iT).Polykpo0, *cellBases(iT).Polyk0, int_monoT_2kp3);

    BPT.bottomRightCorner(cellBases(iT).GolyComplk->dimension(), cellBases(iT).GolyComplkpo->dimension())
      += GramMatrix(T, *cellBases(iT).GolyComplk, *cellBases(iT).GolyComplkpo, int_monoT_2kp3);
  } // if degree() > 0

  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
  
    MonomialScalarBasisFace basis_Pkpo_F(F, degree()+1);
    DecomposePoly dec(F, basis_Pkpo_F);
    boost::multi_array<double, 2> Pkpo0_T_on_F_nodes 
            = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo0, dec.get_nodes());
    auto Pkpo0_T_on_F_family = dec.family(Pkpo0_T_on_F_nodes);
    MonomialFaceIntegralsType int_mono_2kp1_F = IntegrateFaceMonomials(F, 2*degree()+1);
    BPT.block(0, iF * PolynomialSpaceDimension<Face>::Poly(degree()), cellBases(iT).Polykpo0->dimension(), PolynomialSpaceDimension<Face>::Poly(degree()))
      += T.face_orientation(iF) * GramMatrix(F, Pkpo0_T_on_F_family, *faceBases(F).Polyk, int_mono_2kp1_F);
  } // for iF
  
   Eigen::MatrixXd proj_function = MPT.partialPivLu().solve(BPT);

  return LocalOperators(proj_div, proj_function);
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//        Compute matrices for local L2 product on Vdiv
//------------------------------------------------------------------------------

Eigen::MatrixXd VDiv::computeL2Product(
                                      const size_t iT,
                                      const double & penalty_factor,
                                      const Eigen::MatrixXd & mass_Pk3_T,
                                      const IntegralWeight & weight
                                      ) const
{
  const Cell & T = *mesh().cell(iT); 
  
  // create the weighted mass matrix, with simple product if weight is constant
  Eigen::MatrixXd w_mass_Pk3_T;
  MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2);
  if (weight.deg(T)==0){
    // constant weight
    if (mass_Pk3_T.rows()==1){
      // We have to compute the mass matrix
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

  // Faces
  for (size_t iF=0; iF < T.n_faces(); iF++){
    const Face & F = *T.face(iF);
    
    // Tangential projections of Pk(T)^3 on F
    const VectorRd nF = F.normal();
    
    QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2*degree());
    boost::multi_array<double, 2> Pk3T_nF_quad
        = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_F), nF);
    boost::multi_array<double, 2> PkF_quad
        = evaluate_quad<Function>::compute(*faceBases(F).Polyk, quad_2k_F);

    // Gram matrices on F:
    // Polykp(F)-Polykp(F)
    Eigen::MatrixXd MFF = GramMatrix(F, *faceBases(F).Polyk);
    // Polyk(F)-Pk3(T).nF
    Eigen::MatrixXd MFT = compute_gram_matrix(PkF_quad, Pk3T_nF_quad, quad_2k_F);
    // Pk3(T).nF-Pk3(T).nF
    Eigen::MatrixXd MTT = compute_gram_matrix(Pk3T_nF_quad, quad_2k_F);
    
    // The stabilisation on F is obtained developing \int_F (dofF(v) - (pi^k v).nF).(dofF(w) - (pi^k w).nF)
    // (we don't include scaling coming from weight here)
    Eigen::MatrixXd proj_k_F = 
        extendOperator(T, F, Eigen::MatrixXd::Identity(faceBases(iF).Polyk->dimension(),faceBases(iF).Polyk->dimension()));
    Eigen::MatrixXd mixed_product = T.diam() * proj_k_F.transpose() * MFT * cellOperators(iT).proj_function;
    L2T += T.diam() * proj_k_F.transpose() * MFF * proj_k_F;
    L2T -= mixed_product;
    L2T -= mixed_product.transpose();
    L2T += T.diam() *  cellOperators(iT).proj_function.transpose() * MTT *  cellOperators(iT).proj_function;
  }

  // Element term, also developing \int_T (dofT(v) - pi^{c,k+1}_G(pi^k v)).(dofT(w) - pi^{c,k+1}_G(pi^k w))
  Eigen::MatrixXd MPP = GramMatrix(T, *cellBases(iT).Polyk3, int_mono_2kp2);
  Eigen::MatrixXd MPG = GramMatrix(T, *cellBases(iT).Polyk3, *cellBases(iT).GolyComplkpo, int_mono_2kp2);
  Eigen::MatrixXd MGG = GramMatrix(T, *cellBases(iT).GolyComplkpo, int_mono_2kp2);
  Eigen::MatrixXd proj_Pik_times_MPG = cellOperators(iT).proj_function.transpose() * MPG;
  L2T.bottomRightCorner(cellBases(iT).GolyComplkpo->dimension(),cellBases(iT).GolyComplkpo->dimension())
    += MGG;
  L2T.rightCols(cellBases(iT).GolyComplkpo->dimension())
    -= proj_Pik_times_MPG;
  L2T.bottomRows(cellBases(iT).GolyComplkpo->dimension())
    -= proj_Pik_times_MPG.transpose();
  // pi^{c,k+1}(pi^k) is obtained as MGG^{-1} MGP
  L2T += cellOperators(iT).proj_function.transpose() * MPG * MGG.ldlt().solve(MPG.transpose()) *  cellOperators(iT).proj_function;
  
  L2T *= penalty_factor;
  
  // CONSISTENT TERM
  L2T += cellOperators(iT).proj_function.transpose() * w_mass_Pk3_T * cellOperators(iT).proj_function;

  return L2T;  
}

