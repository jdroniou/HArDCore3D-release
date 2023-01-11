#include <cassert>

#include <hhospace.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_face.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------

HHOSpace::HHOSpace(const Mesh & mesh, size_t K, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(mesh,
	     0,
	     0,
	     PolynomialSpaceDimension<Face>::Poly(K),
	     PolynomialSpaceDimension<Cell>::Poly(K)
	     ),
	  m_mesh(mesh),
    m_K(K),
    m_use_threads(use_threads),
    m_output(output),
    m_cell_bases(mesh.n_cells()),
    m_face_bases(mesh.n_faces()),
    m_operators(mesh.n_cells())
{
  m_output << "[HHOSpace] Initializing" << std::endl;
  
  // Construct element bases
  std::function<void(size_t, size_t)> construct_all_cell_bases
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iT = start; iT < end; iT++) {
	        this->m_cell_bases[iT].reset( new CellBases(this->_construct_cell_bases(iT)) );
	      } // for iT
      };

  m_output << "[HHOSpace] Constructing element bases" << std::endl;
  parallel_for(mesh.n_cells(), construct_all_cell_bases, use_threads);
  
  // Construct face bases
  std::function<void(size_t, size_t)> construct_all_face_bases
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iF = start; iF < end; iF++) {
	        this->m_face_bases[iF].reset( new FaceBases(_construct_face_bases(iF)) );
	      } // for iF
      };
  
  m_output << "[HHOSpace] Constructing face bases" << std::endl;
  parallel_for(mesh.n_faces(), construct_all_face_bases, use_threads);

  // Construct gradients, potentials and stabilisation
  std::function<void(size_t, size_t)> construct_all_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_operators[iT].reset( new LocalOperators(_compute_operators(iT)) );
        } // for iT
      };

  m_output << "[HHOSpace] Constructing operators" << std::endl;
  parallel_for(mesh.n_cells(), construct_all_operators, use_threads);

}

//------------------------------------------------------------------------------
// Polynomial bases
//------------------------------------------------------------------------------

HHOSpace::CellBases HHOSpace::_construct_cell_bases(size_t iT)
{
  const Cell & T = *m_mesh.cell(iT);

  CellBases bases_T;
  
  MonomialCellIntegralsType int_monoT_2kp2 = IntegrateCellMonomials(T, 2*(m_K+1));
  
  //------------------------------------------------------------------------------
  // Basis for Pk+1(T), Pk(T)
  //------------------------------------------------------------------------------
  
  MonomialScalarBasisCell basis_Pkpo_T(T, m_K + 1);
  bases_T.Polykpo.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pkpo_T, GramMatrix(T, basis_Pkpo_T, int_monoT_2kp2))) );  

  MonomialScalarBasisCell basis_Pk_T(T, m_K);
  bases_T.Polyk.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pk_T, GramMatrix(T, basis_Pk_T, int_monoT_2kp2))) );  

  // Check that we got the dimensions right
  assert( bases_T.Polykpo->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K + 1) );
  assert( bases_T.Polyk->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K) );
  
  //------------------------------------------------------------------------------
  // Basis for Pk(T)^d
  //------------------------------------------------------------------------------

  bases_T.Polykd.reset( new PolydBasisCellType(*bases_T.Polyk) );
  assert( bases_T.Polykd->dimension() == dimspace * PolynomialSpaceDimension<Cell>::Poly(m_K) );
  
  return bases_T;
}

//------------------------------------------------------------------------------

HHOSpace::FaceBases HHOSpace::_construct_face_bases(size_t iF)
{
  const Face & F = *m_mesh.face(iF);
  
  FaceBases bases_F;

  MonomialFaceIntegralsType int_monoF_2k = IntegrateFaceMonomials(F, 2*m_K);

  //------------------------------------------------------------------------------
  // Basis for Pk(F)
  //------------------------------------------------------------------------------
  MonomialScalarBasisFace basis_Pk_F(F, m_K);
  bases_F.Polyk.reset( new PolyBasisFaceType(l2_orthonormalize(basis_Pk_F, GramMatrix(F, basis_Pk_F, int_monoF_2k))) );
  
  // Check that we got the dimensions right
  assert( bases_F.Polyk->dimension() == PolynomialSpaceDimension<Face>::Poly(m_K) );

  return bases_F;
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd HHOSpace::interpolate(const FunctionType & q, const int doe_cell, const int doe_face) const
{
  Eigen::VectorXd qh = Eigen::VectorXd::Zero(dimension());

  // Degrees of quadrature rules
  size_t dqr_cell = (doe_cell >= 0 ? doe_cell : 2 * degree() + 3);
  size_t dqr_face = (doe_face >= 0 ? doe_face : 2 * degree() + 3);
      
  // Interpolate at faces
  std::function<void(size_t, size_t)> interpolate_faces
    = [this, &qh, q, &dqr_face](size_t start, size_t end)->void
      {
        for (size_t iF = start; iF < end; iF++) {
          const Face & F = *mesh().face(iF);
          QuadratureRule quad_dqr_F = generate_quadrature_rule(F, dqr_face);
          auto basis_Pk_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).Polyk, quad_dqr_F);
          qh.segment(globalOffset(F), PolynomialSpaceDimension<Face>::Poly(degree())) 
            = l2_projection(q, *faceBases(iF).Polyk, quad_dqr_F, basis_Pk_F_quad, GramMatrix(F, *faceBases(iF).Polyk));
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
          auto basis_Pk_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk, quad_dqr_T);
          qh.segment(globalOffset(T), PolynomialSpaceDimension<Cell>::Poly(degree())) 
            = l2_projection(q, *cellBases(iT).Polyk, quad_dqr_T, basis_Pk_T_quad, GramMatrix(T, *cellBases(iT).Polyk));
        } // for iT
      };
  parallel_for(mesh().n_cells(), interpolate_cells, m_use_threads);


  return qh;
}

//------------------------------------------------------------------------------
// Operators
//------------------------------------------------------------------------------

HHOSpace::LocalOperators HHOSpace::_compute_operators(size_t iT)
{
  const Cell & T = *mesh().cell(iT);

  // Dimension
  size_t dim_Pkpo = PolynomialSpaceDimension<Cell>::Poly(degree()+1);

  //------------------------------------------------------------------------------
  // Gradient
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix

  // Compute all integrals of monomial powers to degree 2k+3 and the mass matrix
  MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2);
  Eigen::MatrixXd MGT = GramMatrix(T, *cellBases(iT).Polykd, int_mono_2kp2);

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BGT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polykd->dimension(), dimensionCell(iT));

  // Boundary contribution
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    
    DecomposePoly dec(F, MonomialScalarBasisFace(F, degree()));
    auto PkdT_dot_nTF_nodes = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polykd, dec.get_nodes()), T.face_normal(iF));
    auto PkdT_dot_nTF_family_PkF = dec.family(PkdT_dot_nTF_nodes);
    // PF extracts the face unknowns corresponding to F when read among all the element & face unknowns of T
    Eigen::MatrixXd PF = extendOperator(T, F, Eigen::MatrixXd::Identity(dimensionFace(F), dimensionFace(F)));
    MonomialFaceIntegralsType int_mono_2k_F = IntegrateFaceMonomials(F, 2*degree());
    BGT += GramMatrix(F, PkdT_dot_nTF_family_PkF, *faceBases(F).Polyk, int_mono_2k_F) * PF;

    // Following commented block could replace the block above, without using DecomposePoly (more expensive, but sometimes better rounding errors) -- not tested here yet
    /*
      QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree() );
      auto basis_Pkd_nTF_F_quad = scalar_product(
					         evaluate_quad<Function>::compute(*cellBases(iT).Polykd, quad_2k_F),
					         T.face_normal(iF)
					         );
      auto basis_Pk_F_quad = evaluate_quad<Function>::compute(*faceBases(F).Polyk, quad_2k_F);
      Eigen::MatrixXd PF = extendOperator(T, F, Eigen::MatrixXd::Identity(dimensionFace(F), dimensionFace(F)));
      BGT += compute_gram_matrix(basis_Pkd_nTF_F_quad, basis_Pk_F_quad, quad_2k_F) * PF;
    */
    
  } // for iF

  // Cell contribution
  DivergenceBasis<HHOSpace::PolydBasisCellType> div_Pkd_basis(*cellBases(iT).Polykd);
  BGT.rightCols(numLocalDofsCell()) -= GramMatrix(T, div_Pkd_basis, *cellBases(iT).Polyk, int_mono_2kp2);

  Eigen::MatrixXd GT = MGT.ldlt().solve(BGT);
  
  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------

  // We write \nabla pT = projection on \nabla P^{k+1} of GT, and add the closure relation:
  //    (nabla pT v, nabla w)_T + lambda_T(p_T v,1)_T(w,1)_T = (GT v,\nabla w)_T + lambda_T(v_T,1)_T(w,1)_T

  // LHS
  GradientBasis<HHOSpace::PolyBasisCellType> GradPolykpo = GradientBasis(*cellBases(iT).Polykpo);
  Eigen::MatrixXd StiffT = GramMatrix(T, GradPolykpo, int_mono_2kp2);

  // RHS, starting with volumetric term
  Eigen::MatrixXd RHS_PT = GramMatrix(T, GradPolykpo, *cellBases(iT).Polykd, int_mono_2kp2) * GT;
  
  // Vector LT of (phi_j,1)_T for phi_j up to degree K+1, and LT^t*LT, for the closure relation
  Eigen::VectorXd LT = GramMatrix(T, *cellBases(iT).Polykpo, RestrictedBasis(*cellBases(iT).Polykpo, 1), int_mono_2kp2);
  Eigen::MatrixXd LTtLT = LT * (LT.transpose());
  double scalT = StiffT.trace() / std::pow(LT.norm(), 2);
 
  // Add closure relation and compute PT
  RHS_PT.block(
              0, localOffset(T), 
              dim_Pkpo, numLocalDofsCell() 
              ) += 
      scalT * LTtLT.topLeftCorner( dim_Pkpo, numLocalDofsCell() );
  Eigen::MatrixXd PT = ((StiffT + scalT*LTtLT).ldlt()).solve(RHS_PT);

  //------------------------------------------------------------------------------
  // Stabilisation.
  //    It is based on bilinear forms on the local space, and the operators Id - I_T P_T
  //  The bilinear forms built here are (depending on choice_bilinear):
  //        1) 'components_norm': ||v_T||_T^2 + h_T\sum_F ||v_F||_F^2 with scaling depending on regT
  //        2) 'original_hho': original hho one, h_T \sum_F ||v_F-v_T||_F^2
  //------------------------------------------------------------------------------
  std::string choice_bilinear = "original_hho";

  Eigen::MatrixXd BilinearForm = Eigen::MatrixXd::Zero(dimensionCell(iT), dimensionCell(iT));
  Eigen::MatrixXd Id_minus_ITPT = Eigen::MatrixXd::Identity(dimensionCell(iT), dimensionCell(iT));

  // Element contributions
  Eigen::MatrixXd GramTk = GramMatrix(T, *cellBases(iT).Polyk, int_mono_2kp2);
  Eigen::MatrixXd GramTk_kpo = GramMatrix(T, *cellBases(iT).Polyk, *cellBases(iT).Polykpo, int_mono_2kp2);

  // Element contribution to bilinear form only for the components norm
  if (choice_bilinear == "components_norm"){
    double regT = std::pow(T.diam(), dimspace)/T.measure() * T.n_faces();
    BilinearForm.bottomRightCorner(numLocalDofsCell(), numLocalDofsCell()) +=  regT * GramTk;
  }
  Id_minus_ITPT.bottomRows(numLocalDofsCell()) -= GramTk.ldlt().solve(GramTk_kpo * PT);

  // Face contributions
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    
    DecomposePoly dec(F, MonomialScalarBasisFace(F, degree()+1));
    auto PkpoT_nodes = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, dec.get_nodes());
    auto PkpoT_family_PkF = dec.family(PkpoT_nodes);

    MonomialFaceIntegralsType int_mono_2kpo_F = IntegrateFaceMonomials(F, 2*degree()+1);
    Eigen::MatrixXd GramFk = GramMatrix(F, *faceBases(F).Polyk, int_mono_2kpo_F);
    Eigen::MatrixXd GramFk_Tkpo = GramMatrix(F, *faceBases(F).Polyk, PkpoT_family_PkF, int_mono_2kpo_F);
    Eigen::LDLT<Eigen::MatrixXd> GramFk_inv(GramFk);
    GramFk_inv.compute(GramFk);
    
    // PF extracts the face unknowns corresponding to F when read among all the element & face unknowns of T
    Eigen::MatrixXd PF = extendOperator(T, F, Eigen::MatrixXd::Identity(dimensionFace(F), dimensionFace(F)));

    // For the original hho stabilisation, PF must correspond to v_F - v_T (the latter being projected on P^k(F)^d)
    if (choice_bilinear == "original_hho"){
      auto PkT_nodes = evaluate_quad<Function>::compute(*cellBases(iT).Polyk, dec.get_nodes());
      auto PkT_family_PkF = dec.family(PkT_nodes);
      Eigen::MatrixXd GramFk_Tk = GramMatrix(F, *faceBases(F).Polyk, PkT_family_PkF, int_mono_2kpo_F);

      PF.rightCols(numLocalDofsCell()) -= GramFk_inv.solve(GramFk_Tk);
    }
    
    BilinearForm += T.diam() * PF.transpose() * GramFk * PF;
    Id_minus_ITPT.middleRows(localOffset(T, F), numLocalDofsFace()) -= GramFk_inv.solve(GramFk_Tkpo * PT);

  } // for iF

  // Creation of the stabilisation
  Eigen::MatrixXd ST = std::pow(T.diam(), -2) * Id_minus_ITPT.transpose() * BilinearForm * Id_minus_ITPT;
  

  return LocalOperators(GT, PT, ST);
}




//------------------------------------------------------------------------------
// Norms
//------------------------------------------------------------------------------

std::vector<std::pair<double,double>> HHOSpace::computeNorms( const std::vector<Eigen::VectorXd> & list_dofs ) const
{
  size_t nb_vectors = list_dofs.size();
  std::vector<Eigen::VectorXd> local_L2_sqnorms(nb_vectors, Eigen::VectorXd::Zero(mesh().n_cells()));
  std::vector<Eigen::VectorXd> local_H1_sqnorms(nb_vectors, Eigen::VectorXd::Zero(mesh().n_cells()));

  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &list_dofs, &local_L2_sqnorms, &local_H1_sqnorms, &nb_vectors](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *mesh().cell(iT);

        // Mass matrices
	      MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*degree());
        Eigen::MatrixXd mass_PkT = GramMatrix(T, *cellBases(iT).Polyk, int_mono_2k);
        Eigen::MatrixXd mass_GradPkT = GramMatrix(T, GradientBasis<HHOSpace::PolyBasisCellType>(*cellBases(iT).Polyk), int_mono_2k);

        // Cell contributions
        for (size_t i=0; i<nb_vectors; i++){
          Eigen::VectorXd uT_cell = restrict(T, list_dofs[i]).tail(numLocalDofsCell());

          // L2 norm
          local_L2_sqnorms[i](iT) += uT_cell.transpose() * mass_PkT * uT_cell;
          // H1 norm
          local_H1_sqnorms[i](iT) += uT_cell.transpose() * mass_GradPkT * uT_cell;
        }
    
        // Face contributions
        for (size_t iF=0; iF<T.n_faces(); iF++){
          Face & F = *T.face(iF);

          MonomialFaceIntegralsType int_mono_2kF = IntegrateFaceMonomials(F, 2*degree());
          Eigen::MatrixXd mass_PkF = GramMatrix(F, *faceBases(F).Polyk, int_mono_2kF);

          DecomposePoly dec(F, MonomialScalarBasisFace(F, degree()));
          auto PkT_nodes = evaluate_quad<Function>::compute(*cellBases(iT).Polyk, dec.get_nodes());
          auto PkT_family_PkF = dec.family(PkT_nodes);
          Eigen::MatrixXd mass_PkF_PkT = GramMatrix(F, *faceBases(F).Polyk, PkT_family_PkF, int_mono_2kF);
          
          for (size_t i=0; i<nb_vectors; i++){
            Eigen::VectorXd uT = restrict(T, list_dofs[i]);

            // Decompose u_F-u_T on Polyk basis of F
            Eigen::VectorXd uF_minus_uT = 
                    uT.segment(iF*numLocalDofsFace(), numLocalDofsFace())
                    -
                    mass_PkF.ldlt().solve(mass_PkF_PkT * uT.tail(numLocalDofsCell()));
            // Contribution of gradients, without any weight (no permeability)
            double sqnorm_uF_minus_uT = uF_minus_uT.transpose() * mass_PkF * uF_minus_uT;
            local_H1_sqnorms[i](iT) += sqnorm_uF_minus_uT/T.diam();
          }
        } // for iF
      }
    };
  parallel_for(mesh().n_cells(), compute_local_squarednorms, m_use_threads);
  
  // Vector of outputs
  std::vector<std::pair<double,double>> list_norms(nb_vectors);
  for (size_t i=0; i<nb_vectors; i++){
    list_norms[i].first = std::sqrt(std::abs(local_L2_sqnorms[i].sum()));
    list_norms[i].second = std::sqrt(std::abs(local_H1_sqnorms[i].sum()));
  }
  
  return list_norms;
}


