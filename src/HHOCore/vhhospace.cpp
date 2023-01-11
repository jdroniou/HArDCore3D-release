#include <cassert>

#include <vhhospace.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_face.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------

VHHOSpace::VHHOSpace(const Mesh & mesh, size_t K, const CellSelection & BoundaryStab, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(mesh,
	     0,
	     0,
	     dimspace * PolynomialSpaceDimension<Face>::Poly(K),
	     dimspace * PolynomialSpaceDimension<Cell>::Poly(K)
	     ),
	  m_mesh(mesh),
    m_K(K),
    m_boundary_stab(BoundaryStab),
    m_use_threads(use_threads),
    m_output(output),
    m_cell_bases(mesh.n_cells()),
    m_face_bases(mesh.n_faces()),
    m_operators(mesh.n_cells())
{
  m_output << "[VHHOSpace] Initializing" << std::endl;
  
  // Construct element bases
  std::function<void(size_t, size_t)> construct_all_cell_bases
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iT = start; iT < end; iT++) {
	        this->m_cell_bases[iT].reset( new CellBases(this->_construct_cell_bases(iT)) );
	      } // for iT
      };

  m_output << "[VHHOSpace] Constructing element bases" << std::endl;
  parallel_for(mesh.n_cells(), construct_all_cell_bases, use_threads);
  
  // Construct face bases
  std::function<void(size_t, size_t)> construct_all_face_bases
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iF = start; iF < end; iF++) {
	        this->m_face_bases[iF].reset( new FaceBases(_construct_face_bases(iF)) );
	      } // for iF
      };
  
  m_output << "[VHHOSpace] Constructing face bases" << std::endl;
  parallel_for(mesh.n_faces(), construct_all_face_bases, use_threads);

  // Construct gradients, potentials and stabilisation
  std::function<void(size_t, size_t)> construct_all_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_operators[iT].reset( new LocalOperators(_compute_operators(iT)) );
        } // for iT
      };

  m_output << "[VHHOSpace] Constructing operators" << std::endl;
  parallel_for(mesh.n_cells(), construct_all_operators, use_threads);

}

//------------------------------------------------------------------------------
// Polynomial bases
//------------------------------------------------------------------------------

VHHOSpace::CellBases VHHOSpace::_construct_cell_bases(size_t iT)
{
  const Cell & T = *m_mesh.cell(iT);

  CellBases bases_T;
  
  MonomialCellIntegralsType int_monoT_2kp2 = IntegrateCellMonomials(T, 2*(m_K+1));
  
  //------------------------------------------------------------------------------
  // Basis for Pk+1(T), Pk(T), and vector versions
  //------------------------------------------------------------------------------
  
  MonomialScalarBasisCell basis_Pkpo_T(T, m_K + 1);
  bases_T.Polykpo.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pkpo_T, GramMatrix(T, basis_Pkpo_T, int_monoT_2kp2))) );
  bases_T.Polykpod.reset( new PolydBasisCellType(*bases_T.Polykpo) );  

  MonomialScalarBasisCell basis_Pk_T(T, m_K);
  bases_T.Polyk.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pk_T, GramMatrix(T, basis_Pk_T, int_monoT_2kp2))) );  
  bases_T.Polykd.reset( new PolydBasisCellType(*bases_T.Polyk) );  

  // Check that we got the dimensions right
  assert( bases_T.Polykpo->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K + 1) );
  assert( bases_T.Polykpod->dimension() == dimspace * PolynomialSpaceDimension<Cell>::Poly(m_K + 1) );
  assert( bases_T.Polyk->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K) );
  assert( bases_T.Polykd->dimension() == dimspace * PolynomialSpaceDimension<Cell>::Poly(m_K) );
  
  //------------------------------------------------------------------------------
  // Basis for Pk(T)^{dxd}
  //------------------------------------------------------------------------------

  bases_T.Polykdxd.reset( new PolydxdBasisCellType(*bases_T.Polyk) );
  assert( bases_T.Polykdxd->dimension() == dimspace * dimspace * PolynomialSpaceDimension<Cell>::Poly(m_K) );
  
  return bases_T;
}

//------------------------------------------------------------------------------

VHHOSpace::FaceBases VHHOSpace::_construct_face_bases(size_t iF)
{
  const Face & F = *m_mesh.face(iF);
  
  FaceBases bases_F;

  MonomialFaceIntegralsType int_monoF_2k = IntegrateFaceMonomials(F, 2*m_K);

  //------------------------------------------------------------------------------
  // Basis for Pk(F)
  //------------------------------------------------------------------------------
  MonomialScalarBasisFace basis_Pk_F(F, m_K);
  bases_F.Polyk.reset( new PolyBasisFaceType(l2_orthonormalize(basis_Pk_F, GramMatrix(F, basis_Pk_F, int_monoF_2k))) );

  //------------------------------------------------------------------------------
  // Basis for Pk(F)^d
  //------------------------------------------------------------------------------
  // Old version: tensorized on the canonical basis of R^3 (requires to change the class PolydBasisFaceType)
  // bases_F.Polykd.reset( new PolydBasisFaceType(*bases_F.Polyk) );
  
  // For internal faces, the basis of Pk(F)^d$ is just a tensorized one, that we deal as a (trivial) Family<Tensorized> to match the expected class.
  // For boundary faces, the basis of Pk(F)^d is designed such that, if r=dim Pk(F), the first 2r polynomials are tangent to F and the last r polynomials are orthogonal to F
  if (!F.is_boundary()){
    TensorizedVectorFamily<PolyBasisFaceType, dimspace> basis_tens_PolykdF(*bases_F.Polyk);
    bases_F.Polykd.reset( 
              new PolydBasisFaceType(basis_tens_PolykdF,
                                     Eigen::MatrixXd::Identity(basis_tens_PolykdF.dimension(), basis_tens_PolykdF.dimension()) )
                        );                                     
  }else{
    bases_F.Polykd.reset( 
        new PolydBasisFaceType( 
            GenericTensorization<PolyBasisFaceType, dimspace>(*bases_F.Polyk, std::vector<Eigen::VectorXd> {F.edge(0)->tangent(), F.edge_normal(0), F.normal()})
            )
          );
  }
  
  // Check that we got the dimensions right
  assert( bases_F.Polyk->dimension() == PolynomialSpaceDimension<Face>::Poly(m_K) );
  assert( bases_F.Polykd->dimension() == dimspace * PolynomialSpaceDimension<Face>::Poly(m_K) );

  return bases_F;
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd VHHOSpace::interpolate(const FunctionType & q, const int doe_cell, const int doe_face) const
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
          auto basis_Pkd_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).Polykd, quad_dqr_F);
          qh.segment(globalOffset(F), numLocalDofsFace()) 
            = l2_projection(q, *faceBases(iF).Polykd, quad_dqr_F, basis_Pkd_F_quad, GramMatrix(F, *faceBases(iF).Polykd));
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
          auto basis_Pkd_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykd, quad_dqr_T);
          qh.segment(globalOffset(T), numLocalDofsCell()) 
            = l2_projection(q, *cellBases(iT).Polykd, quad_dqr_T, basis_Pkd_T_quad, GramMatrix(T, *cellBases(iT).Polykd));
        } // for iT
      };
  parallel_for(mesh().n_cells(), interpolate_cells, m_use_threads);


  return qh;
}

//------------------------------------------------------------------------------
// Operators
//------------------------------------------------------------------------------

VHHOSpace::LocalOperators VHHOSpace::_compute_operators(size_t iT)
{
  const Cell & T = *mesh().cell(iT);

  // Dimension
  size_t dim_Pkpod = dimspace * PolynomialSpaceDimension<Cell>::Poly(degree()+1);

  //------------------------------------------------------------------------------
  // Gradient
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix

  // Compute all integrals of monomial powers to degree 2k+3 and the mass matrix
  MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2);
  Eigen::MatrixXd MGT = GramMatrix(T, *cellBases(iT).Polykdxd, int_mono_2kp2);

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BGT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polykdxd->dimension(), dimensionCell(iT));

  // Boundary contribution
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    
    DecomposePoly dec(F, TensorizedVectorFamily<MonomialScalarBasisFace, dimspace>(MonomialScalarBasisFace(F, degree())));
    VectorRd nTF = T.face_normal(iF);
    auto PkdxdT_nTF_nodes = transform_values_quad<VectorRd>(evaluate_quad<Function>::compute(*cellBases(iT).Polykdxd, dec.get_nodes()), [&nTF](const MatrixRd &M)->VectorRd { return M*nTF;});
    auto PkdxdT_nTF_family_PkdF = dec.family(PkdxdT_nTF_nodes);
    // PF extracts the face unknowns corresponding to F when read among all the element & face unknowns of T
    Eigen::MatrixXd PF = extendOperator(T, F, Eigen::MatrixXd::Identity(dimensionFace(F), dimensionFace(F)));
    MonomialFaceIntegralsType int_mono_2k_F = IntegrateFaceMonomials(F, 2*degree());
    BGT += GramMatrix(F, PkdxdT_nTF_family_PkdF, *faceBases(F).Polykd, int_mono_2k_F) * PF;

    // Following commented block could replace the block above, without using DecomposePoly (more expensive, but sometimes better rounding errors) -- not tested here yet
    /*
      QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree() );
      auto basis_Pkdxd_nTF_F_quad = transform_values_quad<VectorRd>(
					         evaluate_quad<Function>::compute(*cellBases(iT).Polykdxd, quad_2k_F),
					         [&nTF](MatrixRd &M)->VectorRd { return M*nTF;}
					         );
      auto basis_Pkd_F_quad = evaluate_quad<Function>::compute(*faceBases(F).Polykd, quad_2k_F);
      Eigen::MatrixXd PF = extendOperator(T, F, Eigen::MatrixXd::Identity(dimensionFace(F), dimensionFace(F)));
      BGT += compute_gram_matrix(basis_Pkdxd_nTF_F_quad, basis_Pkd_F_quad, quad_2k_F) * PF;
    */
    
  } // for iF

  // Cell contribution
  DivergenceBasis<VHHOSpace::PolydxdBasisCellType> div_Pkdxd_basis(*cellBases(iT).Polykdxd);
  BGT.rightCols(numLocalDofsCell()) -= GramMatrix(T, div_Pkdxd_basis, *cellBases(iT).Polykd, int_mono_2kp2);

  Eigen::MatrixXd GT = MGT.ldlt().solve(BGT);
  
  //------------------------------------------------------------------------------
  // HHO Potential
  //------------------------------------------------------------------------------

  // We write \nabla pT = projection on \nabla P^{k+1} of GT, and add the closure relation:
  //    (nabla pT v, nabla w)_T + lambda_T(p_T v,1)_T(w,1)_T = (GT v,\nabla w)_T + lambda_T(v_T,1)_T(w,1)_T

  // LHS
  GradientBasis<VHHOSpace::PolydBasisCellType> GradPolykpod = GradientBasis(*cellBases(iT).Polykpod);
  Eigen::MatrixXd StiffT = GramMatrix(T, GradPolykpod, int_mono_2kp2);

  // RHS, starting with volumetric term
  Eigen::MatrixXd RHS_PT = GramMatrix(T, GradPolykpod, *cellBases(iT).Polykdxd, int_mono_2kp2) * GT;
  
  // Closure matrix (\int phi_i).(\int phi_j) computed as \sum_k \int phi_i.e_k \int phi_j.e_k, with e_k basis of R^d represented by a polynomial basis for P^0(T)^d
  auto basis_P0d = TensorizedVectorFamily<MonomialScalarBasisCell, dimspace>(MonomialScalarBasisCell(T, 0));
  Eigen::MatrixXd Integrate_Pkpod_P0d = GramMatrix(T, *cellBases(iT).Polykpod, basis_P0d, int_mono_2kp2);
  Eigen::MatrixXd Integrate_Pkd_P0d = GramMatrix(T, *cellBases(iT).Polykd, basis_P0d, int_mono_2kp2);

  Eigen::MatrixXd Closure = Integrate_Pkpod_P0d * Integrate_Pkpod_P0d.transpose();
  double scalT = StiffT.trace() / Closure.trace();

  RHS_PT.block(
              0, localOffset(T), 
              dim_Pkpod, numLocalDofsCell() 
              ) += scalT * (Integrate_Pkpod_P0d * Integrate_Pkd_P0d.transpose());

  Eigen::MatrixXd PT = ((StiffT + scalT * Closure).ldlt()).solve(RHS_PT);
  
  //------------------------------------------------------------------------------
  // Potential linked to divergence
  //------------------------------------------------------------------------------
  size_t dim_Pkd = cellBases(iT).Polykd->dimension();
  Eigen::MatrixXd LHS_Pdiv = Eigen::MatrixXd::Zero(dim_Pkd, dim_Pkd);
  Eigen::MatrixXd RHS_Pdiv = Eigen::MatrixXd::Zero(dim_Pkd, dimensionCell(iT));

  // LHS - gradient component

	MonomialCellIntegralsType int_mono_2kpo = IntegrateCellMonomials(T, 2*degree()+1);
	auto basis_Pkpo0 = ShiftedBasis<PolyBasisCellType>(*cellBases(iT).Polykpo, 1);
  auto basis_gradPkpo0 = GradientBasis(basis_Pkpo0);
  LHS_Pdiv.topRows(basis_gradPkpo0.dimension()) = GramMatrix(T, basis_gradPkpo0, *cellBases(iT).Polykd, int_mono_2kpo);
    
  // RHS - gradient component
  auto basis_qId = IsotropicMatrixFamily<ShiftedBasis<PolyBasisCellType>, dimspace>(basis_Pkpo0);
  RHS_Pdiv.topRows(basis_Pkpo0.dimension()) = - GramMatrix(T, basis_qId, *cellBases(iT).Polykdxd, int_mono_2kpo) * GT;
  for (size_t iF = 0; iF < T.n_faces(); iF++){
    const Face & F = *T.face(iF);
    VectorRd nF = F.normal();
    
    // Decompose q n_F on a vector polynomial basis on the face
    DecomposePoly dec(F, TensorizedVectorFamily<MonomialScalarBasisFace, dimspace>(MonomialScalarBasisFace(F, degree()+1)));
    auto Pkpo0T_nF_nodes = transform_values_quad<VectorRd>(
                                    evaluate_quad<Function>::compute(basis_Pkpo0, dec.get_nodes()),
                                    [&nF](const double & q)->VectorRd { return q*nF;}
                                    );
    auto Pkpo0T_nF_family_PkpoF = dec.family(Pkpo0T_nF_nodes);
    
    MonomialFaceIntegralsType int_mono_2kpo_F = IntegrateFaceMonomials(F, 2*degree()+1);
    RHS_Pdiv.block(0, localOffset(T, F), basis_Pkpo0.dimension(), dimensionFace(F))
      += T.face_orientation(iF) * GramMatrix(F, Pkpo0T_nF_family_PkpoF, *faceBases(F).Polykd, int_mono_2kpo_F);
  }
  
  // LHS and RHS - Gck component
  if (degree()>0){
    GolyComplBasisCell basis_Gck_T(T, degree());
    Family<GolyComplBasisCell> on_basis_Gck_T(l2_orthonormalize(basis_Gck_T, GramMatrix(T, basis_Gck_T, int_mono_2kpo)));
  
    LHS_Pdiv.bottomRows(on_basis_Gck_T.dimension()) = GramMatrix(T, on_basis_Gck_T, *cellBases(iT).Polykd, int_mono_2kpo); 

    RHS_Pdiv.block(basis_Pkpo0.dimension(), localOffset(T), on_basis_Gck_T.dimension(), numLocalDofsCell()) 
      = GramMatrix(T, on_basis_Gck_T, *cellBases(iT).Polykd, int_mono_2kpo);
  }
  
  Eigen::MatrixXd PT_div = LHS_Pdiv.partialPivLu().solve(RHS_Pdiv);

  
  //------------------------------------------------------------------------------
  // Stabilisations: HHO and coming from the divergence.
  //    They are both based on bilinear forms on the local space, and operators Id - I_T P_T
  //    (with P_T the HHO or divergence potential)
  //  The bilinear forms for HHO built here are (depending on choice_bilinear):
  //        1) 'components_norm': ||v_T||_T^2 + h_T\sum_F ||v_F||_F^2 with scaling depending on regT
  //        2) 'original_hho': original hho one, h_T \sum_F ||v_F-v_T||_F^2
  //  For the divergence stabilisation, it is always 'components_norm'
  //------------------------------------------------------------------------------
  std::string choice_bilinear = "components_norm";

  Eigen::MatrixXd BilinearForm = Eigen::MatrixXd::Zero(dimensionCell(iT), dimensionCell(iT));
  Eigen::MatrixXd Id_minus_ITPT = Eigen::MatrixXd::Identity(dimensionCell(iT), dimensionCell(iT));

  Eigen::MatrixXd BilinearForm_div = Eigen::MatrixXd::Zero(dimensionCell(iT), dimensionCell(iT));
  Eigen::MatrixXd Id_minus_ITPT_div = Eigen::MatrixXd::Identity(dimensionCell(iT), dimensionCell(iT));

  // Element contributions
  Eigen::MatrixXd GramTkd = GramMatrix(T, *cellBases(iT).Polykd, int_mono_2kp2);
  Eigen::MatrixXd GramTkd_kpod = GramMatrix(T, *cellBases(iT).Polykd, *cellBases(iT).Polykpod, int_mono_2kp2);

  double regT = std::pow(T.diam(), dimspace)/T.measure() * T.n_faces();
  // Element contribution to bilinear form only for the components norm
  if (choice_bilinear == "components_norm"){
    BilinearForm.bottomRightCorner(numLocalDofsCell(), numLocalDofsCell()) +=  regT * GramTkd;
  }
  Id_minus_ITPT.bottomRows(numLocalDofsCell()) -= GramTkd.ldlt().solve(GramTkd_kpod * PT);

  BilinearForm_div.bottomRightCorner(numLocalDofsCell(), numLocalDofsCell()) +=  regT * GramTkd;
  Id_minus_ITPT_div.bottomRows(numLocalDofsCell()) -= PT_div;

  // Face contributions
  for (size_t iF = 0; iF < T.n_faces(); iF++) {
    const Face & F = *T.face(iF);
    
    // For boundary faces, we include the contribution to the inner product only if m_boundary_stab says so (we do not complete
    // the construction of Id_minus_ITPT for these faces because that part will not be used by BilinearForm; so Id_minus_ITPT
    // is not empty on these faces but incorrect, only contain u_F)
    if ( !F.is_boundary() || m_boundary_stab(T) ){    
      DecomposePoly dec(F, TensorizedVectorFamily<MonomialScalarBasisFace, dimspace>(MonomialScalarBasisFace(F, degree()+1)));
      auto PkpodT_nodes = evaluate_quad<Function>::compute(*cellBases(iT).Polykpod, dec.get_nodes());
      auto PkpodT_family_PkdF = dec.family(PkpodT_nodes);
      auto PkdT_nodes = evaluate_quad<Function>::compute(*cellBases(iT).Polykd, dec.get_nodes());
      auto PkdT_family_PkdF = dec.family(PkdT_nodes);

      MonomialFaceIntegralsType int_mono_2kpo_F = IntegrateFaceMonomials(F, 2*degree()+1);
      Eigen::MatrixXd GramFkd = GramMatrix(F, *faceBases(F).Polykd, int_mono_2kpo_F);
      Eigen::MatrixXd GramFkd_Tkpod = GramMatrix(F, *faceBases(F).Polykd, PkpodT_family_PkdF, int_mono_2kpo_F);
      Eigen::MatrixXd GramFkd_Tkd = GramMatrix(F, *faceBases(F).Polykd, PkdT_family_PkdF, int_mono_2kpo_F);
      Eigen::LDLT<Eigen::MatrixXd> GramFkd_inv(GramFkd);
      GramFkd_inv.compute(GramFkd);
      
      // PF extracts the face unknowns corresponding to F when read among all the element & face unknowns of T
      Eigen::MatrixXd PF = extendOperator(T, F, Eigen::MatrixXd::Identity(dimensionFace(F), dimensionFace(F)));

      BilinearForm_div += T.diam() * PF.transpose() * GramFkd * PF;
      Id_minus_ITPT_div.middleRows(localOffset(T, F), numLocalDofsFace()) -= GramFkd_inv.solve(GramFkd_Tkd * PT_div);

      // For the original hho stabilisation, PF must correspond to v_F - v_T (the latter being projected on P^k(F)^d)
      if (choice_bilinear == "original_hho"){
        PF.rightCols(numLocalDofsCell()) -= GramFkd_inv.solve(GramFkd_Tkd);
      }
      
      BilinearForm += T.diam() * PF.transpose() * GramFkd * PF;
      Id_minus_ITPT.middleRows(localOffset(T, F), numLocalDofsFace()) -= GramFkd_inv.solve(GramFkd_Tkpod * PT);

    }    
  } // for iF

  // Creation of the stabilisations
  Eigen::MatrixXd ST = std::pow(T.diam(), -2) * Id_minus_ITPT.transpose() * BilinearForm * Id_minus_ITPT;
  Eigen::MatrixXd ST_div = Id_minus_ITPT_div.transpose() * BilinearForm_div * Id_minus_ITPT_div;
  

  return LocalOperators(GT, PT, ST, PT_div, ST_div);
}




//------------------------------------------------------------------------------
// Norms
//------------------------------------------------------------------------------

std::vector<std::pair<double,double>> VHHOSpace::computeNorms( const std::vector<Eigen::VectorXd> & list_dofs ) const
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
        Eigen::MatrixXd mass_PkdT = GramMatrix(T, *cellBases(iT).Polykd, int_mono_2k);
        Eigen::MatrixXd mass_GradPkdT = GramMatrix(T, GradientBasis<VHHOSpace::PolydBasisCellType>(*cellBases(iT).Polykd), int_mono_2k);

        // Local vectors
        std::vector<Eigen::VectorXd> uT(nb_vectors, Eigen::VectorXd::Zero(dimensionCell(iT)));
        for (size_t i=0; i < nb_vectors; i++){
          uT[i] = restrict(T, list_dofs[i]);
        }

        // Cell contributions
        for (size_t i=0; i<nb_vectors; i++){
          Eigen::VectorXd uT_cell = uT[i].tail(numLocalDofsCell());

          // L2 norm
          local_L2_sqnorms[i](iT) += uT_cell.transpose() * mass_PkdT * uT_cell;
          // H1 norm
          local_H1_sqnorms[i](iT) += uT_cell.transpose() * mass_GradPkdT * uT_cell;
        }
    
        // Face contributions
        for (size_t iF=0; iF<T.n_faces(); iF++){
          Face & F = *T.face(iF);

          MonomialFaceIntegralsType int_mono_2kF = IntegrateFaceMonomials(F, 2*degree());
          Eigen::MatrixXd mass_PkdF = GramMatrix(F, *faceBases(F).Polykd, int_mono_2kF);

          DecomposePoly dec(F, TensorizedVectorFamily<MonomialScalarBasisFace, dimspace>(MonomialScalarBasisFace(F, degree())));
          auto PkdT_nodes = evaluate_quad<Function>::compute(*cellBases(iT).Polykd, dec.get_nodes());
          auto PkdT_family_PkdF = dec.family(PkdT_nodes);
          Eigen::MatrixXd mass_PkdF_PkdT = GramMatrix(F, *faceBases(F).Polykd, PkdT_family_PkdF, int_mono_2kF);
          
          for (size_t i=0; i<nb_vectors; i++){
            // Decompose u_F-u_T on Polyk basis of F
            Eigen::VectorXd uF_minus_uT = 
                    uT[i].segment(iF*numLocalDofsFace(), numLocalDofsFace())
                    -
                    mass_PkdF.ldlt().solve(mass_PkdF_PkdT * uT[i].tail(numLocalDofsCell()));
            // Contribution of gradients, without any weight (no permeability)
            double sqnorm_uF_minus_uT = uF_minus_uT.transpose() * mass_PkdF * uF_minus_uT;
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

//------------------------------------------------------------------------------
// Vertex values 
//------------------------------------------------------------------------------
std::vector<VectorRd> VHHOSpace::computeVertexValues(const Eigen::VectorXd & u) const
{
  std::vector<VectorRd> values(m_mesh.n_vertices(), VectorRd::Zero());
  
  // Value at each vertex obtained averaging the values from all the cells around
  for (Vertex * V : mesh().get_vertices()){
    size_t iV = V->global_index();
 
    for (Cell * T : V->get_cells()){
      size_t iT = T->global_index();
      Eigen::VectorXd pTuT = operators(iT).potential * restrict(*T, u);
      
      for (size_t i=0; i < cellBases(iT).Polykpod->dimension(); i++){
        values[iV] += pTuT(i) * cellBases(iT).Polykpod->function(i, V->coords());
      }
    }
    
    values[iV] /= V->n_cells();

  }
  
  return values;
}                  


