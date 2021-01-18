#include <cassert>

#include <ddrcore.hpp>
#include <parallel_for.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------

DDRCore::DDRCore(const Mesh & mesh, size_t K, bool use_threads, std::ostream & output)
  : m_mesh(mesh),
    m_K(K),
    m_output(output),
    m_cell_bases(mesh.n_cells()),
    m_face_bases(mesh.n_faces()),
    m_edge_bases(mesh.n_edges())
{
  m_output << "[DDRCore-ORTH] Initializing" << std::endl;

  // Construct element bases
  std::function<void(size_t, size_t)> construct_all_cell_bases
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iT = start; iT < end; iT++) {
	        this->m_cell_bases[iT].reset( new CellBases(this->_construct_cell_bases(iT)) );
	      } // for iT
      };

  m_output << "[DDRCore-ORTH] Constructing element bases" << std::endl;
  parallel_for(mesh.n_cells(), construct_all_cell_bases, use_threads);
  
  // Construct face bases
  std::function<void(size_t, size_t)> construct_all_face_bases
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iF = start; iF < end; iF++) {
	        this->m_face_bases[iF].reset( new FaceBases(_construct_face_bases(iF)) );
	      } // for iF
      };
  
  m_output << "[DDRCore-ORTH] Constructing face bases" << std::endl;
  parallel_for(mesh.n_faces(), construct_all_face_bases, use_threads);

  // Construct edge bases
  std::function<void(size_t, size_t)> construct_all_edge_bases   
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iE = start; iE < end; iE++) {
	        this->m_edge_bases[iE].reset( new EdgeBases(_construct_edge_bases(iE)) );
	      } // for iF
      };
  
  m_output << "[DDRCore-ORTH] Constructing edge bases" << std::endl;
  parallel_for(mesh.n_edges(), construct_all_edge_bases, use_threads);
}

//------------------------------------------------------------------------------

DDRCore::CellBases DDRCore::_construct_cell_bases(size_t iT)
{
  const Cell & T = *m_mesh.cell(iT);

  CellBases bases_T;
  
  //------------------------------------------------------------------------------
  // Basis for Pk+1(T)
  //------------------------------------------------------------------------------
  
  MonomialScalarBasisCell basis_Pkpo_T(T, m_K + 1);
  QuadratureRule quad_2kpo_T = generate_quadrature_rule(T, 2 * (m_K + 1));
  auto on_basis_Pkpo_T_quad = evaluate_quad<Function>::compute(basis_Pkpo_T, quad_2kpo_T);
  // Orthonormalize and store
  bases_T.Polykpo.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pkpo_T, quad_2kpo_T, on_basis_Pkpo_T_quad)) );   
  // Check that we got the dimension right
  assert( bases_T.Polykpo->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K + 1) );

  //------------------------------------------------------------------------------
  // Basis for Pk(T) and Pk-1(T)
  //------------------------------------------------------------------------------

  // Given that the basis for Pk+1(T) is hierarchical, bases for Pk(T) and
  // Pk-1(T) can be obtained by restricting the former
  bases_T.Polyk.reset( new RestrictedBasis<PolyBasisCellType>(*bases_T.Polykpo, PolynomialSpaceDimension<Cell>::Poly(m_K)) );  
  if (PolynomialSpaceDimension<Cell>::Poly(m_K - 1) > 0) {
    bases_T.Polykmo.reset( new RestrictedBasis<PolyBasisCellType>(*bases_T.Polykpo, PolynomialSpaceDimension<Cell>::Poly(m_K - 1)) );
  }
    
  //------------------------------------------------------------------------------
  // Basis for Pk(T)^3
  //------------------------------------------------------------------------------
  // Construct an orthonormal scalar basis for Pk(T) first
  MonomialScalarBasisCell basis_Pk_T(T, m_K);
  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * m_K);
  auto on_basis_Pk_T_quad = evaluate_quad<Function>::compute(basis_Pk_T, quad_2k_T);
  // After calling the following function, on_basis_Pk_T_ quad contains the values
  // of the orthonormalized basis at quadrature nodes
  PolyBasisCellType on_basis_Pk_T = l2_orthonormalize(basis_Pk_T, quad_2k_T, on_basis_Pk_T_quad);    
  // We obtain the basis for Pk(T)^3 by tensorization
  bases_T.Polyk3.reset( new Poly3BasisCellType(on_basis_Pk_T) );
  // Check that we got the dimension right
  assert( bases_T.Polyk3->dimension() == 3 * PolynomialSpaceDimension<Cell>::Poly(m_K) );
  
  //------------------------------------------------------------------------------
  // Basis for Gk-1(T)
  //------------------------------------------------------------------------------

  if (PolynomialSpaceDimension<Cell>::Goly(m_K - 1) > 0) {
    GradientBasis<ShiftedBasis<MonomialScalarBasisCell> >
      basis_Gkmo_T(ShiftedBasis<MonomialScalarBasisCell>(MonomialScalarBasisCell(T, m_K), 1));
    QuadratureRule quad_2kmo_T = generate_quadrature_rule(T, 2 * (m_K -1));  
    auto basis_Gkmo_T_quad = evaluate_quad<Function>::compute(basis_Gkmo_T, quad_2kmo_T);
    // Orthonormalize and store the basis
    bases_T.Golykmo.reset( new GolyBasisCellType(l2_orthonormalize(basis_Gkmo_T, quad_2kmo_T, basis_Gkmo_T_quad)) );
    // Check that we got the dimension right
    assert( bases_T.Golykmo->dimension() == PolynomialSpaceDimension<Cell>::Goly(m_K - 1) );
  } // if

  //------------------------------------------------------------------------------
  // Bases for GOk(T), GOk+1(T), and Rk-1(T)
  //------------------------------------------------------------------------------

  GradientBasis<ShiftedBasis<MonomialScalarBasisCell> >
    basis_Gkpo_T(ShiftedBasis<MonomialScalarBasisCell>(MonomialScalarBasisCell(T, m_K + 2), 1));
  auto basis_Gkpo_T_quad = evaluate_quad<Function>::compute(basis_Gkpo_T, quad_2kpo_T);
  GolyBasisCellType on_basis_Gkpo_T = l2_orthonormalize(basis_Gkpo_T, quad_2kpo_T, basis_Gkpo_T_quad);

  // Generate a basis for GOk+1(T) by computing the kernel of MGpo
  Poly3BasisCellType basis_Pkpo3(*bases_T.Polykpo);
  Eigen::MatrixXd MGpo = compute_gram_matrix(basis_Gkpo_T_quad, on_basis_Pkpo_T_quad, quad_2kpo_T);
  bases_T.GolyComplkpo.reset( new GolyComplBasisCellType(basis_Pkpo3, MGpo.fullPivLu().kernel().transpose()));
  // Check that we got the dimension right
  assert( bases_T.GolyComplkpo->dimension() == PolynomialSpaceDimension<Cell>::GolyCompl(m_K + 1) );

  if (PolynomialSpaceDimension<Cell>::GolyCompl(m_K) > 0) {
    // Create a basis for Gk(T)
    GradientBasis<ShiftedBasis<MonomialScalarBasisCell> >
      basis_Gk_T(ShiftedBasis<MonomialScalarBasisCell>(MonomialScalarBasisCell(T, m_K + 1), 1));
    auto basis_Gk_T_quad = evaluate_quad<Function>::compute(basis_Gk_T, quad_2k_T);
    GolyBasisCellType on_basis_Gk_T = l2_orthonormalize(basis_Gk_T, quad_2k_T, basis_Gk_T_quad);
    
    // Generate a basis for GOk(T) by computing the kernel of MG
    Eigen::MatrixXd MG = compute_gram_matrix(basis_Gk_T_quad, on_basis_Pk_T_quad, quad_2k_T);
    bases_T.GolyComplk.reset( new GolyComplBasisCellType(*bases_T.Polyk3, MG.fullPivLu().kernel().transpose()) );
    // Generate a basis for Rk-1(T) by taking the curl of functions in GOk(T)
    bases_T.Rolykmo.reset ( new RolyBasisCellType(*bases_T.GolyComplk) );
    // Check that we got the dimensions right
    assert( bases_T.GolyComplk->dimension() == PolynomialSpaceDimension<Cell>::GolyCompl(m_K) );
    assert( bases_T.Rolykmo->dimension() == PolynomialSpaceDimension<Cell>::Roly(m_K - 1));
  } // if
  
  //------------------------------------------------------------------------------
  // Basis for ROk(T)
  //------------------------------------------------------------------------------

  if (PolynomialSpaceDimension<Cell>::RolyCompl(m_K) > 0) {
    // Create a spanning set for Rk(T)
    CurlBasis<TensorizedVectorFamily<MonomialScalarBasisCell, 3> >
      spanning_set_Rk_T(TensorizedVectorFamily<MonomialScalarBasisCell, 3>(MonomialScalarBasisCell(T, m_K + 1)));
    auto spanning_set_Rk_T_quad = evaluate_quad<Function>::compute(spanning_set_Rk_T, quad_2k_T);
    
    // Obtain a basis for ROk(T) by computing the kernel of the matrix MR
    Eigen::MatrixXd MR = compute_gram_matrix(spanning_set_Rk_T_quad, on_basis_Pk_T_quad, quad_2k_T);
    bases_T.RolyComplk.reset( new RolyComplBasisCellType(*bases_T.Polyk3, MR.fullPivLu().kernel().transpose()) );
    // Check that we got the dimension right
    assert ( bases_T.RolyComplk->dimension() == PolynomialSpaceDimension<Cell>::RolyCompl(m_K) );
  } // if

  return bases_T;
}

//------------------------------------------------------------------------------

DDRCore::FaceBases DDRCore::_construct_face_bases(size_t iF)
{
  const Face & F = *m_mesh.face(iF);
  
  FaceBases bases_F;

  //------------------------------------------------------------------------------
  // Basis for Pk+1(F)
  //------------------------------------------------------------------------------
  
  MonomialScalarBasisFace basis_Pkpo_F(F, m_K + 1);
  QuadratureRule quad_2kpo_F = generate_quadrature_rule(F, 2 * (m_K + 1));
  auto basis_Pkpo_F_quad = evaluate_quad<Function>::compute(basis_Pkpo_F, quad_2kpo_F);
  // Orthonormalize and store the basis
  bases_F.Polykpo.reset( new PolyBasisFaceType(l2_orthonormalize(basis_Pkpo_F, quad_2kpo_F, basis_Pkpo_F_quad)) );
  // Check that we got the dimension right
  assert( bases_F.Polykpo->dimension() == PolynomialSpaceDimension<Face>::Poly(m_K + 1) );

  //------------------------------------------------------------------------------
  // Basis for Pk(F) and Pk-1(F)
  //------------------------------------------------------------------------------

  // Given that the basis for Pk+1(F) is hierarchical, bases for Pk(F) and
  // Pk-1(F) can be obtained by restricting the former
  bases_F.Polyk.reset( new RestrictedBasis<PolyBasisFaceType>(*bases_F.Polykpo, PolynomialSpaceDimension<Face>::Poly(m_K)) );
  if (PolynomialSpaceDimension<Face>::Poly(m_K - 1) > 0) {
    bases_F.Polykmo.reset( new RestrictedBasis<PolyBasisFaceType>(*bases_F.Polykpo, PolynomialSpaceDimension<Face>::Poly(m_K - 1)) );
  }
  //------------------------------------------------------------------------------
  // Basis for Pk(F)^2
  //------------------------------------------------------------------------------

  MonomialScalarBasisFace basis_Pk_F(F, m_K);    
  QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * m_K);
  auto basis_Pk_F_quad = evaluate_quad<Function>::compute(basis_Pk_F, quad_2k_F);
  PolyBasisFaceType on_basis_Pk_F = l2_orthonormalize(basis_Pk_F, quad_2k_F, basis_Pk_F_quad);
  // Vectorize and store the basis
  bases_F.Polyk2.reset( new Poly2BasisFaceType(on_basis_Pk_F, basis_Pk_F.jacobian()) );
  // Check that we got the dimension right
  assert( bases_F.Polyk2->dimension() == 2 * PolynomialSpaceDimension<Face>::Poly(m_K) );
  
  //------------------------------------------------------------------------------
  // Basis for Rk-1(F)
  //------------------------------------------------------------------------------

  if (PolynomialSpaceDimension<Face>::Roly(m_K - 1) > 0) {
    bases_F.Rolykmo.reset( new RolyBasisFaceType(ShiftedBasis<PolyBasisFaceType>(on_basis_Pk_F, 1)) );
    // Check that we got the dimension right
    assert( bases_F.Rolykmo->dimension() == PolynomialSpaceDimension<Face>::Roly(m_K - 1) );
  }
  
  //------------------------------------------------------------------------------
  // Basis fore ROk(F)
  //------------------------------------------------------------------------------

  if (PolynomialSpaceDimension<Face>::RolyCompl(m_K) > 0) {
    auto basis_Pk2_F_quad = evaluate_quad<Function>::compute(*bases_F.Polyk2, quad_2k_F);
    RolyBasisFaceType basis_Rk_F(ShiftedBasis<PolyBasisFaceType>(*bases_F.Polykpo, 1));
    auto basis_Rk_F_quad = evaluate_quad<Function>::compute(basis_Rk_F, quad_2k_F);
    Eigen::MatrixXd MR = compute_gram_matrix(basis_Rk_F_quad, basis_Pk2_F_quad, quad_2k_F);
    bases_F.RolyComplk.reset( new RolyComplBasisFaceType(*bases_F.Polyk2, MR.fullPivLu().kernel().transpose()) );
    // Check that we got the dimension right
    assert ( bases_F.RolyComplk->dimension() == PolynomialSpaceDimension<Face>::RolyCompl(m_K) );
  }

  //------------------------------------------------------------------------------
  // Basis for ROk+2(F)
  //------------------------------------------------------------------------------

  CurlBasis<ShiftedBasis<MonomialScalarBasisFace> >
    spanning_set_Rkp2_F(ShiftedBasis<MonomialScalarBasisFace>(MonomialScalarBasisFace(F, m_K + 3), 1));
  TangentFamily<MonomialScalarBasisFace> basis_Pkp2_2_F(MonomialScalarBasisFace(F, m_K + 2), basis_Pk_F.jacobian());
  QuadratureRule quad_2kp2_F = generate_quadrature_rule(F, 2 * (m_K + 2));
  auto spanning_set_Rkp2_F_quad = evaluate_quad<Function>::compute(spanning_set_Rkp2_F, quad_2kp2_F);
  auto basis_Pkp2_2_F_quad = evaluate_quad<Function>::compute(basis_Pkp2_2_F, quad_2kp2_F);
  Eigen::MatrixXd MRp2 = compute_gram_matrix(spanning_set_Rkp2_F_quad, basis_Pkp2_2_F_quad, quad_2kp2_F);
  Family<TangentFamily<MonomialScalarBasisFace> >
    basis_ROkp2_F(basis_Pkp2_2_F, MRp2.fullPivLu().kernel().transpose());
  auto basis_ROkp2_F_quad = evaluate_quad<Function>::compute(basis_ROkp2_F, quad_2kp2_F);
  bases_F.RolyComplkp2.reset( new PotentialTestFunctionBasisType(l2_orthonormalize(basis_ROkp2_F, quad_2kp2_F, basis_ROkp2_F_quad)) );
  // Check that we got the dimension right
  assert( bases_F.RolyComplkp2->dimension() == PolynomialSpaceDimension<Face>::RolyCompl(m_K + 2) );
  
  return bases_F;
}

//------------------------------------------------------------------------------

DDRCore::EdgeBases DDRCore::_construct_edge_bases(size_t iE)
{
  const Edge & E = *m_mesh.edge(iE);

  EdgeBases bases_E;

  // Basis for Pk+1(E)
  MonomialScalarBasisEdge basis_Pkpo_E(E, m_K + 1);
  QuadratureRule quad_2kpo_E = generate_quadrature_rule(E, 2 * (m_K + 1));
  auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(basis_Pkpo_E, quad_2kpo_E);
  bases_E.Polykpo.reset( new PolyEdgeBasisType(l2_orthonormalize(basis_Pkpo_E, quad_2kpo_E, basis_Pkpo_E_quad)) );

  // Basis for Pk(E)
  bases_E.Polyk.reset( new RestrictedBasis<PolyEdgeBasisType>(*bases_E.Polykpo, PolynomialSpaceDimension<Edge>::Poly(m_K)) );
  
  // Basis for Pk-1(E)
  if (PolynomialSpaceDimension<Edge>::Poly(m_K - 1) > 0) {
    // Given that the basis for Pk+1(E) is hierarchical, a basis for Pk-1(E)
    // can be obtained by restricting the former
    bases_E.Polykmo.reset( new RestrictedBasis<PolyEdgeBasisType>(*bases_E.Polykpo, PolynomialSpaceDimension<Edge>::Poly(m_K - 1)) );
  }

  return bases_E;
}
