#include <cassert>

#include <vemcore.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_face.hpp>
#include <GMpoly_edge.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------

VEMCore::VEMCore(const Mesh & mesh, size_t K, bool use_threads, std::ostream & output)
  : m_mesh(mesh),
    m_K(K),
    m_output(output),
    m_cell_bases(mesh.n_cells()),
    m_face_bases(mesh.n_faces()),
    m_edge_bases(mesh.n_edges())
{
  m_output << "[VEMCore] Initializing" << std::endl;
  
  // Construct element bases
  std::function<void(size_t, size_t)> construct_all_cell_bases
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iT = start; iT < end; iT++) {
	        this->m_cell_bases[iT].reset( new CellBases(this->_construct_cell_bases(iT)) );
	      } // for iT
      };

  m_output << "[VEMCore] Constructing element bases" << std::endl;
  parallel_for(mesh.n_cells(), construct_all_cell_bases, use_threads);
  
  // Construct face bases
  std::function<void(size_t, size_t)> construct_all_face_bases
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iF = start; iF < end; iF++) {
	        this->m_face_bases[iF].reset( new FaceBases(_construct_face_bases(iF)) );
	      } // for iF
      };
  
  m_output << "[VEMCore] Constructing face bases" << std::endl;
  parallel_for(mesh.n_faces(), construct_all_face_bases, use_threads);

  // Construct edge bases
  std::function<void(size_t, size_t)> construct_all_edge_bases   
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iE = start; iE < end; iE++) {
	        this->m_edge_bases[iE].reset( new EdgeBases(_construct_edge_bases(iE)) );
	      } // for iF
      };
  
  m_output << "[VEMCore] Constructing edge bases" << std::endl;
  parallel_for(mesh.n_edges(), construct_all_edge_bases, use_threads);
}

//------------------------------------------------------------------------------

VEMCore::CellBases VEMCore::_construct_cell_bases(size_t iT)
{
  const Cell & T = *m_mesh.cell(iT);

  CellBases bases_T;
  
  MonomialCellIntegralsType int_monoT_2kp4 = IntegrateCellMonomials(T, 2*(m_K+2));
  
  //------------------------------------------------------------------------------
  // Basis for Pk+1(T), Pk+1_0(T), Pk(T), Pk-1(T), Pk0(T) and Pk(T)^3
  //------------------------------------------------------------------------------
  
  MonomialScalarBasisCell basis_Pkpo_T(T, m_K + 1);
  bases_T.Polykpo.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pkpo_T, GramMatrix(T, basis_Pkpo_T, int_monoT_2kp4))) );  

  bases_T.Polykpo0.reset( new ShiftedBasis<PolyBasisCellType>(*bases_T.Polykpo, 1) );

  MonomialScalarBasisCell basis_Pk_T(T, m_K);
  bases_T.Polyk.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pk_T, GramMatrix(T, basis_Pk_T, int_monoT_2kp4))) );  

  // Check that we got the dimensions right
  assert( bases_T.Polykpo->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K + 1) );
  assert( bases_T.Polyk->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K) );

  if (PolynomialSpaceDimension<Cell>::Poly(m_K - 1) > 0) {
    MonomialScalarBasisCell basis_Pkmo_T(T, m_K-1);
    bases_T.Polykmo.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pkmo_T, GramMatrix(T, basis_Pkmo_T, int_monoT_2kp4))) );  
    assert( bases_T.Polykmo->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K-1) );

    bases_T.Polyk0.reset( new ShiftedBasis<PolyBasisCellType>(*bases_T.Polyk, 1) );
  }  

  bases_T.Polyk3.reset( new Poly3BasisCellType(*bases_T.Polyk) );
  assert( bases_T.Polyk3->dimension() == 3 * PolynomialSpaceDimension<Cell>::Poly(m_K) );
    
  //------------------------------------------------------------------------------
  // Basis for Gk-1(T)
  //------------------------------------------------------------------------------

  if (PolynomialSpaceDimension<Cell>::Goly(m_K - 1) > 0) {
    GradientBasis<ShiftedBasis<MonomialScalarBasisCell> >
      basis_Gkmo_T(ShiftedBasis<MonomialScalarBasisCell>(MonomialScalarBasisCell(T, m_K), 1));

    // Orthonormalize and store the basis
    bases_T.Golykmo.reset( new GolyBasisCellType(l2_orthonormalize(basis_Gkmo_T, GramMatrix(T, basis_Gkmo_T, int_monoT_2kp4))) );

    // Check that we got the dimension right
    assert( bases_T.Golykmo->dimension() == PolynomialSpaceDimension<Cell>::Goly(m_K - 1) );
  } // if

  //------------------------------------------------------------------------------
  // Bases for Gck(T), Gck+1(T), and Rk-1(T)
  //------------------------------------------------------------------------------
 
  // Gck+1(T) (orthonormalised)
  GolyComplBasisCell basis_Gckpo_T(T, m_K+1);
  bases_T.GolyComplkpo.reset( new GolyComplBasisCellType(l2_orthonormalize(basis_Gckpo_T, GramMatrix(T, basis_Gckpo_T, int_monoT_2kp4))) );

  // check dimension
  assert( bases_T.GolyComplkpo->dimension() == PolynomialSpaceDimension<Cell>::GolyCompl(m_K + 1) );
 

  if (PolynomialSpaceDimension<Cell>::GolyCompl(m_K) > 0) {
    // Gck(T)
    GolyComplBasisCell basis_Gck_T(T, m_K);
    bases_T.GolyComplk.reset( new GolyComplBasisCellType(l2_orthonormalize(basis_Gck_T, GramMatrix(T, basis_Gck_T, int_monoT_2kp4))) );
    assert( bases_T.GolyComplk->dimension() == PolynomialSpaceDimension<Cell>::GolyCompl(m_K) );

    // Basis for curl Gck. We do not want to restart from bases_T.GolyComplk because it is orthonormalised (so a 
    // Family of other bases); if we started from this one, after orthonormalisation, the basis of Rk-1(T) would be
    // a Family of a Family, for which any evaluation could be quite expensive. 
    CurlBasis<GolyComplBasisCell> basis_curl_Gck_T(basis_Gck_T);
    bases_T.Rolykmo.reset( new RolyBasisCellType(l2_orthonormalize(basis_curl_Gck_T, GramMatrix(T, basis_curl_Gck_T, int_monoT_2kp4))) );   
    assert( bases_T.Rolykmo->dimension() == PolynomialSpaceDimension<Cell>::Roly(m_K - 1));
  } // if


  //------------------------------------------------------------------------------
  // Basis for Rck(T) and Rck+2(T)
  //------------------------------------------------------------------------------
  // Rck+2(T) (orthonormalised)
  RolyComplBasisCell basis_Rckp2_T(T, m_K+2);
  bases_T.RolyComplkp2.reset( new RolyComplBasisCellType(l2_orthonormalize(basis_Rckp2_T, GramMatrix(T, basis_Rckp2_T, int_monoT_2kp4))) );
  assert ( bases_T.RolyComplkp2->dimension() == PolynomialSpaceDimension<Cell>::RolyCompl(m_K+2) );

  // Rck(T) (orthonormalised). Could probably also be obtained as a RestrictedBasis of the previous one, but would
  // need to check if the basis for Rck are indeed hierarchical
  if (PolynomialSpaceDimension<Cell>::RolyCompl(m_K) > 0) { 
    RolyComplBasisCell basis_Rck_T(T, m_K);
    bases_T.RolyComplk.reset( new RolyComplBasisCellType(l2_orthonormalize(basis_Rck_T, GramMatrix(T, basis_Rck_T, int_monoT_2kp4))) );
    assert ( bases_T.RolyComplk->dimension() == PolynomialSpaceDimension<Cell>::RolyCompl(m_K) );
  } // if

  return bases_T;
}

//------------------------------------------------------------------------------

VEMCore::FaceBases VEMCore::_construct_face_bases(size_t iF)
{
  const Face & F = *m_mesh.face(iF);
  
  FaceBases bases_F;

  MonomialFaceIntegralsType int_monoF_2kp4 = IntegrateFaceMonomials(F, 2*(m_K+2));

  //------------------------------------------------------------------------------
  // Basis for Pk+2(F), Pk+1(F), Pk(F), Pk-1(F)
  //------------------------------------------------------------------------------
  MonomialScalarBasisFace basis_Pkp2_F(F, m_K + 2);
  bases_F.Polykp2.reset( new PolyBasisFaceType(l2_orthonormalize(basis_Pkp2_F, GramMatrix(F, basis_Pkp2_F, int_monoF_2kp4))) );

  MonomialScalarBasisFace basis_Pkpo_F(F, m_K + 1);
  bases_F.Polykpo.reset( new PolyBasisFaceType(l2_orthonormalize(basis_Pkpo_F, GramMatrix(F, basis_Pkpo_F, int_monoF_2kp4))) );
  
  MonomialScalarBasisFace basis_Pk_F(F, m_K);
  bases_F.Polyk.reset( new PolyBasisFaceType(l2_orthonormalize(basis_Pk_F, GramMatrix(F, basis_Pk_F, int_monoF_2kp4))) );

  // Check that we got the dimensions right
  assert( bases_F.Polykp2->dimension() == PolynomialSpaceDimension<Face>::Poly(m_K + 2) );
  assert( bases_F.Polykpo->dimension() == PolynomialSpaceDimension<Face>::Poly(m_K + 1) );
  assert( bases_F.Polyk->dimension() == PolynomialSpaceDimension<Face>::Poly(m_K) );

  if (PolynomialSpaceDimension<Face>::Poly(m_K - 1) > 0) {
    MonomialScalarBasisFace basis_Pkmo_F(F, m_K-1);
    bases_F.Polykmo.reset( new PolyBasisFaceType(l2_orthonormalize(basis_Pkmo_F, GramMatrix(F, basis_Pkmo_F, int_monoF_2kp4))) );
    assert( bases_F.Polykmo->dimension() == PolynomialSpaceDimension<Face>::Poly(m_K-1) );

    bases_F.Polyk0.reset( new ShiftedBasis<PolyBasisFaceType>(*bases_F.Polyk, 1) );
  }

  //------------------------------------------------------------------------------
  // Basis for Pk(F)^2, Pk+1(F)^2
  //------------------------------------------------------------------------------

  // Bases of Pk(F)^2 and Pk+1(F)^2 as TangentFamily. We use the system of coordinates of the basis on the face as generators of the face
  bases_F.Polyk2.reset( new Poly2BasisFaceType(*bases_F.Polyk, basis_Pkpo_F.coordinates_system()) );
  bases_F.Polykpo2.reset( new TangentFamily<PolyBasisFaceType>(*bases_F.Polykpo, basis_Pkpo_F.coordinates_system()) );
  // Check dimension
  assert( bases_F.Polyk2->dimension() == 2 * PolynomialSpaceDimension<Face>::Poly(m_K) );
  assert( bases_F.Polykpo2->dimension() == 2 * PolynomialSpaceDimension<Face>::Poly(m_K+1) );
  
  
  //------------------------------------------------------------------------------
  // Basis for Rk-1(F)
  //------------------------------------------------------------------------------

  if (PolynomialSpaceDimension<Face>::Roly(m_K - 1) > 0) {
    // Non-orthonormalised basis of Rk-1(F). 
    MonomialScalarBasisFace basis_Pk_F(F, m_K);
    ShiftedBasis<MonomialScalarBasisFace> basis_Pk0_F(basis_Pk_F,1);
    CurlBasis<ShiftedBasis<MonomialScalarBasisFace>> basis_Rkmo_F(basis_Pk0_F);
    // Orthonormalise, store and check dimension
    bases_F.Rolykmo.reset( new RolyBasisFaceType(l2_orthonormalize(basis_Rkmo_F, GramMatrix(F, basis_Rkmo_F, int_monoF_2kp4))) );
    assert( bases_F.Rolykmo->dimension() == PolynomialSpaceDimension<Face>::Roly(m_K - 1) );
  }
  
  //------------------------------------------------------------------------------
  // Basis for Rck(F)
  //------------------------------------------------------------------------------

  if (PolynomialSpaceDimension<Face>::RolyCompl(m_K) > 0) {
    RolyComplBasisFace basis_Rck_F(F, m_K);
    bases_F.RolyComplk.reset( new RolyComplBasisFaceType(l2_orthonormalize(basis_Rck_F, GramMatrix(F, basis_Rck_F, int_monoF_2kp4))) );
    assert ( bases_F.RolyComplk->dimension() == PolynomialSpaceDimension<Face>::RolyCompl(m_K) );
  }

  //------------------------------------------------------------------------------
  // Basis for Rck-1(F)
  //------------------------------------------------------------------------------

  if (PolynomialSpaceDimension<Face>::RolyCompl(m_K-1) > 0) {
    RolyComplBasisFace basis_Rckmo_F(F, m_K-1);
    bases_F.RolyComplkmo.reset( new RolyComplBasisFaceType(l2_orthonormalize(basis_Rckmo_F, GramMatrix(F, basis_Rckmo_F, int_monoF_2kp4))) );
    assert ( bases_F.RolyComplkmo->dimension() == PolynomialSpaceDimension<Face>::RolyCompl(m_K-1) );
  }
  
  //------------------------------------------------------------------------------
  // Basis for Rck+1(F)
  //------------------------------------------------------------------------------

  RolyComplBasisFace basis_Rckpo_F(F, m_K+1);
  bases_F.RolyComplkpo.reset( new RolyComplBasisFaceType(l2_orthonormalize(basis_Rckpo_F, GramMatrix(F, basis_Rckpo_F, int_monoF_2kp4))) );
  assert ( bases_F.RolyComplkpo->dimension() == PolynomialSpaceDimension<Face>::RolyCompl(m_K+1) );

  return bases_F;
}

//------------------------------------------------------------------------------

VEMCore::EdgeBases VEMCore::_construct_edge_bases(size_t iE)
{
  const Edge & E = *m_mesh.edge(iE);

  EdgeBases bases_E;

  MonomialEdgeIntegralsType int_monoE_2kp2 = IntegrateEdgeMonomials(E, 2*degree()+2);

  // Basis for Pk+1(E)
  MonomialScalarBasisEdge basis_Pkpo_E(E, m_K + 1);
  bases_E.Polykpo.reset( new PolyBasisEdgeType(l2_orthonormalize(basis_Pkpo_E, GramMatrix(E, basis_Pkpo_E, int_monoE_2kp2))) );

  // Basis for Pk(E)
  MonomialScalarBasisEdge basis_Pk_E(E, m_K);
  bases_E.Polyk.reset( new PolyBasisEdgeType(l2_orthonormalize(basis_Pk_E, GramMatrix(E, basis_Pk_E, int_monoE_2kp2))) );

  // Basis for Pk-1(E)
  if (PolynomialSpaceDimension<Edge>::Poly(m_K - 1) > 0) {
    MonomialScalarBasisEdge basis_Pkmo_E(E, m_K - 1);
    bases_E.Polykmo.reset( new PolyBasisEdgeType(l2_orthonormalize(basis_Pkmo_E, GramMatrix(E, basis_Pkmo_E, int_monoE_2kp2))) );
  }

  return bases_E;
}
