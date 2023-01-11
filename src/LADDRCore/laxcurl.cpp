#include <laxcurl.hpp>
#include <parallel_for.hpp>
#include <unsupported/Eigen/KroneckerProduct>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

LAXCurl::LAXCurl(const LieAlgebra & lie_algebra, const XCurl & xcurl, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(
	     xcurl.mesh(),
	     0,
	     PolynomialSpaceDimension<Edge>::Poly(xcurl.degree()) * lie_algebra.dimension(),
	     PolynomialSpaceDimension<Face>::Roly(xcurl.degree() - 1) * lie_algebra.dimension() + PolynomialSpaceDimension<Face>::RolyCompl(xcurl.degree()) * lie_algebra.dimension(),
	     PolynomialSpaceDimension<Cell>::Roly(xcurl.degree() - 1) * lie_algebra.dimension() + PolynomialSpaceDimension<Cell>::RolyCompl(xcurl.degree())	* lie_algebra.dimension()     
	     ),
    m_xcurl(xcurl),
    m_lie_algebra(lie_algebra),
    m_use_threads(use_threads),
    m_output(output),
    m_cell_operators(xcurl.mesh().n_cells()),
    m_face_operators(xcurl.mesh().n_faces())
{

  output << "[LAXCurl] Initializing" << std::endl;
  if (m_use_threads) {
    m_output << "[LAXCurl] Parallel execution" << std::endl;
  } else {
    m_output << "[LAXCurl] Sequential execution" << std::endl;
  }

  // Construct face curls and potentials
  std::function<void(size_t, size_t)> construct_all_face_curls_potentials
   = [this](size_t start, size_t end)->void
     {
       for (size_t iF = start; iF < end; iF++) {
         m_face_operators[iF].reset( new LocalOperators(_compute_face_curl_potential(iF)) );
       } // for iF
     };

  m_output << "[LAXCurl] Constructing face curls and potentials" << std::endl;
  parallel_for(mesh().n_faces(), construct_all_face_curls_potentials, m_use_threads);

  // Construct cell curls and potentials
  std::function<void(size_t, size_t)> construct_all_cell_curls_potentials
   = [this](size_t start, size_t end)->void
     {
       for (size_t iT = start; iT < end; iT++) {
         m_cell_operators[iT].reset( new LocalOperators(_compute_cell_curl_potential(iT)) );
       } // for iT
     };

  m_output << "[LAXCurl] Constructing cell curls and potentials" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_curls_potentials, m_use_threads);  
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd LAXCurl::interpolate(
                                    const LAFunctionType & v, 
                                    const int doe_cell, 
                                    const int doe_face, 
                                    const int doe_edge) const
{
  // The matrix where the rows are storing the interpolate of each function of v in XCurl. This lines up in the columns the DOFs in LAXCurl that are the Lie algebra components of the same polynomial DOF in XCurl. 
  Eigen::MatrixXd v_M = Eigen::MatrixXd(m_lie_algebra.dimension(), m_xcurl.dimension());

  for (size_t i = 0; i < m_lie_algebra.dimension(); i++){
    v_M.row(i) = m_xcurl.interpolate(v[i], doe_cell, doe_face, doe_edge);
  }
  // To construct the element of LAXCurl, we read the matrix v_M column by column.
  Eigen::VectorXd vh(Eigen::Map<Eigen::VectorXd>(v_M.data(), v_M.size()));

  return vh;

}

//------------------------------------------------------------------------------
// Curl and potential reconstruction
//------------------------------------------------------------------------------

 LAXCurl::LocalOperators LAXCurl::_compute_face_curl_potential(size_t iF){

  Eigen::MatrixXd F_curl = m_xcurl.faceOperators(iF).curl;
  Eigen::MatrixXd F_potential = m_xcurl.faceOperators(iF).potential;

  Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());

  return LocalOperators(Eigen::KroneckerProduct(F_curl, id), Eigen::KroneckerProduct(F_potential, id));
 }

 LAXCurl::LocalOperators LAXCurl::_compute_cell_curl_potential(size_t iT){
  
  Eigen::MatrixXd T_curl = m_xcurl.cellOperators(iT).curl;
  Eigen::MatrixXd T_potential = m_xcurl.cellOperators(iT).potential;

  Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());

  return LocalOperators(Eigen::KroneckerProduct(T_curl, id), Eigen::KroneckerProduct(T_potential, id));
 }

//------------------------------------------------------------------------------
// Local L2 products on LAXCurl
//------------------------------------------------------------------------------

Eigen::MatrixXd LAXCurl::computeL2Product(
                                        const size_t iT,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pk3_T,
                                        const IntegralWeight & weight
                                        ) const
{
  
  Eigen::MatrixXd L2_xcurl = m_xcurl.computeL2Product(iT, penalty_factor, mass_Pk3_T, weight);

  return Eigen::KroneckerProduct(L2_xcurl, m_lie_algebra.massMatrix());

}

Eigen::MatrixXd LAXCurl::computeL2ProductGradient(
                                                  const size_t iT, 
                                                  const XGrad & x_grad, 
                                                  const std::string & side, 
                                                  const double & penalty_factor,
                                                  const Eigen::MatrixXd & mass_Pk3_T, 
                                                  const IntegralWeight & weight
                                                  ) const
{

  Eigen::MatrixXd L2_xcurl_grad = m_xcurl.computeL2ProductGradient(iT, x_grad, side, penalty_factor, mass_Pk3_T, weight);

  return Eigen::KroneckerProduct(L2_xcurl_grad, m_lie_algebra.massMatrix());

}
