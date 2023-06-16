#include <lasxdiv.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

LASXDiv::LASXDiv(const LieAlgebra & lie_algebra, const XDiv & xdiv, const SXDiv & sxdiv, bool use_threads, std::ostream & output)
  : LAXDiv(lie_algebra, xdiv, use_threads, output),
  m_sxdiv(sxdiv)
{
  m_output << "[LASXDiv] Initializing" << std::endl;
}


// L2 product with Serendipity curl

Eigen::MatrixXd LASXDiv::computeL2ProductCurl(
                                        const size_t iT,
                                        const SXCurl & sx_curl,
                                        const std::string & side,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pk3_T,
                                        const IntegralWeight & weight
                                        ) const
{
 return Eigen::KroneckerProduct(m_sxdiv.computeL2ProductCurl(iT, sx_curl, side, penalty_factor, mass_Pk3_T, weight), m_lie_algebra.massMatrix());
}

