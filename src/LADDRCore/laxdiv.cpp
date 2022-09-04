#include <laxdiv.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>
#include <unsupported/Eigen/KroneckerProduct>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

LAXDiv::LAXDiv(const LieAlgebra & lie_algebra, const XDiv & xdiv, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(
	     xdiv.mesh(), 
       0, 
       0,
	     PolynomialSpaceDimension<Face>::Poly(xdiv.degree())*lie_algebra.dimension(),
	     PolynomialSpaceDimension<Cell>::Goly(xdiv.degree() - 1)*lie_algebra.dimension() + PolynomialSpaceDimension<Cell>::GolyCompl(xdiv.degree())*lie_algebra.dimension()
	    ),
    m_lie_algebra(lie_algebra),
    m_xdiv(xdiv),
    m_use_threads(use_threads),
    m_output(output)
    // m_cell_operators(ddr_core.mesh().n_cells())
{
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd LAXDiv::interpolate(const LAFunctionType & v,  const int doe_cell, const int doe_face) const
{

  Eigen::MatrixXd v_M(m_lie_algebra.dimension(), m_xdiv.dimension());

  for (size_t i = 0; i < m_lie_algebra.dimension(); i++){
    v_M.row(i) = m_xdiv.interpolate(v[i], doe_cell, doe_face);
  }

  Eigen::VectorXd vh(Eigen::Map<Eigen::VectorXd>(v_M.data(), v_M.size()));

  return vh;

}

//------------------------------------------------------------------------------
// local L2 products on LAXDiv
//------------------------------------------------------------------------------

Eigen::MatrixXd LAXDiv::computeL2Product(
                                        const size_t iT,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pk3_T,
                                        const IntegralWeight & weight
                                        ) const
{

  Eigen::MatrixXd L2_xdiv = m_xdiv.computeL2Product(iT, penalty_factor, mass_Pk3_T, weight);  

  return Eigen::KroneckerProduct(L2_xdiv, m_lie_algebra.massMatrix());

}

Eigen::MatrixXd LAXDiv::computeL2ProductCurl(
                                        const size_t iT,
                                        const XCurl & x_curl,
                                        const std::string & side,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pk3_T,
                                        const IntegralWeight & weight
                                        ) const
{

  Eigen::MatrixXd L2_xdiv_curl = m_xdiv.computeL2ProductCurl(iT, x_curl, side, penalty_factor, mass_Pk3_T, weight);

  return Eigen::KroneckerProduct(L2_xdiv_curl, m_lie_algebra.massMatrix());

}


