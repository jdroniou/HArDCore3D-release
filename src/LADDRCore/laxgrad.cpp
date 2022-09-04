#include <laxgrad.hpp>
#include <unsupported/Eigen/KroneckerProduct>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

LAXGrad::LAXGrad(const LieAlgebra & lie_algebra, const XGrad & xgrad, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(xgrad.mesh(),
	     1*lie_algebra.dimension(),
	     PolynomialSpaceDimension<Edge>::Poly(xgrad.degree() - 1)*lie_algebra.dimension(),
	     PolynomialSpaceDimension<Face>::Poly(xgrad.degree() - 1)*lie_algebra.dimension(),
	     PolynomialSpaceDimension<Cell>::Poly(xgrad.degree() - 1)*lie_algebra.dimension()
	     ),
    m_lie_algebra(lie_algebra),
    m_xgrad(xgrad),
    m_use_threads(use_threads),
    m_output(output)
{
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd LAXGrad::interpolate(const LAFunctionType & q, const int doe_cell, const int doe_face, const int doe_edge) const
{

  Eigen::MatrixXd q_M(m_lie_algebra.dimension(), m_xgrad.dimension());

  for (size_t i = 0; i < m_lie_algebra.dimension(); i++){
    q_M.row(i) = m_xgrad.interpolate(q[i], doe_cell, doe_face, doe_edge);
  }

  Eigen::VectorXd qh(Eigen::Map<Eigen::VectorXd>(q_M.data(), q_M.size()));
  
  return qh;

}

//------------------------------------------------------------------------------
// local L2 products on LAXGrad
//------------------------------------------------------------------------------
Eigen::MatrixXd LAXGrad::computeL2Product(
                                        const size_t iT,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pkpo_T,
                                        const IntegralWeight & weight
                                        ) const
{

  Eigen::MatrixXd L2_xgrad = m_xgrad.computeL2Product(iT, penalty_factor, mass_Pkpo_T, weight);

  return Eigen::KroneckerProduct(L2_xgrad, m_lie_algebra.massMatrix());

}
