#ifndef LASXDIV_HPP
#define LASXDIV_HPP

#include <laxdiv.hpp>
#include <sxdiv.hpp>
#include <sxcurl.hpp>
#include <unsupported/Eigen/KroneckerProduct>

namespace HArDCore3D
{
  /*!
   *	\addtogroup LADDRCore
   * @{
   */

  /// Discrete Serendipity Hdiv space: local operators, L2 product and global interpolator
  /** This is actually just an XDiv space with additional function to compute the L2 product of discrete curl of elements in SXcurl. */

  class LASXDiv : public LAXDiv
  {
  public:

    /// Constructor
    LASXDiv(const LieAlgebra & lie_algebra, const XDiv & xdiv, const SXDiv & sxdiv, bool use_threads = true, std::ostream & output = std::cout);

    /// Return the Lie algebra
    inline const LieAlgebra & lieAlg() const
    {
      return m_lie_algebra;
    }

    /// Compute the matrix of the (weighted) L2-product as 'computeL2Product', with application of the discrete serendipity curl on the left/right/both sides (depending on argument "side").
    Eigen::MatrixXd computeL2ProductCurl(
                                     const size_t iT, ///< index of the cell
                                     const SXCurl & sx_curl, ///< instance of SXCurl to access the serendipity curls
                                     const std::string & side, ///< which side (left, right, both) we apply the gradient to
                                     const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                     const Eigen::MatrixXd & mass_Pk3_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of (P^k(T))^3; if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to constant 1.
                                     ) const;
  private:
    const SXDiv & m_sxdiv;

  };
} // end of namespace HArDCore3D

#endif
