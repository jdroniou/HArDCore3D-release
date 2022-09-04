#ifndef LAXDIV_HPP
#define LAXDIV_HPP

#include <globaldofspace.hpp>
#include <xdiv.hpp>
#include <liealgebra.hpp>

namespace HArDCore3D
{
  /*!
   *	\addtogroup LADDRCore
   * @{
   */

  /// Discrete Hdiv space: local operators, L2 product and global interpolator
  /** On each face, the DOFs correspond to the polynomial bases on the face provided by m_ddr_core. On each element, the DOFs are first those of the \f$\mathcal{G}^{k-1}\f$ component and then of the \f$\mathcal{G}^{c,k}\f$ component, each one of them following the bases of these spaces provided by m_ddr_core. */

  class LAXDiv : public GlobalDOFSpace
  {
  public:

    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> FunctionType;
    typedef std::vector<FunctionType> LAFunctionType;

    /// A structure to store local operators (divergence and potential)
    struct LocalOperators
    {
      LocalOperators(
		     const Eigen::MatrixXd & _divergence,     ///< Divergence operator
		     const Eigen::MatrixXd & _divergence_rhs, ///< RHS for the divergence operator problem
		     const Eigen::MatrixXd & _potential       ///< Potential operator
		     )
	    : divergence(_divergence),
	      divergence_rhs(_divergence_rhs),
	      potential(_potential)
      {
	      // Do nothing
      }

      Eigen::MatrixXd divergence;
      Eigen::MatrixXd divergence_rhs;
      Eigen::MatrixXd potential;
    };

    /// Constructor
    LAXDiv(const LieAlgebra & lie_algebra, const XDiv & xdiv, bool use_threads = true, std::ostream & output = std::cout);
    
    /// Interpolator of a continuous function
    Eigen::VectorXd interpolate(
          const LAFunctionType & v, ///< The Lie algebra valued function to interpolate
          const int doe_cell = -1, ///< The optional degre of cell quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_face = -1 ///< The optional degre of face quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
                ) const;
    
    /// Compute the matrix of the (weighted) L2-product for the cell of index iT.
    // The mass matrix of P^k(T)^3 is the most expensive mass matrix in the calculation of this norm, which
    // is why there's the option of passing it as parameter if it's been already pre-computed when the norm is called.
    Eigen::MatrixXd computeL2Product(
                                     const size_t iT, ///< index of the cell
                                     const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                     const Eigen::MatrixXd & mass_Pk3_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of (P^k(T))^3; if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to 1
                                     ) const;
                                     
    /// Compute the matrix of the (weighted) L2-product as 'computeL2Product', with application of the discrete curl on the left/right/both sides (depending on argument "side").
    Eigen::MatrixXd computeL2ProductCurl(
                                     const size_t iT, ///< index of the cell
                                     const XCurl & x_curl, ///< instance of XCurl to access the full curls
                                     const std::string & side, ///< which side (left, right, both) we apply the gradient to
                                     const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                     const Eigen::MatrixXd & mass_Pk3_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of (P^k(T))^3; if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to constant 1.
                                     ) const;

  private:
    // LocalOperators _compute_cell_divergence_potential(size_t iT);

    const LieAlgebra & m_lie_algebra;
    const XDiv & m_xdiv;
    bool m_use_threads;
    std::ostream & m_output;

    // Containers for local operators
    // std::vector<std::unique_ptr<LocalOperators> > m_cell_operators;
  };
} // end of namespace HArDCore3D

#endif
