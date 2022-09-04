#ifndef LAXGRAD_HPP
#define LAXGRAD_HPP

#include <globaldofspace.hpp>
#include <xgrad.hpp>
#include <liealgebra.hpp>

namespace HArDCore3D
{
  /*!
   *  \addtogroup LADDRCore
   * @{
   */

  /// Discrete H1 space: local operators, L2 product and global interpolator
  /** On each edge/face/element, the DOFs (if any) correspond to the polynomial bases on the edge/face/element provided by m_ddr_core */
  class LAXGrad : public GlobalDOFSpace
  {
  public:

    typedef std::function<double(const Eigen::Vector3d &)> FunctionType;
    typedef std::vector<FunctionType> LAFunctionType;

    // /// A structure to store local operators (gradient and potential)
    // struct LocalOperators
    // {
    //   LocalOperators(
    //                  const Eigen::MatrixXd & _gradient, ///< Gradient operator
    //                  const Eigen::MatrixXd & _potential ///< Potential operator
    //                  )
    //     : gradient(_gradient),
    //       potential(_potential)
    //   {
    //     // Do nothing
    //   }
      
    //   Eigen::MatrixXd gradient;
    //   Eigen::MatrixXd potential;
    // };
    
    /// Constructor
    LAXGrad(const LieAlgebra & lie_algebra, const XGrad & xgrad, bool use_threads = true, std::ostream & output = std::cout);
    
    /// Interpolator of a continuous function
    Eigen::VectorXd interpolate(
          const LAFunctionType & q, ///< The function to interpolate
          const int doe_cell = -1, ///< The optional degre of cell quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_face = -1, ///< The optional degre of face quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_edge = -1 ///< The optional degre of edge quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          ) const;

    /// Compute the matrix of the (weighted) L2-product for the cell of index iT. The stabilisation here is based on cell and face potentials.
    // The mass matrix of P^{k+1}(T) is the most expensive mass matrix in the calculation of this norm, which
    // is why there's the option of passing it as parameter if it's been already pre-computed when the norm is called.
    Eigen::MatrixXd computeL2Product(
                                     const size_t iT, ///< index of the cell
                                     const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                     const Eigen::MatrixXd & mass_Pkpo_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of P^{k+1}(T); if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to 1
                                     ) const;

  private:   

    const LieAlgebra & m_lie_algebra;
    const XGrad & m_xgrad;
    bool m_use_threads;
    std::ostream & m_output;
    
  };

} // end of namespace HArDCore3D

#endif
