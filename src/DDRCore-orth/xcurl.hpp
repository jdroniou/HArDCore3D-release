#ifndef XCURL_HPP
#define XCURL_HPP

#include <ddrcore.hpp>
#include <ddrspace.hpp>
#include <integralweight.hpp>
#include <xgrad.hpp>

namespace HArDCore3D
{
  /*!
   *  \addtogroup DDRcore
   * @{
   */

  /// Discrete Hcurl space
  class XCurl : public DDRSpace
  {
  public:
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> FunctionType;

    /// A structure to store local operators
    struct LocalOperators
    {
      LocalOperators(
                     const Eigen::MatrixXd & _curl,     ///< Curl operator
                     const Eigen::MatrixXd & _potential ///< Potential operator
                     )
        : curl(_curl),
          potential(_potential)
      {
        // Do nothing
      }

      Eigen::MatrixXd curl;
      Eigen::MatrixXd potential;
    };

    /// Constructor
    XCurl(const DDRCore & ddr_core, bool use_threads = true, std::ostream & output = std::cout);

    /// Return the mesh
    const Mesh & mesh() const
    {
      return m_ddr_core.mesh();
    }
    
    /// Return the polynomial degree
    const size_t & degree() const
    {
      return m_ddr_core.degree();
    }

    /// Interpolator
    Eigen::VectorXd interpolate(const FunctionType & v) const;

    /// Return cell operators for the cell of index iT
    inline const LocalOperators & cellOperators(size_t iT) const
    {
      return *m_cell_operators[iT];
    }

    /// Return cell operators for cell T
    inline const LocalOperators & cellOperators(const Cell & T) const
    {
      return * m_cell_operators[T.global_index()];  
    }

    /// Return face operators for the face of index iF
    inline const LocalOperators & faceOperators(size_t iF) const
    {
      return *m_face_operators[iF];
    }

    /// Return face operators for face F
    inline const LocalOperators & faceOperators(const Face & F) const
    {
      return * m_face_operators[F.global_index()];  
    }

    /// Return cell bases for the face of index iT
    inline const DDRCore::CellBases & cellBases(size_t iT) const
    {
      return m_ddr_core.cellBases(iT);
    }

    /// Return cell bases for cell T
    inline const DDRCore::CellBases & cellBases(const Cell & T) const
    {
      return m_ddr_core.cellBases(T.global_index());
    }
    
    /// Return face bases for the face of index iF
    inline const DDRCore::FaceBases & faceBases(size_t iF) const
    {
      return m_ddr_core.faceBases(iF);
    }

    /// Return cell bases for face F
    inline const DDRCore::FaceBases & faceBases(const Face & F) const
    {
      return m_ddr_core.faceBases(F.global_index());
    }
    
    /// Return edge bases for the edge of index iE
    inline const DDRCore::EdgeBases & edgeBases(size_t iE) const
    {
      return m_ddr_core.edgeBases(iE);
    }

    /// Return edge bases for edge E
    inline const DDRCore::EdgeBases & edgeBases(const Edge & E) const
    {
      return m_ddr_core.edgeBases(E.global_index());
    }

    /// Compute the matrix of the (weighted) L2-product for the cell of index iT. The stabilisation here is based on cell and face potentials.
    // The mass matrix of P^k(T)^3 is the most expensive mass matrix in the calculation of this norm, which
    // is why there's the option of passing it as parameter if it's been already pre-computed when the norm is called.
    Eigen::MatrixXd computeL2Product (
                                     const size_t iT, ///< index of the cell
                                     const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                     const Eigen::MatrixXd & mass_Pk3_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of (P^k(T))^3; if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to 1
                                     ) const;
                                     
    /// Compute the matrix of the (weighted) L2-product as 'computeL2Product', with application of the discrete gradient on the left/right/both sides (depending on argument "side").
    Eigen::MatrixXd computeL2ProductGradient(
                                     const size_t iT, ///< index of the cell
                                     const XGrad & x_grad, ///< instance of XGrad to access the full gradients
                                     const std::string side, ///< which side (left, right, both) we apply the gradient to
                                     const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                     const Eigen::MatrixXd & mass_Pk3_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of (P^k(T))^3; if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to constant 1.
                                     ) const;

    /// Compute the matrix of the L2 product, applying leftOp and rightOp to the variables; this is the main function called by computeL2Product and computeL2ProductGradient
    Eigen::MatrixXd computeL2Product_with_Ops(
                                     const size_t iT, ///< index of the cell
                                     std::vector<Eigen::MatrixXd> & leftOp, ///< edge, face and element operators to apply on the left
                                     std::vector<Eigen::MatrixXd> & rightOp, ///< edge, face and element operators to apply on the right
                                     const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                     const Eigen::MatrixXd & mass_Pk3_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of (P^k(T))^3; if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to 1
                                     ) const;
                                

    /// Alternate version of L2 product, using only cell potentials and 'DOFs' (as in Di-Pietro & Droniou, JCP 2020)
    Eigen::MatrixXd computeL2ProductDOFs(
                                     size_t iT, ///< index of the cell
                                     const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                     const Eigen::MatrixXd & mass_Pk3_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of (P^k(T))^3; if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, default to constant equal to 1.
                                     );

////////    /// Compute L2-product for the cell of index iT
////////    Eigen::MatrixXd computeL2Product(size_t iT, const double & penalty_factor = 1.);

////////    /// Compute weighted L2-product for the cell of index iT
////////    Eigen::MatrixXd computeWeightedL2Product(
////////                                             size_t iT,
////////                                             const std::function<double(const Eigen::Vector3d &)> & weight,
////////                                             size_t qr_doe_increase = 0,
////////                                             const double & penalty_factor = 1.
////////                                             );

  private:
    LocalOperators _compute_face_curl_potential(size_t iF);
    LocalOperators _compute_cell_curl_potential(size_t iT);
    
    const DDRCore & m_ddr_core;
    bool m_use_threads;
    std::ostream & m_output;

    // Containers for local operators
    std::vector<std::unique_ptr<LocalOperators> > m_cell_operators;
    std::vector<std::unique_ptr<LocalOperators> > m_face_operators;
  };

} // end of namespace HArDCore3D
#endif
