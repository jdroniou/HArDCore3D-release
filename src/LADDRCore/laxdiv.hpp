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

  /// Discrete Lie algebra valued Hdiv space: local operators, L2 product and global interpolator
  /** Each DOF in the XDiv space corresponds to the dimension of the Lie algebra number of DOFs in the LAXDiv space. These DOFs are stored consecutively with overall ordering given by XDiv; they represent and are locally ordered by the basis in the LieAlgebra */


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

    /// Return the mesh
    const Mesh & mesh() const
    {
      return m_xdiv.mesh();
    }
    
    /// Return the polynomial degree
    const size_t & degree() const
    {
      return m_xdiv.degree();
    }

    /// Interpolator of a continuous function
    Eigen::VectorXd interpolate(
          const LAFunctionType & v, ///< The Lie algebra valued function to interpolate
          const int doe_cell = -1, ///< The optional degre of cell quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_face = -1 ///< The optional degre of face quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
                ) const;
    
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

    /// Return cell bases for the face of index iT
    inline const DDRCore::CellBases & cellBases(size_t iT) const
    {
      return m_xdiv.cellBases(iT);
    }

    /// Return cell bases for cell T
    inline const DDRCore::CellBases & cellBases(const Cell & T) const
    {
      return m_xdiv.cellBases(T.global_index());
    }
    
    /// Return face bases for the face of index iF
    inline const DDRCore::FaceBases & faceBases(size_t iF) const
    {
      return m_xdiv.faceBases(iF);
    }

    /// Return cell bases for face F
    inline const DDRCore::FaceBases & faceBases(const Face & F) const
    {
      return m_xdiv.faceBases(F.global_index());
    }
    
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

  protected:
    LocalOperators _compute_cell_divergence_potential(size_t iT);

    const LieAlgebra & m_lie_algebra;
    const XDiv & m_xdiv;
    bool m_use_threads;
    std::ostream & m_output;

    // Containers for local operators
    std::vector<std::unique_ptr<LocalOperators> > m_cell_operators;
  };
} // end of namespace HArDCore3D

#endif
