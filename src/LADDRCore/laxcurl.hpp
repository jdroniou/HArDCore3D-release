#ifndef LAXCURL_HPP
#define LAXCURL_HPP

#include <globaldofspace.hpp>
#include <xcurl.hpp>
#include <liealgebra.hpp>

namespace HArDCore3D
{
    /*!
   *  \addtogroup LADDRCore
   * @{
   */

  /// Discrete Lie algebra valued Hcurl space
    /** Each DOF in the XCurl space corresponds to the dimension of the Lie algebra number of DOFs in the LAXCurl space. These DOFs are stored consecutively with overall ordering given by XCurl; they represent and are locally ordered by the basis in the LieAlgebra */
  class LAXCurl : public GlobalDOFSpace
  {
  public:
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> FunctionType;
    typedef std::vector<FunctionType> LAFunctionType;

    /// Returns the Lie algebra
    inline const LieAlgebra & lieAlgebra() const
    {
      return m_lie_algebra;
    }
    /// Returns the space XCurl
    inline const XCurl & xCurl() const
    {
      return m_xcurl;
    }

     /// A structure to store the local operators (curl and potential)
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
    LAXCurl(const LieAlgebra & lie_algebra, const XCurl & xcurl, bool use_threads = true, std::ostream & output = std::cout);
    
    /// Return the mesh
    const Mesh & mesh() const
    {
      return m_xcurl.mesh();
    }
    
    /// Return the polynomial degree
    const size_t & degree() const
    {
      return m_xcurl.degree();
    }
    
    /// Interpolator of a continuous Lie algebra valued function decomposed on the basis of the LieAlgebra, given as a vector of functions. 
    Eigen::VectorXd interpolate(
                                const LAFunctionType & v, ///< The function to interpolate
                                const int doe_cell = -1, ///< The optional degre of cell quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
                                const int doe_face = -1, ///< The optional degre of face quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
                                const int doe_edge = -1 ///< The optional degre of edge quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
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
      return m_xcurl.cellBases(iT);
    }

    /// Return cell bases for cell T
    inline const DDRCore::CellBases & cellBases(const Cell & T) const
    {
      return m_xcurl.cellBases(T.global_index());
    }
    
    /// Return face bases for the face of index iF
    inline const DDRCore::FaceBases & faceBases(size_t iF) const
    {
      return m_xcurl.faceBases(iF);
    }

    /// Return cell bases for face F
    inline const DDRCore::FaceBases & faceBases(const Face & F) const
    {
      return m_xcurl.faceBases(F.global_index());
    }
    
    /// Return edge bases for the edge of index iE
    inline const DDRCore::EdgeBases & edgeBases(size_t iE) const
    {
      return m_xcurl.edgeBases(iE);
    }

    /// Return edge bases for edge E
    inline const DDRCore::EdgeBases & edgeBases(const Edge & E) const
    {
      return m_xcurl.edgeBases(E.global_index());
    }

    /// Compute the matrix of the (weighted) L2-product for the cell of index iT.
    /// We required the same arguments as the product functions in XCurl since they are called in the construction.
    /// From the choice of ordering of the DOFs, this matrix is the Kronecker product of the L2-product matrix of iT in XCurl with the mass matrix of the Lie algebra.
    Eigen::MatrixXd computeL2Product(
                                    const size_t iT, ///< index of the cell
                                    const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                    const Eigen::MatrixXd & mass_Pk3_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of (P^k(T))^3; if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                    const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to 1
                                    ) const;

    /// Compute the matrix of the (weighted) L2-product as 'computeL2Product', with application of the discrete gradient on the left/right/both sides (depending on argument "side").
    Eigen::MatrixXd computeL2ProductGradient(
                                      const size_t iT, ///< index of the cell
                                      const XGrad & x_grad, ///< instance of XGrad to access the full gradients
                                      const std::string & side, ///< which side (left, right, both) we apply the gradient to
                                      const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                      const Eigen::MatrixXd & mass_Pk3_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of (P^k(T))^3; if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                      const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to constant 1.
                                      ) const;

  private:
    LocalOperators _compute_face_curl_potential(size_t iF);
    LocalOperators _compute_cell_curl_potential(size_t iT);

    const XCurl & m_xcurl;
    const LieAlgebra & m_lie_algebra;
    bool m_use_threads;
    std::ostream & m_output;

    // Containers for local operators
    std::vector<std::unique_ptr<LocalOperators> > m_cell_operators;
    std::vector<std::unique_ptr<LocalOperators> > m_face_operators;

  };

}// end of namespace HArDCore3D

#endif
