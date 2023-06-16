#ifndef LASXCURL_HPP
#define LASXCURL_HPP

#include <sxcurl.hpp>
#include <liealgebra.hpp>
#include <unsupported/Eigen/KroneckerProduct>

namespace HArDCore3D
{
  /*!
   *  \addtogroup LADDRCore
   * @{
   */

  /// Discrete Lie algebra valued Serendipity Hcurl space: local operators, L2 product and global interpolator
    /** Each DOF in the XCurl space corresponds to the dimension of the Lie algebra number of DOFs in the LAXCurl space. These DOFs are stored consecutively with overall ordering given by XCurl; they represent and are locally ordered by the basis in the LieAlgebra */
  class LASXCurl : public VariableDOFSpace
  {
  public:
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> FunctionType;
    typedef std::vector<FunctionType> LAFunctionType;

    /// A structure to store the serendipity, extension and reduction operators
    /** The extension operators do SXCurl -> Xcurl. The reduction operators do not produce all unknowns, only the ones in the mesh entities they're associated to; so they do XCurl -> \f$R^{k-1}(P) \times R^{c,\ell_P+1}(P)\f$ where P is a face or an element (the only use of the reduction operator is to compute the interpolate, and this version is what is useful). */
    struct TransferOperators
    {
      TransferOperators(
                     const Eigen::MatrixXd & _serendipity,   ///< Serendipity reconstruction
                     const Eigen::MatrixXd & _extension,     ///< Extension SXCurl->XCurl
                     const Eigen::MatrixXd & _reduction      ///< Reduction XCurl->SXCurl
                     )
        : serendipity(_serendipity),
          extension(_extension),
          reduction(_reduction)
      {
        // Do nothing
      }

      Eigen::MatrixXd serendipity;
      Eigen::MatrixXd extension;
      Eigen::MatrixXd reduction;
    };

    /// Constructor
    LASXCurl(const LieAlgebra & lie_algebra, const SXCurl & sx_curl, bool use_threads = true, std::ostream & output = std::cout);

    /// Return the mesh
    const Mesh & mesh() const
    {
      return m_sxcurl.mesh();
    }
    
    /// Return the polynomial degree
    const size_t & degree() const
    {
      return m_sxcurl.degree();
    }

    /// Return the Lie algebra
    inline const LieAlgebra & lieAlg() const
    {
      return m_lie_algebra;
    }

    /// Return the serendipity operators
    const SerendipityProblem & serPro() const
    {
      return m_sxcurl.serPro();
    }

    /// Interpolator of a continuous function
    Eigen::VectorXd interpolate(
          const LAFunctionType & v, ///< The function to interpolate
          const int doe_cell = -1, ///< The optional degre of cell quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_face = -1, ///< The optional degre of face quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_edge = -1 ///< The optional degre of edge quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          ) const;

    //---------------------------------//
    //---- Lie algebra operators-------//
    
    inline const Eigen::VectorXd LaxcurlI(size_t iG, Eigen::VectorXd & v) const
    {
      return Eigen::Map<Eigen::VectorXd,0,Eigen::InnerStride<>>(v.data()+iG, m_sxcurl.dimension(),Eigen::InnerStride<>(m_lie_algebra.dimension()));
    }

    //---------------------------------//
    //---- Face transfer operators-----//
    
    /// Return the serendipity reconstruction for the face of index iF
    inline const Eigen::MatrixXd & ScurlFace(size_t iF) const
    {
      return (*m_face_transfer_operators[iF]).serendipity;
    }

    /// Return the extension for the face of index iF
    inline const Eigen::MatrixXd & EcurlFace(size_t iF) const
    {
      return (*m_face_transfer_operators[iF]).extension;
    }

    /// Return the reduction for the face of index iF
    inline const Eigen::MatrixXd & RcurlFace(size_t iF) const
    {
      return (*m_face_transfer_operators[iF]).reduction;
    }

    /// Return the serendipity reconstruction for face F
    inline const Eigen::MatrixXd & ScurlFace(const Face & F) const
    {
      return ScurlFace(F.global_index());
    }

    /// Return the extension for face F
    inline const Eigen::MatrixXd & EcurlFace(const Face & F) const
    {
      return EcurlFace(F.global_index());
    }

    /// Return cell reduction for cell T
    inline const Eigen::MatrixXd & RcurlFace(const Face & F) const
    {
      return RcurlFace(F.global_index());
    }
    
    //----------------------------------//
    //---- Cell transfer operators -----//
    /// Return the serendipity reconstruction for the cell of index iT
    inline const Eigen::MatrixXd & ScurlCell(size_t iT) const
    {
      return (*m_cell_transfer_operators[iT]).serendipity;
    }

    /// Return the extension for the cell of index iT
    inline const Eigen::MatrixXd & EcurlCell(size_t iT) const
    {
      return (*m_cell_transfer_operators[iT]).extension;
    }

    /// Return the reduction for the cell of index iT
    inline const Eigen::MatrixXd & RcurlCell(size_t iT) const
    {
      return (*m_cell_transfer_operators[iT]).reduction;
    }

    /// Return the serendipity reconstruction for cell T
    inline const Eigen::MatrixXd & ScurlCell(const Cell & T) const
    {
      return ScurlCell(T.global_index());
    }

    /// Return the extension for cell T
    inline const Eigen::MatrixXd & EcurlCell(const Cell & T) const
    {
      return EcurlCell(T.global_index());
    }

    /// Return the reduction for cell T
    inline const Eigen::MatrixXd & RcurlCell(const Cell & T) const
    {
      return RcurlCell(T.global_index());
    }

    //-----------------------------------------------------------------//
    //---- Full curl and potential reconstructions, and L2 product ----//
    
    /// Return the full curl operator on the face of index iF
    inline const Eigen::MatrixXd faceCurl(size_t iF) const
    {
      Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
      return Eigen::KroneckerProduct(m_sxcurl.faceCurl(iF), id);
    }
    
    /// Return the full curl operator on face F
    inline const Eigen::MatrixXd faceCurl(const Face & F) const
    {
      return faceCurl(F.global_index());
    }

    /// Return the potential operator on the face of index iF
    inline const Eigen::MatrixXd facePotential(size_t iF) const
    {
      Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
      return Eigen::KroneckerProduct(m_sxcurl.facePotential(iF), id);
    }

    /// Return the potential operator on face F
    inline const Eigen::MatrixXd facePotential(const Face & F) const
    {
      return facePotential(F.global_index());
    }

    /// Return the full curl operator on the cell of index iT
    inline const Eigen::MatrixXd cellCurl(size_t iT) const
    {
      Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
      return Eigen::KroneckerProduct(m_sxcurl.cellCurl(iT), id);
    }

    /// Return the full curl operator on cell T
    inline const Eigen::MatrixXd cellCurl(const Cell & T) const
    {
      return cellCurl(T.global_index());
    }

    /// Return the potential operator on the cell of index iT
    inline const Eigen::MatrixXd cellPotential(size_t iT) const
    {
      Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
      return Eigen::KroneckerProduct(m_sxcurl.cellPotential(iT), id);
    }

    /// Return the potential operator on cell T
    inline const Eigen::MatrixXd cellPotential(const Cell & T) const
    {
      return cellPotential(T.global_index());
    }

    /// Compute the matrix of the (weighted) L2-product
    Eigen::MatrixXd computeL2Product(
                                     const size_t iT, ///< index of the cell
                                     const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                     const Eigen::MatrixXd & mass_Pk3_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of (P^k(T))^3; if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to 1
                                     ) const
    {
      return Eigen::KroneckerProduct(m_sxcurl.computeL2Product(iT, penalty_factor, mass_Pk3_T, weight), m_lie_algebra.massMatrix());
    }

         
    /// Compute the matrix of the (weighted) L2-product as 'computeL2Product', with application of the discrete gradient on the left/right/both sides (depending on argument "side").
    Eigen::MatrixXd computeL2ProductGradient(
                                     const size_t iT, ///< index of the cell
                                     const SXGrad & sx_grad, ///< instance of SXGrad to access the full gradients
                                     const std::string & side, ///< which side (left, right, both) we apply the gradient to
                                     const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                     const Eigen::MatrixXd & mass_Pk3_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of (P^k(T))^3; if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to constant 1.
                                     ) const;

    //-----------------------//       
    //---- Getters ----------//

    /// Return cell bases for the face of index iT
    inline const DDRCore::CellBases & cellBases(size_t iT) const
    {
      return m_sxcurl.cellBases(iT);
    }

    /// Return cell bases for cell T
    inline const DDRCore::CellBases & cellBases(const Cell & T) const
    {
      return m_sxcurl.cellBases(T.global_index());
    }
    
    /// Return face bases for the face of index iF
    inline const DDRCore::FaceBases & faceBases(size_t iF) const
    {
      return m_sxcurl.faceBases(iF);
    }

    /// Return cell bases for face F
    inline const DDRCore::FaceBases & faceBases(const Face & F) const
    {
      return m_sxcurl.faceBases(F.global_index());
    }
    
    /// Return edge bases for the edge of index iE
    inline const DDRCore::EdgeBases & edgeBases(size_t iE) const
    {
      return m_sxcurl.edgeBases(iE);
    }

    /// Return edge bases for edge E
    inline const DDRCore::EdgeBases & edgeBases(const Edge & E) const
    {
      return m_sxcurl.edgeBases(E.global_index());
    }

  private:
    TransferOperators _compute_face_transfer_operators(size_t iF);
    TransferOperators _compute_cell_transfer_operators(size_t iT);
    
    const LieAlgebra & m_lie_algebra;
    const SXCurl & m_sxcurl;
       
    // Containers for serendipity, extension and reduction operators
    std::vector<std::unique_ptr<TransferOperators> > m_face_transfer_operators;
    std::vector<std::unique_ptr<TransferOperators> > m_cell_transfer_operators;

    bool m_use_threads;
    std::ostream & m_output;

  };

} // end of namespace HArDCore3D
#endif
