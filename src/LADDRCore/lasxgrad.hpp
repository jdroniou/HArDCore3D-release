#ifndef LASXGRAD_HPP
#define LASXGRAD_HPP

#include <sxgrad.hpp>
#include <liealgebra.hpp>
#include <unsupported/Eigen/KroneckerProduct>

namespace HArDCore3D
{
  /*!
   *  \addtogroup LADDRCore
   * @{
   */

  /// Discrete Serendipity Hgrad space: local operators, L2 product and global interpolator
  /** On each edge/face/element, the DOFs (if any) correspond to the polynomial bases on the edge/face/element provided by m_ddr_core */
  class LASXGrad : public VariableDOFSpace
  {
  public:
    typedef std::function<double(const Eigen::Vector3d &)> FunctionType;
    typedef std::vector<FunctionType> LAFunctionType;

    /// A structure to store the serendipity, extension and reduction operators
    /** The extension operators do SXGrad -> XGrad. The reduction operators do not produce all unknowns, only the ones in the mesh entities they're associated to; so they do XGrad -> \f$P^{\ell_P}(P)\f$ where P is a face or an element (the only use of the reduction operator is to compute the interpolate, and this version is what is useful). */
    struct TransferOperators
    {
      TransferOperators(
                     const Eigen::MatrixXd & _serendipity,    ///< Serendipity reconstruction
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
    LASXGrad(const LieAlgebra & lie_algebra, const SXGrad & sxgrad, bool use_threads = true, std::ostream & output = std::cout);

    /// Return the mesh
    const Mesh & mesh() const
    {
      return m_sxgrad.mesh();
    }
    
    /// Return the polynomial degree
    const size_t & degree() const
    {
      return m_sxgrad.degree();
    }

    /// Return the Lie algebra
    inline const LieAlgebra & lieAlg() const
    {
      return m_lie_algebra;
    }

    /// Return the serendipity operators
    const SerendipityProblem & serPro() const
    {
      return m_sxgrad.serPro();
    }

    /// Interpolator of a continuous function
    Eigen::VectorXd interpolate(
          const LAFunctionType & q, ///< The function to interpolate
          const int doe_cell = -1, ///< The optional degre of cell quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_face = -1, ///< The optional degre of face quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_edge = -1 ///< The optional degre of edge quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          ) const;

    //---------------------------------//
    //---- Face transfer operators-----//
    
    /// Return the serendipity reconstruction for the face of index iF
    inline const Eigen::MatrixXd & SgradFace(size_t iF) const
    {
      return (*m_face_transfer_operators[iF]).serendipity;
    }

    /// Return the extension for the face of index iF
    inline const Eigen::MatrixXd & EgradFace(size_t iF) const
    {
      return (*m_face_transfer_operators[iF]).extension;
    }

    /// Return the reduction for the face of index iF
    inline const Eigen::MatrixXd & RgradFace(size_t iF) const
    {
      return (*m_face_transfer_operators[iF]).reduction;
    }

    /// Return the serendipity reconstruction for face F
    inline const Eigen::MatrixXd & SgradFace(const Face & F) const
    {
      return SgradFace(F.global_index());
    }

    /// Return the extension for face F
    inline const Eigen::MatrixXd & EgradFace(const Face & F) const
    {
      return EgradFace(F.global_index());
    }

    /// Return cell reduction for cell T
    inline const Eigen::MatrixXd & RgradFace(const Face & F) const
    {
      return RgradFace(F.global_index());
    }
    
    //---------------------------------//
    //---- Cell transfer operators -----//
    /// Return the serendipity reconstruction for the cell of index iT
    inline const Eigen::MatrixXd & SgradCell(size_t iT) const
    {
      return (*m_cell_transfer_operators[iT]).serendipity;
    }

    /// Return the extension for the cell of index iT
    inline const Eigen::MatrixXd & EgradCell(size_t iT) const
    {
      return (*m_cell_transfer_operators[iT]).extension;
    }

    /// Return the reduction for the cell of index iT
    inline const Eigen::MatrixXd & RgradCell(size_t iT) const
    {
      return (*m_cell_transfer_operators[iT]).reduction;
    }

    /// Return the serendipity reconstruction for cell T
    inline const Eigen::MatrixXd & SgradCell(const Cell & T) const
    {
      return SgradCell(T.global_index());
    }

    /// Return the extension for cell T
    inline const Eigen::MatrixXd & EgradCell(const Cell & T) const
    {
      return EgradCell(T.global_index());
    }

    /// Return the reduction for cell T
    inline const Eigen::MatrixXd & RgradCell(const Cell & T) const
    {
      return RgradCell(T.global_index());
    }

    //---------------------------------------------------------------------//
    //---- Full gradient and potential reconstructions, and L2 product ----//
    
    /// Return the full gradient operator on the edge of index iE
    inline const Eigen::MatrixXd edgeGradient(size_t iE) const
    {
      Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
      return Eigen::KroneckerProduct(m_sxgrad.edgeGradient(iE), id);
    }
    
    /// Return the full gradient operator on edge E
    inline const Eigen::MatrixXd edgeGradient(const Edge & E) const
    {
      return edgeGradient(E.global_index());
    }

    /// Return the potential operator on the edge of index iE
    inline const Eigen::MatrixXd edgePotential(size_t iE) const
    {
      Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
      return Eigen::KroneckerProduct(m_sxgrad.edgePotential(iE), id);
    }

    /// Return the potential operator on edge E
    inline const Eigen::MatrixXd edgePotential(const Edge & E) const
    {
      return edgePotential(E.global_index());
    }

    /// Return the full gradient operator on the face of index iF
    inline const Eigen::MatrixXd faceGradient(size_t iF) const
    {
      Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
      return Eigen::KroneckerProduct(m_sxgrad.faceGradient(iF), id);
    }
    
    /// Return the full gradient operator on face F
    inline const Eigen::MatrixXd faceGradient(const Face & F) const
    {
      return faceGradient(F.global_index());
    }

    /// Return the potential operator on the face of index iF
    inline const Eigen::MatrixXd facePotential(size_t iF) const
    {
      Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
      return Eigen::KroneckerProduct(m_sxgrad.facePotential(iF), id);
    }

    /// Return the potential operator on face F
    inline const Eigen::MatrixXd facePotential(const Face & F) const
    {
      return facePotential(F.global_index());
    }

    /// Return the full gradient operator on the cell of index iT
    inline const Eigen::MatrixXd cellGradient(size_t iT) const
    {
      Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
      return Eigen::KroneckerProduct(m_sxgrad.cellGradient(iT), id);
    }

    /// Return the full gradient operator on cell T
    inline const Eigen::MatrixXd cellGradient(const Cell & T) const
    {
      return cellGradient(T.global_index());
    }

    /// Return the potential operator on the cell of index iT
    inline const Eigen::MatrixXd cellPotential(size_t iT) const
    {
      Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
      return Eigen::KroneckerProduct(m_sxgrad.cellPotential(iT), id);
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
      return Eigen::KroneckerProduct(m_sxgrad.computeL2Product(iT,penalty_factor,mass_Pk3_T, weight), m_lie_algebra.massMatrix()); 
    }

         
    //-----------------------//       
    //---- Getters ----------//

    /// Return cell bases for the face of index iT
    inline const DDRCore::CellBases & cellBases(size_t iT) const
    {
      return m_sxgrad.cellBases(iT);
    }

    /// Return cell bases for cell T
    inline const DDRCore::CellBases & cellBases(const Cell & T) const
    {
      return m_sxgrad.cellBases(T.global_index());
    }
    
    /// Return face bases for the face of index iF
    inline const DDRCore::FaceBases & faceBases(size_t iF) const
    {
      return m_sxgrad.faceBases(iF);
    }

    /// Return cell bases for face F
    inline const DDRCore::FaceBases & faceBases(const Face & F) const
    {
      return m_sxgrad.faceBases(F.global_index());
    }
    
    /// Return edge bases for the edge of index iE
    inline const DDRCore::EdgeBases & edgeBases(size_t iE) const
    {
      return m_sxgrad.edgeBases(iE);
    }

    /// Return edge bases for edge E
    inline const DDRCore::EdgeBases & edgeBases(const Edge & E) const
    {
      return m_sxgrad.edgeBases(E.global_index());
    }

  private:
    TransferOperators _compute_face_transfer_operators(size_t iF);
    TransferOperators _compute_cell_transfer_operators(size_t iT);
    
    const LieAlgebra & m_lie_algebra;
    const SXGrad & m_sxgrad;
       
    // Containers for serendipity, extension and reduction operators
    std::vector<std::unique_ptr<TransferOperators> > m_face_transfer_operators;
    std::vector<std::unique_ptr<TransferOperators> > m_cell_transfer_operators;

    bool m_use_threads;
    std::ostream & m_output;

  };

} // end of namespace HArDCore3D
#endif
