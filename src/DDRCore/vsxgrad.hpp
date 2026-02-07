#ifndef VSXGRAD_HPP
#define VSXGRAD_HPP

#include <sxgrad.hpp>
#include <unsupported/Eigen/KroneckerProduct>

namespace HArDCore3D
{
  /*!
   *  \addtogroup DDRCore
   * @{
   */
  
  /// Vector version of sXgrad, the arbitrary order space with nodal primal unknowns
  /** This space is built by tensorizing R^d with the sXgrad space. With \f$\phi=(\phi_i)_i\f$ the basis of sXgrad (local or global)
  and \f$e=(e_j)_j\f$ the canonical basis of \f$R^3\f$ then the basis of  
  vsXgrad is \f$e\otimes\phi=(e_1\phi_1,e_2\phi_1,e_3\phi_1,e_1\phi_2,\ldots)\f$.

  Given the ordering of DOFs in sXgrad, this ensures that we first list the DOFs for each vertex one by one,
  then the DOFs for each edge one by one, etc. Tensorising the other way means that the DOFs attached, say,
  to one vertex would be scrambled with DOFs attached to edges, etc.
  
  The class provides the matrices of the potentials (on edges, faces and cells) as well as the gradient (on the edges and cells)
  and symmetric gradient (in cells). The bases of the co-domain for these operators are those in Edge/Face/CellBases
  provided by this class (see the types provided by this class). We do not provide the face gradients as we do not,
  at the moment, have a suitable basis for this one (it would require a "MatrixFamily" for 2x3 non-square matrices).
  */
  class VSXGrad : public VariableDOFSpace
  {
  public:
    // Type of continuous functions corresponding to the discrete space VSXgrad
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> FunctionType;

    // Types for element bases
    typedef TensorizedVectorFamily<DDRCore::PolyBasisCellType, dimspace> Poly3BasisCellType;
    typedef MatrixFamily<DDRCore::PolyBasisCellType, dimspace> Poly3x3BasisCellType;
    
    // Types for face bases
    typedef TensorizedVectorFamily<DDRCore::PolyBasisFaceType, dimspace> Poly3BasisFaceType;
  
    // Type for edge basis
    typedef TensorizedVectorFamily<DDRCore::PolyBasisEdgeType, dimspace> Poly3BasisEdgeType;

    // Structure to store vectorial/matrix cell bases (the scalar ones are in m_sxgrad)
    struct CellBases
    {
      std::unique_ptr<Poly3BasisCellType> Polykpo3;
      std::unique_ptr<Poly3x3BasisCellType> Polyk3x3;
      std::unique_ptr<Poly3BasisCellType> Polykmo3;
    };
  
    // Structure to store vectorial/matrix face bases (the scalar ones are in m_sxgrad)
    struct FaceBases
    {
      std::unique_ptr<Poly3BasisFaceType> Polykpo3;
      std::unique_ptr<Poly3BasisFaceType> Polykmo3;
    };
  
    // Structure to store vectorial/matrix edge bases (the scalar ones are in m_sxgrad)
    struct EdgeBases
    {
      std::unique_ptr<Poly3BasisEdgeType> Polyk3;
      std::unique_ptr<Poly3BasisEdgeType> Polykpo3;
    };

    /// Constructor
    VSXGrad(const SerendipityProblem & sp, bool use_threads = true, std::ostream & output = std::cout);

    /// Return the mesh
    const Mesh & mesh() const
    {
      return m_sxgrad.mesh();
    }
    
    /// Return the underlying sXgrad space
    const SXGrad & sXgrad() const
    {
      return m_sxgrad;
    }
    
    /// Return the polynomial degree
    const size_t & degree() const
    {
      return m_sxgrad.degree();
    }

    /// Interpolator of a continuous function
    Eigen::VectorXd interpolate(
          const FunctionType & q, ///< The function to interpolate
          const int doe_cell = -1, ///< The optional degre of cell quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_face = -1, ///< The optional degre of face quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_edge = -1, ///< The optional degre of edge quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const bool boundary_only = false ///< If true, only boundary values will be computed, the other ones will be left to 0
          ) const;

    //---------------------------------------------------------------------//
    //---- Gradients and potential reconstructions, and stabilisation -----//
    
    /* Operators from sXgrad are easily extended into operators on vsXgrad using the Kronecker product
     These extension corresponds to the abovementioned tensorized basis of sXgrad and, for the codomain,
     to the tensorized basis of R^3 \otimes P (for P the polynomial space in which the operator of sXgrad
     lands). We need to do a change of basis of that codomain to represent the operators in the bases
     stored in the structures EdgeBases, FaceBases and CellBases.
    */
    
    /// Return the gradient operator (expressed on Polyk3) on the edge of index iE
    inline const Eigen::MatrixXd edgeGradient(size_t iE) const
    {
      // The codomain for the edge gradient of sXgrad is Polyk, so we need to change from
      // R^3 \otimes Polyk to Polyk3 (which is, by construction, Polyk \otimes R^3)
      return PermuteTensorization(m_sxgrad.edgeBases(iE).Polyk->dimension(), dimspace) 
                * Eigen::KroneckerProduct(m_sxgrad.edgeGradient(iE), Eigen::Matrix3d::Identity());
    }
    
    /// Return the gradient operator (expressed on Polyk3) on edge E
    inline const Eigen::MatrixXd edgeGradient(const Edge & E) const
    {
      return edgeGradient(E.global_index());
    }

    /// Return the potential operator (expressed on Polykpo3) on the edge of index iE
    inline const Eigen::MatrixXd edgePotential(size_t iE) const
    {
      // The codomain for the edge potential of sXgrad is Polykpo, so we need to change from
      // R^3 \otimes Polykpo to Polykpo3 which is, by construction, Polykpo \otimes R^3
      return PermuteTensorization(m_sxgrad.edgeBases(iE).Polykpo->dimension(), dimspace) 
                * Eigen::KroneckerProduct(m_sxgrad.edgePotential(iE), Eigen::Matrix3d::Identity());
    }

    /// Return the potential operator (expressed on Polykpo3) on edge E
    inline const Eigen::MatrixXd edgePotential(const Edge & E) const
    {
      return edgePotential(E.global_index());
    }

// Not implemented yet because we do not have a Matrix basis for non-square matrices
//    /// Return the full gradient operator on the face of index iF
//    inline const Eigen::MatrixXd faceGradient(size_t iF) const
//    {
//      Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
//      return Eigen::KroneckerProduct(m_sxgrad.faceGradient(iF), id);
//    }
    
//    /// Return the full gradient operator on face F
//    inline const Eigen::MatrixXd faceGradient(const Face & F) const
//    {
//      return faceGradient(F.global_index());
//    }

    /// Return the potential operator (expressed on Polykpo3) on the face of index iF
    inline const Eigen::MatrixXd facePotential(size_t iF) const
    {
      // The codomain for the face potential of sXgrad is Polykpo, so we need to change from
      // R^3 \otimes Polykpo to Polykpo3 which is, by construction, Polykpo \otimes R^3
      return PermuteTensorization(m_sxgrad.faceBases(iF).Polykpo->dimension(), dimspace) 
                * Eigen::KroneckerProduct(m_sxgrad.facePotential(iF), Eigen::Matrix3d::Identity());
    }

    /// Return the potential operator (expressed on Polykpo3) on face F
    inline const Eigen::MatrixXd facePotential(const Face & F) const
    {
      return facePotential(F.global_index());
    }

    /// Return the gradient operator (expressed on Polyk3x3) on the cell of index iT
    inline const Eigen::MatrixXd cellGradient(size_t iT) const
    {
      // The codomain for the cell gradient of sXgrad is Polyk3, so we need to change from
      // R^3 \otimes Polyk3 to Polyk3x3. 
      // Here, Polyk3 = Polyk \otimes R^3, where the R^3=R^3_nabla component represents the derivatives.
      // We first permute the codomain of m_sxgrad.cellGradient(iT) to make it R^3_nabla \otimes Polyk.
      Eigen::MatrixXd GradPermutedBasis 
          = PermuteTensorization(dimspace, m_sxgrad.cellBases(iT).Polyk->dimension()) * m_sxgrad.cellGradient(iT);

      // The Kronecker product then gives the cellGradient of vsXgrad in (R^3 \otimes R^3_nabla) \otimes Polyk,
      // where the first copy of R^3 lists the components of the functions in the vectorized vsXgrad.
      // We permute this basis to Polyk \otimes (R^3 \otimes R^3_nabla). In this basis, the first index
      // of R^3 \otimes R^3_nabla goes over the components of the functions, and the second over the coordinates
      // derivatives. 
      // Polyk3x3 is Polyk \otimes M_N(R) with the basis of M_N(R) listing the canonical
      // basis vectors in columns; so identifying the basis of R^3 \otimes R^3_nabla with the basis of M_N(R)
      // organises the functions components in columns and their gradients in rows, which is exactly what
      // we want. Hence, the matrix obtained on Polyk \otimes (R^3 \otimes R^3_nabla) can directly be interpreted
      // as the matrix on Polyk \otimes M_N(R) = Polyk3x3.
      return PermuteTensorization(m_sxgrad.cellBases(iT).Polyk->dimension(), dimspace*dimspace) 
                * Eigen::KroneckerProduct(GradPermutedBasis, Eigen::Matrix3d::Identity());
    }

    /// Return the gradient operator (expressed on Polyk3x3) on cell T
    inline const Eigen::MatrixXd cellGradient(const Cell & T) const
    {
      return cellGradient(T.global_index());
    }

    /// Return the symmetric gradient operator (expressed on Polyk3x3) on the cell of index iT
    inline const Eigen::MatrixXd cellSymGradient(size_t iT) const
    {
      return cellBases(iT).Polyk3x3->symmetriseOperator() * cellGradient(iT);
    }

    /// Return the symmetric gradient operator (expressed on Polyk3x3) on cell T
    inline const Eigen::MatrixXd cellSymGradient(const Cell & T) const
    {
      return cellSymGradient(T.global_index());
    }

    /// Return the divergence operator (expressed on the ancester of Polyk3x3, which should be m_sxgrad.cellBases(iT).Polyk) on the cell of index iT
    inline const Eigen::MatrixXd cellDivergence(size_t iT) const
    {
      return cellBases(iT).Polyk3x3->traceOperator() * cellGradient(iT);
    }

    /// Return the divergence operator (expressed on the ancester of Polyk3x3, which should be m_sxgrad.cellBases(iT).Polyk) on the cell of index iT
    inline const Eigen::MatrixXd cellDivergence(const Cell & T) const
    {
      return cellDivergence(T.global_index());
    }

    /// Return the potential operator (expressed on Polykpo3) on the cell of index iT
    inline const Eigen::MatrixXd cellPotential(size_t iT) const
    {
      // The codomain for the cell potential of sXgrad is Polykpo, so we need to change from
      // R^3 \otimes Polykpo to Polykpo3 which is, by construction, Polykpo \otimes R^3
      return PermuteTensorization(m_sxgrad.cellBases(iT).Polykpo->dimension(), dimspace) 
                * Eigen::KroneckerProduct(m_sxgrad.cellPotential(iT), Eigen::Matrix3d::Identity());
    }

    /// Return the potential operator on cell T
    inline const Eigen::MatrixXd cellPotential(const Cell & T) const
    {
      return cellPotential(T.global_index());
    }

    /// Compute the H^1-stabilisation (based on the L^2 stabilisation implemented in SXgrad)
    Eigen::MatrixXd computeStabilisation(
                                     const size_t iT, ///< index of the cell
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to 1
                                     ) const
    {
      // The H^1-stabilisation is the L^2-stabilisation divided by h_T^2
      double hTm2 = std::pow(mesh().cell(iT)->diam(), -2);
      return hTm2 * Eigen::KroneckerProduct(m_sxgrad.computeStabilisation(iT, weight), Eigen::Matrix3d::Identity()); 
    }

    /// Compute the H^1 discrete norm, with full or symmetric gradient
    std::vector<double> computeH1norms(
                          const std::vector<Eigen::VectorXd> & list_dofs, ///< The list of vectors representing the dofs                          
                          const std::string typegrad="full" ///< computes the norm using the "full" or "symmetric" gradient
                          ) const;
         
    /// Computes the values of the potential reconstruction at the mesh vertices
    std::vector<VectorRd> computeVertexValues(
                  const Eigen::VectorXd & u   ///< DOFs in the discrete space
                  ) const;
         
    //-----------------------//       
    //---- Getters ----------//

    /// Return vector/matrix cell bases for the face of index iT
    inline const CellBases & cellBases(size_t iT) const
    {
      assert( m_cell_bases[iT] );
      return *m_cell_bases[iT].get();
    }

    /// Return vector/matrix cell bases for cell T
    inline const CellBases & cellBases(const Cell & T) const
    {
      return cellBases(T.global_index());
    }
    
    /// Return vector face bases for the face of index iF
    inline const FaceBases & faceBases(size_t iF) const
    {
      assert( m_face_bases[iF] );
      return *m_face_bases[iF].get();
    }

    /// Return vector bases for face F
    inline const FaceBases & faceBases(const Face & F) const
    {
      return faceBases(F.global_index());
    }
    
    /// Return edge bases for the edge of index iE
    inline const EdgeBases & edgeBases(size_t iE) const
    {
      assert( m_edge_bases[iE] );
      return *m_edge_bases[iE].get();
    }

    /// Return edge bases for edge E
    inline const EdgeBases & edgeBases(const Edge & E) const
    {
      return edgeBases(E.global_index());
    }

  private:
    EdgeBases _compute_edge_bases(size_t iE);
    FaceBases _compute_face_bases(size_t iF);
    CellBases _compute_cell_bases(size_t iT);
    
    const SXGrad m_sxgrad; // Not a reference as it is not created outside this class
       
    // Cell bases
    std::vector<std::unique_ptr<CellBases> > m_cell_bases;
    // Face bases
    std::vector<std::unique_ptr<FaceBases> > m_face_bases;
    // Edge bases
    std::vector<std::unique_ptr<EdgeBases> > m_edge_bases;

    bool m_use_threads;
    std::ostream & m_output;

  };

} // end of namespace HArDCore3D
#endif
