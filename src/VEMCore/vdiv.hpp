#ifndef VDIV_HPP
#define VDIV_HPP

#include <globaldofspace.hpp>
#include <vemcore.hpp>
#include <integralweight.hpp>

namespace HArDCore3D
{
  /*!
   *  \addtogroup VEMCore
   * @{
   */

  /// Virtual Hdiv space: local operators, L2 product and global interpolator
    /** The order of the DOFs in each element are, in this order: 
       
      Mesh entity | Type     | Space              | Description
      ------------|----------|--------------------|-------------
      F           | function | \f$P^{k}(F)\f$     | projection of normal component of \f$\mathbf{w}\f$ 
      T           | div      | \f$P^{k}_0(T)\f$   | projection of \f$\mathrm{div}(\mathbf{w})\f$ 
      T           | function | \f$G^{c,k+1}(T)\f$ | projection of \f$\mathbf{w}\f$ 

   */
  class VDiv : public GlobalDOFSpace
  {
  public:
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> FunctionType;
    typedef std::function<double(const Eigen::Vector3d &)> DivergenceType;

    /// A structure to store the local operators (projections of div and function)
    struct LocalOperators
    {
      LocalOperators(
                     const Eigen::MatrixXd & _proj_divergence,     ///< Projecton in \f$P^k\f$ of div (or rot for faces) of the function
                     const Eigen::MatrixXd & _proj_function ///< Projection in \f$P^{k}\f$ (for faces) or \f$P^k\f$ (for elements) of the function
                     )
        : proj_divergence(_proj_divergence),
          proj_function(_proj_function)
      {
        // Do nothing
      }

      Eigen::MatrixXd proj_divergence;
      Eigen::MatrixXd proj_function;
    };

    /// Constructor
    VDiv(const VEMCore & vem_core, bool use_threads = true, std::ostream & output = std::cout);

    /// Return the mesh
    const Mesh & mesh() const
    {
      return m_vem_core.mesh();
    }
    
    /// Return the polynomial degree
    const size_t & degree() const
    {
      return m_vem_core.degree();
    }

    /// Interpolator of a continuous function
    Eigen::VectorXd interpolate(
          const FunctionType & v, ///< The function to interpolate
          const DivergenceType & div_v, ///< Its div
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
    inline const VEMCore::CellBases & cellBases(size_t iT) const
    {
      return m_vem_core.cellBases(iT);
    }

    /// Return cell bases for cell T
    inline const VEMCore::CellBases & cellBases(const Cell & T) const
    {
      return m_vem_core.cellBases(T.global_index());
    }
    
    /// Return face bases for the face of index iF
    inline const VEMCore::FaceBases & faceBases(size_t iF) const
    {
      return m_vem_core.faceBases(iF);
    }

    /// Return cell bases for face F
    inline const VEMCore::FaceBases & faceBases(const Face & F) const
    {
      return m_vem_core.faceBases(F.global_index());
    }
    
    /// Return edge bases for the edge of index iE
    inline const VEMCore::EdgeBases & edgeBases(size_t iE) const
    {
      return m_vem_core.edgeBases(iE);
    }

    /// Return edge bases for edge E
    inline const VEMCore::EdgeBases & edgeBases(const Edge & E) const
    {
      return m_vem_core.edgeBases(E.global_index());
    }

    /// Compute the matrix of the (weighted) L2-product for the cell of index iT. The stabilisation here is based on cell and face potentials.
    // The mass matrix of P^k(T)^3 is the most expensive mass matrix in the calculation of this norm, which
    // is why there's the option of passing it as parameter if it's been already pre-computed when the norm is called.
    Eigen::MatrixXd computeL2Product(
                                     const size_t iT, ///< index of the cell
                                     const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                     const Eigen::MatrixXd & mass_Pk3_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of (P^k(T))^3; if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to 1
                                     ) const;
                                     
  private:
    LocalOperators _compute_cell_operators(size_t iT);
    
    const VEMCore & m_vem_core;
    bool m_use_threads;
    std::ostream & m_output;

    // Containers for local operators
    std::vector<std::unique_ptr<LocalOperators> > m_cell_operators;
    std::vector<std::unique_ptr<LocalOperators> > m_face_operators;
    
  };

} // end of namespace HArDCore3D
#endif
