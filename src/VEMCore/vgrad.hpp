#ifndef VGRAD_HPP
#define VGRAD_HPP

#include <globaldofspace.hpp>
#include <vemcore.hpp>
#include <integralweight.hpp>

namespace HArDCore3D
{
  /*!
   *  \addtogroup VEMCore
   * @{
   */

  /// Virtual H1 space: local operators, L2 product and global interpolator
  /** The order of the DOFs in each element are: 
  
      Mesh entity | Type     | Space              | Description
      ------------|----------|--------------------|-------------
      V           | function | \f$\mathbb{R}\f$   | vertex value of \f$q\f$
      E           | function | \f$P^{k-1}(E)\f$   | projection of \f$q\f$
      F           | gradient | \f$R^{c,k+1}(F)\f$ | projection of \f$\mathrm{grad}_Fq\f$
      T           | gradient | \f$R^{c,k}(T)\f$   | projection of \f$\mathrm{grad}q\f$

   No serendipity here (hence the \f$R^{c,k+1}\f$ on the faces) */
  class VGrad : public GlobalDOFSpace
  {
  public:
    typedef std::function<double(const Eigen::Vector3d &)> FunctionType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> GradientType;

    /// A structure to store local operators (projection of gradient and function, dof of gradient in Vcurl)
    struct LocalOperators
    {
      LocalOperators(
                     const Eigen::MatrixXd & _proj_gradient, ///< Projection of gradient on \f$P^k\f$
                     const Eigen::MatrixXd & _dofs_gradient, ///< DOFs of gradient in Vcurl
                     const Eigen::MatrixXd & _proj_function ///< Projection of function (on \f$P^{k+1}\f$ for edge, \f$P^k\f$ for faces, \f$P^{k-1}\f$ for cells)
                     )
        : proj_gradient(_proj_gradient),
          dofs_gradient(_dofs_gradient),
          proj_function(_proj_function)
      {
        // Do nothing
      }
      
      Eigen::MatrixXd proj_gradient;
      Eigen::MatrixXd dofs_gradient;
      Eigen::MatrixXd proj_function;
    };
    
    /// Constructor
    VGrad(const VEMCore & vem_core, bool use_threads = true, std::ostream & output = std::cout);

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
          const FunctionType & q, ///< The function to interpolate
          const GradientType & grad_q, ///< Its gradient
          const int doe_cell = -1, ///< The optional degre of cell quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_face = -1, ///< The optional degre of face quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_edge = -1 ///< The optional degre of edge quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          ) const;

    /// Return edge operators for the edge of index iE
    inline const LocalOperators & edgeOperators(size_t iE) const
    {
      return *m_edge_operators[iE];
    }

    /// Return edge operators for edge E
    inline const LocalOperators & edgeOperators(const Edge & E) const
    {
      return *m_edge_operators[E.global_index()];
    }
    
    /// Return face operators for the face of index iF
    inline const LocalOperators & faceOperators(size_t iF) const
    {
      return *m_face_operators[iF];
    }

    /// Return face operators for face F
    inline const LocalOperators & faceOperators(const Face & F) const
    {
      return *m_face_operators[F.global_index()];
    }

    /// Return cell operators for the cell of index iT
    inline const LocalOperators & cellOperators(size_t iT) const
    {
      return *m_cell_operators[iT];
    }

    /// Return cell operators for cell T
    inline const LocalOperators & cellOperators(const Cell & T) const
    {
      return *m_cell_operators[T.global_index()];
    }

    /// Return cell bases for the cell of index iT
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

    /// Compute the matrix of the (weighted) L2-product for the cell of index iT. The stabilisation here is based on cell and face potentials. NOT COMPLETE: only computes the consistent component of this L2 product
    // The mass matrix of P^{k+1}(T) is the most expensive mass matrix in the calculation of this norm, which
    // is why there's the option of passing it as parameter if it's been already pre-computed when the norm is called.
    Eigen::MatrixXd computeL2Product(
                                     const size_t iT, ///< index of the cell
                                     const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                     const Eigen::MatrixXd & mass_Pkmo_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of P^{k-1}(T); if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to 1
                                     ) const;

    
  private:    
    LocalOperators _compute_edge_operators(size_t iE);
    LocalOperators _compute_face_operators(size_t iF);
    LocalOperators _compute_cell_operators(size_t iT);    

    const VEMCore & m_vem_core;
    bool m_use_threads;
    std::ostream & m_output;

    // Containers for local operators
    std::vector<std::unique_ptr<LocalOperators> > m_edge_operators;
    std::vector<std::unique_ptr<LocalOperators> > m_face_operators;
    std::vector<std::unique_ptr<LocalOperators> > m_cell_operators;
  };

} // end of namespace HArDCore3D

#endif
