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

  /// Discrete Lie algebra valued H1 space: local operators, L2 product and global interpolator
  /** Each DOF in the XGrad space corresponds to the dimension of the Lie algebra number of DOFs in the LAXGrad space. These DOFs are stored consecutively with overall ordering given by XGrad; they represent and are locally ordered by the basis in the LieAlgebra */
  class LAXGrad : public GlobalDOFSpace
  {
  public:

    typedef std::function<double(const Eigen::Vector3d &)> FunctionType;
    typedef std::vector<FunctionType> LAFunctionType;

     /// A structure to store local operators (gradient and potential)
     struct LocalOperators
     {
       LocalOperators(
                      const Eigen::MatrixXd & _gradient, ///< Gradient operator
                      const Eigen::MatrixXd & _potential ///< Potential operator
                      )
         : gradient(_gradient),
           potential(_potential)
       {
         // Do nothing
       }
      
       Eigen::MatrixXd gradient;
       Eigen::MatrixXd potential;
     };
    
    /// Constructor
    LAXGrad(const LieAlgebra & lie_algebra, const XGrad & xgrad, bool use_threads = true, std::ostream & output = std::cout);
    
    /// Return the mesh
    const Mesh & mesh() const
    {
      return m_xgrad.mesh();
    }
    
    /// Return the polynomial degree
    const size_t & degree() const
    {
      return m_xgrad.degree();
    }
    
    /// Interpolator of a continuous function
    Eigen::VectorXd interpolate(
          const LAFunctionType & q, ///< The function to interpolate
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
    inline const DDRCore::CellBases & cellBases(size_t iT) const
    {
      return m_xgrad.cellBases(iT);
    }

    /// Return cell bases for cell T
    inline const DDRCore::CellBases & cellBases(const Cell & T) const
    {
      return m_xgrad.cellBases(T.global_index());
    }
    
    /// Return face bases for the face of index iF
    inline const DDRCore::FaceBases & faceBases(size_t iF) const
    {
      return m_xgrad.faceBases(iF);
    }

    /// Return cell bases for face F
    inline const DDRCore::FaceBases & faceBases(const Face & F) const
    {
      return m_xgrad.faceBases(F.global_index());
    }

    /// Return edge bases for the edge of index iE
    inline const DDRCore::EdgeBases & edgeBases(size_t iE) const
    {
      return m_xgrad.edgeBases(iE);
    }

    /// Return edge bases for edge E
    inline const DDRCore::EdgeBases & edgeBases(const Edge & E) const
    {
      return m_xgrad.edgeBases(E.global_index());
    }

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
    LocalOperators _compute_edge_gradient_potential(size_t iE);
    LocalOperators _compute_face_gradient_potential(size_t iF);
    LocalOperators _compute_cell_gradient_potential(size_t iT);

    const LieAlgebra & m_lie_algebra;
    const XGrad & m_xgrad;
    bool m_use_threads;
    std::ostream & m_output;
    
    // Containers for local operators
    std::vector<std::unique_ptr<LocalOperators> > m_edge_operators;
    std::vector<std::unique_ptr<LocalOperators> > m_face_operators;
    std::vector<std::unique_ptr<LocalOperators> > m_cell_operators;
  };

} // end of namespace HArDCore3D

#endif
