#ifndef XGRAD_HPP
#define XGRAD_HPP

#include <ddrcore.hpp>
#include <ddrspace.hpp>

namespace HArDCore3D
{
  /*!
   *  \addtogroup DDRCore
   * @{
   */

  /// Discrete H1 space: local operators, L2 product and global interpolator
  class XGrad : public DDRSpace
  {
  public:
    typedef std::function<double(const Eigen::Vector3d &)> FunctionType;

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
    XGrad(const DDRCore & ddr_core, bool use_threads = true, std::ostream & output = std::cout);

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
    Eigen::VectorXd interpolate(const FunctionType & q) const;

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

    /// Build the components of the gradient operator
    Eigen::MatrixXd buildGradientComponentsCell(size_t iT) const;
    
  private:    
    LocalOperators _compute_edge_gradient_potential(size_t iE);
    LocalOperators _compute_face_gradient_potential(size_t iF);
    LocalOperators _compute_cell_gradient_potential(size_t iT);    

    const DDRCore & m_ddr_core;
    bool m_use_threads;
    std::ostream & m_output;

    // Containers for local operators
    std::vector<std::unique_ptr<LocalOperators> > m_edge_operators;
    std::vector<std::unique_ptr<LocalOperators> > m_face_operators;
    std::vector<std::unique_ptr<LocalOperators> > m_cell_operators;
  };

} // end of namespace HArDCore3D

#endif
