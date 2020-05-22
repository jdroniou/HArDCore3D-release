#ifndef XDIV_HPP
#define XDIV_HPP

#include <ddrcore.hpp>
#include <ddrspace.hpp>

namespace HArDCore3D
{
  /*!
   *	\addtogroup DDRcore
   * @{
   */

  /// Discrete Hdiv space
  class XDiv : public DDRSpace
  {
  public:
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> FunctionType;

    /// A structure to store local operators
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
    XDiv(const DDRCore & ddr_core, bool use_threads = true, std::ostream & output = std::cout);

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

  private:
    LocalOperators _compute_cell_divergence_potential(size_t iT);

    const DDRCore & m_ddr_core;
    bool m_use_threads;
    std::ostream & m_output;

    // Containers for local operators
    std::vector<std::unique_ptr<LocalOperators> > m_cell_operators;
  };
} // end of namespace HArDCore3D

#endif
