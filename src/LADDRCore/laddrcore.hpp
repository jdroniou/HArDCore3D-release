#ifndef LADDRCORE_HPP
#define LADDRCORE_HPP

#include <ddrcore.hpp>
#include <liealgebra.hpp>

namespace HArDCore3D
{

  /*!
   *  \addtogroup LADDRCore
   * @{
   */

  //------------------------------------------------------------------------------

  /// Construct the spaces for the LADDR sequence
  class LADDRCore : public DDRCore
  {
  public:
    
    /// Constructor
    LADDRCore(const LieAlgebra & lie_algebra, const Mesh & mesh, size_t K, bool use_threads = true, std::ostream & output = std::cout);
    
    /// Return cell bases for element iT
    inline const DDRCore::Poly3BasisCellType & P2k3(size_t iT) const
    {
      // Make sure that the basis has been created
      assert( m_P2k3_cell_basis[iT] );
      return *m_P2k3_cell_basis[iT].get();
    }
    
    /// Return the Lie algebra
    inline const LieAlgebra & lieAlg() const
    {
      return m_lie_algebra;
    }

  private:
    /// Compute the bases on an element T
    Poly3BasisCellType _construct_P2k3_cell_bases(size_t iT);

    // Lie algebra
    const LieAlgebra & m_lie_algebra;

    // Output stream
    std::ostream & m_output;

    // Cell bases
    std::vector<std::unique_ptr<DDRCore::Poly3BasisCellType>> m_P2k3_cell_basis;
  };
  
} // end of namespace HArDCore3D

#endif // DDRCORE_HPP
