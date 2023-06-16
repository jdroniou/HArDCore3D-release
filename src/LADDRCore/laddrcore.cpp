#include <laddrcore.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>

using namespace HArDCore3D;

LADDRCore::LADDRCore(const LieAlgebra & liealgebra, const Mesh & mesh, size_t K, bool use_threads, std::ostream & output)
  : DDRCore({mesh, K, use_threads, output}),
    m_lie_algebra(liealgebra),
    m_output(output),
    m_P2k3_cell_basis(mesh.n_cells())
{
  m_output << "[LADDRCore] Initializing" << std::endl;

  // Construct element bases
  std::function<void(size_t, size_t)> construct_extra_cell_bases
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iT = start; iT < end; iT++) {
	        this->m_P2k3_cell_basis[iT].reset( new Poly3BasisCellType(this->_construct_P2k3_cell_bases(iT)) );
	      } // for iT
      };

  m_output << "[LADDRCore] Constructing P2k3 element bases" << std::endl;
  parallel_for(this->mesh().n_cells(), construct_extra_cell_bases, use_threads);
}

DDRCore::Poly3BasisCellType LADDRCore::_construct_P2k3_cell_bases(size_t iT)
{
  const Cell & T = *this->mesh().cell(iT);
  
  MonomialCellIntegralsType int_monoT_4k = IntegrateCellMonomials(T, 4*this->degree());
  
  //------------------------------------------------------------------------------
  // Basis for P2k(T)
  //------------------------------------------------------------------------------

  MonomialScalarBasisCell basis_P2k_T(T, 2*this->degree());
  PolyBasisCellType Poly2k = l2_orthonormalize(basis_P2k_T, GramMatrix(T, basis_P2k_T, int_monoT_4k));  

  // Check that we got the dimensions right
  assert( Poly2k.dimension() == PolynomialSpaceDimension<Cell>::Poly(2*this->degree()) );
  
  //------------------------------------------------------------------------------
  // Basis for Pk(T)^3
  //------------------------------------------------------------------------------

  Poly3BasisCellType Poly2k3(Poly2k);
  assert( Poly2k3.dimension() == 3 * PolynomialSpaceDimension<Cell>::Poly(2*this->degree()) );

  return Poly2k3;
}