// Core data structures and methods required to implement the discrete de Rham sequence in 3D
//
// Provides:
//  - Full and partial polynomial spaces on the element, faces, and edges
//  - Generic routines to create quadrature nodes over cells and faces of the mesh
//
// Author: Daniele Di Pietro (daniele.di-pietro@umontpellier.fr)
//

/*
 *
 *      This library was developed around HHO methods, although some parts of it have a more
 * general purpose. If you use this code or part of it in a scientific publication, 
 * please mention the following book as a reference for the underlying principles
 * of HHO schemes:
 *
 * The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications. 
 *  D. A. Di Pietro and J. Droniou. Modeling, Simulation and Applications, vol. 19. 
 *  Springer International Publishing, 2020, xxxi + 525p. doi: 10.1007/978-3-030-37203-3. 
 *  url: https://hal.archives-ouvertes.fr/hal-02151813.
 *
 */

/*
 * The DDR sequence has been designed in
 *
 *  Fully discrete polynomial de Rham sequences of arbitrary degree on polygons and polyhedra.
 *   D. A. Di Pietro, J. Droniou, and F. Rapetti, 33p, 2019. url: https://arxiv.org/abs/1911.03616.
 *
 * If you use this code in a scientific publication, please mention the above article.
 *
 */
 

#ifndef DDRCORE_HPP
#define DDRCORE_HPP

#include <memory>
#include <iostream>

#include <basis.hpp>

/*!	
 * @defgroup DDRcore 
 * @brief Classes providing tools to the Discrete De Rham sequence
 */


namespace HArDCore3D
{

  /*!
   *	\addtogroup DDRcore
   * @{
   */

  /// Basis dimensions for the full and trimmed polynomial spaces
  template<typename GeometricSupport>
  struct PolynomialSpaceDimension
  {
    // Only specializations are relevant
  };
  
  template<>
  struct PolynomialSpaceDimension<Cell>
  {
    /// Dimension of Pk(T)
    static size_t Poly(int k)
    {
      return (k >= 0 ? (k * (k + 1) * (2 * k + 1) + 9 * k * (k + 1) + 12 * (k + 1)) / 12 : 0);
    }
    /// Dimension of Gk(T)
    static size_t Goly(int k)
    {
      return (k >= 0 ? PolynomialSpaceDimension<Cell>::Poly(k + 1) - 1 : 0);
    }
    /// Dimension of GOk(T)
    static size_t GolyOrth(int k)
    {
      return 3 * PolynomialSpaceDimension<Cell>::Poly(k) - PolynomialSpaceDimension<Cell>::Goly(k);
    }
    /// Dimension of Rk(T)
    static size_t Roly(int k)
    {
      return (k >= 0 ? PolynomialSpaceDimension<Cell>::GolyOrth(k + 1) : 0);
    }
    /// Dimension of ROk(T)
    static size_t RolyOrth(int k)
    {
      return 3 * PolynomialSpaceDimension<Cell>::Poly(k) - PolynomialSpaceDimension<Cell>::Roly(k);
    }
  };

  template<>
  struct PolynomialSpaceDimension<Face>
  {
    /// Dimension of Pk(F)
    static size_t Poly(int k)
    {
      return (k >= 0 ? (k + 1) * (k + 2) / 2 : 0);
    }
    /// Dimension of Rk(F)
    static size_t Roly(int k)
    {
      return (k >= 0 ? PolynomialSpaceDimension<Face>::Poly(k + 1) - 1 : 0);
    }
    /// Dimension of ROk(F)
    static size_t RolyOrth(int k)
    {
      return 2 * PolynomialSpaceDimension<Face>::Poly(k) - PolynomialSpaceDimension<Face>::Roly(k);
    }
  };

  template<>
  struct PolynomialSpaceDimension<Edge>
  {
    /// Dimension of Pk(E)
    static size_t Poly(int k)
    {
      return (k >= 0 ? k + 1 : 0);
    }
  };

  //------------------------------------------------------------------------------

  /// Construct the spaces for the DDR sequence
  class DDRCore
  {
  public:
    // Types for element bases
    typedef Family<MonomialScalarBasisCell> PolyCellBasisType;
    typedef TensorizedVectorFamily<PolyCellBasisType, 3> Poly3CellBasisType;
    typedef Family<GradientBasis<ShiftedBasis<MonomialScalarBasisCell> > > GolyCellBasisType;
    typedef Family<Poly3CellBasisType> GolyOrthCellBasisType;
    typedef CurlBasis<GolyOrthCellBasisType> RolyCellBasisType;
    typedef Family<Poly3CellBasisType> RolyOrthCellBasisType;

    // Types for face bases
    typedef Family<MonomialScalarBasisFace> PolyFaceBasisType;
    typedef TangentFamily<PolyFaceBasisType> Poly2FaceBasisType;
    typedef CurlBasis<ShiftedBasis<PolyFaceBasisType> > RolyFaceBasisType;
    typedef Family<Poly2FaceBasisType> RolyOrthFaceBasisType;
    typedef Family<Family<TangentFamily<MonomialScalarBasisFace> > > PotentialTestFunctionBasisType;

    // Type for edge basis
    typedef Family<MonomialScalarBasisEdge> PolyEdgeBasisType;

    /// Structure for element bases
    struct CellBases
    {
      /// Geometric support
      typedef Cell GeometricSupport;

      std::unique_ptr<PolyCellBasisType> Polykpo;
      std::unique_ptr<RestrictedBasis<PolyCellBasisType> > Polyk;
      std::unique_ptr<RestrictedBasis<PolyCellBasisType> > Polykmo;
      std::unique_ptr<Poly3CellBasisType> Polyk3;
      std::unique_ptr<GolyCellBasisType> Golykmo;
      std::unique_ptr<GolyOrthCellBasisType> GolyOrthk;
      std::unique_ptr<GolyOrthCellBasisType> GolyOrthkpo;
      std::unique_ptr<RolyCellBasisType>  Rolykmo;
      std::unique_ptr<RolyOrthCellBasisType> RolyOrthk;
    };

    /// Structure for face bases
    struct FaceBases
    {
      /// Geometric support
      typedef Face GeometricSupport;

      std::unique_ptr<PolyFaceBasisType> Polykpo;
      std::unique_ptr<RestrictedBasis<PolyFaceBasisType> > Polyk;
      std::unique_ptr<RestrictedBasis<PolyFaceBasisType> > Polykmo;
      std::unique_ptr<Poly2FaceBasisType> Polyk2;
      std::unique_ptr<RolyFaceBasisType> Rolykmo;
      std::unique_ptr<RolyOrthFaceBasisType> RolyOrthk;
      std::unique_ptr<PotentialTestFunctionBasisType> RolyOrthkp2;
    };

    /// Structure for edge bases
    struct EdgeBases
    {
      /// Geometric support
      typedef Edge GeometricSupport;

      std::unique_ptr<PolyEdgeBasisType> Polykpo;
      std::unique_ptr<RestrictedBasis<PolyEdgeBasisType> > Polyk;
      std::unique_ptr<RestrictedBasis<PolyEdgeBasisType> > Polykmo;
    };    
    
    /// Constructor
    DDRCore(const Mesh & mesh, size_t K, bool use_threads = true, std::ostream & output = std::cout);
    
    // Return a const reference to the mesh
    const Mesh & mesh() const
    {
      return m_mesh;
    }

    // Return the polynomial degree
    const size_t & degree() const
    {
      return m_K;
    }
    
    /// Return cell bases for element iT
    inline const CellBases & cellBases(size_t iT) const
    {
      // Make sure that the basis has been created
      assert( m_cell_bases[iT] );
      return *m_cell_bases[iT].get();
    }

    /// Return face bases for face iF
    inline const FaceBases & faceBases(size_t iF) const
    {
      // Make sure that the basis has been created
      assert( m_face_bases[iF] );
      return *m_face_bases[iF].get();
    }

    /// Return edge bases for edge iE
    inline const EdgeBases & edgeBases(size_t iE) const
    {
      // Make sure that the basis has been created
      assert( m_edge_bases[iE] );
      return *m_edge_bases[iE].get();
    }

  private:
    /// Compute the bases on an element T
    CellBases _construct_cell_bases(size_t iT);

    /// Compute the bases on a face F
    FaceBases _construct_face_bases(size_t iF);

    /// Compute the bases on an edge E
    EdgeBases _construct_edge_bases(size_t iE);
    
    // Pointer to the mesh
    const Mesh & m_mesh;
    // Degree
    const size_t m_K;
    // Output stream
    std::ostream & m_output;
    
    // Cell bases
    std::vector<std::unique_ptr<CellBases> > m_cell_bases;
    // Face bases
    std::vector<std::unique_ptr<FaceBases> > m_face_bases;
    // Edge bases
    std::vector<std::unique_ptr<EdgeBases> > m_edge_bases;
  };
  
} // end of namespace HArDCore3D

#endif // DDRCORE_HPP
