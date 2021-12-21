// Core data structures and methods required to implement the discrete de Rham sequence in 3D
//
// Provides:
//  - Full and partial polynomial spaces on the element, faces, and edges
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
#include <polynomialspacedimension.hpp>

/*!	
 * @defgroup DDRCore 
 * @brief Classes providing tools to the Discrete De Rham sequence
 */


namespace HArDCore3D
{

  /*!
   *	\addtogroup DDRCore
   * @{
   */


  //------------------------------------------------------------------------------

  /// Construct all polynomial spaces for the DDR sequence
  class DDRCore
  {
  public:
    // Types for element bases
    typedef Family<MonomialScalarBasisCell> PolyBasisCellType;
    typedef TensorizedVectorFamily<PolyBasisCellType, 3> Poly3BasisCellType;
    typedef Family<GradientBasis<ShiftedBasis<MonomialScalarBasisCell> > > GolyBasisCellType;
    typedef Family<GolyComplBasisCell> GolyComplBasisCellType;
    typedef Family<CurlBasis<GolyComplBasisCell>> RolyBasisCellType;
    typedef Family<RolyComplBasisCell> RolyComplBasisCellType;

    // Types for face bases
    typedef Family<MonomialScalarBasisFace> PolyBasisFaceType;
    typedef TangentFamily<PolyBasisFaceType> Poly2BasisFaceType;
    typedef Family<CurlBasis<ShiftedBasis<MonomialScalarBasisFace>>> RolyBasisFaceType;
    typedef Family<RolyComplBasisFace> RolyComplBasisFaceType;

    // Type for edge basis
    typedef Family<MonomialScalarBasisEdge> PolyBasisEdgeType;

    /// Structure to store element bases
    /** 'Poly': basis of polynomial space; 'Goly': gradient basis; 'Roly': curl basis.\n
        'k', 'kmo' (k-1) and 'kpo' (k+1) determines the degree.\n
        'Compl' for the complement of the corresponding 'Goly' or 'Roly' in the 'Poly' space */
    struct CellBases
    {
      /// Geometric support
      typedef Cell GeometricSupport;

      std::unique_ptr<PolyBasisCellType> Polykpo;
      std::unique_ptr<PolyBasisCellType> Polyk;
      std::unique_ptr<PolyBasisCellType> Polykmo;
      std::unique_ptr<Poly3BasisCellType> Polyk3;
      std::unique_ptr<GolyBasisCellType> Golykmo;
      std::unique_ptr<GolyComplBasisCellType> GolyComplk;
      std::unique_ptr<GolyComplBasisCellType> GolyComplkpo;
      std::unique_ptr<RolyBasisCellType>  Rolykmo;
      std::unique_ptr<RolyComplBasisCellType> RolyComplk;
      std::unique_ptr<RolyComplBasisCellType> RolyComplkp2;
    };

    /// Structure to store face bases
    /** See CellBases for details */
    struct FaceBases
    {
      /// Geometric support
      typedef Face GeometricSupport;

      std::unique_ptr<PolyBasisFaceType> Polykpo;
      std::unique_ptr<PolyBasisFaceType> Polyk;
      std::unique_ptr<PolyBasisFaceType> Polykmo;
      std::unique_ptr<Poly2BasisFaceType> Polyk2;
      std::unique_ptr<RolyBasisFaceType> Rolykmo;
      std::unique_ptr<RolyComplBasisFaceType> RolyComplk;
      std::unique_ptr<RolyComplBasisFaceType> RolyComplkp2;
    };

    /// Structure to store edge bases
    /** See CellBases for details */
    struct EdgeBases
    {
      /// Geometric support
      typedef Edge GeometricSupport;

      std::unique_ptr<PolyBasisEdgeType> Polykpo;
      std::unique_ptr<PolyBasisEdgeType> Polyk;
      std::unique_ptr<PolyBasisEdgeType> Polykmo;
    };    
    
    /// Constructor
    DDRCore(const Mesh & mesh, size_t K, bool use_threads = true, std::ostream & output = std::cout);
    
    /// Return a const reference to the mesh
    const Mesh & mesh() const
    {
      return m_mesh;
    }

    /// Return the polynomial degree
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
