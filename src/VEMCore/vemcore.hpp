// Core data structures and methods required to implement the de Rham VEM sequence in 3D
//
// Provides:
//  - Full and partial polynomial spaces on the element, faces, and edges
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//



#ifndef VEMCORE_HPP
#define VEMCORE_HPP

#include <memory>
#include <iostream>

#include <basis.hpp>
#include <polynomialspacedimension.hpp>

/*!	
 * @defgroup VEMCore 
 * @brief Classes providing tools to the de Rham VEM sequence
 */


namespace HArDCore3D
{

  /*!
   *	\addtogroup VEMCore
   * @{
   */


  //------------------------------------------------------------------------------

  /// Construct all polynomial spaces for the VEM sequence
  class VEMCore
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
      std::unique_ptr<ShiftedBasis<PolyBasisCellType> > Polyk0;
      std::unique_ptr<ShiftedBasis<PolyBasisCellType> > Polykpo0;
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

      std::unique_ptr<PolyBasisFaceType> Polykp2;
      std::unique_ptr<PolyBasisFaceType> Polykpo;
      std::unique_ptr<PolyBasisFaceType> Polyk;
      std::unique_ptr<PolyBasisFaceType> Polykmo;
      std::unique_ptr<ShiftedBasis<PolyBasisFaceType> > Polyk0;
      std::unique_ptr<TangentFamily<PolyBasisFaceType>> Polykpo2;
      std::unique_ptr<Poly2BasisFaceType> Polyk2;
      std::unique_ptr<RolyBasisFaceType> Rolykmo;
      std::unique_ptr<RolyComplBasisFaceType> RolyComplkmo;
      std::unique_ptr<RolyComplBasisFaceType> RolyComplk;
      std::unique_ptr<RolyComplBasisFaceType> RolyComplkpo;
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
    VEMCore(const Mesh & mesh, size_t K, bool use_threads = true, std::ostream & output = std::cout);
    
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

#endif // VEMCORE_HPP
