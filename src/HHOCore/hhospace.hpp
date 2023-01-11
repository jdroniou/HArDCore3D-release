// Core data structures and methods required to implement the Hybrid High-Order in 3D
//
// Provides:
//  - Polynomial spaces on the element and faces
//  - Interpolator of smooth functions
//  - Full gradient, potential and stabilisation bilinear form in the elements
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
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


#ifndef HHOSPACE_HPP
#define HHOSPACE_HPP

#include <iostream>
#include <globaldofspace.hpp>
#include <basis.hpp>
#include <polynomialspacedimension.hpp>

/*!	
 * @defgroup HHOSpace
 * @brief Classes defining the HHO method (scalar and vector-valued)
 */


namespace HArDCore3D
{

  /*!
   *	\addtogroup HHOSpace
   * @{
   */


  //------------------------------------------------------------------------------

  /// Class definition: polynomial bases and operators
  class HHOSpace : public GlobalDOFSpace
  {
  public:
    // Types for element bases
    typedef Family<MonomialScalarBasisCell> PolyBasisCellType;
    typedef TensorizedVectorFamily<PolyBasisCellType, dimspace> PolydBasisCellType;

    // Types for face bases
    typedef Family<MonomialScalarBasisFace> PolyBasisFaceType;

    // Types for functions to interpolate
    typedef std::function<double(const VectorRd &)> FunctionType;

    /// Structure to store element bases
    /** 'Poly': basis of polynomial space.\n
        'k' and 'kpo' (k+1) determines the degree.\n
        'd' for vector-valued.
      */
    struct CellBases
    {
      /// Geometric support
      typedef Cell GeometricSupport;

      std::unique_ptr<PolyBasisCellType> Polykpo;
      std::unique_ptr<PolyBasisCellType> Polyk;
      std::unique_ptr<PolydBasisCellType> Polykd;
    };

    /// Structure to store face bases
    /** See CellBases for details */
    struct FaceBases
    {
      /// Geometric support
      typedef Face GeometricSupport;

      std::unique_ptr<PolyBasisFaceType> Polyk;
    };
    
    /// A structure to store local operators (gradient, potential, stabilisation)
    struct LocalOperators
    {
      LocalOperators(
                     const Eigen::MatrixXd & _gradient, ///< Gradient operator
                     const Eigen::MatrixXd & _potential, ///< Potential operator
                     const Eigen::MatrixXd & _stabilisation ///< Stabilisation bilinear form
                     )
        : gradient(_gradient),
          potential(_potential),
          stabilisation(_stabilisation)
      {
        // Do nothing
      }
      
      Eigen::MatrixXd gradient;
      Eigen::MatrixXd potential;
      Eigen::MatrixXd stabilisation;
    };

    /// Constructor
    HHOSpace(const Mesh & mesh, size_t K, bool use_threads = true, std::ostream & output = std::cout);
    
    /// Return a const reference to the mesh
    const Mesh & mesh() const
    {
      return m_mesh;
    }

    /// Return the polynomial degree (common face and elements)
    const size_t & degree() const
    {
      return m_K;
    }

    /// Interpolator of a continuous function
    Eigen::VectorXd interpolate(
          const FunctionType & q, ///< The function to interpolate
          const int doe_cell = -1, ///< The optional degre of cell quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_face = -1 ///< The optional degre of face quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          ) const;
    
    /// Return cell bases for element iT
    inline const CellBases & cellBases(size_t iT) const
    {
      // Make sure that the basis has been created
      assert( m_cell_bases[iT] );
      return *m_cell_bases[iT].get();
    }

    /// Return cell bases for cell T
    inline const CellBases & cellBases(const Cell & T) const
    {
      return cellBases(T.global_index());
    }
    
    /// Return face bases for face iF
    inline const FaceBases & faceBases(size_t iF) const
    {
      // Make sure that the basis has been created
      assert( m_face_bases[iF] );
      return *m_face_bases[iF].get();
    }

    /// Return cell bases for face F
    inline const FaceBases & faceBases(const Face & F) const
    {
      return faceBases(F.global_index());
    }

    /// Return operators for the cell of index iT
    inline const LocalOperators & operators(size_t iT) const
    {
      assert( m_operators[iT] );
      return *m_operators[iT];
    }

    /// Return cell operators for cell T
    inline const LocalOperators & operators(const Cell & T) const
    {
      return operators(T.global_index());
    }
    
    /// Computes the discrete L2 (cell unknowns only) and H1 norms of a list of vectors
    std::vector<std::pair<double,double>> computeNorms(
                   const std::vector<Eigen::VectorXd> & list_dofs   ///< The list of vectors representing the dofs
                  ) const;

    
  private:
    /// Compute the bases on an element T
    CellBases _construct_cell_bases(size_t iT);

    /// Compute the bases on a face F
    FaceBases _construct_face_bases(size_t iF);
    
    /// Compute operators in an element T
    LocalOperators _compute_operators(size_t iT);
    
    // Pointer to the mesh
    const Mesh & m_mesh;
    // Degrees
    const size_t m_K;
    // Parallel or not
    bool m_use_threads;
    // Output stream
    std::ostream & m_output;    
    
    // Cell bases
    std::vector<std::unique_ptr<CellBases> > m_cell_bases;
    // Face bases
    std::vector<std::unique_ptr<FaceBases> > m_face_bases;

    // Local operators
    std::vector<std::unique_ptr<LocalOperators> > m_operators;
        
  };
  
} // end of namespace HArDCore3D

#endif // HHOSPACE_HPP
