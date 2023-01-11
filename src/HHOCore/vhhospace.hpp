// Core data structures and methods required to implement the Hybrid High-Order in 3D, for vector-valued functions
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


#ifndef VHHOSPACE_HPP
#define VHHOSPACE_HPP

#include <iostream>
#include <globaldofspace.hpp>
#include <basis.hpp>
#include <polynomialspacedimension.hpp>


namespace HArDCore3D
{

  /*!
   *	\addtogroup HHOSpace
   * @{
   */


  //------------------------------------------------------------------------------

  // Type for selecting some cells (return 0 or 1 for each cell T)
  typedef std::function<bool(const Cell &)> CellSelection;
  static const CellSelection allcells = [](const Cell &)->bool {return true;};

  /// Class definition: polynomial bases and operators
  class VHHOSpace : public GlobalDOFSpace
  {
  public:
    // Types for element bases
    typedef Family<MonomialScalarBasisCell> PolyBasisCellType;
    typedef TensorizedVectorFamily<PolyBasisCellType, dimspace> PolydBasisCellType;
    typedef MatrixFamily<PolyBasisCellType, dimspace> PolydxdBasisCellType;

    // Types for face bases
    typedef Family<MonomialScalarBasisFace> PolyBasisFaceType;
    typedef Family<TensorizedVectorFamily<PolyBasisFaceType, dimspace>> PolydBasisFaceType;

    // Types for functions to interpolate
    typedef std::function<VectorRd(const VectorRd &)> FunctionType;
    
    /// Structure to store element bases
    /** 'Poly': basis of polynomial space.\n
        'k' and 'kpo' (k+1) determines the degree.\n
        'd' for vector-valued, 'dxd' for matrix-valued
      */
    struct CellBases
    {
      /// Geometric support
      typedef Cell GeometricSupport;

      std::unique_ptr<PolyBasisCellType> Polykpo;
      std::unique_ptr<PolyBasisCellType> Polyk;
      std::unique_ptr<PolydBasisCellType> Polykpod;
      std::unique_ptr<PolydBasisCellType> Polykd;
      std::unique_ptr<PolydxdBasisCellType> Polykdxd;
    };

    /// Structure to store face bases
    /** See CellBases for details */
    struct FaceBases
    {
      /// Geometric support
      typedef Face GeometricSupport;

      std::unique_ptr<PolyBasisFaceType> Polyk;
      std::unique_ptr<PolydBasisFaceType> Polykd;
    };
    
    /// A structure to store local operators (gradient, potential, stabilisation)
    struct LocalOperators
    {
      LocalOperators(
                     const Eigen::MatrixXd & _gradient, ///< Gradient operator
                     const Eigen::MatrixXd & _potential, ///< HHO Potential operator
                     const Eigen::MatrixXd & _stabilisation, ///< H1 Stabilisation bilinear form associated to the HHO potential
                     const Eigen::MatrixXd & _potential_div, ///< Potential reconstructed from divergence
                     const Eigen::MatrixXd & _stabilisation_div ///< L2 Stabilisation associated to _potential_div
                     )
        : gradient(_gradient),
          potential(_potential),
          stabilisation(_stabilisation),
          potential_div(_potential_div),
          stabilisation_div(_stabilisation_div)
      {
        // Do nothing
      }
      
      Eigen::MatrixXd gradient;
      Eigen::MatrixXd potential;
      Eigen::MatrixXd stabilisation;
      Eigen::MatrixXd potential_div;
      Eigen::MatrixXd stabilisation_div;
    };

    /// Constructor (with function to select cells in which boundary faces are used in the stabilisation)
    VHHOSpace(const Mesh & mesh, size_t K, const CellSelection & BoundaryStab, bool use_threads = true, std::ostream & output = std::cout);

    /// Overloaded constructor when the selection of boundary stabilisation is not entered (all boundary faces are then used)
    VHHOSpace(const Mesh & mesh, size_t K, bool use_threads = true, std::ostream & output = std::cout)
            : VHHOSpace(mesh, K, allcells, use_threads, output) {};
    
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

    /// Return the function to select the cells with boundary stabilisation
    const CellSelection & boundaryStab() const
    {
      return m_boundary_stab;
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

    /// Computes the values of the potential reconstruction at the mesh vertices
    std::vector<VectorRd> computeVertexValues(
                  const Eigen::VectorXd & u   ///< DOFs in the discrete space
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
    // Choice of boundary stabilisation
    CellSelection m_boundary_stab;
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

#endif // VHHOSPACE_HPP
