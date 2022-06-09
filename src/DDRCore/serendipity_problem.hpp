// Class to select edges/faces to create serendipity spaces and compute serendipity operator
//
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


#ifndef SERENDIPITY_PROBLEM_HPP
#define SERENDIPITY_PROBLEM_HPP

#include <ddrcore.hpp>

namespace HArDCore3D
{

  /*!
   *	\addtogroup DDRCore
   * @{
   */


  //------------------------------------------------------------------------------

  /// Construct all polynomial spaces for the DDR sequence
  class SerendipityProblem
  {
  public:
    typedef RestrictedBasis<DDRCore::PolyBasisFaceType> PolylBasisFaceType;
    typedef RestrictedBasis<DDRCore::PolyBasisCellType> PolylBasisCellType;
    typedef RestrictedBasis<DDRCore::RolyComplBasisFaceType> RolyCompllpoBasisFaceType;
    typedef RestrictedBasis<DDRCore::RolyComplBasisCellType> RolyCompllpoBasisCellType;
    
    /// Type for inverses of matrix for serendipity problem
    /* PartialPivLU does not seem to always work... */
    typedef Eigen::FullPivLU<Eigen::MatrixXd> InverseProblem;
    
    /// Constructor
    SerendipityProblem(const DDRCore & ddrcore, bool use_threads = true, std::ostream & output = std::cout);
    
    //---------- Serendipity edges and faces, and degrees -------------//
    /// Return the list of serendipity edges in a face
    inline const std::vector<size_t> & serendipityEdges(size_t iF) const
    {
      return m_serendipity_edges[iF];
    }
    
    /// Return the number of serendipity edges in a face
    inline const int n_serendipityEdges(size_t iF) const
    {
      return m_serendipity_edges[iF].size();
    }

    /// Return the serendipity degree ell_F in a face
    inline const int serDegreeFace(size_t iF) const
    {
      return m_ddrcore.degree() + 1 - n_serendipityEdges(iF);
    }

    /// Return the list of serendipity face in a cell
    inline const std::vector<size_t> & serendipityFaces(size_t iT) const
    {
      return m_serendipity_faces[iT];
    }
    
    /// Return the number of serendipity faces in a cell
    inline const int n_serendipityFaces(size_t iT) const
    {
      return m_serendipity_faces[iT].size();
    }

    /// Return the serendipity degree ell_T in a cell
    inline const int serDegreeCell(size_t iT) const
    {
      return m_ddrcore.degree() + 1 - n_serendipityFaces(iT);
    }

    //---------- Serendipity dimensions and bases for grad -------------//
    /// Return the dimension of P^l on face of index iF
    inline size_t dimFacePolyl(size_t iF) const
    {
      return PolynomialSpaceDimension<Face>::Poly(serDegreeFace(iF));
    }

    /// Return the dimension of P^{l+1} on face F
    inline size_t dimFacePolyl(const Face & F) const
    {
      return dimFacePolyl(F.global_index());
    }

    /// Return the basis of P^l on face of index iF
    inline const PolylBasisFaceType & faceBasisPolyl(size_t iF) const
    {
      // Make sure that the basis has been created
      assert( m_face_bases_Polyl[iF] );
      return *m_face_bases_Polyl[iF].get();
 
    }

    /// Return the basis of P^l on face F
    inline const PolylBasisFaceType & faceBasisPolyl(const Face & F) const
    {
      return faceBasisPolyl(F.global_index()); 
    }
    
    /// Return the dimension of P^l on cell of index iT
    inline size_t dimCellPolyl(size_t iT) const
    {
      return PolynomialSpaceDimension<Cell>::Poly(serDegreeCell(iT));
    }

    /// Return the dimension of P^l on cell T
    inline size_t dimCellPolyl(const Cell & T) const
    {
      return dimCellPolyl(T.global_index());
    }

    /// Return the basis of P^l on cell of index iT
    inline const PolylBasisCellType & cellBasisPolyl(size_t iT) const
    {
      // Make sure that the basis has been created
      assert( m_cell_bases_Polyl[iT] );
      return *m_cell_bases_Polyl[iT].get();
 
    }

    /// Return the basis of P^l on cell T
    inline const PolylBasisCellType & cellBasisPolyl(const Cell & T) const
    {
      return cellBasisPolyl(T.global_index()); 
    }

    /// Number of DOFs on faces for serendipity XGrad space
    Eigen::VectorXd nDOFs_faces_SXGrad() const;
    /// Number of DOFs on cells for serendipity XGrad space
    Eigen::VectorXd nDOFs_cells_SXGrad() const;

    //---------- Serendipity dimensions and bases for curl -------------//
    /// Return the dimension of R^{c,l+1} on face of index iF
    inline size_t dimFaceRolyCompllpo(size_t iF) const
    {
      return PolynomialSpaceDimension<Face>::RolyCompl(serDegreeFace(iF)+1);
    }

    /// Return the dimension of R^{c,l+1} on face F
    inline size_t dimFaceRolyCompllpo(const Face & F) const
    {
      return dimFaceRolyCompllpo(F.global_index());
    }

    /// Return the basis of R^{c,l+1} on face of index iF
    inline const RolyCompllpoBasisFaceType & faceBasisRolyCompllpo(size_t iF) const
    {
      // Make sure that the basis has been created
      assert( m_face_bases_RolyCompllpo[iF] );
      return *m_face_bases_RolyCompllpo[iF].get();
 
    }

    /// Return the basis of R^{c,l+1} on face F
    inline const RolyCompllpoBasisFaceType & faceBasisRolyCompllpo(const Face & F) const
    {
      return faceBasisRolyCompllpo(F.global_index());
    }
    
    /// Return the dimension of R^{c,l+1} on cell of index iT
    inline size_t dimCellRolyCompllpo(size_t iT) const
    {
      return PolynomialSpaceDimension<Cell>::RolyCompl(serDegreeCell(iT)+1);
    }

    /// Return the dimension of R^{c,l+1} on cell T
    inline size_t dimCellRolyCompllpo(const Cell & T) const
    {
      return dimCellRolyCompllpo(T.global_index());
    }

    /// Return the basis of R^{c,l+1} on cell of index iT
    inline const RolyCompllpoBasisCellType & cellBasisRolyCompllpo(size_t iT) const
    {
      // Make sure that the basis has been created
      assert( m_cell_bases_RolyCompllpo[iT] );
      return *m_cell_bases_RolyCompllpo[iT].get();
 
    }

    /// Return the basis of R^{c,l+1} on cell T
    inline const RolyCompllpoBasisCellType & cellBasisRolyCompllpo(const Cell & T) const
    {
      return cellBasisRolyCompllpo(T.global_index()); 
    }
    
    /// Number of DOFs on faces for serendipity XCurl space
    Eigen::VectorXd nDOFs_faces_SXCurl() const;
    /// Number of DOFs on cells for serendipity XCurl space
    Eigen::VectorXd nDOFs_cells_SXCurl() const;

    //------------ Serendipity operators ---------------//
    /// Compute the serendipity operator on the face of index iF
    const Eigen::MatrixXd SerendipityOperatorFace(const size_t iF, const Eigen::MatrixXd & LF) const;

    /// Compute the serendipity operator on the face F
    inline const Eigen::MatrixXd SerendipityOperatorFace(const Face & F, const Eigen::MatrixXd & LF) const
    {
      return SerendipityOperatorFace(F.global_index(), LF);
    };

    /// Compute the serendipity operator on the cell of index iT
    const Eigen::MatrixXd SerendipityOperatorCell(const size_t iT, const Eigen::MatrixXd & LT) const;
    
    /// Compute the serendipity operator on the Cell T
    inline const Eigen::MatrixXd SerendipityOperatorCell(const Cell & T, const Eigen::MatrixXd & LT) const
    {
      return SerendipityOperatorCell(T.global_index(), LT);
    };


    /// Return a const reference to the mesh
    inline const Mesh & mesh() const
    {
      return m_ddrcore.mesh();
    }

  private:
    /// Create list of serendipity edge on face iF
    std::vector<size_t> _compute_serendipity_edges(size_t iF);

    /// Create inverse for serendipity problem on face iF
    InverseProblem _compute_inverse_problem_faces(size_t iF);

    /// Create list of serendipity edge on cell iT
    std::vector<size_t> _compute_serendipity_faces(size_t iT);

    /// Create inverse for serendipity problem on cell iT
    InverseProblem _compute_inverse_problem_cells(size_t iT);

    // Pointer do DDRCore
    const DDRCore & m_ddrcore;
    // Output stream
    std::ostream & m_output;

    // Face and cell basis P^l
    std::vector<std::unique_ptr<PolylBasisFaceType> > m_face_bases_Polyl;
    std::vector<std::unique_ptr<PolylBasisCellType> > m_cell_bases_Polyl;

    // Face and cell basis R^{c,l+1}
    std::vector<std::unique_ptr<RolyCompllpoBasisFaceType> > m_face_bases_RolyCompllpo;
    std::vector<std::unique_ptr<RolyCompllpoBasisCellType> > m_cell_bases_RolyCompllpo;

    // List the local indices of edges/faces used for the serendipity on each face/cell
    std::vector<std::vector<size_t>> m_serendipity_edges;  
    std::vector<std::vector<size_t>> m_serendipity_faces;
  
    // Store the inverse of the matrix to compute serendipity operators on faces/cells
    std::vector<InverseProblem> m_inverse_problem_faces;
    std::vector<InverseProblem> m_inverse_problem_cells;
    
  };
  
} // end of namespace HArDCore3D

#endif // SERENDIPITY_PROBLEM_HPP
