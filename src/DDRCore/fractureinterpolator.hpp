#ifndef FRACTUREINTERPOLATOR_HPP
#define FRACTUREINTERPOLATOR_HPP

#pragma once

#include<mesh.hpp>
#include<connected_information_generator.hpp>
#include<vsxgrad.hpp>

namespace HArDCore3D
{
  /*!
   * \addtogroup Common
   * @{
   */

  /// Base class for global DOF spaces. 
  /** The DOFs are organised by increasing geometric entities dimensions: DOFs of vertices, DOFs of edges, DOFs of faces, DOFs of cells.
  This version allows for cases where the number of DOF on each mesh entity varies from one mesh entity to the next */
  
  class FractureInterpolator : public VariableDOFSpace
  {
  public:
    typedef std::function<double(const Eigen::Vector3d &, const Cell *)> ScalarFunctionType;

    /// Constructor
    FractureInterpolator(
             const VSXGrad & vsx,
             const ConnectedInformationGenerator & info_generator,
             bool use_threads = true
             );


    //------------------------------------------------------------------------------
    // Accessors
    //------------------------------------------------------------------------------
    
    /// Return the mesh
    const Mesh & mesh() const
    {
      return m_vsxgrad.mesh();
    } 

    /// Return the polynomial degree
    const size_t & degree() const
    {
      return m_vsxgrad.degree();
    }

    /// Interpolator of a continuous function
    Eigen::VectorXd interpolate(
          const ScalarFunctionType & q, ///< The function to interpolate
          const int doe_cell = -1, ///< The optional degre of cell quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_face = -1, ///< The optional degre of face quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_edge = -1 ///< The optional degre of edge quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          ) const;

    /// Return edge bases for the edge of index iE
    inline const DDRCore::EdgeBases & edgeBases(size_t iE) const
    {
      return m_vsxgrad.sXgrad().edgeBases(iE);
    }

    /// Return edge bases for edge E
    inline const DDRCore::EdgeBases & edgeBases(const Edge & E) const
    {
      return m_vsxgrad.sXgrad().edgeBases(E.global_index());
    }

    /// Return face bases for the face of index iF
    inline const DDRCore::FaceBases & faceBases(size_t iF) const
    {
      return m_vsxgrad.sXgrad().faceBases(iF);
    }

    /// Return face bases for Face F
    inline const DDRCore::FaceBases & faceBases(const Face & F) const
    {
      return m_vsxgrad.sXgrad().faceBases(F.global_index());
    }

    /// Return cell bases for the cell of index iT
    inline const DDRCore::CellBases & cellBases(size_t iT) const
    {
      return m_vsxgrad.sXgrad().cellBases(iT);
    }

    /// Return cell bases for cell T
    inline const DDRCore::CellBases & cellBases(const Cell & C) const
    {
      return m_vsxgrad.sXgrad().cellBases(C.global_index());
    }
    
    
    
    private:
    const VSXGrad & m_vsxgrad;
    const ConnectedInformationGenerator & m_info_generator;
    bool m_use_threads;
    // std::ostream & m_output;

      
  };

} // namespace HArDCore3D

#endif
