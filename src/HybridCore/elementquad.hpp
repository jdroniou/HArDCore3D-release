// Class to create and store values of cell and face basis functions on quadrature points
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
 *
 *	This library was developed around HHO methods, although some parts of it have a more
 * general purpose. If you use this code or part of it in a scientific publication, 
 * please mention the following book as a reference for the underlying principles
 * of HHO schemes:
 *
 * The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications.
 * D. A. Di Pietro and J. Droniou. 2019, 516p. 
 * url: https://hal.archives-ouvertes.fr/hal-02151813.
 *
 */

#ifndef ELEMENTQUAD_HPP
#define ELEMENTQUAD_HPP

//#include <cassert>
//#include <cmath>

//#include <functional>
//#include <memory>
#include <Eigen/Dense>
#include <hybridcore.hpp>
//#include <mesh.hpp>
#include <quadraturerule.hpp>


/*!	
 * @defgroup HybridCore 
 * @brief Classes providing cell and face quadrature rules, and values of basis functions at the nodes
 */

namespace HArDCore3D {

  /*!
   *	\addtogroup HybridCore
   * @{
   */
  // ----------------------------------------------------------------------------
  //                            Class definition
  // ----------------------------------------------------------------------------

  /** The ElementQuad class creates cell and face quadrature rules, and vectors of values of basis functions and
   *  gradients at these points
   **/

  class ElementQuad {

  public:
    ///@brief Class constructor: loads the quadrature rules and values of basis functions/gradients at these points
    ElementQuad(
		const HybridCore& hho, ///< A reference to the hybridcore instance
		const size_t iT, ///< Number of cell
		const size_t doeT, ///< The degree of exactness for cell quadratures 
		const size_t doeF ///< The degree of exactness of face quadratures
		); 

    inline QuadratureRule get_quadT() const;	///< Returns quadrature rules in cell
    inline std::vector<Eigen::ArrayXd> get_phiT_quadT() const;	///< Returns values of cell basis functions at cell quadrature rules in cell
    inline std::vector<Eigen::ArrayXXd> get_dphiT_quadT() const;	///< Returns values of gradients of cell basis functions at cell quadrature rules in cell

    inline QuadratureRule get_quadF(size_t ilF) const;	///< Returns quadrature rules on face with local number ilF
    inline std::vector<Eigen::ArrayXd> get_phiT_quadF(size_t ilF) const;	///< Returns values of cell basis functions at cell quadrature rules on face with local number ilF
    inline std::vector<Eigen::ArrayXd> get_phiF_quadF(size_t ilF) const;	///< Returns values of faces basis functions at face quadrature rules on face with local number ilF
    inline std::vector<Eigen::ArrayXXd> get_dphiT_quadF(size_t ilF) const;	///< Returns values of gradients of cell basis functions at cell quadrature rules on face with local number ilF

    /// Builds on the fly the values of cell basis functions at cell quadrature nodes. The vector basis is obtained by tensorization of the scalar one: (phi_0,0,0), (phi_1,0,0), ..., (phi_N,0,0), (0,phi_0,0), (0,phi_1,0) ... (0,phi_N,0), (0,0,phi_0) ... (0,0,phi_N).
    std::vector<Eigen::ArrayXXd> get_vec_phiT_quadT(
      size_t degree   /// maximal degree of basis functions that is required
    ) const; ///< @returns vec_phiT_quadT such that, for r=0,..,dim-1 and i=0,..,dim_Pcell(degree)-1, vec_phiT_quadT[r*dim_Pcell(degree)-1 + i] has, on its row r, the values of the i-th scalar basis function at the quadrature nodes, and 0 on its other rows. 

    /// Builds on the fly the values of cell basis functions at face quadrature nodes. The vector basis is obtained by tensorization of the scalar one: (phi_0,0,0), (phi_1,0,0), ..., (phi_N,0,0), (0,phi_0,0), (0,phi_1,0) ... (0,phi_N,0), (0,0,phi_0) ... (0,0,phi_N). 
    std::vector<Eigen::ArrayXXd> get_vec_phiT_quadF(
        size_t ilF,   /// local number of face
        size_t degree   /// maximum degree of basis functions required
    ) const; ///< @returns vec_phiT_quadT such that, for r=0,..,dim-1 and i=0,..,dim_Pcell(degree)-1, vec_phiT_quadT[r*dim_Pcell(degree)-1 + i] has, on its row r, the values of the i-th scalar basis function at the quadrature nodes, and 0 on its other rows. 

    /// Builds on the fly the values of face basis functions at face quadrature nodes. The vector basis is obtained by tensorization of the scalar one: (phi_0,0,0), (phi_1,0,0), ..., (phi_N,0,0), (0,phi_0,0), (0,phi_1,0) ... (0,phi_N,0), (0,0,phi_0) ... (0,0,phi_N). 
    std::vector<Eigen::ArrayXXd> get_vec_phiF_quadF(
        size_t ilF,   /// local number of face
        size_t degree   /// required degree of basis function
     ) const; ///< @returns vec_phiT_quadT such that, for r=0,..,dim-1 and i=0,..,dim_Pcell(degree)-1, vec_phiT_quadT[r*dim_Pcell(degree)-1 + i] has, on its row r, the values of the i-th scalar basis function at the quadrature nodes, and 0 on its other rows. 
  
  private:
    /// Mesh, cell, degrees
    const HybridCore& _hho;  // reference to the hybridcore instance
    const size_t _iT; // cell number
    const size_t _doeT; // degree of exactness of cell quadrature rules
    const size_t _doeF; // degree of exactness of face quadrature rules

    /// Quadrature and values of basis functions in cells
    QuadratureRule _quadT;
    std::vector<Eigen::ArrayXd> _phiT_quadT;
    std::vector<Eigen::ArrayXXd> _dphiT_quadT;

    /// Quadratures and values of basis functions on faces
    std::vector<QuadratureRule> _quadF;
    std::vector<std::vector<Eigen::ArrayXd>> _phiT_quadF;
    std::vector<std::vector<Eigen::ArrayXd>> _phiF_quadF;
    std::vector<std::vector<Eigen::ArrayXXd>> _dphiT_quadF;

  };



  // --------------------------------------------------------------------------------------------------
  // ------- Functions that return class elements


  QuadratureRule ElementQuad::get_quadT() const { return _quadT; }
  std::vector<Eigen::ArrayXd> ElementQuad::get_phiT_quadT() const { return _phiT_quadT; }
  std::vector<Eigen::ArrayXXd> ElementQuad::get_dphiT_quadT() const { return _dphiT_quadT; }
  QuadratureRule ElementQuad::get_quadF(size_t ilF) const { return _quadF[ilF]; }
  std::vector<Eigen::ArrayXd> ElementQuad::get_phiT_quadF(size_t ilF) const { return _phiT_quadF[ilF]; }
  std::vector<Eigen::ArrayXd> ElementQuad::get_phiF_quadF(size_t ilF) const { return _phiF_quadF[ilF]; }
  std::vector<Eigen::ArrayXXd> ElementQuad::get_dphiT_quadF(size_t ilF) const { return _dphiT_quadF[ilF]; }

  //@}

}  // end of namespace HArDCore3D

#endif /* ELEMENTQUAD_HPP */
