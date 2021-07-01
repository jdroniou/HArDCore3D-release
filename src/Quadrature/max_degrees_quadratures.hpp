// Maximal degrees of quadratures for cell, faces and edges
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
*
*  This library was developed around HHO methods, although some parts of it have a more
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



#ifndef MAX_DEGREES_QUADRATURES_HPP
#define MAX_DEGREES_QUADRATURES_HPP


namespace HArDCore3D {

/*!
*	@addtogroup Quadratures
* @{
*/

  static constexpr size_t MAX_DOE_CELL = 14; ///< Maximum degree of the cell quadrature rules
  static constexpr size_t MAX_DOE_FACE = 20; ///< Maximum degree of the face quadrature rules
  static constexpr size_t MAX_DOE_EDGE = 20; ///< Maximum degree of the edge quadrature rules

//@}
}

#endif /* MAX_DEGREES_QUADRATURES_HPP */
