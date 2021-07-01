#ifndef _QUADRULE_FACE_BASE_HPP
#define _QUADRULE_FACE_BASE_HPP

#include <cassert>
#include <cmath>
#include <memory>
#include <triangle_dunavant_rule.hpp>

namespace HArDCore3D {

/*!
*	@addtogroup Quadratures
* @{
*/

class QuadRuleTriangle{

public:
    /**
    * @brief Default constructor
    *
    * @param doe degrees of exactness (e.g. how many points for approximating
    *integral
    * @param warn
    */
  QuadRuleTriangle( size_t _doe, bool warn );
  ~QuadRuleTriangle();
   double xq( size_t i );
  double yq( size_t i );
   double zq( size_t i );
  double wq( size_t i );
   size_t nq();
  void setup( double xV[], double yV[], double zV[] );

protected:
   size_t _npts;
  double*  qsx ;
  double*  qsy ;
  double*  qsz ;
  double* qswg ;

    double* _xy;
    double* _w;

  double area( double xV[], double yV[], double zV[] ) const;

};

} // end of namespace HArDCore3D

#endif // _QUADRULE_FACE_BASE_HPP
