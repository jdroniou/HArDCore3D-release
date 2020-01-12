# ifndef _QUADTD_HPP
# define _QUADTD_HPP

/*

  1) Classes implementing some (most used) quadrature formulas on tetrahedra
  
  2) This implementation does not require any specific mesh structures, but
  only arrays of coordinates of triangles vertices given as double[].

*/

#include <cassert>
#include <cmath>

#include <algorithm>
#include <array>
#include <iostream>
#include <vector>

#include "tetra_quad.hpp"

/*
  
  cq0[0..nqn-1],  barycentric coordinates of tetrahedron's vertices V0,V1,V2,V3
  cq1[0..nqn-1],
  cq2[0..nqn-1],
  cq3[0..nqn-1],

  cqw[0..nqn-1],  weights
  
*/


/*!	
*	@defgroup Quadratures 
* @brief Classes providing quadratures on edges and in cells
*/

namespace HArDCore3D {

/*!
*	@addtogroup Quadratures
* @{
*/

class QuadRuleTetra {
protected:
  static constexpr size_t max_doe = 14;
  static constexpr size_t num_rules = 15;

  const size_t rule;
  const size_t nqn;

  double vol_T;

  std::vector<double> xV, yV, zV;
  std::vector<double> cq0, cq1, cq2, cq3, cwq;

  void init();

  /// Compute the minimum rule required to achieve the desired degree of exactness.
  /** If no such rule exists, returns the highest rule. */
  size_t required_rule(size_t doe) const;

public:

  /// Create a quadrature rule with at least the given degree of exactness (if available).
  /** The smallest such quadrature rule will be chosen. If none are available, the highest
   * degree available will be used. */
  QuadRuleTetra(size_t doe, bool warn);

  QuadRuleTetra(size_t _rule = 0);

  inline size_t nq() const;

  inline double xq(int iq) const;

  inline double yq(int iq) const;

  inline double zq(int iq) const;

  inline double wq(int iq) const;


  void setup(double (&_xV)[4], double (&_yV)[4], double (&_zV)[4]);

  void get_quadrule(int iq, double &_xq, double &_yq, double &_zq, double &_wq) const;

};
// -------------------------------------------------------------------------------------------- 

inline size_t QuadRuleTetra::nq() const { return nqn; }

inline double QuadRuleTetra::xq(int iq) const { return cq0[iq] * xV[0] + cq1[iq] * xV[1] + cq2[iq] * xV[2] + cq3[iq] * xV[3]; }

inline double QuadRuleTetra::yq(int iq) const { return cq0[iq] * yV[0] + cq1[iq] * yV[1] + cq2[iq] * yV[2] + cq3[iq] * yV[3]; }

inline double QuadRuleTetra::zq(int iq) const { return cq0[iq] * zV[0] + cq1[iq] * zV[1] + cq2[iq] * zV[2] + cq3[iq] * zV[3]; }

inline double QuadRuleTetra::wq(int iq) const { return cwq[iq] * vol_T; }

/*@}*/
} // end of namespace HArDCore3D

#endif // end of _QUADTD_HPP
