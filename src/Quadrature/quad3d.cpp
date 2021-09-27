
/*

  1) Classes implementing some (most used) quadrature formulas on tetrahedra
  
  2) This implementation does not require any specific mesh structures, but
  only arrays of coordinates of triangles vertices given as double[].

*/

#include <quad3d.hpp>
#include <max_degrees_quadratures.hpp>

/*
  
  cq0[0..nqn-1],  barycentric coordinates of tetrahedron's vertices V0,V1,V2,V3
  cq1[0..nqn-1],
  cq2[0..nqn-1],
  cq3[0..nqn-1],

  cqw[0..nqn-1],  weights
  
*/


using namespace HArDCore3D;

QuadRuleTetra::QuadRuleTetra(size_t doe, bool warn)
      : rule(required_rule(doe)),
        nqn(tetra_unit_size(rule)),
        xV(4), yV(4), zV(4) {
    if (doe > MAX_DOE_CELL && warn)
      std::cerr << "Warning: quadrature rule of degree " << doe <<
                " requested, but the maximum available is degree " << MAX_DOE_CELL << std::endl;
    init();
  }

QuadRuleTetra::QuadRuleTetra(size_t _rule)
      : rule(_rule),
        nqn(tetra_unit_size(rule)),
        xV(4), yV(4), zV(4) {
    init();
  }


void QuadRuleTetra::init() {
    assert(rule < num_rules);
    const size_t size = nqn;

    cq0.resize(size);
    cq1.resize(size);
    cq2.resize(size);
    cq3.resize(size);
    cwq.resize(size);

    tetra_unit_set(rule, size, cq0.data(), cq1.data(), cq2.data(), cwq.data());

    for (size_t iq = 0; iq < size; ++iq) {
      cq3[iq] = 1. - cq0[iq] - cq1[iq] - cq2[iq];
    }
  }


size_t QuadRuleTetra::required_rule(size_t doe) const {
//  static const std::array<size_t, num_rules> exactness = {1, 1, 2, 2, 3, 3, 4, 8};
  static const std::array<size_t, num_rules> exactness = {1, 1, 2, 3, 4, 5, 6, 7, 8, 9 , 10 , 11, 12, 13, 14};
  return std::min(size_t(std::lower_bound(std::begin(exactness), std::end(exactness), doe) - std::begin(exactness)),
                  num_rules - 1);
}


void QuadRuleTetra::setup(double (&_xV)[4], double (&_yV)[4], double (&_zV)[4]) {
  std::copy(std::begin(_xV), std::end(_xV), std::begin(xV));
  std::copy(std::begin(_yV), std::end(_yV), std::begin(yV));
  std::copy(std::begin(_zV), std::end(_zV), std::begin(zV));
  vol_T = tetra_volume(xV.data(), yV.data(), zV.data());
}

void QuadRuleTetra::get_quadrule(int iq, double &_xq, double &_yq, double &_zq, double &_wq) const {
  _xq = xq(iq);
  _yq = yq(iq);
  _zq = zq(iq);
  _wq = wq(iq);
}



