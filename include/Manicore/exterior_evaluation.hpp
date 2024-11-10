#ifndef EXTERIOR_EVALUATION_HPP
#define EXTERIOR_EVALUATION_HPP

#include "exterior.hpp"

namespace Manicore {
  ///------------------------------------------------------------------------------------------------------------------------------
  // Initial basis evaluation
  ///------------------------------------------------------------------------------------------------------------------------------
  /** Contains methods to evaluate $P_r\Lambda^l(\Real^d)$
    */


  template<size_t d_outer, size_t d>
  class Monomial_scalar_basis_linear_ID {
    public:
      Monomial_scalar_basis_linear_ID(double A, Eigen::Vector<double,d> const &C, Eigen::Matrix<double, d, d_outer> const &T,int r):
        _r(r),_A(A),_C(C),_T(T) {static_assert(d < d_outer && "Object dimension should be less than space dimension");}
      Monomial_scalar_basis_linear_ID(): _r(-1) {;}

      double evaluate(const Eigen::Vector<double, d_outer> &x, size_t i) const {
        assert((_r >= 0) && "Basis not initialized");
        assert((i < Dimension::PolyDim(_r,d)) && "Index out of range");
        auto const & power = Monomial_powers<d>::complete(_r)[i];
        Eigen::Vector<double, d> y = _A*_T*x + _C;
        double rv = std::pow(y(0),power[0]);
        for (size_t i = 1; i < d; ++i) {
          rv *= std::pow(y(i),power[i]);
        }
        return rv;
      }

   private:
    int _r;
    double _A;
    Eigen::Vector<double, d> _C;
    Eigen::Matrix<double, d, d_outer> _T;
  };
  
  template<size_t d>
  class Monomial_scalar_basis_linear_ID<d,d> {
    public:
      Monomial_scalar_basis_linear_ID(double A, Eigen::Vector<double,d> const &C,Eigen::Matrix<double, d, d> const &T,int r):
        _r(r),_A(A),_C(C) {;}
      Monomial_scalar_basis_linear_ID(): _r(-1) {;}

      double evaluate(const Eigen::Vector<double, d> &x, size_t i) const {
        assert((_r >= 0) && "Basis not initialized");
        assert((i < Dimension::PolyDim(_r,d)) && "Index out of range");
        auto const & power = Monomial_powers<d>::complete(_r)[i];
        Eigen::Vector<double, d> y = _A*x + _C;
        double rv = std::pow(y(0),power[0]);
        for (size_t i = 1; i < d; ++i) {
          rv *= std::pow(y(i),power[i]);
        } 
        return rv;
      }

   private:
    int _r;
    double _A;
    Eigen::Vector<double, d> _C;
  };


}

#endif

