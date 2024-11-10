// Core data structure for spaces and operators on manifolds.
//
// Author: Marien Hanot (marien-lorenzo.hanot@umontpellier.fr)
//


#ifndef EXTERIOR_HPP
#define EXTERIOR_HPP

#include "exterior_dimension.hpp"

// Std utilities
#include <vector>
#include <array>
#include <unordered_map>
#include <algorithm> 
#include <cstdlib>
#include <cassert>

// Linear algebra utilities
#include <Eigen/Dense> 
#include <unsupported/Eigen/KroneckerProduct> // Used to couple the action on polynomial and on the exterior algebra

namespace Manicore {
   /**
     The methods in this file are meant to compute the action of everything that is independant of the atlas.

     The most useful are:
      Compute_pullback : computes the action of a pullback to the exterior algebra
      Monomial_full : gives a mapping between index and monomial powers
      Koszul_full : gives the matrix of the Koszul operator
      Diff_full : gives the matrix of the Diff operator

     Call Initialize_exterior_module<d>::init(int r) to initialize every module on dimension d
   */
  

  ///------------------------------------------------------------------------------------------------------------------------------
  // Exterior algebra
  ///------------------------------------------------------------------------------------------------------------------------------
  /** Contains the part of methods dealing with the exterior algebra basis
    */

  /// Class to handle the exterior algebra basis
  template<size_t l, size_t d>
  class ExteriorBasis 
  {
    public:
      static void init() noexcept {
        static_assert(0 <= l && l <= d,"Error: Tried to generate the exterior algebra basis outside the range [0,d]");
        if (initialized == 1) return; // initialize only the first time
        initialized = 1;
        size_t acc = 0;
        std::array<size_t, l> cbasis;
        for (size_t i = 0; i < l; ++i) {
          cbasis[i] = i; // fill to first element
        }
        _basis.try_emplace(acc,cbasis);
        while (_find_next_tuple(cbasis)) {
          ++acc;
          _basis.try_emplace(acc,cbasis);
        };
      };

      static const std::array<size_t,l> & expand_basis(size_t i) 
      {
        assert(initialized==1 && "ExteriorBasis was not initialized");
        return _basis.at(i);
      }

      static const size_t index_from_tuple(const std::array<size_t,l> & tuple) {
        auto found = std::find_if(_basis.begin(), _basis.end(), [&tuple](const auto& p)
                                                                {return p.second == tuple;});
        if (found != _basis.end()) {
          return found->first;
        } else {
          throw; // tuple not in range
        }
      }
    
    private:
      // Unordered_map provide faster lookup time than map
      static inline std::unordered_map<size_t, std::array<size_t, l>> _basis;
      static inline int initialized = 0; // We need inline to be able to define it inside the class
      // Helper, find the next element of the basis, return false if this is the last element
      static bool _find_next_tuple(std::array<size_t, l> &basis) noexcept {
        for (size_t i = 0; i < l; ++i) {
          size_t nval = basis[l - i - 1] + 1;
          if (nval < d - i) { // value valid, we must reset everything after
            basis[l - i - 1] = nval;
            for (size_t j = 1; j <= i ; ++j){
              basis[l - i - 1 + j] = nval + j;
            }
            return true;
          } else { // already using the last element available
            continue;
          }
        }
        return false;
      };
  };

  /// Compute the action of Kozsul and Diff on the exterior algebra.
  /** The action is returned as a list of matrix between the exterior algebra basis
   To get the full action, one must take the kronecker product between the action on the polynomial and the action on the exterior algebra
   In the case of Koszul, the action is the multiply by x_i, and in the case of Diff, it is the differentiated by x_i. */
  template<size_t l, size_t d>
  class Koszul_exterior {
    public:
      typedef Eigen::Matrix<double, Dimension::ExtDim(l-1,d), Dimension::ExtDim(l,d)> ExtAlgMatType;

      static const ExtAlgMatType & get_transform(size_t i) {
        return _transforms.at(i);
      }

      static void init() noexcept {
        static_assert(0 < l && l <= d,"Error: Tried to generate Koszul basis outside the range [1,d]");
        if (initialized == 1) return;
        initialized = 1;
        ExteriorBasis<l-1,d>::init(); // ensure that the exterior basis is initialized
        ExteriorBasis<l,d>::init(); // ensure that the exterior basis is initialized
        for (size_t i = 0; i < d; ++i) {
          _transforms[i].setZero();
        }
        if (l == 1) { // special case, ext_basis = \emptyset
          for (size_t i = 0; i < d; ++i) {
            _transforms[i](0,i) = 1;
          }
          return;
        } 
        for (size_t i = 0; i < Dimension::ExtDim(l-1,d); ++i) {
          int sign = 1;
          size_t try_pos = 0;
          const std::array<size_t, l-1> & cbasis = ExteriorBasis<l-1,d>::expand_basis(i);
          std::array<size_t,l> basis_cp;
          std::copy(cbasis.cbegin(),cbasis.cend(),basis_cp.begin()+1);
          // basis_cp contains a copy of cbasis with one more slot
          for (size_t j = 0; j < d; ++j) {
            if (try_pos == l-1 || j < cbasis[try_pos]) { // insert here
              basis_cp[try_pos] = j;
              _transforms[j](i,ExteriorBasis<l,d>::index_from_tuple(basis_cp)) = sign;
            } else { // j == cbasis[try_pos]
              basis_cp[try_pos] = basis_cp[try_pos+1]; 
              ++try_pos;
              sign = -sign;
              continue;
            }
          }
        }
      }; // end init()

    private:
      static inline int initialized = 0;
      static inline std::array<ExtAlgMatType,d> _transforms;
  };

  template<size_t l, size_t d>
  class Diff_exterior {
    public:
      typedef Eigen::Matrix<double, Dimension::ExtDim(l+1,d), Dimension::ExtDim(l,d)> ExtAlgMatType;

      static const ExtAlgMatType & get_transform(size_t i) {
        return _transforms.at(i);
      }

      static void init() noexcept {
        static_assert(l < d,"Error: Tried to generate diff basis outside the range [0,d-1]");
        if (initialized == 1) return;
        initialized = 1;
        Koszul_exterior<l+1,d>::init();
        for (size_t i = 0; i < d; ++i) {
          _transforms[i] = Koszul_exterior<l+1,d>::get_transform(i).transpose().eval();
        }
      };
    private:
      static inline int initialized = 0;
      static inline std::array<ExtAlgMatType,d> _transforms;
  };

  ///------------------------------------------------------------------------------------------------------------------------------
  // Pullback helpers
  ///------------------------------------------------------------------------------------------------------------------------------

  /// Generic determinant computation
  /** The first two arguments should be the list of indexes to use, and the last the matrix
    This function returns the determinant of the partial matrix
    */
  template<typename V, typename Derived>
  double Compute_partial_det(const V& a1, const V& a2, const Eigen::MatrixBase<Derived>& A) {
    constexpr size_t N = std::tuple_size<V>::value;
    if constexpr (N==0) {
      return 0.;
    } else if constexpr (N == 1) {
      return A(a1[0],a2[0]);
    } else {
      double sign = 1.;
      double sum = 0.;
      std::array<typename V::value_type,N-1> b1,b2;
      std::copy(a1.cbegin()+1,a1.cend(),b1.begin());
      std::copy(a2.cbegin()+1,a2.cend(),b2.begin());
      for (size_t i = 0; i < N-1;++i){
        sum += sign*A(a1[i],a2[0])*Compute_partial_det(b1,b2,A);
        sign *= -1.;
        b1[i] = a1[i];
      }
      sum += sign*A(a1[N-1],a2[0])*Compute_partial_det(b1,b2,A);
      return sum;
    }
  }

  /// Generic pullback computation
  /// The matrix A go from the space 1 to the space 2
  /// Specialized for some case, try to avoid the generic definition
  template<size_t l, size_t d1, size_t d2>
  struct Compute_pullback {
    template<typename Derived> 
    static Eigen::Matrix<double, Dimension::ExtDim(l,d1), Dimension::ExtDim(l,d2)> compute(Eigen::MatrixBase<Derived> const & A) {
      Eigen::Matrix<double, Dimension::ExtDim(l,d1), Dimension::ExtDim(l,d2)> rv;
      for (size_t i = 0; i < Dimension::ExtDim(l,d1); ++i) {
        for (size_t j = 0; j < Dimension::ExtDim(l,d2); ++j) {
          rv(i,j) = Compute_partial_det(ExteriorBasis<l,d2>::expand_basis(j),ExteriorBasis<l,d1>::expand_basis(i),A);
        }
      }
      return rv;
    }
  };

  template<size_t d1, size_t d2> 
  struct Compute_pullback<0,d1,d2> {
    template<typename Derived> 
    static Eigen::Matrix<double,1,1> compute(Eigen::MatrixBase<Derived> const & A) {
      return Eigen::Matrix<double,1,1>{1.};
    }
  };

  template<size_t d1, size_t d2> 
  struct Compute_pullback<1,d1,d2> {
    template<typename Derived> 
    static Eigen::Matrix<double,d1,d2> compute(Eigen::MatrixBase<Derived> const & A) {
      return A.transpose();
    }
  };

  template<size_t d> 
  struct Compute_pullback<d,d,d> {
    template<typename Derived> 
    static Eigen::Matrix<double,1,1> compute(Eigen::MatrixBase<Derived> const & A) {
      return Eigen::Matrix<double,1,1>{A.determinant()};
    }
  };

  template<> 
  struct Compute_pullback<1,1,1> {
    template<typename Derived> 
    static Eigen::Matrix<double,1,1> compute(Eigen::MatrixBase<Derived> const & A) {
      return Eigen::Matrix<double,1,1>{A(0,0)};
    }
  };

  template<>
  struct Compute_pullback<2,2,3> {
    template<typename Derived> 
    static Eigen::Matrix<double,1,3> compute(Eigen::MatrixBase<Derived> const & A) {
      return Eigen::Matrix<double,1,3>{{A(0,0)*A(1,1) - A(0,1)*A(1,0), A(0,0)*A(2,1) - A(0,1)*A(2,0), A(1,0)*A(2,1) - A(1,1)*A(2,0)}};
    }
  };
  template<>
  struct Compute_pullback<2,3,2> {
    template<typename Derived> 
    static Eigen::Matrix<double,3,1> compute(Eigen::MatrixBase<Derived> const & A) {
      return Eigen::Matrix<double,3,1>{A(0,0)*A(1,1) - A(0,1)*A(1,0), A(0,0)*A(1,2) - A(0,2)*A(1,0), A(0,1)*A(1,2) - A(0,2)*A(1,1)};
    }
  };
  template<>
  struct Compute_pullback<2,3,3> {
    template<typename Derived> 
    static Eigen::Matrix<double,3,3> compute(Eigen::MatrixBase<Derived> const & A) {
      return Eigen::Matrix<double,3,3>{
        {A(0,0)*A(1,1) - A(0,1)*A(1,0), A(0,0)*A(2,1) - A(0,1)*A(2,0), A(1,0)*A(2,1) - A(1,1)*A(2,0)},
        {A(0,0)*A(1,2) - A(0,2)*A(1,0), A(0,0)*A(2,2) - A(0,2)*A(2,0), A(1,0)*A(2,2) - A(1,2)*A(2,0)},
        {A(0,1)*A(1,2) - A(0,2)*A(1,1), A(0,1)*A(2,2) - A(0,2)*A(2,1), A(1,1)*A(2,2) - A(1,2)*A(2,1)}
      };
    }
  };

  ///------------------------------------------------------------------------------------------------------------------------------
  // Polynomial algebra
  ///------------------------------------------------------------------------------------------------------------------------------
  /** Contains the part of methods dealing with the polynomial basis
    */
  /// Generate a basis of monomial powers of degree r
  template<size_t d> 
  struct Monomial_powers {
    static std::vector<std::array<size_t, d>> const & homogeneous(const int r) {
      assert(r <= _init_deg && "Error: Monomial_powers must be initialized first");
      return _powers[r];
    }

    static std::vector<std::array<size_t,d>> const & complete(const int r) {
      assert(r <= _init_deg && "Error: Monomial_powers must be initialized first");
      return _powers_full[r];
    }

    static void init(const int r) {
      static_assert(d > 0,"Error: Tried to construct a monomial basis in dimension 0");
      assert(r <= _max_deg && "Maximum degree reached, increase _max_deg if needed");
      if (r >_init_deg) {
        for (int h_r = _init_deg + 1; h_r <= r; ++h_r) { // init homogeneous
          std::vector<std::array<size_t,d>> powers;
          std::array<size_t,d> current;
          inner_loop(0,h_r,current,powers);
          _powers[h_r] = powers;
          if (h_r == 0) {
            _powers_full[0] = powers;
          } else {
            _powers_full[h_r] = _powers_full[h_r-1];
            for (size_t j = 0; j < powers.size(); ++j) {
              _powers_full[h_r].emplace_back(powers[j]);
            }
          }
        }
        _init_deg = r;
      }
    }

    private:
      static constexpr int _max_deg = 20;
      static inline int _init_deg = -1;
      static inline std::array<std::vector<std::array<size_t,d>>,_max_deg+1> _powers;
      static inline std::array<std::vector<std::array<size_t,d>>,_max_deg+1> _powers_full;

      static void inner_loop(size_t cindex, int max, std::array<size_t,d> & current, std::vector<std::array<size_t,d>> & powers) {
        if (cindex == d-1) {
          current[cindex] = max; // Last direction to enforce the degree
          powers.emplace_back(current);
        } else {
          for(int i = 0; i <= max; ++i) {
            current[cindex] = i;
            inner_loop(cindex+1,max-i,current,powers);
          }
        }
      }
  };

  /// Generate the matrices for the Koszul operator on homogeneous monomial
  template<size_t d, size_t index>
  struct Koszul_homogeneous_mat {
    static Eigen::MatrixXd get (const int r) {
      static_assert(index < d,"Error: Tried to take the koszul operator on a direction outside the dimension");
      Eigen::MatrixXd M = Eigen::MatrixXd::Zero(Dimension::HDim(r+1,d), Dimension::HDim(r,d));
      for (size_t i = 0; i < Dimension::HDim(r,d); ++i) {
        std::array<size_t, d> current = Monomial_powers<d>::homogeneous(r)[i];
        current.at(index) += 1;
        size_t val = std::find(Monomial_powers<d>::homogeneous(r+1).cbegin(), 
                               Monomial_powers<d>::homogeneous(r+1).cend(),current) 
                    - Monomial_powers<d>::homogeneous(r+1).cbegin();
        M(val,i) = 1.;
      }
      return M;
    }
  };
  /// Generate the matrices for the Differential operator on homogeneous monomial
  template<size_t d, size_t index>
  struct Diff_homogeneous_mat {
    static Eigen::MatrixXd get (const int r) {
      static_assert(index < d,"Error: Tried to take the differential operator on a direction outside the dimension");
      assert(r > 0 && "Error: Cannot generate a matrix for the differential on P_0");
      Eigen::MatrixXd M = Eigen::MatrixXd::Zero(Dimension::HDim(r-1,d), Dimension::HDim(r,d)); 
      for (size_t i = 0; i < Dimension::HDim(r,d); ++i) {
        std::array<size_t, d> current = Monomial_powers<d>::homogeneous(r)[i];
        size_t cval = current.at(index);
        if (cval == 0) continue; // Zero on that col
        current.at(index) -= 1;
        size_t val = std::find(Monomial_powers<d>::homogeneous(r-1).cbegin(), 
                               Monomial_powers<d>::homogeneous(r-1).cend(),current) 
                    - Monomial_powers<d>::homogeneous(r-1).cbegin();
        M(val,i) = cval;
      }
      return M;
    }
  };

  ///------------------------------------------------------------------------------------------------------------------------------
  // Coupling polynomial and exterior
  ///------------------------------------------------------------------------------------------------------------------------------
  /** Couple the part on the exterior algebra and the part on polynomials
    */
  /// Koszul operator from $P_r\Lambda^l(R^d)$ to $P_{r+1}\Lambda^{l-1}(R^d)$.
  template<size_t l, size_t d>
  struct Koszul_full {
    static Eigen::MatrixXd get (const int r) {
      if constexpr (l==0 || l > d) {
        return Eigen::MatrixXd(0,0);
      }
      Eigen::MatrixXd M(Dimension::PLDim(r+1,l-1,d),Dimension::PLDim(r,l,d));
      M.setZero();
      if (Dimension::PLDim(r+1,l-1,d) > 0 && Dimension::PLDim(r,l,d) > 0) {
        Eigen::MatrixXd scalar_part(Dimension::PolyDim(r+1,d),Dimension::PolyDim(r,d));
        init_loop_for<0>(r,scalar_part,M);
      }
      return M;
    }

    private:
      template<size_t i, typename T> static void init_loop_for(int r, T & scalar_part,T & M) {
        if constexpr(i < d) {
          scalar_part.setZero();
          int offset_l = Dimension::HDim(0,d);
          int offset_c = 0;
          for (int s = 0; s <= r; ++s) {
            int inc_l = Dimension::HDim(s+1,d);
            int inc_c = Dimension::HDim(s,d);
            scalar_part.block(offset_l,offset_c,inc_l,inc_c) = Koszul_homogeneous_mat<d,i>::get(s);
            offset_l += inc_l;
            offset_c += inc_c;
          }
          M += Eigen::KroneckerProduct(Koszul_exterior<l,d>::get_transform(i),scalar_part);
          init_loop_for<i+1>(r,scalar_part,M);
        }
      }
  };

  /// Differential operator from $P_r\Lambda^l(R^d)$ to $P_{r-1}\Lambda^{l+1}(R^d)$.
  template<size_t l, size_t d>
  struct Diff_full {
    static Eigen::MatrixXd get (const int r) {
      if (l >= d) {
        return Eigen::MatrixXd(0,0);
      }
      Eigen::MatrixXd M(Dimension::PLDim(r-1,l+1,d),Dimension::PLDim(r,l,d));
      M.setZero();
      if (Dimension::PLDim(r-1,l+1,d) > 0 || Dimension::PLDim(r,l,d) > 0) {
        Eigen::MatrixXd scalar_part(Dimension::PolyDim(r-1,d),Dimension::PolyDim(r,d));
        init_loop_for<0>(r,scalar_part,M);
      }
      return M;
    };

    static Eigen::MatrixXd get_as_degr (const int r) {
      auto inj = Eigen::KroneckerProduct(Eigen::Matrix<double,Dimension::ExtDim(l+1,d),Dimension::ExtDim(l+1,d)>::Identity(),
                                 Eigen::MatrixXd::Identity(Dimension::PolyDim(r,d),Dimension::PolyDim(r-1,d)));
      return inj*get(r);
    };

    private:
      template<size_t i,typename T> static void init_loop_for(int r,T &scalar_part,T &M) {
        if constexpr (i < d) {
          scalar_part.setZero();
          int offset_l = 0;
          int offset_c = Dimension::HDim(0,d);
          for (int s = 0; s < r; ++s) {
            int inc_l = Dimension::HDim(s,d);
            int inc_c = Dimension::HDim(s+1,d);
            scalar_part.block(offset_l,offset_c,inc_l,inc_c) = Diff_homogeneous_mat<d,i>::get(s+1);
            offset_l += inc_l;
            offset_c += inc_c;
          }
          M += Eigen::KroneckerProduct(Diff_exterior<l,d>::get_transform(i),scalar_part);
          init_loop_for<i+1>(r,scalar_part,M);
        }
      }
  };

  /// Initialize every class related to the polynomial degree r
  template<size_t d>
  struct Initialize_exterior_module{
    static void init(int r) 
    {
      init_loop_for<0,1>(r+1);
    }

    private:
      template<size_t l,size_t k> static void init_loop_for(int r) {
        if constexpr(k <= d) {
          if constexpr(l < k) {
            Diff_exterior<l,k>::init();
            init_loop_for<l+1,k>(r);
          } else {
            Monomial_powers<k>::init(r);
            init_loop_for<0,k+1>(r);
          }
        }
      }
  };
  // As a general rule, global classes will call init() but not the local classes (the classes which operates differently on each element)

} // End namespace

#endif
