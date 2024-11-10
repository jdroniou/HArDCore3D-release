#ifndef EXTERIOR_DIMENSION_HPP
#define EXTERIOR_DIMENSION_HPP

#include <cstddef>

namespace Manicore::Dimension {
  constexpr size_t binom(size_t n, size_t k) noexcept
  {
    return
      (        k> n  )? 0 :          // out of range
      (k==0 || k==n  )? 1 :          // edge
      (k==1 || k==n-1)? n :          // first
      (     k+k < n  )?              // recursive:
      (binom(n-1,k-1) * n)/k :       //  path to k=1   is faster
      (binom(n-1,k) * n)/(n-k);      //  path to k=n-1 is faster
  }

  /// Dimension of the exterior algebra
  constexpr size_t ExtDim(size_t l,size_t d) noexcept
  {
    return binom(d,l);
  }

  /// Dimension of P^r(\Real^d)
  constexpr size_t PolyDim(int r,size_t d) noexcept
  {
    return (r >= 0) ? binom(r+d,d):0;
  }

  /// Dimension of homogeneous polynomials
  constexpr size_t HDim(int r, size_t d) noexcept
  {
    return PolyDim(r,d) - PolyDim(r-1,d);
  }

  /// Dimension of $P_r\Lambda^l(R^d)$
  constexpr size_t PLDim(int r, size_t l, size_t d) noexcept
  {
    return ExtDim(l,d)*PolyDim(r,d);
  }

  /// Dimension of the image of Koszul on homogeneous polynomials
  constexpr size_t kHDim(int r, size_t l, size_t d) noexcept
  {
    return
      (l > d) ? 0 :
      (l == 0) ? 0 :
      (r < 0) ? 0 :
      binom(d+r,d-l)*binom(r+l-1,l-1);
  }
  /// Dimension of the image of d on homogeneous polynomials
  constexpr size_t dHDim(int r, size_t l, size_t d) noexcept
  {
    return kHDim(r-1,l+1,d);
  }
  /// Dimension of the image of Koszul on polynomials
  constexpr size_t kPLDim(int r, size_t l, size_t d) noexcept
  {
    size_t sum = 0;
    for (int i = 0; i <= r; ++i) {
      sum += kHDim(i,l,d);
    }
    return sum;
  }

  /// Dimension of the image of d on polynomials
  constexpr size_t dPLDim(int r, size_t l, size_t d) noexcept
  {
    size_t sum = 0;
    for (int i = 0; i <= r; ++i) {
      sum += dHDim(i,l,d);
    }
    return sum;
  }

  // Dimension of trimmed polynomial spaces
  constexpr size_t PLtrimmedDim(int r, size_t l, size_t d) noexcept
  {
    return (l == 0)? 
            ((r >= 0)? kPLDim(r-1,l+1,d) + 1 : 0) // We must add P_0\Lambda^0
            : kPLDim(r-1,l+1,d) + dPLDim(r,l-1,d);
  }

}

#endif

