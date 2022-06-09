# ifndef _GMPOLY_CELL_HPP
# define _GMPOLY_CELL_HPP

/*

  Classes to implement the quadrature-free polynomial integration rules
  
*/

#include <cassert>
#include <cmath>

#include <algorithm>
#include <array>
#include <iostream>
#include <vector>
#include <unordered_map>

#include <mesh.hpp>
#include <quadraturerule.hpp>
#include <basis.hpp>
#include <typeinfo>

namespace HArDCore3D {

/*!
*	@addtogroup Quadratures
* @{
*/


//----------------------------------------------------------------------
//----------------------------------------------------------------------
// INTEGRALS OF MONOMIALS
//----------------------------------------------------------------------
//----------------------------------------------------------------------

/// Hash function for VectorZd type
struct VecHash
{
  size_t operator()(const VectorZd & p) const
  {
    constexpr size_t b = 100; // max degree must be less than the basis b
    return p(0) + p(1) * b + p(2) * b*b;
  }
};

/// Type for list of integrals of monomials
typedef std::unordered_map<VectorZd, double, VecHash> MonomialCellIntegralsType;
    
/// Compute all the integrals, on the edges of a cell, of this cell's monomials up to a max degree
std::vector<MonomialCellIntegralsType> IntegrateCellMonomials_onEdges
          (const Cell & T,      ///< Cell
          const size_t maxdeg           ///< Maximal total degree
          );

/// Compute all the integrals, on the faces of a cell, of this cell's monomials up to a max degree
std::vector<MonomialCellIntegralsType> IntegrateCellMonomials_onFaces
          (const Cell & T,      ///< Cell
          const size_t maxdeg,           ///< Maximal total degree
          std::vector<MonomialCellIntegralsType> & integrals_edges  ///< List of integrals of the monomials over the edges
          );

/// Compute all the integrals of a cell's monomials on the cell
MonomialCellIntegralsType IntegrateCellMonomials
          (const Cell & T,      ///< Cell
          const int maxdeg           ///< Maximal total degree
          );

/// Checks if the degree of an existing list of monomial integrals is sufficient, other re-compute and return a proper list
MonomialCellIntegralsType CheckIntegralsDegree
         (const Cell & T,              ///< Cell
           const size_t degree,         ///< Expected degree
           const MonomialCellIntegralsType & mono_int_map = {}    ///< Existing list, optional
          );

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//  TRANSFORMATIONS OF GRAM MATRICES FOR DERIVED BASES
//----------------------------------------------------------------------
//----------------------------------------------------------------------

/// Transforms a Gram Matrix from an ancestor to a family basis
template<typename BasisType>
inline Eigen::MatrixXd transformGM(const Family<BasisType> & family_basis,     ///< Family 
                               const char RC,                ///< R if transformation applied on rows (left), C if applied on columns (right)
                               const Eigen::MatrixXd & anc_GM          ///< Gram matrix of the ancestor basis
                               )
  {
    if (RC=='R'){
      return family_basis.matrix() * anc_GM;
    }else{
      return anc_GM * (family_basis.matrix()).transpose();
    }
  }

/// Transforms a Gram Matrix from an ancestor to a restricted basis
template<typename BasisType>
inline Eigen::MatrixXd transformGM(const RestrictedBasis<BasisType> & restr_basis,     ///< Restricted basis 
                               const char RC,                ///< R if transformation applied on rows (left), C if applied on columns (right)
                               const Eigen::MatrixXd & anc_GM              ///< Gram matrix of the ancestor basis
                               )
  {
    if (RC=='R'){
      return anc_GM.topRows(restr_basis.dimension());
    }else{
      return anc_GM.leftCols(restr_basis.dimension());
    }
  }

/// Transforms a Gram Matrix from an ancestor to a shifted basis
template<typename BasisType>
inline Eigen::MatrixXd transformGM(const ShiftedBasis<BasisType> & shifted_basis,     ///< Shifted basis 
                               const char RC,                ///< R if transformation applied on rows (left), C if applied on columns (right)
                               const Eigen::MatrixXd & anc_GM              ///< Gram matrix of the ancestor basis
                               )
  {
    if (RC=='R'){
      return anc_GM.bottomRows(shifted_basis.dimension());
    }else{
      return anc_GM.rightCols(shifted_basis.dimension());
    }
  }


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//  FUNCTIONS TO COMPUTE GRAM MATRICES
//----------------------------------------------------------------------
//----------------------------------------------------------------------

//
/* BASIC ONES: Monomial, tensorized, and generic template */

/// Computes the Gram Matrix of a pair of local scalar monomial bases
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                         ///< Cell to which the basis corresponds
                    const MonomialScalarBasisCell & basis1, ///< First basis
                    const MonomialScalarBasisCell & basis2, ///< Second basis
                    MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );
  
/// This overload to simplify the call to GramMatrix in case the two bases are the same
template<typename BasisType>
Eigen::MatrixXd GramMatrix(const Cell& T, const BasisType & basis, MonomialCellIntegralsType mono_int_map = {})
  {
    return GramMatrix(T, basis, basis, mono_int_map);
  };

/// Template to compute the Gram Matrix of any pair of tensorized scalar bases
template<typename BasisType1, typename BasisType2, size_t N>
Eigen::MatrixXd GramMatrix(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const TensorizedVectorFamily<BasisType1, N> & basis1, ///< First basis (rows of the Gram matrix)
                     const TensorizedVectorFamily<BasisType2, N> & basis2,  ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(basis1.dimension(), basis2.dimension());

    Eigen::MatrixXd anc_gm = GramMatrix(T, basis1.ancestor(), basis2.ancestor(), mono_int_map);
    size_t dim1 = anc_gm.rows();
    size_t dim2 = anc_gm.cols();

    for (size_t i=0; i<N; i++){
      gm.block(i*dim1, i*dim2, dim1, dim2) = anc_gm;    
    }
    return gm;
  };
  
/// Computes the Gram Matrix of a pair of RolyCompl bases
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                         ///< Cell to which the basis corresponds
                    const RolyComplBasisCell & basis1,      ///< First basis
                    const RolyComplBasisCell & basis2,      ///< Second basis
                    MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );

/// Template to compute the Gram Matrix of a RolyCompl basis and a tensorized scalar basis
template<typename BasisType1, size_t N>
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                                         ///< Cell to which the basis corresponds
                    const RolyComplBasisCell & rolycompl_basis,             ///< First basis (RolyCompl basis)
                    const TensorizedVectorFamily<BasisType1, N> & tens_family,  ///< Second basis (tensorized basis)
                    MonomialCellIntegralsType mono_int_map = {}                 ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    size_t dim1 = rolycompl_basis.dimension();
    size_t dim2 = tens_family.dimension();
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1, dim2);

    // Integrals of monomials
    size_t totaldegree = rolycompl_basis.max_degree()+tens_family.max_degree();
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

    for (size_t m=0; m<N; m++){
      gm.block(0, m*dim2/N, dim1, dim2/N) = GMRolyComplScalar(T, rolycompl_basis, tens_family.ancestor(), m, intmap);
    }
    return gm;
  };
  
/// Template to compute the Gram Matrix of a tensorized scalar basis and a RolyCompl basis
template<typename BasisType1, size_t N>
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                                         ///< Cell to which the basis corresponds
                    const TensorizedVectorFamily<BasisType1, N> & tens_family,  ///< First basis (tensorized basis)
                    const RolyComplBasisCell & rolycompl_basis,             ///< Second basis (RolyCompl basis)
                    MonomialCellIntegralsType mono_int_map = {}                 ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    return GramMatrix(T, rolycompl_basis, tens_family, mono_int_map).transpose();
  };
  
/// Computes the Gram Matrix of a pair of GolyCompl bases
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                         ///< Cell to which the basis corresponds
                    const GolyComplBasisCell & basis1,      ///< First basis
                    const GolyComplBasisCell & basis2,      ///< Second basis
                    MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );

/// Template to compute the Gram Matrix of a GolyCompl basis and a tensorized scalar basis
template<typename BasisType1, size_t N>
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                                         ///< Cell to which the basis corresponds
                    const GolyComplBasisCell & golycompl_basis,             ///< First basis (GolyCompl basis)
                    const TensorizedVectorFamily<BasisType1, N> & tens_family,  ///< Second basis (tensorized basis)
                    MonomialCellIntegralsType mono_int_map = {}                 ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    size_t dim1 = golycompl_basis.dimension();
    size_t dim2 = tens_family.dimension();
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1, dim2);

    // Integrals of monomials
    size_t totaldegree = golycompl_basis.max_degree()+tens_family.max_degree();
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

    size_t dim_l1 = golycompl_basis.dimPkmo();
    size_t dim_l2 = tens_family.ancestor().dimension();
    gm.block(0, dim_l2, dim_l1, dim_l2) = GMGolyComplScalar(T, golycompl_basis, tens_family.ancestor(), 0, intmap, 2);
    gm.block(0, 2*dim_l2, dim_l1, dim_l2) = -GMGolyComplScalar(T, golycompl_basis, tens_family.ancestor(), 0, intmap, 1);
    gm.block(dim_l1, 0, dim_l1, dim_l2) = -GMGolyComplScalar(T, golycompl_basis, tens_family.ancestor(), 1, intmap, 2);
    gm.block(dim_l1, 2*dim_l2, dim_l1, dim_l2) = GMGolyComplScalar(T, golycompl_basis, tens_family.ancestor(), 1, intmap, 0);
    gm.block(2*dim_l1, 0, dim1-2*dim_l1, dim_l2) = GMGolyComplScalar(T, golycompl_basis, tens_family.ancestor(), 2, intmap, 1);
    gm.block(2*dim_l1, dim_l2, dim1-2*dim_l1, dim_l2) = -GMGolyComplScalar(T, golycompl_basis, tens_family.ancestor(), 2, intmap, 0);
    
    return gm;
  };
  
/// Template to compute the Gram Matrix of a tensorized scalar basis and a GolyCompl basis
template<typename BasisType1, size_t N>
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                                         ///< Cell to which the basis corresponds
                    const TensorizedVectorFamily<BasisType1, N> & tens_family,  ///< First basis (tensorized basis)
                    const GolyComplBasisCell & golycompl_basis,             ///< Second basis (GolyCompl basis)
                    MonomialCellIntegralsType mono_int_map = {}                 ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    return GramMatrix(T, golycompl_basis, tens_family, mono_int_map).transpose();
  };
  
/// Computes the Gram Matrix of the mth component of a RolyCompl Basis and a monomial basis
Eigen::MatrixXd GMRolyComplScalar(
                    const Cell & T, ///< Cell to which the basis corresponds
                    const RolyComplBasisCell & rolycompl_basis, ///< First basis
                    const MonomialScalarBasisCell & mono_basis, ///< Second basis
                    const size_t m, ///< Add one to the power of the mth variable
                    MonomialCellIntegralsType mono_int_map = {} ///< list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );

/// Generic template to compute the Gram Matrix of the mth component of a RolyCompl Basis and any basis
template<typename BasisType>
Eigen::MatrixXd GMRolyComplScalar(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const RolyComplBasisCell & basis1, ///< First basis
                     const BasisType & basis2,  ///< Second basis (columns of the Gram matrix)
                     const size_t m, ///< Differentiate basis1 with respect to the mth variable
                     MonomialCellIntegralsType mono_int_map = {} ///< list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(BasisType::hasAncestor, "No method to compute this Gram matrix of derivatives");
    return transformGM(basis2, 'C', GMRolyComplScalar(T, basis1, basis2.ancestor(), m, mono_int_map) );
  };
  
/// Computes the Gram Matrix of the sth section of a GolyCompl Basis and a monomial basis
Eigen::MatrixXd GMGolyComplScalar(
                    const Cell & T, ///< Cell to which the basis corresponds
                    const GolyComplBasisCell & golycompl_basis, ///< First basis
                    const MonomialScalarBasisCell & mono_basis, ///< Second basis
                    const size_t s, ///< the GolyCompl basis has three sections; use the sth.
                    MonomialCellIntegralsType mono_int_map, ///< list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    const size_t m = 3, ///< Optionally (m < 3), add one to the power of the mth variable
                    const size_t k1 = 3, ///< Optionally (k1 < 3), take the kth derivative of the GolyCompl scalar
                    const size_t k2 = 3 ///< Optionally (k2 < 3), take the kth derivative of the MonomialScalar scalar
                    );
  
/// Generic template to compute the Gram Matrix of the sth section of a GolyCompl Basis and any basis with an extra power of the mth variable
template<typename BasisType>
Eigen::MatrixXd GMGolyComplScalar(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const GolyComplBasisCell & basis1, ///< First basis
                     const BasisType & basis2,  ///< Second basis (columns of the Gram matrix)
                     const size_t s, ///< the GolyCompl basis has three sections; use the sth.
                     MonomialCellIntegralsType mono_int_map, ///< list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     const size_t m = 3, ///< Optionally (m < 3), add one to the power of the mth variable
                     const size_t k1 = 3, ///< Optionally (k1 < 3), take the kth derivative of the GolyCompl scalar
                     const size_t k2 = 3 ///< Optionally (k2 < 3), take the kth derivative of the MonomialScalar scalar
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(BasisType::hasAncestor, "No method to compute this Gram matrix of derivatives");
    return transformGM(basis2, 'C', GMGolyComplScalar(T, basis1, basis2.ancestor(), s, mono_int_map, m, k1, k2) );
  };
  
/// Computes the Gram Matrix of the (optionally k1th derivative of the) s1th section of GolyCompl (optionally multiplied by the m1th variable) and the (optionally k2th derivative of the) s2th section of GolyCompl (optionally multiplied by the m2th variable)
Eigen::MatrixXd GMGolyCompl(
                    const Cell & T, ///< Cell to which the basis corresponds
                    const GolyComplBasisCell & basis1, ///< First basis
                    const GolyComplBasisCell & basis2, ///< Second basis
                    const size_t s1, ///< the first GolyCompl basis has three sections; use the s1th.
                    const size_t s2, ///< the second GolyCompl basis has three sections; use the s2th.
                    MonomialCellIntegralsType mono_int_map, ///< list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    const size_t m1 = 3, ///< Optionally (m1 < 3), add one to the power of the m1th variable
                    const size_t m2 = 3, ///< Optionally (m2 < 3), add one to the power of the m2th variable
                    const size_t k1 = 3, ///< Optionally (k1 < 3), take the k1th derivative of the first GolyCompl scalar
                    const size_t k2 = 3 ///< Optionally (k2 < 3), take the k2th derivative of the second GolyCompl scalar
                    );
  
  
/// Determines if the ancestor of a basis will be used to compute a Gram matrix for this basis
template<typename BasisType>
constexpr bool useAncestor()
  {
    // For a given basis, we only use the ancestor if it has an ancestor and has the same rank as its
    //  ancestor. If it has a different rank than its ancestor, it must be dealt with a specific overload

    if constexpr(BasisType::hasAncestor){
      return (BasisType::tensorRank == BasisType::AncestorType::tensorRank);
    } else {
      return false;
    }
  };

/// Generic template to compute the Gram Matrix of any pair of bases
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(const Cell& T,         ///< Cell to which the basis corresponds
                     const BasisType1 & basis1,   ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,   ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(useAncestor<BasisType1>() || useAncestor<BasisType2>(), "No method to compute this Gram matrix");
      
    if constexpr (!useAncestor<BasisType1>() && useAncestor<BasisType2>()) {
      return transformGM(basis2, 'C', GramMatrix(T, basis1, basis2.ancestor(), mono_int_map) );
    } else if constexpr (useAncestor<BasisType1>() && !useAncestor<BasisType2>()) {
      return transformGM(basis1, 'R', GramMatrix(T, basis1.ancestor(), basis2, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GramMatrix(T, basis1.ancestor(), basis2.ancestor(), mono_int_map) ) );
    }

  };
  
//
/* Gram matrix of scalar basis with one or two derivatives */

/// Computes the Gram Matrix of a pair of local scalar monomial bases, taking a partial derivative of the first (w.r.t. homogeneous coordinates, without scaling)
Eigen::MatrixXd GMScalarDerivative(
                  const Cell & T,  ///< Cell to which the basis corresponds
                  const MonomialScalarBasisCell & basis1, ///< First basis
                  const MonomialScalarBasisCell & basis2, ///< Second basis
                  const size_t m, ///< Differentiate basis1 with respect to the mth variable
                  MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                  );
  
/// Computes the Gram Matrix of a pair of local scalar monomial bases, taking partial derivatives of each of them (w.r.t. homogeneous coordinates, without scaling)
Eigen::MatrixXd GMScalarDerivative(
                  const Cell & T,  ///< Cell to which the basis corresponds
                  const MonomialScalarBasisCell & basis1, ///< First basis
                  const MonomialScalarBasisCell & basis2, ///< Second basis
                  const size_t m, ///< Differentiate basis1 with respect to the mth variable
                  const size_t l, ///< Differentiate basis2 with respect to the lth variable
                  MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                  );

/// Generic template to compute the Gram Matrix of any pair of scalar bases, taking a partial derivative of the first (w.r.t. homogeneous coordinates, without scaling)
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GMScalarDerivative(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const BasisType1 & basis1, ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,  ///< Second basis (columns of the Gram matrix)
                     const size_t m, ///< Differentiate basis1 with respect to the mth variable
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                             )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(BasisType1::hasAncestor || BasisType2::hasAncestor, "No method to compute this Gram matrix of derivatives");
      
    if constexpr (!BasisType1::hasAncestor && BasisType2::hasAncestor) {
      return transformGM(basis2, 'C', GMScalarDerivative(T, basis1, basis2.ancestor(), m, mono_int_map) );
    } else if constexpr (BasisType1::hasAncestor && !BasisType2::hasAncestor) {
      return transformGM(basis1, 'R', GMScalarDerivative(T, basis1.ancestor(), basis2, m, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), m, mono_int_map) ) );
    }
  };
  
/// Generic template to compute the Gram Matrix of any pair of scalar bases, taking partial derivatives of each of them  (w.r.t. homogeneous coordinates, without scaling)
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GMScalarDerivative(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const BasisType1 & basis1, ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,  ///< Second basis (columns of the Gram matrix)
                     const size_t m, ///< Differentiate basis1 with respect to the mth variable
                     const size_t l, ///< Differentiate basis2 with respect to the lth variable
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                             )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(BasisType1::hasAncestor || BasisType2::hasAncestor, "No method to compute this Gram matrix of derivatives");
      
    if constexpr (!BasisType1::hasAncestor && BasisType2::hasAncestor) {
      return transformGM(basis2, 'C', GMScalarDerivative(T, basis1, basis2.ancestor(), m, l, mono_int_map) );
    } else if constexpr (BasisType1::hasAncestor && !BasisType2::hasAncestor) {
      return transformGM(basis1, 'R', GMScalarDerivative(T, basis1.ancestor(), basis2, m, l, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), m, l, mono_int_map) ) );
    }
    
  };
  
//
/* Gram matrices of vector-valued basis with gradient, curl... */

/// Template to compute the Gram Matrix of a gradient basis and a tensorized scalar basis
template<typename BasisType1, typename BasisType2, size_t N>
Eigen::MatrixXd GramMatrix(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const GradientBasis<BasisType1> & grad_basis, ///< First basis (rows of the Gram matrix) - gradient basis
                     const TensorizedVectorFamily<BasisType2, N> & tens_family,  ///< Second basis (columns of the Gram matrix) - tensorized basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = grad_basis.dimension();
    size_t dim2 = tens_family.dimension();
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1, dim2);

    // Integrals of monomials
    size_t totaldegree = grad_basis.ancestor().max_degree()+tens_family.max_degree()-1;
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

    for (size_t m=0; m<N; m++){
      gm.block(0, m*dim2/N, dim1, dim2/N) = GMScalarDerivative(T, grad_basis.ancestor(), tens_family.ancestor(), m, intmap);
    }
    return gm/T.diam();
  };


/// Template to compute the Gram Matrix of a tensorized scalar basis and a gradient basis
template<typename BasisType1, typename BasisType2, size_t N>
Eigen::MatrixXd GramMatrix(
             const Cell& T, ///< Cell to which the basis corresponds
             const TensorizedVectorFamily<BasisType1, N> & tens_family, ///< First basis (rows of the Gram matrix) - gradient basis
             const GradientBasis<BasisType2> & grad_basis,  ///< Second basis (columns of the Gram matrix) - tensorized basis
             MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
             )
  {
    return GramMatrix(T, grad_basis, tens_family, mono_int_map).transpose();
  };


/// Template to compute the Gram Matrix of a gradient basis and another gradient basis
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const GradientBasis<BasisType1> & grad_basis1, ///< First basis (rows of the Gram matrix) - gradient basis
                     const GradientBasis<BasisType2> & grad_basis2,  ///< Second basis (columns of the Gram matrix) - gradient basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = grad_basis1.dimension();
    size_t dim2 = grad_basis2.dimension();
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1, dim2);

    // Integrals of monomials
    size_t totaldegree = grad_basis1.ancestor().max_degree()+grad_basis2.ancestor().max_degree()-2;
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

    for (size_t m=0; m<3; m++){
      gm += GMScalarDerivative(T, grad_basis1.ancestor(), grad_basis2.ancestor(), m, m, intmap);
    }
    return gm/std::pow(T.diam(), 2);
  };
  
/// Generic template to compute the Gram Matrix of Curl of any pair of bases
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(const Cell& T,         ///< Cell to which the basis corresponds
                     const CurlBasis<BasisType1> & basis1,   ///< First basis (rows of the Gram matrix)
                     const CurlBasis<BasisType2> & basis2,   ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    return GramMatrixCurlCurl(T, basis1.ancestor(), basis2.ancestor(), mono_int_map);
  };
  
/// Generic template to compute the Gram Matrix of a Curl basis and any other basis
template<typename BasisType1, typename BasisType2>
typename boost::disable_if<boost::is_same<BasisType2, MonomialCellIntegralsType>, Eigen::MatrixXd>::type GramMatrix(
                     const Cell& T,                          ///< Cell to which the basis corresponds
                     const CurlBasis<BasisType1> & basis1,   ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,              ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    return GramMatrixCurl(T, basis1.ancestor(), basis2, mono_int_map);
  };
  
/// Template to compute the Gram Matrix of any basis and a Curl basis
template<typename BasisType1, typename Basis2>
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                         ///< Cell to which the basis corresponds
                    const BasisType1 & basis1,                  ///< First basis (tensorized basis)
                    const CurlBasis<Basis2> & basis2,       ///< Second basis (Curl basis)
                    MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    return GramMatrixCurl(T, basis2.ancestor(), basis1, mono_int_map).transpose();
  };
  
/// Template to compute the Gram Matrix of the curl of a tensorized basis and the curl of another tensorized basis
template<typename BasisType1, typename BasisType2, size_t N>
Eigen::MatrixXd GramMatrixCurlCurl(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const TensorizedVectorFamily<BasisType1, N> & basis1, ///< First basis (rows of the Gram matrix) - curl basis
                     const TensorizedVectorFamily<BasisType2, N> & basis2,  ///< Second basis (columns of the Gram matrix) - curl basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();
    Eigen::MatrixXd gm(dim1, dim2);

    // Integrals of monomials
    size_t totaldegree = basis1.max_degree()+basis2.max_degree()-2;
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
    
    Eigen::MatrixXd GMxx = GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), 0, 0, intmap);
    Eigen::MatrixXd GMyy = GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), 1, 1, intmap);
    Eigen::MatrixXd GMzz = GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), 2, 2, intmap);
    
    gm.block(0, 0, dim1/N, dim2/N) = GMzz + GMyy;
    gm.block(0, dim2/N, dim1/N, dim2/N) = -GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), 1, 0, intmap);
    gm.block(0, 2*dim2/N, dim1/N, dim2/N) = -GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), 2, 0, intmap);
    gm.block(dim1/N, 0, dim1/N, dim2/N) = -GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), 0, 1, intmap);
    gm.block(dim1/N, dim2/N, dim1/N, dim2/N) = GMzz + GMxx;
    gm.block(dim1/N, 2*dim2/N, dim1/N, dim2/N) = -GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), 2, 1, intmap);
    gm.block(2*dim1/N, 0, dim1/N, dim2/N) = -GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), 0, 2, intmap);
    gm.block(2*dim1/N, dim2/N, dim1/N, dim2/N) = -GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), 1, 2, intmap);
    gm.block(2*dim1/N, 2*dim2/N, dim1/N, dim2/N) = GMyy + GMxx;
    
    return gm/std::pow(T.diam(), 2);
  };

/// Compute the Gram Matrix of the curl of a GolyCompl basis and the curl of another GolyCompl basis
Eigen::MatrixXd GramMatrixCurlCurl(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const GolyComplBasisCell & basis1, ///< First basis (rows of the Gram matrix) - curl basis
                     const GolyComplBasisCell & basis2,  ///< Second basis (columns of the Gram matrix) - curl basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     );

/// Template to compute the Gram Matrix of the curl of a GolyCompl basis and the curl of a tensorized basis
template<typename BasisType2, size_t N>
Eigen::MatrixXd GramMatrixCurlCurl(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const GolyComplBasisCell & basis1, ///< First basis (rows of the Gram matrix) - curl basis
                     const TensorizedVectorFamily<BasisType2, N> & basis2,  ///< Second basis (columns of the Gram matrix) - curl basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();
    Eigen::MatrixXd gm(dim1, dim2);
    
    // Integrals of monomials
    size_t totaldegree = basis1.max_degree()+basis2.max_degree();
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
    
    size_t dim_l1 = basis1.dimPkmo();
    gm.block(0, 0, dim_l1, dim2/N) = GMGolyComplScalar(T, basis1, basis2.ancestor(), 0, intmap, 1, 0, 2) - GMGolyComplScalar(T, basis1, basis2.ancestor(), 0, intmap, 2, 0, 1);
    gm.block(0, dim2/N, dim_l1, dim2/N) = GMGolyComplScalar(T, basis1, basis2.ancestor(), 0, intmap, 2, 0, 0) + GMGolyComplScalar(T, basis1, basis2.ancestor(), 0, intmap, 1, 1, 2)
                                        + GMGolyComplScalar(T, basis1, basis2.ancestor(), 0, intmap, 2, 2, 2) + 2/T.diam()*GMGolyComplScalar(T, basis1, basis2.ancestor(), 0, intmap, 3, 3, 2);
    gm.block(0, 2*dim2/N, dim_l1, dim2/N) = -GMGolyComplScalar(T, basis1, basis2.ancestor(), 0, intmap, 1, 0, 0) - GMGolyComplScalar(T, basis1, basis2.ancestor(), 0, intmap, 1, 1, 1)
                                            -GMGolyComplScalar(T, basis1, basis2.ancestor(), 0, intmap, 2, 2, 1) - 2/T.diam()*GMGolyComplScalar(T, basis1, basis2.ancestor(), 0, intmap, 3, 3, 1);
    gm.block(dim_l1, 0, dim_l1, dim2/N) = -GMGolyComplScalar(T, basis1, basis2.ancestor(), 1, intmap, 2, 1, 1) - GMGolyComplScalar(T, basis1, basis2.ancestor(), 1, intmap, 0, 0, 2)
                                          -GMGolyComplScalar(T, basis1, basis2.ancestor(), 1, intmap, 2, 2, 2) - 2/T.diam()*GMGolyComplScalar(T, basis1, basis2.ancestor(), 1, intmap, 3, 3, 2);
    gm.block(dim_l1, dim2/N, dim_l1, dim2/N) = -GMGolyComplScalar(T, basis1, basis2.ancestor(), 1, intmap, 0, 1, 2) + GMGolyComplScalar(T, basis1, basis2.ancestor(), 1, intmap, 2, 1, 0);
    gm.block(dim_l1, 2*dim2/N, dim_l1, dim2/N) = GMGolyComplScalar(T, basis1, basis2.ancestor(), 1, intmap, 0, 1, 1) + GMGolyComplScalar(T, basis1, basis2.ancestor(), 1, intmap, 0, 0, 0)
                                                +GMGolyComplScalar(T, basis1, basis2.ancestor(), 1, intmap, 2, 2, 0) + 2/T.diam()*GMGolyComplScalar(T, basis1, basis2.ancestor(), 1, intmap, 3, 3, 0);
    gm.block(2*dim_l1, 0, dim1-2*dim_l1, dim2/N) = -GMGolyComplScalar(T, basis1, basis2.ancestor(), 2, intmap, 1, 2, 2) + GMGolyComplScalar(T, basis1, basis2.ancestor(), 2, intmap, 0, 0, 1)
                                                   +GMGolyComplScalar(T, basis1, basis2.ancestor(), 2, intmap, 1, 1, 1) + 2/T.diam()*GMGolyComplScalar(T, basis1, basis2.ancestor(), 2, intmap, 3, 3, 1);
    gm.block(2*dim_l1, dim2/N, dim1-2*dim_l1, dim2/N) = -GMGolyComplScalar(T, basis1, basis2.ancestor(), 2, intmap, 0, 2, 2) - GMGolyComplScalar(T, basis1, basis2.ancestor(), 2, intmap, 0, 0, 0)
                                                        -GMGolyComplScalar(T, basis1, basis2.ancestor(), 2, intmap, 1, 1, 0) - 2/T.diam()*GMGolyComplScalar(T, basis1, basis2.ancestor(), 2, intmap, 3, 3, 0);
    gm.block(2*dim_l1, 2*dim2/N, dim1-2*dim_l1, dim2/N) = GMGolyComplScalar(T, basis1, basis2.ancestor(), 2, intmap, 0, 2, 1) - GMGolyComplScalar(T, basis1, basis2.ancestor(), 2, intmap, 1, 2, 0);
    
    return gm;
  };

/// Template to compute the Gram Matrix of the curl of a tensorized basis and the curl of a GolyCompl basis
template<typename BasisType1, size_t N>
Eigen::MatrixXd GramMatrixCurlCurl(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const TensorizedVectorFamily<BasisType1, N> & basis1, ///< First basis (rows of the Gram matrix) - curl basis
                     const GolyComplBasisCell & basis2,  ///< Second basis (columns of the Gram matrix) - curl basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    return GramMatrixCurlCurl(T, basis2, basis1, mono_int_map).transpose();
  };

/// Template to compute the Gram Matrix of the curl of any basis and the curl of any other basis
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrixCurlCurl(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const BasisType1 & basis1, ///< First basis (rows of the Gram matrix) - curl basis
                     const BasisType2 & basis2,  ///< Second basis (columns of the Gram matrix) - curl basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(useAncestor<BasisType1>() || useAncestor<BasisType2>(), "No method to compute this Gram matrix");
      
    if constexpr (!useAncestor<BasisType1>() && useAncestor<BasisType2>()) {
      return transformGM(basis2, 'C', GramMatrixCurlCurl(T, basis1, basis2.ancestor(), mono_int_map) );
    } else if constexpr (useAncestor<BasisType1>() && !useAncestor<BasisType2>()) {
      return transformGM(basis1, 'R', GramMatrixCurlCurl(T, basis1.ancestor(), basis2, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GramMatrixCurlCurl(T, basis1.ancestor(), basis2.ancestor(), mono_int_map) ) );
    }
  };

/// Template to compute the Gram Matrix of a Curl<Tensorized> basis and a tensorized scalar basis
template<typename BasisType1, typename BasisType2, size_t N>
Eigen::MatrixXd GramMatrixCurl(
                    const Cell & T,                                         ///< Cell to which the basis corresponds
                    const TensorizedVectorFamily<BasisType1, N> & basis1,   ///< First basis (tensorized basis)
                    const TensorizedVectorFamily<BasisType2, N> & basis2,   ///< Second basis (tensorized basis)
                    MonomialCellIntegralsType mono_int_map = {}                 ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1, dim2);

    // Integrals of monomials
    size_t totaldegree = basis1.max_degree()+basis2.max_degree()-1;
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

    Eigen::MatrixXd GMx = GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), 0, intmap);
    Eigen::MatrixXd GMy = GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), 1, intmap);
    Eigen::MatrixXd GMz = GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), 2, intmap);
    
    gm.block(0, dim2/N, dim1/N, dim2/N) = GMz;
    gm.block(0, 2*dim2/N, dim1/N, dim2/N) = -GMy;
    gm.block(dim1/N, 0, dim1/N, dim2/N) = -GMz;
    gm.block(dim1/N, 2*dim2/N, dim1/N, dim2/N) = GMx;
    gm.block(2*dim1/N, 0, dim1/N, dim2/N) = GMy;
    gm.block(2*dim1/N, dim2/N, dim1/N, dim2/N) = -GMx;
    
    return gm/T.diam();
  };
  
/// Template to compute the Gram Matrix of a Curl<GolyCompl> basis and a tensorized scalar basis
template<typename BasisType2, size_t N>
Eigen::MatrixXd GramMatrixCurl(
                    const Cell & T,                                         ///< Cell to which the basis corresponds
                    const GolyComplBasisCell & basis1,                      ///< First basis (GolyCompl basis)
                    const TensorizedVectorFamily<BasisType2, N> & basis2,   ///< Second basis (tensorized basis)
                    MonomialCellIntegralsType mono_int_map = {}                 ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();
    Eigen::MatrixXd gm(dim1, dim2);
    
    // Integrals of monomials
    size_t totaldegree = basis1.max_degree()+basis2.max_degree();
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
    
    size_t dim_l1 = basis1.dimPkmo();
    gm.block(0, 0, dim_l1, dim2/N) = -GMGolyComplScalar(T, basis1, basis2.ancestor(), 0, intmap, 1, 1)
                                     -GMGolyComplScalar(T, basis1, basis2.ancestor(), 0, intmap, 2, 2)
                                     -2*GMGolyComplScalar(T, basis1, basis2.ancestor(), 0, intmap)/T.diam();
    gm.block(0, dim2/N, dim_l1, dim2/N) = GMGolyComplScalar(T, basis1, basis2.ancestor(), 0, intmap, 1, 0);
    gm.block(0, 2*dim2/N, dim_l1, dim2/N) = GMGolyComplScalar(T, basis1, basis2.ancestor(), 0, intmap, 2, 0);
    gm.block(dim_l1, 0, dim_l1, dim2/N) = GMGolyComplScalar(T, basis1, basis2.ancestor(), 1, intmap, 0, 1);
    gm.block(dim_l1, dim2/N, dim_l1, dim2/N) = -GMGolyComplScalar(T, basis1, basis2.ancestor(), 1, intmap, 0, 0)
                                               -GMGolyComplScalar(T, basis1, basis2.ancestor(), 1, intmap, 2, 2)
                                               -2*GMGolyComplScalar(T, basis1, basis2.ancestor(), 1, intmap)/T.diam();
    gm.block(dim_l1, 2*dim2/N, dim_l1, dim2/N) = GMGolyComplScalar(T, basis1, basis2.ancestor(), 1, intmap, 2, 1);
    gm.block(2*dim_l1, 0, dim1-2*dim_l1, dim2/N) = GMGolyComplScalar(T, basis1, basis2.ancestor(), 2, intmap, 0, 2);
    gm.block(2*dim_l1, dim2/N, dim1-2*dim_l1, dim2/N) = GMGolyComplScalar(T, basis1, basis2.ancestor(), 2, intmap, 1, 2);
    gm.block(2*dim_l1, 2*dim2/N, dim1-2*dim_l1, dim2/N) = -GMGolyComplScalar(T, basis1, basis2.ancestor(), 2, intmap, 0, 0)
                                                          -GMGolyComplScalar(T, basis1, basis2.ancestor(), 2, intmap, 1, 1)
                                                          -2*GMGolyComplScalar(T, basis1, basis2.ancestor(), 2, intmap)/T.diam();
    
    return gm;
  };
  
/// Template to compute the Gram Matrix of the curl of any basis and any other basis
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrixCurl(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const BasisType1 & basis1, ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,  ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(useAncestor<BasisType1>() || useAncestor<BasisType2>(), "No method to compute this Gram matrix");
      
    if constexpr (!useAncestor<BasisType1>() && useAncestor<BasisType2>()) {
      return transformGM(basis2, 'C', GramMatrixCurl(T, basis1, basis2.ancestor(), mono_int_map) );
    } else if constexpr (useAncestor<BasisType1>() && !useAncestor<BasisType2>()) {
      return transformGM(basis1, 'R', GramMatrixCurl(T, basis1.ancestor(), basis2, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GramMatrixCurl(T, basis1.ancestor(), basis2.ancestor(), mono_int_map) ) );
    }
  };

/// Template to compute the Gram Matrix of the curl of any two bases when one CurlBasis is at a lower level than the other
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrixCurl(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const BasisType1 & basis1, ///< First basis (rows of the Gram matrix)
                     const CurlBasis<BasisType2> & basis2,  ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    if constexpr (!useAncestor<BasisType1>()) {
      return GramMatrixCurlCurl(T, basis1, basis2.ancestor(), mono_int_map);
    } else {
      return transformGM(basis1, 'R', GramMatrixCurlCurl(T, basis1.ancestor(), basis2.ancestor(), mono_int_map) );
    }
  };
 
/// Generic template to compute the Gram Matrix of Divergence of any pair of bases
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(const Cell& T,         ///< Cell to which the basis corresponds
                     const DivergenceBasis<BasisType1> & basis1,   ///< First basis (rows of the Gram matrix)
                     const DivergenceBasis<BasisType2> & basis2,   ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    return GramMatrixDivDiv(T, basis1.ancestor(), basis2.ancestor(), mono_int_map);
  };
 
/// Gram Matrix of the divergence of a RolyCompl Basis and the kth derivative of a monomial basis
Eigen::MatrixXd GMRolyComplScalarDiv(
                    const Cell & T, ///< Cell to which the basis corresponds
                    const MonomialScalarBasisCell & mono_basis, ///< First basis
                    const RolyComplBasisCell & rolycompl_basis, ///< Second basis
                    const size_t k, ///< Take the kth derivative of the monomial scalar
                    MonomialCellIntegralsType mono_int_map ///< list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );
  
/// Gram Matrix of the divergence of a RolyCompl basis and the k-th derivative of any scalar basis
template<typename BasisType>
Eigen::MatrixXd GMRolyComplScalarDiv(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const BasisType & basis1, ///< First basis
                     const RolyComplBasisCell & basis2, ///< Second basis (columns of the Gram matrix)
                     const size_t k, ///< Take the kth derivative of the scalar basis
                     MonomialCellIntegralsType mono_int_map ///< list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(BasisType::hasAncestor, "No method to compute this Gram matrix of derivatives");
    return transformGM(basis1, 'R', GMRolyComplScalarDiv(T, basis1.ancestor(), basis2, k, mono_int_map) );
  };
  
/// Template to compute the Gram Matrix of the divergence of a tensorized basis and the divergence of another tensorized basis
template<typename BasisType1, typename BasisType2, size_t N>
Eigen::MatrixXd GramMatrixDivDiv(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const TensorizedVectorFamily<BasisType1, N> & basis1, ///< First basis (rows of the Gram matrix)
                     const TensorizedVectorFamily<BasisType2, N> & basis2,  ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();
    Eigen::MatrixXd gm(dim1, dim2);

    // Integrals of monomials
    size_t totaldegree = basis1.max_degree()+basis2.max_degree()-2;
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
    
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        gm.block(i*dim1/N, j*dim2/N, dim1/N, dim2/N) = GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), i, j, intmap);
      }
    }
    
    return gm/std::pow(T.diam(), 2);
  };
 
/// Compute the Gram Matrix of the divergence of a RolyCompl basis and the divergence of another RolyCompl basis
Eigen::MatrixXd GramMatrixDivDiv(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const RolyComplBasisCell & basis1, ///< First basis (rows of the Gram matrix)
                     const RolyComplBasisCell & basis2,  ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     );

/// Template to compute the Gram Matrix of the divergence of a RolyCompl basis and the divergence of a tensorized basis
template<typename BasisType2, size_t N>
Eigen::MatrixXd GramMatrixDivDiv(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const RolyComplBasisCell & basis1, ///< First basis (rows of the Gram matrix)
                     const TensorizedVectorFamily<BasisType2, N> & basis2,  ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    return GramMatrixDivDiv(T, basis2, basis1, mono_int_map).transpose();
  };

/// Template to compute the Gram Matrix of the divergence of a tensorized basis and the divergence of a RolyCompl basis
template<typename BasisType1, size_t N>
Eigen::MatrixXd GramMatrixDivDiv(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const TensorizedVectorFamily<BasisType1, N> & basis1, ///< First basis (rows of the Gram matrix)
                     const RolyComplBasisCell & basis2,  ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();
    Eigen::MatrixXd gm(dim1, dim2);
    
    // Integrals of monomials
    size_t totaldegree = basis1.max_degree()+basis2.max_degree();
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
    
    for (size_t i = 0; i < N; i++) {
      gm.block(i*dim1/N, 0, dim1/N, dim2) = GMRolyComplScalarDiv(T, basis1.ancestor(), basis2, i, mono_int_map);
    }
    
    return gm;
  };

/// Template to compute the Gram Matrix of the divergence of any basis and the divergence of any other basis
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrixDivDiv(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const BasisType1 & basis1, ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,  ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(useAncestor<BasisType1>() || useAncestor<BasisType2>(), "No method to compute this Gram matrix");
      
    if constexpr (!useAncestor<BasisType1>() && useAncestor<BasisType2>()) {
      return transformGM(basis2, 'C', GramMatrixDivDiv(T, basis1, basis2.ancestor(), mono_int_map) );
    } else if constexpr (useAncestor<BasisType1>() && !useAncestor<BasisType2>()) {
      return transformGM(basis1, 'R', GramMatrixDivDiv(T, basis1.ancestor(), basis2, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GramMatrixDivDiv(T, basis1.ancestor(), basis2.ancestor(), mono_int_map) ) );
    }
  };

/// Generic template to compute the Gram Matrix of a Divergence basis and any other basis
template<typename BasisType1, typename BasisType2>
typename boost::disable_if<boost::is_same<BasisType2, MonomialCellIntegralsType>, Eigen::MatrixXd>::type GramMatrix(
                     const Cell& T,                          ///< Cell to which the basis corresponds
                     const DivergenceBasis<BasisType1> & basis1,   ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,              ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    return GramMatrixDiv(T, basis1.ancestor(), basis2, mono_int_map);
  };
  
/// Template to compute the Gram Matrix of any basis and a Divergence basis
template<typename BasisType1, typename Basis2>
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                         ///< Cell to which the basis corresponds
                    const BasisType1 & basis1,                  ///< First basis (scalar basis)
                    const DivergenceBasis<Basis2> & basis2,       ///< Second basis (divergence basis)
                    MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    return GramMatrixDiv(T, basis2.ancestor(), basis1, mono_int_map).transpose();
  };
  
/// Template to compute the Gram Matrix of a Divergence<Tensorized> basis and a monomial scalar basis
template<typename BasisType1, size_t N>
Eigen::MatrixXd GramMatrixDiv(
                    const Cell & T,                                         ///< Cell to which the basis corresponds
                    const TensorizedVectorFamily<BasisType1, N> & basis1,   ///< First basis (divergence basis)
                    const MonomialScalarBasisCell & basis2,                 ///< Second basis (monomial scalar basis)
                    MonomialCellIntegralsType mono_int_map = {}                 ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1, dim2);

    // Integrals of monomials
    size_t totaldegree = basis1.max_degree()+basis2.max_degree()-1;
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

    for (size_t i = 0; i < N; i++) {
      gm.block(i*dim1/N, 0, dim1/N, dim2) = GMScalarDerivative(T, basis1.ancestor(), basis2, i, intmap);
    }
    
    return gm/T.diam();
  };
  
/// Computes the Gram Matrix of a Divergence<RolyCompl> basis and a monomial scalar basis
Eigen::MatrixXd GramMatrixDiv(
                    const Cell & T,                           ///< Cell to which the basis corresponds
                    const RolyComplBasisCell & basis1,        ///< First basis (RolyCompl basis)
                    const MonomialScalarBasisCell & basis2,   ///< Second basis (scalar basis)
                    MonomialCellIntegralsType mono_int_map = {}   ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );
  
/// Template to compute the Gram Matrix of the divergence of any basis and any other basis
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrixDiv(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const BasisType1 & basis1, ///< First basis (vector basis)
                     const BasisType2 & basis2,  ///< Second basis (scalar basis)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(useAncestor<BasisType1>() || useAncestor<BasisType2>(), "No method to compute this Gram matrix");
      
    if constexpr (!useAncestor<BasisType1>() && useAncestor<BasisType2>()) {
      return transformGM(basis2, 'C', GramMatrixDiv(T, basis1, basis2.ancestor(), mono_int_map) );
    } else if constexpr (useAncestor<BasisType1>() && !useAncestor<BasisType2>()) {
      return transformGM(basis1, 'R', GramMatrixDiv(T, basis1.ancestor(), basis2, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GramMatrixDiv(T, basis1.ancestor(), basis2.ancestor(), mono_int_map) ) );
    }
  };

/// Template to compute the Gram Matrix of the divergence of any two bases when one DivergenceBasis is at a lower level than the other
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrixDiv(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const BasisType1 & basis1, ///< First basis (vector basis)
                     const DivergenceBasis<BasisType2> & basis2,  ///< Second basis (divergence basis)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    if constexpr (!useAncestor<BasisType1>()) {
      return GramMatrixDivDiv(T, basis1, basis2.ancestor(), mono_int_map);
    } else {
      return transformGM(basis1, 'R', GramMatrixDivDiv(T, basis1.ancestor(), basis2.ancestor(), mono_int_map) );
    }
  };
 

/*@}*/
} // end of namespace HArDCore3D

#endif // end of _GMPOLY_CELL_HPP
