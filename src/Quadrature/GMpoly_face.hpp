# ifndef _GMPOLY_FACE_HPP
# define _GMPOLY_FACE_HPP

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

#include <GMpoly_cell.hpp>

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
struct VecFaceHash
{
  size_t operator()(const Eigen::Vector2i & p) const
  {
    constexpr size_t b = 100; // max degree must be less than the basis b
    return p(0) + p(1) * b;
  }
};

/// Type for list of face integrals of monomials
typedef std::unordered_map<Eigen::Vector2i, double, VecFaceHash> MonomialFaceIntegralsType;
    
/// Compute all integrals, on the edges of a face, of the face monomials up to a total degree
std::vector<MonomialFaceIntegralsType> IntegrateFaceMonomials_onEdges
          (const Face & F,      ///< Face
          const size_t maxdeg   ///< Maximal total degree
          );

/// Compute all integrals on a face of face monomials up to a total degree
MonomialFaceIntegralsType IntegrateFaceMonomials
          (const Face & F,      ///< Face
          const int maxdeg   ///< Maximal total degree
          );
          
/// Checks if the degree of an existing list of monomial integrals is sufficient, other re-compute and return a proper list
MonomialFaceIntegralsType CheckIntegralsDegree
         (const Face & F,              ///< Face
           const size_t degree,         ///< Expected degree
           const MonomialFaceIntegralsType & mono_int_map = {}    ///< Existing list, optional
          );


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//  FUNCTIONS TO COMPUTE GRAM MATRICES
//----------------------------------------------------------------------
//----------------------------------------------------------------------

//
/* BASIC ONES: Monomial, tensorized, and generic template */

/// Generic template to compute the Gram Matrix of any pair of bases
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(const Face& F,         ///< Face to which the basis corresponds
                     const BasisType1 & basis1,   ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,   ///< Second basis (columns of the Gram matrix)
                     MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(useAncestor<BasisType1>() || useAncestor<BasisType2>(), "No method to compute this Gram matrix");
      
    if constexpr (!useAncestor<BasisType1>() && useAncestor<BasisType2>()) {
      return transformGM(basis2, 'C', GramMatrix(F, basis1, basis2.ancestor(), mono_int_map) );
    } else if constexpr (useAncestor<BasisType1>() && !useAncestor<BasisType2>()) {
      return transformGM(basis1, 'R', GramMatrix(F, basis1.ancestor(), basis2, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GramMatrix(F, basis1.ancestor(), basis2.ancestor(), mono_int_map) ) );
    }

  };

/// Computes the Gram Matrix of a pair of local scalar monomial bases
Eigen::MatrixXd GramMatrix(
                    const Face & F,                         ///< Face to which the basis corresponds
                    const MonomialScalarBasisFace & basis1, ///< First basis
                    const MonomialScalarBasisFace & basis2, ///< Second basis
                    MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );
  
/// This overload to simplify the call to GramMatrix in case the two bases are the same
template<typename BasisType>
Eigen::MatrixXd GramMatrix(const Face& F, const BasisType & basis, MonomialFaceIntegralsType mono_int_map = {})
  {
    return GramMatrix(F, basis, basis, mono_int_map);
  };

/// Template to compute the Gram Matrix of any pair of tensorized scalar bases
template<typename BasisType1, typename BasisType2, size_t N>
Eigen::MatrixXd GramMatrix(
                     const Face& F, ///< Face to which the basis corresponds
                     const TensorizedVectorFamily<BasisType1, N> & basis1, ///< First basis (rows of the Gram matrix)
                     const TensorizedVectorFamily<BasisType2, N> & basis2,  ///< Second basis (columns of the Gram matrix)
                     MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(basis1.dimension(), basis2.dimension());

    Eigen::MatrixXd anc_gm = GramMatrix(F, basis1.ancestor(), basis2.ancestor(), mono_int_map);
    size_t dim1 = anc_gm.rows();
    size_t dim2 = anc_gm.cols();

    for (size_t i=0; i<N; i++){
      gm.block(i*dim1, i*dim2, dim1, dim2) = anc_gm;    
    }
    return gm;
  };
  
/// Computes the Gram Matrix of a pair of tangent bases
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(
                    const Face & F,                             ///< Face to which the basis corresponds
                    const TangentFamily<BasisType1> & basis1,   ///< First basis
                    const TangentFamily<BasisType2> & basis2,   ///< Second basis
                    MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    Eigen::MatrixXd gm(basis1.dimension(), basis2.dimension());
    
    Eigen::MatrixXd anc_gm = GramMatrix(F, basis1.ancestor(), basis2.ancestor(), mono_int_map);
    size_t dim1 = anc_gm.rows();
    size_t dim2 = anc_gm.cols();
    
    Eigen::Matrix<double, 2, 3> generators1 = basis1.generators();
    Eigen::Matrix<double, 2, 3> generators2 = basis2.generators();
    
    for (size_t i=0; i<2; i++) {
      for (size_t j=0; j<2; j++) {
        gm.block(i*dim1, j*dim2, dim1, dim2) = generators1.row(i).dot(generators2.row(j)) * anc_gm;
      }
    }
    
    return gm;
  };

/// Computes the Gram Matrix of a pair of RolyCompl bases
Eigen::MatrixXd GramMatrix(
                    const Face & F,                         ///< Face to which the basis corresponds
                    const RolyComplBasisFace & basis1,      ///< First basis
                    const RolyComplBasisFace & basis2,      ///< Second basis
                    MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );

/// Computes the Gram Matrix of a pair of GolyCompl bases
Eigen::MatrixXd GramMatrix(
                    const Face & F,                         ///< Face to which the basis corresponds
                    const GolyComplBasisFace & basis1,      ///< First basis
                    const GolyComplBasisFace & basis2,      ///< Second basis
                    MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );

/// Computes the Gram Matrix of a RolyCompl basis and a tangent basis
template<typename BasisType>
Eigen::MatrixXd GramMatrix(
                    const Face & F,                             ///< Face to which the basis corresponds
                    const RolyComplBasisFace & basis1,          ///< First basis
                    const TangentFamily<BasisType> & basis2,    ///< Second basis
                    MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension()/2;
    size_t totaldegree = basis1.max_degree()+basis2.max_degree();
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1, 2*dim2);
    
    // Obtain integration data from IntegrateFaceMonomials
    MonomialFaceIntegralsType intmap = CheckIntegralsDegree(F, totaldegree, mono_int_map);
    
    std::vector<Eigen::MatrixXd> gmrs = {GMRolyComplScalar(F, basis1, basis2.ancestor(), 0, intmap), GMRolyComplScalar(F, basis1, basis2.ancestor(), 1, intmap)};
    Eigen::Matrix<double, 2, 3> jacobian1 = basis1.jacobian();
    Eigen::Matrix<double, 2, 3> generators2 = basis2.generators();
    Eigen::Matrix2d W = jacobian1 * generators2.transpose();
    
    for (size_t k=0; k<2; k++) {
      for (size_t j=0; j<2; j++) {
        gm.block(0, j*dim2, dim1, dim2) += W(k,j) * gmrs[k];
      }
    }
    
    return F.diam() * gm;
  };

/// Computes the Gram Matrix of a RolyCompl basis and a tangent basis
template<typename BasisType>
Eigen::MatrixXd GramMatrix(
                    const Face & F,                             ///< Face to which the basis corresponds
                    const TangentFamily<BasisType> & basis1,    ///< First basis
                    const RolyComplBasisFace & basis2,          ///< Second basis
                    MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    return GramMatrix(F, basis2, basis1, mono_int_map).transpose();
  };

/// Computes the Gram Matrix of the scalar part of a RolyCompl Basis and a monomial basis with an extra power on the mth variable
Eigen::MatrixXd GMRolyComplScalar(
                    const Face & F, ///< Face to which the basis corresponds
                    const RolyComplBasisFace & rolycompl_basis, ///< First basis
                    const MonomialScalarBasisFace & mono_basis, ///< Second basis
                    const size_t m, ///< Add one to the power of the mth variable
                    MonomialFaceIntegralsType mono_int_map = {} ///< list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );

/// Generic template to compute the Gram Matrix of the scalar part of a RolyCompl Basis and any basis with an extra power on the mth variable
template<typename BasisType>
Eigen::MatrixXd GMRolyComplScalar(
                     const Face & F, ///< Face to which the basis corresponds
                     const RolyComplBasisFace & basis1, ///< First basis
                     const BasisType & basis2,  ///< Second basis (columns of the Gram matrix)
                     const size_t m, ///< Differentiate basis1 with respect to the mth variable
                     MonomialFaceIntegralsType mono_int_map = {} ///< list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(BasisType::hasAncestor, "No method to compute this Gram matrix of derivatives");
    return transformGM(basis2, 'C', GMRolyComplScalar(F, basis1, basis2.ancestor(), m, mono_int_map) );
  };
  
/// Computes the Gram Matrix of a pair of local scalar monomial bases, taking a partial derivative of the first (w.r.t. homogeneous coordinates on the face, no change of variable)
Eigen::MatrixXd GMScalarDerivative(
                  const Face & F,  ///< Face to which the basis corresponds
                  const MonomialScalarBasisFace & basis1, ///< First basis
                  const MonomialScalarBasisFace & basis2, ///< Second basis
                  const size_t m, ///< Differentiate basis1 with respect to the mth variable
                  MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                  );
  
/// Computes the Gram Matrix of a pair of local scalar monomial bases, taking a partial derivative of each (w.r.t. homogeneous coordinates on the face, no change of variable)
Eigen::MatrixXd GMScalarDerivative(
                  const Face & F,  ///< Face to which the basis corresponds
                  const MonomialScalarBasisFace & basis1, ///< First basis
                  const MonomialScalarBasisFace & basis2, ///< Second basis
                  const size_t m, ///< Differentiate basis1 with respect to the mth variable
                  const size_t l, ///< Differentiate basis2 with respect to the lth variable
                  MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                  );

/// Generic template to compute the Gram Matrix of any pair of scalar bases, taking a partial derivative of the first (w.r.t. homogeneous coordinates on the face, no change of variable)
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GMScalarDerivative(
                     const Face& F, ///< Face to which the basis corresponds
                     const BasisType1 & basis1, ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,  ///< Second basis (columns of the Gram matrix)
                     const size_t m, ///< Differentiate basis1 with respect to the mth variable
                     MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                             )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(BasisType1::hasAncestor || BasisType2::hasAncestor, "No method to compute this Gram matrix of derivatives");
      
    if constexpr (!BasisType1::hasAncestor && BasisType2::hasAncestor) {
      return transformGM(basis2, 'C', GMScalarDerivative(F, basis1, basis2.ancestor(), m, mono_int_map) );
    } else if constexpr (BasisType1::hasAncestor && !BasisType2::hasAncestor) {
      return transformGM(basis1, 'R', GMScalarDerivative(F, basis1.ancestor(), basis2, m, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GMScalarDerivative(F, basis1.ancestor(), basis2.ancestor(), m, mono_int_map) ) );
    }
  };
  
/// Generic template to compute the Gram Matrix of any pair of scalar bases, taking a partial derivative of each (w.r.t. homogeneous coordinates on the face, no change of variable)
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GMScalarDerivative(
                     const Face& F, ///< Face to which the basis corresponds
                     const BasisType1 & basis1, ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,  ///< Second basis (columns of the Gram matrix)
                     const size_t m, ///< Differentiate basis1 with respect to the mth variable
                     const size_t l, ///< Differentiate basis2 with respect to the lth variable
                     MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                             )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(BasisType1::hasAncestor || BasisType2::hasAncestor, "No method to compute this Gram matrix of derivatives");
      
    if constexpr (!BasisType1::hasAncestor && BasisType2::hasAncestor) {
      return transformGM(basis2, 'C', GMScalarDerivative(F, basis1, basis2.ancestor(), m, l, mono_int_map) );
    } else if constexpr (BasisType1::hasAncestor && !BasisType2::hasAncestor) {
      return transformGM(basis1, 'R', GMScalarDerivative(F, basis1.ancestor(), basis2, m, l, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GMScalarDerivative(F, basis1.ancestor(), basis2.ancestor(), m, l, mono_int_map) ) );
    }
    
  };
  
/// Generic template to compute the Gram Matrix of a pair of Curl bases
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(
                    const Face & F,                             ///< Face to which the basis corresponds
                    const CurlBasis<BasisType1> & basis1,       ///< First basis
                    const CurlBasis<BasisType2> & basis2,       ///< Second basis
                    MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    // Dimension of the gram matrix
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();
    size_t totaldegree = basis1.ancestor().max_degree()+basis2.ancestor().max_degree()-2;
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
    
    // Obtain integration data from IntegrateFaceMonomials
    MonomialFaceIntegralsType intmap = CheckIntegralsDegree(F, totaldegree, mono_int_map);
    
    for (size_t m=0; m<2; m++){
      gm += GMScalarDerivative(F, basis1.ancestor(), basis2.ancestor(), m, m, intmap);
    }
    
    return gm/std::pow(F.diam(), 2);
  }

/// Generic template to compute the Gram Matrix of a Curl basis and a Tangent basis
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(
                    const Face & F,                             ///< Face to which the basis corresponds
                    const CurlBasis<BasisType1> & basis1,       ///< First basis
                    const TangentFamily<BasisType2> & basis2,   ///< Second basis
                    MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    // Dimension of the gram matrix
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension()/2;
    size_t totaldegree = basis1.ancestor().max_degree()+basis2.max_degree()-1;
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2*2);
    
    // Obtain integration data from IntegrateFaceMonomials
    MonomialFaceIntegralsType intmap = CheckIntegralsDegree(F, totaldegree, mono_int_map);
    
    std::vector<Eigen::MatrixXd> gmd = {GMScalarDerivative(F, basis1.ancestor(), basis2.ancestor(), 0, intmap), GMScalarDerivative(F, basis1.ancestor(), basis2.ancestor(), 1, intmap)};
    // JF must be the jacobian (change of variable 3D->2D) of a scalar basis on the face (underlying basis1, perhaps after
    // several derived bases). We create one to make sure we capture the correct jacobian
    MonomialScalarBasisFace mono_P0_F(F, 0);
    Eigen::Matrix<double, 2, 3> JF = mono_P0_F.jacobian();
    VectorRd nF = F.normal();
    Eigen::Matrix<double, 2, 3> generators2 = basis2.generators();
    Eigen::Matrix2d W;
    W.col(0) = JF * nF.cross(generators2.row(0));
    W.col(1) = JF * nF.cross(generators2.row(1));
    
    for (size_t k=0; k<2; k++) {
      for (size_t j=0; j<2; j++) {
        gm.block(0, j*dim2, dim1, dim2) += W(k,j) * gmd[k];
      }
    }
    
    return gm;
  }

/// Generic template to compute the Gram Matrix of a Tangent basis and a Curl basis
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(
                    const Face & F,                             ///< Face to which the basis corresponds
                    const TangentFamily<BasisType1> & basis1,   ///< First basis
                    const CurlBasis<BasisType2> & basis2,       ///< Second basis
                    MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    return GramMatrix(F, basis2, basis1, mono_int_map).transpose();
  }


/// Generic template to compute the Gram Matrix of a Divergence basis and any other basis
template<typename BasisType1, typename BasisType2>
typename boost::disable_if<boost::is_same<BasisType2, MonomialFaceIntegralsType>, Eigen::MatrixXd>::type GramMatrix(
                     const Face& F,                          ///< Face to which the basis corresponds
                     const DivergenceBasis<BasisType1> & basis1,   ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,              ///< Second basis (scalar basis, columns of the Gram matrix)
                     MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    return GramMatrixDiv(F, basis1.ancestor(), basis2, mono_int_map);
  };
  
/// Template to compute the Gram Matrix of any basis and a Divergence basis
template<typename BasisType1, typename Basis2>
Eigen::MatrixXd GramMatrix(
                    const Face & F,                         ///< Face to which the basis corresponds
                    const BasisType1 & basis1,                  ///< First basis (scalar basis)
                    const DivergenceBasis<Basis2> & basis2,       ///< Second basis (divergence basis)
                    MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    return GramMatrixDiv(F, basis2.ancestor(), basis1, mono_int_map).transpose();
  };
  
/// Template to compute the Gram Matrix of a Divergence<TangentFamily> basis and a monomial scalar basis
/* This will work only if the generators of the TangentFamily are an orthonormal version of the Jacobian (change of variable 3D->2D) on the face */
template<typename BasisType1>
Eigen::MatrixXd GramMatrixDiv(
                    const Face & F,                                         ///< Face to which the basis corresponds
                    const TangentFamily<BasisType1> & basis1,   ///< First basis (divergence basis)
                    const MonomialScalarBasisFace & basis2,                 ///< Second basis (monomial scalar basis)
                    MonomialFaceIntegralsType mono_int_map = {}                 ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1, dim2);

    // Integrals of monomials
    size_t totaldegree = basis1.max_degree()+basis2.max_degree()-1;
    MonomialFaceIntegralsType intmap = CheckIntegralsDegree(F, totaldegree, mono_int_map);

    for (size_t i = 0; i < 2; i++) {
      gm.block(i*dim1/2, 0, dim1/2, dim2) = GMScalarDerivative(F, basis1.ancestor(), basis2, i, intmap);
    }
    
    return gm/F.diam();
  };
  
/// Computes the Gram Matrix of a Divergence<RolyCompl> basis and a monomial scalar basis
Eigen::MatrixXd GramMatrixDiv(
                    const Face & F,                           ///< Cell to which the basis corresponds
                    const RolyComplBasisFace & basis1,        ///< First basis (RolyCompl basis)
                    const MonomialScalarBasisFace & basis2,   ///< Second basis (scalar basis)
                    MonomialFaceIntegralsType mono_int_map = {}   ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );
  
/// Template to compute the Gram Matrix of the divergence of any basis and any other basis
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrixDiv(
                     const Face& F, ///< Face to which the basis corresponds
                     const BasisType1 & basis1, ///< First basis (vector basis)
                     const BasisType2 & basis2,  ///< Second basis (scalar basis)
                     MonomialFaceIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(useAncestor<BasisType1>() || useAncestor<BasisType2>(), "No method to compute this Gram matrix");
      
    if constexpr (!useAncestor<BasisType1>() && useAncestor<BasisType2>()) {
      return transformGM(basis2, 'C', GramMatrixDiv(F, basis1, basis2.ancestor(), mono_int_map) );
    } else if constexpr (useAncestor<BasisType1>() && !useAncestor<BasisType2>()) {
      return transformGM(basis1, 'R', GramMatrixDiv(F, basis1.ancestor(), basis2, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GramMatrixDiv(F, basis1.ancestor(), basis2.ancestor(), mono_int_map) ) );
    }
  };




/*@}*/
} // end of namespace HArDCore3D

#endif // end of _GMPOLY_FACE_HPP
