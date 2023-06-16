#ifndef INTEGRATE_TRIPLE_PRODUCT_HPP
#define INTEGRATE_TRIPLE_PRODUCT_HPP

#include <GMpoly_cell.hpp>

namespace HArDCore3D {

    /*!
    *	@addtogroup Quadratures
    * @{
    */

typedef boost::multi_array<double, 3> Scalar3Tensor;
typedef Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>> MatrixSlice;
typedef Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Eigen::InnerStride<Eigen::Dynamic>> VectorSlice;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//  FUNCTIONS TO SLICE SCALAR 3_TENSORS
//----------------------------------------------------------------------
//----------------------------------------------------------------------

//
/// Function to slice a 3-tensor with respect to one index (returns a 2-tensor)
/** Slicing is always in order e.g. fixing y=5, returns matrix[i][k] = tensor[i][5][k].
The output is an Eigen::Map of the tensor (so accessing the values of the 2-tensor actually
access the values of the original 3-tensor) */
MatrixSlice slice(
                 Scalar3Tensor & tensor,  ///< Tensor to slice
                 size_t fixed_dim,  ///< Position of the index with respect to which we slice (fix the value)
                 size_t index   ///< Value of the index where the slicing is done
                 )
{
  size_t size_dim1 = tensor.shape()[0];
  size_t size_dim2 = tensor.shape()[1];
  size_t size_dim3 = tensor.shape()[2];

  size_t start;
  size_t num_rows;
  size_t num_cols;
  size_t column_stride; // no. elements between two column entries
  size_t row_stride;  // no. elements between two row entries

  if (fixed_dim == 1 && index < size_dim1){
    start = index*size_dim2*size_dim3;
    num_rows = size_dim2;
    num_cols = size_dim3;
    column_stride = 1;
    row_stride = size_dim3;
  } else if (fixed_dim == 2 && index < size_dim2){
    start = index*size_dim3;
    num_rows = size_dim1;
    num_cols = size_dim3;
    column_stride = 1;
    row_stride = size_dim2*size_dim3;
  } else if (fixed_dim == 3 && index < size_dim3){
    // Note: this map can be slow in calculations possibly due to the the spacing of the data (non-contiguous in either index), try using a copy if tests are slow.
    start = index;
    num_rows = size_dim1;
    num_cols = size_dim2;
    column_stride = size_dim3;
    row_stride = size_dim2*size_dim3;
  } else{
    std::cerr << "[IntegrateTripleProduct] ERROR: Multiarray slicing dimension or index out of range" << std::endl;
    exit(1);
  }
return MatrixSlice(tensor.data() + start, num_rows, num_cols, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(column_stride, row_stride));
};

/// Function to slice a 3-tensor with respect to two indices (returns a 1-tensor)
/** See comment in the other slice function */
VectorSlice slice(
                  Scalar3Tensor & tensor,  ///< Tensor to slice
                  size_t fixed_dim1,  ///< Position of the first index with respect to which we slice (fix the value)
                  size_t index1,      ///< Value of the first index to fix
                  size_t fixed_dim2,  ///< Position of the second index with respect to which we slice (fix the value)
                  size_t index2     ///< Value of the second index to fix
                  )
{
  size_t size_dim1 = tensor.shape()[0];
  size_t size_dim2 = tensor.shape()[1];
  size_t size_dim3 = tensor.shape()[2];

  size_t start;
  size_t num_elements;
  size_t stride;  // no. elements between two entries

  if (fixed_dim1 == 1 && index1 < size_dim1 && fixed_dim2 == 2 && index2 < size_dim2){
    start = index1*size_dim2*size_dim3 + index2*size_dim3;
    num_elements = size_dim3;
    stride = 1;
  } else if (fixed_dim1 == 2 && index1 < size_dim2 && fixed_dim2 == 1 && index2 < size_dim1){
    start = index2*size_dim2*size_dim3 + index1*size_dim3;
    num_elements = size_dim3;
    stride = 1;
  } else if (fixed_dim1 == 1 && index1 < size_dim1 && fixed_dim2 == 3 && index2 < size_dim3){
    start = index1*size_dim2*size_dim3 + index2;
    num_elements = size_dim2;
    stride = size_dim3;
  } else if (fixed_dim1 == 3 && index1 < size_dim3 && fixed_dim2 == 1 && index2 < size_dim1){
    start = index1 + index2*size_dim2*size_dim3;
    num_elements = size_dim2;
    stride = size_dim3;
  } else if (fixed_dim1 == 2 && index1 < size_dim2 && fixed_dim2 == 3 && index2 < size_dim3){
    start = index1*size_dim3 + index2;
    num_elements = size_dim1;
    stride = size_dim2*size_dim3;
  } else if (fixed_dim1 == 3 && index1 < size_dim3 && fixed_dim2 == 2 && index2 < size_dim2) {
    start = index1 + index2*size_dim3;
    num_elements = size_dim1;
    stride = size_dim2*size_dim3;
  }else{
    std::cerr << "[IntegrateTripleProduct] ERROR: Multiarray slicing dimension or index out of range" << std::endl;
    exit(1);
  }
  return VectorSlice(tensor.data()+start, num_elements, Eigen::InnerStride<>(stride));
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//  FUNCTIONS TO COMPUTE TRIPLE INTEGRALS
//----------------------------------------------------------------------
//----------------------------------------------------------------------

//
/* BASIC ONES: Monomial, tensorized, and generic template */

/// Computes the triple integral product of a scalar times the dot product of two vectors - basis1(basis2 . basis2).
template<size_t N>
Scalar3Tensor tripleInt(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const MonomialScalarBasisCell & basis1, ///< First basis 
                     const TensorizedVectorFamily<MonomialScalarBasisCell, N> & basis2,  ///< Second basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();
    size_t anc_dim2 = basis2.ancestor().dimension();
 
    size_t totaldegree = basis1.max_degree()+basis2.ancestor().max_degree()+basis2.ancestor().max_degree();
    Scalar3Tensor scalar_integrals(boost::extents[dim1][anc_dim2][anc_dim2]);

    // Obtain integration data from IntegrateCellMonomials
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

    Scalar3Tensor integrals(boost::extents[dim1][dim2][dim2]);
    std::fill_n(integrals.data(), integrals.num_elements(), 0);

    for (size_t i = 0; i < dim1; i++) {
      for (size_t j = 0; j < anc_dim2; j++) {
        for (size_t k = 0; k <= j; k++) {
          scalar_integrals[i][j][k] = intmap.at(basis1.powers(i) + basis2.ancestor().powers(j) + basis2.ancestor().powers(k));
          scalar_integrals[i][k][j] = scalar_integrals[i][j][k];
        }
      }
    }

    for (size_t i = 0; i < N; i++){
      integrals[boost::indices[boost::multi_array_types::index_range(0,dim1)][boost::multi_array_types::index_range(i*anc_dim2,i*anc_dim2+anc_dim2)][boost::multi_array_types::index_range(i*anc_dim2,i*anc_dim2+anc_dim2)]] = scalar_integrals;    
    }

    return integrals;
  };

/// Computes the triple integral product of a vector basis dot the cross product of the second vector basis - tens_family1 . (tens_family2 x tens_family2).
template<size_t N>
Scalar3Tensor tripleInt(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const TensorizedVectorFamily<MonomialScalarBasisCell, N> & tens_family1, ///< First basis 
                     const TensorizedVectorFamily<MonomialScalarBasisCell, N> & tens_family2,  ///< Second basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = tens_family1.dimension();
    size_t anc_dim1 = tens_family1.ancestor().dimension();
    size_t dim2 = tens_family2.dimension();
    size_t anc_dim2 = tens_family2.ancestor().dimension();
    size_t totaldegree = tens_family1.ancestor().max_degree() + tens_family2.ancestor().max_degree() + tens_family2.ancestor().max_degree();

    // Obtain integration data from IntegrateCellMonomials
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

    Scalar3Tensor integrals(boost::extents[dim1][dim2][dim2]);
    std::fill_n(integrals.data(), integrals.num_elements(), 0);

    Eigen::Matrix3i id = Eigen::Matrix3i::Identity();

    // Loop over the basis vectors: e_i . (e_j x e_k)
    for (size_t i = 0; i < N; i++){
      for (size_t j = 0; j < N; j++){
        for (size_t k = 0; k < j; k++){
          int sign = id.col(i).dot(id.col(j).cross(id.col(k)));
          if (sign){
            // Loop over the monomials
            for (size_t sub_i = 0; sub_i < anc_dim1; sub_i++){
              for (size_t sub_j = 0; sub_j < anc_dim2; sub_j++){
                for (size_t sub_k = 0; sub_k < anc_dim2; sub_k++){
                  integrals[i*anc_dim1 + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k] = 
                        sign * intmap.at(
                                  tens_family1.ancestor().powers(sub_i)
                                  + tens_family2.ancestor().powers(sub_j)
                                  + tens_family2.ancestor().powers(sub_k)
                                        );
                  if (sub_j != sub_k){
                    integrals[i*anc_dim1 + sub_i][j*anc_dim2 + sub_k][k*anc_dim2 + sub_j] = integrals[i*anc_dim1 + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                  }
                  integrals[i*anc_dim1 + sub_i][k*anc_dim2 + sub_j][j*anc_dim2 + sub_k] = -integrals[i*anc_dim1 + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                  integrals[i*anc_dim1 + sub_i][k*anc_dim2 + sub_k][j*anc_dim2 + sub_j] = -integrals[i*anc_dim1 + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                }// for sub_k
              }// for sub_j
            }// for sub_i
          };
        }// for k
      }// for j
    }// for i
    return integrals;
  };

/// Computes the triple integral product of a vector basis dot the cross product of the second vector basis - grad_basis . (tens_family x tens_family).
template<size_t N>
Scalar3Tensor tripleInt(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const GradientBasis<MonomialScalarBasisCell> & grad_basis, ///< First basis 
                     const TensorizedVectorFamily<MonomialScalarBasisCell, N> & tens_family,  ///< Second basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = grad_basis.dimension();
    size_t dim2 = tens_family.dimension();
    size_t anc_dim2 = tens_family.ancestor().dimension();
    size_t totaldegree = grad_basis.ancestor().max_degree()-1 + tens_family.ancestor().max_degree() + tens_family.ancestor().max_degree();

    // Obtain integration data from IntegrateCellMonomials
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

    Scalar3Tensor integrals(boost::extents[dim1][dim2][dim2]);
    std::fill_n(integrals.data(), integrals.num_elements(), 0);

    Eigen::Matrix3i id = Eigen::Matrix3i::Identity();

    // Loop over the basis vectors: e_i . (e_j x e_k)
    for (size_t i = 0; i < N; i++){
      for (size_t j = 0; j < N; j++){
        for (size_t k = 0; k < j; k++){
          int sign = id.col(i).dot(id.col(j).cross(id.col(k)));
          if (sign){
            // Loop over the monomials
            for (size_t sub_i = 0; sub_i < dim1; sub_i++){
              if (grad_basis.ancestor().powers(sub_i)[i]>0){
                for (size_t sub_j = 0; sub_j < anc_dim2; sub_j++){
                  for (size_t sub_k = 0; sub_k < anc_dim2; sub_k++){
                    integrals[sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k] = 
                          sign * grad_basis.ancestor().powers(sub_i)[i] 
                                  * intmap.at(
                                      grad_basis.ancestor().powers(sub_i)-id.col(i) 
                                      + tens_family.ancestor().powers(sub_j) 
                                      + tens_family.ancestor().powers(sub_k)
                                            ) / T.diam();
                    if (sub_j != sub_k){
                      integrals[sub_i][j*anc_dim2 + sub_k][k*anc_dim2 + sub_j] = integrals[sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                    }
                    integrals[sub_i][k*anc_dim2 + sub_j][j*anc_dim2 + sub_k] = -integrals[sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                    integrals[sub_i][k*anc_dim2 + sub_k][j*anc_dim2 + sub_j] = -integrals[sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                  }// for sub_k
                }// for sub_j
              };
            }// for sub_i
          };
        }// for k
      }// for j
    }// for i
    return integrals;
  };

/// Computes the triple integral product of a vector basis dot the cross product of the second vector basis - golycompl_basis . (tens_family x tens_family).
template<size_t N>
Scalar3Tensor tripleInt(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const GolyComplBasisCell & golycompl_basis, ///< First basis 
                     const TensorizedVectorFamily<MonomialScalarBasisCell, N> & tens_family,  ///< Second basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = golycompl_basis.dimension();
    size_t dimPkmo = golycompl_basis.dimPkmo();
    size_t dimPkmo2 = dim1 - 2*dimPkmo;
    size_t dim2 = tens_family.dimension();
    
    size_t anc_dim2 = tens_family.ancestor().dimension();
    size_t totaldegree = golycompl_basis.max_degree()+ tens_family.ancestor().max_degree() + tens_family.ancestor().max_degree();

    // Obtain integration data from IntegrateCellMonomials
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

    Scalar3Tensor integrals(boost::extents[dim1][dim2][dim2]);
    std::fill_n(integrals.data(), integrals.num_elements(), 0);

    Eigen::Matrix3i id = Eigen::Matrix3i::Identity();

    // Loop over the basis vectors: (e_j x e_k)
    for (size_t j = 0; j < N; j++){
      for (size_t k = 0; k < N; k++){
        Eigen::Vector3i direction = id.col(j).cross(id.col(k));
        int sign = direction.sum();
          // Loop over the basis/monomials
          for (size_t sub_j = 0; sub_j < anc_dim2; sub_j++){
            for (size_t sub_k = 0; sub_k <= sub_j; sub_k++){
              if (direction==id.col(0)){
                // Second direction of golycompl
                for (size_t sub_i = 0; sub_i < dimPkmo; sub_i++){
                  integrals[dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k] = -sign * intmap.at(golycompl_basis.powers(dimPkmo + sub_i) + id.col(2) + tens_family.ancestor().powers(sub_j) + tens_family.ancestor().powers(sub_k));
                  integrals[dimPkmo + sub_i][j*anc_dim2 + sub_k][k*anc_dim2 + sub_j] = integrals[dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                  integrals[dimPkmo + sub_i][k*anc_dim2 + sub_j][j*anc_dim2 + sub_k] = -integrals[dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                  integrals[dimPkmo + sub_i][k*anc_dim2 + sub_k][j*anc_dim2 + sub_j] = -integrals[dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                }
                // Third direction of golycompl
                for (size_t sub_i = 0; sub_i < dimPkmo2; sub_i++){
                  integrals[2*dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k] = sign * intmap.at(golycompl_basis.powers(2*dimPkmo + sub_i) + id.col(1) + tens_family.ancestor().powers(sub_j) + tens_family.ancestor().powers(sub_k));
                  integrals[2*dimPkmo + sub_i][j*anc_dim2 + sub_k][k*anc_dim2 + sub_j] = integrals[2*dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                  integrals[2*dimPkmo + sub_i][k*anc_dim2 + sub_j][j*anc_dim2 + sub_k] = -integrals[2*dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                  integrals[2*dimPkmo + sub_i][k*anc_dim2 + sub_k][j*anc_dim2 + sub_j] = -integrals[2*dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                }
              }else if (direction==id.col(1)){
                // First direction of golycompl
                for (size_t sub_i = 0; sub_i < dimPkmo; sub_i++){
                  integrals[sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k] = sign * intmap.at(golycompl_basis.powers(sub_i)+ id.col(2) + tens_family.ancestor().powers(sub_j) + tens_family.ancestor().powers(sub_k));
                  integrals[sub_i][j*anc_dim2 + sub_k][k*anc_dim2 + sub_j] = integrals[sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                  integrals[sub_i][k*anc_dim2 + sub_j][j*anc_dim2 + sub_k] = -integrals[sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                  integrals[sub_i][k*anc_dim2 + sub_k][j*anc_dim2 + sub_j] = -integrals[sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                }
                // Third direction of golycompl
                for (size_t sub_i = 0; sub_i < dimPkmo2; sub_i++){
                  integrals[2*dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k] = -sign * intmap.at(golycompl_basis.powers(2*dimPkmo + sub_i) + id.col(0) + tens_family.ancestor().powers(sub_j) + tens_family.ancestor().powers(sub_k));
                  integrals[2*dimPkmo + sub_i][j*anc_dim2 + sub_k][k*anc_dim2 + sub_j] = integrals[2*dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                  integrals[2*dimPkmo + sub_i][k*anc_dim2 + sub_j][j*anc_dim2 + sub_k] = -integrals[2*dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                  integrals[2*dimPkmo + sub_i][k*anc_dim2 + sub_k][j*anc_dim2 + sub_j] = -integrals[2*dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                }
              }else if (direction==id.col(2)){
              // First direction of golycompl
                for (size_t sub_i = 0; sub_i < dimPkmo; sub_i++){
                  integrals[sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k] = -sign * intmap.at(golycompl_basis.powers(sub_i)+ id.col(1) + tens_family.ancestor().powers(sub_j) + tens_family.ancestor().powers(sub_k));
                  integrals[sub_i][j*anc_dim2 + sub_k][k*anc_dim2 + sub_j] = integrals[sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                  integrals[sub_i][k*anc_dim2 + sub_j][j*anc_dim2 + sub_k] = -integrals[sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                  integrals[sub_i][k*anc_dim2 + sub_k][j*anc_dim2 + sub_j] = -integrals[sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                }
                // Second direction of golycompl
                for (size_t sub_i = 0; sub_i < dimPkmo; sub_i++){
                  integrals[dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k] = sign * intmap.at(golycompl_basis.powers(dimPkmo + sub_i) + id.col(0) + tens_family.ancestor().powers(sub_j) + tens_family.ancestor().powers(sub_k));
                  integrals[dimPkmo + sub_i][j*anc_dim2 + sub_k][k*anc_dim2 + sub_j] = integrals[dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                  integrals[dimPkmo + sub_i][k*anc_dim2 + sub_j][j*anc_dim2 + sub_k] = -integrals[dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                  integrals[dimPkmo + sub_i][k*anc_dim2 + sub_k][j*anc_dim2 + sub_j] = -integrals[dimPkmo + sub_i][j*anc_dim2 + sub_j][k*anc_dim2 + sub_k];
                }
              } //else
            }// for sub_k
          }// for sub_j
      }// for k
    }// for j
    return integrals;
  };

//----------------------------------------------------------------------
// Family, Shift, Tensorized

template<typename ScalarBasisType1, typename ScalarBasisType2, size_t N>
Scalar3Tensor tripleInt(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const GradientBasis<ShiftedBasis<ScalarBasisType1>> & grad_shift_basis, ///< First basis 
                     const TensorizedVectorFamily<Family<ScalarBasisType2>, N> & basis2,  ///< Second basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim2 = basis2.dimension();
    size_t unshifted_dim1 = grad_shift_basis.ancestor().ancestor().dimension();
    size_t shift = grad_shift_basis.ancestor().shift();

    GradientBasis<ScalarBasisType1> grad_basis(grad_shift_basis.ancestor().ancestor());

    Scalar3Tensor anc_integrals = tripleInt(T, grad_basis, basis2, mono_int_map);

    Scalar3Tensor integrals = anc_integrals[boost::indices[boost::multi_array_types::index_range(shift,unshifted_dim1)][boost::multi_array_types::index_range(0,dim2)][boost::multi_array_types::index_range(0,dim2)]];

    return integrals;
  };

template<typename ScalarBasisType, typename BasisType>
Scalar3Tensor tripleInt(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const GradientBasis<ShiftedBasis<ScalarBasisType>> & grad_shift_basis, ///< First basis 
                     const BasisType & basis2,  ///< Second basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim2 = basis2.dimension();
    size_t unshifted_dim1 = grad_shift_basis.ancestor().ancestor().dimension();
    size_t shift = grad_shift_basis.ancestor().shift();

    GradientBasis<ScalarBasisType> grad_basis(grad_shift_basis.ancestor().ancestor());

    Scalar3Tensor anc_integrals = tripleInt(T, grad_basis, basis2, mono_int_map);

    Scalar3Tensor integrals = anc_integrals[boost::indices[boost::multi_array_types::index_range(shift,unshifted_dim1)][boost::multi_array_types::index_range(0,dim2)][boost::multi_array_types::index_range(0,dim2)]];

    return integrals;
  };

template<typename BasisType, typename ScalarBasisType, size_t N>
Scalar3Tensor tripleInt(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const BasisType & basis1, ///< First basis 
                     const TensorizedVectorFamily<Family<ScalarBasisType>, N> & basis2,  ///< Second basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();
    size_t anc_dim2 = basis2.ancestor().dimension();
    size_t anc_anc_dim2 = basis2.ancestor().ancestor().dimension();

    Eigen::MatrixXd M2 = basis2.ancestor().matrix();

    Scalar3Tensor integrals(boost::extents[dim1][dim2][dim2]);
    std::fill_n(integrals.data(), integrals.num_elements(), 0);
    
    TensorizedVectorFamily<ScalarBasisType, N> red_basis2(basis2.ancestor().ancestor());

    Scalar3Tensor tensor_integrals = tripleInt(T, basis1, red_basis2, mono_int_map);

    for (size_t j = 0; j < N; j++){
      for (size_t k = 0; k < N; k++){
        Scalar3Tensor anc_scalar_integrals = tensor_integrals[boost::indices[boost::multi_array_types::index_range(0,dim1)][boost::multi_array_types::index_range(j*anc_anc_dim2,(j+1)*anc_anc_dim2)][boost::multi_array_types::index_range(k*anc_anc_dim2,(k+1)*anc_anc_dim2)]];
        for (size_t i = 0; i < dim1; i++){
          MatrixSlice int_jk = slice(anc_scalar_integrals, 1, i);
          for (size_t sub_j = 0; sub_j < anc_dim2; sub_j++){
            for (size_t sub_k = 0; sub_k <= sub_j; sub_k++){
              integrals[i][anc_dim2*j + sub_j][anc_dim2*k + sub_k] = M2.row(sub_j) * int_jk * M2.row(sub_k).transpose();
              integrals[i][anc_dim2*j + sub_k][anc_dim2*k + sub_j] = integrals[i][anc_dim2*j + sub_j][anc_dim2*k + sub_k];
            }
          }
        }
      }
    }
    return integrals;
  };

template<typename BasisType, typename ScalarBasisType, size_t N>
Scalar3Tensor tripleInt(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const Family<BasisType> & basis1, ///< First basis 
                     const Family<TensorizedVectorFamily<ScalarBasisType, N>> & basis2,  ///< Second basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();
    size_t anc_dim1 = basis1.ancestor().dimension();
    size_t anc_dim2 = basis2.ancestor().dimension();

    Eigen::MatrixXd M1 = basis1.matrix();
    Eigen::MatrixXd M2 = basis2.matrix();

    Scalar3Tensor integrals(boost::extents[dim1][dim2][dim2]);
    std::fill_n(integrals.data(), integrals.num_elements(), 0);

    Scalar3Tensor anc_integrals = tripleInt(T, basis1.ancestor(), basis2.ancestor(), mono_int_map);
    
    for (size_t i = 0; i < dim1; i++){
      Eigen::MatrixXd int_jk = Eigen::MatrixXd::Zero(anc_dim2, anc_dim2);
      for (size_t l = 0; l < anc_dim1; l++){
        int_jk += M1(i, l) * slice(anc_integrals, 1, l);
      }
      
      for (size_t j = 0; j < dim2; j++){
        for (size_t k = 0; k < dim2; k++){
          integrals[i][j][k] = M2.row(j) * int_jk * M2.row(k).transpose();
        }
      }
    }
    return integrals;
  };

template<typename BasisType, size_t N>
Scalar3Tensor tripleInt(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const Family<BasisType> & basis1, ///< First basis 
                     const TensorizedVectorFamily<MonomialScalarBasisCell, N> & basis2,  ///< Second basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();
    size_t anc_dim1 = basis1.ancestor().dimension();

    Eigen::MatrixXd M1 = basis1.matrix();

    Scalar3Tensor integrals(boost::extents[dim1][dim2][dim2]);
    std::fill_n(integrals.data(), integrals.num_elements(), 0);

    Scalar3Tensor scalar_integrals = tripleInt(T, basis1.ancestor(), basis2, mono_int_map);

    for (size_t i = 0; i < dim1; i++){
      MatrixSlice int_jk = slice(integrals, 1, i);
      for (size_t l = 0; l < anc_dim1; l++){
        int_jk += M1(i, l) * slice(scalar_integrals, 1, l);
      }
    }
    return integrals;
  };

template<typename BasisType1, typename BasisType2>
Scalar3Tensor tripleInt(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const BasisType1 & basis1, ///< First basis 
                     const Family<BasisType2> & basis2,  ///< Second basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();

    Eigen::MatrixXd M2 = basis2.matrix();

    Scalar3Tensor integrals(boost::extents[dim1][dim2][dim2]);
    std::fill_n(integrals.data(), integrals.num_elements(), 0);

    Scalar3Tensor scalar_integrals = tripleInt(T, basis1, basis2.ancestor(), mono_int_map);

    for (size_t i = 0; i < dim1; i++){
      MatrixSlice int_jk = slice(scalar_integrals, 1, i);
      for (size_t j = 0; j < dim2; j++){
        for (size_t k = 0; k <= j; k++){
          integrals[i][j][k] = M2.row(j) * int_jk * M2.row(k).transpose();
          integrals[i][k][j] = integrals[i][j][k];
        }
      }
    }
    return integrals;
  };

template<size_t N, typename ScalarBasisType>
Scalar3Tensor tripleInt(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const TensorizedVectorFamily<Family<ScalarBasisType>, N> & tens_family1, ///< First basis 
                     const TensorizedVectorFamily<MonomialScalarBasisCell, N> & tens_family2,  ///< Second basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = tens_family1.dimension();
    size_t anc_dim1 = tens_family1.ancestor().dimension();
    size_t anc_anc_dim1 = tens_family1.ancestor().ancestor().dimension();
    size_t dim2 = tens_family2.dimension();
    size_t anc_dim2 = tens_family2.ancestor().dimension();

    Eigen::MatrixXd M1 = tens_family1.ancestor().matrix();

    Scalar3Tensor integrals(boost::extents[dim1][dim2][dim2]);
    std::fill_n(integrals.data(), integrals.num_elements(), 0);
    
    TensorizedVectorFamily<ScalarBasisType, N> red_basis1(tens_family1.ancestor().ancestor());

    Scalar3Tensor tensor_integrals = tripleInt(T, red_basis1, tens_family2, mono_int_map);
    Eigen::Matrix3i id = Eigen::Matrix3i::Identity();
    // Loop over the basis vectors: e_i . (e_j x e_k)
    for (size_t i = 0; i < N; i++){
      for (size_t j = 0; j < N; j++){
        for (size_t k = 0; k < j; k++){
          int sign = id.col(i).dot(id.col(j).cross(id.col(k)));
          if (sign){
            Scalar3Tensor anc_scalar_integrals = tensor_integrals[boost::indices[boost::multi_array_types::index_range(i*anc_anc_dim1,(i+1)*anc_anc_dim1)][boost::multi_array_types::index_range(j*anc_dim2,(j+1)*anc_dim2)][boost::multi_array_types::index_range(k*anc_dim2,(k+1)*anc_dim2)]];
            for (size_t sub_j = 0; sub_j < anc_dim2; sub_j++){
              Eigen::MatrixXd fam_int_jk = M1 * slice(anc_scalar_integrals, 2, sub_j);
              for (size_t sub_i = 0; sub_i < anc_dim1; sub_i++){
                for (size_t sub_k = 0; sub_k < anc_dim2; sub_k++){
                  integrals[i*anc_dim1+sub_i][j*anc_dim2+sub_j][k*anc_dim2+sub_k] = fam_int_jk(sub_i, sub_k);
                  integrals[i*anc_dim1+sub_i][k*anc_dim2+sub_j][j*anc_dim2+sub_k] = -fam_int_jk(sub_i, sub_k);
                }
              }
            }
          };
        }// for k
      }// for j
    }// for i
  return integrals;
  };
}

#endif
