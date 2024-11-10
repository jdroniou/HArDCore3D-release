#include "exterior_cell.hpp"
#include "exterior_evaluation.hpp"

using namespace HArDCore3D;
using namespace Manicore;

template<size_t space_dim,size_t object_dim>
Cell_basis<space_dim,object_dim>::Cell_basis(const MeshND::MeshObject<space_dim,object_dim> & f, int r)
  : _r(r),
    _scale(1/f.diam()),
    _tr(Cell_basis<space_dim,object_dim>::_init_trace(f)),
    _basis(_scale, -_scale*_tr*f.center_mass(),_tr,r) 
{
  if constexpr(space_dim == object_dim) {
    for (size_t i = 0; i <= object_dim; ++i) {
      _exterior_l2[i] = std::pow(_scale,2*i) * Eigen::MatrixXd::Identity(Dimension::ExtDim(i,object_dim),Dimension::ExtDim(i,object_dim));
    }
  } else {
    _cmp_pb_recurse<0>();
  }
}

template<size_t space_dim,size_t object_dim>
template<size_t l>
void Cell_basis<space_dim,object_dim>::_cmp_pb_recurse () {
  if constexpr (l <= object_dim) {
    auto ext_pb = Compute_pullback<l,space_dim,object_dim>::compute(_tr*_scale);
    _exterior_l2[l] = ext_pb.transpose()*ext_pb;
    _cmp_pb_recurse<l+1>();
  }
}

  ///------------------------------------------------------------------------------------------------------------------------------
  // Specialization of the trace for all dimensions
  ///------------------------------------------------------------------------------------------------------------------------------
template<>
Eigen::Matrix<double,3,3> 
Cell_basis<3,3>::_init_trace(const MeshND::MeshObject<3,3> & f) {
  return Eigen::Matrix<double,3,3>::Identity();
}
template<>
Eigen::Matrix<double,2,3> 
Cell_basis<3,2>::_init_trace(const MeshND::MeshObject<3,2> & f) {
  return f.face_tangent().transpose();
}
template<>
Eigen::Matrix<double,1,3> 
Cell_basis<3,1>::_init_trace(const MeshND::MeshObject<3,1> & f) {
  return f.tangent().transpose();
}
template<>
Eigen::Matrix<double,2,2> 
Cell_basis<2,2>::_init_trace(const MeshND::MeshObject<2,2> & f) {
  return Eigen::Matrix<double,2,2>::Identity();
}
template<>
Eigen::Matrix<double,1,2> 
Cell_basis<2,1>::_init_trace(const MeshND::MeshObject<2,1> & f) {
  return f.tangent().transpose();
}
  ///------------------------------------------------------------------------------------------------------------------------------

// 3D
template class HArDCore3D::Cell_basis<3,3>;
template class HArDCore3D::Cell_basis<3,2>;
template class HArDCore3D::Cell_basis<3,1>;
// 2D
//template class HArDCore3D::Cell_basis<2,2>;
//template class HArDCore3D::Cell_basis<2,1>;

