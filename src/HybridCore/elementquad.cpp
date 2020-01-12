// Class to create and store values of cell and face basis functions on quadrature points
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include "elementquad.hpp"
#include <string>
#include <mesh.hpp>
#include <quadraturerule.hpp>
//#include <vector>
//#include <cell.hpp>
//#include <face.hpp>

using namespace HArDCore3D;

// Creation class

ElementQuad::ElementQuad(const HybridCore& hho, const size_t iT, const size_t doeT, const size_t doeF )
  : _hho(hho),
    _iT(iT),
    _doeT(doeT),
    _doeF(doeF) {
  // Create cell elements
  const Mesh* mesh = _hho.get_mesh_ptr();
  Cell* cell = mesh->cell(_iT);
  _quadT = generate_quadrature_rule(*cell, _doeT);
  _phiT_quadT = _hho.basis_quad("cell", _iT, _quadT, _hho.K()+1);
  _dphiT_quadT = _hho.grad_basis_quad(_iT, _quadT, _hho.K()+1);

  // Create face elements
  size_t nfaces = cell->n_faces();
  _quadF.resize(nfaces);
  _phiT_quadF.resize(nfaces);
  _phiF_quadF.resize(nfaces);
  _dphiT_quadF.resize(nfaces);
  for (size_t ilF = 0; ilF < nfaces ; ilF++){
    Face* face = cell->face(ilF);
    const size_t iF = face->global_index();
    _quadF[ilF] = generate_quadrature_rule(*face, _doeF);
    _phiT_quadF[ilF] = _hho.basis_quad("cell", _iT, _quadF[ilF], _hho.K()+1);
    _phiF_quadF[ilF] = _hho.basis_quad("face", iF, _quadF[ilF], _hho.K());
    _dphiT_quadF[ilF] = _hho.grad_basis_quad(_iT, _quadF[ilF], _hho.K()+1);
  }

}

//ElementQuad::~ElementQuad() {}

//------------------------------------------------------------------------------
// Class methods
//------------------------------------------------------------------------------

std::vector<Eigen::ArrayXXd> ElementQuad::get_vec_phiT_quadT(size_t degree) const
{
  const size_t n_scalar_basis = _hho.dim_Pcell(degree);
  size_t dim = _hho.get_mesh_ptr()->dim(); 
  std::vector<Eigen::ArrayXXd> vec_phiT_quadT(dim * n_scalar_basis, Eigen::ArrayXXd::Zero(dim, _quadT.size()));
  for (size_t r = 0; r < dim; r++) {
    for (size_t i = 0; i < n_scalar_basis; i++) {
      vec_phiT_quadT[r * n_scalar_basis + i].row(r) = _phiT_quadT[i].transpose();
    } // for i
  } // for r
  return vec_phiT_quadT;
}

std::vector<Eigen::ArrayXXd> ElementQuad::get_vec_phiT_quadF(size_t ilF, size_t degree) const
{
  const size_t n_scalar_basis = _hho.dim_Pcell(degree);
  size_t dim = _hho.get_mesh_ptr()->dim(); 
  std::vector<Eigen::ArrayXXd> vec_phiT_quadF(dim * n_scalar_basis, Eigen::ArrayXXd::Zero(dim, _quadF[ilF].size())); 
  for (size_t r = 0; r < dim; r++) {
    for (size_t i = 0; i < n_scalar_basis; i++) {
      vec_phiT_quadF[r * n_scalar_basis + i].row(r) = _phiT_quadF[ilF][i].transpose();
    } // for i
  } // for r
  return vec_phiT_quadF;
}

std::vector<Eigen::ArrayXXd> ElementQuad::get_vec_phiF_quadF(size_t ilF, size_t degree) const
{
  const size_t n_scalar_basis = _hho.dim_Pface(degree);
  size_t dim = _hho.get_mesh_ptr()->dim(); 
  std::vector<Eigen::ArrayXXd> vec_phiF_quadF(dim * n_scalar_basis, Eigen::ArrayXXd::Zero(dim, _quadF[ilF].size()));
  for (size_t r = 0; r < dim; r++) {
    for (size_t i = 0; i < n_scalar_basis; i++) {
      vec_phiF_quadF[r * n_scalar_basis + i].row(r) = _phiF_quadF[ilF][i].transpose();
    } // for i
  } // for r
  return vec_phiF_quadF;
}
