// Core data structures and methods required to implement hybrid schemes in 3D (polynomial unknowns
// in the cells and on the faces, such as Hybrid High-order (HHO) schemes).
//
// Provides:
//  - Hybrid polynomial basis functions (on the cells and faces of the mesh)
//  - Generic routines to create quadrature nodes over cells and faces of the mesh
//  - Interpolation of general functions onto the HHO space
//  - Methods for integrating, evaluating, and computing norms of HHO solutions
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include "hybridcore.hpp"
#include <vertex.hpp>
#include <edge.hpp>
#include <quad3d.hpp>
#include <quad3d_face.hpp>
#include <Eigen/Dense>
#include <iostream>
using namespace HArDCore3D;

// Creation class

HybridCore::HybridCore(const Mesh* mesh_ptr, const size_t K, const size_t L, const std::string choice_basis)
  : _mesh_ptr(mesh_ptr),
    _K(K),
    _L(L),
    _Ldeg(std::max(_L,0)),
    _nlocal_cell_dofs(dim_Pcell(_Ldeg)),
    _nlocal_face_dofs(dim_Pface(K)),
    _nhighorder_dofs(dim_Pcell(K+1)),
    _ngradient_dofs(_nhighorder_dofs - 1),
    _ntotal_cell_dofs(_nlocal_cell_dofs * mesh_ptr->n_cells()),
    _ntotal_face_dofs(_nlocal_face_dofs * mesh_ptr->n_faces()),
    _ninternal_face_dofs(_nlocal_face_dofs * mesh_ptr->n_i_faces()),
    _nboundary_face_dofs(_nlocal_face_dofs * mesh_ptr->n_b_faces()),
    _ntotal_dofs(_ntotal_cell_dofs+_ntotal_face_dofs),
    _choice_basis(choice_basis),
    _cell_monomials(0),
    _cell_monomials_gradients(0),
    _cell_bases(0),
    _cell_gradients(0),
    _face_monomials(0),
    _face_bases(0),
    _M_cell_basis(0),
    _M_face_basis(0),
    _offset_doe(0)        {
  std::cout << "Construct HybridCore" << std::endl;
  // Initialise basis functions on cells
  _M_cell_basis.reserve(mesh_ptr->n_cells());
  _cell_monomials.reserve(mesh_ptr->n_cells());
  _cell_monomials_gradients.reserve(mesh_ptr->n_cells());
  _cell_bases.reserve(mesh_ptr->n_cells());
  _cell_gradients.reserve(mesh_ptr->n_cells());
  for (size_t iT = 0; iT < mesh_ptr->n_cells(); iT++) {
    auto cell_mon = create_cell_monomials(iT);
    _cell_monomials.push_back(cell_mon.first);
    _cell_monomials_gradients.push_back(cell_mon.second);

    auto c_basis = create_basis("cell", iT);
    _cell_bases.push_back(std::get<0>(c_basis));
    _cell_gradients.push_back(std::get<1>(c_basis));
    _M_cell_basis.push_back(std::get<2>(c_basis));

  }

  // Initialise basis function on faces
  _M_face_basis.reserve(mesh_ptr->n_faces());
  _face_monomials.reserve(mesh_ptr->n_faces());
  _face_bases.reserve(mesh_ptr->n_faces());
  for (size_t iF = 0; iF < mesh_ptr->n_faces(); iF++){
    auto face_mon = create_face_monomials(iF);
    _face_monomials.push_back(face_mon);

    auto f_basis = create_basis("face", iF);
    _face_bases.push_back(std::get<0>(f_basis));
    _M_face_basis.push_back(std::get<2>(f_basis));

  }
}

// -------------------------------------------------
// ------- Cell and face basis functions
// -------------------------------------------------

size_t HybridCore::dim_Pcell(const size_t m) const {
  return (m * (m + 1) * (2 * m + 1) + 9 * m * (m + 1) + 12 * (m + 1)) / 12;
}

size_t HybridCore::dim_Pface(const size_t m) const {
  return (m + 1) * (m + 2) / 2;
}

const HybridCore::cell_basis_type& HybridCore::cell_monomial(size_t iT, size_t i) const {
  assert(iT < _mesh_ptr->n_cells());
  assert(i < _cell_monomials[iT].size());
  return _cell_monomials[iT][i];
}

const HybridCore::face_basis_type& HybridCore::face_monomial(size_t iF, size_t i) const {
  assert(iF < _mesh_ptr->n_faces());
  assert(i < _face_monomials[iF].size());
  return _face_monomials[iF][i];
}

const HybridCore::cell_basis_type& HybridCore::cell_basis(size_t iT, size_t i) const {
  assert(iT < _mesh_ptr->n_cells());
  assert(i < _cell_bases[iT].size());
  return _cell_bases[iT][i];
}

const HybridCore::face_basis_type& HybridCore::face_basis(size_t iF, size_t i) const {
  assert(iF < _mesh_ptr->n_faces());
  assert(i < _face_bases[iF].size());
  return _face_bases[iF][i];
}

const HybridCore::cell_gradient_type& HybridCore::cell_monomials_gradient(size_t iT, size_t i) const {
  assert(iT < _mesh_ptr->n_cells());
  assert(i < _cell_gradients[iT].size());
  return _cell_monomials_gradients[iT][i];
}

const HybridCore::cell_gradient_type& HybridCore::cell_gradient(size_t iT, size_t i) const {
  assert(iT < _mesh_ptr->n_cells());
  assert(i < _cell_gradients[iT].size());
  return _cell_gradients[iT][i];
}

std::pair<std::vector<HybridCore::cell_basis_type>,
          std::vector<HybridCore::cell_gradient_type> >
HybridCore::create_cell_monomials(const size_t iT) const {
  std::vector<cell_basis_type> cell_monomials;
  std::vector<cell_gradient_type> cell_monomials_gradient;
  cell_monomials.reserve(dim_Pcell(_K+1));
  cell_monomials_gradient.reserve(dim_Pcell(_K+1));

  Vector3d xT = _mesh_ptr->cell(iT)->center_mass();
  double hT = _mesh_ptr->cell(iT)->diam();

  // Create monomial functions up to degree K+1 since we need a high-order
  // basis for the reconstruction operators
  for (size_t degree = 0; degree <= _K+1; degree++){
    for (size_t i = 0; i <= degree; i++) {
      for (size_t j = 0; i + j <= degree; j++) {
        size_t k = degree - i - j;

        cell_basis_type phi = [xT, i, j, k, hT](double x, double y, double z)->double {
                                return std::pow( (x - xT.x())/hT , i) * std::pow( (y - xT.y())/hT , j) * std::pow( (z - xT.z())/hT , k);
                              };

        cell_gradient_type dphi = [xT, i, j, k, hT](double x, double y, double z) -> Eigen::Vector3d {
                                    Eigen::Vector3d gradient;
                                    gradient(0) = i == 0 ? 0.0 : i * std::pow( (x - xT.x())/hT , i - 1) * std::pow( (y - xT.y())/hT , j) * std::pow( (z - xT.z())/hT , k) / hT;
                                    gradient(1) = j == 0 ? 0.0 : std::pow( (x - xT.x())/hT , i) * j * std::pow( (y - xT.y())/hT , j - 1) * std::pow( (z - xT.z())/hT , k) / hT;
                                    gradient(2) = k == 0 ? 0.0 : std::pow( (x - xT.x())/hT , i) * std::pow( (y - xT.y())/hT , j) * k * std::pow( (z - xT.z())/hT , k - 1) / hT;
                                    return gradient;
                                  };

        // Store basis function and gradient
        cell_monomials.push_back(std::move(phi));
        cell_monomials_gradient.push_back(std::move(dphi));

      }
    }
  }

  return std::make_pair(std::move(cell_monomials), std::move(cell_monomials_gradient));
}

std::vector<HybridCore::face_basis_type> HybridCore::create_face_monomials(const size_t iF) const {
  std::vector<face_basis_type> face_monomials;
  face_monomials.reserve(dim_Pface(_K));

  Face* face = _mesh_ptr->face(iF);
  auto xF = face->center_mass();
  auto hF = face->diam();

  // Create a coordinate system on the face, using the center of the face and the two vertices of its
  // first edge
  auto x0 = face->edge(0)->vertex(0)->coords();
  auto x1 = face->edge(0)->vertex(1)->coords();
  Eigen::Vector3d v0 = (xF - x0).normalized();
  Eigen::Vector3d v1 = (xF - x1).normalized();
  v1 = (v1 - v0 * v1.dot(v0)).normalized();    // Orthonormalise the coordinate system (v0,v1)

  for (size_t degree = 0; degree <= _K; degree++) {
    for (size_t i = 0; i <= degree; i++) {
      size_t j = degree - i;

      face_basis_type phi = [xF, v0, v1, i, j, hF](double x, double y, double z) -> double {
                              Eigen::Vector3d face_vec = Eigen::Vector3d(x, y, z) - xF;
                              double s = face_vec.dot(v0);
                              double t = face_vec.dot(v1);
                              return std::pow(s/hF, i) * std::pow(t/hF, j);
                            };

      face_monomials.push_back(std::move(phi));
    }
  }

  return face_monomials;
}

std::tuple<std::vector<HybridCore::cell_basis_type>, std::vector<HybridCore::cell_gradient_type>, Eigen::MatrixXd> 
HybridCore::create_basis(const std::string cellface, const size_t iTF){

  std::vector<cell_basis_type> basis;
  std::vector<cell_gradient_type> gradients;
  Eigen::MatrixXd B;
        
  if (_choice_basis == "Mon"){
    // If monomial basis functions, we just copy
    if (cellface == "cell"){
      basis = _cell_monomials[iTF];
      gradients = _cell_monomials_gradients[iTF];
      B = Eigen::MatrixXd::Identity(dim_Pcell(_K+1), dim_Pcell(_K+1));
    }else{
      basis = _face_monomials[iTF];
      B = Eigen::MatrixXd::Identity(dim_Pface(_K), dim_Pface(_K));
    }

  }  else if (_choice_basis == "ON"){
    // Orthonormal basis functions
        
    // Grab quadrature rules and existing basis, and compute it at quadrature nodes
    QuadratureRule quad;
    std::vector<cell_basis_type> monomials;

    size_t degree = (cellface == "cell" ? _K+1 : _K);
    size_t nbasis = (cellface == "cell" ? dim_Pcell(degree) : dim_Pface(degree));
    if (cellface == "cell") {
      quad = generate_quadrature_rule(*_mesh_ptr->cell(iTF), 2 * (_K+1));
      monomials = _cell_monomials[iTF];
    }else{
      quad = generate_quadrature_rule(*_mesh_ptr->face(iTF), 2 * _K);
      monomials = _face_monomials[iTF];
    }
    size_t nbq = quad.size();
    std::vector<Eigen::ArrayXd> phi_quad = basis_quad(cellface, iTF, quad, degree, "monomials");

    // Arrange quadrature weights as an array, to simplify formulas
    Eigen::ArrayXd weight = Eigen::ArrayXd::Zero(nbq);
    for (size_t iq = 0; iq < nbq; iq++){
      weight(iq) = quad[iq].w;
    }

    // Orthonormalisation process.
    // We simultaneously compute the coefficients of the i-th ON basis onto the 1st,...,i-th monomials,
    // and the values of the i-th ON basis on the quadrature points (phi_quad is replaced with these values as we go)
    // B is triangular inferior, with row i storing the coefficients of the ON basis functions i on the monomials
    B = Eigen::MatrixXd::Zero(nbasis, nbasis);

    // Normalise first monomial
    double norm = sqrt( (weight * phi_quad[0] * phi_quad[0]).sum() );
    phi_quad[0] /= norm;
    B(0,0) = 1/norm;
    for (size_t ib = 1; ib < nbasis; ib++){
      // coefs represent the coefficients of the ib-th ON function on the ON basis functions 0 to ib-1
      Eigen::RowVectorXd coefs = Eigen::RowVectorXd::Zero(ib);
      for (size_t pb = 0; pb < ib; pb++){
        coefs(pb) = - (weight * phi_quad[ib] * phi_quad[pb]).sum();
      }
      // orthogonalise the ib-th monomial
      for (size_t pb = 0; pb < ib; pb++){
        phi_quad[ib] += coefs(pb) * phi_quad[pb];
      }
      // normalise ib-th basis function
      double norm = sqrt( (weight * phi_quad[ib] * phi_quad[ib]).sum() );
      phi_quad[ib] /= norm;
      coefs /= norm;
      // Compute ib-th row of B.
      // B.topLeftCorner(ib, ib) contains the rows that represent the ON basis functions 0 to ib-1 on
      // the monomials. Multiplying on the left by coefs gives the coefficients of the ON basis
      // functions on the monomials 0 to ib-1.
      B.block(ib, 0, 1, ib) = coefs * B.topLeftCorner(ib, ib);
      B(ib, ib) = 1/norm;
    }

    // Create basis functions (and gradients for cells)
    basis.resize(nbasis);
    for (size_t i=0; i < nbasis; i++){
      Eigen::VectorXd rowi = B.row(i);
      basis[i] = [monomials, iTF, i, rowi](double x, double y, double z)->double {
                   double val = 0;
                   for (size_t j = 0; j < i+1; j++){
                     val += rowi(j) * monomials[j](x, y, z);
                   }
                   return val;
                 };
    }

    if (cellface == "cell") {
      // Gradients
      std::vector<cell_gradient_type> monomials_gradients = _cell_monomials_gradients[iTF];
      gradients.resize(nbasis);
      for (size_t i=0; i < nbasis; i++){
        Eigen::VectorXd rowi = B.row(i);
        gradients[i] = [monomials_gradients, iTF, i, rowi](double x, double y, double z)->Vector3d {
                         Eigen::Vector3d grad = Eigen::Vector3d::Zero();
                         for (size_t j = 0; j < i+1; j++){
                           grad += rowi(j) * monomials_gradients[j](x, y, z);
                         }
                         return grad;
                       };
      }
    }
        
  } else {
    std::cout << "choice_basis unknown: " << _choice_basis << "\n\n";
    exit(EXIT_FAILURE);
  }

  return std::make_tuple(std::move(basis), std::move(gradients), std::move(B));

}

// -------------------------------------------------------------------
// ------  Gram matrices for scalar and vector-valued functions
// -------------------------------------------------------------------

Eigen::MatrixXd HybridCore::gram_matrix(const std::vector<Eigen::ArrayXd>& f_quad, const std::vector<Eigen::ArrayXd>& g_quad, const size_t& nrows, const size_t& ncols, const QuadratureRule& quad, const bool& sym, std::vector<double> L2weight) const {

  Eigen::MatrixXd GGM = Eigen::MatrixXd::Zero(nrows, ncols);

  size_t nbq = quad.size();

  // Recast product of quadrature and L2weight into an Eigen::ArrayXd
  Eigen::ArrayXd quad_L2_weights = Eigen::ArrayXd::Zero(nbq);
  if (L2weight.size() == 0){
    for (size_t iqn = 0; iqn < nbq; iqn++){
      quad_L2_weights(iqn) = quad[iqn].w;
    }
  }else{
    for (size_t iqn = 0; iqn < nbq; iqn++){
      quad_L2_weights(iqn) = quad[iqn].w * L2weight[iqn];
    }
  }

  for (size_t i = 0; i < nrows; i++){
    size_t jcut = 0;
    if (sym) jcut = i;
    for (size_t j = 0; j < jcut; j++){
      GGM(i, j) = GGM(j, i);
    }
    for (size_t j = jcut; j < ncols; j++){
      // Integrate f_i * g_j
      // The products here are component-wise since the terms are Eigen::ArrayXd
      GGM(i, j) = (quad_L2_weights * f_quad[i] * g_quad[j]).sum();
    }
  }

  return GGM;
}

Eigen::MatrixXd HybridCore::gram_matrix(const std::vector<Eigen::ArrayXXd>& F_quad, const std::vector<Eigen::ArrayXXd>& G_quad, const size_t& nrows, const size_t& ncols, const QuadratureRule& quad, const bool& sym, std::vector<Eigen::Matrix3d> L2Weight) const {

  Eigen::MatrixXd GSM = Eigen::MatrixXd::Zero(nrows, ncols);

  size_t nbq = quad.size();
  for (size_t i = 0; i < nrows; i++){
    size_t jcut = 0;
    if (sym) jcut = i;
    for (size_t j = 0; j < jcut; j++){
      GSM(i, j) = GSM(j, i);
    }
    // Multiply F_i by quadrature weights and matrix L2Weight, if required
    Eigen::ArrayXXd WeightsTimesF_quad = Eigen::ArrayXXd::Zero(_mesh_ptr->dim(), nbq);
    if (L2Weight.size() == 0){
      for (size_t iqn = 0; iqn < nbq; iqn++){
        WeightsTimesF_quad.col(iqn) = quad[iqn].w * F_quad[i].col(iqn);
      }
    }else{
      for (size_t iqn = 0; iqn < nbq; iqn++){
        WeightsTimesF_quad.col(iqn) = quad[iqn].w * (L2Weight[iqn] * F_quad[i].col(iqn).matrix());
      }
    }
    for (size_t j = jcut; j < ncols; j++){
      // Integrate F_i * G_j
      GSM(i, j) = (WeightsTimesF_quad * G_quad[j]).sum();
    }
  }

  return GSM;
}

// ----------------------------------------------------------------
// ------- Basis functions and gradients on quadrature nodes
// ----------------------------------------------------------------

const std::vector<Eigen::ArrayXd> HybridCore::basis_quad(const std::string cellface, 
      const size_t iTF, 
      const QuadratureRule quad, 
      const size_t degree, 
      std::string type_basis) const
{
  size_t nbq = quad.size();
  size_t nbasis = (cellface == "cell" ? dim_Pcell(degree) : dim_Pface(degree));
  std::vector<Eigen::ArrayXd> phi_quad(nbasis, Eigen::ArrayXd::Zero(nbq));      
  // We first compute the values of monomials at quadrature points
  for (size_t i = 0; i < nbasis; i++){
    if (cellface == "cell"){
      const auto &phi_i = cell_monomial(iTF, i);
      for (size_t iqn = 0; iqn < nbq; iqn++){
        phi_quad[i](iqn) = phi_i(quad[iqn].x, quad[iqn].y, quad[iqn].z);
      }

    }else{
      const auto &phi_i = face_monomial(iTF, i);
      for (size_t iqn = 0; iqn < nbq; iqn++){
        phi_quad[i](iqn) = phi_i(quad[iqn].x, quad[iqn].y, quad[iqn].z);
      }
    }
  }

  // If we actually wanted the basis functions, we modify these values by using the matrix to change from
  // monomials to basis functions
  if (type_basis == "basis"){
    std::vector<Eigen::ArrayXd> phi_quad_new(nbasis, Eigen::ArrayXd::Zero(nbq));        
    Eigen::MatrixXd B;
    if (cellface == "cell"){
      B = _M_cell_basis[iTF];
    }else{
      B = _M_face_basis[iTF];
    }

    for (size_t i = 0; i < nbasis; i++){
      for (size_t j = 0; j < nbasis; j++){
        phi_quad_new[i] += B(i,j) * phi_quad[j];
      }
    }

    phi_quad = phi_quad_new;
  }

  return phi_quad;
}

const std::vector<Eigen::ArrayXXd> HybridCore::grad_basis_quad(const size_t iT, 
      const QuadratureRule quad, 
      const size_t degree, 
      std::string type_basis) const 
{
  size_t nbq = quad.size();
  size_t nbasis = dim_Pcell(degree);
  std::vector<Eigen::ArrayXXd> dphi_quad(nbasis, Eigen::ArrayXXd::Zero(_mesh_ptr->dim(), nbq)); 

  // No need to consider i=0 since the first basis function is constant
  // We start by computing for gradients of monomials, and if required on basis functions we then change the values
  for (size_t i = 1; i < nbasis; i++){
    auto &dphi_i = cell_monomials_gradient(iT, i);

    for (size_t iqn=0; iqn<nbq; iqn++){
      dphi_quad[i].col(iqn) = dphi_i(quad[iqn].x, quad[iqn].y, quad[iqn].z);
    }
  }             

  if (type_basis == "basis"){
    std::vector<Eigen::ArrayXXd> dphi_quad_new(nbasis, Eigen::ArrayXXd::Zero(_mesh_ptr->dim(), nbq));   
    Eigen::MatrixXd B = _M_cell_basis[iT];

    for (size_t i = 0; i < nbasis; i++){
      for (size_t j = 0; j < nbasis; j++){
        dphi_quad_new[i] += B(i,j) * dphi_quad[j];
      }
    }

    dphi_quad = dphi_quad_new;
  }

  return dphi_quad;
}

// -----------------------------------------------------------------
// ------- Weights to eliminate cell unknows if L=-1
// -----------------------------------------------------------------


Eigen::VectorXd HybridCore::compute_weights(size_t iT) const {
  Cell* iCell = _mesh_ptr->cell(iT);
  size_t nlocal_faces = iCell->n_faces();
  Eigen::VectorXd barycoefT = Eigen::VectorXd::Zero(nlocal_faces);

  // Rule degree 0: all coefficients identical
  //    barycoefT = (1/double(nlocal_faces)) * Eigen::VectorXd::Ones(nlocal_faces);

  // Rule degree 1: coefficient is |F|d_{TF}/(d|T|)
  for (size_t ilF=0; ilF < nlocal_faces; ilF++){
    Eigen::Vector3d normalTF = iCell->face_normal(ilF);
    Eigen::Vector3d xFxT = iCell->face(ilF)->center_mass() - iCell->center_mass();

    double dTF = xFxT.dot(normalTF); 

    barycoefT(ilF) = iCell->face(ilF)->measure() * dTF / (_mesh_ptr->dim() * iCell->measure());

  }

  return barycoefT;
}

// ----------------------------------------------------------
// ---------- Norms of discrete unknowns
// ----------------------------------------------------------

Eigen::VectorXd HybridCore::restr(const Eigen::VectorXd &Xh, size_t iT) const {
        
  Cell* cell = _mesh_ptr->cell(iT);
  size_t nfacesT = cell->n_faces();

  Eigen::VectorXd XTF = Eigen::VectorXd::Zero(_nlocal_cell_dofs + nfacesT * _nlocal_face_dofs);

  XTF.head(_nlocal_cell_dofs) = Xh.segment(iT * _nlocal_cell_dofs, _nlocal_cell_dofs);
  for (size_t ilF = 0; ilF < nfacesT; ilF++){
    size_t offset_F = _ntotal_cell_dofs + cell->face(ilF)->global_index() * _nlocal_face_dofs;
    XTF.segment(_nlocal_cell_dofs + ilF * _nlocal_face_dofs, _nlocal_face_dofs) =
      Xh.segment(offset_F, _nlocal_face_dofs);
  }

  return XTF;
}

double HybridCore::L2norm(const Eigen::VectorXd &Xh) const {
  double value = 0.0;
  for (size_t iT = 0; iT < _mesh_ptr->n_cells(); iT++) {
    // L2 norm computed using the mass matrix
    // Compute cell quadrature nodes and values of cell basis functions at these nodes
    Cell* cell = _mesh_ptr->cell(iT);
    QuadratureRule quadT = generate_quadrature_rule(*cell, 2*(_K+1));
    std::vector<Eigen::ArrayXd> phiT_quadT = basis_quad("cell", iT, quadT, _Ldeg);

    Eigen::MatrixXd MTT = gram_matrix(phiT_quadT, phiT_quadT, _nlocal_cell_dofs, _nlocal_cell_dofs, quadT, true);
    Eigen::VectorXd XT = Xh.segment(iT*_nlocal_cell_dofs,_nlocal_cell_dofs);

    value += XT.dot(MTT*XT);

  }
  return sqrt(value);
}

double HybridCore::H1norm(const Eigen::VectorXd &Xh) const {

  double value = 0.0;
  for (size_t iT = 0; iT < _mesh_ptr-> n_cells(); iT++) {
    Cell* cell = _mesh_ptr->cell(iT);
    size_t nfacesT = cell->n_faces();
    size_t nlocal_dofs = _nlocal_cell_dofs + nfacesT * _nlocal_face_dofs;

    // Local matrix of the bilinear form corresponding to the cell contribution in the discrete H1 norm
    Eigen::MatrixXd H1aT = Eigen::MatrixXd::Zero(nlocal_dofs, nlocal_dofs);

    // Compute cell quadrature nodes and values of gradients of cell basis functions
    QuadratureRule quadT = generate_quadrature_rule(*cell, 2*_Ldeg);
    std::vector<Eigen::ArrayXXd> dphiT_quadT = grad_basis_quad(iT, quadT, _Ldeg);

    // CELL CONTRIBUTION
    //
    // Stiffness matrix
    H1aT.topLeftCorner(_nlocal_cell_dofs, _nlocal_cell_dofs) = 
      gram_matrix(dphiT_quadT, dphiT_quadT, _nlocal_cell_dofs, _nlocal_cell_dofs, quadT, true);


    // FACES CONTRIBUTION               
    for (size_t ilF=0; ilF < nfacesT ; ilF++){
      Face* face = cell->face(ilF);
      size_t iF = face->global_index();
      size_t offset_F = _nlocal_cell_dofs + ilF * _nlocal_face_dofs;
      // Different scalings for boundary terms
      //                double dTF = _mesh_ptr->face(iF)->diam();
      double dTF = cell->measure() / _mesh_ptr->face(iF)->measure();
      // Face quadrature nodes and values of cell and face basis functions at these nodes
      QuadratureRule quadF = generate_quadrature_rule(*face, 2*_K+1);
      std::vector<Eigen::ArrayXd> phiT_quadF = basis_quad("cell", iT, quadF, _Ldeg);
      std::vector<Eigen::ArrayXd> phiF_quadF = basis_quad("face", iF, quadF, _K);

      // Face, face-cell and face-face Gram matrix on F
      Eigen::MatrixXd MFF = gram_matrix(phiF_quadF, phiF_quadF, _nlocal_face_dofs, _nlocal_face_dofs, quadF, true);
      Eigen::MatrixXd MFT = gram_matrix(phiF_quadF, phiT_quadF, _nlocal_face_dofs, _nlocal_cell_dofs, quadF, false);
      Eigen::MatrixXd MTT_on_F = gram_matrix(phiT_quadF, phiT_quadF, _nlocal_cell_dofs, _nlocal_cell_dofs, quadF, true);

      // Contribution of the face to the local bilinear form
      H1aT.block(offset_F, offset_F, _nlocal_face_dofs, _nlocal_face_dofs) += MFF / dTF;
      H1aT.block(offset_F, 0, _nlocal_face_dofs, _nlocal_cell_dofs) -= MFT / dTF;
      H1aT.block(0, offset_F, _nlocal_cell_dofs, _nlocal_face_dofs) -= MFT.transpose() / dTF;
      H1aT.topLeftCorner(_nlocal_cell_dofs, _nlocal_cell_dofs) += MTT_on_F / dTF;
    }

    Eigen::VectorXd XTF = restr(Xh, iT);
    value += XTF.transpose() * H1aT * XTF;
  }

  return sqrt(value);
}


double HybridCore::Linf_face(const Eigen::VectorXd &Xh) const {
  double value = 0.0;

  Eigen::VectorXd XF = Xh.segment(_mesh_ptr->n_cells()*_nlocal_cell_dofs,_mesh_ptr->n_faces()*_nlocal_face_dofs);

  for (size_t iF = _mesh_ptr->n_cells()*_nlocal_cell_dofs; iF < _mesh_ptr->n_faces()*_nlocal_face_dofs; iF++) {
    value = std::max(value, std::abs(XF(iF)));

  }
  return value;
}


// -----------------------------------------------------------
// --------- Evaluate discrete functions in cell or faces
// -----------------------------------------------------------

double HybridCore::evaluate_in_cell(const Eigen::VectorXd XTF, size_t iT, double x, double y, double z) const {
  double value = 0.0;
  for (size_t i = 0; i < _nlocal_cell_dofs; i++) {
    const auto &phi_i = cell_basis(iT, i);
    const size_t index = i + iT * _nlocal_cell_dofs;
    value += XTF(index) * phi_i(x,y,z);
  }
  return value;
}

double HybridCore::evaluate_in_face(const Eigen::VectorXd XTF, size_t iF, double x, double y, double z) const {
  double value = 0.0;
  for (size_t i = 0; i < _nlocal_face_dofs; i++) {
    const auto &phi_i = face_basis(iF, i);
    const size_t index = _ntotal_cell_dofs + iF * _nlocal_face_dofs + i;
    value += XTF(index) * phi_i(x,y,z);
  }
  return value;
}


// ---------------------------------------------------
// ------- Compute vertex values from a hybrid function
// ---------------------------------------------------

Eigen::VectorXd HybridCore::VertexValues(const Eigen::VectorXd Xh, const std::string from_dofs) {
  Eigen::VectorXd function_vertex = Eigen::VectorXd::Zero(_mesh_ptr->n_vertices());

  if (from_dofs == "cell"){
    // compute from the cell polynomials, averaging all values coming from the cells around each vertex
    for (size_t iV = 0; iV < _mesh_ptr->n_vertices(); iV++){
      auto xV = _mesh_ptr->vertex(iV)->coords();
      auto cList = _mesh_ptr->vertex(iV)->get_cells();
      auto weight = cList.size();
      for (size_t ilT = 0; ilT < weight; ilT++){
        auto temp = this->evaluate_in_cell(Xh, cList[ilT]->global_index(), xV.x(), xV.y(), xV.z());
        function_vertex(iV) += temp;
      }
      function_vertex(iV) = function_vertex(iV)/(weight);
    }

  }else{
    // compute from the face polynomials, averaging all values coming from the faces around each vertex
    for (size_t iV = 0; iV < _mesh_ptr->n_vertices(); iV++){
      auto xV = _mesh_ptr->vertex(iV)->coords();
      auto eList = _mesh_ptr->vertex(iV)->get_faces();
      double weight = eList.size();
      for (size_t ilF = 0; ilF < weight; ilF++){
        auto temp = this->evaluate_in_face(Xh, eList[ilF]->global_index(), xV.x(), xV.y(), xV.z());
        function_vertex(iV) += temp;
      }
      function_vertex(iV) = function_vertex(iV)/(weight);
    }
  }

  return function_vertex;
}


