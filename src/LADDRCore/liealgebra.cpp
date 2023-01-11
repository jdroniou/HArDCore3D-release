#include "liealgebra.hpp"
#include <iostream>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructors
//------------------------------------------------------------------------------

LieAlgebra::LieAlgebra() 
{
  std::cout << "[LieAlgebra] Initializing default Lie algebra su(2)" << std::endl;
  // Default Lie algebra su(2)
  Eigen::Matrix2cd e1, e2, e3;
  e1 << 0, 1, 1, 0;
  e2 << 0, std::complex<double>{0, -1}, std::complex<double>{0, 1}, 0;
  e3 << 1, 0, 0, -1;
  std::complex<double> coeff {0, -0.5};
  m_basis = {coeff*e1, coeff*e2, coeff*e3};
  m_product = [](Eigen::MatrixXcd & A, Eigen::MatrixXcd & B) -> double {return -2.*(A*B).trace().real();};
  m_mass_matrix = _compute_mass_matrix();
  m_strucConst = _compute_structure_constants();
}

LieAlgebra::LieAlgebra(std::vector<LieAlgValue> & basis, LieAlgProduct & product)
  : m_basis(basis),
    m_product(product),
    m_mass_matrix(_compute_mass_matrix()),
    m_strucConst(_compute_structure_constants())
{
  std::cout << "[LieAlgebra] Initializing" << std::endl;
  Eigen::MatrixXd id = Eigen::MatrixXd::Identity(dimension(), dimension());
  if ((m_mass_matrix - id).norm() > 1e-12)
  {
    std::cout << "[LieAlgebra] Basis in Lie algebra is not orthonormal" << std::endl;
    exit(1);
  }
}

Eigen::MatrixXd LieAlgebra::_compute_mass_matrix()
{
  Eigen::MatrixXd M(dimension(), dimension());
  for (size_t i = 0; i < dimension(); i++){
    for (size_t j = 0; j <= i; j++){
      M(i, j) = M(j, i) = m_product(m_basis[i], m_basis[j]);
    }
  }
  return M;
}

std::vector<Eigen::MatrixXd> LieAlgebra::_compute_structure_constants()
{
  std::vector<Eigen::MatrixXd> struc_consts;
  for (size_t i = 0; i < dimension(); i++){
    struc_consts.emplace_back(Eigen::MatrixXd::Zero(dimension(), dimension()));
  }
  for (size_t i = 0; i < dimension(); i++){
    for (size_t j = 0; j < i; j++){
      LieAlgValue br_ij = lieBracket(m_basis[i], m_basis[j]);
      for (size_t k = 0; k < dimension(); k++){
        struc_consts[k](i, j) = m_product(br_ij, m_basis[k]);
        struc_consts[k](j, i) = -struc_consts[k](i, j);
      }// for k
    }// for j
  }// for i
  return struc_consts;
}
