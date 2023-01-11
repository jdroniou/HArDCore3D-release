#include <laxgrad.hpp>
#include <parallel_for.hpp>
#include <unsupported/Eigen/KroneckerProduct>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

LAXGrad::LAXGrad(const LieAlgebra & lie_algebra, const XGrad & xgrad, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(xgrad.mesh(),
	     1*lie_algebra.dimension(),
	     PolynomialSpaceDimension<Edge>::Poly(xgrad.degree() - 1)*lie_algebra.dimension(),
	     PolynomialSpaceDimension<Face>::Poly(xgrad.degree() - 1)*lie_algebra.dimension(),
	     PolynomialSpaceDimension<Cell>::Poly(xgrad.degree() - 1)*lie_algebra.dimension()
	     ),
    m_lie_algebra(lie_algebra),
    m_xgrad(xgrad),
    m_use_threads(use_threads),
    m_output(output),
    m_edge_operators(xgrad.mesh().n_edges()),
    m_face_operators(xgrad.mesh().n_faces()),
    m_cell_operators(xgrad.mesh().n_cells())
{

  output << "[LAXGrad] Initializing" << std::endl;
  if (m_use_threads) {
    m_output << "[LAXGrad] Parallel execution" << std::endl;
  } else {
    m_output << "[LAXGrad] Sequential execution" << std::endl;
  }

  // Construct edge gradients and potentials
  std::function<void(size_t, size_t)> construct_all_edge_gradients_potentials
    = [this](size_t start, size_t end)->void
      {
        for (size_t iE = start; iE < end; iE++) {
          m_edge_operators[iE].reset( new LocalOperators(_compute_edge_gradient_potential(iE)) );
        } // for iE
      };

  m_output << "[LAXGrad] Constructing edge gradients and potentials" << std::endl;
  parallel_for(mesh().n_edges(), construct_all_edge_gradients_potentials, m_use_threads);

  // Construct face gradients and potentials
  std::function<void(size_t, size_t)> construct_all_face_gradients_potentials
    = [this](size_t start, size_t end)->void
      {
        for (size_t iF = start; iF < end; iF++) {
          m_face_operators[iF].reset( new LocalOperators(_compute_face_gradient_potential(iF)) );
        } // for iF
      };

  m_output << "[LAXGrad] Constructing face gradients and potentials" << std::endl;
  parallel_for(mesh().n_faces(), construct_all_face_gradients_potentials, m_use_threads);

  // Construct cell gradients and potentials
  std::function<void(size_t, size_t)> construct_all_cell_gradients_potentials
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_cell_operators[iT].reset( new LocalOperators(_compute_cell_gradient_potential(iT)) );
        } // for iT
      };

  m_output << "[LAXGrad] Constructing cell gradients and potentials" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_gradients_potentials, m_use_threads);

}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd LAXGrad::interpolate(const LAFunctionType & q, const int doe_cell, const int doe_face, const int doe_edge) const
{

  Eigen::MatrixXd q_M(m_lie_algebra.dimension(), m_xgrad.dimension());

  for (size_t i = 0; i < m_lie_algebra.dimension(); i++){
    q_M.row(i) = m_xgrad.interpolate(q[i], doe_cell, doe_face, doe_edge);
  }

  Eigen::VectorXd qh(Eigen::Map<Eigen::VectorXd>(q_M.data(), q_M.size()));
  
  return qh;

}

//------------------------------------------------------------------------------
// Gradient and potential reconstructions
//------------------------------------------------------------------------------

LAXGrad::LocalOperators LAXGrad::_compute_edge_gradient_potential(size_t iE)
{

  Eigen::MatrixXd E_gradient = m_xgrad.edgeOperators(iE).gradient;
  Eigen::MatrixXd E_potential = m_xgrad.edgeOperators(iE).potential;
  
  Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
  
  return LocalOperators(Eigen::KroneckerProduct(E_gradient, id), Eigen::KroneckerProduct(E_potential, id));
}

//------------------------------------------------------------------------------

LAXGrad::LocalOperators LAXGrad::_compute_face_gradient_potential(size_t iF)
{
  
  Eigen::MatrixXd F_gradient = m_xgrad.faceOperators(iF).gradient;
  Eigen::MatrixXd F_potential = m_xgrad.faceOperators(iF).potential;
  
  Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
  
  return LocalOperators(Eigen::KroneckerProduct(F_gradient, id), Eigen::KroneckerProduct(F_potential, id));
}

//------------------------------------------------------------------------------

LAXGrad::LocalOperators LAXGrad::_compute_cell_gradient_potential(size_t iT)
{
  Eigen::MatrixXd T_gradient = m_xgrad.cellOperators(iT).gradient;
  Eigen::MatrixXd T_potential = m_xgrad.cellOperators(iT).potential;
  
  Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());   
                                        
  return LocalOperators(Eigen::KroneckerProduct(T_gradient, id), Eigen::KroneckerProduct(T_potential, id));
}

//------------------------------------------------------------------------------
// Local L2 products on LAXGrad
//------------------------------------------------------------------------------

Eigen::MatrixXd LAXGrad::computeL2Product(
                                        const size_t iT,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pkpo_T,
                                        const IntegralWeight & weight
                                        ) const
{

  Eigen::MatrixXd L2_xgrad = m_xgrad.computeL2Product(iT, penalty_factor, mass_Pkpo_T, weight);

  return Eigen::KroneckerProduct(L2_xgrad, m_lie_algebra.massMatrix());

}
