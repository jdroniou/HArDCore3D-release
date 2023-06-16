#include <lasxgrad.hpp>
#include <parallel_for.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

LASXGrad::LASXGrad(const LieAlgebra & lie_algebra, const SXGrad & sxgrad, bool use_threads, std::ostream & output)
  : VariableDOFSpace(
	     sxgrad.mesh(),
	     1*lie_algebra.dimension(),
	     PolynomialSpaceDimension<Edge>::Poly(sxgrad.degree()-1)*lie_algebra.dimension(), 
	     sxgrad.serPro().nDOFs_faces_SXGrad() * lie_algebra.dimension(),
	     sxgrad.serPro().nDOFs_cells_SXGrad() * lie_algebra.dimension()    
	     ),
    m_lie_algebra(lie_algebra),
    m_sxgrad(sxgrad),
    m_face_transfer_operators(sxgrad.mesh().n_faces()),
    m_cell_transfer_operators(sxgrad.mesh().n_cells()),
    m_use_threads(use_threads),
    m_output(output)
{
  output << "[LASXGrad] Initializing" << std::endl;
  if (m_use_threads) {
    m_output << "[LASXGrad] Parallel execution" << std::endl;
  } else {
    m_output << "[LASXGrad] Sequential execution" << std::endl;
  }

  // Construct face operators
  std::function<void(size_t, size_t)> construct_all_face_transfer_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iF = start; iF < end; iF++) {
          m_face_transfer_operators[iF].reset( new TransferOperators(_compute_face_transfer_operators(iF)) );
        } // for iF
      };

  m_output << "[LASXGrad] Constructing face transfer operators" << std::endl;
  parallel_for(mesh().n_faces(), construct_all_face_transfer_operators, m_use_threads);

  // Construct cell operators
  std::function<void(size_t, size_t)> construct_all_cell_extensions_reductions
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_cell_transfer_operators[iT].reset( new TransferOperators(_compute_cell_transfer_operators(iT)) );
        } // for iT
      };

  m_output << "[LASXGrad] Constructing cell transfer operators" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_extensions_reductions, m_use_threads);  
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd LASXGrad::interpolate(const LAFunctionType & q, const int doe_cell, const int doe_face, const int doe_edge) const
{

  Eigen::MatrixXd q_M = Eigen::MatrixXd(m_lie_algebra.dimension(), m_sxgrad.dimension());

  for (size_t i = 0; i < m_lie_algebra.dimension(); i++){
    q_M.row(i) = m_sxgrad.interpolate(q[i], doe_cell, doe_face, doe_edge);
  }

  Eigen::VectorXd qh(Eigen::Map<Eigen::VectorXd>(q_M.data(), q_M.size()));
  
  return qh;
}

//------------------------------------------------------------------------------
// Transfer operators
//------------------------------------------------------------------------------

LASXGrad::TransferOperators LASXGrad::_compute_face_transfer_operators(size_t iF)
{
  Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
  Eigen::MatrixXd SerF = Eigen::KroneckerProduct(m_sxgrad.SgradFace(iF), id);
  Eigen::MatrixXd ExtF = Eigen::KroneckerProduct(m_sxgrad.EgradFace(iF), id);
  Eigen::MatrixXd RedF = Eigen::KroneckerProduct(m_sxgrad.RgradFace(iF), id);
  return TransferOperators(SerF, ExtF, RedF);
}


LASXGrad::TransferOperators LASXGrad::_compute_cell_transfer_operators(size_t iT)
{
  Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
  Eigen::MatrixXd SerT = Eigen::KroneckerProduct(m_sxgrad.SgradCell(iT), id);
  Eigen::MatrixXd ExtT = Eigen::KroneckerProduct(m_sxgrad.EgradCell(iT), id);
  Eigen::MatrixXd RedT = Eigen::KroneckerProduct(m_sxgrad.RgradCell(iT), id);
  return TransferOperators(SerT, ExtT, RedT);
}


