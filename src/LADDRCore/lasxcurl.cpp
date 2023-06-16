#include <lasxcurl.hpp>
#include <parallel_for.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

LASXCurl::LASXCurl(const LieAlgebra & lie_algebra, const SXCurl & sxcurl, bool use_threads, std::ostream & output)
  : VariableDOFSpace(
	     sxcurl.mesh(),
	     0,
	     PolynomialSpaceDimension<Edge>::Poly(sxcurl.degree()) * lie_algebra.dimension(), 
	     sxcurl.serPro().nDOFs_faces_SXCurl() * lie_algebra.dimension(),
	     sxcurl.serPro().nDOFs_cells_SXCurl() * lie_algebra.dimension()    
	     ),
    m_lie_algebra(lie_algebra),
    m_sxcurl(sxcurl),
    m_face_transfer_operators(sxcurl.mesh().n_faces()),
    m_cell_transfer_operators(sxcurl.mesh().n_cells()),
    m_use_threads(use_threads),
    m_output(output)
{
  output << "[LASXCurl] Initializing" << std::endl;
  if (m_use_threads) {
    m_output << "[LASXCurl] Parallel execution" << std::endl;
  } else {
    m_output << "[LASXCurl] Sequential execution" << std::endl;
  }

  // Construct face operators
  std::function<void(size_t, size_t)> construct_all_face_transfer_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iF = start; iF < end; iF++) {
          m_face_transfer_operators[iF].reset( new TransferOperators(_compute_face_transfer_operators(iF)) );
        } // for iF
      };

  m_output << "[LASXCurl] Constructing face transfer operators" << std::endl;
  parallel_for(mesh().n_faces(), construct_all_face_transfer_operators, m_use_threads);

  // Construct cell operators
  std::function<void(size_t, size_t)> construct_all_cell_extensions_reductions
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_cell_transfer_operators[iT].reset( new TransferOperators(_compute_cell_transfer_operators(iT)) );
        } // for iT
      };

  m_output << "[LASXCurl] Constructing cell transfer operators" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_extensions_reductions, m_use_threads);  
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd LASXCurl::interpolate(const LAFunctionType & v, const int doe_cell, const int doe_face, const int doe_edge) const
{
  // The matrix where the rows are storing the interpolate of each function of v in XCurl. This lines up in the columns the DOFs in LAXCurl that are the Lie algebra components of the same polynomial DOF in XCurl. 
  Eigen::MatrixXd v_M = Eigen::MatrixXd(m_lie_algebra.dimension(), m_sxcurl.dimension());

  for (size_t i = 0; i < m_lie_algebra.dimension(); i++){
    v_M.row(i) = m_sxcurl.interpolate(v[i], doe_cell, doe_face, doe_edge);
  }
  // To construct the element of LAXCurl, we read the matrix v_M column by column.
  Eigen::VectorXd vh(Eigen::Map<Eigen::VectorXd>(v_M.data(), v_M.size()));
  
  return vh;
}

//------------------------------------------------------------------------------
// Transfer operators
//------------------------------------------------------------------------------

LASXCurl::TransferOperators LASXCurl::_compute_face_transfer_operators(size_t iF)
{
  Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
  Eigen::MatrixXd SerF = Eigen::KroneckerProduct(m_sxcurl.ScurlFace(iF), id);
  Eigen::MatrixXd ExtF = Eigen::KroneckerProduct(m_sxcurl.EcurlFace(iF), id);
  Eigen::MatrixXd RedF = Eigen::KroneckerProduct(m_sxcurl.RcurlFace(iF), id);
  return TransferOperators(SerF, ExtF, RedF);
}


LASXCurl::TransferOperators LASXCurl::_compute_cell_transfer_operators(size_t iT)
{  
  Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_lie_algebra.dimension(), m_lie_algebra.dimension());
  Eigen::MatrixXd SerT = Eigen::KroneckerProduct(m_sxcurl.ScurlCell(iT), id);
  Eigen::MatrixXd ExtT = Eigen::KroneckerProduct(m_sxcurl.EcurlCell(iT), id);
  Eigen::MatrixXd RedT = Eigen::KroneckerProduct(m_sxcurl.RcurlCell(iT), id);
  return TransferOperators(SerT, ExtT, RedT);
}

//------------- L2 product with gradient ----------------//

Eigen::MatrixXd LASXCurl::computeL2ProductGradient(
                                     const size_t iT,
                                     const SXGrad & sx_grad,
                                     const std::string & side,
                                     const double & penalty_factor,
                                     const Eigen::MatrixXd & mass_Pk3_T,
                                     const IntegralWeight & weight
                                     ) const
{
  Eigen::MatrixXd L2_sxcurl_grad = m_sxcurl.computeL2ProductGradient(iT, sx_grad, side, penalty_factor, mass_Pk3_T, weight);

  return Eigen::KroneckerProduct(L2_sxcurl_grad, m_lie_algebra.massMatrix());
}

