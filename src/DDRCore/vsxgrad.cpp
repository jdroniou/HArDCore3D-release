#include <vsxgrad.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

VSXGrad::VSXGrad(const SerendipityProblem & sp, bool use_threads, std::ostream & output)
  : VariableDOFSpace(
	     sp.mesh(),
	     dimspace,
	     dimspace * PolynomialSpaceDimension<Edge>::Poly(sp.ddrCore().degree()-1), 
	     dimspace * sp.nDOFs_faces_SXGrad(),
	     dimspace * sp.nDOFs_cells_SXGrad()    
	     ),
    m_sxgrad(SXGrad(sp.ddrCore(), sp, use_threads)),
    m_cell_bases(mesh().n_cells()),
    m_face_bases(mesh().n_faces()),
    m_edge_bases(mesh().n_edges()),
    m_use_threads(use_threads),
    m_output(output)
{
  output << "[VSXGrad] Initializing" << std::endl;
  if (m_use_threads) {
    m_output << "[VSXGrad] Parallel execution" << std::endl;
  } else {
    m_output << "[VSXGrad] Sequential execution" << std::endl;
  }

  // Construct edge bases
  std::function<void(size_t, size_t)> construct_all_edge_bases
    = [this](size_t start, size_t end)->void
      {
        for (size_t iE = start; iE < end; iE++) {
          m_edge_bases[iE].reset( new EdgeBases(_compute_edge_bases(iE)) );
        } // for iE
      };

  m_output << "[VSXGrad] Constructing edge bases" << std::endl;
  parallel_for(mesh().n_edges(), construct_all_edge_bases, m_use_threads);

  // Construct face bases
  std::function<void(size_t, size_t)> construct_all_face_bases
    = [this](size_t start, size_t end)->void
      {
        for (size_t iF = start; iF < end; iF++) {
          m_face_bases[iF].reset( new FaceBases(_compute_face_bases(iF)) );
        } // for iF
      };

  m_output << "[VSXGrad] Constructing face bases" << std::endl;
  parallel_for(mesh().n_faces(), construct_all_face_bases, m_use_threads);

  // Construct cell bases
  std::function<void(size_t, size_t)> construct_all_cell_bases
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_cell_bases[iT].reset( new CellBases(_compute_cell_bases(iT)) );
        } // for iT
      };

  m_output << "[VSXGrad] Constructing cell bases" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_bases, m_use_threads);  
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd VSXGrad::interpolate(const FunctionType & q, const int doe_cell, const int doe_face, const int doe_edge, bool boundary_only) const
{
  Eigen::MatrixXd q_M = Eigen::MatrixXd(dimspace, m_sxgrad.dimension());

  for (size_t i = 0; i < dimspace; i++){
    std::function<double(const VectorRd &)> q_i = 
              [&q, i](const VectorRd & x)-> double { return q(x).dot(basisRd[i]); };
    q_M.row(i) = m_sxgrad.interpolate(q_i, doe_cell, doe_face, doe_edge, boundary_only);
  }

  Eigen::VectorXd qh(Eigen::Map<Eigen::VectorXd>(q_M.data(), q_M.size()));
  
  return qh;
}

//------------------------------------------------------------------------------
// Bases
//------------------------------------------------------------------------------
VSXGrad::EdgeBases VSXGrad::_compute_edge_bases(size_t iE)
{
  Eigen::MatrixXd id = Eigen::Matrix3d::Identity();
  
  EdgeBases bases_E;
  bases_E.Polyk3.reset( new Poly3BasisEdgeType(*m_sxgrad.edgeBases(iE).Polyk));
  bases_E.Polykpo3.reset( new Poly3BasisEdgeType(*m_sxgrad.edgeBases(iE).Polykpo));

  return bases_E;
}

VSXGrad::FaceBases VSXGrad::_compute_face_bases(size_t iF)
{
  Eigen::MatrixXd id = Eigen::Matrix3d::Identity();
  
  FaceBases bases_F;
  bases_F.Polykpo3.reset( new Poly3BasisFaceType(*m_sxgrad.faceBases(iF).Polykpo));
  if (m_sxgrad.degree()>0) bases_F.Polykmo3.reset( new Poly3BasisFaceType(*m_sxgrad.faceBases(iF).Polykmo));

  return bases_F;
}

VSXGrad::CellBases VSXGrad::_compute_cell_bases(size_t iT)
{
  Eigen::MatrixXd id = Eigen::Matrix3d::Identity();
  
  CellBases bases_T;
  bases_T.Polykpo3.reset( new Poly3BasisCellType(*m_sxgrad.cellBases(iT).Polykpo));
  bases_T.Polyk3x3.reset( new Poly3x3BasisCellType(*m_sxgrad.cellBases(iT).Polyk));
  if (m_sxgrad.degree()>0) bases_T.Polykmo3.reset( new Poly3BasisCellType(*m_sxgrad.cellBases(iT).Polykmo));

  return bases_T;
}

std::vector<double> VSXGrad::computeH1norms(const std::vector<Eigen::VectorXd> & list_dofs, const std::string typegrad) const
{
  size_t nb_vectors = list_dofs.size();
  std::vector<Eigen::VectorXd> local_sqnorms(nb_vectors, Eigen::VectorXd::Zero(mesh().n_cells()));

  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &list_dofs, &local_sqnorms, &nb_vectors, &typegrad](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *mesh().cell(iT);

        // Full or symmetric gradient
        Eigen::MatrixXd GT = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3x3->dimension(), dimensionCell(iT));
        if (typegrad == "symmetric"){
          GT = cellSymGradient(iT);
        }else{
          GT = cellGradient(iT);
        }

        // Contribution of cell T
	      MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*degree());
        Eigen::MatrixXd mass_Pk3x3_T = GramMatrix(T, *cellBases(iT).Polyk3x3, int_mono_2k);
        Eigen::MatrixXd ST = computeStabilisation(iT);
        for (size_t i=0; i<nb_vectors; i++){
          Eigen::VectorXd uT = restrict(T, list_dofs[i]);

          local_sqnorms[i](iT) += uT.transpose() * (GT.transpose() * mass_Pk3x3_T * GT + ST) * uT;
        }
        
      }
    };
  parallel_for(mesh().n_cells(), compute_local_squarednorms, m_use_threads);
  
  // Vector of outputs
  std::vector<double> list_norms(nb_vectors, 0.);
  for (size_t i=0; i<nb_vectors; i++){
    list_norms[i] = std::sqrt(std::abs(local_sqnorms[i].sum()));
  }
  
  return list_norms;

}


std::vector<VectorRd> VSXGrad::computeVertexValues(const Eigen::VectorXd & u) const
{
  std::vector<VectorRd> values(mesh().n_vertices(), VectorRd::Zero());
  
  // Value at each vertex obtained averaging the values from all the cells around
  for (Vertex * V : mesh().get_vertices()){
    size_t iV = V->global_index();
 
    for (Cell * T : V->get_cells()){
      size_t iT = T->global_index();
      Eigen::VectorXd pTuT = cellPotential(iT) * restrict(*T, u);
      
      for (size_t i=0; i < cellBases(iT).Polykpo3->dimension(); i++){
        values[iV] += pTuT(i) * cellBases(iT).Polykpo3->function(i, V->coords());
      }
    }
    
    values[iV] /= V->n_cells();

  }
  
  return values;
}                  

