#include <serendipity_problem.hpp>
#include <mesh.hpp>
#include <parallel_for.hpp>
#include <GMpoly_face.hpp>
#include <GMpoly_cell.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------

SerendipityProblem::SerendipityProblem(const DDRCore & ddrcore, bool use_threads, std::ostream & output)
  : m_ddrcore(ddrcore),
    m_output(output),
    m_face_bases_Polyl(ddrcore.mesh().n_faces()),
    m_cell_bases_Polyl(ddrcore.mesh().n_cells()),
    m_face_bases_RolyCompllpo(ddrcore.mesh().n_faces()),
    m_cell_bases_RolyCompllpo(ddrcore.mesh().n_cells()),
    m_serendipity_edges(ddrcore.mesh().n_faces()),
    m_serendipity_faces(ddrcore.mesh().n_cells()),
    m_inverse_problem_faces(ddrcore.mesh().n_faces()),
    m_inverse_problem_cells(ddrcore.mesh().n_cells())
{
  m_output << "[SerendipityProblem] Constructing" << std::endl;
 
  // Construct serendipity on faces
  std::function<void(size_t, size_t)> construct_all_faces
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iF = start; iF < end; iF++) {
	        m_serendipity_edges[iF].reserve(mesh().face(iF)->n_edges());
	        m_serendipity_edges[iF] = _compute_serendipity_edges(iF);
	        if (dimFacePolyl(iF)>0){
  	        m_face_bases_Polyl[iF].reset( new RestrictedBasis<DDRCore::PolyBasisFaceType>(*m_ddrcore.faceBases(iF).Polyk, dimFacePolyl(iF)));
  	      }
	        if (dimFaceRolyCompllpo(iF)>0){
  	        m_face_bases_RolyCompllpo[iF].reset( new RestrictedBasis<DDRCore::RolyComplBasisFaceType>(*m_ddrcore.faceBases(iF).RolyComplk, dimFaceRolyCompllpo(iF)));
  	      }
	        m_inverse_problem_faces[iF] = _compute_inverse_problem_faces(iF);
	      } // for iF
      };      
  parallel_for(mesh().n_faces(), construct_all_faces, use_threads);

  // Construct serendipity on cells
  std::function<void(size_t, size_t)> construct_all_cells
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iT = start; iT < end; iT++) {
	        m_serendipity_faces[iT].reserve(mesh().cell(iT)->n_faces());
	        m_serendipity_faces[iT] = _compute_serendipity_faces(iT);
	        if (dimCellPolyl(iT)>0){
  	        m_cell_bases_Polyl[iT].reset( new RestrictedBasis<DDRCore::PolyBasisCellType>(*m_ddrcore.cellBases(iT).Polyk, dimCellPolyl(iT)));
  	      }
	        if (dimCellRolyCompllpo(iT)>0){
  	        m_cell_bases_RolyCompllpo[iT].reset( new RestrictedBasis<DDRCore::RolyComplBasisCellType>(*m_ddrcore.cellBases(iT).RolyComplk, dimCellRolyCompllpo(iT)));
  	      }
	        m_inverse_problem_cells[iT] = _compute_inverse_problem_cells(iT);
	      } // for iT
      };
  parallel_for(mesh().n_cells(), construct_all_cells, use_threads);

}

//------------------------------------------------------------------------------

std::vector<size_t> SerendipityProblem::_compute_serendipity_edges(size_t iF)
{
  const Face & F = *mesh().face(iF);
  std::vector<size_t> ser_edges;
  ser_edges.reserve(F.n_edges());

  // For an edge to be selected, it must not be aligned with previously selected edges,
  // and the face must lie only on one side of the edge. We check that by asking that
  // the signed distances between the edge midpoint and other edges in the list
  // are above a certain threshold
  // The selection of the threshold depends on how the face is skewed
  double skew = F.measure() / std::pow(F.diam(), 2);
  double threshold = 0.1 * skew * F.diam();
  for (size_t iE=0; iE<F.n_edges(); iE++){
    const Edge & E = *F.edge(iE);
    bool use_edge = true;
    
    // Check if not aligned with previously selected edges
    for (size_t iEp : ser_edges){
      const VectorRd outer_normal_Ep = F.edge_orientation(iEp) * F.edge_normal(iEp);
      double dist_iE_iEp = (F.edge(iEp)->center_mass() - E.center_mass()).dot(outer_normal_Ep); 
      if (dist_iE_iEp < threshold){
        use_edge = false;
        break;
      }
    }
    // Check if face on one side of the edge
    if (use_edge){
      const VectorRd outer_normal_E = F.edge_orientation(iE) * F.edge_normal(iE);
      for (Vertex * V : F.get_vertices()){
        if ( (E.center_mass() - V->coords()).dot(outer_normal_E) < -1e-3*threshold){
          use_edge = false;
          break;
       }      
      }
    }
  
    if (use_edge){
      ser_edges.push_back(iE);
    }

  }

  if (ser_edges.size()<2){
    std::cout << "[SerendipityProblem] ERROR: fewer than 2 serendipity edges in a face." << std::endl;
    exit(1);
  }

  return ser_edges;
}

//------------------------------------------------------------------------------

SerendipityProblem::InverseProblem SerendipityProblem::_compute_inverse_problem_faces(size_t iF)
{
  const Face & F = *mesh().face(iF);

  size_t dim_Pk2 = m_ddrcore.faceBases(iF).Polyk2->dimension();
  Eigen::MatrixXd AF = Eigen::MatrixXd::Zero(dim_Pk2 + dimFaceRolyCompllpo(iF), dim_Pk2 + dimFaceRolyCompllpo(iF));

  // Edge contributions
  for (Edge * E : F.get_edges()){
    QuadratureRule quad_2k_E = generate_quadrature_rule(*E, 2*m_ddrcore.degree());
    boost::multi_array<double, 2> Pk2_tE_quad 
              = scalar_product( evaluate_quad<Function>::compute( *m_ddrcore.faceBases(iF).Polyk2, quad_2k_E),
                               E->tangent()
                               );
    AF.topLeftCorner(dim_Pk2, dim_Pk2) += compute_gram_matrix(Pk2_tE_quad, quad_2k_E);
  }
  AF.topLeftCorner(dim_Pk2, dim_Pk2) *= F.diam();

  // rot_F contribution
  auto scalar_rot_Pk2F = ScalarRotFamily<DDRCore::PolyBasisFaceType>(*m_ddrcore.faceBases(iF).Polyk2, F);
  MonomialFaceIntegralsType int_monoF_2k = IntegrateFaceMonomials(F, 2*m_ddrcore.degree());
  //// Currently there's no method for GramMatrix of divbasis/divbasis, will have to implement that later
////  DivergenceBasis<TangentFamily<DDRCore::PolyBasisFaceType>> basis_rot_Pk2F(rotated_basis);
////  AF.topLeftCorner(dim_Pk2, dim_Pk2) += std::pow(F.diam(), 2) * GramMatrixDiv(F, basis_rot_Pk2F, int_monoF_2k);
  auto quad_2k_F = generate_quadrature_rule(F, 2*m_ddrcore.degree());
  AF.topLeftCorner(dim_Pk2, dim_Pk2) += std::pow(F.diam(), 2) * 
                      compute_gram_matrix(evaluate_quad<Function>::compute(scalar_rot_Pk2F, quad_2k_F), quad_2k_F);

  // Lagrange multiplier terms
  if (dimFaceRolyCompllpo(iF)>0){
    Eigen::MatrixXd Gram_Rclpo_Pk2 = GramMatrix(F, faceBasisRolyCompllpo(iF), *m_ddrcore.faceBases(iF).Polyk2, int_monoF_2k);
    AF.bottomLeftCorner(dimFaceRolyCompllpo(iF), dim_Pk2) += Gram_Rclpo_Pk2;
    AF.topRightCorner(dim_Pk2, dimFaceRolyCompllpo(iF)) -= Gram_Rclpo_Pk2.transpose();
  }
  
  return InverseProblem(AF);
}

//------------------------------------------------------------------------------
const Eigen::MatrixXd SerendipityProblem::SerendipityOperatorFace(const size_t iF, const Eigen::MatrixXd & LF) const
{
  return ( m_inverse_problem_faces[iF].solve(LF) ).topRows(m_ddrcore.faceBases(iF).Polyk2->dimension());
}


//------------------------------------------------------------------------------

std::vector<size_t> SerendipityProblem::_compute_serendipity_faces(size_t iT)
{
  const Cell & T = *mesh().cell(iT);
  std::vector<size_t> ser_faces;
  ser_faces.reserve(T.n_faces());

  double skew = T.measure() / std::pow(T.diam(), 3);
  double threshold = 0.1 * skew * T.diam();

  for (size_t iF=0; iF<T.n_faces(); iF++){
    const Face & F = *T.face(iF);
    bool use_face = true;
    
    // Check if not aligned with previously selected faces
    for (size_t iFp : ser_faces){
      double dist_iF_iFp = (T.face(iFp)->center_mass() - F.center_mass()).dot(T.face_normal(iFp)); 
      if (dist_iF_iFp < threshold){
        use_face = false;
        break;
      }
    }
    // Check if cell on one side of the face
    if (use_face){
      for (Vertex * V : T.get_vertices()){
        if ( (F.center_mass() - V->coords()).dot(T.face_normal(iF)) < -1e-3*threshold){
          use_face = false;
          break;
       }      
      }
    }
  
    if (use_face){
      ser_faces.push_back(iF);
    }

  }

  if (ser_faces.size()<2){
    std::cout << "[SerendipityProblem] ERROR: fewer than 2 serendipity faces in a cell." << std::endl;
    exit(1);
  }

  return ser_faces;
}

//------------------------------------------------------------------------------

SerendipityProblem::InverseProblem SerendipityProblem::_compute_inverse_problem_cells(size_t iT)
{
  const Cell & T = *mesh().cell(iT);

  size_t dim_Pk3 = m_ddrcore.cellBases(iT).Polyk3->dimension();
  Eigen::MatrixXd AT = Eigen::MatrixXd::Zero(dim_Pk3 + dimCellRolyCompllpo(iT), dim_Pk3 + dimCellRolyCompllpo(iT));

  // Face contributions
  for (Face * F : T.get_faces()){
    QuadratureRule quad_2k_F = generate_quadrature_rule(*F, 2*m_ddrcore.degree());
    const VectorRd nF = F->normal();
    boost::multi_array<VectorRd, 2> Pk3_tangentF_quad 
              = transform_values_quad<VectorRd>( evaluate_quad<Function>::compute( *m_ddrcore.cellBases(iT).Polyk3, quad_2k_F),
                               [&nF](const VectorRd &z)->VectorRd { return z-(z.dot(nF))*nF;}
                               );
    AT.topLeftCorner(dim_Pk3, dim_Pk3) += compute_gram_matrix(Pk3_tangentF_quad, quad_2k_F);
  }
  AT.topLeftCorner(dim_Pk3, dim_Pk3) *= T.diam();

  // curl curl contribution
  CurlBasis<DDRCore::Poly3BasisCellType> curl_Pk3(*m_ddrcore.cellBases(iT).Polyk3);
  MonomialCellIntegralsType int_monoT_2k = IntegrateCellMonomials(T, 2*m_ddrcore.degree());
  AT.topLeftCorner(dim_Pk3, dim_Pk3) += std::pow(T.diam(), 2) * GramMatrix(T, curl_Pk3, int_monoT_2k);

  // Lagrange multiplier terms
  if (dimCellRolyCompllpo(iT)>0){
    Eigen::MatrixXd Gram_Rclpo_Pk3 = GramMatrix(T, cellBasisRolyCompllpo(iT), *m_ddrcore.cellBases(iT).Polyk3, int_monoT_2k);
    AT.bottomLeftCorner(dimCellRolyCompllpo(iT), dim_Pk3) += Gram_Rclpo_Pk3;
    AT.topRightCorner(dim_Pk3, dimCellRolyCompllpo(iT)) -= Gram_Rclpo_Pk3.transpose();
  }
  
  return InverseProblem(AT);
}

//------------------------------------------------------------------------------
const Eigen::MatrixXd SerendipityProblem::SerendipityOperatorCell(const size_t iT, const Eigen::MatrixXd & LT) const
{
  return ( m_inverse_problem_cells[iT].solve(LT) ).topRows(m_ddrcore.cellBases(iT).Polyk3->dimension());
}


//------------------------------------------------------------------------------

Eigen::VectorXd SerendipityProblem::nDOFs_faces_SXGrad() const
{
  Eigen::VectorXd n_local_face_dofs = Eigen::VectorXd::Zero(mesh().n_faces());
  for (size_t iF = 0; iF < mesh().n_faces(); iF++){
    n_local_face_dofs(iF) = dimFacePolyl(iF);
  }
  
  return n_local_face_dofs;
}

Eigen::VectorXd SerendipityProblem::nDOFs_cells_SXGrad() const
{
  Eigen::VectorXd n_local_cell_dofs = Eigen::VectorXd::Zero(mesh().n_cells());
  for (size_t iT = 0; iT < mesh().n_cells(); iT++){
    n_local_cell_dofs(iT) += dimCellPolyl(iT);
  }
  
  return n_local_cell_dofs;
}

//------------------------------------------------------------------------------

Eigen::VectorXd SerendipityProblem::nDOFs_faces_SXCurl() const
{
  size_t dim_RkmoF = PolynomialSpaceDimension<Face>::Roly(m_ddrcore.degree()-1);
  Eigen::VectorXd n_local_face_dofs = Eigen::VectorXd::LinSpaced(mesh().n_faces(), dim_RkmoF, dim_RkmoF);
  for (size_t iF = 0; iF < mesh().n_faces(); iF++){
    n_local_face_dofs(iF) += dimFaceRolyCompllpo(iF);
  }
  
  return n_local_face_dofs;
}

Eigen::VectorXd SerendipityProblem::nDOFs_cells_SXCurl() const
{
  size_t dim_RkmoT = PolynomialSpaceDimension<Cell>::Roly(m_ddrcore.degree()-1);
  Eigen::VectorXd n_local_cell_dofs = Eigen::VectorXd::LinSpaced(mesh().n_cells(), dim_RkmoT, dim_RkmoT);
  for (size_t iT = 0; iT < mesh().n_cells(); iT++){
    n_local_cell_dofs(iT) += dimCellRolyCompllpo(iT);
  }
  
  return n_local_cell_dofs;
}

