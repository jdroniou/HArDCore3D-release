#include <sxgrad.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_face.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

SXGrad::SXGrad(const DDRCore & ddr_core, const SerendipityProblem & ser_pro, bool use_threads, std::ostream & output)
  : VariableDOFSpace(
	     ddr_core.mesh(),
	     1,
	     PolynomialSpaceDimension<Edge>::Poly(ddr_core.degree()-1), 
	     ser_pro.nDOFs_faces_SXGrad(),
	     ser_pro.nDOFs_cells_SXGrad()    
	     ),
    m_ddr_core(ddr_core),
    m_ser_pro(ser_pro),
    m_xgrad(ddr_core, use_threads, output),
    m_face_transfer_operators(ddr_core.mesh().n_faces()),
    m_cell_transfer_operators(ddr_core.mesh().n_cells()),
    m_use_threads(use_threads),
    m_output(output)
{
  output << "[SXGrad] Initializing" << std::endl;
  if (m_use_threads) {
    m_output << "[SXGrad] Parallel execution" << std::endl;
  } else {
    m_output << "[SXGrad] Sequential execution" << std::endl;
  }

  // Construct face operators
  std::function<void(size_t, size_t)> construct_all_face_transfer_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iF = start; iF < end; iF++) {
          m_face_transfer_operators[iF].reset( new TransferOperators(_compute_face_transfer_operators(iF)) );
        } // for iF
      };

  m_output << "[SXGrad] Constructing face transfer operators" << std::endl;
  parallel_for(mesh().n_faces(), construct_all_face_transfer_operators, m_use_threads);

  // Construct cell operators
  std::function<void(size_t, size_t)> construct_all_cell_extensions_reductions
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_cell_transfer_operators[iT].reset( new TransferOperators(_compute_cell_transfer_operators(iT)) );
        } // for iT
      };

  m_output << "[SXGrad] Constructing cell transfer operators" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_extensions_reductions, m_use_threads);  
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd SXGrad::interpolate(const FunctionType & q, const int doe_cell, const int doe_face, const int doe_edge) const
{
  Eigen::VectorXd qh = Eigen::VectorXd::Zero(dimension());
  
  // We interpolate via xgrad and then apply the reduction.
  Eigen::VectorXd qh_xgrad = m_xgrad.interpolate(q, doe_cell, doe_face, doe_edge);
  
  // No change to vertex and edge DOFs
  qh.head(nDOFs_vertices() + nDOFs_edges()) = qh_xgrad.head(nDOFs_vertices() + nDOFs_edges());
  
  if (degree()>0){
    // Face DOFs
    for (size_t iF=0; iF < mesh().n_faces(); iF++){
      if (numLocalDofsFace(iF)>0) {
        const Face & F = *mesh().face(iF);
        qh.segment(globalOffset(F), numLocalDofsFace(iF)) = RgradFace(iF) * m_xgrad.restrict(F, qh_xgrad);
      }
    }
      
    // Cell DOFs
    for (size_t iT=0; iT < mesh().n_cells(); iT++){
      if (m_ser_pro.dimCellPolyl(iT)>0) {
        const Cell & T = *mesh().cell(iT);
        qh.segment(globalOffset(T), numLocalDofsCell(iT)) = RgradCell(iT) * m_xgrad.restrict(T, qh_xgrad);      
      }
    }
  }
  
  return qh;
}

//------------------------------------------------------------------------------
// Transfer operators
//------------------------------------------------------------------------------

SXGrad::TransferOperators SXGrad::_compute_face_transfer_operators(size_t iF)
{
  const Face & F = *mesh().face(iF);
 
  MonomialFaceIntegralsType int_mono_2k_F = IntegrateFaceMonomials(F, 2*degree());
  
  //------------------------------------------------------------------------------
  // Serendipity operator
  //------------------------------------------------------------------------------
  // Compute LcurlF
  size_t dim_RclpoF = m_ser_pro.dimFaceRolyCompllpo(iF);
  size_t dim_PlF = m_ser_pro.dimFacePolyl(iF);
  
  Eigen::MatrixXd LF = Eigen::MatrixXd::Zero(faceBases(iF).Polyk2->dimension() + dim_RclpoF, dimensionFace(iF));

  // Edges
  for (size_t iE=0; iE < F.n_edges(); iE++){
    const Edge & E = *F.edge(iE);
    QuadratureRule quad_2kpo_E = generate_quadrature_rule(E, 2*degree()+1);

    // contribution q_E'
    boost::multi_array<double, 2> Pk2F_tE_quad = 
          scalar_product( evaluate_quad<Function>::compute(*faceBases(iF).Polyk2, quad_2kpo_E),
                          E.tangent() );
    boost::multi_array<double, 2> PkE_quad = 
           evaluate_quad<Function>::compute(*edgeBases(E.global_index()).Polyk, quad_2kpo_E);
           
    LF.topRows(faceBases(iF).Polyk2->dimension()) +=  
          F.diam() * compute_gram_matrix(Pk2F_tE_quad, PkE_quad, quad_2kpo_E) * extendOperator(F, E, m_xgrad.edgeOperators(E).gradient);

    // contribution q_E with mu.nFE
    if (dim_RclpoF > 0){
      boost::multi_array<double, 2> Rclpo_nFE_quad = 
          scalar_product( evaluate_quad<Function>::compute(m_ser_pro.faceBasisRolyCompllpo(iF), quad_2kpo_E),
                          F.edge_normal(iE) );
      LF.bottomRows(dim_RclpoF) += F.edge_orientation(iE) * compute_gram_matrix(
                                      Rclpo_nFE_quad,
                                      evaluate_quad<Function>::compute(*edgeBases(E).Polykpo, quad_2kpo_E),
                                      quad_2kpo_E 
                                      ) 
                                     * extendOperator(F, E, m_xgrad.edgeOperators(E).potential);
    }
  }

  // Face: q_F
  if (dim_RclpoF > 0){
    LF.block(faceBases(iF).Polyk2->dimension(), localOffset(F), dim_RclpoF, dim_PlF) 
         -= GramMatrix(F, 
                      DivergenceBasis<SerendipityProblem::RolyCompllpoBasisFaceType>(m_ser_pro.faceBasisRolyCompllpo(iF)), 
                      m_ser_pro.faceBasisPolyl(iF), 
                      int_mono_2k_F );
  }
    
  // Serendipity operator
  Eigen::MatrixXd SerF = m_ser_pro.SerendipityOperatorFace(F, LF);


  //------------------------------------------------------------------------------
  // Extension operator
  //------------------------------------------------------------------------------
  Eigen::MatrixXd ExtF = Eigen::MatrixXd::Zero(m_xgrad.dimensionFace(F), dimensionFace(F));

  // Identity on vertices and edges
  size_t dim_PkmoF = PolynomialSpaceDimension<Face>::Poly(degree()-1);
  ExtF.topLeftCorner(m_xgrad.dimensionFace(F) - dim_PkmoF, dimensionFace(F) - dim_PlF) 
      = Eigen::MatrixXd::Identity(m_xgrad.dimensionFace(F) - dim_PkmoF, dimensionFace(F) - dim_PlF);
  
  if (degree()>0){
    // Face extension E_{PF}^{k-1}
    // LHS matrix
    Eigen::MatrixXd MFEpoly = GramMatrix(F, 
                                      DivergenceBasis<DDRCore::RolyComplBasisFaceType>(*faceBases(iF).RolyComplk),
                                      *faceBases(iF).Polykmo,
                                      int_mono_2k_F);
    // RHS matrix, starting with serendipity contribution
    Eigen::MatrixXd BFEpoly = - GramMatrix(F, *faceBases(iF).RolyComplk, *faceBases(iF).Polyk2, int_mono_2k_F) * SerF;
    
    // Edge contributions
    for (size_t iE=0; iE < F.n_edges(); iE++){
      const Edge & E = *F.edge(iE);
      QuadratureRule quad_2kpo_E = generate_quadrature_rule(E, 2*degree()+1);

      boost::multi_array<double, 2> Rck_nFE_quad = 
          scalar_product( evaluate_quad<Function>::compute(*faceBases(iF).RolyComplk, quad_2kpo_E),
                          F.edge_normal(iE) );
      BFEpoly += F.edge_orientation(iE) * compute_gram_matrix(
                                      Rck_nFE_quad,
                                      evaluate_quad<Function>::compute(*edgeBases(E).Polykpo, quad_2kpo_E),
                                      quad_2kpo_E 
                                      ) 
                                     * extendOperator(F, E, m_xgrad.edgeOperators(E).potential);

    }
    
    ExtF.bottomRows(dim_PkmoF) = MFEpoly.fullPivLu().solve(BFEpoly);
  }
  
  //------------------------------------------------------------------------------
  // Reduction operator
  //------------------------------------------------------------------------------
  Eigen::MatrixXd RedF = Eigen::MatrixXd::Zero(dim_PlF, m_xgrad.dimensionFace(F));
  
  // Projection of face component on P^l(F)
  if (dim_PlF>0){
    Eigen::MatrixXd Mass_Pl = GramMatrix(F, m_ser_pro.faceBasisPolyl(iF), int_mono_2k_F);
    Eigen::MatrixXd Gram_Pl_Pkmo = GramMatrix(F, m_ser_pro.faceBasisPolyl(iF), *faceBases(iF).Polykmo, int_mono_2k_F);
    RedF.rightCols(faceBases(iF).Polykmo->dimension()) = Mass_Pl.ldlt().solve( Gram_Pl_Pkmo );
  }
  
  return TransferOperators(SerF, ExtF, RedF);
}


SXGrad::TransferOperators SXGrad::_compute_cell_transfer_operators(size_t iT)
{
  const Cell & T = *mesh().cell(iT);
  
  MonomialCellIntegralsType int_mono_2k_T = IntegrateCellMonomials(T, 2*degree());
  
  //------------------------------------------------------------------------------
  // Serendipity operator
  //------------------------------------------------------------------------------
  // Compute LgradT
  size_t dim_RclpoT = m_ser_pro.dimCellRolyCompllpo(iT);
  
  Eigen::MatrixXd LT = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension() + dim_RclpoT, dimensionCell(iT));

  // Faces contributions
  for (size_t iF=0; iF < T.n_faces(); iF++){
    const Face & F = *T.face(iF);
    const VectorRd nF = F.normal();
    QuadratureRule quad_2kpo_F = generate_quadrature_rule(F, 2*degree()+1);
    
    // Face gradient and tangential component
    boost::multi_array<VectorRd, 2> Pk3T_tangentF_quad 
                      = transform_values_quad<VectorRd>(
                            evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2kpo_F),
                            [&nF](const VectorRd &z)->VectorRd { return z-(z.dot(nF))*nF;}
                          );
    boost::multi_array<VectorRd, 2> Pk2F_quad = 
           evaluate_quad<Function>::compute(*faceBases(F).Polyk2, quad_2kpo_F);
           
    Eigen::MatrixXd GF = extendOperator(T, F, faceGradient(F));           
    LT.topRows(cellBases(iT).Polyk3->dimension()) += T.diam() * compute_gram_matrix(Pk3T_tangentF_quad, Pk2F_quad, quad_2kpo_F) * GF;

    // Lagrange multiplier
    if (dim_RclpoT>0){
     Eigen::MatrixXd PF = extendOperator(T, F, facePotential(F));           
      boost::multi_array<double, 2> RclpoT_dot_nF_quad
                      = scalar_product(
                            evaluate_quad<Function>::compute(m_ser_pro.cellBasisRolyCompllpo(iT), quad_2kpo_F),
                            nF );
      LT.bottomRows(dim_RclpoT) += T.face_orientation(iF) 
                                      * compute_gram_matrix( RclpoT_dot_nF_quad,
                                                             evaluate_quad<Function>::compute(*faceBases(F).Polykpo, quad_2kpo_F),
                                                             quad_2kpo_F )
                                        * PF;
    }  
  }

  // Element: Lagrange multiplier
  if (dim_RclpoT>0){
    // Lagrange multiplier
    LT.block(cellBases(iT).Polyk3->dimension(), localOffset(T), dim_RclpoT, numLocalDofsCell(iT))
      -= GramMatrix(T, 
          DivergenceBasis<SerendipityProblem::RolyCompllpoBasisCellType>(m_ser_pro.cellBasisRolyCompllpo(iT)), 
          m_ser_pro.cellBasisPolyl(iT),
          int_mono_2k_T);
  }
        
  // Serendipity operator
  Eigen::MatrixXd SerT = m_ser_pro.SerendipityOperatorCell(T, LT);


  //------------------------------------------------------------------------------
  // Extension operator
  //------------------------------------------------------------------------------
  Eigen::MatrixXd ExtT = Eigen::MatrixXd::Zero(m_xgrad.dimensionCell(T), dimensionCell(T));
  
  // Identity for vertex and edge dofs
  const size_t n_dofs_vertex_edges = T.n_vertices() * m_xgrad.numLocalDofsVertex() + T.n_edges() * m_xgrad.numLocalDofsEdge();
  ExtT.topLeftCorner(n_dofs_vertex_edges, n_dofs_vertex_edges) = Eigen::MatrixXd::Identity(n_dofs_vertex_edges, n_dofs_vertex_edges);
  
  if (degree()>0){
    // Faces grabbed from EcurlF
    for (size_t iF = 0; iF < T.n_faces(); iF++){
      const Face & F = *T.face(iF);
      ExtT.middleRows(m_xgrad.localOffset(T, F), m_xgrad.numLocalDofsFace()) 
          = extendOperator(T, F, EgradFace(F).bottomRows(m_xgrad.numLocalDofsFace()));  
    }
    
    // Cell
    // LHS matrix
    Eigen::MatrixXd MET = GramMatrix(T,
                                     DivergenceBasis<DDRCore::RolyComplBasisCellType>(*cellBases(iT).RolyComplk),
                                     *cellBases(iT).Polykmo,
                                     int_mono_2k_T );    
    // RHS matrix
    Eigen::MatrixXd BET = - GramMatrix(T, *cellBases(iT).RolyComplk, *cellBases(iT).Polyk3, int_mono_2k_T) * SerT;
    for (size_t iF=0; iF < T.n_faces(); iF++){
      const Face & F = *T.face(iF);
      QuadratureRule quad_2kpo_F = generate_quadrature_rule(F, 2*degree()+1);
      
      boost::multi_array<double, 2> RckT_dot_nF_quad =
            scalar_product( evaluate_quad<Function>::compute(*cellBases(iT).RolyComplk, quad_2kpo_F), F.normal() );
      boost::multi_array<double, 2> PkpoF_quad = evaluate_quad<Function>::compute(*faceBases(F).Polykpo, quad_2kpo_F);

      BET += T.face_orientation(iF) * compute_gram_matrix(RckT_dot_nF_quad, PkpoF_quad, quad_2kpo_F) 
                * extendOperator(T, F, facePotential(F));

    }    
    ExtT.bottomRows(m_xgrad.numLocalDofsCell()) = MET.fullPivLu().solve( BET );
    
  }
  
  //------------------------------------------------------------------------------
  // Reduction operator
  //------------------------------------------------------------------------------
  size_t dim_PlT = m_ser_pro.dimCellPolyl(T);
  Eigen::MatrixXd RedT = Eigen::MatrixXd::Zero(dim_PlT, m_xgrad.dimensionCell(T));
  
  if (dim_PlT>0){
    // RPolyT on P^{l}
    // LHS matrix
    Eigen::MatrixXd MRT = GramMatrix( T, 
                            DivergenceBasis<SerendipityProblem::RolyCompllpoBasisCellType>(m_ser_pro.cellBasisRolyCompllpo(iT)), 
                            m_ser_pro.cellBasisPolyl(iT),
                            int_mono_2k_T );
    // RHS matrix, starting with GT contribution
    Eigen::MatrixXd BRT = - GramMatrix(T, m_ser_pro.cellBasisRolyCompllpo(iT), *cellBases(iT).Polyk3) * m_xgrad.cellOperators(iT).gradient;
    for (size_t iF = 0; iF < T.n_faces(); iF++){
      const Face & F = *T.face(iF);
      QuadratureRule quad_2kpo_F = generate_quadrature_rule(F, 2*degree()+1);
      
      boost::multi_array<double, 2> Rclpo_dot_nF_quad = 
          scalar_product( evaluate_quad<Function>::compute(m_ser_pro.cellBasisRolyCompllpo(iT), quad_2kpo_F), F.normal() );
      boost::multi_array<double, 2> PkpoF_quad = evaluate_quad<Function>::compute(*faceBases(F).Polykpo, quad_2kpo_F);
      
      // Need to complete RgradF with vertex and edge unknowns
      Eigen::MatrixXd CompleteRgradF = Eigen::MatrixXd::Zero(dimensionFace(F), m_xgrad.dimensionFace(F));
      CompleteRgradF.topLeftCorner(localOffset(F), localOffset(F)) = Eigen::MatrixXd::Identity(localOffset(F), localOffset(F));
      CompleteRgradF.bottomRows(numLocalDofsFace(F)) = RgradFace(F);
      Eigen::MatrixXd PERF = m_xgrad.extendOperator(T, F, facePotential(F) * CompleteRgradF);

      BRT += T.face_orientation(iF) * compute_gram_matrix(Rclpo_dot_nF_quad, PkpoF_quad, quad_2kpo_F) * PERF;    
    }
    
    RedT = MRT.fullPivLu().solve(BRT);
  }
  
  return TransferOperators(SerT, ExtT, RedT);
}

//------------------------------------------------------------------------------
// Vertex values 
//------------------------------------------------------------------------------
std::vector<double> SXGrad::computeVertexValues(const Eigen::VectorXd & u) const
{
  std::vector<double> values(mesh().n_vertices(), 0.);
  
  // Value at each vertex obtained averaging the values from all the cells around
  for (Vertex * V : mesh().get_vertices()){
    size_t iV = V->global_index();
 
    for (Cell * T : V->get_cells()){
      size_t iT = T->global_index();
      Eigen::VectorXd pTuT = cellPotential(iT) * restrict(*T, u);
      
      for (size_t i=0; i < cellBases(iT).Polykpo->dimension(); i++){
        values[iV] += pTuT(i) * cellBases(iT).Polykpo->function(i, V->coords());
      }
    }
    
    values[iV] /= V->n_cells();

  }
  
  return values;
}                  

