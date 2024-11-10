#include <sxcurl.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_face.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

SXCurl::SXCurl(const DDRCore & ddr_core, const SerendipityProblem & ser_pro, bool use_threads, std::ostream & output)
  : VariableDOFSpace(
	     ddr_core.mesh(),
	     0,
	     PolynomialSpaceDimension<Edge>::Poly(ddr_core.degree()), 
	     ser_pro.nDOFs_faces_SXCurl(),
	     ser_pro.nDOFs_cells_SXCurl()    
	     ),
    m_ddr_core(ddr_core),
    m_ser_pro(ser_pro),
    m_xcurl(ddr_core, use_threads, output),
    m_face_transfer_operators(ddr_core.mesh().n_faces()),
    m_cell_transfer_operators(ddr_core.mesh().n_cells()),
    m_use_threads(use_threads),
    m_output(output)
{
  output << "[SXCurl] Initializing" << std::endl;
  if (m_use_threads) {
    m_output << "[SXCurl] Parallel execution" << std::endl;
  } else {
    m_output << "[SXCurl] Sequential execution" << std::endl;
  }

  // Construct face operators
  std::function<void(size_t, size_t)> construct_all_face_transfer_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iF = start; iF < end; iF++) {
          m_face_transfer_operators[iF].reset( new TransferOperators(_compute_face_transfer_operators(iF)) );
        } // for iF
      };

  m_output << "[SXCurl] Constructing face transfer operators" << std::endl;
  parallel_for(mesh().n_faces(), construct_all_face_transfer_operators, m_use_threads);

  // Construct cell operators
  std::function<void(size_t, size_t)> construct_all_cell_extensions_reductions
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_cell_transfer_operators[iT].reset( new TransferOperators(_compute_cell_transfer_operators(iT)) );
        } // for iT
      };

  m_output << "[SXCurl] Constructing cell transfer operators" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_extensions_reductions, m_use_threads);  
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd SXCurl::interpolate(const FunctionType & v, const int doe_cell, const int doe_face, const int doe_edge) const
{
  Eigen::VectorXd vh = Eigen::VectorXd::Zero(dimension());
  
  // We interpolate via xcurl and then apply the reduction.
  Eigen::VectorXd vh_xcurl = m_xcurl.interpolate(v, doe_cell, doe_face, doe_edge);
  
  // No change to edge DOFs
  vh.head(nDOFs_edges()) = vh_xcurl.head(nDOFs_edges());
  
  if (degree()>0){
    // Face DOFs
    for (size_t iF=0; iF < mesh().n_faces(); iF++){
      const Face & F = *mesh().face(iF);
      vh.segment(globalOffset(F), numLocalDofsFace(iF)) = RcurlFace(iF) * m_xcurl.restrict(F, vh_xcurl);
    }
    
    // Cell DOFs
    for (size_t iT=0; iT < mesh().n_cells(); iT++){
      const Cell & T = *mesh().cell(iT);
      vh.segment(globalOffset(T), numLocalDofsCell(iT)) = RcurlCell(iT) * m_xcurl.restrict(T, vh_xcurl);      
    }

  }
  
  return vh;
}

//------------------------------------------------------------------------------
// Transfer operators
//------------------------------------------------------------------------------

SXCurl::TransferOperators SXCurl::_compute_face_transfer_operators(size_t iF)
{
  const Face & F = *mesh().face(iF);
  
  MonomialFaceIntegralsType int_mono_2k_F = IntegrateFaceMonomials(F, 2*degree());
  
  //------------------------------------------------------------------------------
  // Serendipity operator
  //------------------------------------------------------------------------------
  // Compute LcurlF
  size_t dim_PkE = PolynomialSpaceDimension<Edge>::Poly(degree());
  size_t dim_RclpoF = m_ser_pro.dimFaceRolyCompllpo(iF);
  
  Eigen::MatrixXd LF = Eigen::MatrixXd::Zero(faceBases(iF).Polyk2->dimension() + dim_RclpoF, dimensionFace(iF));

  // Edges
  for (size_t iE=0; iE < F.n_edges(); iE++){
    const Edge & E = *F.edge(iE);
    QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2*degree());

    // Tangential component
    boost::multi_array<double, 2> Pk2F_tE_quad = 
          scalar_product( evaluate_quad<Function>::compute(*faceBases(iF).Polyk2, quad_2k_E),
                          E.tangent() );
    boost::multi_array<double, 2> PkE_quad = 
           evaluate_quad<Function>::compute(*edgeBases(E.global_index()).Polyk, quad_2k_E);
           
    LF.block(0, localOffset(F, E), faceBases(iF).Polyk2->dimension(), dim_PkE) +=  
          F.diam() * compute_gram_matrix(Pk2F_tE_quad, PkE_quad, quad_2k_E);

    // rot component
    if (degree()>0){
      auto rot_Pk2F = ScalarRotFamily<DDRCore::PolyBasisFaceType>(*faceBases(iF).Polyk2, F);
      LF.block(0, localOffset(F, E), faceBases(iF).Polyk2->dimension(), dim_PkE) 
          -= std::pow(F.diam(), 2) * F.edge_orientation(iE) *
                  compute_gram_matrix(
                      evaluate_quad<Function>::compute(rot_Pk2F, quad_2k_E),
                      evaluate_quad<Function>::compute(*edgeBases(E).Polyk, quad_2k_E),
                      quad_2k_E );
    }
  }
  
  // Face: rot rot
  if (degree()>1){
    QuadratureRule quad_2kmo_F = generate_quadrature_rule(F, 2*degree()-1);
    LF.block(0, localOffset(F), faceBases(iF).Polyk2->dimension(), faceBases(iF).Rolykmo->dimension()) 
         += std::pow(F.diam(), 2) * compute_gram_matrix(
                                         evaluate_quad<CurlCurl>::compute(*faceBases(iF).Polyk2, quad_2kmo_F),
                                         evaluate_quad<Function>::compute(*faceBases(iF).Rolykmo, quad_2kmo_F),
                                         quad_2kmo_F);
  } // degree()>1
    
  // Face: Lagrange multiplier
  if (dim_RclpoF>0){
    LF.block(faceBases(iF).Polyk2->dimension(), localOffset(F) + faceBases(iF).Rolykmo->dimension(), dim_RclpoF, dim_RclpoF)
      += GramMatrix(F, m_ser_pro.faceBasisRolyCompllpo(iF), int_mono_2k_F);
  }
          
  // Serendipity operator
  Eigen::MatrixXd SerF = m_ser_pro.SerendipityOperatorFace(F, LF);


  //------------------------------------------------------------------------------
  // Extension operator
  //------------------------------------------------------------------------------
  Eigen::MatrixXd ExtF = Eigen::MatrixXd::Zero(m_xcurl.dimensionFace(F), dimensionFace(F));
  
  // Identity for edge dofs, and Rk-1
  const size_t dim_RkmoF = PolynomialSpaceDimension<Face>::Roly(degree()-1);
  ExtF.topLeftCorner(localOffset(F) + dim_RkmoF, localOffset(F) + dim_RkmoF)
      = Eigen::MatrixXd::Identity(localOffset(F) + dim_RkmoF, localOffset(F) + dim_RkmoF);

  if (degree()>0){
    // Projection on Rck of Serendipity for the last part
    Eigen::MatrixXd Mass_RckF = GramMatrix(F, *faceBases(iF).RolyComplk, int_mono_2k_F);
    Eigen::MatrixXd Gram_RckF_Pk2F = GramMatrix(F, *faceBases(iF).RolyComplk, *faceBases(iF).Polyk2, int_mono_2k_F);
    ExtF.bottomRows(faceBases(iF).RolyComplk->dimension()) = Mass_RckF.ldlt().solve( Gram_RckF_Pk2F * SerF );
  }
  
  //------------------------------------------------------------------------------
  // Reduction operator
  //------------------------------------------------------------------------------
  Eigen::MatrixXd RedF = Eigen::MatrixXd::Zero(numLocalDofsFace(iF), m_xcurl.dimensionFace(F));
  
  // Identity Rk-1
  RedF.block(0, m_xcurl.localOffset(F), dim_RkmoF, dim_RkmoF) = Eigen::MatrixXd::Identity(dim_RkmoF, dim_RkmoF);

  if (dim_RclpoF>0){
    // Projection on Rclpo for the component on Rck
    Eigen::MatrixXd Mass_Rclpo = GramMatrix(F, m_ser_pro.faceBasisRolyCompllpo(iF), int_mono_2k_F);
    Eigen::MatrixXd Gram_Rclpo_Rck = GramMatrix(F, m_ser_pro.faceBasisRolyCompllpo(iF), *faceBases(iF).RolyComplk, int_mono_2k_F);
    RedF.bottomRightCorner(dim_RclpoF, faceBases(iF).RolyComplk->dimension())
          = Mass_Rclpo.ldlt().solve( Gram_Rclpo_Rck );
  }
  
  return TransferOperators(SerF, ExtF, RedF);
}


SXCurl::TransferOperators SXCurl::_compute_cell_transfer_operators(size_t iT)
{
  const Cell & T = *mesh().cell(iT);
  
  MonomialCellIntegralsType int_mono_2k_T = IntegrateCellMonomials(T, 2*degree());
  
  //------------------------------------------------------------------------------
  // Serendipity operator
  //------------------------------------------------------------------------------
  // Compute LcurlT
  size_t dim_RclpoT = m_ser_pro.dimCellRolyCompllpo(iT);
  
  Eigen::MatrixXd LT = Eigen::MatrixXd::Zero(cellBases(iT).Polyk3->dimension() + dim_RclpoT, dimensionCell(iT));

  // Faces contributions
  for (size_t iF=0; iF < T.n_faces(); iF++){
    const Face & F = *T.face(iF);
    const VectorRd nF = F.normal();
    QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2*degree());
    
    // Tangential component
    boost::multi_array<VectorRd, 2> Pk3T_tangentF_quad 
                      = transform_values_quad<VectorRd>(
                            evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2k_F),
                            [&nF](const VectorRd &z)->VectorRd { return z-(z.dot(nF))*nF;}
                          );
    boost::multi_array<VectorRd, 2> Pk2F_quad = 
           evaluate_quad<Function>::compute(*faceBases(F).Polyk2, quad_2k_F);
           
    Eigen::MatrixXd PF = extendOperator(T, F, facePotential(F));           
    LT.topRows(cellBases(iT).Polyk3->dimension()) += T.diam() * compute_gram_matrix(Pk3T_tangentF_quad, Pk2F_quad, quad_2k_F) * PF;
    
    // Curl term
    if (degree()>0){
      LT.topRows(cellBases(iT).Polyk3->dimension()) 
          += std::pow(T.diam(), 2) * T.face_orientation(iF) *
                compute_gram_matrix(
                      vector_product( evaluate_quad<Curl>::compute(*cellBases(iT).Polyk3, quad_2k_F), F.normal()),
                      evaluate_quad<Function>::compute(*faceBases(F).Polyk2, quad_2k_F),
                      quad_2k_F) * PF;
    }
  }
  
  // Element: curl curl component
  if (degree()>1){
    QuadratureRule quad_2kmt_T = generate_quadrature_rule(T, 2*degree()-2);
    LT.block(0, localOffset(T), cellBases(iT).Polyk3->dimension(), cellBases(iT).Rolykmo->dimension()) 
         += std::pow(T.diam(), 2) * compute_gram_matrix(
                                         evaluate_quad<CurlCurl>::compute(*cellBases(iT).Polyk3, quad_2kmt_T),
                                         evaluate_quad<Function>::compute(*cellBases(iT).Rolykmo, quad_2kmt_T),
                                         quad_2kmt_T);
  }

  // Element: Lagrange multiplier    
  if (dim_RclpoT>0){
    // Lagrange multiplier
    LT.block(cellBases(iT).Polyk3->dimension(), localOffset(T) + cellBases(iT).Rolykmo->dimension(), dim_RclpoT, dim_RclpoT)
      += GramMatrix(T, m_ser_pro.cellBasisRolyCompllpo(iT), int_mono_2k_T);
  }
        
  // Serendipity operator
  Eigen::MatrixXd SerT = m_ser_pro.SerendipityOperatorCell(T, LT);


  //------------------------------------------------------------------------------
  // Extension operator
  //------------------------------------------------------------------------------
  Eigen::MatrixXd ExtT = Eigen::MatrixXd::Zero(m_xcurl.dimensionCell(T), dimensionCell(T));
  
  // Identity for edge dofs
  const size_t n_dofs_edges = T.n_edges() * m_xcurl.numLocalDofsEdge();
  ExtT.topLeftCorner(n_dofs_edges, n_dofs_edges) = Eigen::MatrixXd::Identity(n_dofs_edges, n_dofs_edges);
  
  const size_t dim_RkmoT = PolynomialSpaceDimension<Cell>::Roly(degree()-1);
  if (degree()>0){
    // Faces grabbed from EcurlF
    for (size_t iF = 0; iF < T.n_faces(); iF++){
      const Face & F = *T.face(iF);
      ExtT.middleRows(m_xcurl.localOffset(T, F), m_xcurl.numLocalDofsFace()) 
          = extendOperator(T, F, EcurlFace(F).bottomRows(m_xcurl.numLocalDofsFace()));  
    }
    
    // Cell: identity on Rk-1, projection of serendipity on Rck
    ExtT.block(m_xcurl.localOffset(T), localOffset(T), dim_RkmoT, dim_RkmoT) = Eigen::MatrixXd::Identity(dim_RkmoT, dim_RkmoT);
    
    Eigen::MatrixXd Mass_RckT = GramMatrix(T, *cellBases(iT).RolyComplk, int_mono_2k_T);
    Eigen::MatrixXd Gram_RckT_Pk3T = GramMatrix(T, *cellBases(iT).RolyComplk, *cellBases(iT).Polyk3, int_mono_2k_T);
    ExtT.bottomRows(cellBases(iT).RolyComplk->dimension()) = Mass_RckT.ldlt().solve( Gram_RckT_Pk3T * SerT );
  }
  
  //------------------------------------------------------------------------------
  // Reduction operator
  //------------------------------------------------------------------------------
  Eigen::MatrixXd RedT = Eigen::MatrixXd::Zero(numLocalDofsCell(iT), m_xcurl.dimensionCell(T));
  
  if (degree()>0){
    // RRolyT on R^{k-1}
    // LHS matrix
    CurlBasis<DDRCore::GolyComplBasisCellType> CurlGck_basis(*cellBases(iT).GolyComplk);
    Eigen::MatrixXd MRT = GramMatrix(T, CurlGck_basis, *cellBases(iT).Rolykmo, int_mono_2k_T);
    // RHS matrix, starting with CT contribution
    Eigen::MatrixXd BRT = GramMatrix(T, *cellBases(iT).GolyComplk, *cellBases(iT).Polyk3) * m_xcurl.cellOperators(iT).curl;
    for (size_t iF = 0; iF < T.n_faces(); iF++){
      const Face & F = *T.face(iF);
      QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2*degree());
      boost::multi_array<VectorRd, 2> Gck_cross_nF_quad = 
          vector_product( evaluate_quad<Function>::compute(*cellBases(iT).GolyComplk, quad_2k_F), F.normal() );
      boost::multi_array<VectorRd, 2> Pk2F_quad = evaluate_quad<Function>::compute(*faceBases(F).Polyk2, quad_2k_F);
      // Need to complete RcurlF with edge unknowns
      Eigen::MatrixXd CompleteRcurlF = Eigen::MatrixXd::Zero(dimensionFace(F), m_xcurl.dimensionFace(F));
      CompleteRcurlF.topLeftCorner(localOffset(F), localOffset(F)) = Eigen::MatrixXd::Identity(localOffset(F), localOffset(F));
      CompleteRcurlF.bottomRows(numLocalDofsFace(F)) = RcurlFace(F);
      Eigen::MatrixXd PERF = m_xcurl.extendOperator(T, F, facePotential(F) * CompleteRcurlF);
      BRT -= T.face_orientation(iF) * compute_gram_matrix(Gck_cross_nF_quad, Pk2F_quad, quad_2k_F) * PERF;    
    }
    RedT.topRows(dim_RkmoT) = MRT.fullPivLu().solve(BRT);

    // Projection on R^{c,l+1} of R^{c,k} component
    if (m_ser_pro.dimCellRolyCompllpo(iT)>0){
      Eigen::MatrixXd Mass_Rclpo = GramMatrix(T, m_ser_pro.cellBasisRolyCompllpo(iT), int_mono_2k_T);
      Eigen::MatrixXd Gram_Rclpo_Rck = GramMatrix(T, m_ser_pro.cellBasisRolyCompllpo(iT), *cellBases(iT).RolyComplk, int_mono_2k_T);
      RedT.bottomRightCorner(m_ser_pro.dimCellRolyCompllpo(iT), cellBases(iT).RolyComplk->dimension())
            = Mass_Rclpo.ldlt().solve( Gram_Rclpo_Rck );
    }
  }
  
  return TransferOperators(SerT, ExtT, RedT);
}

//------------- L2 product with gradient ----------------//

Eigen::MatrixXd SXCurl::computeL2ProductGradient(
                                     const size_t iT,
                                     const SXGrad & sx_grad,
                                     const std::string & side,
                                     const double & penalty_factor,
                                     const Eigen::MatrixXd & mass_Pk3_T,
                                     const IntegralWeight & weight
                                     ) const
{
  const Cell & T = *mesh().cell(iT);

  // create the weighted mass matrix, with simple product if weight is constant
  Eigen::MatrixXd w_mass_Pk3_T;
  if (weight.deg(T)==0){
    // constant weight
    if (mass_Pk3_T.rows()==1){
      // We have to compute the mass matrix
      MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2);
      w_mass_Pk3_T = weight.value(T, T.center_mass()) * GramMatrix(T, *cellBases(iT).Polyk3, int_mono_2kp2);
    }else{
      w_mass_Pk3_T = weight.value(T, T.center_mass()) * mass_Pk3_T;
    }
  }else{
    // weight is not constant, we create a weighted mass matrix
    QuadratureRule quad_2kpw_T = generate_quadrature_rule(T, 2 * degree() + weight.deg(T));
    auto basis_Pk3_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk3, quad_2kpw_T);
    std::function<double(const Eigen::Vector3d &)> weight_T 
              = [&T, &weight](const Eigen::Vector3d &x)->double {
                  return weight.value(T, x);
                };
    w_mass_Pk3_T = compute_weighted_gram_matrix(weight_T, basis_Pk3_T_quad, basis_Pk3_T_quad, quad_2kpw_T, "sym");
  }


  // list of full gradients
  std::vector<Eigen::MatrixXd> gradientOp(T.n_edges()+T.n_faces()+1);
  for (size_t iE = 0; iE < T.n_edges(); iE++){
    const Edge & E = *T.edge(iE);
    gradientOp[iE] = sx_grad.extendOperator(T, E, sx_grad.edgeGradient(E));
  }
  for (size_t iF = 0; iF < T.n_faces(); iF++){
    const Face & F = *T.face(iF);
    gradientOp[T.n_edges()+iF] = sx_grad.extendOperator(T, F, sx_grad.faceGradient(F));
  }
  gradientOp[T.n_edges()+T.n_faces()] = sx_grad.cellGradient(iT);
  
  // If we apply the gradient on one side only we'll need the potentials
  if (side != "both"){
    // list of potentials
    std::vector<Eigen::MatrixXd> potentialOp(T.n_edges()+T.n_faces()+1);
    for (size_t iE = 0; iE < T.n_edges(); iE++){
      const Edge & E = *T.edge(iE);
      potentialOp[iE] = extendOperator(T, E, Eigen::MatrixXd::Identity(dimensionEdge(E),dimensionEdge(E)));
    }
    for (size_t iF = 0; iF < T.n_faces(); iF++){
      const Face & F = *T.face(iF);
      potentialOp[T.n_edges()+iF] = extendOperator(T, F, facePotential(F));
    }
    potentialOp[T.n_edges()+T.n_faces()] = cellPotential(iT);
  
    // Depending on side of gradient
    if (side == "left"){
      return m_xcurl.computeL2Product_with_Ops(iT, gradientOp, potentialOp, penalty_factor, w_mass_Pk3_T, weight);
    }else{
      return m_xcurl.computeL2Product_with_Ops(iT, potentialOp, gradientOp, penalty_factor, w_mass_Pk3_T, weight);
    }
    
  }

  // Default: gradient on both sides
  return m_xcurl.computeL2Product_with_Ops(iT, gradientOp, gradientOp, penalty_factor, w_mass_Pk3_T, weight);

}


//------------------------------------------------------------------------------
// Vertex values 
//------------------------------------------------------------------------------
std::vector<VectorRd> SXCurl::computeVertexValues(const Eigen::VectorXd & u) const
{
  std::vector<VectorRd> values(mesh().n_vertices(), VectorRd::Zero());
  
  // Value at each vertex obtained averaging the values from all the cells around
  for (Vertex * V : mesh().get_vertices()){
    size_t iV = V->global_index();
 
    for (Cell * T : V->get_cells()){
      size_t iT = T->global_index();
      Eigen::VectorXd pTuT = cellPotential(iT) * restrict(*T, u);
      
      for (size_t i=0; i < cellBases(iT).Polyk3->dimension(); i++){
        values[iV] += pTuT(i) * cellBases(iT).Polyk3->function(i, V->coords());
      }
    }
    
    values[iV] /= V->n_cells();

  }
  
  return values;
}                  



