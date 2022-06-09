#include <sxdiv.hpp>
#include <GMpoly_cell.hpp>

using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

SXDiv::SXDiv(const DDRCore & ddr_core, bool use_threads, std::ostream & output)
  : XDiv(ddr_core, use_threads, output)
{
  m_output << "[SXDiv] Initializing" << std::endl;
}


// L2 product with Serendipity curl

Eigen::MatrixXd SXDiv::computeL2ProductCurl(
                                        const size_t iT,
                                        const SXCurl & sx_curl,
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
      MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*degree());
      w_mass_Pk3_T = weight.value(T, T.center_mass()) * GramMatrix(T, *cellBases(iT).Polyk3, int_mono_2k);
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

  // list of full curls
  std::vector<Eigen::MatrixXd> curlOp(T.n_faces()+1);
  for (size_t iF = 0; iF < T.n_faces(); iF++){
    const Face & F = *T.face(iF);
    curlOp[iF] = sx_curl.extendOperator(T, F, sx_curl.faceCurl(F));
  }
  curlOp[T.n_faces()] = sx_curl.cellCurl(iT);
  
  // If we apply the curl on one side only we'll need the potentials
  if (side != "both"){
    // list of potentials
    std::vector<Eigen::MatrixXd> potentialOp(T.n_faces()+1);
    for (size_t iF = 0; iF < T.n_faces(); iF++){
      const Face & F = *T.face(iF);
      potentialOp[iF] = extendOperator(T, F, Eigen::MatrixXd::Identity(dimensionFace(F),dimensionFace(F)));
    }
    potentialOp[T.n_faces()] = m_cell_operators[iT]->potential;

  
    // Depending on side of curl
    if (side == "left"){
      return computeL2Product_with_Ops(iT, curlOp, potentialOp, penalty_factor, w_mass_Pk3_T, weight);
    }else{
      return computeL2Product_with_Ops(iT, potentialOp, curlOp, penalty_factor, w_mass_Pk3_T, weight);
    }
    
  }

  // Default: curl on both sides
  return computeL2Product_with_Ops(iT, curlOp, curlOp, penalty_factor, w_mass_Pk3_T, weight);

}

