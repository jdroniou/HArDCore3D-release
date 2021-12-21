#include <GMpoly_cell.hpp>

using namespace HArDCore3D;


// IntegrateCellMonomials_onEdges
std::vector<MonomialCellIntegralsType> HArDCore3D::IntegrateCellMonomials_onEdges(const Cell & T, const size_t maxdeg) {
  // Number of monomials
  const size_t nb_poly = PolynomialSpaceDimension<Cell>::Poly(maxdeg);

  std::vector<MonomialCellIntegralsType> integrals;
  integrals.resize(T.n_edges());
  
  // Dimension assumed
  const size_t d = 3;
  
  // Coordinate transformation data
  VectorRd xT = T.center_mass();
  double hT = T.diam();
  
  // Create powers of all monomials
  std::vector<VectorZd> powers = MonomialPowers<Cell>::complete(maxdeg);
  
  // Loop over the edges
  for (size_t iE = 0; iE < T.n_edges(); iE++){
    VectorRd v1 = T.edge(iE)->vertex(0)->coords();
    VectorRd v2 = T.edge(iE)->vertex(1)->coords();
    VectorRd xE = T.edge(iE)->center_mass();
    double d1 = (v1-xE).norm();
    double d2 = (v2-xE).norm();
    
    VectorRd v1_trans = (v1 - xT)/hT;
    VectorRd v2_trans = (v2 - xT)/hT;
    
    // Loop over all the monomials
    for (size_t m = 0; m < nb_poly; m++){
      // Start with vertex values
      double int_e_mono = d1 * std::pow(v1_trans.x(),powers[m](0)) * std::pow(v1_trans.y(),powers[m](1)) * std::pow(v1_trans.z(),powers[m](2));
      int_e_mono += d2 * std::pow(v2_trans.x(),powers[m](0)) * std::pow(v2_trans.y(),powers[m](1)) * std::pow(v2_trans.z(),powers[m](2));

      // Gradient correction (each partial derivative)
      for (size_t ip = 0; ip < d; ip++)
      {
        if (powers[m](ip) > 0){
          VectorZd powers_diff = VectorZd::Zero(d);
          powers_diff(ip) = -1;
          int_e_mono += (xE(ip)-xT(ip))/hT * powers[m](ip) * integrals[iE][powers[m]+powers_diff];
        }
      }

      // Store
      integrals[iE][powers[m]] = int_e_mono / (d-2+powers[m].sum());

    } // for m
  } // for iE

  return integrals;
}

// IntegrateCellMonomials_onFaces
std::vector<MonomialCellIntegralsType> HArDCore3D::IntegrateCellMonomials_onFaces(const Cell & T, const size_t maxdeg, std::vector<MonomialCellIntegralsType> & integrals_edges) {

  // Number of monomials
  const size_t nb_poly = PolynomialSpaceDimension<Cell>::Poly(maxdeg);

  std::vector<MonomialCellIntegralsType> integrals;
  integrals.resize(T.n_faces());
  
  // Dimension assumed
  const size_t d = 3;
  
  // Coordinate transformation data
  VectorRd xT = T.center_mass();
  double hT = T.diam();
  
  // Create powers of all monomials
  std::vector<VectorZd> powers = MonomialPowers<Cell>::complete(maxdeg);

  // Loop over faces
  for (size_t iF=0; iF<T.n_faces(); iF++){
    Face *F = T.face(iF);
    // Outer unit normal and arbitrary point on face
    VectorRd xF = F->center_mass();

    // Loop over all the monomials
    for (size_t m = 0; m < nb_poly; m++){
      double int_f_mono = 0.;
      // Loop over this face's edges and add their contribution
      for (size_t iE = 0; iE < F->n_edges(); iE++){
        VectorRd nTFE = F->edge_normal(iE); // not necessarily outer normal, so abs() required below - we assume that F is star-shaped w.r.t. xF
        VectorRd xFE = F->edge(iE)->center_mass();
        
        // Find where this edge is in the list of edges on the cell (and thus in integrals_edges)
        size_t iE_inT = T.index_edge(F->edge(iE));

        int_f_mono += abs(nTFE.dot(xFE-xF))*integrals_edges[iE_inT][powers[m]];
      } // for iE
      // Gradient correction: linear combination (partial derivatives) of three face integrals of one lower degree
      for (size_t ip = 0; ip < d; ip++){
        if (powers[m](ip) > 0){
          VectorZd powers_diff = VectorZd::Zero(d);
          powers_diff(ip) = -1;
          int_f_mono += (xF(ip)-xT(ip))/hT * powers[m](ip) * integrals[iF][powers[m]+powers_diff];
        }
      }

      // Store
      integrals[iF][powers[m]] = int_f_mono/(d-1+powers[m].sum());
    } // for m
    
  } // for iF
  
  return integrals;
}

// IntegrateCellMonomials
MonomialCellIntegralsType HArDCore3D::IntegrateCellMonomials(const Cell & T, const size_t maxdeg) {

  // Number of monomials
  const size_t nb_poly = PolynomialSpaceDimension<Cell>::Poly(maxdeg);

  MonomialCellIntegralsType integrals;
  
  // Dimension assumed
  const size_t d = 3;
  
  // Coordinate transformation data
  VectorRd xT = T.center_mass();
  
  // Create powers of all monomials
  std::vector<VectorZd> powers = MonomialPowers<Cell>::complete(maxdeg);

  // All integrals over the faces
  std::vector<MonomialCellIntegralsType> integrals_edges = IntegrateCellMonomials_onEdges(T, maxdeg);
  std::vector<MonomialCellIntegralsType> integrals_faces = IntegrateCellMonomials_onFaces(T, maxdeg, integrals_edges);

  // Loop over all the monomials
  for (size_t m = 0; m < nb_poly; m++){
    double int_mono = 0.;
    // Adding face contributions
    for (size_t iF = 0; iF < T.n_faces(); iF++){
      // Outer unit normal and arbitrary point on face
      VectorRd nTF = T.face_normal(iF);
      VectorRd xF = T.face(iF)->center_mass();
    
      int_mono += nTF.dot(xF-xT)*integrals_faces[iF][powers[m]];
    } // for iF
    
    // Store
    integrals[powers[m]] = int_mono / (d+powers[m].sum());
  } // for m
  
  return integrals;
}

// Check integrals list
MonomialCellIntegralsType HArDCore3D::CheckIntegralsDegree(const Cell & T, const size_t degree, const MonomialCellIntegralsType & mono_int_map){
  if (mono_int_map.size() >= PolynomialSpaceDimension<Cell>::Poly(degree)){
    return mono_int_map;
  }else{
    return IntegrateCellMonomials(T, degree);
  }
}


// GramMatrix
Eigen::MatrixXd HArDCore3D::GramMatrix(const Cell & T, const MonomialScalarBasisCell & basis1, const MonomialScalarBasisCell & basis2, MonomialCellIntegralsType mono_int_map)
{
  // Dimension of the gram matrix
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  size_t totaldegree = basis1.max_degree()+basis2.max_degree();
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  // Obtain integration data from IntegrateCellMonomials
  MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
  
  for (size_t i = 0; i < dim1; i++) {
    for (size_t j = 0; j < dim2; j++) {
      gm(i,j) = intmap.at(basis1.powers(i) + basis2.powers(j));
    } // for j
  }   // for i
  
  return gm;
}

// GramMatrix for RolyCompl
Eigen::MatrixXd HArDCore3D::GramMatrix(const Cell & T, const RolyComplBasisCell & basis1, const RolyComplBasisCell & basis2, MonomialCellIntegralsType mono_int_map)
{
  // Dimension of the gram matrix
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  size_t totaldegree = basis1.max_degree()+basis2.max_degree();
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  // Obtain integration data from IntegrateCellMonomials
  MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
  
  for (size_t k = 0; k < 3; k++) {
    VectorZd powers = VectorZd::Zero(3);
    powers(k) = 2;
    for (size_t i = 0; i < dim1; i++) {
      for (size_t j = 0; j < dim2; j++) {
          gm(i,j) += intmap.at(powers + basis1.powers(i) + basis2.powers(j));
      } // for j
    }   // for i
  }     // for k
  
  return gm;
}

// GramMatrix for GolyCompl
Eigen::MatrixXd HArDCore3D::GramMatrix(const Cell & T, const GolyComplBasisCell & basis1, const GolyComplBasisCell & basis2, MonomialCellIntegralsType mono_int_map)
{
  // Dimension of the gram matrix
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  size_t totaldegree = basis1.max_degree()+basis2.max_degree();
  Eigen::MatrixXd gm(dim1,dim2);
  
  // Obtain integration data from IntegrateCellMonomials
  MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
  
  size_t dim_l1 = basis1.dimPkmo();
  size_t dim_l2 = basis2.dimPkmo();
  
  gm.block(0, 0, dim_l1, dim_l2) = GMGolyCompl(T, basis1, basis2, 0, 0, intmap, 1, 1) + GMGolyCompl(T, basis1, basis2, 0, 0, intmap, 2, 2);
  gm.block(0, dim_l2, dim_l1, dim_l2) = -GMGolyCompl(T, basis1, basis2, 0, 1, intmap, 0, 1);
  gm.block(0, 2*dim_l2, dim_l1, dim2-2*dim_l2) = -GMGolyCompl(T, basis1, basis2, 0, 2, intmap, 0, 2);
  gm.block(dim_l1, 0, dim_l1, dim_l2) = -GMGolyCompl(T, basis1, basis2, 1, 0, intmap, 1, 0);
  gm.block(dim_l1, dim_l2, dim_l1, dim_l2) = GMGolyCompl(T, basis1, basis2, 1, 1, intmap, 0, 0) + GMGolyCompl(T, basis1, basis2, 1, 1, intmap, 2, 2);
  gm.block(dim_l1, 2*dim_l2, dim_l1, dim2-2*dim_l2) = -GMGolyCompl(T, basis1, basis2, 1, 2, intmap, 1, 2);
  gm.block(2*dim_l1, 0, dim1-2*dim_l1, dim_l2) = -GMGolyCompl(T, basis1, basis2, 2, 0, intmap, 2, 0);
  gm.block(2*dim_l1, dim_l2, dim1-2*dim_l1, dim_l2) = -GMGolyCompl(T, basis1, basis2, 2, 1, intmap, 2, 1);
  gm.block(2*dim_l1, 2*dim_l2, dim1-2*dim_l1, dim2-2*dim_l2) = GMGolyCompl(T, basis1, basis2, 2, 2, intmap, 0, 0) + GMGolyCompl(T, basis1, basis2, 2, 2, intmap, 1, 1);
  
  return gm;
}

// GMScalarDerivative, one derivative
Eigen::MatrixXd HArDCore3D::GMScalarDerivative(const Cell & T, const MonomialScalarBasisCell & basis1, const MonomialScalarBasisCell & basis2, const size_t m, MonomialCellIntegralsType mono_int_map)
{
  // Dimension of the gram matrix
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  size_t totaldegree = basis1.max_degree()+basis2.max_degree()-1;
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  // Obtain integration data from IntegrateCellMonomials
  MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
  
  for (size_t i = 0; i < dim1; i++) {
    if(basis1.powers(i)(m) > 0) {
      VectorZd powers1 = basis1.powers(i);
      powers1(m) -= 1;
      for (size_t j = 0; j < dim2; j++) {
        gm(i,j) = basis1.powers(i)(m) * intmap.at(powers1 + basis2.powers(j));
      } // for j
    }   // if
  }     // for i
  
  return gm;
}

// GMScalarDerivative, two derivatives
Eigen::MatrixXd HArDCore3D::GMScalarDerivative(const Cell & T, const MonomialScalarBasisCell & basis1, const MonomialScalarBasisCell & basis2, const size_t m, const size_t l, MonomialCellIntegralsType mono_int_map)
{
  // Dimension of the gram matrix
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  size_t totaldegree = basis1.max_degree()+basis2.max_degree()-2;
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  // Obtain integration data from IntegrateCellMonomials
  MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
  
  for (size_t i = 0; i < dim1; i++) {
    if(basis1.powers(i)(m) > 0) {
      VectorZd powers1 = basis1.powers(i);
      powers1(m) -= 1;
      for (size_t j = 0; j < dim2; j++) {
        if (basis2.powers(j)(l) > 0){
          VectorZd powers2 = basis2.powers(j);
          powers2(l) -= 1;
          gm(i,j) = basis1.powers(i)(m) * basis2.powers(j)(l) * intmap.at(powers1 + powers2);
        }
      } // for j
    }   // if
  }     // for i
  
  return gm;
}
 
// GMRolyScalar, the mth component of RolyCompl and the scalar ancestor of a tensorized basis
Eigen::MatrixXd HArDCore3D::GMRolyComplScalar(const Cell & T, const RolyComplBasisCell & rolycompl_basis, const MonomialScalarBasisCell & mono_basis, const size_t m, MonomialCellIntegralsType mono_int_map)
{
  // Dimension of the gram matrix
  size_t dim1 = rolycompl_basis.dimension();
  size_t dim2 = mono_basis.dimension();
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);

  // Obtain integration data from IntegrateFaceMonomials
  size_t totaldegree = rolycompl_basis.max_degree() + mono_basis.max_degree();
  MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

  VectorZd powers = VectorZd::Zero(3);
  powers(m) = 1;
  for (size_t i = 0; i < dim1; i++) {
    for (size_t j = 0; j < dim2; j++) {
      gm(i,j) = intmap.at(powers + rolycompl_basis.powers(i) + mono_basis.powers(j));
    }   // for j
  }     // for i
  
  return gm;
}

// GMGolyComplScalar, the (optionally k1th derivative of the) sth section of GolyCompl and the (optionally k2th derivative of the) scalar ancestor of a tensorized basis (optionally multiplied by the mth variable)
Eigen::MatrixXd HArDCore3D::GMGolyComplScalar(const Cell & T, const GolyComplBasisCell & golycompl_basis, const MonomialScalarBasisCell & mono_basis, const size_t s, MonomialCellIntegralsType intmap, const size_t m, const size_t k1, const size_t k2)
{
  // Dimension of the gram matrix
  size_t dim_l1 = golycompl_basis.dimPkmo();
  size_t dim1 = dim_l1;
  if (s == 2) dim1 = golycompl_basis.dimension()-2*dim_l1;
  size_t dim2 = mono_basis.dimension();
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  VectorZd powers = VectorZd::Zero(3);
  if (m < 3) powers(m) = 1;
  for (size_t i = 0; i < dim1; i++) {
    if (k1 > 2) {
      for (size_t j = 0; j < dim2; j++) {
        if (k2 > 2) {
          gm(i,j) = intmap.at(powers + golycompl_basis.powers(s*dim_l1+i) + mono_basis.powers(j));
        } else if (mono_basis.powers(j)(k2) > 0) {
          VectorZd powers2 = mono_basis.powers(j);
          powers2(k2) -= 1;
          gm(i,j) = mono_basis.powers(j)(k2) / T.diam() * intmap.at(powers + golycompl_basis.powers(s*dim_l1+i) + powers2);
        }
      } // for j
    } else if (golycompl_basis.powers(s*dim_l1+i)(k1) > 0) {
      VectorZd powers1 = golycompl_basis.powers(s*dim_l1+i);
      powers1(k1) -= 1;
      for (size_t j = 0; j < dim2; j++) {
        if (k2 > 2) {
          gm(i,j) = golycompl_basis.powers(s*dim_l1+i)(k1) / T.diam() * intmap.at(powers + powers1 + mono_basis.powers(j));
        } else if (mono_basis.powers(j)(k2) > 0) {
          VectorZd powers2 = mono_basis.powers(j);
          powers2(k2) -= 1;
          gm(i,j) = golycompl_basis.powers(s*dim_l1+i)(k1) * mono_basis.powers(j)(k2) / std::pow(T.diam(),2) * intmap.at(powers + powers1 + powers2);
        }
      } // for j
    }   // if
  }     // for i
  
  return gm;
}

// GMGolyCompl, the (optionally k1th derivative of the) s1th section of GolyCompl (optionally multiplied by the m1th variable) and the (optionally k2th derivative of the) s2th section of GolyCompl (optionally multiplied by the m2th variable)
Eigen::MatrixXd HArDCore3D::GMGolyCompl(const Cell & T, const GolyComplBasisCell & basis1, const GolyComplBasisCell & basis2, const size_t s1, const size_t s2, MonomialCellIntegralsType intmap, const size_t m1, const size_t m2, const size_t k1, const size_t k2)
{
  // Dimension of the gram matrix
  size_t dim_l1 = basis1.dimPkmo();
  size_t dim_l2 = basis2.dimPkmo();
  size_t dim1 = dim_l1;
  size_t dim2 = dim_l2;
  if (s1 == 2) dim1 = basis1.dimension()-2*dim_l1;
  if (s2 == 2) dim2 = basis2.dimension()-2*dim_l2;
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  VectorZd powers = VectorZd::Zero(3);
  if (m1 < 3) powers(m1) = 1;
  if (m2 < 3) powers(m2) += 1;
  for (size_t i = 0; i < dim1; i++) {
    if (k1 > 2) {
      for (size_t j = 0; j < dim2; j++) {
        if (k2 > 2) {
          gm(i,j) = intmap.at(powers + basis1.powers(s1*dim_l1+i) + basis2.powers(s2*dim_l2+j));
        } else if (basis2.powers(s2*dim_l2+j)(k2) > 0) {
          VectorZd powers2 = basis2.powers(s2*dim_l2+j);
          powers2(k2) -= 1;
          gm(i,j) = basis2.powers(s2*dim_l2+j)(k2) / T.diam() * intmap.at(powers + basis1.powers(s1*dim_l1+i) + powers2);
        } // if
      } // for j
    } else if (basis1.powers(s1*dim_l1+i)(k1) > 0) {
      VectorZd powers1 = basis1.powers(s1*dim_l1+i);
      powers1(k1) -= 1;
      for (size_t j = 0; j < dim2; j++) {
        if (k2 > 2) {
          gm(i,j) = basis1.powers(s1*dim_l1+i)(k1) / T.diam() * intmap.at(powers + powers1 + basis2.powers(s2*dim_l2+j));
        } else if (basis2.powers(s2*dim_l2+j)(k2) > 0) {
          VectorZd powers2 = basis2.powers(s2*dim_l2+j);
          powers2(k2) -= 1;
          gm(i,j) = basis1.powers(s1*dim_l1+i)(k1) * basis2.powers(s2*dim_l2+j)(k2) / std::pow(T.diam(), 2) * intmap.at(powers + powers1 + powers2);
        } // if
      } // for j
    }   // if
  }     // for i
  
  return gm;
}

/// Compute the Gram Matrix of the curl of a GolyCompl basis and the curl of another GolyCompl basis
Eigen::MatrixXd HArDCore3D::GramMatrixCurlCurl(const Cell& T, const GolyComplBasisCell & basis1, const GolyComplBasisCell & basis2, MonomialCellIntegralsType mono_int_map)
{
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  Eigen::MatrixXd gm(dim1, dim2);

  // Integrals of monomials
  size_t totaldegree = basis1.max_degree()+basis2.max_degree()-2;
  MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
  
  size_t dim_l1 = basis1.dimPkmo();
  size_t dim_l2 = basis2.dimPkmo();
  
  gm.block(0, 0, dim_l1, dim_l2) = GMGolyCompl(T, basis1, basis2, 0, 0, intmap, 1, 1, 0, 0) + GMGolyCompl(T, basis1, basis2, 0, 0, intmap, 2, 2, 0, 0)
                                 + GMGolyCompl(T, basis1, basis2, 0, 0, intmap, 1, 1, 1, 1) + GMGolyCompl(T, basis1, basis2, 0, 0, intmap, 1, 2, 1, 2)
                                 + GMGolyCompl(T, basis1, basis2, 0, 0, intmap, 2, 1, 2, 1) + GMGolyCompl(T, basis1, basis2, 0, 0, intmap, 2, 2, 2, 2)
                                 + 2/T.diam() * (GMGolyCompl(T, basis1, basis2, 0, 0, intmap, 1, 3, 1) + GMGolyCompl(T, basis1, basis2, 0, 0, intmap, 3, 1, 3, 1)
                                               + GMGolyCompl(T, basis1, basis2, 0, 0, intmap, 2, 3, 2) + GMGolyCompl(T, basis1, basis2, 0, 0, intmap, 3, 2, 3, 2))
                                 + 4/std::pow(T.diam(),2) * GMGolyCompl(T, basis1, basis2, 0, 0, intmap);
  gm.block(0, dim_l2, dim_l1, dim_l2) = GMGolyCompl(T, basis1, basis2, 0, 1, intmap, 2, 2, 0, 1)
                                      - GMGolyCompl(T, basis1, basis2, 0, 1, intmap, 1, 0, 0, 0) - GMGolyCompl(T, basis1, basis2, 0, 1, intmap, 1, 2, 0, 2)
                                      - GMGolyCompl(T, basis1, basis2, 0, 1, intmap, 1, 0, 1, 1) - GMGolyCompl(T, basis1, basis2, 0, 1, intmap, 2, 0, 2, 1)
                                      - 2/T.diam() * (GMGolyCompl(T, basis1, basis2, 0, 1, intmap, 1, 3, 0) + GMGolyCompl(T, basis1, basis2, 0, 1, intmap, 3, 0, 3, 1));
  gm.block(0, 2*dim_l2, dim_l1, dim2-2*dim_l2) = GMGolyCompl(T, basis1, basis2, 0, 2, intmap, 1, 1, 0, 2)
                                               - GMGolyCompl(T, basis1, basis2, 0, 2, intmap, 2, 0, 0, 0) - GMGolyCompl(T, basis1, basis2, 0, 2, intmap, 2, 1, 0, 1)
                                               - GMGolyCompl(T, basis1, basis2, 0, 2, intmap, 1, 0, 1, 2) - GMGolyCompl(T, basis1, basis2, 0, 2, intmap, 2, 0, 2, 2)
                                               - 2/T.diam() * (GMGolyCompl(T, basis1, basis2, 0, 2, intmap, 2, 3, 0) + GMGolyCompl(T, basis1, basis2, 0, 2, intmap, 3, 0, 3, 2));
  gm.block(dim_l1, 0, dim_l1, dim_l2) = GMGolyCompl(T, basis1, basis2, 1, 0, intmap, 2, 2, 1, 0)
                                      - GMGolyCompl(T, basis1, basis2, 1, 0, intmap, 0, 1, 0, 0) - GMGolyCompl(T, basis1, basis2, 1, 0, intmap, 2, 1, 2, 0)
                                      - GMGolyCompl(T, basis1, basis2, 1, 0, intmap, 0, 1, 1, 1) - GMGolyCompl(T, basis1, basis2, 1, 0, intmap, 0, 2, 1, 2)
                                      - 2/T.diam() * (GMGolyCompl(T, basis1, basis2, 1, 0, intmap, 0, 3, 1) + GMGolyCompl(T, basis1, basis2, 1, 0, intmap, 3, 1, 3, 0));
  gm.block(dim_l1, dim_l2, dim_l1, dim_l2) = GMGolyCompl(T, basis1, basis2, 1, 1, intmap, 0, 0, 1, 1) + GMGolyCompl(T, basis1, basis2, 1, 1, intmap, 2, 2, 1, 1)
                                           + GMGolyCompl(T, basis1, basis2, 1, 1, intmap, 0, 0, 0, 0) + GMGolyCompl(T, basis1, basis2, 1, 1, intmap, 0, 2, 0, 2)
                                           + GMGolyCompl(T, basis1, basis2, 1, 1, intmap, 2, 0, 2, 0) + GMGolyCompl(T, basis1, basis2, 1, 1, intmap, 2, 2, 2, 2)
                                           + 2/T.diam() * (GMGolyCompl(T, basis1, basis2, 1, 1, intmap, 0, 3, 0) + GMGolyCompl(T, basis1, basis2, 1, 1, intmap, 3, 0, 3, 0)
                                                         + GMGolyCompl(T, basis1, basis2, 1, 1, intmap, 2, 3, 2) + GMGolyCompl(T, basis1, basis2, 1, 1, intmap, 3, 2, 3, 2))
                                           + 4/std::pow(T.diam(),2) * GMGolyCompl(T, basis1, basis2, 1, 1, intmap);
  gm.block(dim_l1, 2*dim_l2, dim_l1, dim2-2*dim_l2) = GMGolyCompl(T, basis1, basis2, 1, 2, intmap, 0, 0, 1, 2)
                                                    - GMGolyCompl(T, basis1, basis2, 1, 2, intmap, 0, 1, 0, 2) - GMGolyCompl(T, basis1, basis2, 1, 2, intmap, 2, 1, 2, 2)
                                                    - GMGolyCompl(T, basis1, basis2, 1, 2, intmap, 2, 0, 1, 0) - GMGolyCompl(T, basis1, basis2, 1, 2, intmap, 2, 1, 1, 1)
                                                    - 2/T.diam() * (GMGolyCompl(T, basis1, basis2, 1, 2, intmap, 2, 3, 1) + GMGolyCompl(T, basis1, basis2, 1, 2, intmap, 3, 1, 3, 2));
  gm.block(2*dim_l1, 0, dim1-2*dim_l1, dim_l2) = GMGolyCompl(T, basis1, basis2, 2, 0, intmap, 1, 1, 2, 0)
                                               - GMGolyCompl(T, basis1, basis2, 2, 0, intmap, 0, 2, 0, 0) - GMGolyCompl(T, basis1, basis2, 2, 0, intmap, 1, 2, 1, 0)
                                               - GMGolyCompl(T, basis1, basis2, 2, 0, intmap, 0, 1, 2, 1) - GMGolyCompl(T, basis1, basis2, 2, 0, intmap, 0, 2, 2, 2)
                                               - 2/T.diam() * (GMGolyCompl(T, basis1, basis2, 2, 0, intmap, 0, 3, 2) + GMGolyCompl(T, basis1, basis2, 2, 0, intmap, 3, 2, 3, 0));
  gm.block(2*dim_l1, dim_l2, dim1-2*dim_l1, dim_l2) = GMGolyCompl(T, basis1, basis2, 2, 1, intmap, 0, 0, 2, 1)
                                                    - GMGolyCompl(T, basis1, basis2, 2, 1, intmap, 0, 2, 0, 1) - GMGolyCompl(T, basis1, basis2, 2, 1, intmap, 1, 2, 1, 1)
                                                    - GMGolyCompl(T, basis1, basis2, 2, 1, intmap, 1, 0, 2, 0) - GMGolyCompl(T, basis1, basis2, 2, 1, intmap, 1, 2, 2, 2)
                                                    - 2/T.diam() * (GMGolyCompl(T, basis1, basis2, 2, 1, intmap, 1, 3, 2) + GMGolyCompl(T, basis1, basis2, 2, 1, intmap, 3, 2, 3, 1));
  gm.block(2*dim_l1, 2*dim_l2, dim1-2*dim_l1, dim2-2*dim_l2) = GMGolyCompl(T, basis1, basis2, 2, 2, intmap, 0, 0, 2, 2) + GMGolyCompl(T, basis1, basis2, 2, 2, intmap, 1, 1, 2, 2)
                                                             + GMGolyCompl(T, basis1, basis2, 2, 2, intmap, 0, 0, 0, 0) + GMGolyCompl(T, basis1, basis2, 2, 2, intmap, 0, 1, 0, 1)
                                                             + GMGolyCompl(T, basis1, basis2, 2, 2, intmap, 1, 0, 1, 0) + GMGolyCompl(T, basis1, basis2, 2, 2, intmap, 1, 1, 1, 1)
                                                             + 2/T.diam() * (GMGolyCompl(T, basis1, basis2, 2, 2, intmap, 0, 3, 0) + GMGolyCompl(T, basis1, basis2, 2, 2, intmap, 3, 0, 3, 0)
                                                                           + GMGolyCompl(T, basis1, basis2, 2, 2, intmap, 1, 3, 1) + GMGolyCompl(T, basis1, basis2, 2, 2, intmap, 3, 1, 3, 1))
                                                             + 4/std::pow(T.diam(),2) * GMGolyCompl(T, basis1, basis2, 2, 2, intmap);
  
  return gm;
}

/// Compute the Gram Matrix of the divergence of a RolyCompl basis and the divergence of another RolyCompl basis
Eigen::MatrixXd HArDCore3D::GramMatrixDivDiv(const Cell& T, const RolyComplBasisCell & basis1, const RolyComplBasisCell & basis2, MonomialCellIntegralsType mono_int_map)
{
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  Eigen::MatrixXd gm(dim1, dim2);

  // Integrals of monomials
  size_t totaldegree = basis1.max_degree()+basis2.max_degree()-2;
  MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
  
  for (size_t i = 0; i < dim1; i++) {
    for (size_t j = 0; j < dim2; j++) {
      gm(i,j) = (3 + basis1.powers(i).sum()) * (3 + basis2.powers(j).sum()) * intmap.at(basis1.powers(i) + basis2.powers(j));
    } // for j
  }   // for i
  
  return gm / std::pow(T.diam(),2);
}

/// Computes the Gram Matrix of a Divergence<RolyCompl> basis and a monomial scalar basis
Eigen::MatrixXd HArDCore3D::GramMatrixDiv(const Cell & T, const RolyComplBasisCell & basis1, const MonomialScalarBasisCell & basis2, MonomialCellIntegralsType mono_int_map)
{
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  Eigen::MatrixXd gm(dim1, dim2);
  
  // Integrals of monomials
  size_t totaldegree = basis1.max_degree()+basis2.max_degree()-1;
  MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
  
  for (size_t i = 0; i < dim1; i++) {
    for (size_t j = 0; j < dim2; j++) {
      gm(i,j) = (3 + basis1.powers(i).sum()) * intmap.at(basis1.powers(i) + basis2.powers(j));
    } // for j
  }   // for i
  
  return gm / T.diam();
};
   
/// Computes the Gram Matrix of the divergence of a RolyCompl Basis and the kth derivative of a monomial basis
Eigen::MatrixXd HArDCore3D::GMRolyComplScalarDiv(const Cell & T, const MonomialScalarBasisCell & mono_basis, const RolyComplBasisCell & rolycompl_basis, const size_t k, MonomialCellIntegralsType intmap)
{
  size_t dim1 = mono_basis.dimension();
  size_t dim2 = rolycompl_basis.dimension();
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1, dim2);
  
  for (size_t i = 0; i < dim1; i++) {
    if (mono_basis.powers(i)(k) > 0) {
      VectorZd powers = mono_basis.powers(i);
      powers (k) -= 1;
      for (size_t j = 0; j < dim2; j++) {
        gm(i,j) = (3 + rolycompl_basis.powers(j).sum()) * mono_basis.powers(i)(k) * intmap.at(powers + rolycompl_basis.powers(j));
      } // for j
    } // if
  }   // for i
  
  return gm / std::pow(T.diam(),2);
};
  
