#include <GMpoly_face.hpp>
#include <max_degrees_quadratures.hpp>

using namespace HArDCore3D;

// IntegrateFaceMonomials_onEdges
std::vector<MonomialFaceIntegralsType> HArDCore3D::IntegrateFaceMonomials_onEdges(const Face & F, const size_t maxdeg) 
  {
  // Number of monomials
  const size_t nb_poly = PolynomialSpaceDimension<Face>::Poly(maxdeg);
  
  std::vector<MonomialFaceIntegralsType> integrals;
  integrals.resize(F.n_edges());
    
  // Dimension of the face
  const size_t d = 2;
  
  // Coordinate transformation data
  VectorRd xF = F.center_mass();
  double hF = F.diam();
  Eigen::Matrix<double, 2, 3> JF;
  JF.row(0) = F.edge(0)->tangent();
  JF.row(1) = F.edge_normal(0);
  JF /= hF;
  
  // Create powers and degrees of all monomials
  std::vector<Eigen::Vector2i> powers = MonomialPowers<Face>::complete(maxdeg);
  
  // Loop over the edges
  for (size_t iE = 0; iE < F.n_edges(); iE++)
  {
    VectorRd v1 = F.edge(iE)->vertex(0)->coords();
    VectorRd v2 = F.edge(iE)->vertex(1)->coords();
    VectorRd xE = F.edge(iE)->center_mass();
    double d1 = (v1-xE).norm();
    double d2 = (v2-xE).norm();

    Eigen::Vector2d v1_trans = JF * (v1-xF);
    Eigen::Vector2d v2_trans = JF * (v2-xF);

    // Loop over all the monomials
    for (size_t m = 0; m < nb_poly; m++)
    {
      // Vertex values
      double int_e_mono = d1 * std::pow(v1_trans.x(),powers[m](0)) * std::pow(v1_trans.y(),powers[m](1));
      int_e_mono += d2 * std::pow(v2_trans.x(),powers[m](0)) * std::pow(v2_trans.y(),powers[m](1));

      // Gradient correction (each partial derivative)
      Eigen::Vector2d xE_trans = JF * (xE-xF);
      for (size_t ip = 0; ip < d; ip++)
      {
        if (powers[m](ip) > 0){
          Eigen::Vector2i powers_diff = Eigen::Vector2i::Zero(d);
          powers_diff(ip) = -1;
          int_e_mono += xE_trans(ip) * powers[m](ip) * integrals[iE][powers[m]+powers_diff];
        }
      }
      
      // Store
      integrals[iE][powers[m]] = int_e_mono / (d-1+powers[m].sum());
    } // for m
  } // for iE      
  
  return integrals;
}

// IntegrateFaceMonomials
MonomialFaceIntegralsType HArDCore3D::IntegrateFaceMonomials(const Face & F, const size_t maxdeg) 
  {
  // Number of monomials
  const size_t nb_poly = PolynomialSpaceDimension<Face>::Poly(maxdeg);
  
  MonomialFaceIntegralsType integrals;
  
  // Dimension of the face
  const size_t d = 2;
  
  // Compute integrals on all the edges
  std::vector<MonomialFaceIntegralsType> integrals_edges = IntegrateFaceMonomials_onEdges(F, maxdeg);
  
  VectorRd xF = F.center_mass();
  
  // Create powers and degrees of all monomials
  std::vector<Eigen::Vector2i> powers = MonomialPowers<Face>::complete(maxdeg);

  // Loop over the monomials
  for (size_t m = 0; m<nb_poly; m++){
    double int_mono = 0.; 
    // Adding edge contributions
    for (size_t iE = 0; iE < F.n_edges(); iE++){
      // Outer unit normal and arbitrary point on edge
      VectorRd nFE = F.edge_normal(iE);
      VectorRd xE = F.edge(iE)->center_mass();
      
      int_mono += abs(nFE.dot(xE-xF)) * integrals_edges[iE][powers[m]];
    } // for iE
    
    // Store
    integrals[powers[m]] = int_mono / (d+powers[m].sum());     
  }
    
  return integrals;
}


// GramMatrix
Eigen::MatrixXd HArDCore3D::GramMatrix(const Face & F, const MonomialScalarBasisFace & basis1, const MonomialScalarBasisFace & basis2, MonomialFaceIntegralsType mono_int_map)
{
  // Dimension of the gram matrix
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  size_t totaldegree = basis1.max_degree()+basis2.max_degree();
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  // Obtain integration data from IntegrateFaceMonomials
  MonomialFaceIntegralsType intmap;
  if (mono_int_map.size()>0){
    intmap = mono_int_map;
  }else{
    intmap = IntegrateFaceMonomials(F, totaldegree);
  }
  
  for (size_t i = 0; i < dim1; i++) {
    for (size_t j = 0; j < dim2; j++) {
      gm(i,j) = intmap[basis1.powers(i) + basis2.powers(j)];
    } // for j
  }   // for i
  
  return gm;
}

/// Computes the Gram Matrix of a pair of RolyCompl bases
Eigen::MatrixXd HArDCore3D::GramMatrix(const Face & F, const RolyComplBasisFace & basis1, const RolyComplBasisFace & basis2, MonomialFaceIntegralsType mono_int_map)
{
  // Dimension of the gram matrix
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  size_t totaldegree = basis1.max_degree()+basis2.max_degree();
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  // Obtain integration data from IntegrateFaceMonomials
  MonomialFaceIntegralsType intmap;
  if (mono_int_map.size()>0){
    intmap = mono_int_map;
  }else{
    intmap = IntegrateFaceMonomials(F, totaldegree);
  }
  
  for (size_t k = 0; k < 2; k++) {
    Eigen::Vector2i powers = Eigen::Vector2i::Zero(2);
    powers(k) = 2;
    for (size_t i = 0; i < dim1; i++) {
      for (size_t j = 0; j < dim2; j++) {
          gm(i,j) += intmap[powers + basis1.powers(i) + basis2.powers(j)];
      } // for j
    }   // for i
  }     // for k
  
  return gm;
}

/// Computes the Gram Matrix of a pair of GolyCompl bases
Eigen::MatrixXd HArDCore3D::GramMatrix(const Face & F, const GolyComplBasisFace & basis1, const GolyComplBasisFace & basis2, MonomialFaceIntegralsType mono_int_map)
{
  return GramMatrix(F, *basis1.rck(), *basis2.rck(), mono_int_map);
}

/// Computes the Gram Matrix of the scalar part of a RolyCompl Basis and a monomial basis with an extra power on the mth variable
Eigen::MatrixXd HArDCore3D::GMRolyComplScalar(const Face & F, const RolyComplBasisFace & rolycompl_basis, const MonomialScalarBasisFace & mono_basis, const size_t m, MonomialFaceIntegralsType intmap )
{
  // Dimension of the gram matrix
  size_t dim1 = rolycompl_basis.dimension();
  size_t dim2 = mono_basis.dimension();
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  Eigen::Vector2i powers = Eigen::Vector2i::Zero(2);
  powers(m) = 1;
  for (size_t i = 0; i < dim1; i++) {
    for (size_t j = 0; j < dim2; j++) {
      gm(i,j) = intmap[powers + rolycompl_basis.powers(i) + mono_basis.powers(j)];
    }   // for j
  }     // for i
  
  return gm;
}

// GMScalarDerivative, one derivative
Eigen::MatrixXd HArDCore3D::GMScalarDerivative(const Face & F, const MonomialScalarBasisFace & basis1, const MonomialScalarBasisFace & basis2, const size_t m, MonomialFaceIntegralsType intmap)
{
  // Dimension of the gram matrix
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  for (size_t i = 0; i < dim1; i++) {
    if(basis1.powers(i)(m) > 0) {
      Eigen::Vector2i powers1 = basis1.powers(i);
      powers1(m) -= 1;
      for (size_t j = 0; j < dim2; j++) {
        gm(i,j) = basis1.powers(i)(m) * intmap[powers1 + basis2.powers(j)];
      } // for j
    }   // if
  }     // for i
  
  return gm;
}

// GMScalarDerivative, two derivatives
Eigen::MatrixXd HArDCore3D::GMScalarDerivative(const Face & F, const MonomialScalarBasisFace & basis1, const MonomialScalarBasisFace & basis2, const size_t m, const size_t l, MonomialFaceIntegralsType intmap)
{
  // Dimension of the gram matrix
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  for (size_t i = 0; i < dim1; i++) {
    if(basis1.powers(i)(m) > 0) {
      Eigen::Vector2i powers1 = basis1.powers(i);
      powers1(m) -= 1;
      for (size_t j = 0; j < dim2; j++) {
        if (basis2.powers(j)(l) > 0){
          Eigen::Vector2i powers2 = basis2.powers(j);
          powers2(l) -= 1;
          gm(i,j) = basis1.powers(i)(m) * basis2.powers(j)(l) * intmap[powers1 + powers2];
        }
      } // for j
    }   // if
  }     // for i
  
  return gm;
}



// Gram Matrix of a Divergence<RolyCompl> basis and a monomial scalar basis
Eigen::MatrixXd HArDCore3D::GramMatrixDiv(const Face & F, const RolyComplBasisFace & basis1, const MonomialScalarBasisFace & basis2, MonomialFaceIntegralsType mono_int_map)
{
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  Eigen::MatrixXd gm(dim1, dim2);
  
  // Integrals of monomials
  size_t totaldegree = basis1.max_degree()+basis2.max_degree()-1;
  MonomialFaceIntegralsType intmap;
  if (mono_int_map.size()>0){
    intmap = mono_int_map;
  }else{
    intmap = IntegrateFaceMonomials(F, totaldegree);
  }
  
  for (size_t i = 0; i < dim1; i++) {
    for (size_t j = 0; j < dim2; j++) {
      gm(i,j) = (2 + basis1.powers(i).sum()) * intmap[basis1.powers(i) + basis2.powers(j)];
    } // for j
  }   // for i
  
  return gm / F.diam();
};



