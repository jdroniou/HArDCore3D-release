// Core data structures and methods required to implement the discrete de Rham sequence in 3D
//
// Provides:
//  - Full and partial polynomial spaces on the element, faces, and edges
//  - Generic routines to create quadrature nodes over cells and faces of the mesh
//
// Authors: Daniele Di Pietro (daniele.di-pietro@umontpellier.fr)
//          Jerome Droniou (jerome.droniou@monash.edu)

/*
 *
 *      This library was developed around HHO methods, although some parts of it have a more
 * general purpose. If you use this code or part of it in a scientific publication, 
 * please mention the following book as a reference for the underlying principles
 * of HHO schemes:
 *
 * The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications. 
 *  D. A. Di Pietro and J. Droniou. Modeling, Simulation and Applications, vol. 19. 
 *  Springer International Publishing, 2020, xxxi + 525p. doi: 10.1007/978-3-030-37203-3. 
 *  url: https://hal.archives-ouvertes.fr/hal-02151813.
 *
 */

#ifndef BASIS_HPP
#define BASIS_HPP

#include <boost/multi_array.hpp>

#include <mesh.hpp>
#include <iostream>

#include <polynomialspacedimension.hpp>

#include <quadraturerule.hpp>

namespace HArDCore3D
{

  /*!	
 * @defgroup Basis 
 * @brief Classes and functions for polynomial basis creation and manipulation
 */

  /*!
   *	\addtogroup Basis
   * @{
   */

  /// Dimension, and generic types for vector in correct dimension (makes it easier to translate a code between 2D and 3D)
  constexpr int dimspace = 3;
  typedef Eigen::Matrix3d MatrixRd;

  template <typename T>
  using FType = std::function<T(const VectorRd&)>; ///< type for function of point. T is the return type of the function

  template <typename T>
  using CellFType = std::function<T(const VectorRd&, const Cell *)>; ///< type for function of point. T is the return type of the function

  template <typename T>
  using BasisQuad = boost::multi_array<T, 2>; ///< type for bases evaluated on quadrature nodes

  enum TensorRankE
  {
    Scalar = 0,
    Vector = 1,
    Matrix = 2
  };

  //-----------------------------------------------------------------------------
  //                Powers of monomial basis
  //-----------------------------------------------------------------------------

  /// Compute vectors listing the powers of monomial basis functions (for a cell or face, only specializations are relevant) up to a certain degree
  template<typename GeometricSupport>
  struct MonomialPowers
  {
    // Only specializations are relevant
  };
  
  template<>
  struct MonomialPowers<Cell>
  {
    // Powers of homogeneous polynomials of degree L
    static std::vector<VectorZd> homogeneous(const size_t L) 
    {
      std::vector<VectorZd> powers;
      powers.reserve( PolynomialSpaceDimension<Cell>::Poly(L) - PolynomialSpaceDimension<Cell>::Poly(L-1) );
      for (size_t i = 0; i <= L; i++)
      {
        for (size_t j = 0; i + j <= L; j++)
        {
          powers.push_back(VectorZd(i, j, L - i - j));
        } // for j
      }   // for i
      return powers;
    }
  
    // Powers for degrees from 0 to degree
    static std::vector<VectorZd> complete(const size_t degree)
    {
      std::vector<VectorZd> powers;
      powers.reserve( PolynomialSpaceDimension<Cell>::Poly(degree) );
      for (size_t l = 0; l <= degree; l++)
      {
        std::vector<VectorZd> pow_hom = homogeneous(l);
        for (VectorZd p : pow_hom){
          powers.push_back(p);
        }
      }     // for l
      return powers;
    }
  };
  
  template<>
  struct MonomialPowers<Face>
  {
    // Powers of homogeneous polynomials of degree L
    static std::vector<Eigen::Vector2i> homogeneous(const size_t L) 
    {
      std::vector<Eigen::Vector2i> powers;
      powers.reserve( PolynomialSpaceDimension<Face>::Poly(L) - PolynomialSpaceDimension<Face>::Poly(L-1) );
      for (size_t i = 0; i <= L; i++)
      {
        powers.push_back(Eigen::Vector2i(i, L - i));
      }   // for i
      return powers;
    }

    // Powers for degrees from 0 to degree
    static std::vector<Eigen::Vector2i> complete(size_t degree)
    {
      std::vector<Eigen::Vector2i> powers;
      powers.reserve( PolynomialSpaceDimension<Face>::Poly(degree) );
      for (size_t l = 0; l <= degree; l++)
      {
        std::vector<Eigen::Vector2i> pow_hom = homogeneous(l);
        for (Eigen::Vector2i p : pow_hom){
          powers.push_back(p);
        }
      }   // for l

      return powers;
    }
  };

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  //          SCALAR MONOMIAL BASIS
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------

  /// Scalar monomial basis on a cell
  class MonomialScalarBasisCell
  {
  public:
    typedef double FunctionValue;
    typedef VectorRd GradientValue;
    typedef VectorRd CurlValue;
    typedef double DivergenceValue;
    typedef MatrixRd HessianValue;

    typedef Cell GeometricSupport;

    constexpr static const TensorRankE tensorRank = Scalar;
    constexpr static const bool hasAncestor = false;
    static const bool hasFunction = true;
    static const bool hasGradient = true;
    static const bool hasCurl = false;
    static const bool hasDivergence = false;
    static const bool hasHessian = true;
    static const bool hasCurlCurl = false;

    /// Constructor
    MonomialScalarBasisCell(
        const Cell &T, ///< A mesh cell
        size_t degree  ///< The maximum polynomial degree to be considered
    );

    /// Compute the dimension of the basis
    inline size_t dimension() const
    {
      return (m_degree >= 0 ? (m_degree * (m_degree + 1) * (2 * m_degree + 1) + 9 * m_degree * (m_degree + 1) + 12 * (m_degree + 1)) / 12 : 0);
    }

    /// Evaluate the i-th basis function at point x
    FunctionValue function(size_t i, const VectorRd &x) const;

    /// Evaluate the gradient of the i-th basis function at point x
    GradientValue gradient(size_t i, const VectorRd &x) const;
    
    /// Evaluate the hessian of the i-th basis function at point x
    HessianValue hessian(size_t i, const VectorRd &x) const;

    /// Returns the maximum degree of the basis functions
    inline size_t max_degree() const
    {
      return m_degree;
    }

    /// Returns the powers of the i-th basis function (its degree can be found using powers(i).sum())
    inline VectorZd powers(size_t i) const
    {
      return m_powers[i];
    }
    
  private:
    /// Coordinate transformation
    inline VectorRd _coordinate_transform(const VectorRd &x) const
    {
      return (x - m_xT) / m_hT;
    }

    Cell m_T;
    size_t m_degree;
    VectorRd m_xT;
    double m_hT;
    std::vector<VectorZd> m_powers;
  };


  //------------------------------------------------------------------------------

  /// Scalar monomial basis on a face
  class MonomialScalarBasisFace
  {
  public:
    typedef double FunctionValue;
    typedef VectorRd GradientValue;
    typedef VectorRd CurlValue;
    typedef double DivergenceValue;
    typedef MatrixRd HessianValue;

    typedef Face GeometricSupport;

    typedef Eigen::Matrix<double, 2, dimspace> JacobianType;

    constexpr static const TensorRankE tensorRank = Scalar;
    constexpr static const bool hasAncestor = false;
    static const bool hasFunction = true;
    static const bool hasGradient = true;
    static const bool hasCurl = true;
    static const bool hasDivergence = false;
    static const bool hasHessian = true;
    static const bool hasCurlCurl = false;
    
    /// Constructor
    MonomialScalarBasisFace(
        const Face &F, ///< A mesh face
        size_t degree  ///< The maximum polynomial degree to be considered
    );

    /// Dimension of the basis
    inline size_t dimension() const
    {
      return (m_degree + 1) * (m_degree + 2) / 2;
    }

    /// Evaluate the i-th basis function at point x
    FunctionValue function(size_t i, const VectorRd &x) const;

    /// Evaluate the gradient of the i-th basis function at point x
    GradientValue gradient(size_t i, const VectorRd &x) const;

    /// Evaluate the two-dimensional curl of the i-th basis function at point x
    CurlValue curl(size_t i, const VectorRd &x) const;

    /// Evaluate the Hessian of the i-th basis function at point x
    HessianValue hessian(size_t i, const VectorRd &x) const;

    /// Returns the maximum degree of the basis functions
    inline size_t max_degree() const
    {
      return m_degree;
    }

    /// Return the normal to the face used in the computation of the curl
    inline const VectorRd &normal() const
    {
      return m_nF;
    }

    /// Return the Jacobian of the coordinate system transformation
    inline const JacobianType &jacobian() const
    {
      return m_jacobian;
    }
    
    /// Return the system of coordinates (basis in rows) on the face
    inline const JacobianType coordinates_system() const
    {
      return m_jacobian * m_hF;
    }

    /// Returns the powers of the i-th basis function (its degree can be found using powers(i).sum())
    inline Eigen::Vector2i powers(size_t i) const
    {
      return m_powers[i];
    }
    
  private:
    /// Coordinate transformation
    inline Eigen::Vector2d _coordinate_transform(const VectorRd &x) const
    {
      return m_jacobian * (x - m_xF);
    }

    size_t m_degree;
    VectorRd m_xF;
    double m_hF;
    VectorRd m_nF;
    JacobianType m_jacobian;
    std::vector<Eigen::Vector2i> m_powers;
  };

  //------------------------------------------------------------------------------

  /// Scalar monomial basis on an edge
  class MonomialScalarBasisEdge
  {
  public:
    typedef double FunctionValue;
    typedef VectorRd GradientValue;
    typedef VectorRd CurlValue;
    typedef double DivergenceValue;
    typedef double HessianValue;
    
    typedef Edge GeometricSupport;

    constexpr static const TensorRankE tensorRank = Scalar;
    constexpr static const bool hasAncestor = false;
    static const bool hasFunction = true;
    static const bool hasGradient = true;
    static const bool hasCurl = false;
    static const bool hasDivergence = false;
    static const bool hasHessian = false;
    static const bool hasCurlCurl = false;
    
    /// Constructor
    MonomialScalarBasisEdge(
        const Edge &E, ///< A mesh edge
        size_t degree  ///< The maximum polynomial degree to be considered
    );

    /// Dimension of the basis
    inline size_t dimension() const
    {
      return m_degree + 1;
    }

    /// Evaluate the i-th basis function at point x
    FunctionValue function(size_t i, const VectorRd &x) const;

    /// Evaluate the gradient of the i-th basis function at point x
    GradientValue gradient(size_t i, const VectorRd &x) const;

    /// Returns the maximum degree of the basis functions
    inline size_t max_degree() const
    {
      return m_degree;
    }

  private:
    inline double _coordinate_transform(const VectorRd &x) const
    {
      return (x - m_xE).dot(m_tE) / m_hE;
    }

    size_t m_degree;
    VectorRd m_xE;
    double m_hE;
    VectorRd m_tE;
  };

  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //          DERIVED BASIS
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  
  //----------------------------FAMILY------------------------------------------
  //
  /// Family of functions expressed as linear combination of the functions of a given basis
  /** If \f$(f_1,...,f_r)\f$ is the basis, then the family is \f$(\phi_1,...,\phi_l)\f$ where
       \f$\phi_i = \sum_j M_{ij} f_j\f$.
   The matrix \f$M\f$ is the member m_matrix. */
  template <typename BasisType>
  class Family
  {
  public:
    typedef typename BasisType::FunctionValue FunctionValue;
    typedef typename BasisType::GradientValue GradientValue;
    typedef typename BasisType::CurlValue CurlValue;
    typedef typename BasisType::DivergenceValue DivergenceValue;
    typedef typename BasisType::HessianValue HessianValue;
    
    typedef typename BasisType::GeometricSupport GeometricSupport;

    constexpr static const TensorRankE tensorRank = BasisType::tensorRank;
    constexpr static const bool hasAncestor = true;
    static const bool hasFunction = BasisType::hasFunction;
    static const bool hasGradient = BasisType::hasGradient;
    static const bool hasCurl = BasisType::hasCurl;
    static const bool hasDivergence = BasisType::hasDivergence;
    static const bool hasHessian = BasisType::hasHessian;
    static const bool hasCurlCurl = BasisType::hasCurlCurl;

    typedef BasisType AncestorType;

    /// Constructor
    Family(
        const BasisType &basis,       ///< The basis in which the family is expressed
        const Eigen::MatrixXd &matrix ///< The coefficient matrix whose i-th line contains the coefficient of the expansion of the i-th function of the family in the basis
        )
        : m_basis(basis),
          m_matrix(matrix)
    {
      assert((size_t)matrix.cols() == basis.dimension() || "Inconsistent family initialization");
    }

    /// Dimension of the family. This is actually the number of functions in the family, not necessarily linearly independent
    inline size_t dimension() const
    {
      return m_matrix.rows();
    }

    /// Evaluate the i-th function at point x
    FunctionValue function(size_t i, const VectorRd &x) const
    {
      static_assert(hasFunction, "Call to function() not available");

      FunctionValue f = m_matrix(i, 0) * m_basis.function(0, x);
      for (auto j = 1; j < m_matrix.cols(); j++)
      {
        f += m_matrix(i, j) * m_basis.function(j, x);
      } // for j
      return f;
    }

    /// Evaluate the i-th function at a quadrature point iqn, knowing all the values of ancestor basis functions at the quadrature nodes (provided by eval_quad)
    FunctionValue function(size_t i, size_t iqn, const boost::multi_array<FunctionValue, 2> &ancestor_value_quad) const
    {
      static_assert(hasFunction, "Call to function() not available");

      FunctionValue f = m_matrix(i, 0) * ancestor_value_quad[0][iqn];
      for (auto j = 1; j < m_matrix.cols(); j++)
      {
        f += m_matrix(i, j) * ancestor_value_quad[j][iqn];
      } // for j
      return f;
    }
   
    
    /// Evaluate the gradient of the i-th function at point x
    GradientValue gradient(size_t i, const VectorRd &x) const
    {
      static_assert(hasGradient, "Call to gradient() not available");

      GradientValue G = m_matrix(i, 0) * m_basis.gradient(0, x);
      for (auto j = 1; j < m_matrix.cols(); j++)
      {
        G += m_matrix(i, j) * m_basis.gradient(j, x);
      } // for j
      return G;
    }

    /// Evaluate the gradient of the i-th function at a quadrature point iqn, knowing all the gradients of ancestor basis functions at the quadrature nodes (provided by eval_quad)
    GradientValue gradient(size_t i, size_t iqn, const boost::multi_array<GradientValue, 2> &ancestor_gradient_quad) const
    {
      static_assert(hasGradient, "Call to gradient() not available");

      GradientValue G = m_matrix(i, 0) * ancestor_gradient_quad[0][iqn];
      for (auto j = 1; j < m_matrix.cols(); j++)
      {
        G += m_matrix(i, j) * ancestor_gradient_quad[j][iqn];
      } // for j
      return G;
    }

    /// Evaluate the curl of the i-th function at point x
    CurlValue curl(size_t i, const VectorRd &x) const
    {
      static_assert(hasCurl, "Call to curl() not available");

      CurlValue C = m_matrix(i, 0) * m_basis.curl(0, x);
      for (auto j = 1; j < m_matrix.cols(); j++)
      {
        C += m_matrix(i, j) * m_basis.curl(j, x);
      } // for j
      return C;
    }

    /// Evaluate the curl of the i-th function at a quadrature point iqn, knowing all the curls of ancestor basis functions at the quadrature nodes (provided by eval_quad)
    CurlValue curl(size_t i, size_t iqn, const boost::multi_array<CurlValue, 2> &ancestor_curl_quad) const
    {
      static_assert(hasCurl, "Call to curl() not available");

      CurlValue C = m_matrix(i, 0) * ancestor_curl_quad[0][iqn];
      for (auto j = 1; j < m_matrix.cols(); j++)
      {
        C += m_matrix(i, j) * ancestor_curl_quad[j][iqn];
      } // for j
      return C;
    }

    /// Evaluate the divergence of the i-th function at point x
    DivergenceValue divergence(size_t i, const VectorRd &x) const
    {
      static_assert(hasDivergence, "Call to divergence() not available");

      DivergenceValue D = m_matrix(i, 0) * m_basis.divergence(0, x);
      for (auto j = 1; j < m_matrix.cols(); j++)
      {
        D += m_matrix(i, j) * m_basis.divergence(j, x);
      } // for j
      return D;
    }

    /// Evaluate the divergence of the i-th function at a quadrature point iqn, knowing all the divergences of ancestor basis functions at the quadrature nodes (provided by eval_quad)
    DivergenceValue divergence(size_t i, size_t iqn, const boost::multi_array<DivergenceValue, 2> &ancestor_divergence_quad) const
    {
      static_assert(hasDivergence, "Call to divergence() not available");

      DivergenceValue D = m_matrix(i, 0) * ancestor_divergence_quad[0][iqn];
      for (auto j = 1; j < m_matrix.cols(); j++)
      {
        D += m_matrix(i, j) * ancestor_divergence_quad[j][iqn];
      } // for j
      return D;
    }

    /// Evaluate the hessian of the i-th function at point x
    HessianValue hessian(size_t i, const VectorRd &x) const
    {
      static_assert(hasHessian, "Call to hessian() not available");

      HessianValue H = m_matrix(i, 0) * m_basis.hessian(0, x);
      for (auto j = 1; j < m_matrix.cols(); j++)
      {
        H += m_matrix(i, j) * m_basis.hessian(j, x);
      } // for j
      return H;
    }

    /// Evaluate the hessian of the i-th function at a quadrature point iqn, knowing all the hessian of ancestor basis functions at the quadrature nodes (provided by eval_quad)
    HessianValue hessian(size_t i, size_t iqn, const boost::multi_array<HessianValue, 2> &ancestor_hessian_quad) const
    {
      static_assert(hasHessian, "Call to hessian() not available");

      HessianValue H = m_matrix(i, 0) * ancestor_hessian_quad[0][iqn];
      for (auto j = 1; j < m_matrix.cols(); j++)
      {
        H += m_matrix(i, j) * ancestor_hessian_quad[j][iqn];
      } // for j
      return H;
    }

    /// Return the coefficient matrix
    inline const Eigen::MatrixXd & matrix() const
    {
      return m_matrix;
    }

    /// Return the ancestor
    constexpr inline const BasisType &ancestor() const
    {
      return m_basis;
    }

    /// Returns the maximum degree of the basis functions
    inline size_t max_degree() const
    {
      return m_basis.max_degree();
    }

  private:
    BasisType m_basis;
    Eigen::MatrixXd m_matrix;
  };

  //----------------------TENSORIZED--------------------------------------------------------

  /// Vector family obtained by tensorization of a scalar family
  /** The tensorization is done the following way: if \f$(f_1,...,f_r)\f$ is the family of scalar functions,
   the tensorized family of rank N is given by (where all vectors are columns of size N):
  
     \f$\left(\begin{array}{c}f_1\\0\\\vdots\\0\end{array}\right)\f$;
     \f$\left(\begin{array}{c}f_2\\0\\\vdots\\0\end{array}\right)\f$;...;
     \f$\left(\begin{array}{c}f_r\\0\\\vdots\\0\end{array}\right)\f$;
     \f$\left(\begin{array}{c}0\\f_1\\0\\\vdots\\0\end{array}\right)\f$;
     \f$\left(\begin{array}{c}0\\f_2\\0\\\vdots\\0\end{array}\right)\f$;...;
     \f$\left(\begin{array}{c}0\\f_r\\0\\\vdots\\0\end{array}\right)\f$;...;
     \f$\left(\begin{array}{c}0\\\vdots\\0\\f_1\end{array}\right)\f$;...;
     \f$\left(\begin{array}{c}0\\\vdots\\0\\f_r\end{array}\right)\f$
  
   The gradient values are therefore matrices of size N*r, where the gradients of the scalar functions are put
   in rows:
  
     \f$\left(\begin{array}{c}(\nabla f_1)^t\\0\\\vdots\\0\end{array}\right)\f$;
     \f$\left(\begin{array}{c}(\nabla f_2)^t\\0\\\vdots\\0\end{array}\right)\f$;...;
     \f$\left(\begin{array}{c}(\nabla f_r)^t\\0\\\vdots\\0\end{array}\right)\f$;
     \f$\left(\begin{array}{c}0\\(\nabla f_1)^t\\0\\\vdots\\0\end{array}\right)\f$;
     \f$\left(\begin{array}{c}0\\(\nabla f_2)^t\\0\\\vdots\\0\end{array}\right)\f$;...;
     \f$\left(\begin{array}{c}0\\(\nabla f_r)^t\\0\\\vdots\\0\end{array}\right)\f$;...;
     \f$\left(\begin{array}{c}0\\\vdots\\0\\(\nabla f_1)^t\end{array}\right)\f$;...;
     \f$\left(\begin{array}{c}0\\\vdots\\0\\(\nabla f_r)^t\end{array}\right)\f$
  */
  template <typename ScalarFamilyType, size_t N>
  class TensorizedVectorFamily
  {
  public:
    typedef typename Eigen::Matrix<double, N, 1> FunctionValue;
    typedef typename Eigen::Matrix<double, N, dimspace> GradientValue;
    typedef VectorRd CurlValue;
    typedef double DivergenceValue;
    typedef void HessianValue;

    typedef typename ScalarFamilyType::GeometricSupport GeometricSupport;

    constexpr static const TensorRankE tensorRank = Vector;
    constexpr static const bool hasAncestor = true;
    static const bool hasFunction = ScalarFamilyType::hasFunction;
    static const bool hasGradient = ScalarFamilyType::hasGradient;
    // We know how to compute the curl and divergence if gradient is available, and
    // if we tensorize at the dimension of the space
    static const bool hasDivergence = (ScalarFamilyType::hasGradient && N == dimspace);
    static const bool hasCurl = (ScalarFamilyType::hasGradient && N == dimspace);
    static const bool hasHessian = false;
    static const bool hasCurlCurl = (ScalarFamilyType::hasHessian && N == dimspace);

    typedef ScalarFamilyType AncestorType;

    TensorizedVectorFamily(const ScalarFamilyType &scalar_family)
        : m_scalar_family(scalar_family)
    {
      static_assert(ScalarFamilyType::tensorRank == Scalar,
                    "Vector family can only be constructed from scalar families");
    }

    /// Return the dimension of the family
    inline size_t dimension() const
    {
      return m_scalar_family.dimension() * N;
    }

    /// Evaluate the i-th basis function at point x
    FunctionValue function(size_t i, const VectorRd &x) const
    {
      static_assert(hasFunction, "Call to function() not available");

      FunctionValue ek = Eigen::Matrix<double, N, 1>::Zero();
      ek(i / m_scalar_family.dimension()) = 1.;
      return ek * m_scalar_family.function(i % m_scalar_family.dimension(), x);
    }

    /// Evaluate the i-th basis function at a quadrature point iqn, knowing all the values of ancestor basis functions at the quadrature nodes (provided by eval_quad)
    FunctionValue function(size_t i, size_t iqn, const boost::multi_array<double, 2> &ancestor_value_quad) const
    {
      static_assert(hasFunction, "Call to function() not available");

      FunctionValue ek = Eigen::Matrix<double, N, 1>::Zero();
      ek(i / m_scalar_family.dimension()) = 1.;
      return ek * ancestor_value_quad[i % m_scalar_family.dimension()][iqn];
    }

    /// Evaluate the gradient of the i-th basis function at point x
    GradientValue gradient(size_t i, const VectorRd &x) const
    {
      static_assert(hasGradient, "Call to gradient() not available");

      GradientValue G = Eigen::Matrix<double, N, dimspace>::Zero();
      G.row(i / m_scalar_family.dimension()) = m_scalar_family.gradient(i % m_scalar_family.dimension(), x);
      return G;
    }

    /// Evaluate the gradient of the i-th basis function at a quadrature point iqn, knowing all the gradients of ancestor basis functions at the quadrature nodes (provided by eval_quad)
    GradientValue gradient(size_t i, size_t iqn, const boost::multi_array<VectorRd, 2> &ancestor_gradient_quad) const
    {
      static_assert(hasGradient, "Call to gradient() not available");

      GradientValue G = Eigen::Matrix<double, N, dimspace>::Zero();
      G.row(i / m_scalar_family.dimension()) = ancestor_gradient_quad[i % m_scalar_family.dimension()][iqn];
      return G;
    }

    /// Evaluate the curl of the i-th basis function at point x
    CurlValue curl(size_t i, const VectorRd &x) const
    {
      static_assert(hasCurl, "Call to curl() not available");

      VectorRd ek = VectorRd::Zero();
      ek(i / m_scalar_family.dimension()) = 1.;
      return m_scalar_family.gradient(i % m_scalar_family.dimension(), x).cross(ek);
    }

    /// Evaluate the curl of the i-th basis function at a quadrature point iqn, knowing all the gradients of ancestor basis functions at the quadrature nodes (provided by eval_quad)
    CurlValue curl(size_t i, size_t iqn, const boost::multi_array<VectorRd, 2> &ancestor_gradient_quad) const
    {
      static_assert(hasCurl, "Call to curl() not available");

      VectorRd ek = VectorRd::Zero();
      ek(i / m_scalar_family.dimension()) = 1.;
      return ancestor_gradient_quad[i % m_scalar_family.dimension()][iqn].cross(ek);
    }

    /// Evaluate the divergence of the i-th basis function at point x
    DivergenceValue divergence(size_t i, const VectorRd &x) const
    {
      static_assert(hasDivergence, "Call to divergence() not available");

      return m_scalar_family.gradient(i % m_scalar_family.dimension(), x)(i / m_scalar_family.dimension());
    }
    
    /// Evaluate the divergence of the i-th basis function at a quadrature point iqn, knowing all the gradients of ancestor basis functions at the quadrature nodes (provided by eval_quad)
    DivergenceValue divergence(size_t i, size_t iqn, const boost::multi_array<VectorRd, 2> &ancestor_gradient_quad) const
    {
      static_assert(hasDivergence, "Call to divergence() not available");

      return ancestor_gradient_quad[i % m_scalar_family.dimension()][iqn](i / m_scalar_family.dimension());      
    }

    /// Evaluate the curl curl of the i-th basis function at point x. We use the formula curl curl = - Laplace + grad div
    CurlValue curlcurl(size_t i, const VectorRd &x) const
    {
      static_assert(hasCurlCurl, "Call to curlcurl() not available");

      FunctionValue ek = Eigen::Matrix<double, N, 1>::Zero();
      int p = i / m_scalar_family.dimension();
      ek(p) = 1.;
      typename ScalarFamilyType::HessianValue hess = m_scalar_family.hessian(i % m_scalar_family.dimension(), x);
      return - hess.trace() * ek + hess.col(p);
    }

    /// Evaluate the curl curl of the i-th basis function at a quadrature point iqn, knowing all the values of the Hessian of ancestor basis functions at the quadrature nodes (provided by eval_quad)
    CurlValue curlcurl(size_t i, size_t iqn, const boost::multi_array<typename ScalarFamilyType::HessianValue, 2> &ancestor_hessian_quad) const
    {
      static_assert(hasCurlCurl, "Call to curlcurl() not available");

      FunctionValue ek = Eigen::Matrix<double, N, 1>::Zero();
      int p = i / m_scalar_family.dimension();
      ek(p) = 1.;
      typename ScalarFamilyType::HessianValue hess = ancestor_hessian_quad[i % m_scalar_family.dimension()][iqn];
      return - hess.trace() * ek + hess.col(p);
    }


    /// Return the ancestor (family that has been tensorized)
    constexpr inline const ScalarFamilyType &ancestor() const
    {
      return m_scalar_family;
    }

    /// Returns the maximum degree of the basis functions
    inline size_t max_degree() const
    {
      return m_scalar_family.max_degree();
    }

  private:
    ScalarFamilyType m_scalar_family;
  };

  //----------------------MATRIX FAMILY--------------------------------------------------------

  /// Matrix family obtained from a scalar family
  /** The transformation into matrix is done the following way: if \f$(f_1,\ldots,f_r)\f$ is the family of scalar functions,
      and \f$(E_1,\ldots,E_{N^2})\f$ is the canonical basis of the space of NxN matrices, the matrix family is 
      \f$(f_1E_1,f_2E_1,\ldots,f_rE_1,f_1E_2,\ldots...)\f$.
      The divergence is taken along each row, and \f$E_m\f$ has 1 at position (m%N, m/N)
  */
  // Useful formulas to manipulate this basis (all indices starting from 0):
  //      the i-th element in the basis is f_{i%r} E_{i/r}
  //      the matrix E_m has 1 at position (m%N, m/N)
  //      the element f_aE_b is the i-th in the family with i = ( [b/N]N + [b%N] )r + a
  template<typename ScalarFamilyType, size_t N>
  class MatrixFamily
  {
  public:
    typedef typename Eigen::Matrix<double, N, N> FunctionValue;
    typedef void GradientValue;
    typedef void CurlValue;
    typedef typename Eigen::Matrix<double, N, 1> DivergenceValue;
    typedef void HessianValue;    

    typedef typename ScalarFamilyType::GeometricSupport GeometricSupport;

    constexpr static const TensorRankE tensorRank = Matrix;
    constexpr static const bool hasAncestor = true;
    static const bool hasFunction = ScalarFamilyType::hasFunction;
    static const bool hasGradient = false;
    static const bool hasDivergence = ( ScalarFamilyType::hasGradient && N==dimspace );
    static const bool hasCurl = false;
    static const bool hasHessian = false;
    static const bool hasCurlCurl = false;
   
    typedef ScalarFamilyType AncestorType;

    MatrixFamily(const ScalarFamilyType & scalar_family)
      : m_scalar_family(scalar_family),
        m_E(N*N),
        m_transposeOperator(Eigen::MatrixXd::Zero(scalar_family.dimension()*N*N, scalar_family.dimension()*N*N))
    {
      static_assert(ScalarFamilyType::tensorRank == Scalar,
                    "Vector family can only be constructed from scalar families");
                    
      // Construct the basis for NxN matrices
      for (size_t j = 0; j < N; j++){
        for (size_t i = 0; i < N; i++){
          m_E[j*N + i] = Eigen::Matrix<double, N, N>::Zero();
          m_E[j*N + i](i,j) = 1.;
        }        
      }
      
      // Construct the transpose operator
      for (size_t i = 0; i < dimension(); i++){
        // The i-th basis element is sent to the j-th basis element s.t. j%r=i%r, (j/r)%N=(i/r)/N, (j/r)/N=(i/r)%N
        size_t r = m_scalar_family.dimension();
        size_t j = ( ( (i/r)%N )*N + (i/r)/N ) * r + (i%r);
        m_transposeOperator(i,j) = 1.;
      }
    }

    /// Return the dimension of the family
    inline size_t dimension() const
    {
      return m_scalar_family.dimension() * N * N;
    }

    /// Evaluate the i-th basis function at point x
    FunctionValue function(size_t i, const VectorRd & x) const
    {
      static_assert(hasFunction, "Call to function() not available");
            
      return m_scalar_family.function(i % m_scalar_family.dimension(), x) * m_E[i / m_scalar_family.dimension()];
    }
    
    /// Evaluate the i-th basis function at a quadrature point iqn, knowing all the values of ancestor basis functions at the quadrature nodes (provided by eval_quad)
    FunctionValue function(size_t i, size_t iqn, const boost::multi_array<double, 2> &ancestor_value_quad) const
    {
      static_assert(hasFunction, "Call to function() not available");

      return ancestor_value_quad[i % m_scalar_family.dimension()][iqn] * m_E[i / m_scalar_family.dimension()];
    }

    /// Evaluate the divergence of the i-th basis function at point x
    DivergenceValue divergence(size_t i, const VectorRd & x) const
    {
      static_assert(hasDivergence, "Call to divergence() not available");

      Eigen::Matrix<double, N, 1> V = m_E[i/m_scalar_family.dimension()].col( (i/m_scalar_family.dimension()) /N);
      return m_scalar_family.gradient(i % m_scalar_family.dimension(), x)( (i / m_scalar_family.dimension()) / N) * V;
    }
    
    /// Evaluate the divergence of the i-th basis function at a quadrature point iqn, knowing all the gradients of ancestor basis functions at the quadrature nodes (provided by eval_quad)
    DivergenceValue divergence(size_t i, size_t iqn, const boost::multi_array<VectorRd, 2> &ancestor_gradient_quad) const
    {
      static_assert(hasDivergence, "Call to divergence() not available");

      Eigen::Matrix<double, N, 1> V = m_E[i/m_scalar_family.dimension()].col( (i/m_scalar_family.dimension()) /N);
      return ancestor_gradient_quad[i % m_scalar_family.dimension()][iqn]( (i / m_scalar_family.dimension()) / N) * V;      
    }

    /// Return the ancestor (family that has been tensorized)
    inline const ScalarFamilyType &ancestor() const
    {
      return m_scalar_family;
    }
    
    /// Return the dimension of the matrices in the family
    inline const size_t matrixSize() const
    {
      return N;
    }
    
    /// Return the transpose operator, the rN^2 square matrix that to a given vector of coefficients on the Matrix family associate the vector of coefficients corresponding to the transpose
    inline const Eigen::MatrixXd transposeOperator() const
    {
      return m_transposeOperator;
    }
    
    /// Return the symmetrisation operator, the rN^2 square matrix that to a given vector of coefficients on the Matrix family associate the vector of coefficients corresponding to the symmetrised matrix
    inline const Eigen::MatrixXd symmetriseOperator() const
    {
      return (Eigen::MatrixXd::Identity(dimension(), dimension()) + m_transposeOperator)/2;
    }

    /// Returns the matrix that can be used, in a Family of this MatrixBasis, to create a basis of the subspace of symmetric matrices
    const Eigen::MatrixXd symmetricBasis() const
    {
      size_t dim_scalar =  m_scalar_family.dimension();
      Eigen::MatrixXd symOpe = symmetriseOperator();
      
      Eigen::MatrixXd SB = Eigen::MatrixXd::Zero(dim_scalar * N*(N+1)/2, dimension());
      
      // The matrix consists in selecting certain rows of symmetriseOperator, corresponding to the lower triangle part
      // To select these rows, we loop over the columns of the underlying matrices, and get the rows in symmetriseOperator()
      // corresponding to the lower triangular part of each column
      size_t position = 0;
      for (size_t i = 0; i<N; i++){
        // Skip the first i*dim_scalar elements in the column (upper triangle) in symOpe, take the remaining (N-i)*dim_scalar
        SB.middleRows(position, (N-i)*dim_scalar) = symOpe.middleRows(i*(N+1)*dim_scalar, (N-i)*dim_scalar);
        // update current position
        position += (N-i)*dim_scalar;
      }
      
      return SB;
    }

  private:
    ScalarFamilyType m_scalar_family;
    std::vector<Eigen::Matrix<double, N, N>> m_E;
    Eigen::MatrixXd m_transposeOperator;
  };


  //----------------------TANGENT FAMILY--------------------------------------------------------

  /// Vector family for polynomial functions that are tangent to a certain place (determined by the generators)
  template <typename ScalarFamilyType>
  class TangentFamily
  {
  public:
    typedef VectorRd FunctionValue;
    typedef VectorRd GradientValue;
    typedef VectorRd CurlValue;
    typedef double DivergenceValue;
    typedef Eigen::Matrix2d HessianValue;

    typedef Face GeometricSupport;

    constexpr static const TensorRankE tensorRank = Vector;
    constexpr static const bool hasAncestor = true;
    static const bool hasFunction = true;
    static const bool hasGradient = false;
    static const bool hasCurl = false;
    static const bool hasDivergence = ScalarFamilyType::hasGradient;
    static const bool hasHessian = false;
    static const bool hasCurlCurl = ScalarFamilyType::hasHessian;
    
    typedef ScalarFamilyType AncestorType;

    /// Constructor
    TangentFamily(
        const ScalarFamilyType &scalar_family,               ///< A basis for the scalar space
        const Eigen::Matrix<double, 2, dimspace> &generators ///< Two generators of the plane
        )
        : m_scalar_family(scalar_family),
          m_generators(generators),
          m_normal( ( generators.row(0).cross(generators.row(1)) ).normalized() )
    {
      static_assert(ScalarFamilyType::hasFunction, "Call to function() not available");
      static_assert(std::is_same<typename ScalarFamilyType::GeometricSupport, Face>::value,
                    "Tangent families can only be defined on faces");
    }

    /// Return the dimension of the family
    inline size_t dimension() const
    {
      return m_scalar_family.dimension() * 2;
    }

    /// Evaluate the i-th basis function at point x
    inline FunctionValue function(size_t i, const VectorRd &x) const
    {
      return m_generators.row(i / m_scalar_family.dimension()) * m_scalar_family.function(i % m_scalar_family.dimension(), x);
    }

    /// Evaluate the divergence of the i-th basis function at point x
    inline DivergenceValue divergence(size_t i, const VectorRd &x) const
    {
      return m_generators.row(i / m_scalar_family.dimension()).dot(m_scalar_family.gradient(i % m_scalar_family.dimension(), x));
    }
    
    /// Evaluate the curl curl of the i-th basis function at point x
    inline CurlValue curlcurl(size_t i, const VectorRd &x) const
    {
      return m_normal.cross( 
            m_scalar_family.hessian(i % m_scalar_family.dimension(), x) * 
              (m_normal.cross(m_generators.row(i / m_scalar_family.dimension())) ) );
    }

    /// Return the ancestor (family used for the tangent)
    constexpr inline const ScalarFamilyType &ancestor() const
    {
      return m_scalar_family;
    }

    /// Returns the maximum degree of the basis functions
    inline size_t max_degree() const
    {
      return m_scalar_family.max_degree();
    }

    /// Returns the generators of the basis functions
    inline Eigen::Matrix<double, 2, dimspace> generators() const
    {
      return m_generators;
    }


  private:
    ScalarFamilyType m_scalar_family;
    Eigen::Matrix<double, 2, dimspace> m_generators;
    VectorRd m_normal; // unit normal
  };

  //--------------------SHIFTED BASIS----------------------------------------------------------

  /// Generate a basis where the function indices are shifted. 
  /** Can be used, e.g., to ignore the constant function in a hierarchical scalar basis */
  template <typename BasisType>
  class ShiftedBasis
  {
  public:
    typedef typename BasisType::FunctionValue FunctionValue;
    typedef typename BasisType::GradientValue GradientValue;
    typedef typename BasisType::CurlValue CurlValue;
    typedef typename BasisType::DivergenceValue DivergenceValue;
    typedef typename BasisType::HessianValue HessianValue;

    typedef typename BasisType::GeometricSupport GeometricSupport;

    constexpr static const TensorRankE tensorRank = BasisType::tensorRank;
    constexpr static const bool hasAncestor = true;
    static const bool hasFunction = BasisType::hasFunction;
    static const bool hasGradient = BasisType::hasGradient;
    static const bool hasCurl = BasisType::hasCurl;
    static const bool hasDivergence = BasisType::hasDivergence;
    static const bool hasHessian = BasisType::hasHessian;
    static const bool hasCurlCurl = BasisType::hasCurlCurl;
    
    typedef BasisType AncestorType;

    /// Constructor
    ShiftedBasis(
        const BasisType &basis, ///< A basis
        const int shift         ///< The shift
        )
        : m_basis(basis),
          m_shift(shift)
    {
      // Do nothing
    }

    /// Return the dimension of the basis
    inline size_t dimension() const
    {
      return m_basis.dimension() - m_shift;
    }
    
    /// Return the underlying complete basis
    constexpr inline const BasisType &ancestor() const
    {
      return m_basis;
    }

    /// Returns the maximum degree of the basis functions
    inline size_t max_degree() const
    {
      return m_basis.max_degree();
    }
    
    /// Evaluate the i-th basis function at point x
    inline FunctionValue function(size_t i, const VectorRd &x) const
    {
      static_assert(hasFunction, "Call to function() not available");

      return m_basis.function(i + m_shift, x);
    }

    /// Evaluate the gradient of the i-th basis function at point x
    inline GradientValue gradient(size_t i, const VectorRd &x) const
    {
      static_assert(hasGradient, "Call to gradient() not available");

      return m_basis.gradient(i + m_shift, x);
    }

    /// Evaluate the curl of the i-th basis function at point x
    inline CurlValue curl(size_t i, const VectorRd &x) const
    {
      static_assert(hasCurl, "Call to curl() not available");

      return m_basis.curl(i + m_shift, x);
    }

    /// Evaluate the divergence of the i-th basis function at point x
    inline DivergenceValue divergence(size_t i, const VectorRd &x) const
    {
      static_assert(hasDivergence, "Call to divergence() not available");

      return m_basis.divergence(i + m_shift, x);
    }

    /// Evaluate the hessian of the i-th basis function at point x
    inline HessianValue hessian(size_t i, const VectorRd &x) const
    {
      static_assert(hasHessian, "Call to hessian() not available");

      return m_basis.hessian(i + m_shift, x);
    }

  private:
    BasisType m_basis;
    int m_shift;
  };

  //-----------------------RESTRICTED BASIS-------------------------------------------------------

  /// Generate a basis restricted to the first "dimension" functions.
  /** This can be useful, e.g., to form bases of subspaces of a given space from a hierarchical basis of the latter */
  template <typename BasisType>
  class RestrictedBasis
  {
  public:
    typedef typename BasisType::FunctionValue FunctionValue;
    typedef typename BasisType::GradientValue GradientValue;
    typedef typename BasisType::CurlValue CurlValue;
    typedef typename BasisType::DivergenceValue DivergenceValue;
    typedef typename BasisType::HessianValue HessianValue;

    typedef typename BasisType::GeometricSupport GeometricSupport;

    constexpr static const TensorRankE tensorRank = BasisType::tensorRank;
    constexpr static const bool hasAncestor = true;
    static const bool hasFunction = BasisType::hasFunction;
    static const bool hasGradient = BasisType::hasGradient;
    static const bool hasCurl = BasisType::hasCurl;
    static const bool hasDivergence = BasisType::hasDivergence;
    static const bool hasHessian = BasisType::hasHessian;
    static const bool hasCurlCurl = BasisType::hasCurlCurl;
    
    typedef BasisType AncestorType;

    /// Constructor
    RestrictedBasis(
        const BasisType &basis, ///< A basis
        const size_t &dimension ///< The dimension of the restricted basis
        )
        : m_basis(basis),
          m_dimension(dimension)
    {
      // Make sure that the provided dimension is smaller than the one
      // of the provided basis
      assert(dimension <= basis.dimension());
    }

    /// Return the dimension of the basis
    inline size_t dimension() const
    {
      return m_dimension;
    }
    
    /// Return the underlying complete basis
    constexpr inline const BasisType &ancestor() const
    {
      return m_basis;
    }
    
    /// Returns the maximum degree of the basis functions
    inline size_t max_degree() const
    {
      // We need to find the degree based on the dimension, assumes the basis is hierarchical
      size_t deg = 0;
      while (PolynomialSpaceDimension<GeometricSupport>::Poly(deg) < m_dimension){
        deg++;
      }
      return deg;
    }

    /// Evaluate the i-th basis function at point x
    inline FunctionValue function(size_t i, const VectorRd &x) const
    {
      static_assert(hasFunction, "Call to function() not available");

      return m_basis.function(i, x);
    }

    /// Evaluate the gradient of the i-th basis function at point x
    GradientValue gradient(size_t i, const VectorRd &x) const
    {
      static_assert(hasGradient, "Call to gradient() not available");

      return m_basis.gradient(i, x);
    }

    /// Evaluate the curl of the i-th basis function at point x
    CurlValue curl(size_t i, const VectorRd &x) const
    {
      static_assert(hasCurl, "Call to curl() not available");

      return m_basis.curl(i, x);
    }

    /// Evaluate the divergence of the i-th basis function at point x
    DivergenceValue divergence(size_t i, const VectorRd &x) const
    {
      static_assert(hasDivergence, "Call to divergence() not available");

      return m_basis.divergence(i, x);
    }

    /// Evaluate the hessian of the i-th basis function at point x
    HessianValue hessian(size_t i, const VectorRd &x) const
    {
      static_assert(hasHessian, "Call to hessian() not available");

      return m_basis.hessian(i, x);
    }

  private:
    BasisType m_basis;
    size_t m_dimension;
  };

  //--------------------GRADIENT BASIS----------------------------------------------------------

  /// Basis for the space of gradients of polynomials.
  /** To construct a basis of G^k, this assumes that the scalar basis it is constructed from is an ancestor basis of P^{k+1}/P^0, space of polynomials of degree k+1 without trivial polynomial with zero gradient (e.g. polynomials with zero average, or a hierarchical basis of P^{k+1} in which we have removed the first, constant, polynomial).
  This can also be used to create a family of gradients (not necessarily linearly independent, if the ancestor basis is not a basis of P^{k+1}/P^0) */
  template <typename BasisType>
  class GradientBasis
  {
  public:
    typedef typename BasisType::GradientValue FunctionValue;
    typedef void GradientValue;
    typedef void CurlValue;
    typedef void DivergenceValue;
    typedef void HessianValue;

    typedef typename BasisType::GeometricSupport GeometricSupport;

    constexpr static const TensorRankE tensorRank = (BasisType::tensorRank==Scalar ? Vector : Matrix);
    constexpr static const bool hasAncestor = true;
    static const bool hasFunction = true;
    static const bool hasGradient = false;
    static const bool hasCurl = false;
    static const bool hasDivergence = false;
    static const bool hasHessian = false;
    static const bool hasCurlCurl = false;
    
    typedef BasisType AncestorType;
    
    /// Constructor
    GradientBasis(const BasisType &basis)
        : m_scalar_basis(basis)
    {
      static_assert(BasisType::tensorRank == Scalar || BasisType::tensorRank == Vector,
                    "Gradient basis can only be constructed starting from scalar or vector bases");
      static_assert(BasisType::hasGradient,
                    "Gradient basis requires gradient() for the original basis to be available");
      // Do nothing
    }

    /// Compute the dimension of the basis
    inline size_t dimension() const
    {
      return m_scalar_basis.dimension();
    }

    /// Evaluate the i-th basis function at point x
    inline FunctionValue function(size_t i, const VectorRd &x) const
    {
      return m_scalar_basis.gradient(i, x);
    }

    /// Return the ancestor (basis that the gradient was taken of)
    constexpr inline const BasisType &ancestor() const
    {
      return m_scalar_basis;
    }


  private:
    BasisType m_scalar_basis;
  };

  //-------------------CURL BASIS-----------------------------------------------------------

  /// Basis for the space of curls of polynomials.
  /** To construct a basis of R^k, assumes that the vector basis from which it is constructed is a basis for G^{k+1,c} (in 3D) or P^{k+1}/P^0 (in 2D).
    This can also be used to create a family of curl (not necessarily linearly independent) */
  template <typename BasisType>
  class CurlBasis
  {
  public:
    typedef VectorRd FunctionValue;
    typedef Eigen::Matrix<double, dimspace, dimspace> GradientValue;
    typedef Eigen::Matrix<double, dimspace, dimspace> CurlValue;
    typedef double DivergenceValue;
    typedef Eigen::Matrix<double, dimspace, dimspace*dimspace> HessianValue;

    typedef typename BasisType::GeometricSupport GeometricSupport;

    constexpr static const TensorRankE tensorRank = Vector;
    constexpr static const bool hasAncestor = true;
    static const bool hasFunction = true;
    static const bool hasGradient = false;
    static const bool hasCurl = false;
    static const bool hasDivergence = false;
    static const bool hasHessian = false;
    static const bool hasCurlCurl = false;

    typedef BasisType AncestorType;

    /// Constructor
    CurlBasis(const BasisType &basis)
        : m_basis(basis)
    {
      static_assert((BasisType::tensorRank == Vector && std::is_same<typename BasisType::GeometricSupport, Cell>::value) ||
                        (BasisType::tensorRank == Scalar && std::is_same<typename BasisType::GeometricSupport, Face>::value),
                    "Curl basis can only be constructed starting from vector bases on elements or scalar bases on faces");
      static_assert(BasisType::hasCurl,
                    "Curl basis requires curl() for the original basis to be available");
    }

    /// Compute the dimension of the basis
    inline size_t dimension() const
    {
      return m_basis.dimension();
    }

    /// Evaluate the i-th basis function at point x
    inline FunctionValue function(size_t i, const VectorRd &x) const
    {
      return m_basis.curl(i, x);
    }
    
    /// Return the ancestor (basis that the gradient was taken of)
    constexpr inline const BasisType &ancestor() const
    {
      return m_basis;
    }


  private:
    size_t m_degree;
    BasisType m_basis;
  };
  
  
  //--------------------DIVERGENCE BASIS----------------------------------------------------------

  /// Basis (or rather family) of divergence of an existing basis
  /** This will be a real basis of the range of divergence (which is just P^k) if the ancestor basis is taken as a basis of R^{c,k+1} */
  template <typename BasisType>
  class DivergenceBasis
  {
  public:
    typedef typename BasisType::DivergenceValue FunctionValue;
    typedef void GradientValue;
    typedef void CurlValue;
    typedef void DivergenceValue;
    typedef void HessianValue;

    typedef typename BasisType::GeometricSupport GeometricSupport;

    constexpr static const TensorRankE tensorRank = Scalar;
    constexpr static const bool hasAncestor = true;
    static const bool hasFunction = true;
    static const bool hasGradient = false;
    static const bool hasCurl = false;
    static const bool hasDivergence = false;
    static const bool hasHessian = false;
    static const bool hasCurlCurl = false;

    typedef BasisType AncestorType;

    /// Constructor
    DivergenceBasis(const BasisType &basis)
        : m_vector_basis(basis)
    {
      static_assert(BasisType::tensorRank == Vector || BasisType::tensorRank == Matrix,
                    "Divergence basis can only be constructed starting from vector bases");
      static_assert(BasisType::hasDivergence,
                    "Divergence basis requires divergence() for the original basis to be available");
      // Do nothing
    }

    /// Compute the dimension of the basis
    inline size_t dimension() const
    {
      return m_vector_basis.dimension();
    }

    /// Evaluate the i-th basis function at point x
    inline FunctionValue function(size_t i, const VectorRd &x) const
    {
      return m_vector_basis.divergence(i, x);
    }

    /// Return the ancestor (basis that the gradient was taken of)
    constexpr inline const BasisType &ancestor() const
    {
      return m_vector_basis;
    }


  private:
    BasisType m_vector_basis;
  };


  /// The following function creates the "scalar rot" basis of a TangentFamily on a face
  /* It is created as the DivergenceBasis of a TangentFamily obtained rotating the generators by -pi/2 in the oriented face */
  template <typename ScalarFamilyType>
  DivergenceBasis<TangentFamily<ScalarFamilyType>> ScalarRotFamily(const TangentFamily<ScalarFamilyType> & tangent_family, const Face & F) 
  {
    // Take the generator and rotates them
    Eigen::MatrixXd gen = tangent_family.generators();
    VectorRd nF = F.normal();
    Eigen::MatrixXd rotated_gen = Eigen::MatrixXd::Zero(2, 3);
    rotated_gen.row(0) = VectorRd(gen.row(0)).cross(nF);
    rotated_gen.row(1) = VectorRd(gen.row(1)).cross(nF);
    
    TangentFamily<ScalarFamilyType> rotated_tangent_family(tangent_family.ancestor(), rotated_gen);

    return DivergenceBasis<TangentFamily<ScalarFamilyType>>(rotated_tangent_family);
  };

  //---------------------------------------------------------------------
  //      BASES FOR THE KOSZUL COMPLEMENTS OF G^k, R^k 
  //---------------------------------------------------------------------

  /// Basis for the complement R^{c,k}(T) in P^k(T)^3 of the range of curl
  class RolyComplBasisCell
  {
  public:
    typedef VectorRd FunctionValue;
    typedef Eigen::Matrix<double, dimspace, dimspace> GradientValue;
    typedef Eigen::Matrix<double, dimspace, dimspace> CurlValue;
    typedef double DivergenceValue;
    typedef Eigen::Matrix<double, dimspace, dimspace*dimspace> HessianValue;

    typedef Cell GeometricSupport;

    constexpr static const TensorRankE tensorRank = Vector;
    constexpr static const bool hasAncestor = false;
    static const bool hasFunction = true;
    static const bool hasGradient = false;
    static const bool hasCurl = false;
    static const bool hasDivergence = true;
    static const bool hasHessian = false;
    static const bool hasCurlCurl = false;

    /// Constructor
    RolyComplBasisCell(
        const Cell &T, ///< A mesh cell
        size_t degree  ///< The maximum polynomial degree to be considered
    );

    /// Compute the dimension of the basis
    inline size_t dimension() const
    {
      return PolynomialSpaceDimension<Cell>::RolyCompl(m_degree);
    }

    /// Evaluate the i-th basis function at point x
    FunctionValue function(size_t i, const VectorRd &x) const;

    /// Evaluate the divergence of the i-th basis function at point x
    DivergenceValue divergence(size_t i, const VectorRd &x) const;
    
    /// Returns the maximum degree of the basis functions
    inline size_t max_degree() const
    {
      return m_degree;
    }
    
    /// Returns the powers of the i-th basis function (not including the vector part)
    inline VectorZd powers(size_t i) const
    {
      return m_powers[i];
    }
    
  private:
    /// Coordinate transformation
    inline VectorRd _coordinate_transform(const VectorRd &x) const
    {
      return (x - m_xT) / m_hT;
    }

    size_t m_degree;
    VectorRd m_xT;
    double m_hT;
    std::vector<VectorZd> m_powers;
  };

  /// Basis for the complement G^{c,k}(T) in P^k(T)^3 of the range of grad.
  // This basis is obtained by translation and scaling of the following basis of x \times (P^{k-1}(R^3))^3 
  // (suggested by Daniel Mathews (daniel.mathews@monash.edu)):
  //     {scalar monomials in (x1,x2,x3) of degree <= k-1} * {(0, x3, -x2), (-x3, 0, x1)}
  //     and
  //     {scalar monomials in (x1,x2) of degree <=k-1} * (x2, -x1, 0) 
  // The vectors (0,x3,-x2), (-x3,0,x1) and (x2,-x1,0) are the "directions", the scalar factors in P^{k-1}
  // are given by m_powers
  class GolyComplBasisCell
  {
  public:
    typedef VectorRd FunctionValue;
    typedef Eigen::Matrix<double, dimspace, dimspace> GradientValue;
    typedef VectorRd CurlValue;
    typedef double DivergenceValue;
    typedef Eigen::Matrix<double, dimspace, dimspace*dimspace> HessianValue;

    typedef Cell GeometricSupport;

    constexpr static const TensorRankE tensorRank = Vector;
    constexpr static const bool hasAncestor = false;
    static const bool hasFunction = true;
    static const bool hasGradient = false;
    static const bool hasCurl = true;
    static const bool hasDivergence = false;
    static const bool hasHessian = false;
    static const bool hasCurlCurl = false;

    /// Constructor
    GolyComplBasisCell(
        const Cell &T, ///< A mesh cell
        size_t degree ///< The maximum polynomial degree to be considered
    );

    /// Compute the dimension of the basis
    inline size_t dimension() const
    {
      return PolynomialSpaceDimension<Cell>::GolyCompl(m_degree);;
    }

    /// Evaluate the i-th basis function at point x
    FunctionValue function(size_t i, const VectorRd &x) const;

    /// Evaluate the curl of the i-th basis function at point x
    CurlValue curl(size_t i, const VectorRd &x) const;
    
    /// Returns the maximum degree of the basis functions
    inline size_t max_degree() const
    {
      return m_degree;
    }
    
    /// Returns the powers of the i-th basis function (not including the vector part)
    inline VectorZd powers(size_t i) const
    {
      return m_powers[i];
    }
    
    /// Returns the dimension of P^{k-1}(R^3)
    inline size_t dimPkmo() const
    {
      return m_dimPkmo3D;
    }
    
  private:
    /// Coordinate transformation
    inline VectorRd _coordinate_transform(const VectorRd &x) const
    {
      return (x - m_xT) / m_hT;
    }

    // To evaluate the value and curl of the "direction"
    VectorRd direction_value(size_t i, const VectorRd &x) const;
    VectorRd direction_curl(size_t i, const VectorRd &x) const;

    size_t m_degree; // Degree k of the basis
    VectorRd m_xT; // center of cell
    double m_hT; // diameter of cell
    size_t m_dimPkmo3D; // shorthand for dimension of P^{k-1}(R^3)
    size_t m_dimPkmo2D; // shorthand for dimension of P^{k-1}(R^2)
    std::vector<VectorZd> m_powers; // vector listing the powers of the scalar factors in the basis

  };

  /// Basis for the complement R^{c,k}(F) in P^k(F)^2 of the range of the vectorial rotational on a face
  class RolyComplBasisFace
  {
  public:
    typedef VectorRd FunctionValue;
    typedef Eigen::Matrix<double, dimspace, dimspace> GradientValue;
    typedef Eigen::Matrix<double, dimspace, dimspace> CurlValue;
    typedef double DivergenceValue;
    typedef Eigen::Matrix<double, dimspace, dimspace*dimspace> HessianValue;

    typedef Face GeometricSupport;

    typedef Eigen::Matrix<double, 2, dimspace> JacobianType;

    constexpr static const TensorRankE tensorRank = Vector;
    constexpr static const bool hasAncestor = false;
    static const bool hasFunction = true;
    static const bool hasGradient = false;
    static const bool hasCurl = false;
    static const bool hasDivergence = true;
    static const bool hasHessian = false;
    static const bool hasCurlCurl = false;

    /// Constructor
    RolyComplBasisFace(
        const Face &F, ///< A mesh face
        size_t degree  ///< The maximum polynomial degree to be considered
    );

    /// Dimension of the basis
    inline size_t dimension() const
    {
      return PolynomialSpaceDimension<Face>::RolyCompl(m_degree);;
    }

    /// Evaluate the i-th basis function at point x
    FunctionValue function(size_t i, const VectorRd &x) const;

    /// Evaluate the divergence of the i-th basis function at point x
    DivergenceValue divergence(size_t i, const VectorRd &x) const;

    /// Return the Jacobian of the coordinate system transformation
    inline const JacobianType &jacobian() const
    {
      return m_jacobian;
    }

    /// Returns the maximum degree of the basis functions
    inline size_t max_degree() const
    {
      return m_degree;
    }
    
    /// Returns the powers of the i-th basis function
    inline Eigen::Vector2i powers(size_t i) const
    {
      return m_powers[i];
    }
    
  private:
    /// Coordinate transformation
    inline Eigen::Vector2d _coordinate_transform(const VectorRd &x) const
    {
      return m_jacobian * (x - m_xF);
    }

    size_t m_degree;
    VectorRd m_xF;
    double m_hF;
    JacobianType m_jacobian;
    std::vector<Eigen::Vector2i> m_powers;
  };


  /// Basis for the complement G^{c,k}(F) in P^k(F)^2 of the range of the gradient on a face.
  class GolyComplBasisFace
  {
  public:
    typedef VectorRd FunctionValue;
    typedef Eigen::Matrix<double, dimspace, dimspace> GradientValue;
    typedef Eigen::Matrix<double, dimspace, dimspace> CurlValue;
    typedef double DivergenceValue;
    typedef Eigen::Matrix<double, dimspace, dimspace*dimspace> HessianValue;

    typedef Face GeometricSupport;

    typedef Eigen::Matrix<double, 2, dimspace> JacobianType;

    constexpr static const TensorRankE tensorRank = Vector;
    constexpr static const bool hasAncestor = false;
    static const bool hasFunction = true;
    static const bool hasGradient = false;
    static const bool hasCurl = false;
    static const bool hasDivergence = false;
    static const bool hasHessian = false;
    static const bool hasCurlCurl = false;

    /// Constructor
    GolyComplBasisFace(
        const Face &F, ///< A mesh face
        size_t degree  ///< The maximum polynomial degree to be considered
    );

    /// Dimension of the basis
    inline size_t dimension() const
    {
      return PolynomialSpaceDimension<Face>::GolyCompl(m_degree);;
    }

    /// Evaluate the i-th basis function at point x
    FunctionValue function(size_t i, const VectorRd &x) const;

    /// Return the normal to the face used in the computation of the curl
    inline const VectorRd &normal() const
    {
      return m_nF;
    }

    /// Return the Rck basis
    inline const std::shared_ptr<RolyComplBasisFace> &rck() const
    {
      return m_Rck_basis;
    }

  private:
    size_t m_degree;
    VectorRd m_nF;
    std::shared_ptr<RolyComplBasisFace> m_Rck_basis;
  };

  //------------------------------------------------------------------------------
  //            Free functions
  //------------------------------------------------------------------------------

  enum BasisFunctionE
  {
    Function,
    Gradient,
    Curl,
    Divergence,
    Hessian,
    CurlCurl
  };
  
  /// Takes an array B_quad of values at quadrature nodes and applies the function F to all of them. F must take inValue and return outValue. The function must be called with outValue as template argument: transform_values_quad<outValue>(...)
  template <typename outValue, typename inValue, typename FunctionType>
  boost::multi_array<outValue, 2> transform_values_quad(
          const boost::multi_array<inValue, 2> &B_quad,   ///< Evaluations to be transformed
          const FunctionType &F                          ///< Transformation
          )
  {
    boost::multi_array<outValue, 2> transformed_B_quad( boost::extents[B_quad.shape()[0]][B_quad.shape()[1]] );
    
    std::transform( B_quad.origin(), B_quad.origin() + B_quad.num_elements(), transformed_B_quad.origin(), F);

    return transformed_B_quad;
  };

  /// From a scalar family B=(B_1..B_r) and vectors (v_1..v_k) in R^N, constructs a "Family" of "TensorizedVectorFamily" (built on B, of size N) that represents the family (B_1v_1..B_rv_1 B_1v_2...B_rv_2... B_1v_k..B_rv_k).
  /** Useful to tensorized scalar family while controlling the directions of tensorization (e.g. to identify tangential and normal directions along a surface, etc.). If B and v are orthonormal, then so is the returned family of tensorized basis. */
  template <typename ScalarBasisType, size_t N>  
  Family<TensorizedVectorFamily<ScalarBasisType, N>> GenericTensorization(
                                                          const ScalarBasisType & B,  /// The scalar family
                                                          const std::vector<Eigen::VectorXd> & v /// The vectors for tensorization
                                                          )
  {
    size_t r = B.dimension();
    size_t k = v.size();

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(r*k, r*N);
    
    for (size_t i=0; i < k; i++){
      // Check the dimension of vector v[i]
      assert(v[i].size() == N);

      // Fill in M
      for (size_t j=0; j < N; j++){
        M.block(i*r, j*r, r, r) = v[i](j) * Eigen::MatrixXd::Identity(r,r);
      }
    }
  
    TensorizedVectorFamily<ScalarBasisType, N> tensbasis(B);
    return Family<TensorizedVectorFamily<ScalarBasisType, N>>(tensbasis, M);
  };


  /// Function to symmetrise a matrix (useful together with transform_values_quad)
  inline static std::function<Eigen::MatrixXd(const Eigen::MatrixXd &)> 
  symmetrise_matrix = [](const Eigen::MatrixXd & x)->Eigen::MatrixXd { return 0.5*(x+x.transpose());};

  /// Function to skew-symmetrise a matrix (useful together with transform_values_quad)
  inline static std::function<Eigen::MatrixXd(const Eigen::MatrixXd &)> 
  skew_symmetrise_matrix = [](const Eigen::MatrixXd & x)->Eigen::MatrixXd { return 0.5*(x-x.transpose());};

  /// From a scalar family B, constructs a "Family" of "MatrixFamily" (built on B, of size NxN) that represents the family B Id on the MatrixFamily.
  /** Useful, e.g., to compute Gram matrices involving traces, by writing tr(G) B = G : B Id and computing the Gram matrix of G with the output of this function. */
  template <typename ScalarBasisType, size_t N>  
  Family<MatrixFamily<ScalarBasisType, N>> IsotropicMatrixFamily(
                                                          const ScalarBasisType & B  /// The scalar family
                                                          )
  {
    size_t r = B.dimension();
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(r, r*N*N);
    
    for (size_t k=0; k < N; k++){
      M.block(0, k*r*(N+1), r, r) = Eigen::MatrixXd::Identity(r, r);
    }
  
    MatrixFamily<ScalarBasisType, N> matbasis(B);
    return Family<MatrixFamily<ScalarBasisType, N>>(matbasis, M);
  };

  //------------------------------------------------------------------------------
  //                      BASIS EVALUATIONS
  //------------------------------------------------------------------------------

  //-----------------------DETAILS FOR EVALUATIONS----------------------

  namespace detail
  {
    /// Basis evaluation traits. Only specialization of 'BasisFunction' (=Function, Gradient, Curl or Divergence) are relevant, and determines what kind of value  we want to evaluate.
    /**  Provides a uniform way of evaluating the value, gradient, curl or divergence of functions in a basis. Specializations for TensorizedVectorFamily is also provided as it includes some additional information on the ancestor basis; this information is useful to optimise eval_quad for tensorized bases. */
    template <typename BasisType, BasisFunctionE BasisFunction>
    struct basis_evaluation_traits
    {
    };

    // Evaluate the function value at x
    template <typename BasisType>
    struct basis_evaluation_traits<BasisType, Function>
    {
      static_assert(BasisType::hasFunction, "Call to function not available");
      typedef typename BasisType::FunctionValue ReturnValue;
      static inline ReturnValue evaluate(const BasisType &basis, size_t i, const VectorRd &x)
      {
        return basis.function(i, x);
      }
      
      // Computes function value at quadrature node iqn, knowing values of ancestor basis at quadrature nodes
      static inline ReturnValue evaluate(
                    const BasisType &basis, 
                    size_t i, 
                    size_t iqn,
                    const boost::multi_array<ReturnValue, 2> &ancestor_value_quad
                    )
      {
        return basis.function(i, iqn, ancestor_value_quad);
      }
      
    };

    // Evaluate the gradient value at x
    template <typename BasisType>
    struct basis_evaluation_traits<BasisType, Gradient>
    {
      static_assert(BasisType::hasGradient, "Call to gradient not available");
      typedef typename BasisType::GradientValue ReturnValue;
      static inline ReturnValue evaluate(const BasisType &basis, size_t i, const VectorRd &x)
      {
        return basis.gradient(i, x);
      }
      
      // Computes gradient value at quadrature node iqn, knowing gradients of ancestor basis at quadrature nodes
      static inline ReturnValue evaluate(
                    const BasisType &basis, 
                    size_t i, 
                    size_t iqn,
                    const boost::multi_array<ReturnValue, 2> &ancestor_gradient_quad
                    )
      {
        return basis.gradient(i, iqn, ancestor_gradient_quad);
      }
    };

    // Evaluate the curl value at x
    template <typename BasisType>
    struct basis_evaluation_traits<BasisType, Curl>
    {
      static_assert(BasisType::hasCurl, "Call to curl not available");
      typedef typename BasisType::CurlValue ReturnValue;
      static inline ReturnValue evaluate(const BasisType &basis, size_t i, const VectorRd &x)
      {
        return basis.curl(i, x);
      }
      
      // Computes curl value at quadrature node iqn, knowing curls of ancestor basis at quadrature nodes
      static inline ReturnValue evaluate(
                    const BasisType &basis, 
                    size_t i, 
                    size_t iqn,
                    const boost::multi_array<ReturnValue, 2> &ancestor_curl_quad
                    )
      {
        return basis.curl(i, iqn, ancestor_curl_quad);
      }
    };

    // Evaluate the divergence value at x
    template <typename BasisType>
    struct basis_evaluation_traits<BasisType, Divergence>
    {
      static_assert(BasisType::hasDivergence, "Call to divergence not available");
      typedef typename BasisType::DivergenceValue ReturnValue;
      static inline ReturnValue evaluate(const BasisType &basis, size_t i, const VectorRd &x)
      {
        return basis.divergence(i, x);
      }
            
      // Computes divergence value at quadrature node iqn, knowing divergences of ancestor basis at quadrature nodes
      static inline ReturnValue evaluate(
                    const BasisType &basis, 
                    size_t i, 
                    size_t iqn,
                    const boost::multi_array<ReturnValue, 2> &ancestor_divergence_quad
                    )
      {
        return basis.divergence(i, iqn, ancestor_divergence_quad);
      }

    };
    
    // Evaluate the hessian value at x
    template <typename BasisType>
    struct basis_evaluation_traits<BasisType, Hessian>
    {
      static_assert(BasisType::hasHessian, "Call to hessian not available");
      typedef typename BasisType::HessianValue ReturnValue;
      static inline ReturnValue evaluate(const BasisType &basis, size_t i, const VectorRd &x)
      {
        return basis.hessian(i, x);
      }
            
      // Computes hessian value at quadrature node iqn, knowing hessian of ancestor basis at quadrature nodes
      static inline ReturnValue evaluate(
                    const BasisType &basis, 
                    size_t i, 
                    size_t iqn,
                    const boost::multi_array<ReturnValue, 2> &ancestor_hessian_quad
                    )
      {
        return basis.hessian(i, iqn, ancestor_hessian_quad);
      }

    };

    // Evaluate the curl curl value at x
    template <typename BasisType>
    struct basis_evaluation_traits<BasisType, CurlCurl>
    {
      static_assert(BasisType::hasCurlCurl, "Call to curl curl not available");
      typedef typename BasisType::CurlValue ReturnValue;
      static inline ReturnValue evaluate(const BasisType &basis, size_t i, const VectorRd &x)
      {
        return basis.curlcurl(i, x);
      }
            
      // Computes hessian value at quadrature node iqn, knowing hessian of ancestor basis at quadrature nodes
      static inline ReturnValue evaluate(
                    const BasisType &basis, 
                    size_t i, 
                    size_t iqn,
                    const boost::multi_array<ReturnValue, 2> &ancestor_curlcurl_quad
                    )
      {
        return basis.curlcurl(i, iqn, ancestor_curlcurl_quad);
      }

    };

    //---------- TENSORIZED FAMILY EVALUATION TRAITS -----------------------------------------
    // Evaluate the value at x of a TensorizedVectorFamily (includes information on the ancestor basis, to optimise eval_quad for tensorized bases)
    template <typename ScalarBasisType, size_t N>
    struct basis_evaluation_traits<TensorizedVectorFamily<ScalarBasisType, N>, Function>
    {
      static_assert(TensorizedVectorFamily<ScalarBasisType, N>::hasFunction, "Call to function not available");
      
      typedef typename TensorizedVectorFamily<ScalarBasisType, N>::FunctionValue ReturnValue;
      static const BasisFunctionE AncestorBasisFunction = Function; // Type of values needed from ancestor basis
      typedef double AncestorBasisFunctionValue; 
      
      // Computes function values at x
      static inline ReturnValue evaluate(
                    const TensorizedVectorFamily<ScalarBasisType, N> &basis, 
                    size_t i, 
                    const VectorRd &x
                    )
      {
        return basis.function(i, x);
      }
      
      // Computes function value at quadrature node iqn, knowing ancestor basis at quadrature nodes
      static inline ReturnValue evaluate(
                    const TensorizedVectorFamily<ScalarBasisType, N> &basis, 
                    size_t i, 
                    size_t iqn,
                    const boost::multi_array<double, 2> &ancestor_basis_quad
                    )
      {
        return basis.function(i, iqn, ancestor_basis_quad);
      }
    };

    // Evaluate the gradient at x of a TensorizedVectorFamily (includes information on the ancestor basis, to optimise eval_quad for tensorized bases)
    template <typename ScalarBasisType, size_t N>
    struct basis_evaluation_traits<TensorizedVectorFamily<ScalarBasisType, N>, Gradient>
    {
      static_assert(TensorizedVectorFamily<ScalarBasisType, N>::hasGradient, "Call to gradient not available");
      
      typedef typename TensorizedVectorFamily<ScalarBasisType, N>::GradientValue ReturnValue;
      static const BasisFunctionE AncestorBasisFunction = Gradient; // Type of values needed from ancestor basis
      typedef VectorRd AncestorBasisFunctionValue;
          
      // Computes gradient value at x
      static inline ReturnValue evaluate(
                    const TensorizedVectorFamily<ScalarBasisType, N> &basis, 
                    size_t i, 
                    const VectorRd &x
                    )
      {
        return basis.gradient(i, x);
      }

      // Computes gradient value at quadrature node iqn, knowing ancestor basis at quadrature nodes
      static inline ReturnValue evaluate(
                    const TensorizedVectorFamily<ScalarBasisType, N> &basis, 
                    size_t i, 
                    size_t iqn,
                    const boost::multi_array<VectorRd, 2> &ancestor_basis_quad
                    )
      {
        return basis.gradient(i, iqn, ancestor_basis_quad);
      }
    };

    // Evaluate the curl at x of a TensorizedVectorFamily (includes information on the ancestor basis, to optimise eval_quad for tensorized bases)
    template <typename ScalarBasisType, size_t N>
    struct basis_evaluation_traits<TensorizedVectorFamily<ScalarBasisType, N>, Curl>
    {
      static_assert(TensorizedVectorFamily<ScalarBasisType, N>::hasCurl, "Call to curl not available");
      
      typedef typename TensorizedVectorFamily<ScalarBasisType, N>::CurlValue ReturnValue;
      static const BasisFunctionE AncestorBasisFunction = Gradient; // Type of values needed from ancestor basis
      typedef VectorRd AncestorBasisFunctionValue;

      // Computes curl values at x
      static inline ReturnValue evaluate(
                    const TensorizedVectorFamily<ScalarBasisType, N> &basis, 
                    size_t i, 
                    const VectorRd &x
                    )
      {
        return basis.curl(i, x);
      }

      // Computes curl values at quadrature node iqn, knowing ancestor basis at quadrature nodes
      static inline ReturnValue evaluate(
                    const TensorizedVectorFamily<ScalarBasisType, N> &basis, 
                    size_t i, 
                    size_t iqn,
                    const boost::multi_array<VectorRd, 2> &ancestor_basis_quad
                    )
      {
        return basis.curl(i, iqn, ancestor_basis_quad);
      }
      
    };

    // Evaluate the divergence at x of a TensorizedVectorFamily (includes information on the ancestor basis, to optimise eval_quad for tensorized bases)
    template <typename ScalarBasisType, size_t N>
    struct basis_evaluation_traits<TensorizedVectorFamily<ScalarBasisType, N>, Divergence>
    {
      static_assert(TensorizedVectorFamily<ScalarBasisType, N>::hasDivergence, "Call to divergence not available");
      
      typedef typename TensorizedVectorFamily<ScalarBasisType, N>::DivergenceValue ReturnValue;
      static const BasisFunctionE AncestorBasisFunction = Gradient; // Type of values needed from ancestor basis
      typedef VectorRd AncestorBasisFunctionValue;

      // Computes divergence value at x
      static inline ReturnValue evaluate(
                    const TensorizedVectorFamily<ScalarBasisType, N> &basis, 
                    size_t i, 
                    const VectorRd &x
                    )
      {
        return basis.divergence(i, x);
      }

      // Computes divergence value at quadrature node iqn, knowing ancestor basis at quadrature nodes
      static inline ReturnValue evaluate(
                    const TensorizedVectorFamily<ScalarBasisType, N> &basis, 
                    size_t i, 
                    size_t iqn,
                    const boost::multi_array<VectorRd, 2> &ancestor_basis_quad
                    )
      {
        return basis.divergence(i, iqn, ancestor_basis_quad);
      }
      
    };

    // Evaluate the curl curl at x of a TensorizedVectorFamily (includes information on the ancestor basis, to optimise eval_quad for tensorized bases)
    template <typename ScalarBasisType, size_t N>
    struct basis_evaluation_traits<TensorizedVectorFamily<ScalarBasisType, N>, CurlCurl>
    {
      static_assert(TensorizedVectorFamily<ScalarBasisType, N>::hasCurlCurl, "Call to curl curl not available");
      
      typedef typename TensorizedVectorFamily<ScalarBasisType, N>::CurlValue ReturnValue;
      static const BasisFunctionE AncestorBasisFunction = Hessian; // Type of values needed from ancestor basis
      typedef MatrixRd AncestorBasisFunctionValue;

      // Computes curl curl value at x
      static inline ReturnValue evaluate(
                    const TensorizedVectorFamily<ScalarBasisType, N> &basis, 
                    size_t i, 
                    const VectorRd &x
                    )
      {
        return basis.curlcurl(i, x);
      }

      // Computes curl curl value at quadrature node iqn, knowing ancestor basis at quadrature nodes
      static inline ReturnValue evaluate(
                    const TensorizedVectorFamily<ScalarBasisType, N> &basis, 
                    size_t i, 
                    size_t iqn,
                    const boost::multi_array<MatrixRd, 2> &ancestor_basis_quad
                    )
      {
        return basis.curlcurl(i, iqn, ancestor_basis_quad);
      }
      
    };

    //---------- MATRIX FAMILY EVALUATION TRAITS -----------------------------------------
    // Evaluate the value at x of a MatrixFamily (includes information on the ancestor basis, to optimise eval_quad for matrix families)
    template <typename ScalarBasisType, size_t N>
    struct basis_evaluation_traits<MatrixFamily<ScalarBasisType, N>, Function>
    {
      static_assert(MatrixFamily<ScalarBasisType, N>::hasFunction, "Call to function not available");
      
      typedef typename MatrixFamily<ScalarBasisType, N>::FunctionValue ReturnValue;
      static const BasisFunctionE AncestorBasisFunction = Function; // Type of values needed from ancestor basis
      typedef double AncestorBasisFunctionValue; 
      
      // Computes function values at x
      static inline ReturnValue evaluate(
                                         const MatrixFamily<ScalarBasisType, N> &basis, 
                                         size_t i, 
                                         const VectorRd &x
                                         )
      {
        return basis.function(i, x);
      }
      
      // Computes function value at quadrature node iqn, knowing ancestor basis at quadrature nodes
      static inline ReturnValue evaluate(
                                         const MatrixFamily<ScalarBasisType, N> &basis, 
                                         size_t i, 
                                         size_t iqn,
                                         const boost::multi_array<double, 2> &ancestor_basis_quad
                                         )
      {
        return basis.function(i, iqn, ancestor_basis_quad);
      }
    };

    // Evaluate the divergence at x of a MatrixFamily (includes information on the ancestor basis, to optimise eval_quad for matrix families)
    template <typename ScalarBasisType, size_t N>
    struct basis_evaluation_traits<MatrixFamily<ScalarBasisType, N>, Divergence>
    {
      static_assert(MatrixFamily<ScalarBasisType, N>::hasDivergence, "Call to divergence not available");
      
      typedef typename MatrixFamily<ScalarBasisType, N>::DivergenceValue ReturnValue;
      static const BasisFunctionE AncestorBasisFunction = Gradient; // Type of values needed from ancestor basis
      typedef VectorRd AncestorBasisFunctionValue;

      // Computes divergence value at x
      static inline ReturnValue evaluate(
                                         const MatrixFamily<ScalarBasisType, N> &basis, 
                                         size_t i, 
                                         const VectorRd &x
                                         )
      {
        return basis.divergence(i, x);
      }

      // Computes divergence value at quadrature node iqn, knowing ancestor basis at quadrature nodes
      static inline ReturnValue evaluate(
                                         const MatrixFamily<ScalarBasisType, N> &basis, 
                                         size_t i, 
                                         size_t iqn,
                                         const boost::multi_array<VectorRd, 2> &ancestor_basis_quad
                                         )
      {
        return basis.divergence(i, iqn, ancestor_basis_quad);
      }
      
    };

  } // end of namespace detail

  
  //-----------------------------------EVALUATE_QUAD--------------------------------------

  /// Evaluate a basis at quadrature nodes. 'BasisFunction' (=Function, Gradient, Curl, Divergence or Hessian) determines what kind of value we want to evaluate.
  template <BasisFunctionE BasisFunction>
  struct evaluate_quad
  {
    /// Generic basis evaluation
    template <typename BasisType>
    static boost::multi_array<typename detail::basis_evaluation_traits<BasisType, BasisFunction>::ReturnValue, 2>
    compute(
        const BasisType &basis,    ///< The basis
        const QuadratureRule &quad ///< The quadrature rule
    )
    {
      typedef detail::basis_evaluation_traits<BasisType, BasisFunction> traits;

      boost::multi_array<typename traits::ReturnValue, 2>
          basis_quad(boost::extents[basis.dimension()][quad.size()]);

      for (size_t i = 0; i < basis.dimension(); i++)
      {
        for (size_t iqn = 0; iqn < quad.size(); iqn++)
        {
          basis_quad[i][iqn] = traits::evaluate(basis, i, quad[iqn].vector());
        } // for iqn
      }   // for i

      return basis_quad;
    }

    /// Evaluate a 'Family' of functions at quadrature nodes (optimised compared the generic basis evaluation, to avoid computing several times the ancestor basis at the quadrature nodes)
    template <typename BasisType>
    static boost::multi_array<typename detail::basis_evaluation_traits<Family<BasisType>, BasisFunction>::ReturnValue, 2>
    compute(
        const Family<BasisType> &basis, ///< The family
        const QuadratureRule &quad      ///< The quadrature rule
    )
    {
      typedef detail::basis_evaluation_traits<Family<BasisType>, BasisFunction> traits;

      boost::multi_array<typename traits::ReturnValue, 2>
          basis_quad(boost::extents[basis.dimension()][quad.size()]);

      // Compute values at quadrature note on ancestor basis, and then do a transformation
      const auto &ancestor_basis = basis.ancestor();
      boost::multi_array<typename traits::ReturnValue, 2>
          ancestor_basis_quad = evaluate_quad<BasisFunction>::compute(ancestor_basis, quad);

      for (size_t iqn = 0; iqn < quad.size(); iqn++)
      {
        for (size_t i = 0; i < size_t(basis.matrix().rows()); i++){
          basis_quad[i][iqn] = traits::evaluate(basis, i, iqn, ancestor_basis_quad);
        }
      }
      return basis_quad;
    }
    
    /// Evaluate a tensorized family at quadrature nodes (optimised compared the generic basis evaluation, to avoid computing several times the ancestor basis at the quadrature nodes)
    template <typename BasisType, size_t N>
    static boost::multi_array<typename detail::basis_evaluation_traits<TensorizedVectorFamily<BasisType, N>, BasisFunction>::ReturnValue, 2>
    compute(
        const TensorizedVectorFamily<BasisType, N> &basis, ///< The family
        const QuadratureRule &quad      ///< The quadrature rule
    )
    {
      typedef detail::basis_evaluation_traits<TensorizedVectorFamily<BasisType, N>, BasisFunction> traits;

      boost::multi_array<typename traits::ReturnValue, 2>
          basis_quad(boost::extents[basis.dimension()][quad.size()]);

      const auto &ancestor_basis = basis.ancestor();
      const boost::multi_array<typename traits::AncestorBasisFunctionValue, 2> ancestor_basis_quad 
                = evaluate_quad<traits::AncestorBasisFunction>::compute(ancestor_basis, quad);

      for (size_t i = 0; i < basis.dimension(); i++){
        for (size_t iqn = 0; iqn < quad.size(); iqn++){
          basis_quad[i][iqn] = traits::evaluate(basis, i, iqn, ancestor_basis_quad);
        }
      }
      return basis_quad;
    }

    /// Evaluate a Matrix family at quadrature nodes (optimised compared the generic basis evaluation, to avoid computing several times the ancestor basis at the quadrature nodes)
    template <typename BasisType, size_t N>
    static boost::multi_array<typename detail::basis_evaluation_traits<MatrixFamily<BasisType, N>, BasisFunction>::ReturnValue, 2>
    compute(
            const MatrixFamily<BasisType, N> &basis, ///< The family
            const QuadratureRule &quad      ///< The quadrature rule
            )
    {
      typedef detail::basis_evaluation_traits<MatrixFamily<BasisType, N>, BasisFunction> traits;

      boost::multi_array<typename traits::ReturnValue, 2>
        basis_quad(boost::extents[basis.dimension()][quad.size()]);

      const auto &ancestor_basis = basis.ancestor();
      const boost::multi_array<typename traits::AncestorBasisFunctionValue, 2> ancestor_basis_quad 
        = evaluate_quad<traits::AncestorBasisFunction>::compute(ancestor_basis, quad);

      for (size_t i = 0; i < basis.dimension(); i++){
        for (size_t iqn = 0; iqn < quad.size(); iqn++){
          basis_quad[i][iqn] = traits::evaluate(basis, i, iqn, ancestor_basis_quad);
        }
      }
      return basis_quad;
    }
    
  };


  //------------------------------------------------------------------------------
  //          ORTHONORMALISATION
  //------------------------------------------------------------------------------
  
  /// Gram-Schmidt algorithm to ortonormalize a basis.
  /// The matrix \f$M\f$ returned by this function gives the coefficients in the original basis of the
  /// orthonormalised basis. If \f$(f_1,...,f_r)\f$ is the original basis, the orthonormalised basis
  /// is \f$(\phi_1,...,\phi_r)\f$ where
  ///     \f$\phi_i = \sum_j M_{ij}f_j\f$.
  ///
  /// The function also modifies the variable basis_eval so that it contains the evaluation on at
  /// quadrature nodes of the new orthonormalised basis.
  template <typename T>
  Eigen::MatrixXd gram_schmidt(
      boost::multi_array<T, 2> &basis_eval,                      ///< Evaluations at quadrature nodes of the original basis.
      const std::function<double(size_t, size_t)> &inner_product ///< inner product (of two original basis functions) with respect to which we orthonormalise. This inner product must only depend on the basis functions through their values basis_eval.
  )
  {
    auto induced_norm = [inner_product](size_t i) -> double {
      return std::sqrt(inner_product(i, i));
    };

    // Number of basis functions
    size_t Nb = basis_eval.shape()[0];
    // Number of quadrature nodes
    size_t Nqn = basis_eval.shape()[1];

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(Nb, Nb);

    // Normalise the first element
    double norm = induced_norm(0);
    for (size_t iqn = 0; iqn < Nqn; iqn++)
    {
      basis_eval[0][iqn] /= norm;
    } // for iqn
    B(0, 0) = 1. / norm;

    for (size_t ib = 1; ib < Nb; ib++)
    {
      // 'coeffs' represents the coefficients of the ib-th orthogonal function on the ON basis functions 0 to ib-1
      Eigen::RowVectorXd coeffs = Eigen::RowVectorXd::Zero(ib);
      for (size_t pb = 0; pb < ib; pb++)
      {
        coeffs(pb) = -inner_product(ib, pb);
      } // for pb

      // store the values of the orthogonalised ib-th basis function
      for (size_t pb = 0; pb < ib; pb++)
      {
        for (size_t iqn = 0; iqn < Nqn; iqn++)
        {
          basis_eval[ib][iqn] += coeffs(pb) * basis_eval[pb][iqn];
        } // for iqn
      }   // for pb

      // normalise ib-th basis function
      double norm = induced_norm(ib);
      for (size_t iqn = 0; iqn < Nqn; iqn++)
      {
        basis_eval[ib][iqn] /= norm;
      } // for iqn
      coeffs /= norm;
      // Compute ib-th row of B.
      // B.topLeftCorner(ib, ib) contains the rows that represent the ON basis functions 0 to ib-1 on
      // the original basis. Multiplying on the left by coeffs gives the coefficients of the ON basis
      // functions on the original basis functions from 0 to ib-1.
      B.block(ib, 0, 1, ib) = coeffs * B.topLeftCorner(ib, ib);
      B(ib, ib) = 1. / norm;
    }
    return B;
  }

  //------------------------------------------------------------------------------

  /// Scalar product between two reals
  double scalar_product(const double &x, const double &y);

  /// Scalar product between one real and one 1-dimension Eigen vector
  double scalar_product(const double &x, const Eigen::Matrix<double, 1, 1> &y);

  /// Scalar product between two vectors
  double scalar_product(const VectorRd &x, const VectorRd &y);

  /// Scalar product between two matrices
  template<int N>
  double scalar_product(const Eigen::Matrix<double, N, N> & x, const Eigen::Matrix<double, N, N> & y)
  {
    return (x.transpose() * y).trace();
  }

  /// This overloading of the scalar_product function computes the scalar product between an evaluation of a basis and a constant value; both basis values and constant value must be of type Value
  template <typename Value>
  boost::multi_array<double, 2> scalar_product(
          const boost::multi_array<Value, 2> &basis_quad, ///< The basis evaluation
          const Value &v                                  ///< The vector to take the scalar product with
          )
  {
    boost::multi_array<double, 2> basis_dot_v_quad(boost::extents[basis_quad.shape()[0]][basis_quad.shape()[1]]);
    std::transform(basis_quad.origin(), basis_quad.origin() + basis_quad.num_elements(),
                   basis_dot_v_quad.origin(), [&v](const Value &x) -> double { return scalar_product(x,v); });
    return basis_dot_v_quad;
  }

  
  /// Compute the vector (cross) product between the evaluation of a basis and a constant vector
  boost::multi_array<VectorRd, 2>
  vector_product(
      const boost::multi_array<VectorRd, 2> &basis_quad, ///< The basis evaluation
      const VectorRd &v                                  ///< The vector to take the vector product with
  );

  /// \f$L^2\f$-orthonormalization: simply consists in using gram_schmidt() with the specific l2 inner product
  template <typename BasisType>
  Family<BasisType> l2_orthonormalize(
      const BasisType &basis,                                              ///< basis to orthonormalise
      const QuadratureRule &qr,                                            ///< quadrature rule for computing the l2 inner product
      boost::multi_array<typename BasisType::FunctionValue, 2> &basis_quad ///< values of basis functions at quadrature nodes
  )
  {
    // Check that the basis evaluation and quadrature rule are coherent
    assert(basis.dimension() == basis_quad.shape()[0] && qr.size() == basis_quad.shape()[1]);

    // The inner product between the i-th and j-th basis vectors
    std::function<double(size_t, size_t)> inner_product = [&basis_quad, &qr](size_t i, size_t j) -> double {
      double r = 0.;
      for (size_t iqn = 0; iqn < qr.size(); iqn++)
      {
        r += qr[iqn].w * scalar_product(basis_quad[i][iqn], basis_quad[j][iqn]);
      } // for iqn
      return r;
    };

    Eigen::MatrixXd B = gram_schmidt(basis_quad, inner_product);

    return Family<BasisType>(basis, B);
  }

  /// \f$L^2\f$-orthonormalization: when the Gram Matrix is passed, we use Cholesky. 
  /* This method is more robust that the previous one based on values at quadrature nodes */
  template <typename BasisType>
  Family<BasisType> l2_orthonormalize(
      const BasisType &basis,       ///< basis to orthonormalise
      const Eigen::MatrixXd & GM    ///< Gram matrix of the basis
  )
  {
    // Check that the basis and Gram matrix are coherent
    assert(basis.dimension() == size_t(GM.rows()) && GM.rows() == GM.cols());

    Eigen::MatrixXd L = GM.llt().matrixL();

    return Family<BasisType>(basis, L.inverse());
  }

  //------------------------------------------------------------------------------
  //      Gram matrices
  //------------------------------------------------------------------------------

  /// Compute the Gram-like matrix given the evaluation of two families of
  /// functions at quadrature nodes. This templated function is very generic, and thus not
  /// the most efficient. More efficient overloads are provided for double- or Vector3d-valued families
  template <typename FunctionValue>
  Eigen::MatrixXd compute_gram_matrix(const boost::multi_array<FunctionValue, 2> &B1, ///< First family at quadrature nodes
                                      const boost::multi_array<FunctionValue, 2> &B2, ///< Second family at quadrature nodes
                                      const QuadratureRule &qr,                       ///< Quadrature rule used for evaluation
                                      const size_t nrows,                             ///< Number of rows of the matrix (nb of members of first family to consider)
                                      const size_t ncols,                             ///< Number of rows of the matrix (nn of members of second family to consider)
                                      const std::string sym = "nonsym"                ///< Optional. "sym" to indicate that the matrix is symmetric (B1=B2)
  )
  {
    // Check that the basis evaluation and quadrature rule are coherent
    assert(qr.size() == B1.shape()[1] && qr.size() == B2.shape()[1]);
    // Check that we don't ask for more members of family than available
    assert(nrows <= B1.shape()[0] && ncols <= B2.shape()[0]);

    // Re-cast quadrature weights into ArrayXd to make computations faster
    Eigen::ArrayXd qr_weights = Eigen::ArrayXd::Zero(qr.size());
    for (size_t iqn = 0; iqn < qr.size(); iqn++)
    {
      qr_weights(iqn) = qr[iqn].w;
    }

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(nrows, ncols);
    for (size_t i = 0; i < nrows; i++)
    {
      size_t jcut = 0;
      if (sym == "sym")
        jcut = i;
      for (size_t j = 0; j < jcut; j++)
      {
        M(i, j) = M(j, i);
      }
      for (size_t j = jcut; j < ncols; j++)
      {
        std::vector<double> tmp(B1.shape()[1]);
        // Extract values at quadrature nodes for elements i of B1 and j of B2
        auto B1i = B1[boost::indices[i][boost::multi_array_types::index_range(0, B1.shape()[1])]];
        auto B2j = B2[boost::indices[j][boost::multi_array_types::index_range(0, B1.shape()[1])]];
        // Compute scalar product of B1i and B2j, and recast it as ArrayXd
        std::transform(B1i.begin(), B1i.end(), B2j.begin(), tmp.begin(), [](FunctionValue a, FunctionValue b) -> double { return scalar_product(a, b); });
        Eigen::ArrayXd tmp_array = Eigen::Map<Eigen::ArrayXd, Eigen::Unaligned>(tmp.data(), tmp.size());
        // Multiply by quadrature weights and sum (using .sum() of ArrayXd makes this step faster than a loop)
        M(i, j) = (qr_weights * tmp_array).sum();

        // Simple version with loop (replaces everything above after the for on j)
        /*
        for (size_t iqn = 0; iqn < qr.size(); iqn++) {
            M(i,j) += qr[iqn].w * scalar_product(B1[i][iqn], B2[j][iqn]);
        } // for iqn
  */

      } // for j
    }   // for i
    return M;
  }

  /// Compute the Gram-like matrix given the evaluation of two families of
  /// functions at quadrature nodes. This version calls the generic one with nrows = nb of elements in family B1
  /// and ncols = nb of elements in family B2.
  template <typename FunctionValue>
  Eigen::MatrixXd compute_gram_matrix(const boost::multi_array<FunctionValue, 2> &B1, ///< First family at quadrature nodes
                                      const boost::multi_array<FunctionValue, 2> &B2, ///< Second family at quadrature nodes
                                      const QuadratureRule &qr,                       ///< Quadrature rule used for evaluation
                                      const std::string sym = "nonsym"                ///< Optional. "sym" to indicate that the matrix is symmetric (B1=B2)
  )
  {
    return compute_gram_matrix<FunctionValue>(B1, B2, qr, B1.shape()[0], B2.shape()[0], sym);
  }

  /// Compute the Gram matrix given the evaluation of one family of functions at quadrature nodes. Consists in calling the generic templated version with B1=B2.
  template <typename FunctionValue>
  inline Eigen::MatrixXd compute_gram_matrix(const boost::multi_array<FunctionValue, 2> &B, ///< Family at quadrature nodes
                                             const QuadratureRule &qr                       ///< Quadrature rule used for evaluation
  )
  {
    return compute_gram_matrix<FunctionValue>(B, B, qr, "sym");
  }

  /// Compute the Gram-like matrix given a family of vector-valued and one of scalar-valued functions by tensorizing the latter
  Eigen::MatrixXd compute_gram_matrix(const boost::multi_array<VectorRd, 2> &B1, ///< First family at quadrature nodes
                                      const boost::multi_array<double, 2> &B2,   ///< Second family (to be tensorized) at quadrature nodes
                                      const QuadratureRule &qr                   ///< Quadrature rule used for evaluation
  );

  /// Compute the Gram-like matrix given the evaluation of two families of
  /// functions at quadrature nodes. This version is an overload for double-valued families, more efficient
  /// than the generic templated version.
  Eigen::MatrixXd compute_gram_matrix(const boost::multi_array<double, 2> &B1, ///< First family at quadrature nodes
                                      const boost::multi_array<double, 2> &B2, ///< Second family at quadrature nodes
                                      const QuadratureRule &qr,                ///< Quadrature rule used for evaluation
                                      const size_t nrows,                      ///< Number of rows of the matrix (nb of members of first family to consider)
                                      const size_t ncols,                      ///< Number of rows of the matrix (nn of members of second family to consider)
                                      const std::string sym = "nonsym"         ///< Optional. "sym" to indicate that the matrix is symmetric (B1=B2)
  );

  /// Compute the Gram-like matrix given the evaluation of two families of
  /// functions at quadrature nodes. Consists in calling the double-valued version with nrows = nb of elements in B1,
  /// ncols = nb of elements in B2
  Eigen::MatrixXd compute_gram_matrix(const boost::multi_array<double, 2> &B1, ///< First family at quadrature nodes
                                      const boost::multi_array<double, 2> &B2, ///< Second family at quadrature nodes
                                      const QuadratureRule &qr,                ///< Quadrature rule used for evaluation
                                      const std::string sym = "nonsym"         ///< Optional. "sym" to indicate that the matrix is symmetric (B1=B2)
  );

  /// Compute the Gram-like matrix given the evaluation of two families of
  /// functions at quadrature nodes. This version is an overload for Vector3d-valued families, more efficient
  /// than the generic templated version.
  Eigen::MatrixXd compute_gram_matrix(const boost::multi_array<VectorRd, 2> &B1, ///< First family at quadrature nodes
                                      const boost::multi_array<VectorRd, 2> &B2, ///< Second family at quadrature nodes
                                      const QuadratureRule &qr,                  ///< Quadrature rule used for evaluation
                                      const size_t nrows,                        ///< Optional. Number of rows of the matrix (nb of members of first family to consider).
                                      const size_t ncols,                        ///< Optional. Number of rows of the matrix (nb of members of second family to consider).
                                      const std::string sym = "nonsym"           ///< Optional. "sym" to indicate that the matrix is symmetric (B1=B2)
  );

  /// Compute the Gram-like matrix given the evaluation of two families of functions at quadrature nodes. Consists in calling the Vector3d-valued version with nrows = nb of elements in B1, ncols = nb of elements in B2
  Eigen::MatrixXd compute_gram_matrix(const boost::multi_array<VectorRd, 2> &B1, ///< First family at quadrature nodes
                                      const boost::multi_array<VectorRd, 2> &B2, ///< Second family at quadrature nodes
                                      const QuadratureRule &qr,                  ///< Quadrature rule used for evaluation
                                      const std::string sym = "nonsym"           ///< Optional. "sym" to indicate that the matrix is symmetric (B1=B2)
  );

  /// Compute the Gram-like matrix for a MatrixFamily. This overload is more efficient than the generic function as it only computes the gram matrix of the underlying scalar family, and then creates the bloc-diagonal gram matrix of the MatrixFamily (which is indeed bloc diagonal given the choice of m_E elements in this class).
  template<typename ScalarFamilyType, size_t N>
  Eigen::MatrixXd compute_gram_matrix(const MatrixFamily<ScalarFamilyType, N> & MatFam, ///< Matrix family
                                      const boost::multi_array<double, 2> & scalar_family_quad, ///< values of the underlying scalar family at quadrature nodes
                                      const QuadratureRule & qr        ///< Quadrature rule used for evaluation
                                      )
  {
    Eigen::MatrixXd Gram = Eigen::MatrixXd::Zero(MatFam.dimension(), MatFam.dimension());
        
    // Gram matrix of the scalar family
    Eigen::MatrixXd scalarGram = compute_gram_matrix<double>(scalar_family_quad, qr);
    for (size_t i = 0; i < MatFam.dimension(); i = i+scalarGram.rows()){
      Gram.block(i, i, scalarGram.rows(), scalarGram.cols()) = scalarGram;
    }   
    
    return Gram;
  }

  /// Computes the vector of integrals (f, phi_i)
  template <typename T>
  Eigen::VectorXd integrate(
      const FType<T> &f,        ///< Function to be integrated. T is any type that has a scalar_product
      const BasisQuad<T> &B,    ///< Family of basis functions at quadrature nodes. Possible types of T are double and VectorRd
      const QuadratureRule &qr, ///< Quadrature rule
      size_t n_rows = 0         ///< Optional argument for number of basis functions to be integrated. Default integrates all in the family.
  )
  {
    // If default, set n_rows to size of family
    if (n_rows == 0)
    {
      n_rows = B.shape()[0];
    }

    // Number of quadrature nodes
    const size_t num_quads = qr.size();

    // Check number of quadrature nodes is compatible with B
    assert(num_quads == B.shape()[1]);

    // Check that we don't ask for more members of family than available
    assert(n_rows <= B.shape()[0]);

    Eigen::VectorXd V = Eigen::VectorXd::Zero(n_rows);

    for (size_t iqn = 0; iqn < num_quads; iqn++)
    {
      double qr_weight = qr[iqn].w;
      T f_on_qr = f(qr[iqn].vector());
      for (size_t i = 0; i < n_rows; i++)
      {
        V(i) += qr_weight * scalar_product(B[i][iqn], f_on_qr);
      }
    }

    return V;
  }

  /// Computes the Gram-like matrix of integrals (f phi_i, phi_j)
  template <typename T, typename U>
  Eigen::MatrixXd compute_weighted_gram_matrix(
      const FType<U> &f,               ///< Weight function. Possible types of U are MatrixRd and double - must be compatible with T
      const BasisQuad<T> &B1,          ///< First family of basis functions at quadrature nodes. Possible types of T are VectorRd and double
      const BasisQuad<T> &B2,          ///< Second family of basis functions at quadrature nodes. Possible types of T are VectorRd and double
      const QuadratureRule &qr,        ///< Quadrature rule
      size_t n_rows = 0,               ///< Optional argument for number of functions from first family to be integrated against. Default integrates whole family
      size_t n_cols = 0,               ///< Optional argument for number of functions from second family to be integrated against. Default integrates whole family
      const std::string sym = "nonsym" ///< Optional argument if matrix is symmetric to increase efficiency
  )
  {
    // If default, set n_rows and n_cols to size of families
    if (n_rows == 0 && n_cols == 0)
    {
      n_rows = B1.shape()[0];
      n_cols = B2.shape()[0];
    }

    // Number of quadrature nodes
    const size_t num_quads = qr.size();
    // Check number of quadrature nodes is compatible with B1 and B2
    assert(num_quads == B1.shape()[1] && num_quads == B2.shape()[1]);
    // Check that we don't ask for more members of family than available
    assert(n_rows <= B1.shape()[0] && n_cols <= B2.shape()[0]);

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n_rows, n_cols);
    for (size_t iqn = 0; iqn < num_quads; iqn++)
    {
      double qr_weight = qr[iqn].w;
      U f_on_qr = f(qr[iqn].vector());
      for (size_t i = 0; i < n_rows; i++)
      {
        T f_B1 = f_on_qr * B1[i][iqn];
        size_t jcut = 0;
        if (sym == "sym")
          jcut = i;
        for (size_t j = 0; j < jcut; j++)
        {
          M(i, j) = M(j, i);
        }
        for (size_t j = jcut; j < n_cols; j++)
        {
          M(i, j) += qr_weight * scalar_product(f_B1, B2[j][iqn]);
        }
      }
    }
    return M;
  }

  /// Computes the Gram-like matrix of integrals (f phi_i, phi_j)
  template <typename T, typename U>
  Eigen::MatrixXd compute_weighted_gram_matrix(
      const FType<U> &f,        ///< Weight function. Possible types of U are MatrixRd and double - must be compatible with T
      const BasisQuad<T> &B1,   ///< First family of basis functions at quadrature nodes. Possible types of T are VectorRd and double
      const BasisQuad<T> &B2,   ///< Second family of basis functions at quadrature nodes. Possible types of T are VectorRd and double
      const QuadratureRule &qr, ///< Quadrature rule
      const std::string sym     ///< Argument if matrix is symmetric to increase efficiency
  )
  {
    return compute_weighted_gram_matrix(f, B1, B2, qr, B1.shape()[0], B2.shape()[0], sym);
  }

  /// Computes the Gram-like matrix of integrals (f dot phi_i, phi_j)
  Eigen::MatrixXd compute_weighted_gram_matrix(
      const FType<VectorRd> &f,      ///< Weight function
      const BasisQuad<VectorRd> &B1, ///< Family of vector basis functions at quadrature nodes
      const BasisQuad<double> &B2,   ///< Family of scalar basis functions at quadrature nodes
      const QuadratureRule &qr,      ///< Quadrature rule
      size_t n_rows = 0,             ///< Optional argument for number of functions from vector family to be integrated against. Default integrates whole family
      size_t n_cols = 0              ///< Optional argument for number of functions from scalar family to be integrated against. Default integrates whole family
  );

  /// Computes the Gram-like matrix of integrals (phi_i, f dot phi_j)
  Eigen::MatrixXd compute_weighted_gram_matrix(
      const FType<VectorRd> &f,      ///< Weight function
      const BasisQuad<double> &B1,   ///< Family of scalar basis functions at quadrature nodes
      const BasisQuad<VectorRd> &B2, ///< Family of vector basis functions at quadrature nodes
      const QuadratureRule &qr,      ///< Quadrature rule
      size_t n_rows = 0,             ///< Optional argument for number of functions from scalar family to be integrated against. Default integrates whole family
      size_t n_cols = 0              ///< Optional argument for number of functions from vector family to be integrated against. Default integrates whole family
  );

  //------------------------------------------------------------------------------
  //        L2 projection of a function
  //------------------------------------------------------------------------------

  /// Compute the L2-projection of a function
  template <typename BasisType>
  Eigen::VectorXd l2_projection(
      const std::function<typename BasisType::FunctionValue(const VectorRd &)> &f, ///< Function to project
      const BasisType &basis,                                                      ///< Basis for the space on which we project
      QuadratureRule &quad,                                                        ///< Quadrature rule
      const boost::multi_array<typename BasisType::FunctionValue, 2> &basis_quad,  ///< Evaluation of the basis at quadrature nodes
      const Eigen::MatrixXd &mass_basis = Eigen::MatrixXd::Zero(1,1)               ///< Optional: the mass matrix of the basis, if previously calculated
  )
  {
    Eigen::MatrixXd Mass = mass_basis;
    if (Mass.norm() < 1e-13){
      Mass = compute_gram_matrix(basis_quad, quad);
    }
    Eigen::LDLT<Eigen::MatrixXd> cholesky_mass(Mass);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(basis.dimension());
    for (size_t i = 0; i < basis.dimension(); i++)
    {
      for (size_t iqn = 0; iqn < quad.size(); iqn++)
      {
        b(i) += quad[iqn].w * scalar_product(f(quad[iqn].vector()), basis_quad[i][iqn]);
      } // for iqn
    }   // for i
    return cholesky_mass.solve(b);
  }

  //------------------------------------------------------------------------------
  //        Decomposition on a polynomial basis
  //------------------------------------------------------------------------------

  /// Structure to decompose a set of polynomials on a basis on a face or a cell
  /** The main interest of this structure for faces is to give a decomposition of traces (scalar, normal, etc.) of
  cell polynomials on a basis of face polynomials. This enables the usage of these traces in a face GramMatrix.
  However, the decomposition leads to increased rouding errors compared to a straight usage of compute_gram_matrix.
  Some tests however seem to indicate that this version of DecomposePoly is relatively robust up to degrees ~4 or 5
  */
  template <typename BasisType>
  struct DecomposePoly
  {
    /// Constructor for face
    /** The basis provided here must be of same rank (scalar, vector) as the polynomials to decompose, and of course
    of degree equal to or higher to those. The polynomials are decomposed on this basis.
    The method seems much more stable if we start from a straight monomial basis (possibly
    vectorialised, but not a family of some basis).
    */
    DecomposePoly(
            const Face &F,    ///< Face on which the polynomial functions are defined
            const BasisType &basis  ///< Basis for face polynomials
            ):
          m_dim(basis.dimension()),
          m_basis(basis),
          m_on_basis(basis, Eigen::MatrixXd::Identity(m_dim,m_dim))  // Need to initialise before we compute the real ON basis
      {
        /* The decomposition will be performed using an orthonormalised version of m_basis, for stability.
         The scalar product for which we orthonormalise is a minimal one (much less expensive than integrating
         over the face): \sum_i phi(x_i) psi(x_i) where (x_i)_i are nodes that are sufficient to determine
         entirely polynomials up to the considered degree.
         The nodes we build here correspond to P^k nodes on a triangle on a face (the triangle being centered
         at the face center) */
        
        // Create nodes
        std::vector<Eigen::Vector2i> indices = MonomialPowers<Face>::complete(basis.max_degree());
        VectorRd tF = F.edge(0)->tangent();
        VectorRd nF = F.edge_normal(0);
        VectorRd x0 = F.center_mass() - F.diam()*(tF+nF)/3.0;
        
        // Nodes 
        m_nb_nodes = indices.size();
        m_nodes.reserve(m_nb_nodes);
        for (size_t i=0; i<m_nb_nodes; i++){
          VectorRd xi = x0 + F.diam()*(double(indices[i](0)) * tF + double(indices[i](1)) * nF)/std::max(1.,double(basis.max_degree()));
          m_nodes.emplace_back(xi.x(), xi.y(), xi.z(), 1.); 
        }

        // We then orthonormalise Orthonormalise basis for simple scalar product
        m_on_basis_nodes.resize(boost::extents[m_dim][m_nb_nodes]);
        m_on_basis_nodes = evaluate_quad<Function>::compute(m_basis, m_nodes);
        m_on_basis = Family<BasisType>(l2_orthonormalize(m_basis, m_nodes, m_on_basis_nodes));
      };
      
    /// Constructor for cell
    DecomposePoly(
            const Cell &T,    ///< Cell on which the polynomial functions are defined
            const BasisType &basis  ///< Basis for face polynomials
            ):
          m_dim(basis.dimension()),
          m_basis(basis),
          m_on_basis(basis, Eigen::MatrixXd::Identity(m_dim,m_dim))  // Need to initialise before we compute the real ON basis
      {
        /* The decomposition will be performed using an orthonormalised version of m_basis, for stability.
         The scalar product for which we orthonormalise is a minimal one (much less expensive than integrating
         over the cell): \sum_i phi(x_i) psi(x_i) where (x_i)_i are nodes that are sufficient to determine
         entirely polynomials up to the considered degree.
         The nodes we build here correspond to P^k nodes on a tetrahedron */
        
        // Create nodes
        std::vector<Eigen::Vector3i> indices = MonomialPowers<Cell>::complete(basis.max_degree());
        VectorRd e0(1., 0., 0.);
        VectorRd e1(0., 1., 0.);
        VectorRd e2(0., 0., 1.);
        VectorRd xT = T.center_mass() - T.diam()*(e0+e1+e2)/3.0;
        
        // Nodes 
        m_nb_nodes = indices.size();
        m_nodes.reserve(m_nb_nodes);
        for (size_t i=0; i<m_nb_nodes; i++){
          VectorRd xi = xT + T.diam()*(double(indices[i](0)) * e0 + double(indices[i](1)) * e1 + double(indices[i](2)) * e2)/std::max(1.,double(basis.max_degree()));
          m_nodes.emplace_back(xi.x(), xi.y(), xi.z(), 1.); 
        }

        // We then orthonormalise Orthonormalise basis for simple scalar product
        m_on_basis_nodes.resize(boost::extents[m_dim][m_nb_nodes]);
        m_on_basis_nodes = evaluate_quad<Function>::compute(m_basis, m_nodes);
        m_on_basis = Family<BasisType>(l2_orthonormalize(m_basis, m_nodes, m_on_basis_nodes));
      };

    /// Return the set of nodes (useful to compute value of polynomial to decompose via evaluate_quad)
    inline QuadratureRule get_nodes() const 
    {
      return m_nodes;    
    };

    /// Returns the decomposed polynomials as a Family of the provided basis
    Family<BasisType> family(
            boost::multi_array<typename BasisType::FunctionValue, 2> &values  ///< Values of polynomials at the nodes (e.g. obtained via evaluate_quad)
            )
    {
      // Gram matrix of the polynomials and the ON basis
      Eigen::MatrixXd Gram_P_onbasis = compute_gram_matrix(values, m_on_basis_nodes, m_nodes);
      // L=matrix of P as a family of the ON basis. In theory, Gram_onbasis is identity, but due to rounding errors it 
      // can be a bit different. Computing L as if it's not the case improves stability 
      Eigen::MatrixXd Gram_onbasis = compute_gram_matrix(m_on_basis_nodes, m_nodes);
      Eigen::MatrixXd L = Gram_onbasis.ldlt().solve(Gram_P_onbasis.transpose());
      return Family<BasisType>(m_basis, L.transpose() * m_on_basis.matrix());
    };


    // Member variables
    size_t m_dim;   // dimension of basis
    BasisType m_basis; // Basis on which we decompose
    Family<BasisType> m_on_basis; // Orthonormalised basis (for simple scalar product)
    boost::multi_array<typename BasisType::FunctionValue, 2> m_on_basis_nodes; // Values of ON basis at the nodes
    size_t m_nb_nodes;    // nb of nodes
    QuadratureRule m_nodes; ///< Nodes for the interpolation
    
  };
  


  //@}

} // end of namespace HArDCore3D

#endif
