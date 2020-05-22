namespace HArDCore3D {
  
  namespace Tests {
    
    /// Compute the L2-orthogonal projection of a function
    template<typename BasisType>
    class L2Projection
    {
    public:
      typedef typename BasisType::FunctionValue FunctionValue;
      typedef boost::multi_array<FunctionValue, 2> BasisEvaluation;
      typedef std::function<typename BasisType::FunctionValue(const Eigen::Vector3d &)> FunctionType;
  
      /// Constructor
      L2Projection(
		   const BasisType & basis,           ///< Basis for the space on which we project
		   QuadratureRule & quad,             ///< Quadrature rule
		   const BasisEvaluation & basis_quad ///< Evaluation of the basis at quadrature nodes
		   )
	: m_basis(basis),
	  m_quad(quad),
	  m_basis_quad(basis_quad),
	  m_cholesky_mass(compute_gram_matrix(basis_quad, basis_quad, quad))
      {
	// Check that the basis evaluation and quadrature rule are coherent
	assert( basis.dimension() == basis_quad.shape()[0] && quad.size() == basis_quad.shape()[1] );
      }

      /// Compute the projection of a function
      Eigen::VectorXd compute(const FunctionType & f) const
      {
	Eigen::VectorXd b = Eigen::VectorXd::Zero(m_basis.dimension());
	for (size_t i = 0; i < m_basis.dimension(); i++) {
	  for (size_t iqn = 0; iqn < m_quad.size(); iqn++) {
	    FunctionValue f_iqn = f(m_quad[iqn].vector());
	    b(i) += m_quad[iqn].w * scalar_product(f_iqn, m_basis_quad[i][iqn]);
	  } // for iqn
	} // for i
	return m_cholesky_mass.solve(b);
      }

      /// Compute the squared L2-error between the L2-orthogonal projection
      /// and the function
      double squared_error(const FunctionType & f, const Eigen::VectorXd & f_dofs) const
      {
	double err = 0.;
	for (size_t iqn = 0; iqn < m_quad.size(); iqn++) {
	  FunctionValue f_iqn = f_dofs(0) * m_basis_quad[0][iqn];
	  for (size_t i = 1; i < m_basis.dimension(); i++) {
	    f_iqn += f_dofs(i) * m_basis_quad[i][iqn];
	  } // for i
	  FunctionValue diff_f_iqn = f_iqn - f(m_quad[iqn].vector());
	  err += m_quad[iqn].w * scalar_product(diff_f_iqn, diff_f_iqn);
	} // for iqn
	return err;
      }

      /// Compute the error between the L2-orthogonal projection and the function
      inline double error(const FunctionType & f, const Eigen::VectorXd & f_dofs) const
      {
	return std::sqrt( squared_error(f, f_dofs) );
      }

      /// Return the basis evaluation
      inline const BasisEvaluation & basisQuad()
      {
	return m_basis_quad;
      }
    
    private:
      BasisType m_basis;
      QuadratureRule m_quad;
      BasisEvaluation m_basis_quad;
      Eigen::LDLT<Eigen::MatrixXd> m_cholesky_mass;
    };

    /// Compute the squared L2-norm of a discrete function
    template<typename T>
    double squared_l2_norm(
			   const Eigen::VectorXd & f_dofs,             ///< Vector of DOFs of the discrete function
			   const QuadratureRule & quad,                ///< Quadrature rule
			   const boost::multi_array<T, 2> & basis_quad ///< Evaluation of the basis at quadrature nodes
			   )
    {
      double norm = 0.;
      for (size_t iqn = 0; iqn < quad.size(); iqn++) {
	T f_iqn = f_dofs(0) * basis_quad[0][iqn];
	for (size_t i = 1; i < basis_quad.shape()[0]; i++) {
	  f_iqn += f_dofs(i) * basis_quad[i][iqn];
	} // for i
	
	norm += quad[iqn].w * scalar_product(f_iqn, f_iqn);
      } // for iqn
      return norm;
    }
  
    /// Compute the L2-norm of a discrete function
    template<typename T>
    double l2_norm(
		   const Eigen::VectorXd & f_dofs,             ///< Vector of DOFs of the discrete function
		   const QuadratureRule & quad,                ///< Quadrature rule
		   const boost::multi_array<T, 2> & basis_quad ///< Evaluation of the basis at quadrature nodes
		   )
    {
      return std::sqrt( squared_l2_norm(f_dofs, quad, basis_quad) );
    }
    
  } // end of namespace Tests
  
} // end of namespace HArDCore3D
