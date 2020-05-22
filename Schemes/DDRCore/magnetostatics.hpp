// Implementation of the discrete de Rham sequence for the magnetostatic problem.
//
// Author: Daniele Di Pietro (daniele.di-pietro@umontpellier.fr)
//

/*
 * The DDR scheme is described in
 *
 *  An arbitrary-order method for magnetostatics on polyhedral meshes based on a discrete de Rham
sequence. 
 *   D. A. Di Pietro and J. Droniou, 31p, 2020. url: https://arxiv.org/abs/2005.06890.
 *
 * If you use this code in a scientific publication, please mention the above article.
 *
 */


#ifndef MAGNETOSTATICS_HPP
#define MAGNETOSTATICS_HPP

#include <iostream>

#include <boost/math/constants/constants.hpp>


/*!
 * @defgroup DDR_magnetostatic
 * @brief Implementation of the DDR scheme for the magnetostatic problem
 */

namespace HArDCore3D
{

  /*!
   *  @addtogroup DDR_magnetostatic
   * @{
   */

  /// Assemble a magnetostatic problem 
  struct Magnetostatics
  {
    typedef Eigen::SparseMatrix<double> SystemMatrixType;
    
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> ForcingTermType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> SolutionPotentialType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> SolutionCurlType;

    /// Constructor
    Magnetostatics(
		   const DDRCore & ddrcore,           ///< Core for the DDR space sequence
		   bool use_threads,                  ///< True for parallel execution, false for sequential execution
		   std::ostream & output = std::cout  ///< Output stream to print status messages
		   );

    /// Assemble the global system    
    void assembleLinearSystem(
			      const ForcingTermType & f,      ///< Forcing term
			      const SolutionPotentialType & u ///< Boundary condition
			      );

    /// Returns the global problem dimension
    inline size_t dimension() const
    {
      return m_xcurl.dimension() + m_xdiv.dimension();
    }

    /// Returns the space XCurl
    inline const XCurl & xCurl() const
    {
      return m_xcurl;
    }

    /// Returns the space XDiv
    inline const XDiv & xDiv() const
    {
      return m_xdiv;
    }

    /// Returns the linear system matrix
    inline const SystemMatrixType & systemMatrix() const {
      return m_A;
    }

    /// Returns the linear system matrix
    inline SystemMatrixType & systemMatrix() {
      return m_A;
    }

    /// Returns the linear system right-hand side vector
    inline const Eigen::VectorXd & systemVector() const {
      return m_b;
    }

    /// Returns the linear system right-hand side vector
    inline Eigen::VectorXd & systemVector() {
      return m_b;
    }

    /// Returns the stabilization parameter
    inline const double & stabilizationParameter() const {
      return m_stab_par;
    }

    /// Returns the stabilization parameter
    inline double & stabilizationParameter() {
      return m_stab_par;
    }


  private:
    /// Compute the local contribution for the element of index iT
    std::pair<Eigen::MatrixXd, Eigen::VectorXd>
    _compute_local_contribution(
				size_t iT,         ///< Element index
				const ForcingTermType & f, ///< Forcing term
				const SolutionPotentialType & u    ///< Boundary condition
				);

    /// Assemble the local contribution for the element of index iT into the global system
    void _assemble_local_contribution(
				      size_t iT,                                               ///< Element index
				      const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT, ///< Local contribution
				      std::list<Eigen::Triplet<double> > & A,                  ///< List of triplets
				      Eigen::VectorXd & b                                      ///< Vector
				      );

    const DDRCore & m_ddrcore;
    bool m_use_threads;
    std::ostream & m_output;
    XCurl m_xcurl;
    XDiv m_xdiv;
    SystemMatrixType m_A;
    Eigen::VectorXd m_b;
    double m_stab_par;
  };

  //------------------------------------------------------------------------------
  // Exact solutions
  //------------------------------------------------------------------------------

  static const double PI = boost::math::constants::pi<double>();
  using std::sin;
  using std::cos;

  //------------------------------------------------------------------------------
  
  static Magnetostatics::SolutionPotentialType
  constant_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
		 return Eigen::Vector3d(1., 2., 3.);
	       };

  static Magnetostatics::SolutionCurlType  
  constant_sigma = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
		     return Eigen::Vector3d::Zero();
		   };

  static Magnetostatics::ForcingTermType
  constant_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
		 return Eigen::Vector3d::Zero();
	       };

  //------------------------------------------------------------------------------

  static Magnetostatics::SolutionPotentialType
  linear_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
	       return Eigen::Vector3d(
				      x(0) + x(1) + x(2),
				      2.*(x(0)-x(1)+x(2)),
				      -x(0) - x(1) + x(2)
				      );
	     };

  static Magnetostatics::SolutionCurlType  
  linear_sigma = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
		   return Eigen::Vector3d(-3., 2., 1.);
		 };

  static Magnetostatics::ForcingTermType
  linear_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
	       return Eigen::Vector3d::Zero();
	     };
  
  //------------------------------------------------------------------------------

  static Magnetostatics::SolutionPotentialType
  trigonometric_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
		      return Eigen::Vector3d(
					     cos(PI*x(0)) * sin(PI*x(1)) * sin(PI*x(2)),
					     -2. * sin(PI*x(0)) * cos(PI*x(1)) * sin(PI*x(2)),
					     sin(PI*x(0)) * sin(PI*x(1)) * cos(PI*x(2))
					     );
		    };

  static Magnetostatics::SolutionCurlType  
  trigonometric_sigma = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
			  return 3. * PI * Eigen::Vector3d(
							   sin(PI*x(0)) * cos(PI*x(1)) * cos(PI*x(2)),
							   0.,
							   -cos(PI*x(0)) * cos(PI*x(1)) * sin(PI*x(2))
							   );
			};

  static Magnetostatics::ForcingTermType
  trigonometric_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
		      return 3. * std::pow(PI, 2) * Eigen::Vector3d(
								    cos(PI*x(0)) * sin(PI*x(1)) * sin(PI*x(2)),
								    -2. * sin(PI*x(0)) * cos(PI*x(1)) * sin(PI*x(2)),
								    sin(PI*x(0)) * sin(PI*x(1)) * cos(PI*x(2))
								    );
		    };

} // end of namespace HArDCore3D

#endif
