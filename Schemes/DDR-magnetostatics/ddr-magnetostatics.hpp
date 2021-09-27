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

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include <mesh.hpp>
#include <mesh_builder.hpp>

#include <xdiv.hpp>
#include <xcurl.hpp>

/*!
 * @defgroup DDR_magnetostatic
 * @brief Implementation of the DDR scheme for the magnetostatic problem
 */

namespace HArDCore3D
{

  /*!
   * @addtogroup DDR_magnetostatic
   * @{
   */

  /// Assemble a magnetostatic problem 
  struct Magnetostatics
  {
    typedef Eigen::SparseMatrix<double> SystemMatrixType;
    
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> ForcingTermType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> SolutionPotentialType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> SolutionCurlType;
    typedef IntegralWeight PermeabilityType;

    /// Constructor
    Magnetostatics(
                   const DDRCore & ddrcore,           ///< Core for the DDR space sequence
                   bool use_threads,                  ///< True for parallel execution, false for sequential execution
                   std::ostream & output = std::cout  ///< Output stream to print status messages
                   );

    /// Assemble the global system    
    void assembleLinearSystem(
                              const ForcingTermType & f,      ///< Forcing term
                              const PermeabilityType & mu,    ///< Permeability
                              const SolutionPotentialType & u ///< Boundary condition
                              );


    /// Returns the dimension of the magnetic field + potential space
    inline size_t dimensionSpace() const
    {
      return m_xcurl.dimension() + m_xdiv.dimension();
    }

    /// Returns the number of statically condensed DOFs (here, the cell magnetic field DOFs)
    inline size_t nbSCDOFs() const
    {
      return m_ddrcore.mesh().n_cells() * 
          (PolynomialSpaceDimension<Cell>::Roly(m_ddrcore.degree() - 1) + PolynomialSpaceDimension<Cell>::RolyCompl(m_ddrcore.degree()));
    }

    /// Returns the size of the statically condensed system
    inline size_t sizeSystem() const
    {
      return dimensionSpace() - nbSCDOFs();
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

    /// Returns the static condensation recovery operator
    inline const SystemMatrixType & scMatrix() const {
      return m_sc_A;
    }

    /// Returns the static condensation rhs
    inline Eigen::VectorXd & scVector() {
      return m_sc_b;
    }

    /// Compute the discrete Hcurl \times Hdiv norm of a family of vectors representing the dofs
    std::vector<double> computeNorms(
                       const std::vector<Eigen::VectorXd> & list_dofs ///< The list of vectors representing the dofs
                      ) const;

  private:
    /// Compute the local contribution for the element of index iT
    std::pair<Eigen::MatrixXd, Eigen::VectorXd>
    _compute_local_contribution(
                                size_t iT,                      ///< Element index
                                const ForcingTermType & f,      ///< Forcing term
                                const PermeabilityType & mu,    ///< Permeability
                                const SolutionPotentialType & u ///< Boundary condition
                                );

    /// Assemble the local contribution for the element of index iT into the global system
    void _assemble_local_contribution(
                                      size_t iT,                                               ///< Element index
                                      const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT, ///< Local contribution
                                      std::list<Eigen::Triplet<double> > & A1,                 ///< List of triplets for system
                                      Eigen::VectorXd & b1,                                    ///< Vector for RHS for sysem
                                      std::list<Eigen::Triplet<double> > & A2,                 ///< List of triplets for sc
                                      Eigen::VectorXd & b2                                     ///< Vector for RHS for sc
                                      );
    
    const DDRCore & m_ddrcore;
    bool m_use_threads;
    std::ostream & m_output;
    XCurl m_xcurl;
    XDiv m_xdiv;
    SystemMatrixType m_A;   // Matrix and RHS of statically condensed system
    Eigen::VectorXd m_b;
    SystemMatrixType m_sc_A; // Static condensation operator and RHS (to recover statically condensed DOFs)
    Eigen::VectorXd m_sc_b;
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

  static Magnetostatics::PermeabilityType
  constant_mu = Magnetostatics::PermeabilityType(1.);
  
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

  static Magnetostatics::PermeabilityType
  linear_mu = Magnetostatics::PermeabilityType(1.);
  
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

  static Magnetostatics::PermeabilityType
  trigonometric_mu = Magnetostatics::PermeabilityType(1.);

  //------------------------------------------------------------------------------
  

  static Magnetostatics::SolutionPotentialType
  variable_permeability_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                     cos(PI*x(0)) * sin(PI*x(1)) * sin(PI*x(2)),
                                                     -2. * sin(PI*x(0)) * cos(PI*x(1)) * sin(PI*x(2)),
                                                     sin(PI*x(0)) * sin(PI*x(1)) * cos(PI*x(2))
                                                     );
                            };

  static Magnetostatics::SolutionCurlType  
  variable_permeability_sigma = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                                  return 3 * PI / (1. + x(0) + x(1) + x(2))
                                    * Eigen::Vector3d(
                                                      sin(PI*x(0)) * cos(PI*x(1)) * cos(PI*x(2)),
                                                      0.,
                                                      -cos(PI*x(0)) * cos(PI*x(1)) * sin(PI*x(2))
                                                      );
                                };

  static Magnetostatics::ForcingTermType
  variable_permeability_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              double a = 3 * PI / std::pow(1. + x(0) + x(1) + x(2), 2);
                              double b = 3 * std::pow(PI, 2) / (1. + x(0) + x(1) + x(2));
                              double dy_sigmax = -a * sin(PI*x(0)) * cos(PI*x(1)) * cos(PI*x(2))
                                - b * sin(PI*x(0)) * sin(PI*x(1)) * cos(PI*x(2));
                              double dz_sigmax = -a * sin(PI*x(0)) * cos(PI*x(1)) * cos(PI*x(2))
                                - b * sin(PI*x(0)) * cos(PI*x(1)) * sin(PI*x(2));
                              ;
                              double dx_sigmaz = a * cos(PI*x(0)) * cos(PI*x(1)) * sin(PI*x(2))
                                + b * sin(PI*x(0)) * cos(PI*x(1)) * sin(PI*x(2));

                              double dy_sigmaz = a * cos(PI*x(0)) * cos(PI*x(1)) * sin(PI*x(2))
                                + b * cos(PI*x(0)) * sin(PI*x(1)) * sin(PI*x(2));

                              return Eigen::Vector3d(
                                                     dy_sigmaz,
                                                     dz_sigmax - dx_sigmaz,
                                                     -dy_sigmax
                                                     );
                            };

  static Magnetostatics::PermeabilityType
  variable_permeability_mu = Magnetostatics::PermeabilityType(
                        [](const Cell & T, const Eigen::Vector3d & x) -> double { return 1. + x(0) + x(1) + x(2); },
                        [](const Cell & T) -> size_t { return 1; }
                        );

} // end of namespace HArDCore3D

#endif
