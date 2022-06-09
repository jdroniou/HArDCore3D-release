// Implementation of the Virtual Element complex for the Stokes problem in curl-curl formulation.
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#ifndef STOKES_HPP
#define STOKES_HPP

#include <iostream>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include <mesh.hpp>
#include <mesh_builder.hpp>
#include <local_static_condensation.hpp>

#include <vdiv.hpp>
#include <vgrad.hpp>
#include <vcurl.hpp>

/*!
 * @defgroup VEM_stokes
 * @brief Implementation of the VEM scheme for the %Stokes problem in curl-curl formulation
 */

namespace HArDCore3D
{

  /*!
   * @addtogroup VEM_stokes
   * @{
   */
   
  /// Structure to store norm components (for velocity, its curl, pressure, and its gradient), as well as Hcurl and Hgrad norms of velocity and pressure
  struct StokesNorms
  {
    /// Constructor
    StokesNorms(double norm_u, double norm_curl_u, double norm_p, double norm_grad_p):
      u(norm_u),
      curl_u(norm_curl_u),
      p(norm_p),
      grad_p(norm_grad_p),
      hcurl_u( std::sqrt( std::pow(norm_u, 2) + std::pow(norm_curl_u, 2) ) ),
      hgrad_p( std::sqrt( std::pow(norm_p, 2) + std::pow(norm_grad_p, 2) ) )
      {
        // do nothing
      };
    
    double u; ///< Norm of velocity
    double curl_u; ///< Norm of curl of velocity (vorticity)
    double p; ///< Norm of pressure
    double grad_p; ///< Norm of grad of pressure
    double hcurl_u; ///< Hcurl norm of u
    double hgrad_p; ///< Hgrad norm of p
  
  };

  /// Assemble a Stokes problem 
  struct Stokes
  {
    typedef Eigen::SparseMatrix<double> SystemMatrixType;
    
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> ForcingTermType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> VelocityType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> VorticityType;
    typedef std::function<double(const Eigen::Vector3d &)> PressureType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> PressureGradientType;
    typedef IntegralWeight ViscosityType;

    /// Constructor
    Stokes(
           const VEMCore & vemcore,           ///< Core for the VEM space sequence
           bool use_threads,                  ///< True for parallel execution, false for sequential execution
           std::ostream & output = std::cout  ///< Output stream to print status messages
         );

    /// Assemble the global system    
    void assembleLinearSystem(
                              const ForcingTermType & f,      ///< Forcing term
                              const ForcingTermType & curl_f,      ///< Curl of forcing term
                              const VelocityType & u, ///< Boundary condition on velocity
                              const VorticityType & omega, ///< Boundary condition on vorticity
                              const ViscosityType & nu    ///< Viscosity
                              );

    /// Returns the global problem dimension (without Lagrange multiplier, just velocity & pressure)
    inline size_t dimension() const
    {
      return m_vcurl.dimension() + m_vgrad.dimension();
    }

    /// Returns the number of statically condensed DOFs (velocity)
    inline size_t nbSCDOFs_u() const
    {
      return m_vemcore.mesh().n_cells() * m_nloc_sc_u;
    }

    /// Returns the number of statically condensed DOFs (pressure)
    inline size_t nbSCDOFs_p() const
    {
      return m_vemcore.mesh().n_cells() * m_nloc_sc_p;
    }

    /// Returns the number of statically condensed DOFs (both velocity and pressure)
    inline size_t nbSCDOFs() const
    {
      return m_vemcore.mesh().n_cells() * (m_nloc_sc_u + m_nloc_sc_p);
    }

    /// Returns the size of the statically condensed system (with Lagrange multiplier)
    inline size_t sizeSystem() const
    {
      return dimension() + 1 - nbSCDOFs();
    }
    /// Returns the space VGrad
    inline const VGrad & vGrad() const
    {
      return m_vgrad;
    }

    /// Returns the space VCurl
    inline const VCurl & vCurl() const
    {
      return m_vcurl;
    }

    /// Returns the space VDiv
    inline const VDiv & vDiv() const
    {
      return m_vdiv;
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

    /// Returns the static condensation recovery operator
    inline const SystemMatrixType & scMatrix() const {
      return m_sc_A;
    }

    /// Returns the static condensation rhs
    inline Eigen::VectorXd & scVector() {
      return m_sc_b;
    }

    /// Returns the stabilization parameter
    inline const double & stabilizationParameter() const {
      return m_stab_par;
    }

    /// Returns the stabilization parameter
    inline double & stabilizationParameter() {
      return m_stab_par;
    }

    /// Compute the discrete L2 norms, for a family of Eigen::VectorXd representing velocities & pressures, of u, curl u, p and grad p
    std::vector<StokesNorms> computeStokesNorms(
                const std::vector<Eigen::VectorXd> & list_dofs ///< List of vectors representing velocities and pressures
                ) const;

    /// Compute the continuous Hcurl errors in u and Hgrad error in p (1st component), and same with continuous norms (second component)
    std::pair<StokesNorms, StokesNorms> computeContinuousErrorsNorms(
                const Eigen::VectorXd & v, ///< The vector, contains both velocity and pressure
                const VelocityType & u, ///< Continuous velocity 
                const VorticityType & curl_u, ///< Curl of continuous velocity
                const PressureType & p, ///< Continuous pressure
                const PressureGradientType & grad_p ///< Grad of continuous pressure
                ) const;


  private:
    /// Compute the local contribution for the element of index iT
    std::pair<Eigen::MatrixXd, Eigen::VectorXd>
    _compute_local_contribution(
                                size_t iT,                      ///< Element index
                                const Eigen::VectorXd & interp_f,  ///< Interpolate on vcurl of forcing term
                                const ViscosityType & nu    ///< Permeability
                                );

    /// Creates the permutation matrix and the global DOFs for the local static condensation
    LocalStaticCondensation _compute_static_condensation(const size_t & iT) const;

    /// Assemble the local contribution for the element of index iT into the global system
    void _assemble_local_contribution(
                                      size_t iT,                                               ///< Element index
                                      const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT, ///< Local contribution
                                      std::list<Eigen::Triplet<double> > & A1,                 ///< List of triplets for system
                                      Eigen::VectorXd & b1,                                    ///< Vector for RHS for sysem
                                      std::list<Eigen::Triplet<double> > & A2,                 ///< List of triplets for sc
                                      Eigen::VectorXd & b2                                     ///< Vector for RHS for sc
                                      );
    
    const VEMCore & m_vemcore;
    bool m_use_threads;
    std::ostream & m_output;
    VGrad m_vgrad;
    VCurl m_vcurl;
    VDiv m_vdiv;
    const size_t m_nloc_sc_u; // Nb of statically condensed DOFs for velocity in each cell (cell unknowns)
    const size_t m_nloc_sc_p; // Nb of statically condensed DOFs for pressure in each cell (cell unknowns)
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
  using std::pow;

  //------------------------------------------------------------------------------
  // Trigonometric solution
  //------------------------------------------------------------------------------

  double pressure_scaling = 1.;
  static Stokes::VelocityType
  trigonometric_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                             0.5 * sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                             0.5 * cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                             -cos(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2))
                                             );
                    };

  static Stokes::VorticityType
  trigonometric_curl_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return 3. * PI * Eigen::Vector3d(
                                                      cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                                      -sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                                      0.
                                                      );
                         };

  static Stokes::PressureType
  trigonometric_p = [](const Eigen::Vector3d & x) -> double {
                      return pressure_scaling * sin(2. * PI * x(0)) * sin(2. * PI * x(1)) * sin(2. * PI * x(2));
                    };

  static Stokes::PressureGradientType
  trigonometric_grad_p = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return pressure_scaling * 2. * PI * Eigen::Vector3d(
                                                      cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                                      sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                                      sin(2. * PI * x(0)) * sin(2. * PI * x(1)) * cos(2. * PI * x(2))
                                                      );
                         };

  static Stokes::ForcingTermType
  trigonometric_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return 6. * std::pow(PI, 2) * Eigen::Vector3d(
                                                                    sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                                                    cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                                                    -2. * cos(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2))
                                                                    )
                        + trigonometric_grad_p(x);
                    };

  static Stokes::ForcingTermType
  trigonometric_curl_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return 36. * std::pow(PI, 3) * Eigen::Vector3d(
                                                                    cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                                                    -sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                                                    0.
                                                                    );
                    };

  static Stokes::ViscosityType
  trigonometric_nu = Stokes::ViscosityType(1.);

  //------------------------------------------------------------------------------
  // Linear solution
  //------------------------------------------------------------------------------

  static Stokes::VelocityType
  linear_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(x(0), -x(1), 0.);
             };

  static Stokes::VorticityType
  linear_curl_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d::Zero();
                  };

  static Stokes::PressureType
  linear_p = [](const Eigen::Vector3d & x) -> double {
               return pressure_scaling * (x(0) - x(1));
             };

  static Stokes::PressureGradientType
  linear_grad_p = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return pressure_scaling * Eigen::Vector3d(1., -1., 0.);
                  };

  static Stokes::ForcingTermType
  linear_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return linear_grad_p(x);
             };

  static Stokes::ForcingTermType
  linear_curl_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d::Zero();
             };

  static Stokes::ViscosityType
  linear_nu = Stokes::ViscosityType(1.);


  //------------------------------------------------------------------------------
  // Solution derived from a field
  //------------------------------------------------------------------------------
  // Setting f(x,y,z)= (sin(pi x) sin(pi y) sin(pi z))^3, we take 
  //    U=(-d_y f, d_x f, 0), p=0
  // Calculations are made symbolically in matlab

  static Stokes::VelocityType
  field_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
               double ux = PI*cos(PI*x(1))*pow(sin(PI*x(0)),3.0)*pow(sin(PI*x(1)),2.0)*pow(sin(PI*x(2)),3.0)*3.0;
               double uy = PI*cos(PI*x(0))*pow(sin(PI*x(0)),2.0)*pow(sin(PI*x(1)),3.0)*pow(sin(PI*x(2)),3.0)*(-3.0);
               return Eigen::Vector3d(ux, uy, 0.);
             };

  static Stokes::VorticityType
  field_curl_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                double cux = (PI*PI)*cos(PI*x(0))*cos(PI*x(2))*pow(sin(PI*x(0)),2.0)*pow(sin(PI*x(1)),3.0)*pow(sin(PI*x(2)),2.0)*9.0;
                double cuy = (PI*PI)*cos(PI*x(1))*cos(PI*x(2))*pow(sin(PI*x(0)),3.0)*pow(sin(PI*x(1)),2.0)*pow(sin(PI*x(2)),2.0)*9.0;
                double cuz = (PI*PI)*pow(sin(PI*x(0)),3.0)*pow(sin(PI*x(1)),3.0)*pow(sin(PI*x(2)),3.0)*6.0-(PI*PI)*pow(cos(PI*x(0)),2.0)*sin(PI*x(0))*pow(sin(PI*x(1)),3.0)*pow(sin(PI*x(2)),3.0)*6.0-(PI*PI)*pow(cos(PI*x(1)),2.0)*pow(sin(PI*x(0)),3.0)*sin(PI*x(1))*pow(sin(PI*x(2)),3.0)*6.0;

                return Eigen::Vector3d(cux, cuy, cuz);
                };

  static Stokes::PressureType
  field_p = [](const Eigen::Vector3d & x) -> double {
               return 0.;
             };

  static Stokes::PressureGradientType
  field_grad_p = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(0., 0., 0.);
                  };

  static Stokes::ForcingTermType
  field_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
              double fx = (PI*PI*PI)*pow(cos(PI*x(1)),3.0)*pow(sin(PI*x(0)),3.0)*pow(sin(PI*x(2)),3.0)*-6.0+(PI*PI*PI)*cos(PI*x(1))*pow(sin(PI*x(0)),3.0)*pow(sin(PI*x(1)),2.0)*pow(sin(PI*x(2)),3.0)*3.9E+1-(PI*PI*PI)*pow(cos(PI*x(0)),2.0)*cos(PI*x(1))*sin(PI*x(0))*pow(sin(PI*x(1)),2.0)*pow(sin(PI*x(2)),3.0)*1.8E+1-(PI*PI*PI)*cos(PI*x(1))*pow(cos(PI*x(2)),2.0)*pow(sin(PI*x(0)),3.0)*pow(sin(PI*x(1)),2.0)*sin(PI*x(2))*1.8E+1;

              double fy = (PI*PI*PI)*pow(cos(PI*x(0)),3.0)*pow(sin(PI*x(1)),3.0)*pow(sin(PI*x(2)),3.0)*6.0-(PI*PI*PI)*cos(PI*x(0))*pow(sin(PI*x(0)),2.0)*pow(sin(PI*x(1)),3.0)*pow(sin(PI*x(2)),3.0)*3.9E+1+(PI*PI*PI)*cos(PI*x(0))*pow(cos(PI*x(1)),2.0)*pow(sin(PI*x(0)),2.0)*sin(PI*x(1))*pow(sin(PI*x(2)),3.0)*1.8E+1+(PI*PI*PI)*cos(PI*x(0))*pow(cos(PI*x(2)),2.0)*pow(sin(PI*x(0)),2.0)*pow(sin(PI*x(1)),3.0)*sin(PI*x(2))*1.8E+1;

              double fz = 0.;
               return Eigen::Vector3d(fx, fy, fz);
             };

  static Stokes::ForcingTermType
  field_curl_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
              double curl_fx = (PI*PI*PI*PI)*cos(PI*x(0))*pow(cos(PI*x(2)),3.0)*pow(sin(PI*x(0)),2.0)*pow(sin(PI*x(1)),3.0)*-1.8E+1-(PI*PI*PI*PI)*pow(cos(PI*x(0)),3.0)*cos(PI*x(2))*pow(sin(PI*x(1)),3.0)*pow(sin(PI*x(2)),2.0)*1.8E+1+(PI*PI*PI*PI)*cos(PI*x(0))*cos(PI*x(2))*pow(sin(PI*x(0)),2.0)*pow(sin(PI*x(1)),3.0)*pow(sin(PI*x(2)),2.0)*1.53E+2-(PI*PI*PI*PI)*cos(PI*x(0))*pow(cos(PI*x(1)),2.0)*cos(PI*x(2))*pow(sin(PI*x(0)),2.0)*sin(PI*x(1))*pow(sin(PI*x(2)),2.0)*5.4E+1;
;

              double curl_fy = (PI*PI*PI*PI)*cos(PI*x(1))*pow(cos(PI*x(2)),3.0)*pow(sin(PI*x(0)),3.0)*pow(sin(PI*x(1)),2.0)*-1.8E+1-(PI*PI*PI*PI)*pow(cos(PI*x(1)),3.0)*cos(PI*x(2))*pow(sin(PI*x(0)),3.0)*pow(sin(PI*x(2)),2.0)*1.8E+1+(PI*PI*PI*PI)*cos(PI*x(1))*cos(PI*x(2))*pow(sin(PI*x(0)),3.0)*pow(sin(PI*x(1)),2.0)*pow(sin(PI*x(2)),2.0)*1.53E+2-(PI*PI*PI*PI)*pow(cos(PI*x(0)),2.0)*cos(PI*x(1))*cos(PI*x(2))*sin(PI*x(0))*pow(sin(PI*x(1)),2.0)*pow(sin(PI*x(2)),2.0)*5.4E+1;
;

              double curl_fz = (PI*PI*PI*PI)*pow(sin(PI*x(0)),3.0)*pow(sin(PI*x(1)),3.0)*pow(sin(PI*x(2)),3.0)*7.8E+1-(PI*PI*PI*PI)*pow(cos(PI*x(0)),2.0)*sin(PI*x(0))*pow(sin(PI*x(1)),3.0)*pow(sin(PI*x(2)),3.0)*1.14E+2-(PI*PI*PI*PI)*pow(cos(PI*x(1)),2.0)*pow(sin(PI*x(0)),3.0)*sin(PI*x(1))*pow(sin(PI*x(2)),3.0)*1.14E+2-(PI*PI*PI*PI)*pow(cos(PI*x(2)),2.0)*pow(sin(PI*x(0)),3.0)*pow(sin(PI*x(1)),3.0)*sin(PI*x(2))*3.6E+1+(PI*PI*PI*PI)*pow(cos(PI*x(0)),2.0)*pow(cos(PI*x(1)),2.0)*sin(PI*x(0))*sin(PI*x(1))*pow(sin(PI*x(2)),3.0)*7.2E+1+(PI*PI*PI*PI)*pow(cos(PI*x(0)),2.0)*pow(cos(PI*x(2)),2.0)*sin(PI*x(0))*pow(sin(PI*x(1)),3.0)*sin(PI*x(2))*3.6E+1+(PI*PI*PI*PI)*pow(cos(PI*x(1)),2.0)*pow(cos(PI*x(2)),2.0)*pow(sin(PI*x(0)),3.0)*sin(PI*x(1))*sin(PI*x(2))*3.6E+1;
;
               return Eigen::Vector3d(curl_fx, curl_fy, curl_fz);
             };

  static Stokes::ViscosityType
  field_nu = Stokes::ViscosityType(1.);

  //------------------------------------------------------------------------------
  

} // end of namespace HArDCore3D

#endif
