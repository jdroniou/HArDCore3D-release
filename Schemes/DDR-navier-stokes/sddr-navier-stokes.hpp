// Implementation of the discrete de Rham sequence for the Navier-Stokes problem in curl-curl formulation.
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#ifndef SDDR_STOKES_HPP
#define SDDR_STOKES_HPP

#include <iostream>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include <mesh.hpp>
#include <mesh_builder.hpp>
#include <local_static_condensation.hpp>
#include <linearsolver.hpp>

#include <BoundaryConditions/BoundaryConditions.hpp>
#include <BoundaryConditions/BChandlers.hpp>   

#include <sxgrad.hpp>
#include <sxcurl.hpp>
#include <sxdiv.hpp>

/*!
 * @defgroup DDR_navier_stokes
 *  
 *  -# <i>Model</i>: solve the equations with the following boundary conditions:\n
 *    - natural: impose u.n and curl u x n on the boundary.\n
 *    - essential: impose the tangential component of u and p on the boundary\n
 * \n
 *  -# <i>Nomenclature</i>: we use the following terms. \n
 *    - DOF (degree of freedom): all the components of the velocity and pressure.\n
 *    - UKN (unknowns): only the non-Dirichlet components (components in the system before static condensation).\n
 *    - SCUKN: statically condensed (eliminated) unknowns.\n
 *    - GLUKN: globally coupled unknowns (non-statically condensed).
 *
 *  All the local calculations are made on all the DOFs (including the Lagrange multiplier). It's only when assembling the
 *  local contributions that we select the UKN (based on the map DOFtoUKN).
 * \n
 *  -# <i>Global vector</i>: the global vector we manipulate has the DOFs of u first, then the DOFs of p, then the Lagrange multiplier. 
 *    For each u,p we first put the vertices, then edge, then face, then cell DOFs. In each slice (vertex/edge/face) of the DOFs, the
 *    Dirichlet ones come first (we have re-numbered the geometric entities to ensure this).
 *
 * @brief Implementation of the serendipity DDR scheme for the %Navier-%Stokes problem in curl-curl formulation
 */

namespace HArDCore3D
{

  /*!
   * @addtogroup DDR_navier_stokes
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

    // Constructor without arguments
    StokesNorms() : StokesNorms(0., 0., 0., 0.) {};
    
    double u; ///< Norm of velocity
    double curl_u; ///< Norm of curl of velocity (vorticity)
    double p; ///< Norm of pressure
    double grad_p; ///< Norm of grad of pressure
    double hcurl_u; ///< Hcurl norm of u
    double hgrad_p; ///< Hgrad norm of p
  
  };

  /// Class to assemble and solve a Navier-Stokes problem 
  struct NavierStokes
  {
    typedef Eigen::SparseMatrix<double> SystemMatrixType;
    
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> ForcingTermType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> VelocityType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> VorticityType;
    typedef std::function<double(const Eigen::Vector3d &)> PressureType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> PressureGradientType;
    typedef IntegralWeight ViscosityType;

    /// Constructor
    NavierStokes(
           const DDRCore & ddrcore,           ///< Core for the DDR space sequence
           const BoundaryConditions & BC,     ///< Type of BC
           bool use_threads,                  ///< True for parallel execution, false for sequential execution
           std::ostream & output = std::cout  ///< Output stream to print status messages
         );

    //--------------------------------------------------
    //    Methods to assemble and solve
    //--------------------------------------------------    

    /// Create a vector of Dirichlet boundary conditions
    /** The vector has dimension velocity-pressure and zero everywhere except for Dirichlet DOFs on the velocity and pressure
    */
    Eigen::VectorXd vectorDirBC(
                                const VelocityType & u,    ///< BC on velocity
                                const PressureType & p     ///< BC on pressure
                               );

    /// Create the linear part of the system ( curl-curl velocity matrix, vel-grad pressure matrix, and RHS). Needs to be called before assembleLinearSystem
    void computeLinearComponents(
                                const ForcingTermType & f,      ///< Forcing term
                                const VelocityType & u,   ///< Boundary condition on velocity
                                const VorticityType & omega,  ///< Boundary condition on vorticity
                                const double & reynolds     ///< Reynolds number for BC on vorticity
                                );

    /// Assemble the global system at each Newton iteration, and returns the residual (norm of F(Xn)-b)
    double assembleLinearSystem(
                              const Eigen::VectorXd & Xn,  ///< State in the Newton iteration
                              const double & reac         ///< Coefficient of mass matrix for pseudo-temporal iterations
                              );

    /// Computes the norm of a vector (in Xcurl-Xgrad space) from the element-wise contributions to that vector, excluding the Dirichlet components (only the UKN are taken into account in the norm)
    double normGlobalRHS(
                        const std::vector<Eigen::VectorXd> & locRHS  ///< List of local contributions
                        );

    /// Solve the linear system and computes the Newton increment (taking into account relaxation)
    Eigen::VectorXd newtonIncrement( 
                                   LinearSolver<SystemMatrixType> & solver,           ///< Linear solver to use
                                   const double & navier_scaling,   ///< Scaling of the Navier term
                                   const double & relax,            ///< Relaxation parameter
                                   double & norm_dX                 ///< Norm of the increment dX (before applying relaxation)
                                   );

    /// Compute the discrete L2 norms, for a family of Eigen::VectorXd representing velocities & pressures, of u, curl u, p and grad p
    std::vector<StokesNorms> computeStokesNorms(
                const std::vector<Eigen::VectorXd> & list_dofs ///< List of vectors representing velocities and pressures
                ) const;

    /// Compute the continuous Hcurl errors in u and Hgrad error in p (1st component), and same with continuous norms (second component)
    std::pair<StokesNorms, StokesNorms> computeContinuousErrorsNorms(
                const Eigen::VectorXd & v,            ///< The vector, contains both velocity and pressure
                const VelocityType & u,               ///< Continuous velocity 
                const VorticityType & curl_u,         ///< Curl of continuous velocity
                const PressureType & p,               ///< Continuous pressure
                const PressureGradientType & grad_p   ///< Grad of continuous pressure
                ) const;

    /// Compute the flux (based on the polynomial reconstructions) for a discrete velocity across a flat convex hull
    double computeFlux(
                  const Eigen::VectorXd & u,            ///< Discrete vector (in Xcurl) of the velocity
                  const Hull & surface                  ///< Surface, described as a convex hull, across which the flux is computed
                  ) const;

    //------------------------------------------------------------------
    //    Getters: dimension of space and numbers of unknowns of various types
    //------------------------------------------------------------------
    
    /// Number of DOFs for velocity and pressure
    inline size_t nDOFs_up() const
    {
      return m_sxcurl.dimension() + m_sxgrad.dimension();
    }

    // ------ Dirichlet DOFs ------
    /// Number of DOFs on Dirichlet edges for the velocity
    inline size_t nDOFs_dir_edge_u() const
    {
      return m_nDOFs_dir_edge_u;
    }

    /// Number of DOFs on Dirichlet faces for the velocity
    inline size_t nDOFs_dir_face_u() const
    {
      return m_nDOFs_dir_face_u;
    }

    /// Total number of Dirichlet DOFs for the velocity
    inline size_t nDOFs_dir_u() const
    {
      return m_nDOFs_dir_edge_u + m_nDOFs_dir_face_u;
    }

    /// Number of DOFs on Dirichlet vertices for the pressure
    inline size_t nDOFs_dir_vertex_p() const
    {
      return m_nDOFs_dir_vertex_p;
    }

    /// Number of DOFs on Dirichlet edges for the pressure
    inline size_t nDOFs_dir_edge_p() const
    {
      return m_nDOFs_dir_edge_p;
    }

    /// Number of DOFs on Dirichlet faces for the pressure
    inline size_t nDOFs_dir_face_p() const
    {
      return m_nDOFs_dir_face_p;
    }

    /// Total number of Dirichlet DOFs for the pressure
    inline size_t nDOFs_dir_p() const
    {
      return m_nDOFs_dir_vertex_p + m_nDOFs_dir_edge_p + m_nDOFs_dir_face_p;
    }

    // ------ statically condensed and global UKN -----
    /// Number of statically condensed unknowns (velocity)
    inline size_t nSCUKN_u() const
    {
      return m_nloc_sc_u.sum();
    }

    /// Number of statically condensed unknowns (pressure)
    inline size_t nSCUKN_p() const
    {
      return m_nloc_sc_p.sum();
    }

    /// Number of statically condensed unknowns (both velocity and pressure)
    inline size_t nSCUKN() const
    {
      return nSCUKN_u() + nSCUKN_p();
    }
    
    /// Number of velocity unknowns (statically condensed and global, but no Dirichlet)
    inline size_t nUKN_u() const
    {
      return m_sxcurl.dimension() - nDOFs_dir_u();
    }
    
    /// Number of pressure unknowns (statically condensed and global)
    inline size_t nUKN_p() const
    {
      return m_sxgrad.dimension() - nDOFs_dir_p();
    }

    /// Number of Lagrange multiplier unknowns (1 if BC=N to impose \int p=0, 0 otherwise)
    inline size_t nUKN_lambda() const
    {
      return m_nUKN_lambda;
    }

    /// Number of unknowns (velocity-pressure-Lagrange multiplier, including statically condensed and globally coupled unknowns)
    inline size_t nUKN_uplambda() const
    {
      return nUKN_u() + nUKN_p() + nUKN_lambda();
    }

    /// Number of globally coupled unknowns for the velocity (statically condensed unknowns removed)
    inline size_t nGLUKN_u() const
    {
      return nUKN_u() - nSCUKN_u();
    }
    
    /// Number of globally coupled unknowns for the pressure (statically condensed unknowns removed)
    inline size_t nGLUKN_p() const
    {
      return nUKN_p() - nSCUKN_p();
    }

    /// Number of globally coupled unknowns for velocity, pressure, Lagrange multiplier
    inline size_t nGLUKN() const
    {
      return nUKN_uplambda() - nSCUKN();
    }

    //--------------------------------------------------------
    //  Getters: underlying SDDR spaces
    //--------------------------------------------------------
    
    /// Returns the space SXGrad
    inline const SXGrad & sxGrad() const
    {
      return m_sxgrad;
    }

    /// Returns the space SXCurl
    inline const SXCurl & sxCurl() const
    {
      return m_sxcurl;
    }

    /// Returns the space XDiv
    inline const SXDiv & sxDiv() const
    {
      return m_sxdiv;
    }

    //-------------------------------------------------
    //    Matrices and RHS
    //-------------------------------------------------
    
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

    /// Returns local RHS for Stokes problem (useful to compute relative residual)
    inline std::vector<Eigen::VectorXd> & rhsSourceNaturalBC() {
      return m_rhs_source_naturalBC;
    }

    //-------------------------------------------
    //  Getters: various parameters
    //-------------------------------------------
    
    /// Returns the stabilization parameter
    inline const double & stabilizationParameter() const {
      return m_stab_par;
    }

    /// Returns the stabilization parameter
    inline double & stabilizationParameter() {
      return m_stab_par;
    }

    /// Returns the selection of iterative solver
    inline const bool & iterativeSolver() const {
      return m_itsolver;
    }

    /// Returns the selection of iterative solver
    inline bool & iterativeSolver() {
      return m_itsolver;
    }

    /// Returns the reynolds number
    inline double & Reynolds() {
      return m_reynolds;
    }

  private:
    /// Compute the local linear contributions for the element of index iT
    /* First matrix: (CT ., CT .)_div
       Second matrix: (., GT .)_curl
       RHS: of size velocity + pressure + nb lagrange multiplier (0 or 1), and includes natural BCs
    */
    std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd>
    _compute_linear_contributions(
                                const size_t iT,                      ///< Element index
                                const Eigen::VectorXd & interp_f,     ///< Interpolate on xcurl of forcing term
                                const ForcingTermType & f             ///< Only required to test the non-pressure robust version
                                );

    /// From the local linear contributions, assemble the local Stokes matrix
    Eigen::MatrixXd _compute_Stokes_matrix(
                                        const size_t iT           ///< Element index
                                        );

    /// Compute the local contribution for the Jacobian matrix and the RHS in the Newton iterations
    std::pair<Eigen::MatrixXd, Eigen::VectorXd>
    _compute_jacobian_rhs(
                          const size_t iT,   ///< Element index
                          const Eigen::VectorXd & Xn,     ///< Current Newton state
                          const double & reac       ///< Coefficient of mass matrix for pseudo-temporal iterations
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
    
    // Constructor parameter
    const DDRCore & m_ddrcore;
    bool m_use_threads;
    std::ostream & m_output;
    SerendipityProblem m_ser_pro;
    const SXGrad m_sxgrad;
    const SXCurl m_sxcurl;
    const SXDiv m_sxdiv;
    const Eigen::VectorXi m_nloc_sc_u; // Nb of statically condensed DOFs for velocity in each cell (cell unknowns)
    const Eigen::VectorXi m_nloc_sc_p; // Nb of statically condensed DOFs for pressure in each cell (cell unknowns)

    // Boundary conditions
    BoundaryConditions m_BC;
    size_t m_nDOFs_dir_edge_u;
    size_t m_nDOFs_dir_face_u;
    size_t m_nDOFs_dir_vertex_p;
    size_t m_nDOFs_dir_edge_p;
    size_t m_nDOFs_dir_face_p;
    size_t m_nUKN_lambda;
    std::vector<std::pair<size_t,size_t>> m_locSCUKN;    // Location, among the DOFs, of the statically condensed unknowns (removing Dirichlet and non-SC unknowns); the unknowns are between m_locSCUCK[i].first and m_locSCUCK[i].first+m_locSCUKN[i].second for all i
    std::vector<std::pair<size_t,size_t>> m_locGLUKN;    // Location, among the DOFs, of the globally coupled unknowns (removing Dirichlet and non-global unknowns); the unknowns are between m_locGLUCK[i].first and m_locGLUCK[i].first+m_locGLUKN[i].second for all i
    Eigen::VectorXi m_DOFtoSCUKN;     // Maps a dof i to its SC unknown in the recovery operator, or to -1 if the DOF is not a SCUKN
    Eigen::VectorXi m_DOFtoGLUKN;     // Maps a dof i to its GL unknown in the global system, or to -1 if the DOF is a not a GLUKN

    // Scheme and solver parameters    
    double m_stab_par;
    bool m_itsolver;    // true if we use the iterative solver, false uses a direct solver

    double m_reynolds = 1.;   // Initialised at 1 by default, has to be adjusted after loading the value

    // Local matrices and RHS (including BCs) for the linear part of the model. The boolean to flag that these quantities have been computed
    std::vector<Eigen::MatrixXd> m_CTuCTu;    // curl - curl matrix
    std::vector<Eigen::MatrixXd> m_PTuGTp;    // velocity - grad pressure matrix
    std::vector<Eigen::VectorXd> m_rhs_source_naturalBC;   // RHS from the source and natural BCs
    std::vector<Eigen::MatrixXd> m_massT;     // mass-matrix for pseudo-temporal iterations
    bool m_linear_computed;
    // Global matrix and RHS of statically condensed system, at each Newton iteration
    SystemMatrixType m_A;   
    Eigen::VectorXd m_b;
    // Global static condensation operator and RHS (to recover statically condensed DOFs), at each Newton iteration
    SystemMatrixType m_sc_A; 
    Eigen::VectorXd m_sc_b;

  };

  //------------------------------------------------------------------------------
  // Exact solutions
  //------------------------------------------------------------------------------

  static const double PI = boost::math::constants::pi<double>();
  using std::sin;
  using std::cos;
  using std::pow;

  // To scale the pressure, the viscosity and the nonlinear term (curl u) x u
  double pressure_scaling = 1.;
  double reynolds = 1.;
  double navier_scaling = 1.;
  
  //------------------------------------------------------------------------------
  // Trigonometric solution
  //------------------------------------------------------------------------------

  static NavierStokes::VelocityType
  trigonometric_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                             0.5 * sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                             0.5 * cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                             -cos(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2))
                                             );
                    };

  static NavierStokes::VorticityType
  trigonometric_curl_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return 3. * PI * Eigen::Vector3d(
                                                      cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                                      -sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                                      0.
                                                      );
                         };

  static NavierStokes::PressureType
  trigonometric_p = [](const Eigen::Vector3d & x) -> double {
                      return pressure_scaling * sin(2. * PI * x(0)) * sin(2. * PI * x(1)) * sin(2. * PI * x(2));
                    };

  static NavierStokes::PressureGradientType
  trigonometric_grad_p = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return pressure_scaling * 2. * PI * Eigen::Vector3d(
                                                      cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                                      sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                                      sin(2. * PI * x(0)) * sin(2. * PI * x(1)) * cos(2. * PI * x(2))
                                                      );
                         };

  static NavierStokes::ForcingTermType
  trigonometric_curl_u_cross_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d (3.0*PI*sin(2*PI*x.x())*pow(sin(2*PI*x.z()), 2)*cos(2*PI*x.x())*pow(cos(2*PI*x.y()), 2),
                                              3.0*PI*sin(2*PI*x.y())*pow(sin(2*PI*x.z()), 2)*pow(cos(2*PI*x.x()), 2)*cos(2*PI*x.y()),
                                              1.5*PI*pow(sin(2*PI*x.x()), 2)*sin(2*PI*x.z())*pow(cos(2*PI*x.y()), 2)*cos(2*PI*x.z()) + 1.5*PI*pow(sin(2*PI*x.y()), 2)*sin(2*PI*x.z())*pow(cos(2*PI*x.x()), 2)*cos(2*PI*x.z())
                                              );
                    };

  static NavierStokes::ForcingTermType
  trigonometric_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return (1./reynolds) * 6. * std::pow(PI, 2) * Eigen::Vector3d(
                                                                    sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                                                    cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                                                    -2. * cos(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2))
                                                                    )
                        + navier_scaling * trigonometric_curl_u_cross_u(x)
                        + trigonometric_grad_p(x);
                    };

  //------------------------------------------------------------------------------
  // Constant solution
  //------------------------------------------------------------------------------

  static NavierStokes::VelocityType
  constant_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(1., 1., 1.);
             };

  static NavierStokes::VorticityType
  constant_curl_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d::Zero();
                  };

  static NavierStokes::ForcingTermType
  constant_curl_u_cross_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d::Zero();
                  };

  static NavierStokes::PressureType
  constant_p = [](const Eigen::Vector3d & x) -> double {
               return 0.;
             };

  static NavierStokes::PressureGradientType
  constant_grad_p = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d::Zero();
                  };

  static NavierStokes::ForcingTermType
  constant_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d::Zero();
             };

  //------------------------------------------------------------------------------
  // Linear solution
  //------------------------------------------------------------------------------

  static NavierStokes::VelocityType
  linear_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return 1e3*Eigen::Vector3d(x(0), -x(1), 0.);
             };

  static NavierStokes::VorticityType
  linear_curl_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d::Zero();
                  };

  static NavierStokes::ForcingTermType
  linear_curl_u_cross_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d::Zero();
                  };

  static NavierStokes::PressureType
  linear_p = [](const Eigen::Vector3d & x) -> double {
               return pressure_scaling * (x(0) - x(1));
             };

  static NavierStokes::PressureGradientType
  linear_grad_p = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return pressure_scaling * Eigen::Vector3d(1., -1., 0.);
                  };

  static NavierStokes::ForcingTermType
  linear_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return linear_grad_p(x);
             };

  //------------------------------------------------------------------------------
  // Vertical flow
  //------------------------------------------------------------------------------

  constexpr double scal_u = 1.;
  static NavierStokes::VelocityType
  vertical_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return scal_u * Eigen::Vector3d( cos(PI*x.y()), cos(PI*x.x()), 1.);
                    };

  static NavierStokes::VorticityType
  vertical_curl_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return scal_u * Eigen::Vector3d(0, 0, -PI*sin(PI*x.x()) + PI*sin(PI*x.y()) );
                         };

  static NavierStokes::PressureType
  vertical_p = trigonometric_p;

  static NavierStokes::PressureGradientType
  vertical_grad_p = trigonometric_grad_p;

  static NavierStokes::ForcingTermType
  vertical_curl_u_cross_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return std::pow(scal_u, 2) 
                           * Eigen::Vector3d (
                                             -(-PI*sin(PI*x.x()) + PI*sin(PI*x.y()))*cos(PI*x.x()),
                                             (-PI*sin(PI*x.x()) + PI*sin(PI*x.y()))*cos(PI*x.y()), 
                                             0.
                                             );
                    };

  static NavierStokes::ForcingTermType
  vertical_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return (1./reynolds) * scal_u * Eigen::Vector3d(pow(PI, 2)*cos(PI*x.y()), pow(PI, 2)*cos(PI*x.x()), 0.)
                        + navier_scaling * vertical_curl_u_cross_u(x)
                        + vertical_grad_p(x);
//                        return Eigen::Vector3d::Zero();
                    };


  //------------------------------------------------------------------------------
  // Impose BC: p=1 at entrance x=0, u.n=1 at exit x=1 on lower bottom corner
  //------------------------------------------------------------------------------
  static NavierStokes::VelocityType
  pressflux_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return ( lower_bottom_corner_x_one.is_in(x) ? VectorRd(1., 0., 0.) : VectorRd::Zero());
                    };

  static NavierStokes::VorticityType
  pressflux_curl_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                          return Eigen::Vector3d::Zero();
                        };

  static NavierStokes::PressureType
  pressflux_p = [](const Eigen::Vector3d & x) -> double {
                          // return (x.x()<1e-8 ? 1. : 0.);
                          return  pressure_scaling * (-x.z());
                        };

  static NavierStokes::PressureGradientType
  pressflux_grad_p = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                         return Eigen::Vector3d::Zero();
                      };

  static NavierStokes::ForcingTermType
  pressflux_curl_u_cross_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d::Zero();
                    };

  static NavierStokes::ForcingTermType
  pressflux_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d::Zero();
                    };


  /**@}*/

} // end of namespace HArDCore3D

#endif
