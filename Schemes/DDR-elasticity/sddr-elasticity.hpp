// Implementation of vsXgrad scheme for the elasticity equation
//
// Author: Jerome Droniou (jerome.droniou@umontpellier.fr)
//


#ifndef SDDR_ELASTICITY_HPP
#define SDDR_ELASTICITY_HPP

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

#include <vsxgrad.hpp>

/*!
 * @defgroup DDR_elasticity
 *  
 *  -# <i>Model</i>: solve the elasticity equation with the following boundary conditions:\n
 *    - essential: impose u on the boundary.\n
 *    - natural: impose eps(u)n on the boundary\n
 * \n
 *  -# <i>Nomenclature</i>: we use the following terms. \n
 *    - DOF (degree of freedom): all the components of the displacement.\n
 *    - UKN (unknowns): only the non-Dirichlet components (components in the system before static condensation).\n
 *    - SCUKN: statically condensed (eliminated) unknowns.\n
 *    - GLUKN: globally coupled unknowns (non-statically condensed).
 *
 *  All the local calculations are made on all the DOFs. It's only when assembling the
 *  local contributions that we select the UKN (based on the map DOFtoUKN).
 * \n
 *  -# <i>Global vector</i>: the global vector we manipulate has the DOFs of u. 
 *    We first put the vertices, then edge, then face, then cell DOFs. In each slice (vertex/edge/face) of the DOFs, the
 *    Dirichlet ones come first (we have re-numbered the geometric entities to ensure this).
 *
 * @brief Implementation of the serendipity Xgrad scheme for various elasticity problems in primal formulation, with or without fractures and with or without bubble DOFs.
 */

namespace HArDCore3D
{

  /*!
   * @addtogroup DDR_elasticity
   * @{
   */ 
   
  /// Class to assemble and solve a Navier-Stokes problem 
  struct Elasticity
  {
    typedef Eigen::SparseMatrix<double> SystemMatrixType;
    
    typedef std::function<Eigen::Vector3d(const double &, const double &, const Eigen::Vector3d &)> ForcingTermType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> DisplacementType;
    typedef std::function<Eigen::Matrix3d(const Eigen::Vector3d &)> GradientDisplacementType;

    /// Constructor
    Elasticity(
               const DDRCore & ddrcore,           ///< Core for the DDR space sequence
               const double & mu,                 ///< First Lame coefficient
               const double & lambda,             ///< Second Lame coefficient
               const BoundaryConditions & BC,     ///< Boundary conditions
               bool use_threads,                  ///< True for parallel execution, false for sequential execution
               std::ostream & output = std::cout  ///< Output stream to print status messages
               );

    //--------------------------------------------------
    //    Methods to assemble and solve
    //--------------------------------------------------    

    /// Create a vector of Dirichlet boundary conditions
    /** The vector has dimension vsXgrad and zero everywhere except for Dirichlet DOFs
    */
    Eigen::VectorXd vectorDirBC(
                                const DisplacementType & u    ///< BC on displacement
                               );

    /// Assemble the global system    
    void assembleLinearSystem(
                 const ForcingTermType & f,      ///< Forcing term
                 const GradientDisplacementType & symgrad_u, ///< Symmetric gradient of u for natural BCs
                 const Eigen::VectorXd & UDir        ///< Vector filled in by Dirichlets BCs (zero elsewhere) for essential BCs
                 );

    /// From dirichlet BCs and the solution of the globally coupled system, re-create the full vector of DOFs (with SC unknowns)
    inline Eigen::VectorXd reAssembleSolution (
                              const Eigen::VectorXd & UDir,    ///< Vector of Dirichlet BCs (of size nDOFs(), 0 outside Dirichlet DOFs)
                              const Eigen::VectorXd & u_gl     ///< Vector of globally coupled DOFs (of size nGLUKN()), solution of SC system
                              )
    {
      // Re-create statically condensed DOFs
      Eigen::VectorXd u_sc = scVector() + scMatrix() * u_gl;
      
      // Distribute u_gl and u_sc in global vector of DOFs, starting from UDir
      return replaceSectionsVector(
                               replaceSectionsVector(UDir, u_gl, m_locGLUKN), 
                               u_sc,
                               m_locSCUKN
                               );
    }
    


    /// Compute the discrete energy norm
//    std::vector<double> computeEnergyNorms(
//                const std::vector<Eigen::VectorXd> & list_dofs ///< List of vectors representing displacements
//                ) const;
                
    //------------------------------------------------------------------
    //    Getters: dimension of space and numbers of unknowns of various types
    //------------------------------------------------------------------
    
    /// Number of DOFs for displacement
    inline size_t nDOFs() const
    {
      return m_vsxgrad.dimension();
    }

    // ------ Dirichlet DOFs ------
    /// Number of DOFs on Dirichlet vertices
    inline size_t nDOFs_dir_vertices() const
    {
      return m_nDOFs_dir_vertices;
    }

    /// Number of DOFs on Dirichlet edges
    inline size_t nDOFs_dir_edges() const
    {
      return m_nDOFs_dir_edges;
    }

    /// Number of DOFs on Dirichlet faces
    inline size_t nDOFs_dir_faces() const
    {
      return m_nDOFs_dir_faces;
    }

    /// Total number of Dirichlet DOFs
    inline size_t nDOFs_dir() const
    {
      return m_nDOFs_dir_vertices + m_nDOFs_dir_edges + m_nDOFs_dir_faces;
    }

    // ------ statically condensed and global UKN -----
    /// Number of statically condensed unknowns
    inline size_t nSCUKN() const
    {
      return m_nloc_sc.sum();
    }
    
    /// Number of displacement unknowns (statically condensed and global, but no Dirichlet)
    inline size_t nUKN() const
    {
      return m_vsxgrad.dimension() - nDOFs_dir();
    }
    
    /// Number of globally coupled unknowns for the displacement (statically condensed unknowns removed)
    inline size_t nGLUKN() const
    {
      return nUKN() - nSCUKN();
    }
    
        
    //--------------------------------------------------------
    //  Getters: underlying SDDR spaces
    //--------------------------------------------------------
    
    /// Returns the space vsxgrad
    inline const VSXGrad & vsxGrad() const
    {
      return m_vsxgrad;
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

  private:
    /// Compute the local contribution for the element of index iT
    std::pair<Eigen::MatrixXd, Eigen::VectorXd>
    _compute_local_contribution (
                      const size_t iT,            ///< Element index
                      const ForcingTermType & f,  ///< Forcing term
                      const GradientDisplacementType & symgrad_u, ///< Symmetric gradient of u for natural BCs
                      const Eigen::VectorXd & UDir ///< Boundary condition
                      ) const ;

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
    const VSXGrad m_vsxgrad;
    const Eigen::VectorXi m_nloc_sc; // Nb of statically condensed DOFs for velocity in each cell (cell unknowns)
    
    // Physical parameters
    double m_mu;
    double m_lambda;

    // Boundary conditions
    BoundaryConditions m_BC;
    size_t m_nDOFs_dir_vertices;
    size_t m_nDOFs_dir_edges;
    size_t m_nDOFs_dir_faces;
    std::vector<std::pair<size_t,size_t>> m_locSCUKN;    // Location, among the DOFs, of the statically condensed unknowns (removing Dirichlet and non-SC unknowns); the unknowns are between m_locSCUCK[i].first and m_locSCUCK[i].first+m_locSCUKN[i].second for all i
    std::vector<std::pair<size_t,size_t>> m_locGLUKN;    // Location, among the DOFs, of the globally coupled unknowns (removing Dirichlet and non-global unknowns); the unknowns are between m_locGLUCK[i].first and m_locGLUCK[i].first+m_locGLUKN[i].second for all i
    Eigen::VectorXi m_DOFtoSCUKN;     // Maps a dof i to its SC unknown in the recovery operator, or to -1 if the DOF is not a SCUKN
    Eigen::VectorXi m_DOFtoGLUKN;     // Maps a dof i to its GL unknown in the global system, or to -1 if the DOF is a not a GLUKN

    // Scheme and solver parameters    
    double m_stab_par;
    
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

  
  //------------------------------------------------------------------------------
  // Trigonometric solution
  //------------------------------------------------------------------------------

  static Elasticity::DisplacementType
  trigonometric_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                             sin(2*PI*x.x())*cos(2*PI*x.y())*cos(2*PI*x.z()),
                                             0.5*sin(2*PI*x.y())*cos(2*PI*x.x())*cos(2*PI*x.z()),
                                             -sin(2*PI*x.z())*cos(2*PI*x.x())*cos(2*PI*x.y())
                                             );
                    };

  static Elasticity::GradientDisplacementType
  trigonometric_grad_u = [](const Eigen::Vector3d & x) -> Eigen::Matrix3d {
                        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << 2*PI*cos(2*PI*x.x())*cos(2*PI*x.y())*cos(2*PI*x.z()),
                                   -2*PI*sin(2*PI*x.x())*sin(2*PI*x.y())*cos(2*PI*x.z()), -2*PI*sin(2*PI*x.x())*sin(2*PI*x.z())*cos(2*PI*x.y());
                        M.row(1) << -1.0*PI*sin(2*PI*x.x())*sin(2*PI*x.y())*cos(2*PI*x.z()),
                                     1.0*PI*cos(2*PI*x.x())*cos(2*PI*x.y())*cos(2*PI*x.z()), -1.0*PI*sin(2*PI*x.y())*sin(2*PI*x.z())*cos(2*PI*x.x());
                        M.row(2) << 2*PI*sin(2*PI*x.x())*sin(2*PI*x.z())*cos(2*PI*x.y()), 
                                    2*PI*sin(2*PI*x.y())*sin(2*PI*x.z())*cos(2*PI*x.x()), -2*PI*cos(2*PI*x.x())*cos(2*PI*x.y())*cos(2*PI*x.z());
                        return M;
                         };

  static Elasticity::ForcingTermType
  trigonometric_div_eps_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                -7.0*pow(PI, 2)*sin(2*PI*x.x())*cos(2*PI*x.y())*cos(2*PI*x.z()),
                                -4.0*pow(PI, 2)*sin(2*PI*x.y())*cos(2*PI*x.x())*cos(2*PI*x.z()),
                                5.0*pow(PI, 2)*sin(2*PI*x.z())*cos(2*PI*x.x())*cos(2*PI*x.y())
                                ); 
                    };

  static Elasticity::ForcingTermType
  trigonometric_grad_div_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                -2.0*pow(PI, 2)*sin(2*PI*x.x())*cos(2*PI*x.y())*cos(2*PI*x.z()),
                                -2.0*pow(PI, 2)*sin(2*PI*x.y())*cos(2*PI*x.x())*cos(2*PI*x.z()),
                                -2.0*pow(PI, 2)*sin(2*PI*x.z())*cos(2*PI*x.x())*cos(2*PI*x.y())
                                );  
                    };

  static Elasticity::ForcingTermType
  trigonometric_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return -2 * mu * trigonometric_div_eps_u(mu,lambda,x) - lambda * trigonometric_grad_div_u(mu,lambda,x); 
                    };

  //------------------------------------------------------------------------------
  // Constant solution
  //------------------------------------------------------------------------------
  static Elasticity::DisplacementType
  constant_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(1., 1., 1.);
                    };

  static Elasticity::GradientDisplacementType
  constant_grad_u = [](const Eigen::Vector3d & x) -> Eigen::Matrix3d {
                        return Eigen::Matrix3d::Zero();
                         };

  static Elasticity::ForcingTermType
  constant_div_eps_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d::Zero(); 
                    };

  static Elasticity::ForcingTermType
  constant_grad_div_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d::Zero();  
                    };

  static Elasticity::ForcingTermType
  constant_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return -2 * mu * constant_div_eps_u(mu,lambda,x) - lambda * constant_grad_div_u(mu,lambda,x); 
                    };

  //------------------------------------------------------------------------------
  // Linear solution
  //------------------------------------------------------------------------------
  static Elasticity::DisplacementType
  linear_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(x.x(), -x.y(), x.z());
                    };

  static Elasticity::GradientDisplacementType
  linear_grad_u = [](const Eigen::Vector3d & x) -> Eigen::Matrix3d {
                        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << 1, 0, 0;
                        M.row(1) << 0, -1, 0;
                        M.row(2) << 0, 0, 1;
                        
                        return M;
                         };

  static Elasticity::ForcingTermType
  linear_div_eps_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d::Zero(); 
                    };

  static Elasticity::ForcingTermType
  linear_grad_div_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d::Zero();  
                    };

  static Elasticity::ForcingTermType
  linear_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return -2 * mu * linear_div_eps_u(mu,lambda,x) - lambda * linear_grad_div_u(mu,lambda,x); 
                    };

  //------------------------------------------------------------------------------
  // Quadratic solution
  //------------------------------------------------------------------------------
  static Elasticity::DisplacementType
  quadratic_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(pow(x.x(), 2), x.y()*x.z(), x.x()*x.z() + pow(x.z(), 2));
                    };

  static Elasticity::GradientDisplacementType
  quadratic_grad_u = [](const Eigen::Vector3d & x) -> Eigen::Matrix3d {
                        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << 2*x.x(), 0, 0;
                        M.row(1) << 0, x.z(), x.y();
                        M.row(2) << x.z(), 0, x.x() + 2*x.z();
                        
                        return M;
                         };

  static Elasticity::ForcingTermType
  quadratic_div_eps_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(2.5, 0, 2.5); 
                    };

  static Elasticity::ForcingTermType
  quadratic_grad_div_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(3, 0, 3);  
                    };

  static Elasticity::ForcingTermType
  quadratic_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return -2 * mu * quadratic_div_eps_u(mu,lambda,x) - lambda * quadratic_grad_div_u(mu,lambda,x); 
                    };

//  //------------------------------------------------------------------------------
//  // Vertical flow
//  //------------------------------------------------------------------------------
//
//  constexpr double scal_u = 1.;
//  static NavierStokes::VelocityType
//  vertical_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
//                      return scal_u * Eigen::Vector3d( cos(PI*x.y()), cos(PI*x.x()), 1.);
//                    };
//
//  static NavierStokes::VorticityType
//  vertical_curl_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
//                      return scal_u * Eigen::Vector3d(0, 0, -PI*sin(PI*x.x()) + PI*sin(PI*x.y()) );
//                         };
//
//  static NavierStokes::PressureType
//  vertical_p = trigonometric_p;
//
//  static NavierStokes::PressureGradientType
//  vertical_grad_p = trigonometric_grad_p;
//
//  static NavierStokes::ForcingTermType
//  vertical_curl_u_cross_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
//                      return std::pow(scal_u, 2) 
//                           * Eigen::Vector3d (
//                                             -(-PI*sin(PI*x.x()) + PI*sin(PI*x.y()))*cos(PI*x.x()),
//                                             (-PI*sin(PI*x.x()) + PI*sin(PI*x.y()))*cos(PI*x.y()), 
//                                             0.
//                                             );
//                    };
//
//  static NavierStokes::ForcingTermType
//  vertical_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
//                      return (1./reynolds) * scal_u * Eigen::Vector3d(pow(PI, 2)*cos(PI*x.y()), pow(PI, 2)*cos(PI*x.x()), 0.)
//                        + navier_scaling * vertical_curl_u_cross_u(x)
//                        + vertical_grad_p(x);
////                        return Eigen::Vector3d::Zero();
//                    };
//
//
//  //------------------------------------------------------------------------------
//  // Impose BC: p=1 at entrance x=0, u.n=1 at exit x=1 on lower bottom corner
//  //------------------------------------------------------------------------------
//  static NavierStokes::VelocityType
//  pressflux_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
//                      return ( lower_bottom_corner_x_one.is_in(x) ? VectorRd(1., 0., 0.) : VectorRd::Zero());
//                    };
//
//  static NavierStokes::VorticityType
//  pressflux_curl_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
//                          return Eigen::Vector3d::Zero();
//                        };
//
//  static NavierStokes::PressureType
//  pressflux_p = [](const Eigen::Vector3d & x) -> double {
//                          // return (x.x()<1e-8 ? 1. : 0.);
//                          return  pressure_scaling * (-x.z());
//                        };
//
//  static NavierStokes::PressureGradientType
//  pressflux_grad_p = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
//                         return Eigen::Vector3d::Zero();
//                      };
//
//  static NavierStokes::ForcingTermType
//  pressflux_curl_u_cross_u = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
//                      return Eigen::Vector3d::Zero();
//                    };
//
//  static NavierStokes::ForcingTermType
//  pressflux_f = [](const Eigen::Vector3d & x) -> Eigen::Vector3d {
//                      return Eigen::Vector3d::Zero();
//                    };


  /**@}*/

} // end of namespace HArDCore3D

#endif
