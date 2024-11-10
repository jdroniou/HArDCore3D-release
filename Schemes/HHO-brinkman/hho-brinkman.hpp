// Implementation of the HHO scheme, with full gradient, for the Brinkman equations -div(mu grad u) + nu u + grad p=f, div u = g.
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
 *
 *  This implementation of HHO was developped following the principles described in 
 * Appendix B of the book
 *
 * The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications. 
 *  D. A. Di Pietro and J. Droniou. Modeling, Simulation and Applications, vol. 19. 
 *  Springer International Publishing, 2020, xxxi + 525p. doi: 10.1007/978-3-030-37203-3. 
 *  url: https://hal.archives-ouvertes.fr/hal-02151813.
 *
 * If you use this code or part of it for a scientific publication, please cite the book
 * above as a reference for the implementation.
 *
 */


#ifndef HHOBRINKMAN_HPP
#define HHOBRINKMAN_HPP

#include <iostream>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include <mesh.hpp>
#include <mesh_builder.hpp>
#include <local_static_condensation.hpp>
#include <linearsolver.hpp>

#include <vhhospace.hpp>
#include <BoundaryConditions/BoundaryConditions.hpp>
#include <integralweight.hpp>

/*!
 * @defgroup HHO_brinkman
 * @brief Implementation of the HHO scheme for the Stokes and Brinkman equations
 */

namespace HArDCore3D
{

  /*!
   * @addtogroup HHO_brinkman
   * @{
   */

  //------------------------------------------------------------------------------
  // Brinkman model
  //------------------------------------------------------------------------------

  /// Structure to store physical parameter (and compute derived parameters)
  struct BrinkmanParameters
  {
    typedef IntegralWeight ViscosityType;
    typedef IntegralWeight PermeabilityInvType;

    /// Constructor
    BrinkmanParameters(
      const ViscosityType & _mu,   ///< Viscosity
      const PermeabilityInvType & _kappainv    ///< Inverse of permeability
      )
      : mu(_mu),
        kappainv(_kappainv)
      {
        // do nothing
      }
   
    /// Returns the friction coefficient
    double Cf(const Cell & T) const
    {
      double muT = mu.value(T, T.center_mass());
      double kappainvT = kappainv.value(T, T.center_mass());
      
      // Thres is the level at which the parameter is put at 0; note that Cf is also bounded below by thres
      // to avoid singularities when taking Cf^{-1}.
      double thres = 1e-14;
      return (muT > thres ? std::max(thres, kappainvT * std::pow(T.diam(), 2) / muT) : 1./thres);
    }
   
    ViscosityType mu;
    PermeabilityInvType kappainv;
  };

  /// Structure to store norm components (energy for velocity, L2 pressure, and global)
  struct BrinkmanNorms
  {
    /// Constructor
    BrinkmanNorms(double norm_u, double norm_p):
      u(norm_u),
      p(norm_p),
      energy( std::sqrt( std::pow(norm_u, 2) + std::pow(norm_p, 2) ) )
      {
        // do nothing
      };
    
    double u; ///< Energy norm of velocity
    double p; ///< L2 norm of p
    double energy; ///< Global energy norm  
  };


  /// Structure for Brinkman model
  /** The global unknowns are in the following order:
        - Dirichlets dofs of the velocity (Dirichlet faces must be numbered first in the mesh)
        - Face and then element unknowns of the velocity
        - Element unknowns of the pressure
        - Lagrange multiplier
     The statically condenseds dofs for the velocity are the element unknowns of the velocity, and all dofs but
     the first one in each element for the pressure.
  */
  struct Brinkman
  {
    typedef Eigen::SparseMatrix<double> SystemMatrixType;
    
    typedef std::function<VectorRd(const VectorRd &)> MomentumForcingTermType;
    typedef std::function<double(const VectorRd &)> CompressibilityForcingTermType;
    typedef std::function<VectorRd(const VectorRd &)> VelocityType;
    typedef std::function<MatrixRd(const VectorRd &)> VelocityGradientType;
    typedef std::function<double(const VectorRd &)> PressureType;
    typedef std::function<VectorRd(const VectorRd &)> PressureGradientType;

    /// Constructor
    Brinkman(
           const VHHOSpace & vhho_space,      ///< Velocity space and operators
           const GlobalDOFSpace & p_space,    ///< Space for the pressure
           const BoundaryConditions & BC,     ///< Boundary conditions
           bool use_threads,                  ///< True for parallel execution, false for sequential execution
           std::ostream & output = std::cout  ///< Output stream to print status messages
           );

    /// Assemble the global system    
    void assembleLinearSystem(
                    const MomentumForcingTermType & f,      ///< Forcing term
                    const CompressibilityForcingTermType & g,      ///< Forcing term for the divergence equation
                    const BrinkmanParameters & para,       ///< Viscosity and permeability
                    const VelocityType & u,         ///< Exact solution for boundary condition
                    Eigen::VectorXd & UDir          ///< Vector filled in by Dirichlets BCs
                    );
                    
    /// Returns the local number of velocity statically condensed DOFs
    inline size_t nloc_sc_u() const
    {
      return m_nloc_sc_u;
    }

    /// Returns the local number of pressure statically condensed DOFs
    inline size_t nloc_sc_p() const
    {
      return m_nloc_sc_p;
    }

    /// Returns the number of velocity statically condensed DOFs
    inline size_t numSCDOFs_u() const
    {
      return mesh().n_cells() * m_nloc_sc_u;
    }

    /// Returns the number of pressure statically condensed DOFs
    inline size_t numSCDOFs_p() const
    {
      return mesh().n_cells() * m_nloc_sc_p;
    }

    /// Returns the number of statically condensed DOFs
    inline size_t numSCDOFs() const
    {
      return numSCDOFs_u() + numSCDOFs_p();
    }

    /// Returns the number of Dirichlet DOFs
    inline size_t numDirDOFs() const
    {
      return m_BC.n_dir_faces()*m_vhhospace.numLocalDofsFace();
    }
    
    /// Returns the dimension of velocity space
    inline size_t dimVelocity() const
    {
      return m_vhhospace.dimension();
    }

    /// Returns the dimension of pressure space
    inline size_t dimPressure() const
    {
      return m_pspace.dimension();
    }
    
    /// Returns the number of DOFs after SC and with Lagrange multiplier, but before eliminating Dirichlet DOFs
    inline size_t numNonSCDOFs() const
    {
      return dimVelocity() + dimPressure() - numSCDOFs() + 1;
    }

    /// Returns the size of the final system with Lagrange multiplier, after application of SC and removal of Dirichlet BCs
    inline size_t sizeSystem() const
    {
      return numNonSCDOFs() - numDirDOFs();
    }

    /// Returns the velocity space
    inline const VHHOSpace & vhhospace() const
    {
      return m_vhhospace;
    }
    
    /// Returns the pressure space
    inline const GlobalDOFSpace & pspace() const
    {
      return m_pspace;
    }

    /// Returns the mesh
    inline const Mesh & mesh() const
    {
      return m_vhhospace.mesh();
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

    /// Returns the stabilization parameter (scaling)
    inline const std::pair<double,double> & stabilizationParameter() const {
      return m_stab_par;
    }

    /// Returns the stabilization parameter
    inline std::pair<double,double> & stabilizationParameter() {
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

    /// Interpolates velocity and pressure
    Eigen::VectorXd interpolate(
          const VelocityType & u, ///< Continuous velocity to interpolate
          const PressureType & p, ///< Continuous pressure to interpolate
          const int doe_cell = -1, ///< The optional degre of cell quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_face = -1 ///< The optional degre of face quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          ) const;    

    /// Compute the discrete energy norm of a family of vectors representing the dofs
    std::vector<BrinkmanNorms> computeBrinkmanNorms(
                       const std::vector<Eigen::VectorXd> & list_dofs,  ///< The list of vectors representing the dofs
                       const BrinkmanParameters & para                  ///< Viscosity and permeability
                      ) const;

    /// Create vertex values for the pressure (from the element values), for plotting
    std::vector<double> pressureVertexValues(
                       const Eigen::VectorXd & p    ///< Vector of pressure DOFs
                       ) const;

    /// Computes the condition number of the matrix
    std::pair<double, double> computeConditionNum() const;

    /// Compute the velocity flux (integral of u.n) across a given surface, and the area of the surface, assumed to be star-shaped with respect the barycenter of its vertices
    std::pair<double,double> computeFlux(
                   const Eigen::VectorXd & u,    ///< Vector of velocity DOFs
                   const std::pair<std::vector<VectorRd>, VectorRd> & surf   ///< Surface described by coplanar vertices (in order around the boundary) and one unit normal vector to indicate the flux to compute
                   ) const;

  private:
    /// Compute the local contribution for the element of index iT
    std::pair<Eigen::MatrixXd, Eigen::VectorXd>
          _compute_local_contribution(
                                size_t iT,                            ///< Element index
                                const MomentumForcingTermType & f,      ///< Forcing term
                                const CompressibilityForcingTermType & g,    ///< Forcing term for divergence equation
                                const BrinkmanParameters & para              ///< Viscosity and permeability
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
    
    /// Spaces for velocity and pressure. The pressure one uses the polynomial bases provided by vhhospace
    const VHHOSpace & m_vhhospace;
    const GlobalDOFSpace & m_pspace; 

    const BoundaryConditions & m_BC;
    bool m_use_threads;
    std::ostream & m_output;
    const size_t m_nloc_sc_u;   // Number of statically condensed velocity DOFs in each cell
    const size_t m_nloc_sc_p;   // Number of statically condensed pressure DOFs in each cell
    SystemMatrixType m_A;   // Matrix and RHS of statically condensed system
    Eigen::VectorXd m_b;
    SystemMatrixType m_sc_A; // Static condensation operator and RHS (to recover statically condensed DOFs)
    Eigen::VectorXd m_sc_b;
    std::pair<double,double> m_stab_par; // Stabilisation parameters: first for Stokes, second for Darcy
  };

  
  //------------------------------------------------------------------------------
  // Exact solutions
  //------------------------------------------------------------------------------

  static const double PI = boost::math::constants::pi<double>();
  using std::sin;
  using std::cos;

  double scaling_mu = 1.;
  double scaling_kappainv = 1.;
  
  // Threshold for checking that quantities are zero or not
  double constexpr eps=1e-14;

  //---------------- Linear -----------------------------------------------------------
  
  
  static BrinkmanParameters::ViscosityType
  linear_mu = BrinkmanParameters::ViscosityType(1.);
  
  static BrinkmanParameters::PermeabilityInvType
  linear_kappainv = BrinkmanParameters::PermeabilityInvType(1.);

  static const MatrixRd mat_u = (MatrixRd() << 1.,2.,3., 0.,-1.,1. , 0.,1.,1.).finished();
  static const VectorRd vec_p = VectorRd(-1., 2., -5.);
    
  static Brinkman::VelocityType
  linear_u = [](const VectorRd & x) -> VectorRd {
                 return mat_u * x;
               };
               
  static Brinkman::PressureType
  linear_p = [](const VectorRd & x) -> double {
                 return vec_p.dot(x - VectorRd(.5,.5,.5));  // Assuming domain (0,1)^3 to have zero average
               };

  static Brinkman::VelocityGradientType  
  linear_gradu = [](const VectorRd & x) -> MatrixRd {
                     return mat_u;
                   };

  static Brinkman::PressureGradientType  
  linear_gradp = [](const VectorRd & x) -> VectorRd {
                     return vec_p;
                   };

  static Brinkman::MomentumForcingTermType
  linear_f = [](const VectorRd & x) -> VectorRd {
                 return VectorRd::Zero() + scaling_kappainv*linear_u(x) + linear_gradp(x);
               };

  static Brinkman::CompressibilityForcingTermType 
  linear_g = [](const VectorRd & x)->double { return linear_gradu(x).trace();};
  

  //---------------- Quadratic ------------------------------------------------
  static BrinkmanParameters::ViscosityType
  quadratic_mu = BrinkmanParameters::ViscosityType(1.);

  static BrinkmanParameters::PermeabilityInvType
  quadratic_kappainv = BrinkmanParameters::PermeabilityInvType(1.);
  
  static Brinkman::VelocityType
  quadratic_u = [](const VectorRd & x) -> VectorRd {
                      return VectorRd(
                                     std::pow(x(0), 2),
                                     std::pow(x(1), 2),
                                     std::pow(x(2), 2)
                                     );
                    };

  static Brinkman::VelocityGradientType
  quadratic_gradu = [](const VectorRd & x) -> MatrixRd {
                      MatrixRd M = MatrixRd::Zero(); 
                      M.row(0) << 2. * x(0),
                                    0.,
                                    0.;
                      M.row(1) << 0.,
                                  2 * x(1),
                                  0.;
                      M.row(2) << 0.,
                                  0.,
                                  2 * x(2);
                      
                      return M;
                    };

  static Brinkman::PressureType
  quadratic_p = [](const VectorRd & x) -> double {
                      return  std::pow(x(0), 2) + x(1)*x(2) - 1./3. - 1./4.;
                    };

  static Brinkman::PressureGradientType
  quadratic_gradp = [](const VectorRd & x) -> VectorRd {
                      return  VectorRd( 2 * x(0),
                                        x(2),
                                        x(1)
                                        );
                         };

  static Brinkman::MomentumForcingTermType
  quadratic_f = [](const VectorRd & x) -> VectorRd {
                      return - scaling_mu * VectorRd( 2., 2., 2.) + scaling_kappainv * quadratic_u(x) + quadratic_gradp(x);
                    };

  static Brinkman::CompressibilityForcingTermType 
  quadratic_g = [](const VectorRd & x)->double { return quadratic_gradu(x).trace();};


  //---------------- Cubic ------------------------------------------------
  static BrinkmanParameters::ViscosityType
  cubic_mu = BrinkmanParameters::ViscosityType(1.);

  static BrinkmanParameters::PermeabilityInvType
  cubic_kappainv = BrinkmanParameters::PermeabilityInvType(1.);
  
  static Brinkman::VelocityType
  cubic_u = [](const VectorRd & x) -> VectorRd {
                      return VectorRd(
                                     std::pow(x(0), 3) + std::pow(x(0), 2)*x(1),
                                     std::pow(x(1), 2) + x(0) * x(1) * x(2),
                                     std::pow(x(2), 3)
                                     );
                    };

  static Brinkman::VelocityGradientType
  cubic_gradu = [](const VectorRd & x) -> MatrixRd {
                      MatrixRd M = MatrixRd::Zero(); 
                      M.row(0) << 3. * std::pow(x(0), 2) + 2 * x(0) * x(1),
                                    std::pow(x(0), 2),
                                    0.;
                      M.row(1) << x(1) * x(2),
                                  2 * x(1) + x(0) * x(2),
                                  x(0) * x(1);
                      M.row(2) << 0.,
                                  0.,
                                  3 * std::pow(x(2), 2);
                      
                      return M;
                    };

  static Brinkman::PressureType
  cubic_p = [](const VectorRd & x) -> double {
                      return  std::pow(x(0), 3) + x(1)*x(2) - 1./4. - 1./4.;
                    };

  static Brinkman::PressureGradientType
  cubic_gradp = [](const VectorRd & x) -> VectorRd {
                      return  VectorRd( 3 * std::pow(x(0), 2),
                                        x(2),
                                        x(1)
                                        );
                         };

  static Brinkman::MomentumForcingTermType
  cubic_f = [](const VectorRd & x) -> VectorRd {
                      return - scaling_mu * VectorRd( 6 * x(0) + 2 * x(1), 2., 6 * x(2)) + scaling_kappainv * cubic_u(x) + cubic_gradp(x);
                    };

  static Brinkman::CompressibilityForcingTermType 
  cubic_g = [](const VectorRd & x)->double { return cubic_gradu(x).trace();};


  //---------------- Trigonometric ------------------------------------------------
  static BrinkmanParameters::ViscosityType
  trigonometric_mu = BrinkmanParameters::ViscosityType(1.);

  static BrinkmanParameters::PermeabilityInvType
  trigonometric_kappainv = BrinkmanParameters::PermeabilityInvType(1.);

  double constexpr pressure_scaling = 1.;
  
  static Brinkman::VelocityType
  trigonometric_u = [](const VectorRd & x) -> VectorRd {
                      return VectorRd(
                                     0.5 * sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                     0.5 * cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                     -cos(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2))
                                     );
                    };

  static Brinkman::VelocityGradientType
  trigonometric_gradu = [](const VectorRd & x) -> MatrixRd {
                      MatrixRd M = MatrixRd::Zero(); 
                      M.row(0) << PI * cos(2 * PI * x(0)) * cos(2 * PI * x(1)) * cos(2 * PI * x(2)),
                                    -PI * sin(2. * PI * x(0)) * sin(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                    -PI * sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2));
                      M.row(1) << - PI * sin(2. * PI * x(0)) * sin(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                    PI * cos(2. * PI * x(0)) * cos(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                    -PI * cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * sin(2. * PI * x(2));
                      M.row(2) << 2 * PI * sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                    2 * PI * cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                    - 2 * PI * cos(2. * PI * x(0)) * cos(2. * PI * x(1)) * cos(2. * PI * x(2));
                      
                      return M;
                    };

  static Brinkman::PressureType
  trigonometric_p = [](const VectorRd & x) -> double {
                      return pressure_scaling * sin(2. * PI * x(0)) * sin(2. * PI * x(1)) * sin(2. * PI * x(2));
                    };

  static Brinkman::PressureGradientType
  trigonometric_gradp = [](const VectorRd & x) -> VectorRd {
                      return pressure_scaling * 2. * PI * VectorRd(
                                                      cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                                      sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                                      sin(2. * PI * x(0)) * sin(2. * PI * x(1)) * cos(2. * PI * x(2))
                                                      );
                         };

  static Brinkman::MomentumForcingTermType
  trigonometric_f = [](const VectorRd & x) -> VectorRd {
                      return scaling_mu * 6. * std::pow(PI, 2) * VectorRd(
                                                            sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                                            cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                                            -2. * cos(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2))
                                                            )
                        + scaling_kappainv*trigonometric_u(x) + trigonometric_gradp(x);
                    };

  static Brinkman::CompressibilityForcingTermType 
  trigonometric_g = [](const VectorRd & x)->double { return trigonometric_gradu(x).trace();};



  //---------------- Various regimes parametrised by global friction coefficient --------------------------------------------

  static BrinkmanParameters::ViscosityType
  regimes_mu = BrinkmanParameters::ViscosityType(1.);

  static BrinkmanParameters::PermeabilityInvType
  regimes_kappainv = BrinkmanParameters::PermeabilityInvType(1.);

  static Brinkman::PressureType
  regimes_p = trigonometric_p;

  static Brinkman::PressureGradientType
  regimes_gradp = trigonometric_gradp;

  // u_D = - kappainv^{-1} grad p
  static Brinkman::VelocityType
  regimes_uD = [](const VectorRd & x) -> VectorRd {
                      VectorRd uDvalue = VectorRd::Zero();
                      if (scaling_kappainv > eps) {
                        uDvalue = - std::pow(scaling_kappainv, -1) * regimes_gradp(x) ;
                      }
                      return uDvalue;
                     };
                       
  static Brinkman::VelocityGradientType
  regimes_graduD = [](const VectorRd & x) -> MatrixRd {
                  MatrixRd M = MatrixRd::Zero(); 
                  
                  if (scaling_kappainv>eps){
                      M.row(0) << sin(2. * PI * x(0)) * sin(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                    -cos(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                    -cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * cos(2. * PI * x(2));
                      M.row(1) <<  -cos(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                   sin(2. * PI * x(0)) * sin(2. * PI * x(1)) * sin(2. * PI * x(2)),
                                   -sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * cos(2. * PI * x(2));
                      M.row(2) <<  -cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                    -sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                    sin(2. * PI * x(0)) * sin(2. * PI * x(1)) * sin(2. * PI * x(2));
                      M *= 4. * std::pow(PI, 2) / scaling_kappainv;
                   }

                   return M;
                };

  // u_S = trigonometric u
  static Brinkman::VelocityType regimes_uS = trigonometric_u;
  static Brinkman::VelocityGradientType regimes_graduS = trigonometric_gradu;

  static std::function<double(const VectorRd &)> 
    XiS = [](const VectorRd & x) -> double {
              return (scaling_mu > eps ? std::exp(-scaling_kappainv/scaling_mu) : 0.);
          };
              
  static Brinkman::VelocityType 
    regimes_u = [](const VectorRd & x) -> VectorRd { return XiS(x) * regimes_uS(x) + (1.-XiS(x)) * regimes_uD(x);};

  static Brinkman::VelocityGradientType 
    regimes_gradu = [](const VectorRd & x) -> MatrixRd { return XiS(x) * regimes_graduS(x) + (1.-XiS(x)) * regimes_graduD(x);};

  static Brinkman::MomentumForcingTermType
  regimes_f = [](const VectorRd & x) -> VectorRd {
                      VectorRd fvalue = XiS(x) * scaling_mu * 6. * std::pow(PI, 2) * VectorRd(
                                                            sin(2. * PI * x(0)) * cos(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                                            cos(2. * PI * x(0)) * sin(2. * PI * x(1)) * cos(2. * PI * x(2)),
                                                            -2. * cos(2. * PI * x(0)) * cos(2. * PI * x(1)) * sin(2. * PI * x(2))
                                                            )
                                        + regimes_gradp(x);
                      
                      if (scaling_kappainv > eps){
                        fvalue += scaling_kappainv * regimes_u(x)
                                  - scaling_mu * (1-XiS(x))/scaling_kappainv * 12. * std::pow(PI, 2) * regimes_gradp(x);
                      }

                      return fvalue;
                    };

  static Brinkman::CompressibilityForcingTermType 
  regimes_g = [](const VectorRd & x)->double { return regimes_gradu(x).trace() ;};


  //---------------- V-cracked boundary conditions and data (no exact solution here)------------------------------------

  // Returns 1 in the crack, 0 outside
  static std::function<double(const VectorRd &)>
      XiV = [](const VectorRd & x)->double { return (x.z() > 2*std::abs(x.x()-1.)+0.5-eps); };

  // Viscosity only in crack
  static BrinkmanParameters::ViscosityType
  vcracked_mu = BrinkmanParameters::ViscosityType([](const VectorRd & x)->double { return 0.01 * XiV(x);}); 

  // Permeability only outside crack
  static BrinkmanParameters::PermeabilityInvType
  vcracked_kappainv = BrinkmanParameters::PermeabilityInvType([](const VectorRd & x)->double { return 1e-5*(1.-XiV(x)) + 0.001*XiV(x);});

  // u: entry on the left (surface (0,0,0)-(0,0.2,0)-(0,0.2,1)-(0,0,1), 
  //    and "wind" on the top above the crack  (surface (0.75,0,1)-(1.25,0,1)-(1.25,0.2,1)-(0.75,0.2,1))
  static Brinkman::VelocityType
  vcracked_u = [](const VectorRd & x) -> VectorRd {
                      VectorRd value = VectorRd::Zero();
                      if ( x.x() < .2) {
                        // Entry on the left
                        value = VectorRd(.3, 0., 0.);
                      }else if ( (x.z() > 0.8) && (x.x() > 0.75-eps) && (x.x() < 1.25+eps) ){
                        // Wind top
                        value = VectorRd(1, 0., 0.);
                      }
                      
                      return value;
                     };


  static Brinkman::VelocityGradientType
  vcracked_gradu = [](const VectorRd & x) -> MatrixRd {
                      return MatrixRd::Zero();
                    };

  static Brinkman::PressureType
  vcracked_p = [](const VectorRd & x) -> double {
                      return  0.;
                    };

  static Brinkman::PressureGradientType
  vcracked_gradp = [](const VectorRd & x) -> VectorRd {
                      return  VectorRd::Zero();
                         };

  // Source is gravity
  static Brinkman::MomentumForcingTermType
  vcracked_f = [](const VectorRd & x) -> VectorRd {
                      return VectorRd(0., 0., -.98);
                    };

  static Brinkman::CompressibilityForcingTermType 
  vcracked_g = [](const VectorRd & x) -> double { 
                      return 0.;
                    };

  //---------------- BCs for the cavity in porous medium ------------------------------------

  // Characteristic functions of cavity, wedge, and surrounding box
  static std::function<double(const VectorRd &)>
      XiCavity = [](const VectorRd & x)->double { 
            return ( (x.x()>0) && (x.x()<1) && (x.y()>0) && (x.y()<1) && (x.z()>-1) ? 1. : 0.); 
          };
  static std::function<double(const VectorRd &)>
      XiWedge = [](const VectorRd & x)->double { 
            return ( (x.x()>=1) && (x.x()<2) && (x.y()>0) && (x.y()<1) && (x.z()> -.75 + .25*(x.x()-1)) ? 1. : 0.); 
          };
  static std::function<double(const VectorRd &)>
      XiBox = [](const VectorRd & x)->double { return 1. - XiCavity(x) - XiWedge(x); };

  // Viscosity in cavity
  constexpr double viscosity_in_cavity = 0.01;
  static BrinkmanParameters::ViscosityType cavity_mu = viscosity_in_cavity * BrinkmanParameters::ViscosityType(XiCavity); 

  // Permeability in wedge and box
  constexpr double permeability_in_wedge = 1e-2;
  constexpr double permeability_in_box = 1e-7;
  static BrinkmanParameters::PermeabilityInvType cavity_kappainv 
        = std::pow(permeability_in_wedge, -1) * BrinkmanParameters::ViscosityType(XiWedge)
            + std::pow(permeability_in_box, -1) * BrinkmanParameters::ViscosityType(XiBox);

  // u: wind above cavity (only for BCs anyway).
  static Brinkman::VelocityType
  cavity_u = [](const VectorRd & x) -> VectorRd {
                      return XiCavity(x) * VectorRd(x.x()*(1.-x.x()), 0., 0.);
                     };

  static Brinkman::VelocityGradientType
  cavity_gradu = [](const VectorRd & x) -> MatrixRd {
                      return MatrixRd::Zero();
                    };

  static Brinkman::PressureType
  cavity_p = [](const VectorRd & x) -> double {
                      return  0.;
                    };

  static Brinkman::PressureGradientType
  cavity_gradp = [](const VectorRd & x) -> VectorRd {
                      return  VectorRd::Zero();
                         };

  // Source is gravity
  static Brinkman::MomentumForcingTermType
  cavity_f = [](const VectorRd & x) -> VectorRd {
                      return VectorRd(0., 0., -.98);
                    };

  static Brinkman::CompressibilityForcingTermType 
  cavity_g = [](const VectorRd & x) -> double { 
                      return 0.;
                    };


  //---------------- Two regions ------------------------------------------------
  // Two regions with different parameters, separated by x=1/2 (region S is x<1/2, D is x>1/2)
  //    u = alpha(y,z) + cos(pi x) beta_i(y,z) where i=S or D
  // This ensures that u is continuous, and that nu \nabla u n is continuous across x=1/2

  // Characteristic functions
  static std::function<double(const VectorRd &)>
      Xi_S = [](const VectorRd & x)->double { 
            return ( (x.x()<0.5) ? 1. : 0.); 
          };
  static std::function<double(const VectorRd &)>
      Xi_D = [](const VectorRd & x)->double { 
            return 1.-Xi_S(x); 
          };

  // Viscosity and permeability
  constexpr double viscosity_S = 1.;
  constexpr double viscosity_D = 0.;
  static BrinkmanParameters::ViscosityType tworegions_mu 
         = viscosity_S * BrinkmanParameters::ViscosityType(Xi_S)
            + viscosity_D * BrinkmanParameters::ViscosityType(Xi_D); 

  constexpr double permeability_S = 1e-7;
  constexpr double permeability_D = 1e-2;
  static BrinkmanParameters::PermeabilityInvType tworegions_kappainv 
        = std::pow(permeability_S, -1) * BrinkmanParameters::ViscosityType(Xi_S)
            + std::pow(permeability_D, -1) * BrinkmanParameters::ViscosityType(Xi_D);

  // Velocity (alpha and beta only depend on y,z)
  static std::function<VectorRd(const VectorRd &)>
      alpha = [](const VectorRd & x)->VectorRd {
        return VectorRd(
                        exp(-x.y()-x.z()),
                        sin(PI * x.y()) * sin(PI * x.z()),
                        x.y() * x.z()
                       );
      };

  static std::function<VectorRd(const VectorRd &)>
      beta_S = [](const VectorRd & x)->VectorRd {
        return cos(PI * x.x()) * (x.x() - 0.5)
                    * VectorRd(
                          x.y() + x.z(),
                          x.y() + cos(PI * x.z()),
                          sin(PI * x.y())
                         );
      };

  static std::function<VectorRd(const VectorRd &)>
      beta_D = [](const VectorRd & x)->VectorRd {
        return cos(PI * x.x()) * (x.x() - 0.5)
                    * VectorRd(
                        sin(PI * x.y()) + cos(PI * x.z()),
                        std::pow(x.z(), 3),
                        std::pow(x.y(), 2) * std::pow(x.z(), 2)
                       );
      }; 

  static Brinkman::VelocityType
  tworegions_u = [](const VectorRd & x) -> VectorRd {
                      return alpha(x) + Xi_S(x) * beta_S(x) + Xi_D(x) * beta_D(x);
                    };

  static Brinkman::VelocityGradientType
  grad_alpha = [](const VectorRd & x) -> MatrixRd {
                      MatrixRd M = MatrixRd::Zero(); 
                      M.row(0) << 0, -exp(-x.y() - x.z()), -exp(-x.y() - x.z());
                      M.row(1) << 0, PI*sin(PI*x.z())*cos(PI*x.y()), PI*sin(PI*x.y())*cos(PI*x.z());
                      M.row(2) << 0, x.z(), x.y();
                      
                      return M;
                    };

  static Brinkman::VelocityGradientType
  grad_beta_S = [](const VectorRd & x) -> MatrixRd {
                      MatrixRd M = MatrixRd::Zero(); 
                      M.row(0) << -PI*(x.x() - 0.5)*(x.y() + x.z())*sin(PI*x.x()) + (x.y() + x.z())*cos(PI*x.x()),
                                  (x.x() - 0.5)*cos(PI*x.x()), 
                                  (x.x() - 0.5)*cos(PI*x.x());
                      M.row(1) << -PI*(x.x() - 0.5)*(x.y() + cos(PI*x.z()))*sin(PI*x.x()) + (x.y() + cos(PI*x.z()))*cos(PI*x.x()),
                                  (x.x() - 0.5)*cos(PI*x.x()),
                                  -PI*(x.x() - 0.5)*sin(PI*x.z())*cos(PI*x.x());
                      M.row(2) << -PI*(x.x() - 0.5)*sin(PI*x.x())*sin(PI*x.y()) + sin(PI*x.y())*cos(PI*x.x()), 
                                  PI*(x.x() - 0.5)*cos(PI*x.x())*cos(PI*x.y()),
                                  0;

                      return M;
                    };

  static Brinkman::VelocityGradientType
  grad_beta_D = [](const VectorRd & x) -> MatrixRd {
                      MatrixRd M = MatrixRd::Zero(); 
                      M.row(0) << -PI*(x.x() - 0.5)*(sin(PI*x.y()) + cos(PI*x.z()))*sin(PI*x.x()) + (sin(PI*x.y()) + cos(PI*x.z()))*cos(PI*x.x()),
                                   PI*(x.x() - 0.5)*cos(PI*x.x())*cos(PI*x.y()), 
                                   -PI*(x.x() - 0.5)*sin(PI*x.z())*cos(PI*x.x());
                      M.row(1) << -PI*pow(x.z(), 3)*(x.x() - 0.5)*sin(PI*x.x()) + pow(x.z(), 3)*cos(PI*x.x()), 
                                  0, 
                                  3*pow(x.z(), 2)*(x.x() - 0.5)*cos(PI*x.x());
                      M.row(2) << -PI*pow(x.y(), 2)*pow(x.z(), 2)*(x.x() - 0.5)*sin(PI*x.x()) + pow(x.y(), 2)*pow(x.z(), 2)*cos(PI*x.x()),
                                  2*x.y()*pow(x.z(), 2)*(x.x() - 0.5)*cos(PI*x.x()), 
                                  2*pow(x.y(), 2)*x.z()*(x.x() - 0.5)*cos(PI*x.x());
                      
                      return M;
                    };

  static Brinkman::VelocityGradientType
  tworegions_gradu = [](const VectorRd & x) -> MatrixRd {
                      return grad_alpha(x) + Xi_S(x) * grad_beta_S(x) + Xi_D(x) * grad_beta_D(x);
                    };

  static Brinkman::PressureType
  tworegions_p = trigonometric_p;

  static Brinkman::PressureGradientType
  tworegions_gradp = trigonometric_gradp;

//  static Brinkman::PressureType
//  tworegions_p = [](const VectorRd & x) -> double {
//                      return  (x.x()-0.5) * x.y() * x.z();
//                    };

//  static Brinkman::PressureGradientType
//  tworegions_gradp = [](const VectorRd & x) -> VectorRd {
//                      return  VectorRd( x.y() * x.z(), (x.x()-0.5) * x.z(), (x.x()-0.5) * x.y());
//                    };

  static Brinkman::MomentumForcingTermType
  DivGrad_alpha = [](const VectorRd & x) -> VectorRd {
         return VectorRd(
                      2*exp(-x.y() - x.z()),
                  	 -2*pow(PI, 2)*sin(PI*x.y())*sin(PI*x.z()),
	                    0 
                     );
          };
                     
  static Brinkman::MomentumForcingTermType
  DivGrad_beta_S = [](const VectorRd & x) -> VectorRd {
         return VectorRd(
                      -pow(PI, 2)*(x.x() - 0.5)*(x.y() + x.z())*cos(PI*x.x()) - 2*PI*(x.y() + x.z())*sin(PI*x.x()), 
                    	-pow(PI, 2)*(x.x() - 0.5)*(x.y() + cos(PI*x.z()))*cos(PI*x.x()) - pow(PI, 2)*(x.x() - 0.5)*cos(PI*x.x())*cos(PI*x.z()) - 2*PI*(x.y() + cos(PI*x.z()))*sin(PI*x.x()),
                      -2*pow(PI, 2)*(x.x() - 0.5)*sin(PI*x.y())*cos(PI*x.x()) - 2*PI*sin(PI*x.x())*sin(PI*x.y()) 
                     );
 
          };

  static Brinkman::MomentumForcingTermType
  DivGrad_beta_D = [](const VectorRd & x) -> VectorRd {
         return VectorRd(
                       -pow(PI, 2)*(x.x() - 0.5)*(sin(PI*x.y()) + cos(PI*x.z()))*cos(PI*x.x()) - pow(PI, 2)*(x.x() - 0.5)*sin(PI*x.y())*cos(PI*x.x()) - pow(PI, 2)*(x.x() - 0.5)*cos(PI*x.x())*cos(PI*x.z()) - 2*PI*(sin(PI*x.y()) + cos(PI*x.z()))*sin(PI*x.x()),
                       -pow(PI, 2)*pow(x.z(), 3)*(x.x() - 0.5)*cos(PI*x.x()) - 2*PI*pow(x.z(), 3)*sin(PI*x.x()) + 6*x.z()*(x.x() - 0.5)*cos(PI*x.x()),
                       -pow(PI, 2)*pow(x.y(), 2)*pow(x.z(), 2)*(x.x() - 0.5)*cos(PI*x.x()) - 2*PI*pow(x.y(), 2)*pow(x.z(), 2)*sin(PI*x.x()) + 2*pow(x.y(), 2)*(x.x() - 0.5)*cos(PI*x.x()) + 2*pow(x.z(), 2)*(x.x() - 0.5)*cos(PI*x.x())
                     );
          };

  static Brinkman::MomentumForcingTermType
  tworegions_f = [](const VectorRd & x) -> VectorRd {
         return - scaling_mu * ( viscosity_S * Xi_S(x) * (DivGrad_alpha(x) + DivGrad_beta_S(x)) + 
                                  viscosity_D * Xi_D(x) * (DivGrad_alpha(x) + DivGrad_beta_D(x)) )
                 + scaling_kappainv * 
                      ( std::pow(permeability_S, -1) * Xi_S(x) + std::pow(permeability_D, -1) * Xi_D(x) ) * tworegions_u(x)
                 + tworegions_gradp(x);
                    };

  static Brinkman::CompressibilityForcingTermType 
  tworegions_g = [](const VectorRd & x)->double { return tworegions_gradu(x).trace();};


} // end of namespace HArDCore3D

#endif
