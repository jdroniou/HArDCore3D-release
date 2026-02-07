// Implementation of fracture scheme for the elasticity equation
//
// Author: Jerome Droniou (jerome.droniou@umontpellier.fr)
//


#ifndef DDR_ELASTICITY_FRACTURE_HPP
#define DDR_ELASTICITY_FRACTURE_HPP

#include <iostream>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include <mesh.hpp>
#include <mesh_builder.hpp>
#include <mesh_reader_fracture.hpp>
#include <local_static_condensation.hpp>
#include <linearsolver.hpp>


#include <BoundaryConditions/BoundaryConditions.hpp>
#include <BoundaryConditions/BChandlers.hpp>   

#include <vsxgrad.hpp>
#include <connected_information_generator.hpp>
#include<fracturedofspace.hpp>
#include <fractureinterpolator.hpp>
#include <vtu_facewriter.hpp>

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
    
    // typedef std::function<Eigen::Vector3d(const double &, const double &, const Eigen::Vector3d &)> ForcingTermType;
    // typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> DisplacementType;
    // typedef std::function<Eigen::Matrix3d(const Eigen::Vector3d &)> GradientDisplacementType;

    typedef std::function<Eigen::Vector3d(const double &, const double &, const Eigen::Vector3d &, const Cell *)> ForcingTermType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &, const Cell *, const double & lambda)> DisplacementType;
    typedef std::function<double(const Face *, const double & lambda)> FrictionType;
    typedef std::function<Eigen::Matrix3d(const Eigen::Vector3d &, const Cell *, const double & lambda)> GradientDisplacementType;
    typedef std::function<Eigen::Matrix3d(const double &, const double &, const Eigen::Vector3d &, const Cell *)> GradientDisplacementType1;

    /// Constructor
    Elasticity(
               const DDRCore & ddrcore,           ///< Core for the DDR space sequence
               const ConnectedInformationGenerator & connectedinfo,
               const double & mu,                 ///< First Lame coefficient
               const double & lambda,             ///< Second Lame coefficient
               const double & beta,               ///< Newton parameter
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
                
    /// Initialising a zero vector for Newton's method
    // Eigen::VectorXd vectorInitialGuess;
    // Eigen::VectorXd vectorDiscreteSolved;
    
    void assembleLinearSystem(
                 const ForcingTermType & f,      ///< Forcing term
                 const FrictionType & g,      ///< Friction
                 const DisplacementType & u,     ///< u for natural BCs
                 const GradientDisplacementType1 & sigma_u, ///< sigma_u for natural BCs
                 const GradientDisplacementType & symgrad_u, ///< Symmetric gradient of u for natural BCs
                 const Eigen::VectorXd & UDir,        ///< Vector filled in by Dirichlets BCs (zero elsewhere) for essential BCs
                 const Eigen::VectorXd & U_Iguess
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
      return m_fractureinfo.dimension();
    }

    /// Number of DOFs for fracture faces
    inline size_t fDOFs() const
    {
      return m_fractureinfo.fracture_dofs();
    }

    /// Number of DOFs for displacement
    inline size_t tDOFs() const
    {
      return m_fractureinfo.dimension()+m_fractureinfo.fracture_dofs();
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
      // return m_fractureinfo.dimension() - nDOFs_dir();
      return m_fractureinfo.dimension() + m_fractureinfo.fracture_dofs() - nDOFs_dir();
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

    //--------------------------------------------------------
    //  Getters: underlying Fracture DOF Information
    //--------------------------------------------------------
    
    /// Returns the space fracturedof
    inline const FractureDOFSpace & fractureinfo() const
    {
      return m_fractureinfo;
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
                      const FrictionType & g,  ///< Friction
                      const DisplacementType & u, ///< Displacement term
                      const GradientDisplacementType1 & sigma_u, ///< sigma_u for natural BCs
                      const GradientDisplacementType & symgrad_u, ///< Symmetric gradient of u for natural BCs
                      const Eigen::VectorXd & UDir, ///< Boundary condition
                      const Eigen::VectorXd & U_Iguess
                      ) const ;

    /// Creates the permutation matrix and the global DOFs for the local static condensation
    LocalStaticCondensation _compute_static_condensation(const size_t & iT) const;

    /// Assemble the local contribution for the element of index iT into the global system
    void _assemble_local_contribution(
                                      size_t iT,                                               ///< Element index
                                      const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT, ///< Local contribution
                                      std::list<Eigen::Triplet<double> > & A1,                 ///< List of triplets for system
                                      Eigen::VectorXd & b1,                                    ///< Vector for RHS for sysem
                                      // std::list<Eigen::Triplet<double> > & A2,                  ///< List of triplets for BC
                                      std::list<Eigen::Triplet<double> > & A2,                 ///< List of triplets for sc
                                      Eigen::VectorXd & b2                                     ///< Vector for RHS for sc
                                      );
    
    
    // Constructor parameter
    const DDRCore & m_ddrcore;
    const ConnectedInformationGenerator & m_connectedinfo;
    bool m_use_threads;
    std::ostream & m_output;
    SerendipityProblem m_ser_pro;
    const VSXGrad m_vsxgrad;
    const FractureInterpolator m_fracture_interpolator;
    const FractureDOFSpace m_fractureinfo;
    const Eigen::VectorXi m_nloc_sc; // Nb of statically condensed DOFs for velocity in each cell (cell unknowns)
    
    // Physical parameters
    double m_mu;
    double m_lambda;
    
    // Newton parameter
    double m_beta;

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
    // Eigen::MatrixXd m_fA;   
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
  trigonometric_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                             sin(2*PI*x.x())*cos(2*PI*x.y())*cos(2*PI*x.z()),
                                             0.5*sin(2*PI*x.y())*cos(2*PI*x.x())*cos(2*PI*x.z()),
                                             -sin(2*PI*x.z())*cos(2*PI*x.x())*cos(2*PI*x.y())
                                             );
                    };

  static Elasticity::GradientDisplacementType
  trigonometric_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
                        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << 2*PI*cos(2*PI*x.x())*cos(2*PI*x.y())*cos(2*PI*x.z()),
                                   -2*PI*sin(2*PI*x.x())*sin(2*PI*x.y())*cos(2*PI*x.z()), -2*PI*sin(2*PI*x.x())*sin(2*PI*x.z())*cos(2*PI*x.y());
                        M.row(1) << -1.0*PI*sin(2*PI*x.x())*sin(2*PI*x.y())*cos(2*PI*x.z()),
                                     1.0*PI*cos(2*PI*x.x())*cos(2*PI*x.y())*cos(2*PI*x.z()), -1.0*PI*sin(2*PI*x.y())*sin(2*PI*x.z())*cos(2*PI*x.x());
                        M.row(2) << 2*PI*sin(2*PI*x.x())*sin(2*PI*x.z())*cos(2*PI*x.y()), 
                                    2*PI*sin(2*PI*x.y())*sin(2*PI*x.z())*cos(2*PI*x.x()), -2*PI*cos(2*PI*x.x())*cos(2*PI*x.y())*cos(2*PI*x.z());
                        return M;
                         };
  
  static Elasticity::GradientDisplacementType1
  trigonometric_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {
        
        Eigen::Matrix3d grad = Eigen::Matrix3d::Zero();

        grad.row(0) <<  2*PI*cos(2*PI*x.x())*cos(2*PI*x.y())*cos(2*PI*x.z()),
                        -PI*sin(2*PI*x.x())*sin(2*PI*x.y())*cos(2*PI*x.z())-0.5*PI*sin(2*PI*x.x())*sin(2*PI*x.y())*cos(2*PI*x.z()), 
                        -PI*sin(2*PI*x.x())*sin(2*PI*x.z())*cos(2*PI*x.y())+PI*sin(2*PI*x.x())*sin(2*PI*x.z())*cos(2*PI*x.y());
        
        grad.row(1) <<  -PI*sin(2*PI*x.x())*sin(2*PI*x.y())*cos(2*PI*x.z())-0.5*PI*sin(2*PI*x.x())*sin(2*PI*x.y())*cos(2*PI*x.z()), 
                        1.0*PI*cos(2*PI*x.x())*cos(2*PI*x.y())*cos(2*PI*x.z()),            
                        -0.5*PI*sin(2*PI*x.y())*sin(2*PI*x.z())*cos(2*PI*x.x())+PI*sin(2*PI*x.y())*sin(2*PI*x.z())*cos(2*PI*x.x());
                  
        grad.row(2) <<  -PI*sin(2*PI*x.x())*sin(2*PI*x.z())*cos(2*PI*x.y())+PI*sin(2*PI*x.x())*sin(2*PI*x.z())*cos(2*PI*x.y()),
                        -0.5*PI*sin(2*PI*x.y())*sin(2*PI*x.z())*cos(2*PI*x.x())+PI*sin(2*PI*x.y())*sin(2*PI*x.z())*cos(2*PI*x.x()),
                        -2*PI*cos(2*PI*x.x())*cos(2*PI*x.y())*cos(2*PI*x.z());

        // All other entries remain zero since other components are 0

        Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
        trace_grad(0,0)=grad(0,0)+grad(1,1)+grad(2,2);  
        trace_grad(1,1)=grad(0,0)+grad(1,1)+grad(2,2);  
        trace_grad(2,2)=grad(0,0)+grad(1,1)+grad(2,2);  
     
        return 2*mu*grad + lambda*trace_grad;
  };                       

  static Elasticity::ForcingTermType
  trigonometric_div_eps_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                -7.0*pow(PI, 2)*sin(2*PI*x.x())*cos(2*PI*x.y())*cos(2*PI*x.z()),
                                -4.0*pow(PI, 2)*sin(2*PI*x.y())*cos(2*PI*x.x())*cos(2*PI*x.z()),
                                5.0*pow(PI, 2)*sin(2*PI*x.z())*cos(2*PI*x.x())*cos(2*PI*x.y())
                                ); 
                    };

  static Elasticity::ForcingTermType
  trigonometric_grad_div_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                -2.0*pow(PI, 2)*sin(2*PI*x.x())*cos(2*PI*x.y())*cos(2*PI*x.z()),
                                -2.0*pow(PI, 2)*sin(2*PI*x.y())*cos(2*PI*x.x())*cos(2*PI*x.z()),
                                -2.0*pow(PI, 2)*sin(2*PI*x.z())*cos(2*PI*x.x())*cos(2*PI*x.y())
                                );  
                    };

  static Elasticity::ForcingTermType
  trigonometric_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                    return -2 * mu * trigonometric_div_eps_u(mu,lambda,x,T) - lambda * trigonometric_grad_div_u(mu,lambda,x,T); 
                    };

  //------------------------------------------------------------------------------
  // Fracture solution
  //------------------------------------------------------------------------------

  static Elasticity::DisplacementType
  fracture_test_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
    if (T->center_mass().x() > 0.5 &&
        T->center_mass().y() >= 0.0 && T->center_mass().y() <= 0.5 &&
        T->center_mass().z() >= 0.5 && T->center_mass().z() <= 1.0) {
        
        double val = std::pow(x.x(), 2)*std::pow(x.y() - 0.5, 2)*std::pow(x.z() - 0.5, 2);
        return Eigen::Vector3d(val, 0.0, 0.0);
    } 
    else {
        return Eigen::Vector3d::Zero();
    }
  };

  static Elasticity::GradientDisplacementType
  fracture_test_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
    if (T->center_mass().x() > 0.5 &&
        T->center_mass().y() >= 0.0 && T->center_mass().y() <= 0.5 &&
        T->center_mass().z() >= 0.5 && T->center_mass().z() <= 1.0) {
        
        Eigen::Matrix3d grad = Eigen::Matrix3d::Zero();

        grad(0, 0) = 2.0 * x.x() * std::pow(x.y() - 0.5, 2) * std::pow(x.z() - 0.5, 2);         // ∂f₁/∂x
        grad(0, 1) = 2.0 * std::pow(x.x(), 2) * (x.y() - 0.5) * std::pow(x.z() - 0.5, 2);         // ∂f₁/∂y
        grad(0, 2) = 2.0 * std::pow(x.x(), 2) * std::pow(x.y() - 0.5, 2) * (x.z() - 0.5);         // ∂f₁/∂z

        // All other entries remain zero since other components are 0
        
        return grad;
    } 
    else {
        return Eigen::Matrix3d::Zero();
    }
  };

  static Elasticity::GradientDisplacementType1
  fracture_test_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {
    if (T->center_mass().x() > 0.5 &&
        T->center_mass().y() >= 0.0 && T->center_mass().y() <= 0.5 &&
        T->center_mass().z() >= 0.5 && T->center_mass().z() <= 1.0) {
        
        Eigen::Matrix3d grad = Eigen::Matrix3d::Zero();

        grad(0, 0) = 2.0 * x.x() * std::pow(x.y() - 0.5, 2) * std::pow(x.z() - 0.5, 2);  
        grad(0, 1) = std::pow(x.x(), 2) * (x.y() - 0.5) * std::pow(x.z() - 0.5, 2);         
        grad(0, 2) = std::pow(x.x(), 2) * std::pow(x.y() - 0.5, 2) * (x.z() - 0.5);         
        grad(1, 0) = std::pow(x.x(), 2) * (x.y() - 0.5) * std::pow(x.z() - 0.5, 2);         
        grad(2, 0) = std::pow(x.x(), 2) * std::pow(x.y() - 0.5, 2) * (x.z() - 0.5);         

        // All other entries remain zero since other components are 0

        Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
        trace_grad(0,0)=2.0 * x.x() * std::pow(x.y() - 0.5, 2) * std::pow(x.z() - 0.5, 2);  
        trace_grad(1,1)=2.0 * x.x() * std::pow(x.y() - 0.5, 2) * std::pow(x.z() - 0.5, 2);  
        trace_grad(2,2)=2.0 * x.x() * std::pow(x.y() - 0.5, 2) * std::pow(x.z() - 0.5, 2);  

        
        return 2*mu*grad + lambda*trace_grad;
    } 
    else {
        return Eigen::Matrix3d::Zero();
    }
  };

  static Elasticity::ForcingTermType
  fracture_test_div_eps_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
    if (T->center_mass().x() > 0.5 && T->center_mass().y() >= 0 && T->center_mass().y() <= 0.5 && T->center_mass().z() >= 0.5 && T->center_mass().z() <= 1.0) {
        double comp1 = 2*std::pow(x.y() - 0.5, 2) * std::pow(x.z() - 0.5, 2)
                     + x.x() * x.x() * std::pow(x.z() - 0.5, 2)
                     + x.x() * x.x() * std::pow(x.y() - 0.5, 2);
        double comp2 = 2 * x.x() * (x.y() - 0.5) * std::pow(x.z() - 0.5, 2);
        double comp3 = 2 * x.x() * std::pow(x.y() - 0.5, 2) * (x.z() - 0.5);
        return Eigen::Vector3d(comp1, comp2, comp3);
    } else {
        return Eigen::Vector3d(0.0, 0.0, 0.0);
    }
};

  static Elasticity::ForcingTermType
  fracture_test_grad_div_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
    if (T->center_mass().x() > 0.5 && T->center_mass().y() >= 0 && T->center_mass().y() <= 0.5 && T->center_mass().z() >= 0.5 && T->center_mass().z() <= 1.0) {
        return Eigen::Vector3d{
            2 * std::pow(x.y() - 0.5, 2) * std::pow(x.z() - 0.5, 2),
            4 * x.x() * (x.y() - 0.5) * std::pow(x.z() - 0.5, 2),
            4 * x.x() * std::pow(x.y() - 0.5, 2) * (x.z() - 0.5)
        };
    } else {
        return Eigen::Vector3d::Zero();
    }
};

  static Elasticity::ForcingTermType
  fracture_test_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                    return -2 * mu * fracture_test_div_eps_u(mu,lambda,x,T) - lambda * fracture_test_grad_div_u(mu,lambda,x,T); 
                    };

  //------------------------------------------------------------------------------
  // continuous trignometric solution with 0 traction and tangential traction component
  //------------------------------------------------------------------------------

  static Elasticity::DisplacementType
  trigonometric_tangential_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                             -sin(PI*x.x()),
                                             (1/(2*PI))*sin(2*PI*x.y())*cos(2*PI*x.z()),
                                             -(1/(2*PI))*cos(2*PI*x.y())*sin(2*PI*x.z()));
                    };   
  
  static Elasticity::GradientDisplacementType
  trigonometric_tangential_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
                        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << -PI*cos(PI*x.x()), 0, 0;
                        M.row(1) << 0, cos(2*PI*x.y())*cos(2*PI*x.z()), -sin(2*PI*x.y())*sin(2*PI*x.z());
                        M.row(2) << 0, sin(2*PI*x.y())*sin(2*PI*x.z()), -cos(2*PI*x.y())*cos(2*PI*x.z());
                        return M;
                         };         
                         
  static Elasticity::GradientDisplacementType1
  trigonometric_tangential_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {
        
        Eigen::Matrix3d grad = Eigen::Matrix3d::Zero();

        // grad.row(0) <<  -PI*cos(PI*x.x()), 0, 0;
        grad.row(0) <<  -PI*cos(PI*x.x()), 0, 0;
        
        grad.row(1) <<  0, cos(2*PI*x.y())*cos(2*PI*x.z()), 0;
                  
        grad.row(2) <<  0, 0, -cos(2*PI*x.y())*cos(2*PI*x.z());

        // All other entries remain zero since other components are 0

        Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
        trace_grad(0,0)=grad(0,0)+grad(1,1)+grad(2,2);  
        trace_grad(1,1)=grad(0,0)+grad(1,1)+grad(2,2);  
        trace_grad(2,2)=grad(0,0)+grad(1,1)+grad(2,2);  
     
        return 2*mu*grad + lambda*trace_grad;
  };                       


  static Elasticity::ForcingTermType
  trigonometric_tangential_div_eps_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                      return Eigen::Vector3d(std::pow(PI,2)*sin(PI*x.x()),
                              -(2*PI)*sin(2*PI*x.y())*cos(2*PI*x.z()),
                              (2*PI)*cos(2*PI*x.y())*sin(2*PI*x.z()));
                    };

  static Elasticity::ForcingTermType
  trigonometric_tangential_grad_div_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                      return Eigen::Vector3d(std::pow(PI,2)*sin(PI*x.x()),0,0);  
                    };

  static Elasticity::ForcingTermType
  trigonometric_tangential_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                    return -2 * mu * trigonometric_tangential_div_eps_u(mu,lambda,x,T) - lambda * trigonometric_tangential_grad_div_u(mu,lambda,x,T); 
                    };                     
  
  //------------------------------------------------------------------------------
  // continuous trignometric solution with 0 tangential component and negative traction normal
  //------------------------------------------------------------------------------

  static Elasticity::DisplacementType
  trigonometric_new_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                             -sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()))*pow(x.z(),2),
                                             pow(x.z(),2),
                                             pow(x.x()-0.5,2)*pow(x.z(),2));
                    };   
  
  static Elasticity::GradientDisplacementType
  trigonometric_new_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
                        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << -(PI/2)*cos((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()))*pow(x.z(),2),
                        (PI/2)*sin((PI/2)*(x.x()-0.5))*sin((PI/2)*(x.y()))*pow(x.z(),2),-2*sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()))*x.z();
                        M.row(1) << 0, 0, 2*x.z();
                        M.row(2) << 2*(x.x()-0.5)*pow(x.z(),2), 0, 2*pow(x.x()-0.5,2)*x.z();
                        return M;
                         };         
                         
  static Elasticity::GradientDisplacementType1
  trigonometric_new_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {
        
        Eigen::Matrix3d grad = Eigen::Matrix3d::Zero();

        grad.row(0) <<  -(PI/2)*cos((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()))*pow(x.z(),2),
                        (PI/4)*sin((PI/2)*(x.x()-0.5))*sin((PI/2)*(x.y()))*pow(x.z(),2),
                        -sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()))*x.z()+(x.x()-0.5)*pow(x.z(),2);
        
        grad.row(1) <<  (PI/4)*sin((PI/2)*(x.x()-0.5))*sin((PI/2)*(x.y()))*pow(x.z(),2),
                        0,
                        x.z();
                  
        grad.row(2) <<  -sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()))*x.z()+(x.x()-0.5)*pow(x.z(),2),
                        x.z(),
                        2*pow(x.x()-0.5,2)*x.z();

        // All other entries remain zero since other components are 0

        Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
        trace_grad(0,0)=grad(0,0)+grad(1,1)+grad(2,2);  
        trace_grad(1,1)=grad(0,0)+grad(1,1)+grad(2,2);  
        trace_grad(2,2)=grad(0,0)+grad(1,1)+grad(2,2);  
     
        return 2*mu*grad + lambda*trace_grad;
  };                       


  static Elasticity::ForcingTermType
  trigonometric_new_div_eps_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                (0.25)*std::pow(PI,2)*sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()))*pow(x.z(),2)+(0.125)*(std::pow(PI,2))*sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()))*pow(x.z(),2)-sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()))+2*(x.x()-0.5)*x.z(),
                                (0.125)*(std::pow(PI,2))*cos((PI/2)*(x.x()-0.5))*sin((PI/2)*(x.y()))*pow(x.z(),2)+1,
                                -(PI/2)*cos((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()))*x.z()+std::pow(x.z(),2)+2*std::pow(x.x()-0.5,2)
                                ); 
                    };

  static Elasticity::ForcingTermType
  trigonometric_new_grad_div_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                (0.25)*std::pow(PI,2)*sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()))*pow(x.z(),2)+4*(x.x()-0.5)*x.z(),
                                (0.25)*std::pow(PI,2)*cos((PI/2)*(x.x()-0.5))*sin((PI/2)*(x.y()))*pow(x.z(),2),
                                -PI*cos((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()))*x.z()+2*std::pow(x.x()-0.5,2)
                                );  
                    };

  static Elasticity::ForcingTermType
  trigonometric_new_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                    return -2 * mu * trigonometric_new_div_eps_u(mu,lambda,x,T) - lambda * trigonometric_new_grad_div_u(mu,lambda,x,T); 
                    };  
                    
  //------------------------------------------------------------------------------
  // fracture trignometric solution with 0 traction
  //------------------------------------------------------------------------------

  static Elasticity::DisplacementType
  fracture_trigonometric_test_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
    Eigen::Vector3d u_base(cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),4),
          -(8/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),4),
          -(8/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),3));

    if (T->center_mass().x() > 0.5 &&
        T->center_mass().y() >= 0.0 && T->center_mass().y() <= 0.5 &&
        T->center_mass().z() >= 0.5 && T->center_mass().z() <= 1.0) {
        return 2*u_base;
    } 
    else {
        return u_base;
    }
  }; 

  static Elasticity::GradientDisplacementType
  fracture_trigonometric_test_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
    Eigen::Matrix3d gradu_base;
    gradu_base(0,0) = -(PI/2)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),4);
    gradu_base(0,1) = 4*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),4); 
    gradu_base(0,2) = 4*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),3);
    gradu_base(1,0) = -4*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),4);
    gradu_base(1,1) = -(24/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),2)*std::pow((x.z()-0.5),4);
    gradu_base(1,2) = -(32/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),3);
    gradu_base(2,0) = -4*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),3);
    gradu_base(2,1) = -(32/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),3);
    gradu_base(2,2) = -(24/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),2);

    if (T->center_mass().x() > 0.5 &&
        T->center_mass().y() >= 0.0 && T->center_mass().y() <= 0.5 &&
        T->center_mass().z() >= 0.5 && T->center_mass().z() <= 1.0) {
        
        return 2*gradu_base;
    } 
    else {
        return gradu_base;
    }
  };       

  static Elasticity::GradientDisplacementType1
  fracture_trigonometric_test_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {
    Eigen::Matrix3d epsu_base;
    epsu_base(0,0) = -(PI/2)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),4);
    epsu_base(0,1) = 0; 
    epsu_base(0,2) = 0;
    epsu_base(1,0) = 0;
    epsu_base(1,1) = -(24/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),2)*std::pow((x.z()-0.5),4);
    epsu_base(1,2) = -(32/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),3);
    epsu_base(2,0) = 0;
    epsu_base(2,1) = -(32/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),3);
    epsu_base(2,2) = -(24/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),2);

    Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
        trace_grad(0,0)=epsu_base(0,0)+epsu_base(1,1)+epsu_base(2,2);  
        trace_grad(1,1)=epsu_base(0,0)+epsu_base(1,1)+epsu_base(2,2);  
        trace_grad(2,2)=epsu_base(0,0)+epsu_base(1,1)+epsu_base(2,2);  

    Eigen::Matrix3d sigmau_base = 2*mu*epsu_base + lambda*trace_grad;

    if (T->center_mass().x() > 0.5 &&
        T->center_mass().y() >= 0.0 && T->center_mass().y() <= 0.5 &&
        T->center_mass().z() >= 0.5 && T->center_mass().z() <= 1.0) {

        return 2*sigmau_base;
    } 
    else {
        return sigmau_base;
    }
  };                                            

  static Elasticity::ForcingTermType
  fracture_trigonometric_test_div_eps_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
    Eigen::Vector3d div_epsu_base(-(0.25)*std::pow(PI,2)*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),4),
          -(48/PI)*sin((PI/2)*(x.x()-0.5))*(x.y()-0.5)*std::pow((x.z()-0.5),4)-(96/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),2),
          -(96/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),2)*std::pow((x.z()-0.5),3)-(48/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*(x.z()-0.5));

    if (T->center_mass().x() > 0.5 && T->center_mass().y() >= 0 && T->center_mass().y() <= 0.5 && T->center_mass().z() >= 0.5 && T->center_mass().z() <= 1.0) {
        return 2*div_epsu_base;
    } else {
        return div_epsu_base;
    }
};

  static Elasticity::ForcingTermType
  fracture_trigonometric_test_grad_div_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
    Eigen::Vector3d grad_divu_base(-(0.25)*std::pow(PI,2)*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),4)-(12)*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),2)*std::pow((x.z()-0.5),4)-(12)*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),2),
            -(2)*PI*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),4)-(48/PI)*sin((PI/2)*(x.x()-0.5))*(x.y()-0.5)*std::pow((x.z()-0.5),4)-(96/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),2),
          -(2)*PI*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),3)-(96/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),2)*std::pow((x.z()-0.5),3)-(48/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*(x.z()-0.5)
          );
    if (T->center_mass().x() > 0.5 && T->center_mass().y() >= 0 && T->center_mass().y() <= 0.5 && T->center_mass().z() >= 0.5 && T->center_mass().z() <= 1.0) {
        return 2*grad_divu_base;
    } else {
        return grad_divu_base;
    }
};

  static Elasticity::ForcingTermType
  fracture_trigonometric_test_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                    return -2 * mu * fracture_trigonometric_test_div_eps_u(mu,lambda,x,T) - lambda * fracture_trigonometric_test_grad_div_u(mu,lambda,x,T); 
                    };   
                    
  //------------------------------------------------------------------------------
  // frictionless static contact problem
  //------------------------------------------------------------------------------

  static Elasticity::DisplacementType
  frictionless_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
    Eigen::Vector3d u_base(cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),4),
          -(8/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),4),
          -(8/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),3));

    if (T->center_mass().x() > 0.5 &&
        T->center_mass().y() >= 0.0 && T->center_mass().y() <= 0.5 &&
        T->center_mass().z() >= 0.5 && T->center_mass().z() <= 1.0) {
        return 2*u_base;
    } 
    else if (T->center_mass().z() < 0.5){
      return Eigen::Vector3d(
                            -sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()-0.5))*pow(x.z()-0.5,2),
                            pow(x.z()-0.5,2),
                            pow(x.x()-0.5,2)*pow(x.z()-0.5,2));
    }
    else {
        return u_base;
    }
  }; 

  static Elasticity::FrictionType
  frictionless_friction = [](const Face * F, const double & lambda) -> double{
    return 0.0;
  };

  static Elasticity::GradientDisplacementType
  frictionless_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
    Eigen::Matrix3d gradu_base;
    gradu_base(0,0) = -(PI/2)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),4);
    gradu_base(0,1) = 4*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),4); 
    gradu_base(0,2) = 4*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),3);
    gradu_base(1,0) = -4*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),4);
    gradu_base(1,1) = -(24/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),2)*std::pow((x.z()-0.5),4);
    gradu_base(1,2) = -(32/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),3);
    gradu_base(2,0) = -4*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),3);
    gradu_base(2,1) = -(32/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),3);
    gradu_base(2,2) = -(24/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),2);

    if (T->center_mass().x() > 0.5 &&
        T->center_mass().y() >= 0.0 && T->center_mass().y() <= 0.5 &&
        T->center_mass().z() >= 0.5 && T->center_mass().z() <= 1.0) {
        
        return 2*gradu_base;
    } 
    else if (T->center_mass().z() < 0.5) {
      Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << -(PI/2)*cos((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()-0.5))*pow(x.z()-0.5,2),
                        (PI/2)*sin((PI/2)*(x.x()-0.5))*sin((PI/2)*(x.y()-0.5))*pow(x.z()-0.5,2),-2*sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()-0.5))*(x.z()-0.5);
                        M.row(1) << 0, 0, 2*(x.z()-0.5);
                        M.row(2) << 2*(x.x()-0.5)*pow(x.z()-0.5,2), 0, 2*pow(x.x()-0.5,2)*(x.z()-0.5);
                        return M;
    }
    else {
        return gradu_base;
    }
  };       

  static Elasticity::GradientDisplacementType1
  frictionless_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {
    Eigen::Matrix3d epsu_base;
    epsu_base(0,0) = -(PI/2)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),4);
    epsu_base(0,1) = 0; 
    epsu_base(0,2) = 0;
    epsu_base(1,0) = 0;
    epsu_base(1,1) = -(24/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),2)*std::pow((x.z()-0.5),4);
    epsu_base(1,2) = -(32/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),3);
    epsu_base(2,0) = 0;
    epsu_base(2,1) = -(32/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),3);
    epsu_base(2,2) = -(24/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),2);

    Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
        trace_grad(0,0)=epsu_base(0,0)+epsu_base(1,1)+epsu_base(2,2);  
        trace_grad(1,1)=epsu_base(0,0)+epsu_base(1,1)+epsu_base(2,2);  
        trace_grad(2,2)=epsu_base(0,0)+epsu_base(1,1)+epsu_base(2,2);  

    Eigen::Matrix3d sigmau_base = 2*mu*epsu_base + lambda*trace_grad;

    if (T->center_mass().x() > 0.5 &&
        T->center_mass().y() >= 0.0 && T->center_mass().y() <= 0.5 &&
        T->center_mass().z() >= 0.5 && T->center_mass().z() <= 1.0) {

        return 2*sigmau_base;
    } 
    else if (T->center_mass().z() < 0.5){
        Eigen::Matrix3d grad = Eigen::Matrix3d::Zero();

        grad.row(0) <<  -(PI/2)*cos((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()-0.5))*pow(x.z()-0.5,2),
                        (PI/4)*sin((PI/2)*(x.x()-0.5))*sin((PI/2)*(x.y()-0.5))*pow(x.z()-0.5,2),
                        -sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()-0.5))*(x.z()-0.5)+(x.x()-0.5)*pow(x.z()-0.5,2);
        
        grad.row(1) <<  (PI/4)*sin((PI/2)*(x.x()-0.5))*sin((PI/2)*(x.y()-0.5))*pow(x.z()-0.5,2),
                        0,
                        (x.z()-0.5);
                  
        grad.row(2) <<  -sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()-0.5))*(x.z()-0.5)+(x.x()-0.5)*pow(x.z()-0.5,2),
                        (x.z()-0.5),
                        2*pow(x.x()-0.5,2)*(x.z()-0.5);

        // All other entries remain zero since other components are 0

        Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
        trace_grad(0,0)=grad(0,0)+grad(1,1)+grad(2,2);  
        trace_grad(1,1)=grad(0,0)+grad(1,1)+grad(2,2);  
        trace_grad(2,2)=grad(0,0)+grad(1,1)+grad(2,2);  
     
        return 2*mu*grad + lambda*trace_grad;
    }
    else {
        return sigmau_base;
    }
  };                                            

  static Elasticity::ForcingTermType
  frictionless_div_eps_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
    Eigen::Vector3d div_epsu_base(-(0.25)*std::pow(PI,2)*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),4),
          -(48/PI)*sin((PI/2)*(x.x()-0.5))*(x.y()-0.5)*std::pow((x.z()-0.5),4)-(96/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),2),
          -(96/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),2)*std::pow((x.z()-0.5),3)-(48/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*(x.z()-0.5));

    if (T->center_mass().x() > 0.5 && T->center_mass().y() >= 0 && T->center_mass().y() <= 0.5 && T->center_mass().z() >= 0.5 && T->center_mass().z() <= 1.0) {
        return 2*div_epsu_base;
    } 
    else if (T->center_mass().z() < 0.5){
      return Eigen::Vector3d(
                                (0.25)*std::pow(PI,2)*sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()-0.5))*pow(x.z()-0.5,2)+(0.125)*(std::pow(PI,2))*sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()-0.5))*pow(x.z()-0.5,2)-sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()-0.5))+2*(x.x()-0.5)*(x.z()-0.5),
                                (0.125)*(std::pow(PI,2))*cos((PI/2)*(x.x()-0.5))*sin((PI/2)*(x.y()-0.5))*pow(x.z()-0.5,2)+1,
                                -(PI/2)*cos((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()-0.5))*(x.z()-0.5)+std::pow(x.z()-0.5,2)+2*std::pow(x.x()-0.5,2)
                                ); 
    }
    else {
        return div_epsu_base;
    }
};

  static Elasticity::ForcingTermType
  frictionless_grad_div_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
    Eigen::Vector3d grad_divu_base(-(0.25)*std::pow(PI,2)*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),4)-(12)*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),2)*std::pow((x.z()-0.5),4)-(12)*cos((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),2),
            -(2)*PI*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),4)-(48/PI)*sin((PI/2)*(x.x()-0.5))*(x.y()-0.5)*std::pow((x.z()-0.5),4)-(96/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),3)*std::pow((x.z()-0.5),2),
          -(2)*PI*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*std::pow((x.z()-0.5),3)-(96/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),2)*std::pow((x.z()-0.5),3)-(48/PI)*sin((PI/2)*(x.x()-0.5))*std::pow((x.y()-0.5),4)*(x.z()-0.5)
          );
    if (T->center_mass().x() > 0.5 && T->center_mass().y() >= 0 && T->center_mass().y() <= 0.5 && T->center_mass().z() >= 0.5 && T->center_mass().z() <= 1.0) {
        return 2*grad_divu_base;
    } 
    else if (T->center_mass().z() < 0.5) {
      return Eigen::Vector3d(
                                (0.25)*std::pow(PI,2)*sin((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()-0.5))*pow(x.z()-0.5,2)+4*(x.x()-0.5)*(x.z()-0.5),
                                (0.25)*std::pow(PI,2)*cos((PI/2)*(x.x()-0.5))*sin((PI/2)*(x.y()-0.5))*pow(x.z()-0.5,2),
                                -PI*cos((PI/2)*(x.x()-0.5))*cos((PI/2)*(x.y()-0.5))*(x.z()-0.5)+2*std::pow(x.x()-0.5,2)
                                ); 
    }
    else {
        return grad_divu_base;
    }
};

  static Elasticity::ForcingTermType
  frictionless_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                    return -2 * mu * frictionless_div_eps_u(mu,lambda,x,T) - lambda * frictionless_grad_div_u(mu,lambda,x,T); 
                    };

  //------------------------------------------------------------------------------
  // frictionless static contact problem
  //------------------------------------------------------------------------------

  static Elasticity::DisplacementType
  frictionless_new_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
    Eigen::Vector3d u_base(cos((PI/2)*x.x())*std::pow(x.z(),4),
          4*cos((PI/2)*x.x())*std::pow(x.z(),3),
          -(8/PI)*sin((PI/2)*x.x())*std::pow(x.z(),3));

    if (T->center_mass().x() < 0 &&
        T->center_mass().z() < 0) {
        return u_base;
    } 
    else if (T->center_mass().z() >= 0){
      return Eigen::Vector3d(
                            -sin((PI/2)*x.x())*cos((PI/2)*x.y())*pow(x.z(),2),
                            pow(x.z(),2),
                            pow(x.x(),2)*pow(x.z(),2));
    }
    else {
        return 2*u_base;
    }
  }; 

  static Elasticity::FrictionType
  frictionless_new_friction = [](const Face * F, const double & lambda) -> double{
    return 0.0;
  };

  static Elasticity::GradientDisplacementType
  frictionless_new_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
    Eigen::Matrix3d gradu_base;
    gradu_base(0,0) = -0.5*PI*pow(x.z(), 4)*sin(0.5*PI*x.x());
    gradu_base(0,1) = 0; 
    gradu_base(0,2) = 4*pow(x.z(), 3)*cos(0.5*PI*x.x());
    gradu_base(1,0) = -2.0*PI*pow(x.z(), 3)*sin(0.5*PI*x.x());
    gradu_base(1,1) = 0;
    gradu_base(1,2) = 12*pow(x.z(), 2)*cos(0.5*PI*x.x());
    gradu_base(2,0) = -4.0*pow(x.z(), 3)*cos(0.5*PI*x.x());
    gradu_base(2,1) = 0;
    gradu_base(2,2) = -24*pow(x.z(), 2)*sin(0.5*PI*x.x())/PI;

    if (T->center_mass().x() < 0 &&
        T->center_mass().z() < 0) {
        
        return gradu_base;
    } 
    else if (T->center_mass().z() >= 0) {
      Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << -0.5*PI*pow(x.z(), 2)*cos(0.5*PI*x.x())*cos(0.5*PI*x.y()),
                         0.5*PI*pow(x.z(), 2)*sin(0.5*PI*x.x())*sin(0.5*PI*x.y()),-2*x.z()*sin(0.5*PI*x.x())*cos(0.5*PI*x.y());
                        M.row(1) << 0, 0, 2*x.z();
                        M.row(2) << 2*x.x()*pow(x.z(), 2), 0, 2*pow(x.x(), 2)*x.z();
                        return M;
    }
    else {
        return 2*gradu_base;
    }
  };       

  static Elasticity::GradientDisplacementType1
  frictionless_new_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {
    Eigen::Matrix3d grad = frictionless_new_grad_u(x,T,lambda);
    Eigen::Matrix3d eps = (0.5)*(grad+grad.transpose());
    Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
    trace_grad(0,0)=grad(0,0)+grad(1,1)+grad(2,2);  
    trace_grad(1,1)=grad(0,0)+grad(1,1)+grad(2,2);  
    trace_grad(2,2)=grad(0,0)+grad(1,1)+grad(2,2);  
    return 2*mu*eps + lambda*trace_grad;
  };

                                            

  static Elasticity::ForcingTermType
  frictionless_new_div_eps_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
    Eigen::Vector3d div_epsu_base(-0.25*pow(PI, 2)*pow(x.z(), 4)*cos(0.5*PI*x.x()),
     -0.5*pow(PI, 2)*pow(x.z(), 3)*cos(0.5*PI*x.x()) + 12.0*x.z()*cos(0.5*PI*x.x()),
     -48.0*x.z()*sin(0.5*PI*x.x())/PI);

    if (T->center_mass().x() < 0 && T->center_mass().z() < 0) {
        return div_epsu_base;
    } 
    else if (T->center_mass().z() >= 0){
      return Eigen::Vector3d(
                                0.375*pow(PI, 2)*pow(x.z(), 2)*sin(0.5*PI*x.x())*cos(0.5*PI*x.y()) + 2.0*x.x()*x.z() - 1.0*sin(0.5*PI*x.x())*cos(0.5*PI*x.y()),
                                 0.125*pow(PI, 2)*pow(x.z(), 2)*sin(0.5*PI*x.y())*cos(0.5*PI*x.x()) + 1.0,
                                 -0.5*PI*x.z()*cos(0.5*PI*x.x())*cos(0.5*PI*x.y()) + 2.0*pow(x.x(), 2) + 1.0*pow(x.z(), 2)
                                ); 
    }
    else {
        return 2*div_epsu_base;
    }
};

  static Elasticity::ForcingTermType
  frictionless_new_grad_div_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
    Eigen::Vector3d grad_divu_base(-0.25*pow(PI, 2)*pow(x.z(), 4)*cos(0.5*PI*x.x()) - 12.0*pow(x.z(), 2)*cos(0.5*PI*x.x()),
     0, -2.0*PI*pow(x.z(), 3)*sin(0.5*PI*x.x()) - 48*x.z()*sin(0.5*PI*x.x())/PI);
    if (T->center_mass().x() < 0 && T->center_mass().z() < 0) {
        return grad_divu_base;
    } 
    else if (T->center_mass().z() >= 0) {
      return Eigen::Vector3d(
                                0.25*pow(PI, 2)*pow(x.z(), 2)*sin(0.5*PI*x.x())*cos(0.5*PI*x.y()) + 4*x.x()*x.z(),
                                 0.25*pow(PI, 2)*pow(x.z(), 2)*sin(0.5*PI*x.y())*cos(0.5*PI*x.x()),
                                 -1.0*PI*x.z()*cos(0.5*PI*x.x())*cos(0.5*PI*x.y()) + 2*pow(x.x(), 2)
                                ); 
    }
    else {
        return 2*grad_divu_base;
    }
};

  static Elasticity::ForcingTermType
  frictionless_new_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                    return -2 * mu * frictionless_new_div_eps_u(mu,lambda,x,T) - lambda * frictionless_new_grad_div_u(mu,lambda,x,T); 
                    };

//------------------------------------------------------------------------------
  // frictionless incompressible static contact problem
  //------------------------------------------------------------------------------

  static Elasticity::DisplacementType
  frictionless_incom_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
    if (T->center_mass().x() < 0 &&
        T->center_mass().z() < 0) {
        return Eigen::Vector3d(pow(x.x(), 3)*(sin(x.z()) + cos(x.y())) + pow(x.z(), 4)*cos((1.0/2.0)*PI*x.x())/lambda,
         -3*pow(x.x(), 2)*sin(x.y()) + 4*pow(x.z(), 3)*cos((1.0/2.0)*PI*x.x())/lambda,
         3*pow(x.x(), 2)*cos(x.z()) - 8*pow(x.z(), 3)*sin((1.0/2.0)*PI*x.x())/(PI*lambda));
    } 
    else if (T->center_mass().z() >= 0){
      return Eigen::Vector3d(
                            pow(x.x(), 3)*(sin(x.z()) + cos(x.y())) - pow(x.z(), 2)*sin(0.5*PI*x.x())*cos(0.5*PI*x.y())/lambda,
                             -3*pow(x.x(), 2)*sin(x.y()) + pow(x.z(), 2)/lambda,
                             3*pow(x.x(), 2)*cos(x.z()) + pow(x.x(), 2)*pow(x.z(), 2)/lambda);
    }
    else {
        return Eigen::Vector3d(pow(x.x(), 3)*(sin(x.z()) + cos(x.y())) + 2*pow(x.z(), 4)*cos((1.0/2.0)*PI*x.x())/lambda,
         -3*pow(x.x(), 2)*sin(x.y()) + 8*pow(x.z(), 3)*cos((1.0/2.0)*PI*x.x())/lambda,
         3*pow(x.x(), 2)*cos(x.z()) - 16*pow(x.z(), 3)*sin((1.0/2.0)*PI*x.x())/(PI*lambda));
    }
  }; 

  static Elasticity::FrictionType
  frictionless_incom_friction = [](const Face * F, const double & lambda) -> double{
    return 0.0;
  };

  static Elasticity::GradientDisplacementType
  frictionless_incom_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
    if (T->center_mass().x() < 0 &&
        T->center_mass().z() < 0) {
        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << -1.0/2.0*PI*pow(x.z(), 4)*sin((1.0/2.0)*PI*x.x())/lambda + 3*pow(x.x(), 2)*(sin(x.z()) + cos(x.y())),
                         -pow(x.x(), 3)*sin(x.y()), pow(x.x(), 3)*cos(x.z()) + 4*pow(x.z(), 3)*cos((1.0/2.0)*PI*x.x())/lambda;
                        M.row(1) << -2*PI*pow(x.z(), 3)*sin((1.0/2.0)*PI*x.x())/lambda - 6*x.x()*sin(x.y()),
                         -3*pow(x.x(), 2)*cos(x.y()), 12*pow(x.z(), 2)*cos((1.0/2.0)*PI*x.x())/lambda;
                        M.row(2) << 6*x.x()*cos(x.z()) - 4*pow(x.z(), 3)*cos((1.0/2.0)*PI*x.x())/lambda,
                         0, -3*pow(x.x(), 2)*sin(x.z()) - 24*pow(x.z(), 2)*sin((1.0/2.0)*PI*x.x())/(PI*lambda);
                        return M;
    } 
    else if (T->center_mass().z() >= 0) {
      Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << -0.5*PI*pow(x.z(), 2)*cos(0.5*PI*x.x())*cos(0.5*PI*x.y())/lambda + 3*pow(x.x(), 2)*(sin(x.z()) + cos(x.y())),
                         0.5*PI*pow(x.z(), 2)*sin(0.5*PI*x.x())*sin(0.5*PI*x.y())/lambda - pow(x.x(), 3)*sin(x.y()),
                          pow(x.x(), 3)*cos(x.z()) - 2*x.z()*sin(0.5*PI*x.x())*cos(0.5*PI*x.y())/lambda;
                        M.row(1) << -6*x.x()*sin(x.y()), -3*pow(x.x(), 2)*cos(x.y()), 2*x.z()/lambda;
                        M.row(2) << 6*x.x()*cos(x.z()) + 2*x.x()*pow(x.z(), 2)/lambda, 0, -3*pow(x.x(), 2)*sin(x.z()) + 2*pow(x.x(), 2)*x.z()/lambda;
                        return M;
    }
    else {
        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << -PI*pow(x.z(), 4)*sin((1.0/2.0)*PI*x.x())/lambda + 3*pow(x.x(), 2)*(sin(x.z()) + cos(x.y())),
                         -pow(x.x(), 3)*sin(x.y()), pow(x.x(), 3)*cos(x.z()) + 8*pow(x.z(), 3)*cos((1.0/2.0)*PI*x.x())/lambda;
                        M.row(1) << -4*PI*pow(x.z(), 3)*sin((1.0/2.0)*PI*x.x())/lambda - 6*x.x()*sin(x.y()),
                         -3*pow(x.x(), 2)*cos(x.y()), 24*pow(x.z(), 2)*cos((1.0/2.0)*PI*x.x())/lambda;
                        M.row(2) << 6*x.x()*cos(x.z()) - 8*pow(x.z(), 3)*cos((1.0/2.0)*PI*x.x())/lambda,
                         0, -3*pow(x.x(), 2)*sin(x.z()) - 48*pow(x.z(), 2)*sin((1.0/2.0)*PI*x.x())/(PI*lambda);
                        return M;
    }
  };       

  static Elasticity::GradientDisplacementType1
  frictionless_incom_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {
    Eigen::Matrix3d grad = frictionless_incom_grad_u(x,T,lambda);
    Eigen::Matrix3d eps = (0.5)*(grad+grad.transpose());
    Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
    trace_grad(0,0)=grad(0,0)+grad(1,1)+grad(2,2);  
    trace_grad(1,1)=grad(0,0)+grad(1,1)+grad(2,2);  
    trace_grad(2,2)=grad(0,0)+grad(1,1)+grad(2,2);  
    return 2*mu*eps + lambda*trace_grad;
  };

                                            

  static Elasticity::ForcingTermType
  frictionless_incom_div_eps_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {

    if (T->center_mass().x() < 0 && T->center_mass().z() < 0) {
        return Eigen::Vector3d(-0.25*pow(PI, 2)*pow(x.z(), 4)*cos((1.0/2.0)*PI*x.x())/lambda - 0.5*pow(x.x(), 3)*sin(x.z()) - 0.5*pow(x.x(), 3)*cos(x.y()) + 6.0*x.x()*(sin(x.z()) + cos(x.y())) - 3.0*x.x()*sin(x.z()) - 3.0*x.x()*cos(x.y()),
         -0.5*pow(PI, 2)*pow(x.z(), 3)*cos((1.0/2.0)*PI*x.x())/lambda + 1.5*pow(x.x(), 2)*sin(x.y()) - 3.0*sin(x.y()) + 12.0*x.z()*cos((1.0/2.0)*PI*x.x())/lambda,
         -1.5*pow(x.x(), 2)*cos(x.z()) + 3.0*cos(x.z()) - 48.0*x.z()*sin((1.0/2.0)*PI*x.x())/(PI*lambda));
    } 
    else if (T->center_mass().z() >= 0){
      return Eigen::Vector3d(0.375*pow(PI, 2)*pow(x.z(), 2)*sin(0.5*PI*x.x())*cos(0.5*PI*x.y())/lambda - 0.5*pow(x.x(), 3)*sin(x.z()) - 0.5*pow(x.x(), 3)*cos(x.y()) + 6.0*x.x()*(sin(x.z()) + cos(x.y())) - 3.0*x.x()*sin(x.z()) - 3.0*x.x()*cos(x.y()) + 2.0*x.x()*x.z()/lambda - 1.0*sin(0.5*PI*x.x())*cos(0.5*PI*x.y())/lambda,
       0.125*pow(PI, 2)*pow(x.z(), 2)*sin(0.5*PI*x.y())*cos(0.5*PI*x.x())/lambda + 1.5*pow(x.x(), 2)*sin(x.y()) - 3.0*sin(x.y()) + 1.0/lambda,
        -0.5*PI*x.z()*cos(0.5*PI*x.x())*cos(0.5*PI*x.y())/lambda - 1.5*pow(x.x(), 2)*cos(x.z()) + 3.0*cos(x.z()) + 2.0*pow(x.x(), 2)/lambda + 1.0*pow(x.z(), 2)/lambda); 
    }
    else {
        return Eigen::Vector3d(-0.5*pow(PI, 2)*pow(x.z(), 4)*cos((1.0/2.0)*PI*x.x())/lambda - 0.5*pow(x.x(), 3)*sin(x.z()) - 0.5*pow(x.x(), 3)*cos(x.y()) + 6.0*x.x()*(sin(x.z()) + cos(x.y())) - 3.0*x.x()*sin(x.z()) - 3.0*x.x()*cos(x.y()),
         -1.0*pow(PI, 2)*pow(x.z(), 3)*cos((1.0/2.0)*PI*x.x())/lambda + 1.5*pow(x.x(), 2)*sin(x.y()) - 3.0*sin(x.y()) + 24.0*x.z()*cos((1.0/2.0)*PI*x.x())/lambda,
         -1.5*pow(x.x(), 2)*cos(x.z()) + 3.0*cos(x.z()) - 96.0*x.z()*sin((1.0/2.0)*PI*x.x())/(PI*lambda)
        );
    }
};

  static Elasticity::ForcingTermType
  frictionless_incom_grad_div_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
    if (T->center_mass().x() < 0 && T->center_mass().z() < 0) {
        return Eigen::Vector3d(-1.0/4.0*pow(PI, 2)*pow(x.z(), 4)*cos((1.0/2.0)*PI*x.x())/lambda + 6*x.x()*(sin(x.z()) + cos(x.y())) - 6*x.x()*sin(x.z()) - 6*x.x()*cos(x.y()) - 12*pow(x.z(), 2)*cos((1.0/2.0)*PI*x.x())/lambda,
         0, -2*PI*pow(x.z(), 3)*sin((1.0/2.0)*PI*x.x())/lambda - 48*x.z()*sin((1.0/2.0)*PI*x.x())/(PI*lambda));
    } 
    else if (T->center_mass().z() >= 0) {
      return Eigen::Vector3d(0.25*pow(PI, 2)*pow(x.z(), 2)*sin(0.5*PI*x.x())*cos(0.5*PI*x.y())/lambda + 6*x.x()*(sin(x.z()) + cos(x.y())) - 6*x.x()*sin(x.z()) - 6*x.x()*cos(x.y()) + 4*x.x()*x.z()/lambda,
       0.25*pow(PI, 2)*pow(x.z(), 2)*sin(0.5*PI*x.y())*cos(0.5*PI*x.x())/lambda,
       -1.0*PI*x.z()*cos(0.5*PI*x.x())*cos(0.5*PI*x.y())/lambda + 2*pow(x.x(), 2)/lambda); 
    }
    else {
        return Eigen::Vector3d(-1.0/2.0*pow(PI, 2)*pow(x.z(), 4)*cos((1.0/2.0)*PI*x.x())/lambda + 6*x.x()*(sin(x.z()) + cos(x.y())) - 6*x.x()*sin(x.z()) - 6*x.x()*cos(x.y()) - 24*pow(x.z(), 2)*cos((1.0/2.0)*PI*x.x())/lambda,
         0, -4*PI*pow(x.z(), 3)*sin((1.0/2.0)*PI*x.x())/lambda - 96*x.z()*sin((1.0/2.0)*PI*x.x())/(PI*lambda));
    }
};

  static Elasticity::ForcingTermType
  frictionless_incom_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                    return -2 * mu * frictionless_incom_div_eps_u(mu,lambda,x,T) - lambda * frictionless_incom_grad_div_u(mu,lambda,x,T); 
                    };

//   //------------------------------------------------------------------------------
//   // frictionless incompressible static contact problem
//   //------------------------------------------------------------------------------

//   static Elasticity::DisplacementType
//   frictionless_incom_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
//     if (T->center_mass().x() < 0 &&
//         T->center_mass().z() < 0) {
//         return Eigen::Vector3d(cos((PI/2)*x.x())*std::pow(x.z(),4)/lambda,
//           4*cos((PI/2)*x.x())*std::pow(x.z(),3)/lambda+3*std::pow(x.x(),2)*sin(x.z()),
//           -(8/PI)*sin((PI/2)*x.x())*std::pow(x.z(),3)/lambda);
//     } 
//     else if (T->center_mass().z() >= 0){
//       return Eigen::Vector3d(
//                             -sin((PI/2)*x.x())*cos((PI/2)*x.y())*pow(x.z(),2)/lambda,
//                             pow(x.z(),2)/lambda+3*std::pow(x.x(),2)*sin(x.z()),
//                             pow(x.x(),2)*pow(x.z(),2)/lambda);
//     }
//     else {
//         return Eigen::Vector3d(2*cos((PI/2)*x.x())*std::pow(x.z(),4)/lambda,
//           8*cos((PI/2)*x.x())*std::pow(x.z(),3)/lambda+3*std::pow(x.x(),2)*sin(x.z()),
//           -(16/PI)*sin((PI/2)*x.x())*std::pow(x.z(),3)/lambda);
//     }
//   }; 

//   static Elasticity::FrictionType
//   frictionless_incom_friction = [](const Face * F, const double & lambda) -> double{
//     return 0.0;
//   };

//   static Elasticity::GradientDisplacementType
//   frictionless_incom_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
//     if (T->center_mass().x() < 0 &&
//         T->center_mass().z() < 0) {
//         Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
//                         M.row(0) << -1.0/2.0*PI*pow(x.z(), 4)*sin((1.0/2.0)*PI*x.x())/lambda,
//                          0, 4*pow(x.z(), 3)*cos((1.0/2.0)*PI*x.x())/lambda;
//                         M.row(1) << -2*PI*pow(x.z(), 3)*sin((1.0/2.0)*PI*x.x())/lambda + 6*x.x()*sin(x.z()),
//                          0, 3*pow(x.x(), 2)*cos(x.z()) + 12*pow(x.z(), 2)*cos((1.0/2.0)*PI*x.x())/lambda;
//                         M.row(2) << -4*pow(x.z(), 3)*cos((1.0/2.0)*PI*x.x())/lambda,
//                          0, -24*pow(x.z(), 2)*sin((1.0/2.0)*PI*x.x())/(PI*lambda);
//                         return M;
//     } 
//     else if (T->center_mass().z() >= 0) {
//       Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
//                         M.row(0) << -0.5*PI*pow(x.z(), 2)*cos(0.5*PI*x.x())*cos(0.5*PI*x.y())/lambda, 
//                         0.5*PI*pow(x.z(), 2)*sin(0.5*PI*x.x())*sin(0.5*PI*x.y())/lambda, -2*x.z()*sin(0.5*PI*x.x())*cos(0.5*PI*x.y())/lambda;
//                         M.row(1) << 6*x.x()*sin(x.z()), 0, 3*pow(x.x(), 2)*cos(x.z()) + 2*x.z()/lambda;
//                         M.row(2) << 2*x.x()*pow(x.z(), 2)/lambda, 0, 2*pow(x.x(), 2)*x.z()/lambda;
//                         return M;
//     }
//     else {
//         Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
//                         M.row(0) << -PI*pow(x.z(), 4)*sin((1.0/2.0)*PI*x.x())/lambda,
//                          0, 8*pow(x.z(), 3)*cos((1.0/2.0)*PI*x.x())/lambda;
//                         M.row(1) << -4*PI*pow(x.z(), 3)*sin((1.0/2.0)*PI*x.x())/lambda + 6*x.x()*sin(x.z()),
//                          0, 3*pow(x.x(), 2)*cos(x.z()) + 24*pow(x.z(), 2)*cos((1.0/2.0)*PI*x.x())/lambda;
//                         M.row(2) << -8*pow(x.z(), 3)*cos((1.0/2.0)*PI*x.x())/lambda,
//                          0, -48*pow(x.z(), 2)*sin((1.0/2.0)*PI*x.x())/(PI*lambda);
//                         return M;
//     }
//   };       

//   static Elasticity::GradientDisplacementType1
//   frictionless_incom_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {
//     Eigen::Matrix3d grad = frictionless_incom_grad_u(x,T,lambda);
//     Eigen::Matrix3d eps = (0.5)*(grad+grad.transpose());
//     Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
//     trace_grad(0,0)=grad(0,0)+grad(1,1)+grad(2,2);  
//     trace_grad(1,1)=grad(0,0)+grad(1,1)+grad(2,2);  
//     trace_grad(2,2)=grad(0,0)+grad(1,1)+grad(2,2);  
//     return 2*mu*eps + lambda*trace_grad;
//   };

                                            

//   static Elasticity::ForcingTermType
//   frictionless_incom_div_eps_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {

//     if (T->center_mass().x() < 0 && T->center_mass().z() < 0) {
//         return Eigen::Vector3d(-0.25*pow(PI, 2)*pow(x.z(), 4)*cos((1.0/2.0)*PI*x.x())/lambda,
//          -0.5*pow(PI, 2)*pow(x.z(), 3)*cos((1.0/2.0)*PI*x.x())/lambda - 1.5*pow(x.x(), 2)*sin(x.z()) + 3.0*sin(x.z()) + 12.0*x.z()*cos((1.0/2.0)*PI*x.x())/lambda,
//          -48.0*x.z()*sin((1.0/2.0)*PI*x.x())/(PI*lambda));
//     } 
//     else if (T->center_mass().z() >= 0){
//       return Eigen::Vector3d(0.375*pow(PI, 2)*pow(x.z(), 2)*sin(0.5*PI*x.x())*cos(0.5*PI*x.y())/lambda + 2.0*x.x()*x.z()/lambda - 1.0*sin(0.5*PI*x.x())*cos(0.5*PI*x.y())/lambda,
//        0.125*pow(PI, 2)*pow(x.z(), 2)*sin(0.5*PI*x.y())*cos(0.5*PI*x.x())/lambda - 1.5*pow(x.x(), 2)*sin(x.z()) + 3.0*sin(x.z()) + 1.0/lambda,
//        -0.5*PI*x.z()*cos(0.5*PI*x.x())*cos(0.5*PI*x.y())/lambda + 2.0*pow(x.x(), 2)/lambda + 1.0*pow(x.z(), 2)/lambda
//                                 ); 
//     }
//     else {
//         return Eigen::Vector3d(-0.5*pow(PI, 2)*pow(x.z(), 4)*cos((1.0/2.0)*PI*x.x())/lambda,
//          -1.0*pow(PI, 2)*pow(x.z(), 3)*cos((1.0/2.0)*PI*x.x())/lambda - 1.5*pow(x.x(), 2)*sin(x.z()) + 3.0*sin(x.z()) + 24.0*x.z()*cos((1.0/2.0)*PI*x.x())/lambda,
//          -96.0*x.z()*sin((1.0/2.0)*PI*x.x())/(PI*lambda)
//         );
//     }
// };

//   static Elasticity::ForcingTermType
//   frictionless_incom_grad_div_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
//     if (T->center_mass().x() < 0 && T->center_mass().z() < 0) {
//         return Eigen::Vector3d(-1.0/4.0*pow(PI, 2)*pow(x.z(), 4)*cos((1.0/2.0)*PI*x.x())/lambda - 12*pow(x.z(), 2)*cos((1.0/2.0)*PI*x.x())/lambda,
//          0, -2*PI*pow(x.z(), 3)*sin((1.0/2.0)*PI*x.x())/lambda - 48*x.z()*sin((1.0/2.0)*PI*x.x())/(PI*lambda));
//     } 
//     else if (T->center_mass().z() >= 0) {
//       return Eigen::Vector3d(0.25*pow(PI, 2)*pow(x.z(), 2)*sin(0.5*PI*x.x())*cos(0.5*PI*x.y())/lambda + 4*x.x()*x.z()/lambda,
//        0.25*pow(PI, 2)*pow(x.z(), 2)*sin(0.5*PI*x.y())*cos(0.5*PI*x.x())/lambda,
//        -1.0*PI*x.z()*cos(0.5*PI*x.x())*cos(0.5*PI*x.y())/lambda + 2*pow(x.x(), 2)/lambda
//                                 ); 
//     }
//     else {
//         return Eigen::Vector3d(-1.0/2.0*pow(PI, 2)*pow(x.z(), 4)*cos((1.0/2.0)*PI*x.x())/lambda - 24*pow(x.z(), 2)*cos((1.0/2.0)*PI*x.x())/lambda,
//          0,-4*PI*pow(x.z(), 3)*sin((1.0/2.0)*PI*x.x())/lambda - 96*x.z()*sin((1.0/2.0)*PI*x.x())/(PI*lambda)
//         );
//     }
// };

//   static Elasticity::ForcingTermType
//   frictionless_incom_f = 
//           [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
//                     return -2 * mu * frictionless_incom_div_eps_u(mu,lambda,x,T) - lambda * frictionless_incom_grad_div_u(mu,lambda,x,T); 
//                     };

  //------------------------------------------------------------------------------
  // friction static contact problem
  //------------------------------------------------------------------------------

  static Elasticity::FrictionType
  friction_friction = [](const Face * F, const double & lambda) -> double{
    return 1.0;
  };

  static Elasticity::DisplacementType
  friction_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
    if (T->center_mass().z() >= 0){
      return Eigen::Vector3d(
                            -sin(x.x())*cos(x.y())*pow(x.z(),2)-x.y(),
                            pow(x.z(),2),
                            pow(x.x(),2)*pow(x.z(),2));
    }
    else if (T->center_mass().x() >= 0 &&
        T->center_mass().z() < 0) {
        return Eigen::Vector3d(
                            (-0.25)*sin(x.x())*cos(x.y())*pow(x.z(),2)-x.y(),
                            (0.25)*pow(x.z(),2),
                            (0.25)*pow(x.x(),2)*pow(x.z(),2));
    }
    else {
        return Eigen::Vector3d(
                            (-0.25)*sin(x.x())*cos(x.y())*pow(x.z(),2)-x.y(),
                            (0.50)*pow(x.z(),2),
                            (0.25)*pow(x.x(),2)*pow(x.z(),2));
    }
  }; 

  static Elasticity::GradientDisplacementType
  friction_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
    if (T->center_mass().z() >= 0) {
      Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << -pow(x.z(), 2)*cos(x.x())*cos(x.y()), 
                        pow(x.z(), 2)*sin(x.x())*sin(x.y()) - 1, -(2*x.z())*sin(x.x())*cos(x.y());
                        M.row(1) << 0, 0, 2*(x.z());
                        M.row(2) << (2*x.x())*pow(x.z(), 2), 0, pow(x.x(), 2)*(2*x.z());
                        return M;
    }
    else if (T->center_mass().z() < 0 && T->center_mass().x() >= 0) {
        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << -0.25*pow(x.z(), 2)*cos(x.x())*cos(x.y()), 
                        0.25*pow(x.z(), 2)*sin(x.x())*sin(x.y()) - 1, -0.25*(2*x.z())*sin(x.x())*cos(x.y());
                        M.row(1) << 0, 0, (0.5)*(x.z());
                        M.row(2) << 0.25*(2*x.x())*pow(x.z(), 2), 0, 0.25*pow(x.x(), 2)*(2*x.z());
                        return M;
    }
    else {
        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << -0.25*pow(x.z(), 2)*cos(x.x())*cos(x.y()), 
                        0.25*pow(x.z(), 2)*sin(x.x())*sin(x.y()) - 1, -0.25*(2*x.z())*sin(x.x())*cos(x.y());
                        M.row(1) << 0, 0, (x.z());
                        M.row(2) << 0.25*(2*x.x())*pow(x.z(), 2), 0, 0.25*pow(x.x(), 2)*(2*x.z());
                        return M;
    }
  };       

  static Elasticity::GradientDisplacementType1
  friction_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {
    Eigen::Matrix3d grad = friction_grad_u(x,T,lambda);
    Eigen::Matrix3d eps = (0.5)*(grad+grad.transpose());
    Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
    trace_grad(0,0)=grad(0,0)+grad(1,1)+grad(2,2);  
    trace_grad(1,1)=grad(0,0)+grad(1,1)+grad(2,2);  
    trace_grad(2,2)=grad(0,0)+grad(1,1)+grad(2,2);  
    return 2*mu*eps + lambda*trace_grad;
  };                                            

  static Elasticity::ForcingTermType
  friction_div_eps_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
    if (T->center_mass().z() >= 0){
          return Eigen::Vector3d(
                                0.5*(2*x.x())*(2*x.z()) + 1.5*pow(x.z(), 2)*sin(x.x())*cos(x.y()) - 1.0*sin(x.x())*cos(x.y()), 
                                0.5*pow(x.z(), 2)*sin(x.y())*cos(x.x()) + 1.0, 
                                2.0*pow(x.x(), 2) + 1.0*pow(x.z(), 2) - 0.5*(2*x.z())*cos(x.x())*cos(x.y())
                                ); 
    }
    else if (T->center_mass().x() >= 0 &&
        T->center_mass().z() < 0) {
          return Eigen::Vector3d(
                                0.125*(2*x.x())*(2*x.z()) + 0.375*pow(x.z(), 2)*sin(x.x())*cos(x.y()) - 0.25*sin(x.x())*cos(x.y()), 
                                0.125*pow(x.z(), 2)*sin(x.y())*cos(x.x()) + 0.25, 
                                0.5*pow(x.x(), 2) + 0.25*pow(x.z(), 2) - 0.125*(2*x.z())*cos(x.x())*cos(x.y())
                                );
    }
    else {
        return Eigen::Vector3d(
                                0.125*(2*x.x())*(2*x.z()) + 0.375*pow(x.z(), 2)*sin(x.x())*cos(x.y()) - 0.25*sin(x.x())*cos(x.y()), 
                                0.125*pow(x.z(), 2)*sin(x.y())*cos(x.x()) + 0.5, 
                                0.5*pow(x.x(), 2) + 0.25*pow(x.z(), 2) - 0.125*(2*x.z())*cos(x.x())*cos(x.y())
                                );
    }
};

  static Elasticity::ForcingTermType
  friction_grad_div_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
    if (T->center_mass().z() >= 0) {
      return Eigen::Vector3d(
                                (2*x.x())*(2*x.z()) + pow(x.z(), 2)*sin(x.x())*cos(x.y()), 
                                pow(x.z(), 2)*sin(x.y())*cos(x.x()), 
                                2*pow(x.x(), 2) - (2*x.z())*cos(x.x())*cos(x.y())
                                ); 
    }
    else if (T->center_mass().x() >= 0 &&
        T->center_mass().z() < 0) {
        return Eigen::Vector3d(
                                0.25*(2*x.x())*(2*x.z()) + 0.25*pow(x.z(), 2)*sin(x.x())*cos(x.y()), 
                                0.25*pow(x.z(), 2)*sin(x.y())*cos(x.x()), 
                                0.5*pow(x.x(), 2) - 0.25*(2*x.z())*cos(x.x())*cos(x.y())
                                );
    }
    else {
        return Eigen::Vector3d(
                                0.25*(2*x.x())*(2*x.z()) + 0.25*pow(x.z(), 2)*sin(x.x())*cos(x.y()), 
                                0.25*pow(x.z(), 2)*sin(x.y())*cos(x.x()), 
                                0.5*pow(x.x(), 2) - 0.25*(2*x.z())*cos(x.x())*cos(x.y())
                                );
    }
};

  static Elasticity::ForcingTermType
  friction_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                    return -2 * mu * friction_div_eps_u(mu,lambda,x,T) - lambda * friction_grad_div_u(mu,lambda,x,T); 
                    };
                    
  //------------------------------------------------------------------------------
  // friction inompressible static contact problem
  //------------------------------------------------------------------------------

  static Elasticity::FrictionType
  friction_incom_friction = [](const Face * F, const double & lambda) -> double{
    return 1.0/lambda;
  };

  static Elasticity::DisplacementType
  friction_incom_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
    if (T->center_mass().z() >= 0){
      return Eigen::Vector3d(
                            pow(x.x(), 3)*(sin(x.z()) + cos(x.y())) + (-x.y() - pow(x.z(), 2)*sin(x.x())*cos(x.y()))/lambda,
                             -3*pow(x.x(), 2)*sin(x.y()) + pow(x.z(), 2)/lambda,
                              3*pow(x.x(), 2)*cos(x.z()) + pow(x.x(), 2)*pow(x.z(), 2)/lambda);
    }
    else if (T->center_mass().x() >= 0 &&
        T->center_mass().z() < 0) {
        return Eigen::Vector3d(pow(x.x(), 3)*(sin(x.z()) + cos(x.y())) + (-x.y() - 0.25*pow(x.z(), 2)*sin(x.x())*cos(x.y()))/lambda,
         -3*pow(x.x(), 2)*sin(x.y()) + 0.25*pow(x.z(), 2)/lambda,
         3*pow(x.x(), 2)*cos(x.z()) + 0.25*pow(x.x(), 2)*pow(x.z(), 2)/lambda);
    }
    else {
        return Eigen::Vector3d(pow(x.x(), 3)*(sin(x.z()) + cos(x.y())) + (-x.y() - 0.25*pow(x.z(), 2)*sin(x.x())*cos(x.y()))/lambda,
         -3*pow(x.x(), 2)*sin(x.y()) + 0.5*pow(x.z(), 2)/lambda,
         3*pow(x.x(), 2)*cos(x.z()) + 0.25*pow(x.x(), 2)*pow(x.z(), 2)/lambda);
    }
  }; 

  static Elasticity::GradientDisplacementType
  friction_incom_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
    if (T->center_mass().z() >= 0) {
      Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << 3*pow(x.x(), 2)*(sin(x.z()) + cos(x.y())) - pow(x.z(), 2)*cos(x.x())*cos(x.y())/lambda,
                         -pow(x.x(), 3)*sin(x.y()) + (pow(x.z(), 2)*sin(x.x())*sin(x.y()) - 1)/lambda,
                          pow(x.x(), 3)*cos(x.z()) - 2*x.z()*sin(x.x())*cos(x.y())/lambda;
                        M.row(1) << -6*x.x()*sin(x.y()), -3*pow(x.x(), 2)*cos(x.y()), 2*x.z()/lambda;
                        M.row(2) << 6*x.x()*cos(x.z()) + 2*x.x()*pow(x.z(), 2)/lambda, 0, -3*pow(x.x(), 2)*sin(x.z()) + 2*pow(x.x(), 2)*x.z()/lambda;
                        return M;
    }
    else if (T->center_mass().z() < 0 && T->center_mass().x() >= 0) {
        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << 3*pow(x.x(), 2)*(sin(x.z()) + cos(x.y())) - 0.25*pow(x.z(), 2)*cos(x.x())*cos(x.y())/lambda,
                         -pow(x.x(), 3)*sin(x.y()) + (0.25*pow(x.z(), 2)*sin(x.x())*sin(x.y()) - 1)/lambda,
                         pow(x.x(), 3)*cos(x.z()) - 0.5*x.z()*sin(x.x())*cos(x.y())/lambda;
                        M.row(1) << -6*x.x()*sin(x.y()), -3*pow(x.x(), 2)*cos(x.y()), 0.5*x.z()/lambda;
                        M.row(2) << 6*x.x()*cos(x.z()) + 0.5*x.x()*pow(x.z(), 2)/lambda, 0, -3*pow(x.x(), 2)*sin(x.z()) + 0.5*pow(x.x(), 2)*x.z()/lambda;
                        return M;
    }
    else {
        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << 3*pow(x.x(), 2)*(sin(x.z()) + cos(x.y())) - 0.25*pow(x.z(), 2)*cos(x.x())*cos(x.y())/lambda,
                         -pow(x.x(), 3)*sin(x.y()) + (0.25*pow(x.z(), 2)*sin(x.x())*sin(x.y()) - 1)/lambda,
                          pow(x.x(), 3)*cos(x.z()) - 0.5*x.z()*sin(x.x())*cos(x.y())/lambda;
                        M.row(1) << -6*x.x()*sin(x.y()), -3*pow(x.x(), 2)*cos(x.y()), 1.0*x.z()/lambda;
                        M.row(2) << 6*x.x()*cos(x.z()) + 0.5*x.x()*pow(x.z(), 2)/lambda, 0, -3*pow(x.x(), 2)*sin(x.z()) + 0.5*pow(x.x(), 2)*x.z()/lambda;
                        return M;
    }
  };       

  static Elasticity::GradientDisplacementType1
  friction_incom_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {
    Eigen::Matrix3d grad = friction_incom_grad_u(x,T,lambda);
    Eigen::Matrix3d eps = (0.5)*(grad+grad.transpose());
    Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
    trace_grad(0,0)=grad(0,0)+grad(1,1)+grad(2,2);  
    trace_grad(1,1)=grad(0,0)+grad(1,1)+grad(2,2);  
    trace_grad(2,2)=grad(0,0)+grad(1,1)+grad(2,2);  
    return 2*mu*eps + lambda*trace_grad;
  };                                            

  static Elasticity::ForcingTermType
  friction_incom_div_eps_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
    if (T->center_mass().z() >= 0){
          return Eigen::Vector3d(-0.5*pow(x.x(), 3)*sin(x.z()) - 0.5*pow(x.x(), 3)*cos(x.y()) + 6.0*x.x()*(sin(x.z()) + cos(x.y())) - 3.0*x.x()*sin(x.z()) - 3.0*x.x()*cos(x.y()) + 2.0*x.x()*x.z()/lambda + 1.5*pow(x.z(), 2)*sin(x.x())*cos(x.y())/lambda - 1.0*sin(x.x())*cos(x.y())/lambda,
           1.5*pow(x.x(), 2)*sin(x.y()) - 3.0*sin(x.y()) + 0.5*pow(x.z(), 2)*sin(x.y())*cos(x.x())/lambda + 1.0/lambda,
            -1.5*pow(x.x(), 2)*cos(x.z()) + 3.0*cos(x.z()) + 2.0*pow(x.x(), 2)/lambda + 1.0*pow(x.z(), 2)/lambda - 1.0*x.z()*cos(x.x())*cos(x.y())/lambda); 
    }
    else if (T->center_mass().x() >= 0 &&
        T->center_mass().z() < 0) {
          return Eigen::Vector3d(-0.5*pow(x.x(), 3)*sin(x.z()) - 0.5*pow(x.x(), 3)*cos(x.y()) + 6.0*x.x()*(sin(x.z()) + cos(x.y())) - 3.0*x.x()*sin(x.z()) - 3.0*x.x()*cos(x.y()) + 0.5*x.x()*x.z()/lambda + 0.375*pow(x.z(), 2)*sin(x.x())*cos(x.y())/lambda - 0.25*sin(x.x())*cos(x.y())/lambda,
           1.5*pow(x.x(), 2)*sin(x.y()) - 3.0*sin(x.y()) + 0.125*pow(x.z(), 2)*sin(x.y())*cos(x.x())/lambda + 0.25/lambda,
            -1.5*pow(x.x(), 2)*cos(x.z()) + 3.0*cos(x.z()) + 0.5*pow(x.x(), 2)/lambda + 0.25*pow(x.z(), 2)/lambda - 0.25*x.z()*cos(x.x())*cos(x.y())/lambda);
    }
    else {
        return Eigen::Vector3d(-0.5*pow(x.x(), 3)*sin(x.z()) - 0.5*pow(x.x(), 3)*cos(x.y()) + 6.0*x.x()*(sin(x.z()) + cos(x.y())) - 3.0*x.x()*sin(x.z()) - 3.0*x.x()*cos(x.y()) + 0.5*x.x()*x.z()/lambda + 0.375*pow(x.z(), 2)*sin(x.x())*cos(x.y())/lambda - 0.25*sin(x.x())*cos(x.y())/lambda,
         1.5*pow(x.x(), 2)*sin(x.y()) - 3.0*sin(x.y()) + 0.125*pow(x.z(), 2)*sin(x.y())*cos(x.x())/lambda + 0.5/lambda,
          -1.5*pow(x.x(), 2)*cos(x.z()) + 3.0*cos(x.z()) + 0.5*pow(x.x(), 2)/lambda + 0.25*pow(x.z(), 2)/lambda - 0.25*x.z()*cos(x.x())*cos(x.y())/lambda);
    }
};

  static Elasticity::ForcingTermType
  friction_incom_grad_div_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
    if (T->center_mass().z() >= 0) {
      return Eigen::Vector3d(6*x.x()*(sin(x.z()) + cos(x.y())) - 6*x.x()*sin(x.z()) - 6*x.x()*cos(x.y()) + 4*x.x()*x.z()/lambda + pow(x.z(), 2)*sin(x.x())*cos(x.y())/lambda,
       pow(x.z(), 2)*sin(x.y())*cos(x.x())/lambda,
       2*pow(x.x(), 2)/lambda - 2*x.z()*cos(x.x())*cos(x.y())/lambda); 
    }
    else if (T->center_mass().x() >= 0 &&
        T->center_mass().z() < 0) {
        return Eigen::Vector3d(6*x.x()*(sin(x.z()) + cos(x.y())) - 6*x.x()*sin(x.z()) - 6*x.x()*cos(x.y()) + 1.0*x.x()*x.z()/lambda + 0.25*pow(x.z(), 2)*sin(x.x())*cos(x.y())/lambda,
         0.25*pow(x.z(), 2)*sin(x.y())*cos(x.x())/lambda,
          0.5*pow(x.x(), 2)/lambda - 0.5*x.z()*cos(x.x())*cos(x.y())/lambda);
    }
    else {
        return Eigen::Vector3d(6*x.x()*(sin(x.z()) + cos(x.y())) - 6*x.x()*sin(x.z()) - 6*x.x()*cos(x.y()) + 1.0*x.x()*x.z()/lambda + 0.25*pow(x.z(), 2)*sin(x.x())*cos(x.y())/lambda, 
        0.25*pow(x.z(), 2)*sin(x.y())*cos(x.x())/lambda,
         0.5*pow(x.x(), 2)/lambda - 0.5*x.z()*cos(x.x())*cos(x.y())/lambda);
    }
};

  static Elasticity::ForcingTermType
  friction_incom_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                    return -2 * mu * friction_incom_div_eps_u(mu,lambda,x,T) - lambda * friction_incom_grad_div_u(mu,lambda,x,T); 
                    };                    

//   //------------------------------------------------------------------------------
//   // friction inompressible static contact problem
//   //------------------------------------------------------------------------------

//   static Elasticity::FrictionType
//   friction_incom_friction = [](const Face * F, const double & lambda) -> double{
//     return 1.0/lambda;
//   };

//   static Elasticity::DisplacementType
//   friction_incom_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
//     if (T->center_mass().z() >= 0){
//       return Eigen::Vector3d(
//                             (-sin(x.x())*cos(x.y())*pow(x.z(),2)-x.y())/lambda,
//                             pow(x.z(),2)/lambda+3*std::pow(x.x(),2)*sin(x.z()),
//                             pow(x.x(),2)*pow(x.z(),2)/lambda);
//     }
//     else if (T->center_mass().x() >= 0 &&
//         T->center_mass().z() < 0) {
//         return Eigen::Vector3d(
//                             ((-0.25)*sin(x.x())*cos(x.y())*pow(x.z(),2)-x.y())/lambda,
//                             (0.25)*pow(x.z(),2)/lambda+3*std::pow(x.x(),2)*sin(x.z()),
//                             (0.25)*pow(x.x(),2)*pow(x.z(),2)/lambda);
//     }
//     else {
//         return Eigen::Vector3d(
//                             ((-0.25)*sin(x.x())*cos(x.y())*pow(x.z(),2)-x.y())/lambda,
//                             (0.50)*pow(x.z(),2)/lambda+3*std::pow(x.x(),2)*sin(x.z()),
//                             (0.25)*pow(x.x(),2)*pow(x.z(),2)/lambda);
//     }
//   }; 

//   static Elasticity::GradientDisplacementType
//   friction_incom_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
//     if (T->center_mass().z() >= 0) {
//       Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
//                         M.row(0) << -pow(x.z(), 2)*cos(x.x())*cos(x.y())/lambda,
//                          (pow(x.z(), 2)*sin(x.x())*sin(x.y()) - 1)/lambda, -2*x.z()*sin(x.x())*cos(x.y())/lambda;
//                         M.row(1) << 6*x.x()*sin(x.z()), 0, 3*pow(x.x(), 2)*cos(x.z()) + 2*x.z()/lambda;
//                         M.row(2) << 2*x.x()*pow(x.z(), 2)/lambda, 0, 2*pow(x.x(), 2)*x.z()/lambda;
//                         return M;
//     }
//     else if (T->center_mass().z() < 0 && T->center_mass().x() >= 0) {
//         Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
//                         M.row(0) << -0.25*pow(x.z(), 2)*cos(x.x())*cos(x.y())/lambda,
//                          (0.25*pow(x.z(), 2)*sin(x.x())*sin(x.y()) - 1)/lambda, -0.5*x.z()*sin(x.x())*cos(x.y())/lambda;
//                         M.row(1) << 6*x.x()*sin(x.z()), 0, 3*pow(x.x(), 2)*cos(x.z()) + 0.5*x.z()/lambda;
//                         M.row(2) << 0.5*x.x()*pow(x.z(), 2)/lambda, 0, 0.5*pow(x.x(), 2)*x.z()/lambda;
//                         return M;
//     }
//     else {
//         Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
//                         M.row(0) << -0.25*pow(x.z(), 2)*cos(x.x())*cos(x.y())/lambda,
//                          (0.25*pow(x.z(), 2)*sin(x.x())*sin(x.y()) - 1)/lambda, -0.5*x.z()*sin(x.x())*cos(x.y())/lambda;
//                         M.row(1) << 6*x.x()*sin(x.z()), 0, 3*pow(x.x(), 2)*cos(x.z()) + 1.0*x.z()/lambda;
//                         M.row(2) << 0.5*x.x()*pow(x.z(), 2)/lambda, 0, 0.5*pow(x.x(), 2)*x.z()/lambda;
//                         return M;
//     }
//   };       

//   static Elasticity::GradientDisplacementType1
//   friction_incom_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {
//     Eigen::Matrix3d grad = friction_incom_grad_u(x,T,lambda);
//     Eigen::Matrix3d eps = (0.5)*(grad+grad.transpose());
//     Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
//     trace_grad(0,0)=grad(0,0)+grad(1,1)+grad(2,2);  
//     trace_grad(1,1)=grad(0,0)+grad(1,1)+grad(2,2);  
//     trace_grad(2,2)=grad(0,0)+grad(1,1)+grad(2,2);  
//     return 2*mu*eps + lambda*trace_grad;
//   };                                            

//   static Elasticity::ForcingTermType
//   friction_incom_div_eps_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
//     if (T->center_mass().z() >= 0){
//           return Eigen::Vector3d(
//                                 2.0*x.x()*x.z()/lambda + 1.5*pow(x.z(), 2)*sin(x.x())*cos(x.y())/lambda - 1.0*sin(x.x())*cos(x.y())/lambda,
//                                  -1.5*pow(x.x(), 2)*sin(x.z()) + 3.0*sin(x.z()) + 0.5*pow(x.z(), 2)*sin(x.y())*cos(x.x())/lambda + 1.0/lambda,
//                                  2.0*pow(x.x(), 2)/lambda + 1.0*pow(x.z(), 2)/lambda - 1.0*x.z()*cos(x.x())*cos(x.y())/lambda
//                                 ); 
//     }
//     else if (T->center_mass().x() >= 0 &&
//         T->center_mass().z() < 0) {
//           return Eigen::Vector3d(
//                                 0.5*x.x()*x.z()/lambda + 0.375*pow(x.z(), 2)*sin(x.x())*cos(x.y())/lambda - 0.25*sin(x.x())*cos(x.y())/lambda,
//                                  -1.5*pow(x.x(), 2)*sin(x.z()) + 3.0*sin(x.z()) + 0.125*pow(x.z(), 2)*sin(x.y())*cos(x.x())/lambda + 0.25/lambda,
//                                   0.5*pow(x.x(), 2)/lambda + 0.25*pow(x.z(), 2)/lambda - 0.25*x.z()*cos(x.x())*cos(x.y())/lambda
//                                 );
//     }
//     else {
//         return Eigen::Vector3d(
//                                 0.5*x.x()*x.z()/lambda + 0.375*pow(x.z(), 2)*sin(x.x())*cos(x.y())/lambda - 0.25*sin(x.x())*cos(x.y())/lambda,
//                                  -1.5*pow(x.x(), 2)*sin(x.z()) + 3.0*sin(x.z()) + 0.125*pow(x.z(), 2)*sin(x.y())*cos(x.x())/lambda + 0.5/lambda,
//                                   0.5*pow(x.x(), 2)/lambda + 0.25*pow(x.z(), 2)/lambda - 0.25*x.z()*cos(x.x())*cos(x.y())/lambda
//                                 );
//     }
// };

//   static Elasticity::ForcingTermType
//   friction_incom_grad_div_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
//     if (T->center_mass().z() >= 0) {
//       return Eigen::Vector3d(4*x.x()*x.z()/lambda + pow(x.z(), 2)*sin(x.x())*cos(x.y())/lambda,
//        pow(x.z(), 2)*sin(x.y())*cos(x.x())/lambda,
//        2*pow(x.x(), 2)/lambda - 2*x.z()*cos(x.x())*cos(x.y())/lambda
//                                 ); 
//     }
//     else if (T->center_mass().x() >= 0 &&
//         T->center_mass().z() < 0) {
//         return Eigen::Vector3d(1.0*x.x()*x.z()/lambda + 0.25*pow(x.z(), 2)*sin(x.x())*cos(x.y())/lambda,
//          0.25*pow(x.z(), 2)*sin(x.y())*cos(x.x())/lambda,
//           0.5*pow(x.x(), 2)/lambda - 0.5*x.z()*cos(x.x())*cos(x.y())/lambda
//                                 );
//     }
//     else {
//         return Eigen::Vector3d(1.0*x.x()*x.z()/lambda + 0.25*pow(x.z(), 2)*sin(x.x())*cos(x.y())/lambda,
//          0.25*pow(x.z(), 2)*sin(x.y())*cos(x.x())/lambda,
//           0.5*pow(x.x(), 2)/lambda - 0.5*x.z()*cos(x.x())*cos(x.y())/lambda
//                                 );
//     }
// };

//   static Elasticity::ForcingTermType
//   friction_incom_f = 
//           [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
//                     return -2 * mu * friction_incom_div_eps_u(mu,lambda,x,T) - lambda * friction_incom_grad_div_u(mu,lambda,x,T); 
//                     };

  //------------------------------------------------------------------------------
  // friction unknown solution problem
  //------------------------------------------------------------------------------

  static Elasticity::FrictionType
  unknown_friction = [](const Face * F, const double & lambda) -> double{
    // return 5e6;
    return 4e6;
  };

  static Elasticity::DisplacementType
  unknown_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
    if (T->center_mass().z() > 0.5){
        // return Eigen::Vector3d(0.0, 0.005, -0.002);
        return Eigen::Vector3d(0.0015, -0.0015, -0.002);
    }
    else {
        return Eigen::Vector3d(0.0, 0.0, 0.0);
    }
  }; 

  static Elasticity::GradientDisplacementType
  unknown_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
      Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        return M;
  };       

  static Elasticity::GradientDisplacementType1
  unknown_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {
    Eigen::Matrix3d grad = unknown_grad_u(x,T,lambda);
    Eigen::Matrix3d eps = (0.5)*(grad+grad.transpose());
    Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
    trace_grad(0,0)=grad(0,0)+grad(1,1)+grad(2,2);  
    trace_grad(1,1)=grad(0,0)+grad(1,1)+grad(2,2);  
    trace_grad(2,2)=grad(0,0)+grad(1,1)+grad(2,2);  
    return 2*mu*eps + lambda*trace_grad;
  };                                            

  static Elasticity::ForcingTermType
  unknown_div_eps_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
    return Eigen::Vector3d(0.0, 0.0, 0.0); 
  };

  static Elasticity::ForcingTermType
  unknown_grad_div_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
    return Eigen::Vector3d(0.0, 0.0, 0.0); 
  };

  static Elasticity::ForcingTermType
  unknown_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                    return -2 * mu * unknown_div_eps_u(mu,lambda,x,T) - lambda * unknown_grad_div_u(mu,lambda,x,T); 
                    };                    

  //------------------------------------------------------------------------------
  // polynomial solution with 0 tangential component of \sigma
  //------------------------------------------------------------------------------

  static Elasticity::DisplacementType
  polynomial_new_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
                      return Eigen::Vector3d(-pow(x.x()-0.5,2),0,0);
                    };   
  
  static Elasticity::GradientDisplacementType
  polynomial_new_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
                        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << -2*(x.x()-0.5), 0, 0;
                        M.row(1) << 0, 0, 0;
                        M.row(2) << 0, 0, 0;
                        return M;
                         };         
                         
  static Elasticity::GradientDisplacementType1
  polynomial_new_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {
        
        Eigen::Matrix3d grad = Eigen::Matrix3d::Zero();

        grad.row(0) << -2*(x.x()-0.5), 0, 0;
        grad.row(1) << 0, 0, 0;
        grad.row(2) << 0, 0, 0;

        // All other entries remain zero since other components are 0

        Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
        trace_grad(0,0)=grad(0,0)+grad(1,1)+grad(2,2);  
        trace_grad(1,1)=grad(0,0)+grad(1,1)+grad(2,2);  
        trace_grad(2,2)=grad(0,0)+grad(1,1)+grad(2,2);  
     
        return 2*mu*grad + lambda*trace_grad;
  };                       


  static Elasticity::ForcingTermType
  polynomial_new_div_eps_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                      return Eigen::Vector3d(-2,0,0); 
                    };

  static Elasticity::ForcingTermType
  polynomial_new_grad_div_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                      return Eigen::Vector3d(-2,0,0);  
                    };

  static Elasticity::ForcingTermType
  polynomial_new_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                    return -2 * mu * polynomial_new_div_eps_u(mu,lambda,x,T) - lambda * polynomial_new_grad_div_u(mu,lambda,x,T); 
                    };                          


  //------------------------------------------------------------------------------
  // Constant solution
  //------------------------------------------------------------------------------
  static Elasticity::FrictionType
  constant_friction = [](const Face * F, const double & lambda) -> double{
    return 0.0;
  };

  static Elasticity::DisplacementType
  constant_u =   [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
                      return Eigen::Vector3d(1., 1., 1.);
                    };

  static Elasticity::GradientDisplacementType
  constant_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
                        return Eigen::Matrix3d::Zero();
                         };

  static Elasticity::GradientDisplacementType1
  constant_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {    
                        return Eigen::Matrix3d::Zero();
                         };                                                

  static Elasticity::ForcingTermType
  constant_div_eps_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                      return Eigen::Vector3d::Zero(); 
                    };

  static Elasticity::ForcingTermType
  constant_grad_div_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                      return Eigen::Vector3d::Zero();  
                    };

  static Elasticity::ForcingTermType
  constant_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                    return -2 * mu * constant_div_eps_u(mu,lambda,x,T) - lambda * constant_grad_div_u(mu,lambda,x,T); 
                    };

  //------------------------------------------------------------------------------
  // Linear solution
  //------------------------------------------------------------------------------
  static Elasticity::DisplacementType
  linear_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Vector3d {
                      return Eigen::Vector3d(x.x(), -x.y(), x.z());
                    };

  static Elasticity::GradientDisplacementType
  linear_grad_u = [](const Eigen::Vector3d & x, const Cell * T, const double & lambda) -> Eigen::Matrix3d {
                        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
                        M.row(0) << 1, 0, 0;
                        M.row(1) << 0, -1, 0;
                        M.row(2) << 0, 0, 1;
                        
                        return M;
                         };
  static Elasticity::GradientDisplacementType1
  linear_sigma_u = [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) ->  Eigen::Matrix3d {    
                        Eigen::Matrix3d grad = Eigen::Matrix3d::Zero();
                        grad.row(0) << 1, 0, 0;
                        grad.row(1) << 0, -1, 0;
                        grad.row(2) << 0, 0, 1;         

        // All other entries remain zero since other components are 0

                        Eigen::Matrix3d trace_grad = Eigen::Matrix3d::Zero();
                        trace_grad.row(0) << 1, 0, 0;
                        trace_grad.row(1) << 0, 1, 0;
                        trace_grad.row(2) << 0, 0, 1;  

        
                        return 2*mu*grad + lambda*trace_grad;
                        };                      

  static Elasticity::ForcingTermType
  linear_div_eps_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                      return Eigen::Vector3d::Zero(); 
                    };

  static Elasticity::ForcingTermType
  linear_grad_div_u = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                      return Eigen::Vector3d::Zero();  
                    };

  static Elasticity::ForcingTermType
  linear_f = 
          [](const double & mu, const double & lambda, const Eigen::Vector3d & x, const Cell * T) -> Eigen::Vector3d {
                    return -2 * mu * linear_div_eps_u(mu,lambda,x,T) - lambda * linear_grad_div_u(mu,lambda,x,T); 
                    };

} // end of namespace HArDCore3D

#endif
