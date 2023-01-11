#ifndef YANGMILLS_LM_HPP
#define YANGMILLS_LM_HPP

#include <iostream>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include <mesh.hpp>
#include <mesh_builder.hpp>

#include <laxgrad.hpp>
#include <laxcurl.hpp>
#include <laxdiv.hpp>

 /*!
  * @defgroup DDR_yangmills
  * @brief Implementation of the DDR scheme for the %Yang-Mills problem 
  */

namespace HArDCore3D
{

   /*!
    * @addtogroup DDR_yangmills
    * @{
    */
   
  /// Structure to store norm components
 struct YangMillsNorms
 {
   /// Constructor
   YangMillsNorms(double norm_E, double norm_A, double norm_lambda):
    E(norm_E),
    A(norm_A),
    lambda(norm_lambda)
     {
       // do nothing
     };
   
  double E; ///< Norm of electric field
  double A; ///< Norm of potential
  double lambda; ///< Norm of Lagrange multiplier
 };

  /// Assemble a YangMills problem 
  struct YangMills
  {
    typedef Eigen::SparseMatrix<double> SystemMatrixType;
    /// Base functions 
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> ForcingTermType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> ElectricFieldType; 
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> MagneticFieldType;
    /// Base functions with time component
    typedef std::function<Eigen::Vector3d(const double &, const Eigen::Vector3d &)> TForcingTermType;
    typedef std::function<Eigen::Vector3d(const double &, const Eigen::Vector3d &)> TElectricFieldType; 
    typedef std::function<Eigen::Vector3d(const double &, const Eigen::Vector3d &)> TMagneticFieldType;
    /// Lie algebra valued functions (vector of dim Lie algebra)
    typedef std::vector<ForcingTermType> LAForcingTermType;
    typedef std::vector<ElectricFieldType> LAElectricFieldType; 
    typedef std::vector<MagneticFieldType> LAMagneticFieldType;
    /// Lie algebra valued functions with time component
    typedef std::vector<TForcingTermType> TLAForcingTermType;
    typedef std::vector<TElectricFieldType> TLAElectricFieldType;
    typedef std::vector<TMagneticFieldType> TLAMagneticFieldType;

    /// Constructor
    YangMills(
           const DDRCore & ddrcore,           ///< Core for the DDR space sequence
           const LieAlgebra & liealgebra,     ///< Lie algebra of the problem
           bool use_threads,                  ///< True for parallel execution, false for sequential execution
           std::ostream & output = std::cout  ///< Output stream to print status messages
         );

    /// Assemble the global system    
    void assembleLinearSystem(
                              double dt ///< Time step
                              );
    /// Sets system vector for the nonlinear problem
    void setSystemVector(
                          const Eigen::VectorXd & interp_f,  ///< Forcing term of first equation
                          const Eigen::VectorXd & interp_dE, ///< Used in forcing term of second equation
                          const Eigen::VectorXd & interp_A,  ///< Used in forcing term of second equation
                          const Eigen::VectorXd & E_old,    ///< Electric field at previous time
                          const Eigen::VectorXd & A_old,    ///< Potential at previous time
                          double dt,                        ///< Time step
                          double theta,                     ///< Theta scheme
                          double nonlinear_coeff            ///< Scaling of the nonlinear terms (0 for Maxwell)
                        );
    /// Assembles the system for Newton iterations
    void assembleSystemNewton(
                          const Eigen::VectorXd & E_i,    ///< Electric field at previous time
                          const Eigen::VectorXd & A_i,    ///< Potential at previous time
                          const Eigen::VectorXd & Elambdac_k, ///< Solution position after previous Newton iteration
                          const Eigen::VectorXd & sysVec,   ///< System vector of nonlinear problem
                          double dt,                        ///< Time step
                          double theta,                     ///< Theta scheme
                          double nonlinear_coeff            ///< Scaling of the nonlinear terms (0 for Maxwell)
                        );
    /// Stops once changes become small
    double stoppingCrit(
                      const Eigen::VectorXd & v, 
                      const Eigen::VectorXd & u  
                    );
    /// Stops once close enough to system vector for nonlinear problem
    double stoppingCrit2(
                        const Eigen::VectorXd & Elambda_k, ///< Solution in Newton iterations
                        const Eigen::VectorXd & E_i,       ///< Electric field at previous time
                        const Eigen::VectorXd & A_i,       ///< Potential at previous time
                        const Eigen::VectorXd & sysVec,    ///< System vector of nonlinear problem
                        double dt,                         ///< Time step
                        double theta,                      ///< Theta scheme
                        double nonlinear_coeff             ///< Scaling of the nonlinear terms (0 for Maxwell)
                        );
    /// Adds boundary condition for chosen solution to the system vector
    void addBoundaryConditions(
                              const LAMagneticFieldType & curl_A, ///< Curl of the potential
                              double dt                           ///< Time step
                              );
    /// Computes the projected initial conditions that satisfy the constraint
    Eigen::VectorXd computeInitialConditions(const Eigen::MatrixXd & Eh,   ///< Interpolate of the initial electric field
                                             const Eigen::MatrixXd & Ah,   ///< Interpolate of the initial potential field
                                             const double nonlinear_coeff, ///< Scaling of the nonlinear terms (0 for Maxwell)
                                             const size_t solver           ///< Solver to use
                                             );

    /// Calculates the matrix of Hstar[v,.] in element index iT
    Eigen::MatrixXd epsBkt_v(size_t iT, const Eigen::VectorXd & v) const;

    /// Wrapper to plug vector into L2 integral product: int (Pcurl 1,[Pcurl 2, Pgrad 3])
    Eigen::MatrixXd L2v_Bkt(size_t iT,                                    ///< Element index
                            boost::multi_array<double, 3> & intPciPcjPgk, ///< L2 integral product (trilinear form)
                            const Eigen::VectorXd & v,                    ///< Vector to plug in
                            const size_t & entry) const;

    Eigen::MatrixXd L2v_epsBkt(size_t iT,                     ///< Element index
                               const Eigen::VectorXd & v_T,   ///< Local vector in T 
                               const Eigen::MatrixXd & L2prod ///< L2 product to be combined with bracket
                               ) const;

    /// Returns the global problem dimension 
    inline size_t dimension() const
    {
      return m_laxcurl.dimension() + m_laxgrad.dimension();
    }

    /// Returns the size of the system
    inline size_t sizeSystem() const
    {
      return dimension(); 
    }

    /// Returns the Lie algebra
    inline const LieAlgebra & lieAlg() const
    {
      return m_liealg;
    }

    /// Returns the space LAXGrad
    inline const LAXGrad & laXGrad() const
    {
      return m_laxgrad;
    }

    /// Returns the space LAXCurl
    inline const LAXCurl & laXCurl() const
    {
      return m_laxcurl;
    }

    /// Returns the space LAXCurl
    inline const LAXDiv & laXDiv() const
    {
      return m_laxdiv;
    }

    /// Returns the space XDiv
    inline const XGrad & xGrad() const
    {
      return m_xgrad;
    }

    /// Returns the space XDiv
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

    /// Compute the discrete L2 norms, for a family of Eigen::VectorXd representing the electric field and potential
    std::vector<YangMillsNorms> computeYangMillsNorms(
                                                  const std::vector<Eigen::VectorXd> & list_dofs ///< List of vectors representing the electric field and potential
                                                  ) const;

    /// Takes a two parameter function and a value and returns a function with the first parameter fixed to value
    template<typename outValue, typename TFct>
    std::function<outValue(const Eigen::Vector3d &)> contractPara(const TFct &F, double t) const
    {
      std::function<outValue(const Eigen::Vector3d &)> f = [&F, t](const Eigen::Vector3d &x) -> outValue {return F(t, x);};
      return f;
    }
    /// Takes a vector of two parameter functions and a value and returns a function with the first parameter fixed to value
    template<typename outValue, typename Fct, typename TFct>
    std::vector<Fct> contractParaLA(const std::vector<TFct> &TF, double t) const
    { 
      std::vector<Fct> f;
      for (size_t i = 0; i < TF.size(); i++){
        f.emplace_back(contractPara<outValue>(TF[i], t));
      }
      return f;
    }

    /// Calculates residual with current matrix and rhs from given x (m_A*x = m_b)
    double computeResidual(const Eigen::VectorXd & x) const;

    /// Computes the constraint [E, grad P'] +  int <E, [A, P']> for the DOFs
    Eigen::VectorXd computeConstraint(const Eigen::VectorXd & E, const Eigen::VectorXd & A, const double nonlinear_coeff);

    /// Solves a system to calculate the norm of the constraint (in the dual space of LaXgrad)
    double computeConstraintNorm(const Eigen::VectorXd & Ch, const size_t itersolver) const;

		std::pair<double, double> computeConditionNum() const;

  private:
    /// Compute the local contribution for the element of index iT
    std::vector<Eigen::MatrixXd>
    _compute_local_contribution(
                                size_t iT, ///< Element index
                                double dt  ///< Time step
                                );

    /// Assemble the local contribution for the element of index iT into the global system
    void _assemble_local_contribution(
                                      size_t iT,    ///< Element index
                                      const std::vector<Eigen::MatrixXd> & lsT,   ///< Local contributions
                                      std::vector<std::list<Eigen::Triplet<double>>> & triplets,    ///< List of triplets for the system
                                      std::vector<Eigen::VectorXd> & vecs   //< List of vectors for the system
                                      );

    std::vector<Eigen::MatrixXd> _epsBkt();

    /// Computes the Lie bracket of two (constant) vectors contracted with epsilon (antisymmetric symbol)
    Eigen::MatrixXd _compute_det_bracket_face(
                                      const size_t iF ///< index of the face
                                      ) const;

    /// Computes the (constant) face tangent functions on the basis of F
    Eigen::MatrixXd _compute_face_tangent(const size_t iF) const;

    /// Computes integrals of three basis function potentials (Xcurl)x(Xcurl)x(xGrad)
    boost::multi_array<double, 3> _integral_ijk(size_t iT) const;

    /// Creates L2 integral product with bracket: int (P_LaXcurl,[P_LaXcurl, P_LaXgrad])
    boost::multi_array<double, 3> _int_Pci_bktPcjPgk(size_t iT);

    const DDRCore & m_ddrcore;
    bool m_use_threads;
    std::ostream & m_output;
    const XGrad m_xgrad;
    const XCurl m_xcurl;
    const XDiv m_xdiv;
    const LieAlgebra m_liealg;
    const LAXGrad m_laxgrad;
    const LAXCurl m_laxcurl;
    const LAXDiv m_laxdiv;
    SystemMatrixType m_A;   // Matrix and RHS system
    Eigen::VectorXd m_b;
    SystemMatrixType m_laxcurl_L2;   // Global laxcurl L2Product matrix
    SystemMatrixType m_laxdiv_L2;   // Global laxdiv L2Product matrix
    SystemMatrixType m_laxgrad_L2;   // Global laxgrad L2Product matrix
    std::vector<Eigen::MatrixXd> m_F_ebkt;  // Matrix operator on the face to recover (LA) face value of ebkt
    std::vector<boost::multi_array<double, 3>> m_int_PciPcjPgk;   // L2 integral product with bracket (trilinear form)
    double m_stab_par;
  };

    /// Template to evaluate a vector of functions 
    template<typename outValue, typename Fct>
    std::vector<Fct> sumLA(const std::vector<Fct> & LAF, const std::vector<Fct> & LAG, double lambda)
    { 
      std::vector<Fct> values;
      for (size_t i = 0; i < LAF.size(); i++){
        values.emplace_back([&, i, lambda](double t, const Eigen::Vector3d & x)-> outValue {return LAF[i](t, x) + lambda * LAG[i](t, x);});
      }
      return values;
    }

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
  static YangMills::TElectricFieldType
  trigonometric_E1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              -0.5*sin(t)*sin(PI*x(0))*cos(PI*x(1))*cos(PI*x(2)), 
                                              sin(t)*sin(PI*x(1))*cos(PI*x(0))*cos(PI*x(2)), 
                                              -0.5*sin(t)*sin(PI*x(2))*cos(PI*x(0))*cos(PI*x(1))
                                              );
                      };
  static YangMills::TElectricFieldType
  trigonometric_E2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                               0.5*sin(PI*x(0))*cos(t)*cos(PI*x(1))*cos(PI*x(2)), 
                                               -sin(PI*x(1))*cos(t)*cos(PI*x(0))*cos(PI*x(2)), 
                                               0.5*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1))
                                              );
                      }; 

  static YangMills::TElectricFieldType
  trigonometric_E3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                                0.5*pow(sin(PI*x(1)), 2)*cos(t), 
                                               sin(t)*pow(cos(PI*x(2)), 2), 
                                               0.5*cos(t)*pow(cos(PI*x(0)), 2)
                                               );
                      };

  static YangMills::TLAElectricFieldType
  trigonometric_E = {trigonometric_E1, trigonometric_E2, trigonometric_E3};

  static YangMills::TElectricFieldType
  trigonometric_A1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              -0.5*sin(PI*x(0))*cos(t)*cos(PI*x(1))*cos(PI*x(2)), 
                                              sin(PI*x(1))*cos(t)*cos(PI*x(0))*cos(PI*x(2)), 
                                              -0.5*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1))
                                              );
                      };
                      
  static YangMills::TElectricFieldType
  trigonometric_A2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              -0.5*sin(t)*sin(PI*x(0))*cos(PI*x(1))*cos(PI*x(2)), 
                                              sin(t)*sin(PI*x(1))*cos(PI*x(0))*cos(PI*x(2)), 
                                              -0.5*sin(t)*sin(PI*x(2))*cos(PI*x(0))*cos(PI*x(1))
                                              );
                      };

  static YangMills::TElectricFieldType
  trigonometric_A3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(-0.5*sin(t)*pow(sin(PI*x(1)), 2), 
                                                cos(t)*pow(cos(PI*x(2)), 2), 
                                                -0.5*sin(t)*pow(cos(PI*x(0)), 2)
                                                );
                      }; 

  static YangMills::TLAElectricFieldType
  trigonometric_A = {trigonometric_A1, trigonometric_A2, trigonometric_A3};

  static YangMills::TMagneticFieldType
  trigonometric_B1_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    1.5*PI*sin(PI*x(1))*sin(PI*x(2))*cos(t)*cos(PI*x(0)), 
                                                    0, 
                                                    -1.5*PI*sin(PI*x(0))*sin(PI*x(1))*cos(t)*cos(PI*x(2))
                                                    );
                            }; 

  static YangMills::TMagneticFieldType
  trigonometric_B2_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
	                                              		1.5*PI*sin(t)*sin(PI*x(1))*sin(PI*x(2))*cos(PI*x(0)), 
                                                    0, 
                                                    -1.5*PI*sin(t)*sin(PI*x(0))*sin(PI*x(1))*cos(PI*x(2))
	                                              		);
                            }; 

  static YangMills::TMagneticFieldType
  trigonometric_B3_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    2*PI*sin(PI*x(2))*cos(t)*cos(PI*x(2)), 
                                                    -1.0*PI*sin(t)*sin(PI*x(0))*cos(PI*x(0)), 
                                                    1.0*PI*sin(t)*sin(PI*x(1))*cos(PI*x(1))
                                                     );
                            }; 
              
  static YangMills::TLAMagneticFieldType
  trigonometric_B_linear = {trigonometric_B1_linear, trigonometric_B2_linear, trigonometric_B3_linear};

  static YangMills::TMagneticFieldType
  trigonometric_B1_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    -0.5*pow(sin(t), 2)*sin(PI*x(1))*pow(cos(PI*x(0)), 3)*cos(PI*x(2)) + 0.5*sin(t)*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1))*pow(cos(PI*x(2)), 2), 
                                                    -0.25*pow(sin(t), 2)*sin(PI*x(0))*pow(cos(PI*x(0)), 2)*cos(PI*x(1))*cos(PI*x(2)) + 0.25*pow(sin(t), 2)*pow(sin(PI*x(1)), 2)*sin(PI*x(2))*cos(PI*x(0))*cos(PI*x(1)), 
                                                    0.5*pow(sin(t), 2)*pow(sin(PI*x(1)), 3)*cos(PI*x(0))*cos(PI*x(2)) - 0.5*sin(t)*sin(PI*x(0))*cos(t)*cos(PI*x(1))*pow(cos(PI*x(2)), 3)
                                                    );
                            }; 

  static YangMills::TMagneticFieldType
  trigonometric_B2_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
	                                              		0.5*sin(t)*sin(PI*x(1))*cos(t)*pow(cos(PI*x(0)), 3)*cos(PI*x(2)) - 0.5*sin(PI*x(2))*pow(cos(t), 2)*cos(PI*x(0))*cos(PI*x(1))*pow(cos(PI*x(2)), 2), 
                                                    0.25*sin(t)*sin(PI*x(0))*cos(t)*pow(cos(PI*x(0)), 2)*cos(PI*x(1))*cos(PI*x(2)) - 0.25*sin(t)*pow(sin(PI*x(1)), 2)*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1)), 
                                                    -0.5*sin(t)*pow(sin(PI*x(1)), 3)*cos(t)*cos(PI*x(0))*cos(PI*x(2)) + 0.5*sin(PI*x(0))*pow(cos(t), 2)*cos(PI*x(1))*pow(cos(PI*x(2)), 3)
	                                              		);
                            }; 

  static YangMills::TMagneticFieldType
  trigonometric_B3_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                     );
                            }; 
              
  static YangMills::TLAMagneticFieldType
  trigonometric_B_nonlinear = {trigonometric_B1_nonlinear, trigonometric_B2_nonlinear, trigonometric_B3_nonlinear};

  static YangMills::TElectricFieldType
  trigonometric_dtE1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              -0.5*sin(PI*x(0))*cos(t)*cos(PI*x(1))*cos(PI*x(2)), 
                                              sin(PI*x(1))*cos(t)*cos(PI*x(0))*cos(PI*x(2)), 
                                              -0.5*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1))
                                              );
                      };
  static YangMills::TElectricFieldType
  trigonometric_dtE2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                               -0.5*sin(t)*sin(PI*x(0))*cos(PI*x(1))*cos(PI*x(2)), 
                                               sin(t)*sin(PI*x(1))*cos(PI*x(0))*cos(PI*x(2)), 
                                               -0.5*sin(t)*sin(PI*x(2))*cos(PI*x(0))*cos(PI*x(1))
                                              );
                      }; 

  static YangMills::TElectricFieldType
  trigonometric_dtE3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                               -0.5*sin(t)*pow(sin(PI*x(1)), 2), 
                                               cos(t)*pow(cos(PI*x(2)), 2), 
                                               -0.5*sin(t)*pow(cos(PI*x(0)), 2)
                                               );
                      };

  static YangMills::TLAElectricFieldType
  trigonometric_dtE = {trigonometric_dtE1, trigonometric_dtE2, trigonometric_dtE3};

  static YangMills::TForcingTermType
  trigonometric_f1_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                               -0.5*sin(PI*x(0))*cos(t)*cos(PI*x(1))*cos(PI*x(2)) + 1.5*pow(PI, 2)*sin(PI*x(0))*cos(t)*cos(PI*x(1))*cos(PI*x(2)), 
                                               -3.0*pow(PI, 2)*sin(PI*x(1))*cos(t)*cos(PI*x(0))*cos(PI*x(2)) + sin(PI*x(1))*cos(t)*cos(PI*x(0))*cos(PI*x(2)), 
                                               -0.5*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1)) + 1.5*pow(PI, 2)*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1))
                                               );
                      };

  static YangMills::TForcingTermType
  trigonometric_f2_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              -0.5*sin(t)*sin(PI*x(0))*cos(PI*x(1))*cos(PI*x(2)) + 1.5*pow(PI, 2)*sin(t)*sin(PI*x(0))*cos(PI*x(1))*cos(PI*x(2)), 
                                              -3.0*pow(PI, 2)*sin(t)*sin(PI*x(1))*cos(PI*x(0))*cos(PI*x(2)) + sin(t)*sin(PI*x(1))*cos(PI*x(0))*cos(PI*x(2)), 
                                              -0.5*sin(t)*sin(PI*x(2))*cos(PI*x(0))*cos(PI*x(1)) + 1.5*pow(PI, 2)*sin(t)*sin(PI*x(2))*cos(PI*x(0))*cos(PI*x(1))
                                              );
                      };

  static YangMills::TForcingTermType
  trigonometric_f3_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                               -0.5*sin(t)*pow(sin(PI*x(1)), 2) + 1.0*pow(PI, 2)*sin(t)*pow(sin(PI*x(1)), 2) - 1.0*pow(PI, 2)*sin(t)*pow(cos(PI*x(1)), 2), 
                                               2*pow(PI, 2)*pow(sin(PI*x(2)), 2)*cos(t) - 2*pow(PI, 2)*cos(t)*pow(cos(PI*x(2)), 2) + cos(t)*pow(cos(PI*x(2)), 2), 
                                               -1.0*pow(PI, 2)*sin(t)*pow(sin(PI*x(0)), 2) - 0.5*sin(t)*pow(cos(PI*x(0)), 2) + 1.0*pow(PI, 2)*sin(t)*pow(cos(PI*x(0)), 2)
                                               );
                      };

  static YangMills::TLAForcingTermType
  trigonometric_f_linear = {trigonometric_f1_linear, trigonometric_f2_linear, trigonometric_f3_linear};

  static YangMills::TForcingTermType
  trigonometric_f1_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              0.5*(0.25*sin(t)*sin(PI*x(0))*cos(t)*pow(cos(PI*x(0)), 2)*cos(PI*x(1))*cos(PI*x(2)) - 0.25*sin(t)*pow(sin(PI*x(1)), 2)*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1)))*sin(t)*pow(cos(PI*x(0)), 2) + (-1.5*PI*sin(t)*sin(PI*x(0))*sin(PI*x(1))*cos(PI*x(2)) - 0.5*sin(t)*pow(sin(PI*x(1)), 3)*cos(t)*cos(PI*x(0))*cos(PI*x(2)) + 0.5*sin(PI*x(0))*pow(cos(t), 2)*cos(PI*x(1))*pow(cos(PI*x(2)), 3))*cos(t)*pow(cos(PI*x(2)), 2) + 0.75*PI*pow(sin(t), 2)*sin(PI*x(0))*sin(PI*x(2))*pow(cos(PI*x(0)), 2)*cos(PI*x(1)) - 2.25*PI*pow(sin(t), 2)*pow(sin(PI*x(1)), 2)*cos(PI*x(0))*cos(PI*x(1))*cos(PI*x(2)) - 0.5*PI*sin(t)*sin(PI*x(0))*sin(PI*x(1))*cos(t)*pow(cos(PI*x(2)), 3), 
                                              0.5*(-1.5*PI*sin(t)*sin(PI*x(0))*sin(PI*x(1))*cos(PI*x(2)) - 0.5*sin(t)*pow(sin(PI*x(1)), 3)*cos(t)*cos(PI*x(0))*cos(PI*x(2)) + 0.5*sin(PI*x(0))*pow(cos(t), 2)*cos(PI*x(1))*pow(cos(PI*x(2)), 3))*sin(t)*pow(sin(PI*x(1)), 2) - 0.5*(1.5*PI*sin(t)*sin(PI*x(1))*sin(PI*x(2))*cos(PI*x(0)) + 0.5*sin(t)*sin(PI*x(1))*cos(t)*pow(cos(PI*x(0)), 3)*cos(PI*x(2)) - 0.5*sin(PI*x(2))*pow(cos(t), 2)*cos(PI*x(0))*cos(PI*x(1))*pow(cos(PI*x(2)), 2))*sin(t)*pow(cos(PI*x(0)), 2) - 0.5*PI*pow(sin(t), 2)*sin(PI*x(0))*pow(sin(PI*x(1)), 3)*cos(PI*x(2)) - 0.5*PI*pow(sin(t), 2)*sin(PI*x(0))*sin(PI*x(1))*pow(cos(PI*x(1)), 2)*cos(PI*x(2)) - 0.5*PI*pow(sin(t), 2)*sin(PI*x(1))*sin(PI*x(2))*pow(cos(PI*x(0)), 3) + 2.0*PI*sin(t)*pow(sin(PI*x(2)), 2)*cos(t)*cos(PI*x(0))*cos(PI*x(1))*cos(PI*x(2)) - 1.0*PI*sin(t)*cos(t)*cos(PI*x(0))*cos(PI*x(1))*pow(cos(PI*x(2)), 3), 
                                              -0.5*(0.25*sin(t)*sin(PI*x(0))*cos(t)*pow(cos(PI*x(0)), 2)*cos(PI*x(1))*cos(PI*x(2)) - 0.25*sin(t)*pow(sin(PI*x(1)), 2)*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1)))*sin(t)*pow(sin(PI*x(1)), 2) - (1.5*PI*sin(t)*sin(PI*x(1))*sin(PI*x(2))*cos(PI*x(0)) + 0.5*sin(t)*sin(PI*x(1))*cos(t)*pow(cos(PI*x(0)), 3)*cos(PI*x(2)) - 0.5*sin(PI*x(2))*pow(cos(t), 2)*cos(PI*x(0))*cos(PI*x(1))*pow(cos(PI*x(2)), 2))*cos(t)*pow(cos(PI*x(2)), 2) - 1.0*PI*pow(sin(t), 2)*pow(sin(PI*x(0)), 2)*cos(PI*x(0))*cos(PI*x(1))*cos(PI*x(2)) + 0.25*PI*pow(sin(t), 2)*sin(PI*x(0))*pow(sin(PI*x(1)), 2)*sin(PI*x(2))*cos(PI*x(1)) - 0.25*PI*pow(sin(t), 2)*pow(cos(PI*x(0)), 3)*cos(PI*x(1))*cos(PI*x(2)) + 1.5*PI*sin(t)*sin(PI*x(1))*sin(PI*x(2))*cos(t)*cos(PI*x(0))*pow(cos(PI*x(2)), 2)
                                              );
                      };

  static YangMills::TForcingTermType
  trigonometric_f2_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              -0.5*(-0.25*pow(sin(t), 2)*sin(PI*x(0))*pow(cos(PI*x(0)), 2)*cos(PI*x(1))*cos(PI*x(2)) + 0.25*pow(sin(t), 2)*pow(sin(PI*x(1)), 2)*sin(PI*x(2))*cos(PI*x(0))*cos(PI*x(1)))*sin(t)*pow(cos(PI*x(0)), 2) - (0.5*pow(sin(t), 2)*pow(sin(PI*x(1)), 3)*cos(PI*x(0))*cos(PI*x(2)) - 0.5*sin(t)*sin(PI*x(0))*cos(t)*cos(PI*x(1))*pow(cos(PI*x(2)), 3) - 1.5*PI*sin(PI*x(0))*sin(PI*x(1))*cos(t)*cos(PI*x(2)))*cos(t)*pow(cos(PI*x(2)), 2) - 0.75*PI*sin(t)*sin(PI*x(0))*sin(PI*x(2))*cos(t)*pow(cos(PI*x(0)), 2)*cos(PI*x(1)) + 2.25*PI*sin(t)*pow(sin(PI*x(1)), 2)*cos(t)*cos(PI*x(0))*cos(PI*x(1))*cos(PI*x(2)) + 0.5*PI*sin(PI*x(0))*sin(PI*x(1))*pow(cos(t), 2)*pow(cos(PI*x(2)), 3), 
                                              0.5*(-0.5*pow(sin(t), 2)*sin(PI*x(1))*pow(cos(PI*x(0)), 3)*cos(PI*x(2)) + 0.5*sin(t)*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1))*pow(cos(PI*x(2)), 2) + 1.5*PI*sin(PI*x(1))*sin(PI*x(2))*cos(t)*cos(PI*x(0)))*sin(t)*pow(cos(PI*x(0)), 2) - 0.5*(0.5*pow(sin(t), 2)*pow(sin(PI*x(1)), 3)*cos(PI*x(0))*cos(PI*x(2)) - 0.5*sin(t)*sin(PI*x(0))*cos(t)*cos(PI*x(1))*pow(cos(PI*x(2)), 3) - 1.5*PI*sin(PI*x(0))*sin(PI*x(1))*cos(t)*cos(PI*x(2)))*sin(t)*pow(sin(PI*x(1)), 2) + 0.5*PI*sin(t)*sin(PI*x(0))*pow(sin(PI*x(1)), 3)*cos(t)*cos(PI*x(2)) + 0.5*PI*sin(t)*sin(PI*x(0))*sin(PI*x(1))*cos(t)*pow(cos(PI*x(1)), 2)*cos(PI*x(2)) + 0.5*PI*sin(t)*sin(PI*x(1))*sin(PI*x(2))*cos(t)*pow(cos(PI*x(0)), 3) - 2.0*PI*pow(sin(PI*x(2)), 2)*pow(cos(t), 2)*cos(PI*x(0))*cos(PI*x(1))*cos(PI*x(2)) + 1.0*PI*pow(cos(t), 2)*cos(PI*x(0))*cos(PI*x(1))*pow(cos(PI*x(2)), 3), 
                                              0.5*(-0.25*pow(sin(t), 2)*sin(PI*x(0))*pow(cos(PI*x(0)), 2)*cos(PI*x(1))*cos(PI*x(2)) + 0.25*pow(sin(t), 2)*pow(sin(PI*x(1)), 2)*sin(PI*x(2))*cos(PI*x(0))*cos(PI*x(1)))*sin(t)*pow(sin(PI*x(1)), 2) + (-0.5*pow(sin(t), 2)*sin(PI*x(1))*pow(cos(PI*x(0)), 3)*cos(PI*x(2)) + 0.5*sin(t)*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1))*pow(cos(PI*x(2)), 2) + 1.5*PI*sin(PI*x(1))*sin(PI*x(2))*cos(t)*cos(PI*x(0)))*cos(t)*pow(cos(PI*x(2)), 2) + 1.0*PI*sin(t)*pow(sin(PI*x(0)), 2)*cos(t)*cos(PI*x(0))*cos(PI*x(1))*cos(PI*x(2)) - 0.25*PI*sin(t)*sin(PI*x(0))*pow(sin(PI*x(1)), 2)*sin(PI*x(2))*cos(t)*cos(PI*x(1)) + 0.25*PI*sin(t)*cos(t)*pow(cos(PI*x(0)), 3)*cos(PI*x(1))*cos(PI*x(2)) - 1.5*PI*sin(PI*x(1))*sin(PI*x(2))*pow(cos(t), 2)*cos(PI*x(0))*pow(cos(PI*x(2)), 2)
                                              );
                      };

  static YangMills::TForcingTermType
  trigonometric_f3_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                               0.5*(-0.25*pow(sin(t), 2)*sin(PI*x(0))*pow(cos(PI*x(0)), 2)*cos(PI*x(1))*cos(PI*x(2)) + 0.25*pow(sin(t), 2)*pow(sin(PI*x(1)), 2)*sin(PI*x(2))*cos(PI*x(0))*cos(PI*x(1)))*sin(t)*sin(PI*x(2))*cos(PI*x(0))*cos(PI*x(1)) - 0.5*(0.25*sin(t)*sin(PI*x(0))*cos(t)*pow(cos(PI*x(0)), 2)*cos(PI*x(1))*cos(PI*x(2)) - 0.25*sin(t)*pow(sin(PI*x(1)), 2)*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1)))*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1)) + (0.5*pow(sin(t), 2)*pow(sin(PI*x(1)), 3)*cos(PI*x(0))*cos(PI*x(2)) - 0.5*sin(t)*sin(PI*x(0))*cos(t)*cos(PI*x(1))*pow(cos(PI*x(2)), 3) - 1.5*PI*sin(PI*x(0))*sin(PI*x(1))*cos(t)*cos(PI*x(2)))*sin(t)*sin(PI*x(1))*cos(PI*x(0))*cos(PI*x(2)) - (-1.5*PI*sin(t)*sin(PI*x(0))*sin(PI*x(1))*cos(PI*x(2)) - 0.5*sin(t)*pow(sin(PI*x(1)), 3)*cos(t)*cos(PI*x(0))*cos(PI*x(2)) + 0.5*sin(PI*x(0))*pow(cos(t), 2)*cos(PI*x(1))*pow(cos(PI*x(2)), 3))*sin(PI*x(1))*cos(t)*cos(PI*x(0))*cos(PI*x(2)), 
                                               -0.5*(-0.5*pow(sin(t), 2)*sin(PI*x(1))*pow(cos(PI*x(0)), 3)*cos(PI*x(2)) + 0.5*sin(t)*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1))*pow(cos(PI*x(2)), 2) + 1.5*PI*sin(PI*x(1))*sin(PI*x(2))*cos(t)*cos(PI*x(0)))*sin(t)*sin(PI*x(2))*cos(PI*x(0))*cos(PI*x(1)) + 0.5*(0.5*pow(sin(t), 2)*pow(sin(PI*x(1)), 3)*cos(PI*x(0))*cos(PI*x(2)) - 0.5*sin(t)*sin(PI*x(0))*cos(t)*cos(PI*x(1))*pow(cos(PI*x(2)), 3) - 1.5*PI*sin(PI*x(0))*sin(PI*x(1))*cos(t)*cos(PI*x(2)))*sin(t)*sin(PI*x(0))*cos(PI*x(1))*cos(PI*x(2)) - 0.5*(-1.5*PI*sin(t)*sin(PI*x(0))*sin(PI*x(1))*cos(PI*x(2)) - 0.5*sin(t)*pow(sin(PI*x(1)), 3)*cos(t)*cos(PI*x(0))*cos(PI*x(2)) + 0.5*sin(PI*x(0))*pow(cos(t), 2)*cos(PI*x(1))*pow(cos(PI*x(2)), 3))*sin(PI*x(0))*cos(t)*cos(PI*x(1))*cos(PI*x(2)) + 0.5*(1.5*PI*sin(t)*sin(PI*x(1))*sin(PI*x(2))*cos(PI*x(0)) + 0.5*sin(t)*sin(PI*x(1))*cos(t)*pow(cos(PI*x(0)), 3)*cos(PI*x(2)) - 0.5*sin(PI*x(2))*pow(cos(t), 2)*cos(PI*x(0))*cos(PI*x(1))*pow(cos(PI*x(2)), 2))*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1)), 
                                               -0.5*(-0.25*pow(sin(t), 2)*sin(PI*x(0))*pow(cos(PI*x(0)), 2)*cos(PI*x(1))*cos(PI*x(2)) + 0.25*pow(sin(t), 2)*pow(sin(PI*x(1)), 2)*sin(PI*x(2))*cos(PI*x(0))*cos(PI*x(1)))*sin(t)*sin(PI*x(0))*cos(PI*x(1))*cos(PI*x(2)) + 0.5*(0.25*sin(t)*sin(PI*x(0))*cos(t)*pow(cos(PI*x(0)), 2)*cos(PI*x(1))*cos(PI*x(2)) - 0.25*sin(t)*pow(sin(PI*x(1)), 2)*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1)))*sin(PI*x(0))*cos(t)*cos(PI*x(1))*cos(PI*x(2)) - (-0.5*pow(sin(t), 2)*sin(PI*x(1))*pow(cos(PI*x(0)), 3)*cos(PI*x(2)) + 0.5*sin(t)*sin(PI*x(2))*cos(t)*cos(PI*x(0))*cos(PI*x(1))*pow(cos(PI*x(2)), 2) + 1.5*PI*sin(PI*x(1))*sin(PI*x(2))*cos(t)*cos(PI*x(0)))*sin(t)*sin(PI*x(1))*cos(PI*x(0))*cos(PI*x(2)) + (1.5*PI*sin(t)*sin(PI*x(1))*sin(PI*x(2))*cos(PI*x(0)) + 0.5*sin(t)*sin(PI*x(1))*cos(t)*pow(cos(PI*x(0)), 3)*cos(PI*x(2)) - 0.5*sin(PI*x(2))*pow(cos(t), 2)*cos(PI*x(0))*cos(PI*x(1))*pow(cos(PI*x(2)), 2))*sin(PI*x(1))*cos(t)*cos(PI*x(0))*cos(PI*x(2)));
                      };

  static YangMills::TLAForcingTermType
  trigonometric_f_nonlinear = {trigonometric_f1_nonlinear, trigonometric_f2_nonlinear, trigonometric_f3_nonlinear};

  //------------------------------------------------------------------------------
  // Linear solution
  //------------------------------------------------------------------------------

  static YangMills::TElectricFieldType
  linear_E1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(-x(2), -x(1), -x(0));
             };

  static YangMills::TElectricFieldType
  linear_E2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(-x(2), -x(0), -x(1));
             };

  static YangMills::TElectricFieldType
  linear_E3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(-x(1) - x(2), -x(0) - x(2), -x(0) - x(1));
             };

  static YangMills::TLAElectricFieldType
  linear_E = {linear_E1, linear_E2, linear_E3};

  static YangMills::TElectricFieldType
  linear_A1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(t*x(2), t*x(1), t*x(0));
                  };

  static YangMills::TElectricFieldType
  linear_A2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(t*x(2), t*x(0), t*x(1));
                  };

  static YangMills::TElectricFieldType
  linear_A3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(t*x(1) + t*x(2), t*x(0) + t*x(2), t*x(0) + t*x(1));
                  };

  static YangMills::TLAElectricFieldType
  linear_A = {linear_A1, linear_A2, linear_A3};

  static YangMills::TMagneticFieldType
  linear_B1_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(
                                      0, 
                                      0, 
                                      0
                                      );
             };

  static YangMills::TMagneticFieldType
  linear_B2_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(
                                      t, 
                                      t, 
                                      t
                                      );
             };

  static YangMills::TMagneticFieldType
  linear_B3_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(
                                      0,
                                      0,
                                      0
                                      );
             };

  static YangMills::TLAMagneticFieldType
  linear_B_linear = {linear_B1_linear, linear_B2_linear, linear_B3_linear};

  static YangMills::TMagneticFieldType
  linear_B1_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(
                                      1.0*t*x(0)*(t*x(0) + t*x(1)) - 1.0*t*x(1)*(t*x(0) + t*x(2)), 
                                      1.0*t*x(1)*(t*x(1) + t*x(2)) - 1.0*t*x(2)*(t*x(0) + t*x(1)), 
                                      -1.0*t*x(0)*(t*x(1) + t*x(2)) + 1.0*t*x(2)*(t*x(0) + t*x(2))
                                      );
             };

  static YangMills::TMagneticFieldType
  linear_B2_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(
                                      1.0*t*x(0)*(t*x(0) + t*x(2)) - 1.0*t*x(1)*(t*x(0) + t*x(1)), 
                                      -1.0*t*x(0)*(t*x(1) + t*x(2)) + 1.0*t*x(2)*(t*x(0) + t*x(1)), 
                                      1.0*t*x(1)*(t*x(1) + t*x(2)) - 1.0*t*x(2)*(t*x(0) + t*x(2))
                                      );
             };

  static YangMills::TMagneticFieldType
  linear_B3_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(
                                      -1.0*pow(t, 2)*pow(x(0), 2) + 1.0*pow(t, 2)*pow(x(1), 2), 
                                      1.0*pow(t, 2)*x(0)*x(2) - 1.0*pow(t, 2)*x(1)*x(2), 
                                      1.0*pow(t, 2)*x(0)*x(2) - 1.0*pow(t, 2)*x(1)*x(2)
                                      );
             };

  static YangMills::TLAMagneticFieldType
  linear_B_nonlinear = {linear_B1_nonlinear, linear_B2_nonlinear, linear_B3_nonlinear};

  static YangMills::TElectricFieldType
  linear_dtE1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, 0, 0);
                      };
  static YangMills::TElectricFieldType
  linear_dtE2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, 0, 0);
                      }; 

  static YangMills::TElectricFieldType
  linear_dtE3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, 0, 0);
                      };

  static YangMills::TLAElectricFieldType
  linear_dtE = {linear_dtE1, linear_dtE2, linear_dtE3};

  static YangMills::TForcingTermType
  linear_f1_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(0, 0, 0);
                  };

  static YangMills::TForcingTermType
  linear_f2_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(0, 0, 0);
                  };

  static YangMills::TForcingTermType
  linear_f3_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(0, 0, 0);
                  };

  static YangMills::TLAForcingTermType
  linear_f_linear = {linear_f1_linear, linear_f2_linear, linear_f3_linear};

  static YangMills::TForcingTermType
  linear_f1_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(1.0*pow(t, 2)*x(0) + 1.0*pow(t, 2)*x(1) - t*x(0)*(1.0*pow(t, 2)*x(0)*x(2) - 1.0*pow(t, 2)*x(1)*x(2)) + t*x(1)*(1.0*pow(t, 2)*x(0)*x(2) - 1.0*pow(t, 2)*x(1)*x(2)) - 1.0*t*(t*x(0) + t*x(1)) - (t*x(0) + t*x(1))*(-1.0*t*x(0)*(t*x(1) + t*x(2)) + 1.0*t*x(2)*(t*x(0) + t*x(1)) + t) + (t*x(0) + t*x(2))*(1.0*t*x(1)*(t*x(1) + t*x(2)) - 1.0*t*x(2)*(t*x(0) + t*x(2)) + t), 1.0*pow(t, 2)*x(1) + 1.0*pow(t, 2)*x(2) - t*x(1)*(-1.0*pow(t, 2)*pow(x(0), 2) + 1.0*pow(t, 2)*pow(x(1), 2)) + t*x(2)*(1.0*pow(t, 2)*x(0)*x(2) - 1.0*pow(t, 2)*x(1)*x(2)) - 1.0*t*(t*x(1) + t*x(2)) + (t*x(0) + t*x(1))*(1.0*t*x(0)*(t*x(0) + t*x(2)) - 1.0*t*x(1)*(t*x(0) + t*x(1)) + t) - (t*x(1) + t*x(2))*(1.0*t*x(1)*(t*x(1) + t*x(2)) - 1.0*t*x(2)*(t*x(0) + t*x(2)) + t), 1.0*pow(t, 2)*x(0) + 1.0*pow(t, 2)*x(2) + t*x(0)*(-1.0*pow(t, 2)*pow(x(0), 2) + 1.0*pow(t, 2)*pow(x(1), 2)) - t*x(2)*(1.0*pow(t, 2)*x(0)*x(2) - 1.0*pow(t, 2)*x(1)*x(2)) - 1.0*t*(t*x(0) + t*x(2)) - (t*x(0) + t*x(2))*(1.0*t*x(0)*(t*x(0) + t*x(2)) - 1.0*t*x(1)*(t*x(0) + t*x(1)) + t) + (t*x(1) + t*x(2))*(-1.0*t*x(0)*(t*x(1) + t*x(2)) + 1.0*t*x(2)*(t*x(0) + t*x(1)) + t)
                                          );
                  };

  static YangMills::TForcingTermType
  linear_f2_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(-1.0*pow(t, 2)*x(0) - 1.0*pow(t, 2)*x(1) - t*x(0)*(1.0*pow(t, 2)*x(0)*x(2) - 1.0*pow(t, 2)*x(1)*x(2)) + t*x(1)*(1.0*pow(t, 2)*x(0)*x(2) - 1.0*pow(t, 2)*x(1)*x(2)) + 1.0*t*(t*x(0) + t*x(1)) - 1.0*t*(t*x(1) + t*x(2)) + (t*x(0) + t*x(1))*(1.0*t*x(1)*(t*x(1) + t*x(2)) - 1.0*t*x(2)*(t*x(0) + t*x(1))) - (t*x(0) + t*x(2))*(-1.0*t*x(0)*(t*x(1) + t*x(2)) + 1.0*t*x(2)*(t*x(0) + t*x(2))), -1.0*pow(t, 2)*x(0) - 1.0*pow(t, 2)*x(2) + t*x(0)*(-1.0*pow(t, 2)*pow(x(0), 2) + 1.0*pow(t, 2)*pow(x(1), 2)) - t*x(2)*(1.0*pow(t, 2)*x(0)*x(2) - 1.0*pow(t, 2)*x(1)*x(2)) - (t*x(0) + t*x(1))*(1.0*t*x(0)*(t*x(0) + t*x(1)) - 1.0*t*x(1)*(t*x(0) + t*x(2))) + (t*x(1) + t*x(2))*(-1.0*t*x(0)*(t*x(1) + t*x(2)) + 1.0*t*x(2)*(t*x(0) + t*x(2))), -1.0*pow(t, 2)*x(1) - 1.0*pow(t, 2)*x(2) - t*x(1)*(-1.0*pow(t, 2)*pow(x(0), 2) + 1.0*pow(t, 2)*pow(x(1), 2)) + t*x(2)*(1.0*pow(t, 2)*x(0)*x(2) - 1.0*pow(t, 2)*x(1)*x(2)) - 1.0*t*(t*x(0) + t*x(1)) + 1.0*t*(t*x(1) + t*x(2)) + (t*x(0) + t*x(2))*(1.0*t*x(0)*(t*x(0) + t*x(1)) - 1.0*t*x(1)*(t*x(0) + t*x(2))) - (t*x(1) + t*x(2))*(1.0*t*x(1)*(t*x(1) + t*x(2)) - 1.0*t*x(2)*(t*x(0) + t*x(1)))
                                           );
                  };

  static YangMills::TForcingTermType
  linear_f3_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(1.0*pow(t, 2)*x(0) - 1.0*pow(t, 2)*x(1) + 1.0*pow(t, 2)*x(2) + t*x(0)*(-1.0*t*x(0)*(t*x(1) + t*x(2)) + 1.0*t*x(2)*(t*x(0) + t*x(2))) + t*x(0)*(-1.0*t*x(0)*(t*x(1) + t*x(2)) + 1.0*t*x(2)*(t*x(0) + t*x(1)) + t) - t*x(1)*(1.0*t*x(1)*(t*x(1) + t*x(2)) - 1.0*t*x(2)*(t*x(0) + t*x(1))) - t*x(1)*(1.0*t*x(1)*(t*x(1) + t*x(2)) - 1.0*t*x(2)*(t*x(0) + t*x(2)) + t), 1.0*pow(t, 2)*x(2) - t*x(0)*(1.0*t*x(0)*(t*x(0) + t*x(2)) - 1.0*t*x(1)*(t*x(0) + t*x(1)) + t) + t*x(1)*(1.0*t*x(0)*(t*x(0) + t*x(1)) - 1.0*t*x(1)*(t*x(0) + t*x(2))) - t*x(2)*(-1.0*t*x(0)*(t*x(1) + t*x(2)) + 1.0*t*x(2)*(t*x(0) + t*x(2))) + t*x(2)*(1.0*t*x(1)*(t*x(1) + t*x(2)) - 1.0*t*x(2)*(t*x(0) + t*x(2)) + t), 2.0*pow(t, 2)*x(1) - 1.0*pow(t, 2)*x(2) - t*x(0)*(1.0*t*x(0)*(t*x(0) + t*x(1)) - 1.0*t*x(1)*(t*x(0) + t*x(2))) + t*x(1)*(1.0*t*x(0)*(t*x(0) + t*x(2)) - 1.0*t*x(1)*(t*x(0) + t*x(1)) + t) + t*x(2)*(1.0*t*x(1)*(t*x(1) + t*x(2)) - 1.0*t*x(2)*(t*x(0) + t*x(1))) - t*x(2)*(-1.0*t*x(0)*(t*x(1) + t*x(2)) + 1.0*t*x(2)*(t*x(0) + t*x(1)) + t)
                                           );
                  };

  static YangMills::TLAForcingTermType
  linear_f_nonlinear = {linear_f1_nonlinear, linear_f2_nonlinear, linear_f3_nonlinear};

  //------------------------------------------------------------------------------
  // Constant solution
  //------------------------------------------------------------------------------

  static YangMills::TElectricFieldType
  const_E1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d::Zero();
             };

  static YangMills::TElectricFieldType
  const_E2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d::Zero();
             };

  static YangMills::TElectricFieldType
  const_E3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d::Zero();
             };

  static YangMills::TLAElectricFieldType
  const_E = {const_E1, const_E2, const_E3};

  static YangMills::TElectricFieldType
  const_A1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(0.2, 0.4, 0.6);
                  };

  static YangMills::TElectricFieldType
  const_A2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(0.8, 1.0, 1.2);
                  };

  static YangMills::TElectricFieldType
  const_A3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(1.4, 1.6, 1.8);
                  };

  static YangMills::TLAElectricFieldType
  const_A = {const_A1, const_A2, const_A3};

  static YangMills::TMagneticFieldType
  const_B1_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(0, 0, 0);
             };

  static YangMills::TMagneticFieldType
  const_B2_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(0, 0, 0);
             };

  static YangMills::TMagneticFieldType
  const_B3_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(0, 0, 0);
             };

  static YangMills::TLAMagneticFieldType
  const_B_linear = {const_B1_linear, const_B2_linear, const_B3_linear};

  static YangMills::TMagneticFieldType
  const_B1_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(-0.12, 0.24, -0.12);
             };

  static YangMills::TMagneticFieldType
  const_B2_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(0.24, -0.48, 0.24);
             };

  static YangMills::TMagneticFieldType
  const_B3_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
               return Eigen::Vector3d(-0.12, 0.24, -0.12);
             };

  static YangMills::TLAMagneticFieldType
  const_B_nonlinear = {const_B1_nonlinear, const_B2_nonlinear, const_B3_nonlinear};

  static YangMills::TElectricFieldType
  const_dtE1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, 0, 0);
                      };
  static YangMills::TElectricFieldType
  const_dtE2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, 0, 0);
                      }; 

  static YangMills::TElectricFieldType
  const_dtE3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, 0, 0);
                      };

  static YangMills::TLAElectricFieldType
  const_dtE = {const_dtE1, const_dtE2, const_dtE3};

  static YangMills::TForcingTermType
  const_f1_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(0, 0, 0);
                  };

  static YangMills::TForcingTermType
  const_f2_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(0, 0, 0);
                  };

  static YangMills::TForcingTermType
  const_f3_linear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(0, 0, 0);
                  };

  static YangMills::TLAForcingTermType
  const_f_linear = {const_f1_linear, const_f2_linear, const_f3_linear};

  static YangMills::TForcingTermType
  const_f1_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(1.656, 0.144, -1.368);
                  };

  static YangMills::TForcingTermType
  const_f2_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(0.432, 0, -0.432);
                  };

  static YangMills::TForcingTermType
  const_f3_nonlinear = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                    return Eigen::Vector3d(-0.792, -0.144, 0.504);
                  };

  static YangMills::TLAForcingTermType
  const_f_nonlinear = {const_f1_nonlinear, const_f2_nonlinear, const_f3_nonlinear};

  //------------------------------------------------------------------------------
  /// Struct to store the systems that need assembling
  struct AssembleSystems
  {
    /// Constructor
    AssembleSystems(std::vector<Eigen::MatrixXd> sys, std::vector<Eigen::VectorXd> vec):
      systems(sys),
      vectors(vec)
      {
      };
    
    std::vector<Eigen::MatrixXd> systems;
    std::vector<Eigen::VectorXd> vectors;
  };
  /// Struct to store the assembled systems
  struct AssembledSystems
  {
    /// Constructor
    AssembledSystems(std::vector<Eigen::SparseMatrix<double>> sys, std::vector<Eigen::VectorXd> vec):
      systems(sys),
      vectors(vec)
      {
      };
    
    std::vector<Eigen::SparseMatrix<double>> systems;
    std::vector<Eigen::VectorXd> vectors;
  };

  /// Function to assemble global matrices from a procedure that compute local triplets
  static inline AssembledSystems
         parallel_assembly_system(
            size_t nb_elements,    //< nb of elements over which the threading will be done
            std::vector<std::pair<size_t, size_t>> size_systems,    //< sizes of each system to assemble
            std::vector<size_t> size_vectors,   //< sizes of each vector to assemble
            std::function<void(size_t start, size_t end, std::vector<std::list<Eigen::Triplet<double>>> * triplets, std::vector<Eigen::VectorXd> * vecs)> batch_local_assembly, //< procedure to compute all the local contributions to the matrices (put in the triplets) and vectors between start and end points (e.g. elements indices) 
            bool use_threads = true   //< determine if threaded process is used or not
            )
  {
    // Matrices and vectors
    std::vector<Eigen::SparseMatrix<double>> systems;
    systems.reserve(size_systems.size());
    std::vector<Eigen::VectorXd> vectors;
    for (auto size : size_vectors){
      vectors.emplace_back(Eigen::VectorXd::Zero(size));
    }

    if (use_threads) {
      // Select the number of threads
      unsigned nb_threads_hint = std::thread::hardware_concurrency();
      unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);

      // Compute the batch size and the remainder
      unsigned batch_size = nb_elements / nb_threads;
      unsigned batch_remainder = nb_elements % nb_threads;

      // Create vectors of triplets and vectors
      std::vector<std::vector<std::list<Eigen::Triplet<double> > > > triplets(nb_threads + 1, std::vector<std::list<Eigen::Triplet<double> > >(size_systems.size()));
      std::vector<std::vector<Eigen::VectorXd>> vecs(nb_threads + 1, vectors);

      // Assign a task to each thread
      std::vector<std::thread> my_threads(nb_threads);
      for (unsigned i = 0; i < nb_threads; ++i) {
        int start = i * batch_size;
        my_threads[i] = std::thread(batch_local_assembly, start, start + batch_size, &triplets[i], &vecs[i]);
      }

      // Execute the elements left
      int start = nb_threads * batch_size;
      batch_local_assembly(start, start + batch_remainder, &triplets[nb_threads], &vecs[nb_threads]);

      // Wait for the other threads to finish their task
      std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));

      // Create systems from triplets
      for (size_t i = 0; i < size_systems.size(); i++){
        systems.emplace_back(Eigen::SparseMatrix<double>(size_systems[i].first, size_systems[i].second));
        size_t n_triplets = 0;
        for (auto triplets_thread : triplets) {
          n_triplets += triplets_thread[i].size();
        }
        std::vector<Eigen::Triplet<double> > all_triplets(n_triplets);
        auto triplet_index = all_triplets.begin();
        for (auto triplets_thread : triplets) {
          triplet_index = std::copy(triplets_thread[i].begin(), triplets_thread[i].end(), triplet_index);
        }
        systems[i].setFromTriplets(all_triplets.begin(), all_triplets.end());
      }

      for (size_t i = 0; i < size_vectors.size(); i++){
        for (auto vec_thread : vecs){
          vectors[i] += vec_thread[i];
        }
      }

    } else {
      std::vector<std::list<Eigen::Triplet<double> > > triplets(size_systems.size());
      batch_local_assembly(0, nb_elements, &triplets, &vectors);
      for (size_t i = 0; i < size_systems.size(); i++){
        systems.emplace_back(Eigen::SparseMatrix<double>(size_systems[i].first, size_systems[i].second));
        systems[i].setFromTriplets(triplets[i].begin(), triplets[i].end());
      }
    }
    return AssembledSystems(systems, vectors);  
  }
  
} // end of namespace HArDCore3D

#endif
