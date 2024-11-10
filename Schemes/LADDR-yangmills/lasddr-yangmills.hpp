#ifndef SDDR_YANGMILLS_HPP
#define SDDR_YANGMILLS_HPP

#include <iostream>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <local_static_condensation.hpp>

#include <mesh.hpp>
#include <mesh_builder.hpp>
#include <GMpoly_cell.hpp>
#include <parallel_for.hpp>

#include <laddrcore.hpp>
#include <lasxgrad.hpp>
#include <lasxcurl.hpp>
#include <lasxdiv.hpp>

 /*!
  * @defgroup SDDR_yangmills
  * @brief Implementation of the arbitrary order Serendipity DDR scheme for the %Yang-Mills problem 
  */

namespace HArDCore3D
{

   /*!
    * @addtogroup SDDR_yangmills
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
           const LADDRCore & laddrcore,       ///< Has extra P2k3_T space for nonlinear operator
           const LieAlgebra & liealgebra,     ///< Lie algebra of the problem
           size_t nonlinear_discretisation,   /// Discretisation of nonlinear terms (1: discrete xdiv bracket, 2: (C_h,C_h)_xdiv + integral bracket)
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
                          const Eigen::VectorXd & E_i,    ///< Electric field at previous time
                          const Eigen::VectorXd & A_i,    ///< Potential at previous time
                          double dt,                        ///< Time step
                          double theta,                     ///< Theta scheme
                          double nonlinear_coeff            ///< Scaling of the nonlinear terms (0 for Maxwell)
                        );
    /// Assembles the system for Newton iterations
    void assembleSystemNewton(
                          const Eigen::VectorXd & E_i,    ///< Electric field at previous time
                          const Eigen::VectorXd & A_i,    ///< Potential at previous time
                          const Eigen::VectorXd & Elambda_k, ///< Solution position after previous Newton iteration
                          double dt,                        ///< Time step
                          double theta,                     ///< Theta scheme
                          double nonlinear_coeff            ///< Scaling of the nonlinear terms (0 for Maxwell)
                        );
    /// Stops once changes become small
    double stoppingCrit(
                      const Eigen::VectorXd & v, 
                      const Eigen::VectorXd & u  
                    );
    /// Sets the nonlinear residual vector b-F(x_k) in m_b_k
    void setNonlinearRes(
                        const Eigen::VectorXd & Elambda_k, ///< Solution in Newton iterations
                        const Eigen::VectorXd & E_i,       ///< Electric field at previous time
                        const Eigen::VectorXd & A_i,       ///< Potential at previous time
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

    /// Does one of two things depending on the bracket discretisation chosen (m_nl_par)
    /// 0: Calculates the local matrix of *[v,.]^div. (i,j), i represents laxdiv index, j represents xcurl index of entry
    /// 0: Calculates the matrix of *[v,.] (i,j), i represents LaP2k3 basis index, j represents xcurl index of entry
    std::vector<Eigen::MatrixXd> epsBkt_v(size_t iT,                               ///< Element index
                             boost::multi_array<double, 3> & crossij_Pot_T,  ///<  The output of _compute_crossij_Pot_T representing the DDR nonlinear operator in the element
                             const std::vector<Eigen::VectorXd> & vec_list                ///< Vector to plug in
                             ) const;             

    /// Wrapper to plug vector into L2 integral product: int (Pcurl 1,[Pcurl 2, Pgrad 3])
    std::vector<Eigen::MatrixXd> L2v_Bkt(size_t iT,                                    ///< Element index
                            boost::multi_array<double, 3> & intPciPcjPgk, ///< L2 integral product (trilinear form)
                            const std::vector<Eigen::VectorXd> & vec_list, ///< Vector to plug in
                            const size_t & entry                          ///< Position to put vector (1, 2, 3)
                            ) const;
    
    /// Does one of two things depending on the bracket discretisation chosen (m_nl_par)
    /// 0: Calculates the bracket inside the L2prod with v: (v,*[.,.]^div)_L2, adds the contribution from the faces, then the element
    /// 1: Calculates the bracket inside the L2prod with v: \int_T <v,*[.,.]>
    std::vector<Eigen::MatrixXd> L2v_epsBkt(size_t iT,                     ///< Element index
                               boost::multi_array<double, 3> & crossij_Pot_T,  ///< The output of _compute_crossij_Pot_T representing the DDR nonlinear operator in the element
                               const std::vector<Eigen::VectorXd> & v_L2prod   ///< Local vector in T (because we often only have the local part of the vector e.g. (*[A,A]^div)_T), representing (v,.)_L2 or \int_T <v,.> depending on m_nl_par
                               ) const;

    /// Returns the global problem dimension 
    inline size_t dimensionSpace() const
    {
      return m_lasxcurl.dimension() + m_lasxgrad.dimension();
    }

    /// Returns the number of statically condensed DOFs (Electric field)
    inline size_t nbSCDOFs_E() const
    {
      return m_nloc_sc_E.sum();
    }

    /// Returns the number of statically condensed DOFs (lambda)
    inline size_t nbSCDOFs_lambda() const
    {
      return m_nloc_sc_lambda.sum();
    }

    /// Returns the number of statically condensed DOFs (both velocity and pressure)
    inline size_t nbSCDOFs() const
    {
      return nbSCDOFs_E() + nbSCDOFs_lambda();
    }

    /// Returns the size of the system
    inline size_t sizeSystem() const
    {
      return dimensionSpace() - nbSCDOFs(); 
    }

    /// Returns the Lie algebra
    inline const LieAlgebra & lieAlg() const
    {
      return m_liealg;
    }

    /// Returns the space LAXGrad
    inline const LASXGrad & laSXGrad() const
    {
      return m_lasxgrad;
    }

    /// Returns the space LAXCurl
    inline const LASXCurl & laSXCurl() const
    {
      return m_lasxcurl;
    }

    /// Returns the space LAXCurl
    inline const LASXDiv & laSXDiv() const
    {
      return m_lasxdiv;
    }

    /// Returns the space XDiv
    inline const SXGrad & xSGrad() const
    {
      return m_sxgrad;
    }

    /// Returns the space XDiv
    inline const SXCurl & xSCurl() const
    {
      return m_sxcurl;
    }

    /// Returns the space XDiv
    inline const SXDiv & xSDiv() const
    {
      return m_sxdiv;
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
    inline const Eigen::VectorXd & systemVectorNewton() const {
      return m_b;
    }

    /// Returns the linear system right-hand side vector
    inline Eigen::VectorXd & systemVectorNewton() {
      return m_b;
    }

    /// Returns the right-hand side to the nonlinear problem
    inline const Eigen::VectorXd & systemVector() const {
      return m_b_i;
    }

    /// Returns the right-hand side to the nonlinear problem
    inline Eigen::VectorXd & systemVector() {
      return m_b_i;
    }

    /// Returns the right-hand side to the nonlinear problem
    inline const Eigen::VectorXd & nonlinearRes() const {
      return m_b_k;
    }

    /// Returns the right-hand side to the nonlinear problem
    inline Eigen::VectorXd & nonlinearRes() {
      return m_b_k;
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

    Eigen::VectorXd _compute_local_vec(size_t iT, ///< Element index
                          const Eigen::VectorXd & interp_f,  ///< Forcing term of first equation
                          const Eigen::VectorXd & interp_dE, ///< Used in forcing term of second equation
                          const Eigen::VectorXd & interp_A,  ///< Used in forcing term of second equation
                          const Eigen::VectorXd & E_i,    ///< Electric field at previous time
                          const Eigen::VectorXd & A_i,    ///< Potential at previous time
                          double dt,                        ///< Time step
                          double theta,                     ///< Theta scheme
                          double nonlinear_coeff            ///< Scaling of the nonlinear terms (0 for Maxwell)
                          );

    void _assemble_local_vec(size_t iT, 
                             const Eigen::VectorXd & v,    ///< Local vector
                             std::vector<std::list<Eigen::Triplet<double>>> & triplets_sys,  ///< List of triplets for the system
                             std::vector<Eigen::VectorXd> & vecs   //< List of vectors for the system
                             );

    /// There is a switch to calculate either the nonlinear residual vector (F(E_k)), or the derivative matrix DF which uses the same calculations (and a bit extra)
    SystemVectors<Eigen::MatrixXd> 
    _compute_local_newton(
                          size_t iT, ///< Element index
                          size_t option_1FEk_2DF,     /// 1 to calculate F(E_k) vector, 2 for newton derivative matrix DF
                          const Eigen::VectorXd & E_i,    ///< Electric field at previous time
                          const Eigen::VectorXd & A_i,    ///< Potential at previous time
                          const Eigen::VectorXd & Elambda_k, ///< Solution position after previous Newton iteration
                          double dt,                        ///< Time step
                          double theta,                     ///< Theta scheme
                          double nonlinear_coeff            ///< Scaling of the nonlinear terms (0 for Maxwell)
                          );

    /// Creates the permutation matrix and the global DOFs for the local static condensation
    LocalStaticCondensation _compute_static_condensation(const size_t & iT) const;

    void _assemble_local_newton(
                                size_t iT,   ///< Element index
                                const SystemVectors<Eigen::MatrixXd> & lsT,    ///< Local contributions
                                std::vector<std::list<Eigen::Triplet<double> > > & triplets_sys,  ///< List of triplets for the system
                                std::vector<Eigen::VectorXd> & vecs   //< List of vectors for the system
                                );

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

    /// Builds off of _compute_detnij_PkF. Calculates the projection onto Pk_F of (P_curl e_i x P_curl e_j).nF (e_ basis of xcurl_F) and stores the coefficients (on Polyk) of the function in the index k
    boost::multi_array<double, 3> _compute_detnij_Pot_PkF(size_t iF) const;

    /// Calculates the projection onto Pk_F of (phi_i x phi_j).nF (phi_ basis of Polyk2) and stores the coefficients (on Polyk) of the function in the index k
    boost::multi_array<double, 3> _compute_detnij_PkF(size_t iF) const;

    /// Does one of two things depending on the bracket discretisation chosen (m_nl_par)
    /// 0: Builds off of _compute_crossij_T. Calculates the projection onto Pk_T of (P_curl e_i x P_curl e_j) (e_ basis of xcurl_T) and stores the coefficients (first Golykmo, then GolyComplk) of the function in the index k
    /// 1: Calculates the triple integral (P_curl e_i x P_curl e_j).psi_k (e_ basis of xcurl_T, psi_k basis of Poly2k3)
    boost::multi_array<double, 3> _compute_crossij_Pot_T(size_t iT) const; 

    /// Calulates the projection onto Golykmo and GolyComplk of (phi_i x phi_j) (phi_ basis of Polyk3) and stores the coefficients (first Golykmo, then GolyComplk) of the function in the index k
    boost::multi_array<double, 3> _compute_crossij_T(size_t iT) const;

    /// Computes integrals of three basis function potentials (Xcurl . Xcurl)(Xgrad)
    boost::multi_array<double, 3> _integral_Pot_ijk(size_t iT) const;

    const DDRCore & m_ddrcore;
    const LADDRCore & m_laddrcore;
    bool m_use_threads;
    std::ostream & m_output;
    SerendipityProblem m_ser_pro;
    const XDiv m_xdiv;
    const SXGrad m_sxgrad;
    const SXCurl m_sxcurl;
    const SXDiv m_sxdiv;
    const LieAlgebra m_liealg;
    const LASXGrad m_lasxgrad;
    const LASXCurl m_lasxcurl;
    const LASXDiv m_lasxdiv;
    const Eigen::VectorXi m_nloc_sc_E; // Nb of statically condensed DOFs for electric field in each cell (cell unknowns)
    const Eigen::VectorXi m_nloc_sc_lambda; // Nb of statically condensed DOFs for lambda in each cell (cell unknowns)
    SystemMatrixType m_A;   // Matrix and RHS system
    Eigen::VectorXd m_b;    // RHS to Newton problem after static condensation
    Eigen::VectorXd m_b_i;  // RHS to nonlinear problem (F(x)=b)
    Eigen::VectorXd m_b_k;  // RHS/Residual to Newton problem  (DF(dx)=b-F(x_(k-1))
    SystemMatrixType m_sc_A; // Static condensation operator and RHS (to recover statically condensed DOFs)
    Eigen::VectorXd m_sc_b;
    SystemMatrixType m_laxgrad_L2;   // Global laxgrad L2Product matrix
    SystemMatrixType m_laxcurl_L2;   // Global laxcurl L2Product matrix
    SystemMatrixType m_laxdiv_L2;   // Global laxdiv L2Product matrix
    int m_nl_par;
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
 
} // end of namespace HArDCore3D

#endif
