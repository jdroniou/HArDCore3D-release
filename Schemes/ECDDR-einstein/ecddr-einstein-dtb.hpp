// Implementation of the discrete exterior calculus de Rham sequence for Einstein's equations in differential form formulation.
//
// Author: Jia Jia Qian(jiajia.qian1@monash.edu)
//
/*
 *
 *   If using this scheme in scientific publications, please cite the following article:
 *
 *    A polytopal discrete de Rham scheme for the exterior calculus Einstein's equations
 *    Todd A. Oliynyk and Jia Jia Qian
 *    arxiv: https://arxiv.org/abs/2505.00286
 */


#ifndef ECDDR_EINSTEIN_DTB_HPP
#define ECDDR_EINSTEIN_DTB_HPP

#include <iostream>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include <mesh.hpp>
#include <mesh_builder.hpp>
#include <local_static_condensation.hpp>
#include <linearsolver.hpp>
#include <parallel_for.hpp>

#include "ddr_spaces.hpp"
#include "ddr_pec.hpp"
#include "ddr_spaces.hpp"

namespace HArDCore3D
{
  /// Structure to store norm components (for B_\alpha, D_\alpha, theta^\alpha)
  struct EinsteinNorms
  {
    /// Constructor
    EinsteinNorms(Eigen::VectorXd norm_D,  Eigen::VectorXd norm_theta,  Eigen::VectorXd norm_B):
      D(norm_D),
      theta(norm_theta),
      B(norm_B)
      {
        // do nothing
      };

    // Constructor without arguments
    EinsteinNorms() : EinsteinNorms(Eigen::VectorXd::Zero(3), Eigen::VectorXd::Zero(3), Eigen::VectorXd::Zero(3)){};
    
    Eigen::VectorXd D;     ///< Norms of D_alpha
    Eigen::VectorXd theta; ///< Norms of theta^\alpha
    Eigen::VectorXd B;     ///< Norms of B_alpha
  };

  /// Class to assemble and solve an Einstein problem 
  struct Einstein
  {
    typedef Eigen::SparseMatrix<double> SystemMatrixType;
    
    typedef std::function<double(const Eigen::Vector3d &)> ZeroFormType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> OneFormType;
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> TwoFormType;
    typedef std::function<double(const Eigen::Vector3d &)> ThreeFormType;
    typedef TwoFormType ForcingTermType;

    typedef std::vector<ZeroFormType> Collection0Forms;
    typedef std::vector<OneFormType> Collection1Forms;
    typedef std::vector<TwoFormType> Collection2Forms;
    typedef std::vector<ThreeFormType> Collection3Forms;

    typedef std::function<double(const double t, const Eigen::Vector3d &)> TZeroFormType;
    typedef std::function<Eigen::Vector3d(const double t, const Eigen::Vector3d &)> TOneFormType;
    typedef std::function<Eigen::Vector3d(const double t, const Eigen::Vector3d &)> TTwoFormType;
    typedef std::function<double(const double t, const Eigen::Vector3d &)> TThreeFormType;
    typedef TTwoFormType TForcingTermType;

    typedef std::vector<ZeroFormType> TCollection0Forms;
    typedef std::vector<TOneFormType> TCollection1Forms;
    typedef std::vector<TTwoFormType> TCollection2Forms;
    typedef std::vector<TThreeFormType> TCollection3Forms;

    typedef std::function<Eigen::Vector3d(const Eigen::MatrixXd & E, const Eigen::MatrixXd & B)> CompositeHType;
    typedef std::vector<CompositeHType> CollectionHforms;
    typedef std::function<Eigen::Vector3d(const Eigen::MatrixXd & D)> CompositeEType;
    typedef std::vector<CompositeEType> CollectionEforms;
    typedef std::function<Eigen::Vector3d(const Eigen::VectorXd & E0, const Eigen::VectorXd & D0, const Eigen::MatrixXd & H, const Eigen::MatrixXd & E, const Eigen::MatrixXd & D, const Eigen::MatrixXd & B)> CompositeUType;
    typedef std::vector<CompositeUType> CollectionUforms; 

    /// Constructor
    Einstein(
           const DDR_Spaces & ddrspaces,           ///< Core for the ECDDR space sequence
           bool use_threads,                  ///< True for parallel execution, false for sequential execution
           double stab_par,
           std::ostream & output = std::cout  ///< Output stream to print status messages
         );

    //--------------------------------------------------
    //    Methods to assemble and solve
    //--------------------------------------------------    

    Eigen::VectorXd updateLapse(double dt, Eigen::VectorXd & lapse_n,  Eigen::VectorXd & starDtB);
    /// 
    void precomputeL2();

    void assembleLinearSystem(
                              double dt,
                              const ZeroFormType & lapse,
                              const OneFormType & E0,
                              const Eigen::VectorXd & starDtB
                              );

    void computeBCRHS(  
                      double dt,
                      const ZeroFormType & lapse_n,
                      const std::vector<TwoFormType> & D,
                      const std::vector<OneFormType> & theta,
                      const std::vector<TwoFormType> & B,
                      const std::vector<TwoFormType> & dtD,
                      const std::vector<OneFormType> & E,
                      const std::vector<OneFormType> & H,
                      const std::vector<TwoFormType> & dH
                      );


    void updateThetaRHS(  
                        const Eigen::VectorXd & starD
                        );

    /// Compute the residual
    double computeResidual(
                          const Eigen::VectorXd & starDB
                          ) const;

    /// Compute the constraints
    Eigen::VectorXd computeConstraints(
                              const ZeroFormType & lapse_n,
                              const std::vector<OneFormType> & theta0,
                              const std::vector<OneFormType> & thetan,
                              const Eigen::VectorXd & starDB0,
                              const Eigen::VectorXd & starDBn,
                              LinearSolver<Einstein::SystemMatrixType> & solver
                              ) const;

    /// Compute the discrete L2 norms, for a family of Eigen::VectorXd 
    std::vector<EinsteinNorms> computeEinsteinNorms(
                const std::vector<Eigen::VectorXd> & list_dofs ///< List of vectors representing velocities and pressures
                ) const;

    std::pair<EinsteinNorms, EinsteinNorms> computeContinuousErrorsNorms(
                const Eigen::VectorXd & starDt, ///< List of vectors representing velocities and pressures
                const std::vector<OneFormType> & starD, 
                const std::vector<TwoFormType> & starTheta,
                const std::vector<TwoFormType> & starB
                ) const;

    /// Change of basis for 2 forms
    Eigen::MatrixXd basisChange(
                                size_t k,                      ///< form degree 1 or 2
                                const Eigen::MatrixXd & forms, ///< Forms either 1-form [A_1, A_2, A_3] or 2-form [B_12, B_13, B_23] in column vectors
                                const Eigen::MatrixXd & basis  ///< Basis vectors by column
                                ) const; 

    /// Contract collection index with first form index (A_i)^{ji}
    Eigen::VectorXd contractIndex(
                                size_t index,                      ///< first or both
                                const Eigen::MatrixXd & forms ///< Forms either 1-form [A_1, A_2, A_3] or 2-form [B_12, B_13, B_23] in column vectors
                                ) const; 

    /// Contract collection index with first form index (A_i)^{ji}
    Eigen::MatrixXd contractIndex(
                                size_t index,                      ///< first or both
                                const Eigen::MatrixXd & twoForm1, ///< Forms either 1-form [A_1, A_2, A_3] or 2-form [B_12, B_13, B_23] in column vectors
                                const Eigen::MatrixXd & twoForm2  ///< Forms either 1-form [A_1, A_2, A_3] or 2-form [B_12, B_13, B_23] in column vectors
                                ) const; 

    /// wedge a 1-form and a 1-form to return 2-form as [A12, A13, A23]
    Eigen::Vector3d wedge11( 
                            const Eigen::Vector3d & form1,
                            const Eigen::Vector3d & form2
                            ) const; 

    /// expand a 2 form vector into a 3x3 matrix
    Eigen::Matrix3d expand2Form( 
                                const Eigen::Vector3d & form ///< 2 form as [A_12, A_13, A_23]
                                ) const; 

    //--------------------------------------------------------
    //  Nonlinear calculations
    //--------------------------------------------------------

    Einstein::CompositeHType
    composite_orth_H0 = [](const Eigen::VectorXd & orth_E0, const Eigen::Matrix3d & orth_B) -> Eigen::Vector3d {
                          return Eigen::Vector3d(
                                                0,
                                                0,
                                                0
                                                );
                        };

    Einstein::CompositeHType 
    composite_orth_H1 = [](const Eigen::VectorXd & orth_E0, const Eigen::Matrix3d & orth_B) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              -0.5*(-orth_B(1,1)+orth_B(0,2))+0.5*orth_B(2,0),
                                              orth_E0(2)+orth_B(2,1),
                                              -orth_E0(1)+orth_B(2,2)
                                              );
                  };

    Einstein::CompositeHType
    composite_orth_H2 = [](const Eigen::VectorXd & orth_E0, const Eigen::Matrix3d & orth_B) -> Eigen::Vector3d {
                          return Eigen::Vector3d(
                                                -orth_E0(2)-orth_B(1,0),
                                                -0.5*(orth_B(0,2)+orth_B(2,0))-0.5*orth_B(1,1),
                                                orth_E0(0)-orth_B(1,2)
                                                );
                        };

    Einstein::CompositeHType
    composite_orth_H3 = [](const Eigen::VectorXd & orth_E0, const Eigen::Matrix3d & orth_B) -> Eigen::Vector3d {
                          return Eigen::Vector3d(
                                                orth_E0(1)+orth_B(0,0),
                                                -orth_E0(0)+orth_B(0,1),
                                                -0.5*(orth_B(2,0)-orth_B(1,1))+0.5*orth_B(0,2)
                                                );
                        };

  Einstein::CollectionHforms composite_orth_H = {composite_orth_H0, composite_orth_H1, composite_orth_H2, composite_orth_H3};

    Einstein::CompositeEType
    composite_orth_E0 = [](const Eigen::Matrix3d & orth_D) -> Eigen::Vector3d {
                          return Eigen::Vector3d(
                                                0,
                                                0,
                                                0
                                                );
                        };

    Einstein::CompositeEType 
    composite_orth_E1 = [](const Eigen::Matrix3d & orth_D) -> Eigen::Vector3d {
                          Eigen::Matrix3d starD = Eigen::Matrix<double,3,3>{{0.,0.,1.},{0.,-1.,0.},{1.,0.,0.}} * orth_D;
                          return Eigen::Vector3d(
                                                -0.5*starD(0,0)+0.5*(starD(1,1)+starD(2,2)),
                                                -starD(0,1),
                                                -starD(0,2)
                                                );
                  };

    Einstein::CompositeEType
    composite_orth_E2 = [](const Eigen::Matrix3d & orth_D) -> Eigen::Vector3d {
                          Eigen::Matrix3d starD = Eigen::Matrix<double,3,3>{{0.,0.,1.},{0.,-1.,0.},{1.,0.,0.}} * orth_D;
                          return Eigen::Vector3d(
                                                -starD(1,0),
                                                -0.5*starD(1,1)+0.5*starD(0,0)+0.5*starD(2,2),
                                                -starD(1,2)
                                                );
                        };

    Einstein::CompositeEType
    composite_orth_E3 = [](const Eigen::Matrix3d & orth_D) -> Eigen::Vector3d {
                          Eigen::Matrix3d starD = Eigen::Matrix<double,3,3>{{0.,0.,1.},{0.,-1.,0.},{1.,0.,0.}} * orth_D;
                          return Eigen::Vector3d(
                                                -starD(2,0),
                                                -starD(2,1),
                                                -0.5*starD(2,2)+0.5*starD(0,0)+0.5*starD(1,1)
                                                );
                        };

  Einstein::CollectionEforms composite_orth_E = {composite_orth_E0, composite_orth_E1, composite_orth_E2, composite_orth_E3};

  Einstein::CompositeUType
  composite_orth_U0 = [](const Eigen::Vector3d & E0, const Eigen::Vector3d & D0, const Eigen::Matrix3d & H, const Eigen::Matrix3d & E, const Eigen::Matrix3d & D, const Eigen::Matrix3d & B) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              0, 
                                              0, 
                                              0
                                              );
                      };

  Einstein::CompositeUType
  composite_orth_U1 = [](const Eigen::Vector3d & E0, const Eigen::Vector3d & D0, const Eigen::Matrix3d & H, const Eigen::Matrix3d & E, const Eigen::Matrix3d & D, const Eigen::Matrix3d & B) -> Eigen::Vector3d {
                          Eigen::Matrix<double,3,3> hstar{{0.,0.,1.},{0.,-1.,0.},{1.,0.,0.}};
                          Eigen::Vector3d starD0 = hstar * D0;
                          Eigen::Matrix3d starD = hstar * D;
                          Eigen::Matrix3d starB = hstar * B;
                          Eigen::MatrixXd starU = Eigen::MatrixXd::Identity(3,3) 
                          * 0.5 * (E0.transpose()*starD0 + (E.transpose()*starD).trace() 
                                                                - (H.transpose()*starB).trace())
                        - starD0*E0.transpose() - starD * E.transpose() + starB * H.transpose();
                        return (hstar*starU).col(0);
                      };

  Einstein::CompositeUType
  composite_orth_U2 = [](const Eigen::Vector3d & E0, const Eigen::Vector3d & D0, const Eigen::Matrix3d & H, const Eigen::Matrix3d & E, const Eigen::Matrix3d & D, const Eigen::Matrix3d & B) -> Eigen::Vector3d {
                        Eigen::Matrix<double,3,3> hstar{{0.,0.,1.},{0.,-1.,0.},{1.,0.,0.}};
                        Eigen::Vector3d starD0 = hstar * D0;
                        Eigen::Matrix3d starD = hstar * D;
                        Eigen::Matrix3d starB = hstar * B;
                        Eigen::MatrixXd starU = Eigen::MatrixXd::Identity(3,3) 
                        * 0.5 * (E0.transpose()*starD0 + (E.transpose()*starD).trace() 
                                                              - (H.transpose()*starB).trace())
                      - starD0*E0.transpose() - starD * E.transpose() + starB * H.transpose();
                      return (hstar*starU).col(1);
                      };

  Einstein::CompositeUType
  composite_orth_U3 = [](const Eigen::Vector3d & E0, const Eigen::Vector3d & D0, const Eigen::Matrix3d & H, const Eigen::Matrix3d & E, const Eigen::Matrix3d & D, const Eigen::Matrix3d & B) -> Eigen::Vector3d {
                        Eigen::Matrix<double,3,3> hstar{{0.,0.,1.},{0.,-1.,0.},{1.,0.,0.}};
                        Eigen::Vector3d starD0 = hstar * D0;
                        Eigen::Matrix3d starD = hstar * D;
                        Eigen::Matrix3d starB = hstar * B;
                        Eigen::MatrixXd starU = Eigen::MatrixXd::Identity(3,3) 
                        * 0.5 * (E0.transpose()*starD0 + (E.transpose()*starD).trace() 
                                                              - (H.transpose()*starB).trace())
                      - starD0*E0.transpose() - starD * E.transpose() + starB * H.transpose();
                      return (hstar*starU).col(2);
                      };

  Einstein::CollectionUforms composite_orth_U = {composite_orth_U0, composite_orth_U1, composite_orth_U2, composite_orth_U3};

    //------------------------------------------------------------------
    //    Getters: dimension of space and numbers of unknowns of various types
    //------------------------------------------------------------------
    
    /// Returns the global problem dimension
    inline size_t dimensionSpace() const
    {
      return 3*(m_ddrspaces.dofspace(1).dimension() + m_ddrspaces.dofspace(2).dimension() + m_ddrspaces.dofspace(1).dimension());
    }

    //--------------------------------------------------------
    //  Getters: underlying ECDDR sequence
    //--------------------------------------------------------
    
    /// Returns the ECDDR sequence
    inline const DDR_Spaces & ddrSpaces() const
    {
      return m_ddrspaces;
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

    /// Returns the L of the LLT decomposition of the L2 product of X1
    inline const SystemMatrixType & L2_X1() const {
      return m_L2_X1;
    }

    /// Returns the L of the LLT decomposition of the L2 product of X1
    inline SystemMatrixType & L2_X1() {
      return m_L2_X1;
    }

    /// Returns the L of the LLT decomposition of the L2 product of X2
    inline const SystemMatrixType & L2_X2() const {
      return m_L2_X2;
    }

    /// Returns the L of the LLT decomposition of the L2 product of X2
    inline SystemMatrixType & L2_X2() {
      return m_L2_X2;
    }

    /// Returns the lhs of theta equation
    inline const SystemMatrixType & lhsTheta() const {
      return m_A_theta;
    }

    /// Returns the lhs of theta equation
    inline SystemMatrixType & lhsTheta() {
      return m_A_theta;
    }

    /// Returns the rhs of theta equation
    inline const Eigen::VectorXd & rhsTheta() const {
      return m_b_theta;
    }

    /// Returns the rhs of theta equation
    inline Eigen::VectorXd & rhsTheta() {
      return m_b_theta;
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

    inline double & min_detTheta() {
      return m_min_detTheta;
    }

    //-------------------------------------------
    //  Manipulate functions
    //-------------------------------------------

    /// Takes a two parameter function and a value and returns a function with the first parameter fixed to value
    template<typename outValue, typename TFct>
    std::function<outValue(const Eigen::Vector3d &)> contractPara(const TFct &F, double t) const
    {
      std::function<outValue(const Eigen::Vector3d &)> f = [&F, t](const Eigen::Vector3d &x) -> outValue {return F(t, x);};
      return f;
    }
    /// Takes a vector of two parameter functions and a value and returns a function with the first parameter fixed to value
    template<typename outValue, typename Fct, typename TFct>
    std::vector<Fct> contractParaVec(const std::vector<TFct> &TF, double t) const
    { 
      std::vector<Fct> f;
      for (size_t i = 0; i < TF.size(); i++){
        f.emplace_back(contractPara<outValue>(TF[i], t));
      }
      return f;
    }

  private:   

    /// Assemble the local contribution for the element of index iT into the global system
    void _assemble_local_contribution(
                                      size_t iT,
                                      const SystemVectors<Eigen::MatrixXd> & lsT,
                                      std::vector<std::list<Eigen::Triplet<double>>> & triplets_sys,
                                      std::vector<Eigen::VectorXd> & vec
                                      );   

    SystemVectors<Eigen::MatrixXd>
    _compute_local_contribution(  
                                const double dt,                      ///< Time increment
                                const size_t iT,                      ///< Element index
                                const ZeroFormType & lapse,
                                const OneFormType & E0,
                                const Eigen::VectorXd & starDtB
                               );

    std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd, Eigen::MatrixXd, double>
    _compute_nonlinear_relations(
                                const size_t iT,                      ///< Element index
                                const ZeroFormType & lapse,
                                const OneFormType & E0,
                                const Eigen::VectorXd & starDtB
                                ); 

    // Constructor parameter
    const DDR_Spaces & m_ddrspaces;
    bool m_use_threads;
    std::ostream & m_output;

    // Product matrix and RHS
    SystemMatrixType m_A;
    Eigen::VectorXd m_b;
    SystemMatrixType m_L2_X1;
    SystemMatrixType m_L2_X2;
    SystemMatrixType m_A_theta;
    Eigen::VectorXd m_b_theta;
    double m_stab_par;
    double m_min_detTheta;
  };

  static const double PI = boost::math::constants::pi<double>();
  using std::sin;
  using std::cos;
  using std::pow;
  using std::cyl_bessel_j;

  //------------------------------------------------------------------------------
  // Hodge star of functions
  //------------------------------------------------------------------------------

  // Hodge star of a continuous function with space and time variables
  std::function<Eigen::Vector3d(double t, const Eigen::Vector3d &)> star(const std::function<Eigen::Vector3d(double t, const Eigen::Vector3d &)>  &F)
  {
    std::function<Eigen::Vector3d(double t, const Eigen::Vector3d &)> f = [&F](double t, const Eigen::Vector3d &x) -> Eigen::Vector3d {return Eigen::Matrix<double,3,3>{{0.,0.,1.},{0.,-1.,0.},{1.,0.,0.}}*F(t, x);};
    return f;
  }
  // Hodge star of a continuous function with space variable
  std::function<Eigen::Vector3d(const Eigen::Vector3d &)> star(const std::function<Eigen::Vector3d(const Eigen::Vector3d &)>  &F)
  {
    std::function<Eigen::Vector3d(const Eigen::Vector3d &)> f = [&F](const Eigen::Vector3d &x) -> Eigen::Vector3d {return Eigen::Matrix<double,3,3>{{0.,0.,1.},{0.,-1.,0.},{1.,0.,0.}}*F(x);};
    return f;
  }
  // Hodge star of a collection of functions in a vector
  template<typename Fct>
  std::vector<Fct> starCol(const std::vector<Fct> &F)
  { 
    std::vector<Fct> f;
    for (size_t i = 0; i < F.size(); i++){
      f.emplace_back(star(F[i]));
    }
    return f;
  }
  // Random noise function
  std::function<Eigen::Vector3d(const Eigen::Vector3d &)> noise(const std::function<Eigen::Vector3d(const Eigen::Vector3d &)> &f)
  { 
    std::function<Eigen::Vector3d(const Eigen::Vector3d &)> fnoise = [f](const Eigen::Vector3d &x) -> Eigen::Vector3d {return f(x)+Eigen::Vector3d::Random()*1e-10;};

    return fnoise;
  }


  //------------------------------------------------------------------------------
  // Minkowski solution (D, theta, B)
  //------------------------------------------------------------------------------

  static Einstein::TZeroFormType
  minkowski_lapse = [](const double t, const Eigen::Vector3d & x) -> double {
                        return 1.0;
                      };

  static Einstein::TOneFormType
  minkowski_theta0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              0, 
                                              0, 
                                              0
                                              );
                      };

  static Einstein::TOneFormType
  minkowski_theta1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              1., 
                                              0, 
                                              0
                                              );
                      };

  static Einstein::TOneFormType
  minkowski_theta2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              0, 
                                              1., 
                                              0
                                              );
                      };

  static Einstein::TOneFormType
  minkowski_theta3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              0, 
                                              0, 
                                              1.
                                              );
                      };

  static Einstein::TCollection1Forms minkowski_theta = {minkowski_theta0, minkowski_theta1, minkowski_theta2, minkowski_theta3};

  static Einstein::TOneFormType
  minkowski_E0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TOneFormType
  minkowski_E1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TOneFormType
  minkowski_E2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TOneFormType
  minkowski_E3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TCollection1Forms minkowski_E = {minkowski_E0, minkowski_E1, minkowski_E2, minkowski_E3};

  static Einstein::TTwoFormType
  minkowski_dE0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_dE1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_dE2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_dE3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TCollection1Forms minkowski_dE = {minkowski_dE0, minkowski_dE1, minkowski_dE2, minkowski_dE3};

  static Einstein::TTwoFormType
  minkowski_B0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_B1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_B2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_B3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TCollection2Forms minkowski_B = {minkowski_B0, minkowski_B1, minkowski_B2, minkowski_B3};

  static Einstein::TTwoFormType
  minkowski_dtB0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_dtB1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_dtB2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_dtB3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TCollection2Forms minkowski_dtB = {minkowski_dtB0, minkowski_dtB1, minkowski_dtB2, minkowski_dtB3};

  static Einstein::TTwoFormType
  minkowski_D0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_D1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_D2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_D3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TCollection2Forms minkowski_D = {minkowski_D0, minkowski_D1, minkowski_D2, minkowski_D3};

  static Einstein::TTwoFormType
  minkowski_dtD0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_dtD1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_dtD2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_dtD3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TCollection2Forms minkowski_dtD = {minkowski_dtD0, minkowski_dtD1, minkowski_dtD2, minkowski_dtD3};

  static Einstein::TOneFormType
  minkowski_H0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TOneFormType
  minkowski_H1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TOneFormType
  minkowski_H2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TOneFormType
  minkowski_H3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TCollection1Forms minkowski_H = {minkowski_H0, minkowski_H1, minkowski_H2, minkowski_H3};

  static Einstein::TTwoFormType
  minkowski_dH0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_dH1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_dH2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  minkowski_dH3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                    );
                            }; 

  static Einstein::TCollection1Forms minkowski_dH = {minkowski_dH0, minkowski_dH1, minkowski_dH2, minkowski_dH3};

  //------------------------------------------------------------------------------
  // Kasner solution (D, theta, B) (1/2, 1/4(1-sqrt(5)), 1/4(1+sqrt(5)))
  //------------------------------------------------------------------------------

  static Einstein::TZeroFormType
  kasner_lapse = [](const double t, const Eigen::Vector3d & x) -> double {
                        return pow(t, 1.0);
                      };

  static Einstein::TOneFormType
  kasner_theta0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, 0, 0);
                      };

  static Einstein::TOneFormType
  kasner_theta1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(sqrt(t), 0, 0);
                      };

  static Einstein::TOneFormType
  kasner_theta2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, pow(t, -0.30901699437494745), 0);
                      };

  static Einstein::TOneFormType
  kasner_theta3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, 0, pow(t, 0.80901699437494745));
                      };

  static Einstein::TCollection1Forms kasner_theta = {kasner_theta0, kasner_theta1, kasner_theta2, kasner_theta3};

  //  This is E0 = 1/N dN
  static Einstein::TOneFormType
  kasner_E0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  kasner_E1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0.5*pow(t, -0.5), 0, 0);
                            }; 

  static Einstein::TOneFormType
  kasner_E2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, -0.30901699437494745*pow(t, -1.3090169943749475), 0);
                            }; 

  static Einstein::TOneFormType
  kasner_E3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0.80901699437494745*pow(t, -0.19098300562505255));
                            }; 

  static Einstein::TCollection1Forms kasner_E = {kasner_E0, kasner_E1, kasner_E2, kasner_E3};

//------------------------------------------------------------------------------
// Det h E
  static Einstein::TOneFormType
  deth_kasner_E0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  deth_kasner_E1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0.5*pow(t, -1.5), 0, 0);
                            }; 

  static Einstein::TOneFormType
  deth_kasner_E2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, -0.30901699437494745*pow(t, -2.3090169943749475), 0);
                            }; 

  static Einstein::TOneFormType
  deth_kasner_E3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0.80901699437494745*pow(t, -1.1909830056250525));
                            }; 

  static Einstein::TCollection1Forms deth_kasner_E = {deth_kasner_E0, deth_kasner_E1, deth_kasner_E2, deth_kasner_E3};

//------------------------------------------------------------------------------

  static Einstein::TTwoFormType
  kasner_D0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  kasner_D1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0.5*pow(t, -0.5));
                            }; 

  static Einstein::TTwoFormType
  kasner_D2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, -1.3090169943749475*pow(t, 0.30901699437494745), 0);
                            }; 

  static Einstein::TTwoFormType
  kasner_D3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0.19098300562505255*pow(t, -0.80901699437494745), 0, 0);
                            }; 

  static Einstein::TCollection2Forms kasner_D = {kasner_D0, kasner_D1, kasner_D2, kasner_D3};

//------------------------------------------------------------------------------
// Det h D

  static Einstein::TTwoFormType
  deth_kasner_D0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  deth_kasner_D1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0.5*pow(t, -1.5));
                            }; 

  static Einstein::TTwoFormType
  deth_kasner_D2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, -1.3090169943749475*pow(t, -0.69098300562505255), 0);
                            }; 

  static Einstein::TTwoFormType
  deth_kasner_D3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0.19098300562505255*pow(t, -1.8090169943749475), 0, 0);
                            }; 

  static Einstein::TCollection2Forms deth_kasner_D = {deth_kasner_D0, deth_kasner_D1, deth_kasner_D2, deth_kasner_D3};

//------------------------------------------------------------------------------

  static Einstein::TTwoFormType
  kasner_dtD0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  kasner_dtD1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, -0.25*pow(t, -1.5));
                            }; 

  static Einstein::TTwoFormType
  kasner_dtD2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, -0.40450849718747378*pow(t, -0.69098300562505255), 0);
                            }; 

  static Einstein::TTwoFormType
  kasner_dtD3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(-0.1545084971874737*pow(t, -1.8090169943749475), 0, 0);
                            }; 

  static Einstein::TCollection2Forms kasner_dtD = {kasner_dtD0, kasner_dtD1, kasner_dtD2, kasner_dtD3};

//------------------------------------------------------------------------------
// Det h dtD

  static Einstein::TTwoFormType
  deth_kasner_dtD0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 
 
  static Einstein::TTwoFormType
  deth_kasner_dtD1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, -0.75*pow(t, -2.5));
                            }; 

  static Einstein::TTwoFormType
  deth_kasner_dtD2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0.90450849718747373*pow(t, -1.6909830056250525), 0);
                            }; 

  static Einstein::TTwoFormType
  deth_kasner_dtD3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(-0.34549150281252622*pow(t, -2.8090169943749475), 0, 0);
                            }; 

  static Einstein::TCollection2Forms deth_kasner_dtD = {deth_kasner_dtD0, deth_kasner_dtD1, deth_kasner_dtD2, deth_kasner_dtD3};

//------------------------------------------------------------------------------

  static Einstein::TTwoFormType
  kasner_B0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  kasner_B1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  kasner_B2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  kasner_B3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection2Forms kasner_B = {kasner_B0, kasner_B1, kasner_B2, kasner_B3};

//------------------------------------------------------------------------------

  static Einstein::TOneFormType
  kasner_H0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  kasner_H1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  kasner_H2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  kasner_H3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection1Forms kasner_H = {kasner_H0, kasner_H1, kasner_H2, kasner_H3};

//------------------------------------------------------------------------------
// Det h H

  static Einstein::TOneFormType
  deth_kasner_H0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  deth_kasner_H1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  deth_kasner_H2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  deth_kasner_H3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection1Forms deth_kasner_H = {deth_kasner_H0, deth_kasner_H1, deth_kasner_H2, deth_kasner_H3};

//------------------------------------------------------------------------------

  static Einstein::TTwoFormType
  kasner_dH0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  kasner_dH1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  kasner_dH2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  kasner_dH3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection2Forms kasner_dH = {kasner_dH0, kasner_dH1, kasner_dH2, kasner_dH3};

//------------------------------------------------------------------------------
// Det h dH

  static Einstein::TTwoFormType
  deth_kasner_dH0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  deth_kasner_dH1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  deth_kasner_dH2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  deth_kasner_dH3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection2Forms deth_kasner_dH = {deth_kasner_dH0, deth_kasner_dH1, deth_kasner_dH2, deth_kasner_dH3};

  //------------------------------------------------------------------------------
  // Gowdy solution (D, theta, B)
  //------------------------------------------------------------------------------

  static Einstein::TZeroFormType
  gowdy_lapse = [](const double t, const Eigen::Vector3d & x) -> double {
                        return pow(t, -0.25)*exp(0.5*pow(M_PI, 2)*pow(t, 2)*(pow(cyl_bessel_j(0, 2*M_PI*t), 2) + pow(cyl_bessel_j(1, 2*M_PI*t), 2)) - 0.5*M_PI*t*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t) - 0.5*pow(M_PI, 2)*(pow(cyl_bessel_j(1, 2*M_PI), 2) + pow(cyl_bessel_j(0, 2*M_PI), 2)) + 0.25*M_PI*cyl_bessel_j(0, 2*M_PI)*cyl_bessel_j(1, 2*M_PI));
                      };

  static Einstein::TOneFormType
  gowdy_theta0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              0, 
                                              0, 
                                              0
                                              );
                      };

  static Einstein::TOneFormType
  gowdy_theta1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              sqrt(t)*exp(0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t)), 
                                              0, 
                                              0
                                              );
                      };

  static Einstein::TOneFormType
  gowdy_theta2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              0, 
                                              sqrt(t)*exp(-0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t)), 
                                              0
                                              );
                      };

  static Einstein::TOneFormType
  gowdy_theta3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              0, 
                                              0, 
                                              pow(t, -0.25)*exp(M_PI*( 0.5*M_PI*pow(t, 2)*(pow(cyl_bessel_j(0, 2*M_PI*t), 2) + pow(cyl_bessel_j(1, 2*M_PI*t), 2)) 
                                                                      -0.5*t*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t) 
                                                                      -0.5*M_PI*(pow(cyl_bessel_j(1, 2*M_PI), 2) + pow(cyl_bessel_j(0, 2*M_PI), 2)) 
                                                                      + 0.25*cyl_bessel_j(0, 2*M_PI)*cyl_bessel_j(1, 2*M_PI)))
                                              );
                      };

  static Einstein::TCollection1Forms gowdy_theta = {gowdy_theta0, gowdy_theta1, gowdy_theta2, gowdy_theta3};

  static Einstein::TOneFormType
  gowdy_E0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                                0, 
                                                0, 
                                                -1.0*pow(M_PI, 2)*pow(t, 1.0)*sin(4*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t)
                                                );
                      };

  static Einstein::TOneFormType
  gowdy_E1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              pow(t, -0.25)*(-1.0*M_PI*pow(t, 1.0)*cos(2*M_PI*x(2))*cyl_bessel_j(1, 2*M_PI*t) + 0.5)*exp(M_PI*(-0.5*M_PI*pow(t, 2)*(pow(cyl_bessel_j(0, 2*M_PI*t), 2) + pow(cyl_bessel_j(1, 2*M_PI*t), 2)) + 0.5*t*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t) - 0.25*cyl_bessel_j(0, 2*M_PI)*cyl_bessel_j(1, 2*M_PI) + 0.5*M_PI*(pow(cyl_bessel_j(1, 2*M_PI), 2) + pow(cyl_bessel_j(0, 2*M_PI), 2))) + 0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t)), 
                                              0, 
                                              0
                                              );
                      };

  static Einstein::TOneFormType
  gowdy_E2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                               0, 
                                               pow(t, -0.25)*(1.0*M_PI*pow(t, 1.0)*cos(2*M_PI*x(2))*cyl_bessel_j(1, 2*M_PI*t) + 0.5)*exp(-0.5*pow(M_PI, 2)*pow(t, 2)*(pow(cyl_bessel_j(0, 2*M_PI*t), 2) + pow(cyl_bessel_j(1, 2*M_PI*t), 2)) + 0.5*M_PI*t*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t) - 0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t) - 0.25*M_PI*cyl_bessel_j(0, 2*M_PI)*cyl_bessel_j(1, 2*M_PI) + 0.5*pow(M_PI, 2)*(pow(cyl_bessel_j(1, 2*M_PI), 2) + pow(cyl_bessel_j(0, 2*M_PI), 2))), 
                                               0
                                              );
                      }; 

  static Einstein::TOneFormType
  gowdy_E3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                               0, 
                                               0, 
                                               1.0/t*(-1.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cos(2*M_PI*x(2)), 2)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) + 1.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cos(2*M_PI*x(2)), 2)*pow(cyl_bessel_j(1, 2*M_PI*t), 2) + 1.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) - 0.25)
                                               );
                      };

  static Einstein::TCollection1Forms gowdy_E = {gowdy_E0, gowdy_E1, gowdy_E2, gowdy_E3};

  static Einstein::TTwoFormType
  gowdy_dE0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              0, 
                                              0, 
                                              0
                                              );
                      };

  static Einstein::TTwoFormType
  gowdy_dE1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                               0, 
                                               0, 
                                               1.0/t*(-1.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cos(2*M_PI*x(2)), 2)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) + 1.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cos(2*M_PI*x(2)), 2)*pow(cyl_bessel_j(1, 2*M_PI*t), 2) + 1.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) - 0.25)
                                              );
                      };

  static Einstein::TTwoFormType
  gowdy_dE2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                               0, 
                                               0, 
                                               1.0*M_PI*pow(t, -0.25)*(2.0*M_PI*pow(t, 1.0)*cyl_bessel_j(1, 2*M_PI*t) + 2.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*pow(cyl_bessel_j(1, 2*M_PI*t), 2) - 0.5*cyl_bessel_j(0, 2*M_PI*t))*exp(-0.5*pow(M_PI, 2)*pow(t, 2)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) - 0.5*pow(M_PI, 2)*pow(t, 2)*pow(cyl_bessel_j(1, 2*M_PI*t), 2) + 0.5*M_PI*t*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t) - 0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t) - 0.25*M_PI*cyl_bessel_j(0, 2*M_PI)*cyl_bessel_j(1, 2*M_PI) + 0.5*pow(M_PI, 2)*pow(cyl_bessel_j(1, 2*M_PI), 2) + 0.5*pow(M_PI, 2)*pow(cyl_bessel_j(0, 2*M_PI), 2))*sin(2*M_PI*x(2))
                                              );
                      }; 

  static Einstein::TTwoFormType
  gowdy_dE3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                               0, 
                                               0, 
                                               0
                                               );
                      };

  static Einstein::TCollection1Forms gowdy_dE = {gowdy_dE0, gowdy_dE1, gowdy_dE2, gowdy_dE3};

  static Einstein::TOneFormType
  gowdy_H0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              0, 
                                              0, 
                                              0
                                              );
                      };

  static Einstein::TOneFormType
  gowdy_H1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              0, 
                                              M_PI*sqrt(t)*(-1.0*pow(t, 0.25) - 2.0*M_PI*pow(t, 1.25)*cos(2*M_PI*x(2))*cyl_bessel_j(1, 2*M_PI*t))*exp(-0.5*pow(M_PI, 2)*pow(t, 2)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) - 0.5*pow(M_PI, 2)*pow(t, 2)*pow(cyl_bessel_j(1, 2*M_PI*t), 2) + 0.5*M_PI*t*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t) - 0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t) - 0.25*M_PI*cyl_bessel_j(0, 2*M_PI)*cyl_bessel_j(1, 2*M_PI) + 0.5*pow(M_PI, 2)*pow(cyl_bessel_j(1, 2*M_PI), 2) + 0.5*pow(M_PI, 2)*pow(cyl_bessel_j(0, 2*M_PI), 2))*sin(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t), 
                                              0
                                              );
                      };

  static Einstein::TOneFormType
  gowdy_H2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              M_PI*sqrt(t)*(-1.0*pow(t, 0.25) + 2.0*M_PI*pow(t, 1.25)*cos(2*M_PI*x(2))*cyl_bessel_j(1, 2*M_PI*t))*exp(-0.5*pow(M_PI, 2)*pow(t, 2)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) - 0.5*pow(M_PI, 2)*pow(t, 2)*pow(cyl_bessel_j(1, 2*M_PI*t), 2) + 0.5*M_PI*t*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t) + 0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t) - 0.25*M_PI*cyl_bessel_j(0, 2*M_PI)*cyl_bessel_j(1, 2*M_PI) + 0.5*pow(M_PI, 2)*pow(cyl_bessel_j(1, 2*M_PI), 2) + 0.5*pow(M_PI, 2)*pow(cyl_bessel_j(0, 2*M_PI), 2))*sin(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t), 
                                              0, 
                                              0
                                              );
                      }; 

  static Einstein::TOneFormType
  gowdy_H3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              0, 
                                              0, 
                                              0
                                               );
                      };

  static Einstein::TCollection1Forms gowdy_H = {gowdy_H0, gowdy_H1, gowdy_H2, gowdy_H3};

  static Einstein::TTwoFormType
  gowdy_dH0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              0, 
                                              0, 
                                              0
                                              );
                      };

  static Einstein::TTwoFormType
  gowdy_dH1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              0, 
                                              0, 
                                              1.0*pow(M_PI, 2)*(1.0*pow(t, 0.75)*pow(sin(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t) + 2.0*pow(t, 0.75)*cos(2*M_PI*x(2)) - 4.0*M_PI*pow(t, 1.75)*pow(sin(2*M_PI*x(2)), 2)*cyl_bessel_j(1, 2*M_PI*t) + 4.0*M_PI*pow(t, 1.75)*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(1, 2*M_PI*t) - 4.0*pow(M_PI, 2)*pow(t, 2.75)*pow(sin(2*M_PI*x(2)), 2)*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*pow(cyl_bessel_j(1, 2*M_PI*t), 2))*exp(-0.5*pow(M_PI, 2)*pow(t, 2)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) - 0.5*pow(M_PI, 2)*pow(t, 2)*pow(cyl_bessel_j(1, 2*M_PI*t), 2) + 0.5*M_PI*t*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t) - 0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t) - 0.25*M_PI*cyl_bessel_j(0, 2*M_PI)*cyl_bessel_j(1, 2*M_PI) + 0.5*pow(M_PI, 2)*pow(cyl_bessel_j(1, 2*M_PI), 2) + 0.5*pow(M_PI, 2)*pow(cyl_bessel_j(0, 2*M_PI), 2))*cyl_bessel_j(0, 2*M_PI*t)
                                              );
                      };

  static Einstein::TTwoFormType
  gowdy_dH2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                               0, 
                                               1.0*pow(M_PI, 2)*(-1.0*pow(t, 0.75)*pow(sin(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t) + 2.0*pow(t, 0.75)*cos(2*M_PI*x(2)) + 4.0*M_PI*pow(t, 1.75)*pow(sin(2*M_PI*x(2)), 2)*cyl_bessel_j(1, 2*M_PI*t) - 4.0*M_PI*pow(t, 1.75)*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(1, 2*M_PI*t) + 4.0*pow(M_PI, 2)*pow(t, 2.75)*pow(sin(2*M_PI*x(2)), 2)*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*pow(cyl_bessel_j(1, 2*M_PI*t), 2))*exp(-0.5*pow(M_PI, 2)*pow(t, 2)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) - 0.5*pow(M_PI, 2)*pow(t, 2)*pow(cyl_bessel_j(1, 2*M_PI*t), 2) + 0.5*M_PI*t*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t) + 0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t) - 0.25*M_PI*cyl_bessel_j(0, 2*M_PI)*cyl_bessel_j(1, 2*M_PI) + 0.5*pow(M_PI, 2)*pow(cyl_bessel_j(1, 2*M_PI), 2) + 0.5*pow(M_PI, 2)*pow(cyl_bessel_j(0, 2*M_PI), 2))*cyl_bessel_j(0, 2*M_PI*t), 
                                               0
                                               );
                      }; 

  static Einstein::TTwoFormType
 gowdy_dH3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(
                                              0, 
                                              0, 
                                              0
                                              );
                      };

  static Einstein::TCollection1Forms gowdy_dH = {gowdy_dH0, gowdy_dH1, gowdy_dH2, gowdy_dH3};

  static Einstein::TTwoFormType
  gowdy_D0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                             0,
                                             0,
                                             0
                                              );
                         };

  static Einstein::TTwoFormType
  gowdy_D1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                            0, 
                                            0, 
                                            pow(t, -0.5)*(1.0*M_PI*pow(t, 1.0)*cos(2*M_PI*x(2))*cyl_bessel_j(1, 2*M_PI*t) - 1.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cos(2*M_PI*x(2)), 2)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) + 1.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cos(2*M_PI*x(2)), 2)*pow(cyl_bessel_j(1, 2*M_PI*t), 2) + 1.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) + 0.25)*exp(-0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t))
                                           );
                         };

  static Einstein::TTwoFormType
  gowdy_D2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                            0, 
                                            pow(t, -0.5)*(1.0*M_PI*pow(t, 1.0)*cos(2*M_PI*x(2))*cyl_bessel_j(1, 2*M_PI*t) + 1.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cos(2*M_PI*x(2)), 2)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) - 1.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cos(2*M_PI*x(2)), 2)*pow(cyl_bessel_j(1, 2*M_PI*t), 2) - 1.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) - 0.25)*exp(0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t)), 
                                            0
                                            );
                         };

  static Einstein::TTwoFormType
  gowdy_D3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                            1.0*pow(t, 0.25)*exp(M_PI*(-0.5*M_PI*pow(t, 2)*(pow(cyl_bessel_j(0, 2*M_PI*t), 2) + pow(cyl_bessel_j(1, 2*M_PI*t), 2)) + 0.5*t*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t) - 0.25*cyl_bessel_j(0, 2*M_PI)*cyl_bessel_j(1, 2*M_PI) + 0.5*M_PI*(pow(cyl_bessel_j(1, 2*M_PI), 2) + pow(cyl_bessel_j(0, 2*M_PI), 2)))), 
                                            0, 
                                            0
                                            );
                         };

  static Einstein::TCollection2Forms gowdy_D = {gowdy_D0, gowdy_D1, gowdy_D2, gowdy_D3};

  static Einstein::TTwoFormType
  gowdy_dtD0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                             0,
                                             0,
                                             0
                                              );
                         };

  static Einstein::TTwoFormType
  gowdy_dtD1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                            0, 
                                            0, 
                                            pow(t, -2.0)*(-0.125*sqrt(t) - 0.25*M_PI*pow(t, 1.5)*cos(2*M_PI*x(2))*cyl_bessel_j(1, 2*M_PI*t) - 1.5*pow(M_PI, 2)*pow(t, 2.5)*pow(cos(2*M_PI*x(2)), 2)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) + 0.5*pow(M_PI, 2)*pow(t, 2.5)*pow(cos(2*M_PI*x(2)), 2)*pow(cyl_bessel_j(1, 2*M_PI*t), 2) + 2.0*pow(M_PI, 2)*pow(t, 2.5)*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t) + 1.5*pow(M_PI, 2)*pow(t, 2.5)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) - 1.0*pow(M_PI, 3)*pow(t, 3.5)*pow(cos(2*M_PI*x(2)), 3)*pow(cyl_bessel_j(0, 2*M_PI*t), 2)*cyl_bessel_j(1, 2*M_PI*t) + 1.0*pow(M_PI, 3)*pow(t, 3.5)*pow(cos(2*M_PI*x(2)), 3)*pow(cyl_bessel_j(1, 2*M_PI*t), 3) + 8.0*pow(M_PI, 3)*pow(t, 3.5)*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t) + 1.0*pow(M_PI, 3)*pow(t, 3.5)*cos(2*M_PI*x(2))*pow(cyl_bessel_j(0, 2*M_PI*t), 2)*cyl_bessel_j(1, 2*M_PI*t) - 4.0*pow(M_PI, 3)*pow(t, 3.5)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t))*exp(-0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t))
                                            );
                         };

  static Einstein::TTwoFormType
  gowdy_dtD2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                            0, 
                                            pow(t, -2.0)*(0.125*sqrt(t) - 0.25*M_PI*pow(t, 1.5)*cos(2*M_PI*x(2))*cyl_bessel_j(1, 2*M_PI*t) + 1.5*pow(M_PI, 2)*pow(t, 2.5)*pow(cos(2*M_PI*x(2)), 2)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) - 0.5*pow(M_PI, 2)*pow(t, 2.5)*pow(cos(2*M_PI*x(2)), 2)*pow(cyl_bessel_j(1, 2*M_PI*t), 2) + 2.0*pow(M_PI, 2)*pow(t, 2.5)*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t) - 1.5*pow(M_PI, 2)*pow(t, 2.5)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) - 1.0*pow(M_PI, 3)*pow(t, 3.5)*pow(cos(2*M_PI*x(2)), 3)*pow(cyl_bessel_j(0, 2*M_PI*t), 2)*cyl_bessel_j(1, 2*M_PI*t) + 1.0*pow(M_PI, 3)*pow(t, 3.5)*pow(cos(2*M_PI*x(2)), 3)*pow(cyl_bessel_j(1, 2*M_PI*t), 3) - 8.0*pow(M_PI, 3)*pow(t, 3.5)*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t) + 1.0*pow(M_PI, 3)*pow(t, 3.5)*cos(2*M_PI*x(2))*pow(cyl_bessel_j(0, 2*M_PI*t), 2)*cyl_bessel_j(1, 2*M_PI*t) + 4.0*pow(M_PI, 3)*pow(t, 3.5)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t))*exp(0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t)), 
                                            0
                                            );
                         };

  static Einstein::TTwoFormType
  gowdy_dtD3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                            pow(t, -0.75)*(1.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cos(2*M_PI*x(2)), 2)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) - 1.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cos(2*M_PI*x(2)), 2)*pow(cyl_bessel_j(1, 2*M_PI*t), 2) - 1.0*pow(M_PI, 2)*pow(t, 2.0)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) + 0.25)*exp(M_PI*(-0.5*M_PI*pow(t, 2)*pow(cyl_bessel_j(0, 2*M_PI*t), 2) - 0.5*M_PI*pow(t, 2)*pow(cyl_bessel_j(1, 2*M_PI*t), 2) + 0.5*t*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t) - 0.25*cyl_bessel_j(0, 2*M_PI)*cyl_bessel_j(1, 2*M_PI) + 0.5*M_PI*pow(cyl_bessel_j(1, 2*M_PI), 2) + 0.5*M_PI*pow(cyl_bessel_j(0, 2*M_PI), 2))), 0, 
                                            0
                                            );
                         };

  static Einstein::TCollection2Forms gowdy_dtD = {gowdy_dtD0, gowdy_dtD1, gowdy_dtD2, gowdy_dtD3};

  static Einstein::TTwoFormType
  gowdy_B0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                            0, 
                                            0, 
                                            0
                                            );
                         };

  static Einstein::TTwoFormType
  gowdy_B1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    1.0*M_PI*sqrt(t)*exp(0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t))*sin(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t), 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  gowdy_B2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
	                                              		0, 
                                                    0, 
                                                    -1.0*M_PI*sqrt(t)*exp(-0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t))*sin(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t)
	                                              		);
                            }; 

  static Einstein::TTwoFormType
  gowdy_B3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    0, 
                                                    0
                                                     );
                            }; 

  static Einstein::TCollection2Forms gowdy_B = {gowdy_B0, gowdy_B1, gowdy_B2, gowdy_B3};

  static Einstein::TTwoFormType
  gowdy_dtB0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                      return Eigen::Vector3d(
                                            0, 
                                            0, 
                                            0
                                            );
                         };

  static Einstein::TTwoFormType
  gowdy_dtB1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                    0, 
                                                    M_PI*pow(t, -0.5)*(-M_PI*pow(t, 1.0)*(1.0*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t) + 2.0)*cyl_bessel_j(1, 2*M_PI*t) + 0.5*cyl_bessel_j(0, 2*M_PI*t))*exp(0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t))*sin(2*M_PI*x(2)), 
                                                    0
                                                    );
                            }; 

  static Einstein::TTwoFormType
  gowdy_dtB2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
	                                              		0, 
                                                    0, 
                                                    M_PI*pow(t, -0.5)*(M_PI*pow(t, 1.0)*(-1.0*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t) + 2.0)*cyl_bessel_j(1, 2*M_PI*t) - 0.5*cyl_bessel_j(0, 2*M_PI*t))*exp(-0.5*cos(2*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t))*sin(2*M_PI*x(2))
	                                              		);
                            }; 

  static Einstein::TTwoFormType
  gowdy_dtB3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(
                                                     0, 
                                                     0, 
                                                     0
                                                     );
                            }; 

  static Einstein::TCollection2Forms gowdy_dtB = {gowdy_dtB0, gowdy_dtB1, gowdy_dtB2, gowdy_dtB3};

  //------------------------------------------------------------------------------
  // Linear solution (D, theta, B)
  //------------------------------------------------------------------------------

  static Einstein::TZeroFormType
  linear_lapse = [](const double t, const Eigen::Vector3d & x) -> double {
                        return pow(t, 3)*pow(x(0), 2)*x(1) + pow(t, 3)*pow(x(0), 2)*x(2) + pow(t, 3)*x(0)*pow(x(1), 2) + 2*pow(t, 3)*x(0)*x(1)*x(2) + pow(t, 3)*x(0)*pow(x(2), 2) + pow(t, 3)*pow(x(1), 2)*x(2) + pow(t, 3)*x(1)*pow(x(2), 2) + pow(t, 2)*pow(x(0), 2) + 3*pow(t, 2)*x(0)*x(1) + 3*pow(t, 2)*x(0)*x(2) + pow(t, 2)*pow(x(1), 2) + 3*pow(t, 2)*x(1)*x(2) + pow(t, 2)*pow(x(2), 2) + 2*t*x(0) + 2*t*x(1) + 2*t*x(2) + 1;
                      };

  static Einstein::TOneFormType
  linear_theta0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, 0, 0);
                      };

  static Einstein::TOneFormType
  linear_theta1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(t*(x(0) + x(1)) + 1, 0, 0);
                      };

  static Einstein::TOneFormType
  linear_theta2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, t*(x(1) + x(2)) + 1, 0);
                      };

  static Einstein::TOneFormType
  linear_theta3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, 0, t*(x(0) + x(2)) + 1);
                      };

  static Einstein::TCollection1Forms linear_theta = {linear_theta0, linear_theta1, linear_theta2, linear_theta3};

  /// This is E0 = 1/N dN
  static Einstein::TOneFormType
  linear_E0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  linear_E1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(1.0*x(0) + 1.0*x(1), 0, 0);
                            }; 

  static Einstein::TOneFormType
  linear_E2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 1.0*x(1) + 1.0*x(2), 0);
                            }; 

  static Einstein::TOneFormType
  linear_E3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 1.0*x(0) + 1.0*x(2));
                            }; 

  static Einstein::TCollection1Forms linear_E = {minkowski_E0, linear_E1, linear_E2, linear_E3};

  //------------------------------------------------------------------------------
  /// This is E0 = 1/N dN
  static Einstein::TOneFormType
  deth_linear_E0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(t*(-2*t*x(0) - t*x(1) - t*x(2) - 2)/(pow(t, 2)*pow(x(0), 2) + pow(t, 2)*x(0)*x(1) + pow(t, 2)*x(0)*x(2) + pow(t, 2)*x(1)*x(2) + 2*t*x(0) + t*x(1) + t*x(2) + 1), t*(-t*x(0) - 2*t*x(1) - t*x(2) - 2)/(pow(t, 2)*x(0)*x(1) + pow(t, 2)*x(0)*x(2) + pow(t, 2)*pow(x(1), 2) + pow(t, 2)*x(1)*x(2) + t*x(0) + 2*t*x(1) + t*x(2) + 1), t*(-t*x(0) - t*x(1) - 2*t*x(2) - 2)/(pow(t, 2)*x(0)*x(1) + pow(t, 2)*x(0)*x(2) + pow(t, 2)*x(1)*x(2) + pow(t, 2)*pow(x(2), 2) + t*x(0) + t*x(1) + 2*t*x(2) + 1));
                            }; 

  static Einstein::TOneFormType
  deth_linear_E1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d((x(0) + x(1))/(pow(t, 3)*pow(x(0), 2)*x(1) + pow(t, 3)*pow(x(0), 2)*x(2) + pow(t, 3)*x(0)*pow(x(1), 2) + 2*pow(t, 3)*x(0)*x(1)*x(2) + pow(t, 3)*x(0)*pow(x(2), 2) + pow(t, 3)*pow(x(1), 2)*x(2) + pow(t, 3)*x(1)*pow(x(2), 2) + pow(t, 2)*pow(x(0), 2) + 3*pow(t, 2)*x(0)*x(1) + 3*pow(t, 2)*x(0)*x(2) + pow(t, 2)*pow(x(1), 2) + 3*pow(t, 2)*x(1)*x(2) + pow(t, 2)*pow(x(2), 2) + 2*t*x(0) + 2*t*x(1) + 2*t*x(2) + 1), 0, 0);
                            }; 

  static Einstein::TOneFormType
  deth_linear_E2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, (x(1) + x(2))/(pow(t, 3)*pow(x(0), 2)*x(1) + pow(t, 3)*pow(x(0), 2)*x(2) + pow(t, 3)*x(0)*pow(x(1), 2) + 2*pow(t, 3)*x(0)*x(1)*x(2) + pow(t, 3)*x(0)*pow(x(2), 2) + pow(t, 3)*pow(x(1), 2)*x(2) + pow(t, 3)*x(1)*pow(x(2), 2) + pow(t, 2)*pow(x(0), 2) + 3*pow(t, 2)*x(0)*x(1) + 3*pow(t, 2)*x(0)*x(2) + pow(t, 2)*pow(x(1), 2) + 3*pow(t, 2)*x(1)*x(2) + pow(t, 2)*pow(x(2), 2) + 2*t*x(0) + 2*t*x(1) + 2*t*x(2) + 1), 0);
                            }; 

  static Einstein::TOneFormType
  deth_linear_E3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, (x(0) + x(2))/(pow(t, 3)*pow(x(0), 2)*x(1) + pow(t, 3)*pow(x(0), 2)*x(2) + pow(t, 3)*x(0)*pow(x(1), 2) + 2*pow(t, 3)*x(0)*x(1)*x(2) + pow(t, 3)*x(0)*pow(x(2), 2) + pow(t, 3)*pow(x(1), 2)*x(2) + pow(t, 3)*x(1)*pow(x(2), 2) + pow(t, 2)*pow(x(0), 2) + 3*pow(t, 2)*x(0)*x(1) + 3*pow(t, 2)*x(0)*x(2) + pow(t, 2)*pow(x(1), 2) + 3*pow(t, 2)*x(1)*x(2) + pow(t, 2)*pow(x(2), 2) + 2*t*x(0) + 2*t*x(1) + 2*t*x(2) + 1));
                            }; 

  static Einstein::TCollection1Forms deth_linear_E = {deth_linear_E0, deth_linear_E1, deth_linear_E2, deth_linear_E3};

  //------------------------------------------------------------------------------

  static Einstein::TOneFormType
  linear_dE0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  linear_dE1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(-1.0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  linear_dE2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, -1.0);
                            }; 

  static Einstein::TOneFormType
  linear_dE3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 1.0, 0);
                            }; 

  static Einstein::TCollection1Forms linear_dE = {linear_dE0, linear_dE1, linear_dE2, linear_dE3};

  //------------------------------------------------------------------------------

  static Einstein::TTwoFormType
  linear_D0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(t*(t*x(0) + t*x(1) + 1)/(t*x(0) + t*x(2) + 1), t*(-t*x(0) - t*x(2) - 1)/(t*x(1) + t*x(2) + 1), t*(t*x(1) + t*x(2) + 1)/(t*x(0) + t*x(1) + 1));
                            }; 

  static Einstein::TTwoFormType
  linear_D1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 2.0*t*x(0)*x(1) + 2.0*t*x(0)*x(2) + 2.0*t*x(1)*x(2) + 2.0*t*pow(x(2), 2) + 1.0*x(0) + 1.0*x(1) + 2.0*x(2));
                            }; 

  static Einstein::TTwoFormType
  linear_D2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, -2.0*t*pow(x(0), 2) - 2.0*t*x(0)*x(1) - 2.0*t*x(0)*x(2) - 2.0*t*x(1)*x(2) - 2.0*x(0) - 1.0*x(1) - 1.0*x(2), 0);
                            }; 

  static Einstein::TTwoFormType
  linear_D3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(2.0*t*x(0)*x(1) + 2.0*t*x(0)*x(2) + 2.0*t*pow(x(1), 2) + 2.0*t*x(1)*x(2) + 1.0*x(0) + 2.0*x(1) + 1.0*x(2), 0, 0);
                            }; 

  static Einstein::TCollection2Forms linear_D = {linear_D0, linear_D1, linear_D2, linear_D3};

  //------------------------------------------------------------------------------

  /// Not calculated since not needed
  static Einstein::TTwoFormType
  linear_dtD0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  linear_dtD1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 2.0*x(0)*x(1) + 2.0*x(0)*x(2) + 2.0*x(1)*x(2) + 2.0*pow(x(2), 2));
                            }; 

  static Einstein::TTwoFormType
  linear_dtD2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, -2.0*pow(x(0), 2) - 2.0*x(0)*x(1) - 2.0*x(0)*x(2) - 2.0*x(1)*x(2), 0);
                            }; 

  static Einstein::TTwoFormType
  linear_dtD3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(2.0*x(0)*x(1) + 2.0*x(0)*x(2) + 2.0*pow(x(1), 2) + 2.0*x(1)*x(2), 0, 0);
                            }; 

  static Einstein::TCollection2Forms linear_dtD = {linear_dtD0, linear_dtD1, linear_dtD2, linear_dtD3};

//------------------------------------------------------------------------------

  static Einstein::TTwoFormType
  linear_B0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  linear_B1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(-t, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  linear_B2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, -t);
                            }; 

  static Einstein::TTwoFormType
  linear_B3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, t, 0);
                            }; 

  static Einstein::TCollection2Forms linear_B = {linear_B0, linear_B1, linear_B2, linear_B3};

  static Einstein::TTwoFormType
  linear_dtB0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  linear_dtB1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(-1, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  linear_dtB2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, -1);
                            }; 

  static Einstein::TTwoFormType
  linear_dtB3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 1, 0);
                            }; 

  static Einstein::TCollection2Forms linear_dtB = {linear_dtB0, linear_dtB1, linear_dtB2, linear_dtB3};

  static Einstein::TOneFormType
  linear_H0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  linear_H1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, -t/(t*x(0) + t*x(2) + 1), 0);
                            }; 

  static Einstein::TOneFormType
  linear_H2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, -t/(t*x(0) + t*x(1) + 1));
                            }; 

  static Einstein::TOneFormType
  linear_H3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(-t/(t*x(1) + t*x(2) + 1), 0, 0);
                            }; 

  static Einstein::TCollection1Forms linear_H = {linear_H0, linear_H1, linear_H2, linear_H3};

  static Einstein::TTwoFormType
  linear_dH0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  linear_dH1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(pow(t, 2)/pow(t*x(0) + t*x(2) + 1, 2), 0, -pow(t, 2)/(pow(t, 2)*pow(x(0), 2) + 2*pow(t, 2)*x(0)*x(2) + pow(t, 2)*pow(x(2), 2) + 2*t*x(0) + 2*t*x(2) + 1));
                            }; 

  static Einstein::TTwoFormType
  linear_dH2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, pow(t, 2)/(pow(t, 2)*pow(x(0), 2) + 2*pow(t, 2)*x(0)*x(1) + pow(t, 2)*pow(x(1), 2) + 2*t*x(0) + 2*t*x(1) + 1), pow(t, 2)/pow(t*x(0) + t*x(1) + 1, 2));
                            }; 

  static Einstein::TTwoFormType
  linear_dH3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(-pow(t, 2)/(pow(t, 2)*pow(x(1), 2) + 2*pow(t, 2)*x(1)*x(2) + pow(t, 2)*pow(x(2), 2) + 2*t*x(1) + 2*t*x(2) + 1), pow(t, 2)*(-t*(x(0) + x(1)) - 1)/((t*x(0) + t*x(1) + 1)*pow(t*x(1) + t*x(2) + 1, 2)), 0);
                            }; 

  static Einstein::TCollection2Forms linear_dH = {linear_dH0, linear_dH1, linear_dH2, linear_dH3};

  //------------------------------------------------------------------------------
  // Constant solution (D, theta, B)
  //------------------------------------------------------------------------------

  static Einstein::TZeroFormType
  constant_lapse = [](const double t, const Eigen::Vector3d & x) -> double {
                        return 6*pow(t, 3) + 26*pow(t, 2) + 26*t + 6;
                      };

  static Einstein::TOneFormType
  constant_theta0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, 0, 0);
                      };

  static Einstein::TOneFormType
  constant_theta1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, 0, 3*t+1);
                      };

  static Einstein::TOneFormType
  constant_theta2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(2*t+2, 0, 0);
                      };

  static Einstein::TOneFormType
  constant_theta3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, t+3, 0);
                      };

  static Einstein::TCollection1Forms constant_theta = {constant_theta0, constant_theta1, constant_theta2, constant_theta3};

  //  This is E0 = 1/N dN
  static Einstein::TOneFormType
  constant_E0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  constant_E1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 3.0);
                            }; 

  static Einstein::TOneFormType
  constant_E2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(2.0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  constant_E3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 1.0, 0);
                            }; 

  static Einstein::TCollection1Forms constant_E = {constant_E0, constant_E1, constant_E2, constant_E3};

  static Einstein::TOneFormType
  constant_dE0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  constant_dE1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  constant_dE2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  constant_dE3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection1Forms constant_dE = {constant_dE0, constant_dE1, constant_dE2, constant_dE3};

  static Einstein::TTwoFormType
  constant_D0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  constant_D1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(4.0*t + 8.0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  constant_D2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 6.0*t + 10.0);
                            }; 

  static Einstein::TTwoFormType
  constant_D3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, -12.0*t - 8.0, 0);
                            }; 

  static Einstein::TCollection2Forms constant_D = {constant_D0, constant_D1, constant_D2, constant_D3};

  static Einstein::TTwoFormType
  constant_dtD0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  constant_dtD1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(4.0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  constant_dtD2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 6.0);
                            }; 

  static Einstein::TTwoFormType
  constant_dtD3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, -12.0, 0);
                            }; 

  static Einstein::TCollection2Forms constant_dtD = {constant_dtD0, constant_dtD1, constant_dtD2, constant_dtD3};

  static Einstein::TTwoFormType
  constant_B0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  constant_B1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  constant_B2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  constant_B3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection2Forms constant_B = {constant_B0, constant_B1, constant_B2, constant_B3};

//------------------------------------------------------------------------------

  static Einstein::TTwoFormType
  constant_dtB0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  constant_dtB1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  constant_dtB2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  constant_dtB3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection2Forms constant_dtB = {constant_dtB0, constant_dtB1, constant_dtB2, constant_dtB3};

//------------------------------------------------------------------------------

  static Einstein::TOneFormType
  constant_H0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  constant_H1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  constant_H2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  constant_H3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection1Forms constant_H = {constant_H0, constant_H1, constant_H2, constant_H3};

  static Einstein::TTwoFormType
  constant_dH0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  constant_dH1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  constant_dH2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  constant_dH3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection2Forms constant_dH = {constant_dH0, constant_dH1, constant_dH2, constant_dH3};

  //------------------------------------------------------------------------------
  // Test solution (D, theta, B)
  //------------------------------------------------------------------------------

  static Einstein::TZeroFormType
  test_lapse = [](const double t, const Eigen::Vector3d & x) -> double {
                        return pow(t, -0.25)*exp(0.5*pow(M_PI, 2)*pow(t, 2)*(pow(cyl_bessel_j(0, 2*M_PI*t), 2) + pow(cyl_bessel_j(1, 2*M_PI*t), 2)) - 0.5*M_PI*t*pow(cos(2*M_PI*x(2)), 2)*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t) - 0.5*pow(M_PI, 2)*(pow(cyl_bessel_j(1, 2*M_PI), 2) + pow(cyl_bessel_j(0, 2*M_PI), 2)) + 0.25*M_PI*cyl_bessel_j(0, 2*M_PI)*cyl_bessel_j(1, 2*M_PI));
                      };

  static Einstein::TOneFormType
  test_theta0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, 0, 0);
                      };

  static Einstein::TOneFormType
  test_theta1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(1, 0, 0);
                      };

  static Einstein::TOneFormType
  test_theta2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, 1, 0);
                      };

  static Einstein::TOneFormType
  test_theta3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                        return Eigen::Vector3d(0, 0, 1);
                      };

  static Einstein::TCollection1Forms test_theta = {test_theta0, test_theta1, test_theta2, test_theta3};

  //  This is E0 = 1/N dN
  static Einstein::TOneFormType
  test_E0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, -1.0*pow(M_PI, 2)*pow(t, 1.0)*sin(4*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t));
                            }; 

  static Einstein::TOneFormType
  test_E1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  test_E2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  test_E3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection1Forms test_E = {test_E0, test_E1, test_E2, test_E3};

  static Einstein::TOneFormType
  test_dE0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  test_dE1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  test_dE2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  test_dE3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection1Forms test_dE = {test_dE0, test_dE1, test_dE2, test_dE3};

  static Einstein::TTwoFormType
  test_D0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                             return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  test_D1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  test_D2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  test_D3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection2Forms test_D = {test_D0, test_D1, test_D2, test_D3};

  static Einstein::TTwoFormType
  test_dtD0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  test_dtD1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  test_dtD2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  test_dtD3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection2Forms test_dtD = {test_dtD0, test_dtD1, test_dtD2, test_dtD3};

  static Einstein::TTwoFormType
  test_B0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  test_B1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  test_B2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  test_B3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection2Forms test_B = {test_B0, test_B1, test_B2, test_B3};

  static Einstein::TTwoFormType
  test_dtB0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  test_dtB1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  test_dtB2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  test_dtB3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection2Forms test_dtB = {test_dtB0, test_dtB1, test_dtB2, test_dtB3};

//------------------------------------------------------------------------------

  static Einstein::TOneFormType
  test_H0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TOneFormType
  test_H1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, -1.0*pow(M_PI, 2)*pow(t, 1.0)*sin(4*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t), 0);
                            }; 

  static Einstein::TOneFormType
  test_H2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(1.0*pow(M_PI, 2)*pow(t, 1.0)*sin(4*M_PI*x(2))*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t), 0, 0);
                            }; 

  static Einstein::TOneFormType
  test_H3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection1Forms test_H = {test_H0, test_H1, test_H2, test_H3};

  static Einstein::TTwoFormType
  test_dH0 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TTwoFormType
  test_dH1 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 4.0*pow(M_PI, 3)*pow(t, 1.0)*(-pow(sin(2*M_PI*x(2)), 2) + pow(cos(2*M_PI*x(2)), 2))*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t));
                            }; 

  static Einstein::TTwoFormType
  test_dH2 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 4.0*pow(M_PI, 3)*pow(t, 1.0)*(pow(sin(2*M_PI*x(2)), 2) - pow(cos(2*M_PI*x(2)), 2))*cyl_bessel_j(0, 2*M_PI*t)*cyl_bessel_j(1, 2*M_PI*t), 0);
                            }; 

  static Einstein::TTwoFormType
  test_dH3 = [](const double t, const Eigen::Vector3d & x) -> Eigen::Vector3d {
                              return Eigen::Vector3d(0, 0, 0);
                            }; 

  static Einstein::TCollection2Forms test_dH = {test_dH0, test_dH1, test_dH2, test_dH3};

//------------------------------------------------------------------------------

} // end of namespace HArDCore3D

#endif
