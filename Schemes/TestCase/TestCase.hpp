// Class to provide various test cases (diffusion, exact solution, and their derivatives)
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
*
*	This library was developed around HHO methods, although some parts of it have a more
* general purpose. If you use this code or part of it in a scientific publication, 
* please mention the following book as a reference for the underlying principles
* of HHO schemes:
*
* The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications.
* D. A. Di Pietro and J. Droniou. 2019, 516p. 
* url: https://hal.archives-ouvertes.fr/hal-02151813.
*
*/

#ifndef _TEST_CASE_HPP
#define _TEST_CASE_HPP

#include "basis.hpp"


using namespace HArDCore3D;

/*!
* @defgroup TestCases
*	@brief Defines test cases (exact solutions, source terms etc.)
*/

// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

// @addtogroup TestCases
//@{



/// The TestCase class provides definition of test cases
class TestCase {

public:
  /// Initialise data
    TestCase(
        std::vector<int> iTC ///< The vector id of the test case: (id of solution, id of diffusion)
    );

    /// Returns the exact solution
    FType<double> sol();

    /// Returns the gradient of the exact solution
    CellFType<VectorRd> grad_sol();

    /// Returns the Hessian of the exact solution
    CellFType<MatrixRd> hess_sol();

    /// Returns the diffusion matrix
    CellFType<MatrixRd> diff();

    /// Returns the divergence by row of the diffusion matrix
    CellFType<VectorRd> div_diff();

    /// Returns the reaction term
    CellFType<double> reac();

    /// Returns the advection term
    CellFType<VectorRd> advec();

    /// Returns the divergence of the advection
    CellFType<double> div_advec();

    /// Returns the diffusion source term
    CellFType<double> diff_source();

    /// Returns the reaction-diffusion source term
    CellFType<double> diff_advec_reac_source();

    /// Returns the value of the parameter lambda
    inline double get_lambda();

    /// Check if the provided test cases are valid (within range, and combination of solution/diffusion valid)
    void validate();

    /// Returns the degree of the diffusion tensor (useful to set up quadrature rules of proper degree)
    inline size_t get_deg_diff();

    /// Returns a flag to check if reaction is constant
    inline bool is_reac_const();

    /// Returns a flag to check if divergence of advection is zero
    inline bool is_div_advec_zero();

    /// Returns a flag to check if divergence of advection is constant
    inline bool is_div_advec_const();

private:
  // Parameters: id of test case, pi
  std::vector<int> m_iTC;
	const double pi = acos(-1);
  const double eps = 1e-5;

	size_t _deg_diff;

  bool reac_const = true;
  bool div_advec_zero = true;
  bool div_advec_const = true;
};

inline size_t TestCase::get_deg_diff() { return _deg_diff; };

inline bool TestCase::is_reac_const() { return reac_const; };
inline bool TestCase::is_div_advec_zero() { return div_advec_zero; };
inline bool TestCase::is_div_advec_const() { return div_advec_const; };

//@}

#endif //_TEST_CASE_HPP
