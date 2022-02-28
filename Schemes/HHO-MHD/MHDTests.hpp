// Class to provide various test cases for the Stokes problem
//
//
// Author: Liam Yemm (liam.yemm@monash.edu)
//

#include <mesh.hpp>
// #include <Math.hpp>
#include "basis.hpp"

using Math::PI;

#ifndef _MHD_TESTS_HPP
#define _MHD_TESTS_HPP


namespace HArDCore3D
{

  /*!
   * @addtogroup HHO_MHD
   * @{
   */
   
    class MHDTests
    {
    public:
        /// Initialise data
        MHDTests(
            size_t u_id, ///< The id of the test case for the velocity
            size_t b_id,  ///< The id of the test case for the magnetic field
            size_t p_id,  ///< The id of the test case for the pressure
            double visc,  ///< The viscocity coefficient
            double diff  ///< The magnetic diffusivity coefficient
        );

        FType<VectorRd> velocity();

        FType<MatrixRd> grad_velocity();

        FType<VectorRd> laplace_velocity();

        FType<VectorRd> magnetic();

        FType<MatrixRd> grad_magnetic();

        FType<VectorRd> laplace_magnetic();

        FType<VectorRd> u_dot_grad_u();
        FType<VectorRd> u_dot_grad_b();
        FType<VectorRd> b_dot_grad_u();
        FType<VectorRd> b_dot_grad_b();

        FType<double> pressure();

        FType<VectorRd> grad_pressure();

        FType<VectorRd> magnetic_source();
        FType<VectorRd> velocity_source();

    private:
        size_t m_u_id;
        size_t m_b_id;
        size_t m_p_id;
        double m_visc;
        double m_diff;

        /// Check if the provided test cases are valid (within range, and combination of solution/diffusion valid)
        void validate();
    };
}

//@}

#endif
