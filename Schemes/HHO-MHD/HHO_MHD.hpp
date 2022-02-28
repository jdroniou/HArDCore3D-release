
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "MHDTests.hpp"
#include "mesh_builder.hpp"
#include "vtu_writer.hpp"
#include <basis.hpp>
#include <chrono>
#include <elementquad.hpp>
#include <hybridcore.hpp>
#include <mesh.hpp>
#include <parallel_for.hpp>

#include <BoundaryConditions/BoundaryConditions.hpp>
#include <TestCase/TestCase.hpp>

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

#ifdef WITH_MKL
#include <Eigen/PardisoSupport>
#endif

/*!
 * @defgroup HHO_MHD
 * @brief Implementation of the HHO scheme for the MHD problem 
 */
 
 
using namespace HArDCore3D;

#ifndef _HHO_MHD_HPP
#define _HHO_MHD_HPP

  /*!
   * @addtogroup HHO_MHD
   * @{
   */


struct SolutionVector
{
public:
    // SolutionVector() : m_cell_deg(0), m_face_deg(0)
    // {
    // }

    SolutionVector(Eigen::VectorXd values, ///< values of the vector
                   const Mesh *mesh_ptr,   ///< reference to the mesh
                   const size_t cell_deg,  ///< polynomial degrees in cell
                   const size_t face_deg   ///< polynomial degrees on face
                   )
        : m_values(values), m_mesh_ptr(mesh_ptr), m_cell_deg(cell_deg),
          m_face_deg(face_deg)
    {
        assert((size_t)(m_values.size()) == n_total_dofs);
    }

    /// Return the values as an Eigen vector
    inline Eigen::VectorXd asVectorXd() const { return m_values; }
    /// Return the values as an Eigen vector
    inline void set_values(Eigen::VectorXd values) const { m_values = values; }

    /// Return the cell degree
    inline const size_t get_cell_deg() const { return m_cell_deg; }

    /// Return the face degree
    inline const size_t get_face_deg() const { return m_face_deg; }

    Eigen::VectorXd velocity_values() const
    {
        return m_values.head(n_total_cell_velocity_dofs +
                             n_total_face_velocity_dofs);
    }

    Eigen::VectorXd pressure_values() const
    {
        return m_values.segment(n_total_cell_velocity_dofs + n_total_face_velocity_dofs, n_total_pressure_dofs);
    }

    Eigen::VectorXd magnetic_values() const
    {
        return m_values.segment( n_total_cell_velocity_dofs + n_total_face_velocity_dofs + n_total_pressure_dofs, n_total_cell_magnetic_dofs + n_total_face_magnetic_dofs);
    }

    Eigen::VectorXd lagrange_values() const
    {
        return m_values.tail(n_total_lagrange_dofs);
    }

    Eigen::VectorXd velocity_restr(size_t iT) const
    {
        Eigen::VectorXd XTF = Eigen::VectorXd::Zero(n_local_velocity_dofs(iT));

        XTF.head(n_local_cell_velocity_dofs) = m_values.segment( iT * n_local_cell_velocity_dofs, n_local_cell_velocity_dofs);
        for (size_t iTF = 0; iTF < m_mesh_ptr->cell(iT)->n_faces(); iTF++)
        {
            size_t global_offset = n_total_cell_velocity_dofs + m_mesh_ptr->cell(iT)->face(iTF)->global_index() * n_local_face_velocity_dofs;
            size_t local_offset = n_local_cell_velocity_dofs + iTF * n_local_face_velocity_dofs;
            XTF.segment(local_offset, n_local_face_velocity_dofs) = m_values.segment(global_offset, n_local_face_velocity_dofs);
        }

        return XTF;
    }

    Eigen::VectorXd pressure_restr(size_t iT) const
    {
        size_t offset = n_total_cell_velocity_dofs + n_total_face_velocity_dofs + iT * n_local_pressure_dofs;

        return m_values.segment(offset, n_local_pressure_dofs);
    }

    Eigen::VectorXd magnetic_restr(size_t iT) const
    {
        Eigen::VectorXd XTF = Eigen::VectorXd::Zero(n_local_magnetic_dofs(iT));

        size_t offset = n_total_cell_velocity_dofs + n_total_face_velocity_dofs + n_total_pressure_dofs;

        XTF.head(n_local_cell_magnetic_dofs) = m_values.segment( offset + iT * n_local_cell_magnetic_dofs, n_local_cell_magnetic_dofs);
        for (size_t iTF = 0; iTF < m_mesh_ptr->cell(iT)->n_faces(); iTF++)
        {
            size_t global_offset = offset + n_total_cell_magnetic_dofs + m_mesh_ptr->cell(iT)->face(iTF)->global_index() * n_local_face_magnetic_dofs;
            size_t local_offset = n_local_cell_magnetic_dofs + iTF * n_local_face_magnetic_dofs;
            XTF.segment(local_offset, n_local_face_magnetic_dofs) = m_values.segment(global_offset, n_local_face_magnetic_dofs);
        }

        return XTF;
    }

    Eigen::VectorXd lagrange_restr(size_t iT) const
    {
        size_t offset = n_total_cell_velocity_dofs + n_total_face_velocity_dofs + n_total_pressure_dofs + n_total_cell_magnetic_dofs + n_total_face_magnetic_dofs + iT * n_local_lagrange_dofs;

        return m_values.segment(offset, n_local_lagrange_dofs);
    }

    Eigen::VectorXd restr(size_t iT) const
    {
        Eigen::VectorXd XTF = Eigen::VectorXd::Zero( n_local_velocity_dofs(iT) + n_local_pressure_dofs + n_local_magnetic_dofs(iT) + n_local_lagrange_dofs);

        XTF.head(n_local_velocity_dofs(iT)) = this->velocity_restr(iT);
        XTF.segment(n_local_velocity_dofs(iT), n_local_pressure_dofs) = this->pressure_restr(iT);
        XTF.segment(n_local_velocity_dofs(iT) + n_local_pressure_dofs, n_local_magnetic_dofs(iT)) = this->magnetic_restr(iT);
        XTF.tail(n_local_lagrange_dofs) = this->lagrange_restr(iT);

        return XTF;
    }

    /// Overloads the addition: adds the coefficients
    SolutionVector operator+(const SolutionVector &b) const
    {
        assert(m_cell_deg == b.get_cell_deg() || m_face_deg == b.get_face_deg());
        return SolutionVector(m_values + b.asVectorXd(), m_mesh_ptr, m_cell_deg, m_face_deg);
    }

    /// Overloads the subtraction: subtracts the coefficients
    SolutionVector operator-(const SolutionVector &b) const
    {
        assert(m_cell_deg == b.get_cell_deg() || m_face_deg == b.get_face_deg());
        return SolutionVector(m_values - b.asVectorXd(), m_mesh_ptr, m_cell_deg, m_face_deg);
    }

    /// Overloads the (): returns the corresponding coefficient
    double operator()(size_t index) const { return m_values(index); }

    // //Assignment operator
    // SolutionVector &operator=(const SolutionVector &a);

private:
    mutable Eigen::VectorXd m_values;
    const Mesh *m_mesh_ptr;

    const size_t m_cell_deg;
    const size_t m_face_deg;

    static const size_t dim = 3;

    const size_t n_cells = m_mesh_ptr->n_cells();
    const size_t n_faces = m_mesh_ptr->n_faces();

    const size_t n_local_cell_velocity_dofs = dim * DimPoly<Cell>(m_cell_deg);
    const size_t n_local_cell_magnetic_dofs = dim * DimPoly<Cell>(m_cell_deg);
    const size_t n_local_face_velocity_dofs = dim * DimPoly<Face>(m_face_deg);
    const size_t n_local_face_magnetic_dofs = dim * DimPoly<Face>(m_face_deg);
    const size_t n_local_pressure_dofs = DimPoly<Cell>(m_cell_deg);
    const size_t n_local_lagrange_dofs = DimPoly<Cell>(m_cell_deg);

    const size_t n_total_cell_velocity_dofs = n_cells * n_local_cell_velocity_dofs;
    const size_t n_total_cell_magnetic_dofs = n_cells * n_local_cell_magnetic_dofs;
    const size_t n_total_face_velocity_dofs = n_faces * n_local_face_velocity_dofs;
    const size_t n_total_face_magnetic_dofs = n_faces * n_local_face_magnetic_dofs;
    const size_t n_total_pressure_dofs = n_cells * n_local_pressure_dofs;
    const size_t n_total_lagrange_dofs = n_cells * n_local_lagrange_dofs;

    const size_t n_total_dofs = n_total_cell_velocity_dofs + n_total_face_velocity_dofs + n_total_pressure_dofs + n_total_cell_magnetic_dofs + n_total_face_magnetic_dofs + n_total_lagrange_dofs;

    inline size_t n_local_velocity_dofs(const size_t iT) const
    {
        return n_local_cell_velocity_dofs + m_mesh_ptr->cell(iT)->n_faces() * n_local_face_velocity_dofs;
    }

    inline size_t n_local_magnetic_dofs(const size_t iT) const
    {
        return n_local_cell_magnetic_dofs + m_mesh_ptr->cell(iT)->n_faces() * n_local_face_magnetic_dofs;
    }
};

class MHDModel
{
public:
    MHDModel(HybridCore &, size_t, size_t, char, char);

    SolutionVector solve_with_static_cond(FType<VectorRd> f_source, FType<VectorRd> g_source, double visc, double diff, double tol, bool threading = true);
    SolutionVector global_interpolant(FType<VectorRd> velocity, FType<double> pressure, FType<VectorRd> magnetic);
    std::vector<double> compute_errors(SolutionVector intepolant, SolutionVector discrete, double visc, double diff);

private:
    Eigen::MatrixXd gradient_term(Cell *cell, ElementQuad &elquad);
    Eigen::MatrixXd stabilisation_term(Cell *cell, ElementQuad &elquad);
    Eigen::MatrixXd divergence_term(Cell *cell, ElementQuad &elquad);
    Eigen::VectorXd source_term(FType<VectorRd> f_vec, Cell *cell, ElementQuad &elquad);
    Eigen::MatrixXd convective_term(Cell *cell, std::vector<Eigen::MatrixXd> tplus_k, std::vector<Eigen::MatrixXd> tminus_k, Eigen::VectorXd u);
    Eigen::MatrixXd th_ijk_plus_th_ikj(Cell *cell, ElementQuad &elquad, size_t k);
    Eigen::MatrixXd th_ijk_minus_th_ikj(Cell *cell, ElementQuad &elquad, size_t k);

    Eigen::VectorXd local_interpolant(FType<VectorRd> func, Cell *cell, ElementQuad &elquad);
    Eigen::VectorXd scalar_l2_projection(FType<double> func, Cell *cell, ElementQuad &elquad);

    HybridCore &m_hho;
    size_t m_L;
    size_t m_K;
    char m_u_bc;
    char m_b_bc;

    static const size_t dim = 3;

    const Mesh *mesh_ptr = m_hho.get_mesh();

    const size_t n_cells = mesh_ptr->n_cells();
    const size_t n_faces = mesh_ptr->n_faces();
    const size_t n_bdry_faces = mesh_ptr->n_b_faces();

    const size_t n_local_cell_scalar_dofs = DimPoly<Cell>(m_L);
    const size_t n_local_face_scalar_dofs = DimPoly<Face>(m_K);
    const size_t n_local_highorder_scalar_dofs = DimPoly<Cell>(m_K + 1);

    const size_t n_local_cell_vector_dofs = dim * DimPoly<Cell>(m_L);
    const size_t n_local_face_vector_dofs = dim * DimPoly<Face>(m_K);
    const size_t n_local_highorder_vector_dofs = dim * DimPoly<Cell>(m_K + 1);

    const size_t n_total_cell_scalar_dofs = n_local_cell_scalar_dofs * n_cells;
    // const size_t n_total_face_scalar_dofs = n_local_face_scalar_dofs * n_faces;

    const size_t n_total_cell_vector_dofs = n_local_cell_vector_dofs * n_cells;
    const size_t n_total_face_vector_dofs = n_local_face_vector_dofs * n_faces;
    // const size_t n_total_vector_dofs = n_total_cell_vector_dofs +
    // n_total_face_vector_dofs;

    size_t n_unknowns = 0;

    Eigen::SparseMatrix<double> inv_LTll_LTlg;
    Eigen::SparseMatrix<double> LTgl;
    Eigen::SparseMatrix<double> LTgg;

    Eigen::VectorXd inv_LTll_rTl;

    std::vector<Eigen::MatrixXd> AT;
    std::vector<Eigen::MatrixXd> RT;
};

std::vector<double> MHDModel::compute_errors(SolutionVector intepolant, SolutionVector discrete, double visc, double diff)
{

    double vel_energy_difference = 0.0;
    double vel_energy_exact = 0.0;

    double mag_energy_difference = 0.0;
    double mag_energy_exact = 0.0;

    double pressure_energy_difference = 0.0;
    double pressure_energy_exact = 0.0;

    double lagrange_energy_difference = 0.0;
    double lagrange_energy_exact = 0.0;

    double vel_l2_difference = 0.0;
    double vel_l2_exact = 0.0;

    double mag_l2_difference = 0.0;
    double mag_l2_exact = 0.0;

    // Construct the local matrices using multithreading if use_threads is true
    // std::function<void(size_t, size_t)> compute_local_errors = [&](size_t
    // start, size_t end) -> void
    // {
    //     for (size_t iT = start; iT < end; ++iT)
    //     {
    for (size_t iT = 0; iT < n_cells; ++iT)
    {
        Eigen::VectorXd local_vel_difference = (intepolant - discrete).velocity_restr(iT);
        Eigen::VectorXd local_vel_exact = intepolant.velocity_restr(iT);
        Eigen::VectorXd local_mag_difference = (intepolant - discrete).magnetic_restr(iT);
        Eigen::VectorXd local_mag_exact = intepolant.magnetic_restr(iT);

        vel_energy_difference += local_vel_difference.transpose() * AT[iT] * local_vel_difference;
        mag_energy_difference += local_mag_difference.transpose() * AT[iT] * local_mag_difference;
        vel_energy_exact += local_vel_exact.transpose() * AT[iT] * local_vel_exact;
        mag_energy_exact += local_mag_exact.transpose() * AT[iT] * local_mag_exact;

        ElementQuad elquad(m_hho, iT, m_L + m_L, m_K + m_K);

        QuadratureRule quadT = elquad.get_quadT();

        BasisQuad<double> phiT_quadT = elquad.get_phiT_quadT();
        BasisQuad<VectorRd> dphiT_quadT = elquad.get_dphiT_quadT();

        Eigen::MatrixXd MTT = compute_gram_matrix(phiT_quadT, phiT_quadT, quadT, n_local_cell_scalar_dofs, n_local_cell_scalar_dofs, "sym");

        Eigen::MatrixXd vec_MTT = Eigen::MatrixXd::Zero(n_local_cell_vector_dofs, n_local_cell_vector_dofs);

        for (size_t d = 0; d < dim; d++)
        {
            for (size_t i = 0; i < n_local_cell_scalar_dofs; i++)
            {
                for (size_t j = 0; j < n_local_cell_scalar_dofs; j++)
                {
                    vec_MTT(dim * i + d, dim * j + d) = MTT(i, j);
                }
            }
        }

        vel_l2_difference += local_vel_difference.head(n_local_cell_vector_dofs).transpose() * vec_MTT * local_vel_difference.head(n_local_cell_vector_dofs);
        mag_l2_difference += local_mag_difference.head(n_local_cell_vector_dofs).transpose() * vec_MTT * local_mag_difference.head(n_local_cell_vector_dofs);
        vel_l2_exact += local_vel_exact.head(n_local_cell_vector_dofs).transpose() * vec_MTT * local_vel_exact.head(n_local_cell_vector_dofs);
        mag_l2_exact += local_mag_exact.head(n_local_cell_vector_dofs).transpose() * vec_MTT * local_mag_exact.head(n_local_cell_vector_dofs);

        Cell *cell = mesh_ptr->cell(iT);

        size_t n_local_bdry_vector_dofs = cell->n_faces() * n_local_face_vector_dofs;
        Eigen::MatrixXd M_BDRY_BDRY = Eigen::MatrixXd::Zero(n_local_bdry_vector_dofs, n_local_bdry_vector_dofs);

        for (size_t iTF = 0; iTF < cell->n_faces(); iTF++)
        {
            const size_t offset_F = iTF * n_local_face_vector_dofs;

            QuadratureRule quadF = elquad.get_quadF(iTF);
            BasisQuad<double> phiF_quadF = elquad.get_phiF_quadF(iTF);

            BasisQuad<VectorRd> vec_phiF_quadF(
                boost::extents[n_local_face_vector_dofs][quadF.size()]);

            for (size_t d = 0; d < dim; d++)
            {
                VectorRd eF = VectorRd::Zero();

                // Sort face unknowns into perpendicular, then paralell components
                if (d == 0)
                {
                    eF = cell->face(iTF)->normal();
                }
                else if (d == 1)
                {
                    eF = cell->face(iTF)->edge(0)->tangent();
                }
                else if (d == 2)
                {
                    // eF = mesh_ptr->cell(iT)->face(iTF)->edge(1)->tangent();
                    eF = (cell->face(iTF)->normal()).cross(cell->face(iTF)->edge(0)->tangent());
                }
                else
                {
                    assert(false); // should never reach here
                }

                for (size_t iqn = 0; iqn < quadF.size(); iqn++)
                {
                    for (size_t i = 0; i < n_local_face_scalar_dofs; i++)
                    {
                        vec_phiF_quadF[dim * i + d][iqn] = phiF_quadF[i][iqn] * eF;
                    }
                }
            }

            M_BDRY_BDRY.block(offset_F, offset_F, n_local_face_vector_dofs, n_local_face_vector_dofs) = compute_gram_matrix(vec_phiF_quadF, vec_phiF_quadF, quadF, n_local_face_vector_dofs, n_local_face_vector_dofs, "sym");
        }

        double hT = cell->diam();
        vel_l2_difference += hT * local_vel_difference.tail(n_local_bdry_vector_dofs).transpose() * M_BDRY_BDRY * local_vel_difference.tail(n_local_bdry_vector_dofs);
        mag_l2_difference += hT * local_mag_difference.tail(n_local_bdry_vector_dofs).transpose() * M_BDRY_BDRY * local_mag_difference.tail(n_local_bdry_vector_dofs);
        vel_l2_exact += hT * local_vel_exact.tail(n_local_bdry_vector_dofs).transpose() * M_BDRY_BDRY * local_vel_exact.tail(n_local_bdry_vector_dofs);
        mag_l2_exact += hT * local_mag_exact.tail(n_local_bdry_vector_dofs).transpose() * M_BDRY_BDRY * local_mag_exact.tail(n_local_bdry_vector_dofs);

        Eigen::VectorXd local_pressure_difference = (intepolant - discrete).pressure_restr(iT);
        Eigen::VectorXd local_pressure_exact = intepolant.pressure_restr(iT);
        Eigen::VectorXd local_lagrange_difference = (intepolant - discrete).lagrange_restr(iT);
        Eigen::VectorXd local_lagrange_exact = intepolant.lagrange_restr(iT);

        pressure_energy_difference += local_pressure_difference.transpose() * MTT * local_pressure_difference;
        lagrange_energy_difference += local_lagrange_difference.transpose() * MTT * local_lagrange_difference;
        pressure_energy_exact += local_pressure_exact.transpose() * MTT * local_pressure_exact;
        lagrange_energy_exact += local_lagrange_exact.transpose() * MTT * local_lagrange_exact;
    }
    // };

    // // Running the local constructions in parallel
    // parallel_for(n_cells, compute_local_errors, threading);

    // If exact is zero, give absolute rather than relative errors
    if (vel_energy_exact < 1E-12)
    {
        vel_energy_exact = 1.0;
    }
    if (mag_energy_exact < 1E-12)
    {
        mag_energy_exact = 1.0;
    }
    if (pressure_energy_exact < 1E-12)
    {
        pressure_energy_exact = 1.0;
    }
    if (lagrange_energy_exact < 1E-12)
    {
        lagrange_energy_exact = 1.0;
    }
    if (vel_l2_exact < 1E-12)
    {
        vel_l2_exact = 1.0;
    }
    if (mag_l2_exact < 1E-12)
    {
        mag_l2_exact = 1.0;
    }

    // double velocity_energy_error = std::pow(visc, 0.5) * std::sqrt(vel_energy_difference / vel_energy_exact);
    // double magnetic_energy_error = std::pow(diff, 0.5) * std::sqrt(mag_energy_difference / mag_energy_exact);

    double velocity_energy_error = std::sqrt(visc * vel_energy_difference / vel_energy_exact);
    double magnetic_energy_error = std::sqrt(diff * mag_energy_difference / mag_energy_exact);
    double pressure_energy_error = std::sqrt(pressure_energy_difference / pressure_energy_exact);
    double lagrange_energy_error = std::sqrt(lagrange_energy_difference / lagrange_energy_exact);
    double velocity_l2_error = std::sqrt(vel_l2_difference / vel_l2_exact);
    double magnetic_l2_error = std::sqrt(mag_l2_difference / mag_l2_exact);

    return std::vector<double>{velocity_energy_error, magnetic_energy_error, pressure_energy_error, lagrange_energy_error, velocity_l2_error, magnetic_l2_error};
    // return std::vector<double>{velocity_energy_error, magnetic_energy_error,
    // pressure_energy_error, lagrange_energy_error};
}

SolutionVector MHDModel::global_interpolant(FType<VectorRd> velocity, FType<double> pressure, FType<VectorRd> magnetic)
{
    Eigen::VectorXd values = Eigen::VectorXd::Zero(2 * (n_total_cell_vector_dofs + n_total_face_vector_dofs + n_total_cell_scalar_dofs));
    for (size_t iT = 0; iT < n_cells; ++iT)
    {
        ElementQuad elquad(m_hho, iT, 14, 14);
        Cell *cell = mesh_ptr->cell(iT);

        Eigen::VectorXd local_Ikv = local_interpolant(velocity, cell, elquad);
        Eigen::VectorXd local_Ikm = local_interpolant(magnetic, cell, elquad);
        Eigen::VectorXd local_pikp = scalar_l2_projection(pressure, cell, elquad);

        values.segment(iT * n_local_cell_vector_dofs, n_local_cell_vector_dofs) = local_Ikv.head(n_local_cell_vector_dofs);
        for (size_t iTF = 0; iTF < cell->n_faces(); ++iTF)
        {
            values.segment(n_total_cell_vector_dofs + cell->face(iTF)->global_index() * n_local_face_vector_dofs, n_local_face_vector_dofs) = local_Ikv.segment(n_local_cell_vector_dofs + iTF * n_local_face_vector_dofs, n_local_face_vector_dofs); // inefficient -- rewriting on internal // faces
        }
        values.segment(n_total_cell_vector_dofs + n_total_face_vector_dofs + iT * n_local_cell_scalar_dofs, n_local_cell_scalar_dofs) = local_pikp;

        size_t magnetic_offset = n_total_cell_vector_dofs + n_total_face_vector_dofs + n_total_cell_scalar_dofs;
        values.segment(magnetic_offset + iT * n_local_cell_vector_dofs, n_local_cell_vector_dofs) = local_Ikm.head(n_local_cell_vector_dofs);
        for (size_t iTF = 0; iTF < cell->n_faces(); ++iTF)
        {
            values.segment(magnetic_offset + n_total_cell_vector_dofs + cell->face(iTF)->global_index() * n_local_face_vector_dofs, n_local_face_vector_dofs) = local_Ikm.segment(n_local_cell_vector_dofs + iTF * n_local_face_vector_dofs, n_local_face_vector_dofs); // inefficient -- rewriting on internal // faces
        }
    }
    return SolutionVector(values, mesh_ptr, m_L, m_K);
}

MHDModel::MHDModel(HybridCore &hho, size_t L, size_t K, char u_bc, char b_bc) : m_hho(hho), m_L(L), m_K(K), m_u_bc(u_bc), m_b_bc(b_bc)
{
    AT.resize(n_cells);
    RT.resize(n_cells);
}

Eigen::VectorXd MHDModel::scalar_l2_projection(FType<double> func, Cell *cell, ElementQuad &elquad)
{
    QuadratureRule quadT = elquad.get_quadT();

    BasisQuad<double> phiT_quadT = elquad.get_phiT_quadT();

    Eigen::MatrixXd MTT = compute_gram_matrix(phiT_quadT, phiT_quadT, quadT, n_local_cell_scalar_dofs, n_local_cell_scalar_dofs, "sym");

    Eigen::VectorXd RHS = integrate(func, phiT_quadT, quadT, n_local_cell_scalar_dofs);

    return (MTT.ldlt()).solve(RHS);
}

Eigen::VectorXd MHDModel::local_interpolant(FType<VectorRd> func, Cell *cell, ElementQuad &elquad)
{
    const size_t n_local_faces = cell->n_faces();
    const size_t n_local_vector_dofs = n_local_cell_vector_dofs + n_local_faces * n_local_face_vector_dofs;

    Eigen::VectorXd ITk = Eigen::VectorXd::Zero(n_local_vector_dofs);

    FType<double> func1 = [&](const VectorRd x) -> double { return func(x)(0); };
    FType<double> func2 = [&](const VectorRd x) -> double { return func(x)(1); };
    FType<double> func3 = [&](const VectorRd x) -> double { return func(x)(2); };

    Eigen::VectorXd piTk_func1 = scalar_l2_projection(func1, cell, elquad);
    Eigen::VectorXd piTk_func2 = scalar_l2_projection(func2, cell, elquad);
    Eigen::VectorXd piTk_func3 = scalar_l2_projection(func3, cell, elquad);

    for (size_t i = 0; i < n_local_cell_scalar_dofs; i++)
    {
        ITk(dim * i + 0) = piTk_func1(i);
        ITk(dim * i + 1) = piTk_func2(i);
        ITk(dim * i + 2) = piTk_func3(i);
    }

    for (size_t iTF = 0; iTF < n_local_faces; iTF++)
    {
        const size_t offset_F = n_local_cell_vector_dofs + iTF * n_local_face_vector_dofs;

        QuadratureRule quadF = elquad.get_quadF(iTF);
        size_t n_quadF = quadF.size();

        BasisQuad<double> phiF_quadF = elquad.get_phiF_quadF(iTF);

        BasisQuad<VectorRd> vec_phiF_quadF(boost::extents[n_local_face_vector_dofs][n_quadF]);

        for (size_t k = 0; k < dim; k++)
        {
            VectorRd eF = VectorRd::Zero();

            // Sort face unknowns into perpendicular, then paralell components
            if (k == 0)
            {
                eF = cell->face(iTF)->normal();
            }
            else if (k == 1)
            {
                eF = cell->face(iTF)->edge(0)->tangent();
            }
            else if (k == 2)
            {
                // eF = mesh_ptr->cell(iT)->face(iTF)->edge(1)->tangent();
                eF = (cell->face(iTF)->normal()).cross(cell->face(iTF)->edge(0)->tangent());
            }
            else
            {
                assert(false); // should never reach here
            }

            for (size_t iqn = 0; iqn < n_quadF; iqn++)
            {
                for (size_t i = 0; i < n_local_face_scalar_dofs; i++)
                {
                    vec_phiF_quadF[dim * i + k][iqn] = phiF_quadF[i][iqn] * eF;
                }
            }
        }

        Eigen::MatrixXd MFF = compute_gram_matrix(vec_phiF_quadF, vec_phiF_quadF, quadF, n_local_face_vector_dofs, n_local_face_vector_dofs, "sym");
        Eigen::VectorXd RHS = integrate(func, vec_phiF_quadF, quadF, n_local_face_vector_dofs);

        ITk.segment(offset_F, n_local_face_vector_dofs) = (MFF.ldlt()).solve(RHS);
    }

    return ITk;
}

Eigen::MatrixXd MHDModel::divergence_term(Cell *cell, ElementQuad &elquad)
{
    const size_t n_local_faces = cell->n_faces();
    const size_t n_local_vector_dofs = n_local_cell_vector_dofs + n_local_faces * n_local_face_vector_dofs;

    QuadratureRule quadT = elquad.get_quadT();
    size_t n_quadT = quadT.size();

    BasisQuad<double> phiT_quadT = elquad.get_phiT_quadT();
    BasisQuad<VectorRd> dphiT_quadT = elquad.get_dphiT_quadT();

    BasisQuad<double> div_vec_phiT_quadT(boost::extents[n_local_cell_vector_dofs][n_quadT]);

    for (size_t k = 0; k < dim; k++)
    {
        // VectorRd ek = VectorRd::Zero();
        // ek(k) = 1;
        for (size_t i = 0; i < n_local_cell_scalar_dofs; i++)
        {
            for (size_t iqn = 0; iqn < n_quadT; iqn++)
            {
                div_vec_phiT_quadT[dim * i + k][iqn] = dphiT_quadT[i][iqn](k);
            }
        }
    }

    Eigen::MatrixXd DT_RHS = Eigen::MatrixXd::Zero(n_local_cell_scalar_dofs, n_local_vector_dofs);

    DT_RHS.topLeftCorner(n_local_cell_scalar_dofs, n_local_cell_vector_dofs) = compute_gram_matrix(phiT_quadT, div_vec_phiT_quadT, quadT, n_local_cell_scalar_dofs, n_local_cell_vector_dofs);

    for (size_t iTF = 0; iTF < n_local_faces; iTF++)
    {
        const size_t offset_F = n_local_cell_vector_dofs + iTF * n_local_face_vector_dofs;
        const VectorRd &nTF = cell->face_normal(iTF);

        QuadratureRule quadF = elquad.get_quadF(iTF);
        size_t n_quadF = quadF.size();

        BasisQuad<double> phiF_quadF = elquad.get_phiF_quadF(iTF);
        BasisQuad<double> phiT_quadF = elquad.get_phiT_quadF(iTF);

        BasisQuad<double> vec_phiF_quadF_dot_nTF(boost::extents[n_local_face_vector_dofs][n_quadF]);
        BasisQuad<double> vec_phiT_quadF_dot_nTF(boost::extents[n_local_cell_vector_dofs][n_quadF]);

        for (size_t k = 0; k < dim; k++)
        {
            VectorRd eT = VectorRd::Zero();
            eT(k) = 1;

            VectorRd eF = VectorRd::Zero();

            // Sort face unknowns into perpendicular, then paralell components
            if (k == 0)
            {
                eF = cell->face(iTF)->normal();
            }
            else if (k == 1)
            {
                eF = cell->face(iTF)->edge(0)->tangent();
            }
            else if (k == 2)
            {
                // eF = mesh_ptr->cell(iT)->face(iTF)->edge(1)->tangent();
                eF = (cell->face(iTF)->normal()).cross(cell->face(iTF)->edge(0)->tangent());
            }
            else
            {
                assert(false); // should never reach here
            }

            for (size_t iqn = 0; iqn < n_quadF; iqn++)
            {
                for (size_t i = 0; i < n_local_face_scalar_dofs; i++)
                {
                    vec_phiF_quadF_dot_nTF[dim * i + k][iqn] = phiF_quadF[i][iqn] * eF.dot(nTF);
                }
                for (size_t i = 0; i < n_local_cell_scalar_dofs; i++)
                {
                    vec_phiT_quadF_dot_nTF[dim * i + k][iqn] = phiT_quadF[i][iqn] * eT.dot(nTF);
                }
            }
        }

        DT_RHS.topLeftCorner(n_local_cell_scalar_dofs, n_local_cell_vector_dofs) -= compute_gram_matrix(phiT_quadF, vec_phiT_quadF_dot_nTF, quadF, n_local_cell_scalar_dofs, n_local_cell_vector_dofs);

        DT_RHS.block(0, offset_F, n_local_cell_scalar_dofs, n_local_face_vector_dofs) = compute_gram_matrix(phiT_quadF, vec_phiF_quadF_dot_nTF, quadF, n_local_cell_scalar_dofs, n_local_face_vector_dofs);
    }

    return -DT_RHS;
}

Eigen::MatrixXd MHDModel::gradient_term(Cell *cell, ElementQuad &elquad)
{
    const size_t n_local_faces = cell->n_faces();
    const size_t n_local_vector_dofs = n_local_cell_vector_dofs + n_local_faces * n_local_face_vector_dofs;

    QuadratureRule quadT = elquad.get_quadT();

    BasisQuad<double> phiT_quadT = elquad.get_phiT_quadT();
    BasisQuad<VectorRd> dphiT_quadT = elquad.get_dphiT_quadT();

    Eigen::MatrixXd RT_RHS = Eigen::MatrixXd::Zero(n_local_highorder_vector_dofs, n_local_vector_dofs);

    Eigen::MatrixXd MTT = compute_gram_matrix(phiT_quadT, phiT_quadT, quadT, n_local_cell_scalar_dofs, n_local_highorder_scalar_dofs, "sym");
    Eigen::MatrixXd StiffT = compute_gram_matrix(dphiT_quadT, dphiT_quadT, quadT, "sym");

    Eigen::MatrixXd vec_MTT = Eigen::MatrixXd::Zero(n_local_cell_vector_dofs, n_local_highorder_vector_dofs);

    Eigen::MatrixXd vec_StiffT = Eigen::MatrixXd::Zero(n_local_highorder_vector_dofs, n_local_highorder_vector_dofs);

    for (size_t k = 0; k < dim; k++)
    {
        for (size_t i = 0; i < n_local_highorder_scalar_dofs; i++)
        {
            for (size_t j = 0; j < n_local_highorder_scalar_dofs; j++)
            {
                if (i < n_local_cell_scalar_dofs)
                {
                    vec_MTT(dim * i + k, dim * j + k) = MTT(i, j);
                }
                vec_StiffT(dim * i + k, dim * j + k) = StiffT(i, j);
            }
        }
    }

    Eigen::MatrixXd LT = (vec_MTT.topRightCorner(dim, n_local_highorder_vector_dofs)).transpose(); // need both components of the constant term
    Eigen::MatrixXd LTtLT = LT * (LT.transpose());

    double scalT = vec_StiffT.norm() / LTtLT.norm();

    RT_RHS.topLeftCorner(n_local_highorder_vector_dofs, n_local_cell_vector_dofs) = vec_StiffT.topLeftCorner(n_local_highorder_vector_dofs, n_local_cell_vector_dofs) + scalT * LTtLT.topLeftCorner(n_local_highorder_vector_dofs, n_local_cell_vector_dofs);

    for (size_t iTF = 0; iTF < n_local_faces; iTF++)
    {
        const size_t offset_F = n_local_cell_vector_dofs + iTF * n_local_face_vector_dofs;
        const VectorRd &nTF = cell->face_normal(iTF);

        QuadratureRule quadF = elquad.get_quadF(iTF);
        size_t n_quadF = quadF.size();

        BasisQuad<double> phiF_quadF = elquad.get_phiF_quadF(iTF);
        BasisQuad<double> phiT_quadF = elquad.get_phiT_quadF(iTF);
        BasisQuad<VectorRd> dphiT_quadF = elquad.get_dphiT_quadF(iTF);

        BasisQuad<VectorRd> vec_phiF_quadF(boost::extents[n_local_face_vector_dofs][n_quadF]);
        BasisQuad<VectorRd> vec_phiT_quadF(boost::extents[n_local_cell_vector_dofs][n_quadF]);

        BasisQuad<VectorRd> nTF_doT_grad_vec_phiT_quadF(boost::extents[n_local_highorder_vector_dofs][n_quadF]);

        for (size_t k = 0; k < dim; k++)
        {
            VectorRd eT = VectorRd::Zero();
            eT(k) = 1;

            VectorRd eF = VectorRd::Zero();

            // Sort face unknowns into perpendicular, then paralell components
            if (k == 0)
            {
                eF = cell->face(iTF)->normal();
            }
            else if (k == 1)
            {
                eF = cell->face(iTF)->edge(0)->tangent();
            }
            else if (k == 2)
            {
                // eF = mesh_ptr->cell(iT)->face(iTF)->edge(1)->tangent();
                eF = (cell->face(iTF)->normal()).cross(cell->face(iTF)->edge(0)->tangent());
            }
            else
            {
                assert(false); // should never reach here
            }

            for (size_t iqn = 0; iqn < n_quadF; iqn++)
            {
                for (size_t i = 0; i < n_local_face_scalar_dofs; i++)
                {
                    vec_phiF_quadF[dim * i + k][iqn] = phiF_quadF[i][iqn] * eF;
                }
                for (size_t i = 0; i < n_local_highorder_scalar_dofs; i++)
                {
                    if (i < n_local_cell_scalar_dofs)
                    {
                        vec_phiT_quadF[dim * i + k][iqn] = phiT_quadF[i][iqn] * eT;
                    }
                    nTF_doT_grad_vec_phiT_quadF[dim * i + k][iqn] = (dphiT_quadF[i][iqn].dot(nTF)) * eT;
                }
            }
        }

        RT_RHS.topLeftCorner(n_local_highorder_vector_dofs, n_local_cell_vector_dofs) -= compute_gram_matrix(nTF_doT_grad_vec_phiT_quadF, vec_phiT_quadF, quadF, n_local_highorder_vector_dofs, n_local_cell_vector_dofs);

        RT_RHS.block(0, offset_F, n_local_highorder_vector_dofs, n_local_face_vector_dofs) = compute_gram_matrix(nTF_doT_grad_vec_phiT_quadF, vec_phiF_quadF, quadF, n_local_highorder_vector_dofs, n_local_face_vector_dofs);
    }

    RT[cell->global_index()] = (vec_StiffT + scalT * LTtLT).ldlt().solve(RT_RHS);

    return RT[cell->global_index()].transpose() * vec_StiffT * RT[cell->global_index()];
}

Eigen::MatrixXd MHDModel::th_ijk_plus_th_ikj(Cell *cell, ElementQuad &elquad, size_t k)
{
    const size_t n_local_faces = cell->n_faces();
    const size_t n_local_vector_dofs = n_local_cell_vector_dofs + n_local_faces * n_local_face_vector_dofs;

    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(n_local_vector_dofs, n_local_vector_dofs);

    assert(k < 2 * n_local_vector_dofs);
    if (k >= n_local_vector_dofs)
    {
        k -= n_local_vector_dofs;
    }

    if (k < n_local_cell_vector_dofs)
    {
        QuadratureRule quadT = elquad.get_quadT();

        size_t k_row = (k % dim);
        size_t k_scalar = (k / dim); // integer division

        VectorRd eK = VectorRd::Zero();
        eK(k_row) = 1.0;

        BasisQuad<double> phiT_quadT = elquad.get_phiT_quadT();
        BasisQuad<VectorRd> dphiT_quadT = elquad.get_dphiT_quadT();

        BasisQuad<VectorRd> vec_phiT_quadT(boost::extents[n_local_cell_vector_dofs][quadT.size()]);
        BasisQuad<VectorRd> phiT_quadT_dot_grad_basis_k(boost::extents[n_local_cell_vector_dofs][quadT.size()]);
        BasisQuad<VectorRd> basis_k_dot_grad_phiT_quadT(boost::extents[n_local_cell_vector_dofs][quadT.size()]);
        BasisQuad<MatrixRd> grad_phiT_quadT(boost::extents[n_local_cell_vector_dofs][quadT.size()]);
        BasisQuad<MatrixRd> basis_k_tensor_phiT_quadT(boost::extents[n_local_cell_vector_dofs][quadT.size()]);

        for (size_t d = 0; d < dim; ++d)
        {
            VectorRd eT = VectorRd::Zero(); 
            eT(d) = 1.0;
            for (size_t i = 0; i < n_local_cell_scalar_dofs; i++)
            {
                for (size_t iqn = 0; iqn < quadT.size(); iqn++)
                {
                    vec_phiT_quadT[dim * i + d][iqn] = phiT_quadT[i][iqn] * eT;
                    basis_k_dot_grad_phiT_quadT[dim * i + d][iqn] = phiT_quadT[k_scalar][iqn] * eK.dot(dphiT_quadT[i][iqn]) * eT;
                    phiT_quadT_dot_grad_basis_k[dim * i + d][iqn] = phiT_quadT[i][iqn] * eT.dot(dphiT_quadT[k_scalar][iqn]) * eK;
                    grad_phiT_quadT[dim * i + d][iqn] = MatrixRd::Zero();
                    grad_phiT_quadT[dim * i + d][iqn].row(d) = dphiT_quadT[i][iqn].transpose();
                    // grad_phiT_quadT[dim * i + d][iqn].col(d) = dphiT_quadT[i][iqn];
                    basis_k_tensor_phiT_quadT[dim * i + d][iqn] = MatrixRd::Zero();
                    basis_k_tensor_phiT_quadT[dim * i + d][iqn](k_row, d) = phiT_quadT[k_scalar][iqn] * phiT_quadT[i][iqn];
                    // basis_k_tensor_phiT_quadT[dim * i + d][iqn](d, k_row) =
                    // phiT_quadT[k_scalar][iqn] * phiT_quadT[i][iqn];
                }
            }
        }

        Eigen::MatrixXd uT_dot_grad_wT_dot_vT = compute_gram_matrix(vec_phiT_quadT, basis_k_dot_grad_phiT_quadT, quadT, n_local_cell_vector_dofs, n_local_cell_vector_dofs);

        T.topLeftCorner(n_local_cell_vector_dofs, n_local_cell_vector_dofs) = compute_gram_matrix(vec_phiT_quadT, phiT_quadT_dot_grad_basis_k, quadT, n_local_cell_vector_dofs, n_local_cell_vector_dofs) + uT_dot_grad_wT_dot_vT - compute_gram_matrix(grad_phiT_quadT, basis_k_tensor_phiT_quadT, quadT, n_local_cell_vector_dofs, n_local_cell_vector_dofs) - uT_dot_grad_wT_dot_vT.transpose();
        // T.topLeftCorner(n_local_cell_vector_dofs, n_local_cell_vector_dofs) =
        // compute_gram_matrix(vec_phiT_quadT, phiT_quadT_dot_grad_basis_k,
        // quadT).transpose() + uT_dot_grad_wT_dot_vT.transpose() -
        // compute_gram_matrix(basis_k_tensor_phiT_quadT, grad_phiT_quadT,
        // quadT).transpose() - uT_dot_grad_wT_dot_vT;

        for (size_t iTF = 0; iTF < n_local_faces; ++iTF)
        {
            const size_t offset_F = n_local_cell_vector_dofs + iTF * n_local_face_vector_dofs;
            VectorRd nTF = cell->face_normal(iTF);

            QuadratureRule quadF = elquad.get_quadF(iTF);

            BasisQuad<double> phiF_quadF = elquad.get_phiF_quadF(iTF);
            BasisQuad<double> phiT_quadF = elquad.get_phiT_quadF(iTF);

            BasisQuad<VectorRd> basis_k_dot_nTF_phiT_quadF(boost::extents[n_local_cell_vector_dofs][quadF.size()]);
            BasisQuad<double> phiT_quadF_dot_nTF(boost::extents[n_local_cell_vector_dofs][quadF.size()]);

            BasisQuad<double> basis_k_dot_phiF_quadF(boost::extents[n_local_face_vector_dofs][quadF.size()]);
            BasisQuad<VectorRd> vec_phiF_quadF(boost::extents[n_local_face_vector_dofs][quadF.size()]);

            for (size_t d = 0; d < dim; ++d)
            {
                // VectorRd eF = ((d == 0) ? cell->face(iTF)->normal() :
                // cell->face(iTF)->tangent());

                VectorRd eF = VectorRd::Zero();

                // Sort face unknowns into perpendicular, then paralell components
                if (d == 0)
                {
                    eF = cell->face(iTF)->normal();
                }
                else if (d == 1)
                {
                    eF = cell->face(iTF)->edge(0)->tangent();
                }
                else if (d == 2)
                {
                    // eF = mesh_ptr->cell(iT)->face(iTF)->edge(1)->tangent();
                    eF = (cell->face(iTF)->normal()).cross(cell->face(iTF)->edge(0)->tangent());
                }
                else
                {
                    assert(false); // should never reach here
                }

                VectorRd eT = VectorRd::Zero();
                eT(d) = 1.0;

                for (size_t iqn = 0; iqn < quadF.size(); iqn++)
                {
                    for (size_t i = 0; i < n_local_cell_scalar_dofs; i++)
                    {
                        basis_k_dot_nTF_phiT_quadF[dim * i + d][iqn] = phiT_quadF[k_scalar][iqn] * eK.dot(nTF) * phiT_quadF[i][iqn] * eT;
                        phiT_quadF_dot_nTF[dim * i + d][iqn] = phiT_quadF[i][iqn] * eT.dot(nTF);
                    }
                    for (size_t i = 0; i < n_local_face_scalar_dofs; i++)
                    {
                        basis_k_dot_phiF_quadF[dim * i + d][iqn] = phiT_quadF[k_scalar][iqn] * phiF_quadF[i][iqn] * eK.dot(eF);
                        vec_phiF_quadF[dim * i + d][iqn] = phiF_quadF[i][iqn] * eF;
                    }
                }
            }

            Eigen::MatrixXd uT_dot_nTF_vT_dot_wF = compute_gram_matrix(basis_k_dot_nTF_phiT_quadF, vec_phiF_quadF, quadF);

            T.block(0, offset_F, n_local_cell_vector_dofs, n_local_face_vector_dofs) = uT_dot_nTF_vT_dot_wF;
            T.block(offset_F, 0, n_local_face_vector_dofs, n_local_cell_vector_dofs) = -(compute_gram_matrix(basis_k_dot_phiF_quadF, phiT_quadF_dot_nTF, quadF) + uT_dot_nTF_vT_dot_wF.transpose());
        }

        return 0.5 * T;
    }

    // const size_t iT = cell->global_index();
    // const size_t n_local_faces = cell->n_faces();
    // const size_t n_local_vector_dofs = n_local_cell_vector_dofs + n_local_faces
    // * n_local_face_vector_dofs;

    assert(k < n_local_vector_dofs);

    size_t iTF = (k - n_local_cell_vector_dofs) / n_local_face_vector_dofs; // integer division
    size_t k_basis = (k - n_local_cell_vector_dofs) % n_local_face_vector_dofs;
    size_t k_scalar = k_basis / dim; // integer division

    VectorRd eF = VectorRd::Zero();

    // Sort face unknowns into perpendicular, then paralell components
    if (k_basis % dim == 0)
    {
        eF = cell->face(iTF)->normal();
    }
    else if (k_basis % dim == 1)
    {
        eF = cell->face(iTF)->edge(0)->tangent();
    }
    else if (k_basis % dim == 2)
    {
        // eF = mesh_ptr->cell(iT)->face(iTF)->edge(1)->tangent();
        eF = (cell->face(iTF)->normal()).cross(cell->face(iTF)->edge(0)->tangent());
    }
    else
    {
        assert(false); // should never reach here
    }

    VectorRd nTF = cell->face_normal(iTF);

    QuadratureRule quadF = elquad.get_quadF(iTF);

    BasisQuad<double> phiT_quadF = elquad.get_phiT_quadF(iTF);
    BasisQuad<double> phiF_quadF = elquad.get_phiF_quadF(iTF);

    BasisQuad<double> phiT_quadF_dot_nTF(boost::extents[n_local_cell_vector_dofs][quadF.size()]);
    BasisQuad<double> phiT_quadF_dot_basis_k(boost::extents[n_local_cell_vector_dofs][quadF.size()]);

    for (size_t d = 0; d < dim; ++d)
    {
        for (size_t i = 0; i < n_local_cell_scalar_dofs; i++)
        {
            for (size_t iqn = 0; iqn < quadF.size(); iqn++)
            {
                phiT_quadF_dot_nTF[dim * i + d][iqn] = phiT_quadF[i][iqn] * nTF(d); // eT.dot(nTF)
                phiT_quadF_dot_basis_k[dim * i + d][iqn] = phiT_quadF[i][iqn] * phiF_quadF[k_scalar][iqn] * eF(d); // eT.dot(eF)
            }
        }
    }

    // Eigen::MatrixXd T = Eigen::MatrixXd::Zero(n_local_vector_dofs,
    // n_local_vector_dofs);
    T.topLeftCorner(n_local_cell_vector_dofs, n_local_cell_vector_dofs) = 0.5 * compute_gram_matrix(phiT_quadF_dot_basis_k, phiT_quadF_dot_nTF, quadF, n_local_cell_vector_dofs, n_local_cell_vector_dofs);
    return T;

    // return 0.5 * compute_gram_matrix(phiT_quadF_dot_basis_k,
    // phiT_quadF_dot_nTF, quadF, n_local_cell_vector_dofs,
    // n_local_cell_vector_dofs);
}

Eigen::MatrixXd MHDModel::th_ijk_minus_th_ikj(Cell *cell, ElementQuad &elquad, size_t k)
{
    const size_t n_local_faces = cell->n_faces();
    const size_t n_local_vector_dofs = n_local_cell_vector_dofs + n_local_faces * n_local_face_vector_dofs;

    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(n_local_vector_dofs, n_local_vector_dofs);

    assert(k < 2 * n_local_vector_dofs);
    if (k >= n_local_vector_dofs)
    {
        k -= n_local_vector_dofs;
    }

    if (k < n_local_cell_vector_dofs)
    {
        QuadratureRule quadT = elquad.get_quadT();

        size_t k_row = (k % dim);
        size_t k_scalar = (k / dim); // integer division

        VectorRd eK = VectorRd::Zero();
        eK(k_row) = 1.0;

        BasisQuad<double> phiT_quadT = elquad.get_phiT_quadT();
        BasisQuad<VectorRd> dphiT_quadT = elquad.get_dphiT_quadT();

        BasisQuad<VectorRd> vec_phiT_quadT( boost::extents[n_local_cell_vector_dofs][quadT.size()]);
        BasisQuad<VectorRd> phiT_quadT_dot_grad_basis_k( boost::extents[n_local_cell_vector_dofs][quadT.size()]);
        BasisQuad<VectorRd> basis_k_dot_grad_phiT_quadT( boost::extents[n_local_cell_vector_dofs][quadT.size()]);
        BasisQuad<MatrixRd> grad_phiT_quadT( boost::extents[n_local_cell_vector_dofs][quadT.size()]);
        BasisQuad<MatrixRd> basis_k_tensor_phiT_quadT( boost::extents[n_local_cell_vector_dofs][quadT.size()]);

        for (size_t d = 0; d < dim; ++d)
        {
            VectorRd eT = VectorRd::Zero();
            eT(d) = 1.0;
            for (size_t i = 0; i < n_local_cell_scalar_dofs; i++)
            {
                for (size_t iqn = 0; iqn < quadT.size(); iqn++)
                {
                    vec_phiT_quadT[dim * i + d][iqn] = phiT_quadT[i][iqn] * eT;
                    basis_k_dot_grad_phiT_quadT[dim * i + d][iqn] = phiT_quadT[k_scalar][iqn] * eK.dot(dphiT_quadT[i][iqn]) * eT;
                    phiT_quadT_dot_grad_basis_k[dim * i + d][iqn] = phiT_quadT[i][iqn] * eT.dot(dphiT_quadT[k_scalar][iqn]) * eK;
                    grad_phiT_quadT[dim * i + d][iqn] = MatrixRd::Zero();
                    grad_phiT_quadT[dim * i + d][iqn].row(d) = dphiT_quadT[i][iqn].transpose();
                    basis_k_tensor_phiT_quadT[dim * i + d][iqn] = MatrixRd::Zero();
                    basis_k_tensor_phiT_quadT[dim * i + d][iqn](k_row, d) = phiT_quadT[k_scalar][iqn] * phiT_quadT[i][iqn];
                }
            }
        }

        Eigen::MatrixXd uT_dot_grad_wT_dot_vT = compute_gram_matrix(vec_phiT_quadT, basis_k_dot_grad_phiT_quadT, quadT);

        T.topLeftCorner(n_local_cell_vector_dofs, n_local_cell_vector_dofs) = -compute_gram_matrix(vec_phiT_quadT, phiT_quadT_dot_grad_basis_k, quadT) + uT_dot_grad_wT_dot_vT + compute_gram_matrix(grad_phiT_quadT, basis_k_tensor_phiT_quadT, quadT) - uT_dot_grad_wT_dot_vT.transpose();

        for (size_t iTF = 0; iTF < n_local_faces; ++iTF)
        {
            const size_t offset_F = n_local_cell_vector_dofs + iTF * n_local_face_vector_dofs;
            VectorRd nTF = cell->face_normal(iTF);

            QuadratureRule quadF = elquad.get_quadF(iTF);

            BasisQuad<double> phiF_quadF = elquad.get_phiF_quadF(iTF);
            BasisQuad<double> phiT_quadF = elquad.get_phiT_quadF(iTF);

            BasisQuad<VectorRd> basis_k_dot_nTF_phiT_quadF( boost::extents[n_local_cell_vector_dofs][quadF.size()]);
            BasisQuad<double> phiT_quadF_dot_nTF( boost::extents[n_local_cell_vector_dofs][quadF.size()]);

            BasisQuad<double> basis_k_dot_phiF_quadF( boost::extents[n_local_face_vector_dofs][quadF.size()]);
            BasisQuad<VectorRd> vec_phiF_quadF( boost::extents[n_local_face_vector_dofs][quadF.size()]);

            for (size_t d = 0; d < dim; ++d)
            {
                VectorRd eF = VectorRd::Zero();

                // Sort face unknowns into perpendicular, then paralell components
                if (d == 0)
                {
                    eF = cell->face(iTF)->normal();
                }
                else if (d == 1)
                {
                    eF = cell->face(iTF)->edge(0)->tangent();
                }
                else if (d == 2)
                {
                    // eF = mesh_ptr->cell(iT)->face(iTF)->edge(1)->tangent();
                    eF = (cell->face(iTF)->normal()).cross(cell->face(iTF)->edge(0)->tangent());
                }
                else
                {
                    assert(false); // should never reach here
                }

                VectorRd eT = VectorRd::Zero();
                eT(d) = 1.0;

                for (size_t iqn = 0; iqn < quadF.size(); iqn++)
                {
                    for (size_t i = 0; i < n_local_cell_scalar_dofs; i++)
                    {
                        basis_k_dot_nTF_phiT_quadF[dim * i + d][iqn] = phiT_quadF[k_scalar][iqn] * eK.dot(nTF) * phiT_quadF[i][iqn] * eT;
                        phiT_quadF_dot_nTF[dim * i + d][iqn] = phiT_quadF[i][iqn] * eT.dot(nTF);
                    }
                    for (size_t i = 0; i < n_local_face_scalar_dofs; i++)
                    {
                        basis_k_dot_phiF_quadF[dim * i + d][iqn] = phiT_quadF[k_scalar][iqn] * phiF_quadF[i][iqn] * eK.dot(eF);
                        vec_phiF_quadF[dim * i + d][iqn] = phiF_quadF[i][iqn] * eF;
                    }
                }
            }

            Eigen::MatrixXd uT_dot_nTF_vT_dot_wF = compute_gram_matrix( basis_k_dot_nTF_phiT_quadF, vec_phiF_quadF, quadF);

            T.block(0, offset_F, n_local_cell_vector_dofs, n_local_face_vector_dofs) = uT_dot_nTF_vT_dot_wF;
            T.block(offset_F, 0, n_local_face_vector_dofs, n_local_cell_vector_dofs) = -uT_dot_nTF_vT_dot_wF.transpose() + compute_gram_matrix(basis_k_dot_phiF_quadF, phiT_quadF_dot_nTF, quadF);
        }

        return 0.5 * T;
    }

    assert(k < n_local_vector_dofs);

    size_t iTF = (k - n_local_cell_vector_dofs) / n_local_face_vector_dofs; // integer division
    size_t k_basis = (k - n_local_cell_vector_dofs) % n_local_face_vector_dofs;
    size_t k_scalar = k_basis / dim; // integer division

    VectorRd eF = VectorRd::Zero();

    // Sort face unknowns into perpendicular, then paralell components
    if (k_basis % dim == 0)
    {
        eF = cell->face(iTF)->normal();
    }
    else if (k_basis % dim == 1)
    {
        eF = cell->face(iTF)->edge(0)->tangent();
    }
    else if (k_basis % dim == 2)
    {
        // eF = mesh_ptr->cell(iT)->face(iTF)->edge(1)->tangent();
        eF = (cell->face(iTF)->normal()).cross(cell->face(iTF)->edge(0)->tangent());
    }
    else
    {
        assert(false); // should never reach here
    }
    VectorRd nTF = cell->face_normal(iTF);

    QuadratureRule quadF = elquad.get_quadF(iTF);

    BasisQuad<double> phiT_quadF = elquad.get_phiT_quadF(iTF);
    BasisQuad<double> phiF_quadF = elquad.get_phiF_quadF(iTF);

    BasisQuad<double> phiT_quadF_dot_nTF( boost::extents[n_local_cell_vector_dofs][quadF.size()]);
    BasisQuad<double> phiT_quadF_dot_basis_k( boost::extents[n_local_cell_vector_dofs][quadF.size()]);

    for (size_t d = 0; d < dim; ++d)
    {
        for (size_t i = 0; i < n_local_cell_scalar_dofs; i++)
        {
            for (size_t iqn = 0; iqn < quadF.size(); iqn++)
            {
                phiT_quadF_dot_nTF[dim * i + d][iqn] = phiT_quadF[i][iqn] * nTF(d); // eT.dot(nTF)
                phiT_quadF_dot_basis_k[dim * i + d][iqn] = phiT_quadF[i][iqn] * phiF_quadF[k_scalar][iqn] * eF(d); // eT.dot(eF)
            }
        }
    }

    // Eigen::MatrixXd T = Eigen::MatrixXd::Zero(n_local_vector_dofs,
    // n_local_vector_dofs);
    T.topLeftCorner(n_local_cell_vector_dofs, n_local_cell_vector_dofs) = -0.5 * compute_gram_matrix(phiT_quadF_dot_basis_k, phiT_quadF_dot_nTF, quadF, n_local_cell_vector_dofs, n_local_cell_vector_dofs);
    return T;
}

Eigen::MatrixXd MHDModel::convective_term(Cell *cell, std::vector<Eigen::MatrixXd> tplus_k, std::vector<Eigen::MatrixXd> tminus_k, Eigen::VectorXd u)
{
    const size_t n_local_faces = cell->n_faces();
    const size_t n_local_vector_dofs = n_local_cell_vector_dofs + n_local_faces * n_local_face_vector_dofs;
    const size_t n_local_dofs = 2 * n_local_vector_dofs;
    assert((size_t)u.size() == n_local_dofs);

    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(n_local_dofs, n_local_dofs);

    for (size_t k = 0; k < n_local_vector_dofs; ++k)
    {
        T.topLeftCorner(n_local_vector_dofs, n_local_vector_dofs) += (u(k) + u(k + n_local_vector_dofs)) * tplus_k[k];
        T.topRightCorner(n_local_vector_dofs, n_local_vector_dofs) += (u(k) + u(k + n_local_vector_dofs)) * -tplus_k[k];
        T.bottomRightCorner(n_local_vector_dofs, n_local_vector_dofs) += (u(k) + u(k + n_local_vector_dofs)) * tminus_k[k];
        T.bottomLeftCorner(n_local_vector_dofs, n_local_vector_dofs) += (u(k) + u(k + n_local_vector_dofs)) * -tminus_k[k];
    }

    return T;
}

Eigen::MatrixXd MHDModel::stabilisation_term(Cell *cell, ElementQuad &elquad)
{
    const size_t n_local_faces = cell->n_faces();
    const size_t n_local_bdry_vector_dofs = n_local_faces * n_local_face_vector_dofs;

    QuadratureRule quadT = elquad.get_quadT();

    BasisQuad<double> phiT_quadT = elquad.get_phiT_quadT();
    BasisQuad<VectorRd> dphiT_quadT = elquad.get_dphiT_quadT();

    Eigen::MatrixXd MTT = compute_gram_matrix(phiT_quadT, phiT_quadT, quadT, n_local_cell_scalar_dofs, n_local_highorder_scalar_dofs, "sym");
    Eigen::MatrixXd StiffT = compute_gram_matrix(dphiT_quadT, dphiT_quadT, quadT, n_local_cell_scalar_dofs, n_local_cell_scalar_dofs, "sym");

    Eigen::MatrixXd vec_MTT = Eigen::MatrixXd::Zero(n_local_cell_vector_dofs, n_local_highorder_vector_dofs);
    Eigen::MatrixXd vec_StiffT = Eigen::MatrixXd::Zero(n_local_cell_vector_dofs, n_local_cell_vector_dofs);

    for (size_t k = 0; k < dim; k++)
    {
        for (size_t i = 0; i < n_local_cell_scalar_dofs; i++)
        {
            for (size_t j = 0; j < n_local_highorder_scalar_dofs; j++)
            {
                if (j < n_local_cell_scalar_dofs)
                {
                    vec_StiffT(dim * i + k, dim * j + k) = StiffT(i, j);
                }
                vec_MTT(dim * i + k, dim * j + k) = MTT(i, j);
            }
        }
    }

    Eigen::MatrixXd diffT = (vec_MTT.topLeftCorner(n_local_cell_vector_dofs, n_local_cell_vector_dofs).ldlt()).solve(vec_MTT) * RT[cell->global_index()];
    diffT.topLeftCorner(n_local_cell_vector_dofs, n_local_cell_vector_dofs) -= Eigen::MatrixXd::Identity(n_local_cell_vector_dofs, n_local_cell_vector_dofs);

    std::vector<Eigen::Triplet<double>> triplets_M_BDRY_BDRY;
    std::vector<Eigen::Triplet<double>> triplets_M_BDRY_T;
    std::vector<Eigen::Triplet<double>> triplets_SCALE_MAT;

    Eigen::SparseMatrix<double> M_BDRY_BDRY(n_local_bdry_vector_dofs, n_local_bdry_vector_dofs);
    Eigen::SparseMatrix<double> M_BDRY_T(n_local_bdry_vector_dofs, n_local_highorder_vector_dofs);
    Eigen::SparseMatrix<double> BDRY_SCALE_MAT(n_local_bdry_vector_dofs, n_local_bdry_vector_dofs);

    for (size_t iTF = 0; iTF < n_local_faces; iTF++)
    {
        const size_t offset_F = iTF * n_local_face_vector_dofs;

        QuadratureRule quadF = elquad.get_quadF(iTF);
        size_t n_quadF = quadF.size();

        BasisQuad<double> phiF_quadF = elquad.get_phiF_quadF(iTF);
        BasisQuad<double> phiT_quadF = elquad.get_phiT_quadF(iTF);

        BasisQuad<VectorRd> vec_phiF_quadF( boost::extents[n_local_face_vector_dofs][n_quadF]);
        BasisQuad<VectorRd> vec_phiT_quadF( boost::extents[n_local_highorder_vector_dofs][n_quadF]);

        for (size_t d = 0; d < dim; d++)
        {
            VectorRd eT = VectorRd::Zero();
            eT(d) = 1;

            VectorRd eF = VectorRd::Zero();

            // Sort face unknowns into perpendicular, then paralell components
            if (d == 0)
            {
                eF = cell->face(iTF)->normal();
            }
            else if (d == 1)
            {
                eF = cell->face(iTF)->edge(0)->tangent();
            }
            else if (d == 2)
            {
                // eF = mesh_ptr->cell(iT)->face(iTF)->edge(1)->tangent();
                eF = (cell->face(iTF)->normal()).cross(cell->face(iTF)->edge(0)->tangent());
            }
            else
            {
                assert(false); // should never reach here
            }

            for (size_t iqn = 0; iqn < n_quadF; iqn++)
            {
                for (size_t i = 0; i < n_local_face_scalar_dofs; i++)
                {
                    vec_phiF_quadF[dim * i + d][iqn] = phiF_quadF[i][iqn] * eF;
                }
                for (size_t i = 0; i < n_local_highorder_scalar_dofs; i++)
                {
                    vec_phiT_quadF[dim * i + d][iqn] = phiT_quadF[i][iqn] * eT;
                }
            }
        }

        Eigen::MatrixXd vec_MFF = compute_gram_matrix(vec_phiF_quadF, vec_phiF_quadF, quadF, "sym");
        Eigen::MatrixXd vec_MFT = compute_gram_matrix(vec_phiF_quadF, vec_phiT_quadF, quadF, n_local_face_vector_dofs, n_local_highorder_vector_dofs);

        Eigen::MatrixXd LOCAL_SCALING = (1.0 / (cell->diam())) * Eigen::MatrixXd::Identity(n_local_face_vector_dofs, n_local_face_vector_dofs);

        for (size_t i = 0; i < n_local_face_vector_dofs; ++i)
        {
            for (size_t j = 0; j < n_local_face_vector_dofs; ++j)
            {
                triplets_M_BDRY_BDRY.emplace_back(offset_F + i, offset_F + j, vec_MFF(i, j));
                triplets_SCALE_MAT.emplace_back(offset_F + i, offset_F + j, LOCAL_SCALING(i, j));
            }
            for (size_t k = 0; k < n_local_highorder_vector_dofs; ++k)
            {
                triplets_M_BDRY_T.emplace_back(offset_F + i, k, vec_MFT(i, k));
            }
        }
    }

    M_BDRY_BDRY.setFromTriplets(std::begin(triplets_M_BDRY_BDRY), std::end(triplets_M_BDRY_BDRY));
    M_BDRY_T.setFromTriplets(std::begin(triplets_M_BDRY_T), std::end(triplets_M_BDRY_T));
    BDRY_SCALE_MAT.setFromTriplets(std::begin(triplets_SCALE_MAT), std::end(triplets_SCALE_MAT));

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(M_BDRY_BDRY);

    Eigen::MatrixXd inv_M_BDRY_BDRY_M_BDRY_T = solver.solve(M_BDRY_T);

    // Eigen::MatrixXd diffBDRY = solver.solve(M_BDRY_T) * RT[cell->global_index()];
    Eigen::MatrixXd diffBDRY = inv_M_BDRY_BDRY_M_BDRY_T * RT[cell->global_index()];
    diffBDRY.bottomRightCorner(n_local_bdry_vector_dofs, n_local_bdry_vector_dofs) -= Eigen::MatrixXd::Identity(n_local_bdry_vector_dofs, n_local_bdry_vector_dofs);

    return diffT.transpose() * vec_StiffT * diffT + diffBDRY.transpose() * BDRY_SCALE_MAT * M_BDRY_BDRY * diffBDRY;
    // return diffBDRY.transpose() * BDRY_SCALE_MAT * M_BDRY_BDRY * diffBDRY;

    // Eigen::MatrixXd diffT_on_bdry = inv_M_BDRY_BDRY_M_BDRY_T.topLeftCorner(n_local_bdry_vector_dofs, n_local_cell_vector_dofs)  * diffT;

    // return (diffT_on_bdry - diffBDRY).transpose() * BDRY_SCALE_MAT * M_BDRY_BDRY * (diffT_on_bdry - diffBDRY);
}

Eigen::VectorXd MHDModel::source_term(FType<VectorRd> f_vec, Cell *cell, ElementQuad &elquad)
{
    const size_t n_local_faces = cell->n_faces();
    const size_t n_local_vector_dofs = n_local_cell_vector_dofs + n_local_faces * n_local_face_vector_dofs;

    // get quadrature rule and basis function on quadrature nodes
    QuadratureRule quadT = elquad.get_quadT();
    BasisQuad<double> phiT_quadT = elquad.get_phiT_quadT();
    BasisQuad<VectorRd> vec_phiT_quadT(boost::extents[n_local_highorder_vector_dofs][quadT.size()]);

    Eigen::VectorXd source_vec = Eigen::VectorXd::Zero(n_local_vector_dofs);

    for (size_t k = 0; k < dim; k++)
    {
        VectorRd eT = VectorRd::Zero();
        eT(k) = 1;

        for (size_t iqn = 0; iqn < quadT.size(); iqn++)
        {
            for (size_t i = 0; i < n_local_cell_scalar_dofs; i++)
            {
                vec_phiT_quadT[dim * i + k][iqn] = phiT_quadT[i][iqn] * eT;
            }
        }
    }

    // Integrate source term with basis functions
    source_vec.head(n_local_cell_vector_dofs) = integrate(f_vec, vec_phiT_quadT, quadT, n_local_cell_vector_dofs);
    return source_vec;
}

SolutionVector MHDModel::solve_with_static_cond(FType<VectorRd> f_source, FType<VectorRd> g_source, double visc, double diff, double tol, bool threading)
{

    std::vector<Eigen::VectorXd> LTf(n_cells);
    std::vector<Eigen::VectorXd> LTg(n_cells);
    std::vector<Eigen::MatrixXd> DT(n_cells);
    // std::vector<std::vector<Eigen::MatrixXd>> TTk(n_cells);
    std::vector<std::vector<Eigen::MatrixXd>> Tplus_k(n_cells);
    std::vector<std::vector<Eigen::MatrixXd>> Tminus_k(n_cells);

    std::cout << "     Assembling local matrices\n";
    boost::timer::cpu_timer assembly_timer;
    assembly_timer.start();

    // Construct the local matrices using multithreading if use_threads is true
    std::function<void(size_t, size_t)> construct_all_local_contributions = [&](size_t start, size_t end) -> void
    {
        for (size_t iT = start; iT < end; ++iT)
        {
            // ElementQuad elquad(m_hho, iT, m_K + 1 + m_K + 1 + 2, m_K + m_K + 1);

            ElementQuad elquad(m_hho, iT, std::max(3 * m_L, m_L + m_K + 1) + 1, std::max(3 * m_K, m_K + m_K + 1) + 1);
            // ElementQuad elquad(m_hho, iT, 15, 15);
            Cell *cell = mesh_ptr->cell(iT);

            AT[iT] = gradient_term(cell, elquad);
            AT[iT] = AT[iT] + stabilisation_term(cell, elquad);
            LTf[iT] = (source_term(f_source, cell, elquad).head(n_local_cell_vector_dofs));
            LTg[iT] = (source_term(g_source, cell, elquad).head(n_local_cell_vector_dofs));
            DT[iT] = (divergence_term(cell, elquad));
            size_t n_local_vector_dofs = n_local_cell_vector_dofs + cell->n_faces() * n_local_face_vector_dofs;
            // for (size_t k = 0; k < 2 * n_local_vector_dofs; ++k)
            for (size_t k = 0; k < n_local_vector_dofs; ++k)
            {
                // TTk[iT].push_back(T_ijk_plus_T_ikj(cell, elquad, k));
                Tplus_k[iT].push_back(th_ijk_plus_th_ikj(cell, elquad, k));
                Tminus_k[iT].push_back(th_ijk_minus_th_ikj(cell, elquad, k));
            }
        }
    };

    // Running the local constructions in parallel
    parallel_for(n_cells, construct_all_local_contributions, threading);

    double exact_rhs_norm = 0.0;
    for (size_t iT = 0; iT < n_cells; ++iT)
    {
        exact_rhs_norm += LTf[iT].squaredNorm() + LTg[iT].squaredNorm();
    }
    exact_rhs_norm = std::sqrt(exact_rhs_norm);

    assembly_timer.stop();
    std::cout << "     Assembly time = " << double(assembly_timer.elapsed().wall) * 1E-9 << "\n";

    size_t n_magnetic_dofs = n_total_cell_vector_dofs + n_total_face_vector_dofs + n_total_cell_scalar_dofs;
    size_t n_velocity_dofs = n_total_cell_vector_dofs + n_total_face_vector_dofs + n_total_cell_scalar_dofs;

    size_t n_total_dofs = n_magnetic_dofs + n_velocity_dofs;

    size_t n_velocity_unknowns = n_total_face_vector_dofs + n_cells; // face dofs plus one scalar dof per cell
    size_t n_magnetic_unknowns = n_total_face_vector_dofs + n_cells; // face dofs plus one scalar dof per cell

    size_t n_implicit_velocity_dofs = n_total_cell_vector_dofs + n_cells * (n_local_cell_scalar_dofs - 1);
    size_t n_implicit_magnetic_dofs = n_total_cell_vector_dofs + n_cells * (n_local_cell_scalar_dofs - 1);

    size_t n_unknowns = n_velocity_unknowns + n_magnetic_unknowns;
    size_t n_implicit_dofs = n_implicit_velocity_dofs + n_implicit_magnetic_dofs;

    SolutionVector SOL(Eigen::VectorXd::Zero(n_total_dofs), mesh_ptr, m_L, m_K);

    std::cout << "     Solving scheme using Newton method\n";

    boost::timer::cpu_timer solving_timer;
    solving_timer.start();

    double residual_error = 0.0;
    int i = 0;

    do
    {
        std::vector<Eigen::VectorXd> local_RHS(n_cells);

        std::vector<Eigen::MatrixXd> local_inv_LTll_LTlg(n_cells);
        std::vector<Eigen::MatrixXd> local_LTgl(n_cells);
        std::vector<Eigen::MatrixXd> local_LTgg(n_cells);

        Eigen::VectorXd inv_LTll_rTl = Eigen::VectorXd::Zero(n_implicit_dofs);

        // Construct the local matrices using multithreading if use_threads is true
        std::function<void(size_t, size_t)> build_local_mats = [&](size_t start, size_t end) -> void
        {
            for (size_t iT = start; iT < end; ++iT)
            {
                Cell *cell = mesh_ptr->cell(iT);
                size_t n_local_faces = cell->n_faces();
                size_t n_local_bdry_vector_dofs = n_local_faces * n_local_face_vector_dofs;
                size_t n_local_vector_dofs = n_local_cell_vector_dofs + n_local_bdry_vector_dofs;
                size_t n_local_unknowns = 2 * (n_local_vector_dofs + n_local_cell_scalar_dofs);
                size_t magnetic_offset = n_local_vector_dofs + n_local_cell_scalar_dofs;

                Eigen::VectorXd uT = SOL.restr(iT);
                Eigen::VectorXd vector_terms = Eigen::VectorXd::Zero(2 * n_local_vector_dofs);
                vector_terms.head(n_local_vector_dofs) = uT.head(n_local_vector_dofs);
                vector_terms.tail(n_local_vector_dofs) = uT.segment(magnetic_offset, n_local_vector_dofs);

                Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n_local_unknowns, n_local_unknowns);

                A.block(0, 0, n_local_vector_dofs, n_local_vector_dofs) = visc * AT[iT];
                A.block(0, n_local_vector_dofs, n_local_vector_dofs, n_local_cell_scalar_dofs) = DT[iT].transpose();
                A.block(n_local_vector_dofs, 0, n_local_cell_scalar_dofs, n_local_vector_dofs) = -DT[iT];

                A.block(magnetic_offset, magnetic_offset, n_local_vector_dofs, n_local_vector_dofs) = diff * AT[iT];
                A.block(magnetic_offset, magnetic_offset + n_local_vector_dofs, n_local_vector_dofs, n_local_cell_scalar_dofs) = DT[iT].transpose();
                A.block(magnetic_offset + n_local_vector_dofs, magnetic_offset, n_local_cell_scalar_dofs, n_local_vector_dofs) = -DT[iT];

                // Eigen::MatrixXd T = convective_term(cell, TTk[iT], vector_terms);
                Eigen::MatrixXd T = convective_term(cell, Tplus_k[iT], Tminus_k[iT], vector_terms);
                Eigen::MatrixXd BMAT = Eigen::MatrixXd::Zero(n_local_unknowns, n_local_unknowns);
                BMAT.block(0, 0, n_local_vector_dofs, n_local_vector_dofs) = T.topLeftCorner(n_local_vector_dofs, n_local_vector_dofs);
                BMAT.block(magnetic_offset, 0, n_local_vector_dofs, n_local_vector_dofs) = T.bottomLeftCorner(n_local_vector_dofs, n_local_vector_dofs);
                BMAT.block(0, magnetic_offset, n_local_vector_dofs, n_local_vector_dofs) = T.topRightCorner(n_local_vector_dofs, n_local_vector_dofs);
                BMAT.block(magnetic_offset, magnetic_offset, n_local_vector_dofs, n_local_vector_dofs) = T.bottomRightCorner(n_local_vector_dofs, n_local_vector_dofs);

                Eigen::MatrixXd local_DG = A + BMAT;
                Eigen::MatrixXd GMAT = A + 0.5 * BMAT;

                Eigen::VectorXd G = GMAT * uT;
                G.segment(0, n_local_cell_vector_dofs) -= LTf[iT];
                G.segment(magnetic_offset, n_local_cell_vector_dofs) -= LTg[iT];

                local_RHS[iT] = -G;

                // local_DG.block(0, 0, n_local_cell_vector_dofs, n_local_cell_vector_dofs) += std::pow(deltaN, -1) * vec_MTT[iT];
                // local_DG.block(magnetic_offset, magnetic_offset, n_local_cell_vector_dofs, n_local_cell_vector_dofs) += std::pow(deltaN, -1) * vec_MTT[iT];

                // local_DG.block(n_local_vector_dofs, n_local_vector_dofs, n_local_cell_scalar_dofs, n_local_cell_scalar_dofs) += std::pow(deltaN, -1) * MTT[iT];
                // local_DG.block(n_local_vector_dofs + magnetic_offset, n_local_vector_dofs + magnetic_offset, n_local_cell_scalar_dofs, n_local_cell_scalar_dofs) += std::pow(deltaN, -1) * MTT[iT];

                size_t n_local_condensed_velocity_dofs = n_local_bdry_vector_dofs + 1;
                size_t n_local_condensed_magnetic_dofs = n_local_bdry_vector_dofs + 1;

                size_t n_local_condensed_dofs = n_local_condensed_velocity_dofs + n_local_condensed_magnetic_dofs;

                size_t n_local_implicit_velocity_dofs = n_local_cell_vector_dofs + n_local_cell_scalar_dofs - 1;
                size_t n_local_implicit_magnetic_dofs = n_local_cell_vector_dofs + n_local_cell_scalar_dofs - 1;

                size_t n_local_implicit_dofs = n_local_implicit_velocity_dofs + n_local_implicit_magnetic_dofs;

                size_t local_magnetic_offset = n_local_vector_dofs + n_local_cell_scalar_dofs;
                size_t local_condensed_magnetic_offset = n_local_bdry_vector_dofs + 1;
                size_t local_implicit_magnetic_offset = n_local_cell_vector_dofs + n_local_cell_scalar_dofs - 1;

                Eigen::MatrixXd local_LTll = Eigen::MatrixXd::Zero(n_local_implicit_dofs, n_local_implicit_dofs);

                // velocity-velocity dofs
                local_LTll.block(0, 0, n_local_cell_vector_dofs, n_local_cell_vector_dofs) = local_DG.block(0, 0, n_local_cell_vector_dofs, n_local_cell_vector_dofs);
                local_LTll.block(0, n_local_cell_vector_dofs, n_local_cell_vector_dofs, n_local_cell_scalar_dofs - 1) = local_DG.block(0, n_local_vector_dofs + 1, n_local_cell_vector_dofs, n_local_cell_scalar_dofs - 1);
                local_LTll.block(n_local_cell_vector_dofs, 0, n_local_cell_scalar_dofs - 1, n_local_cell_vector_dofs) = local_DG.block(n_local_vector_dofs + 1, 0, n_local_cell_scalar_dofs - 1, n_local_cell_vector_dofs);

                // velocity pressure-pressure dofs (relaxation)
                // local_LTll.block(n_local_cell_vector_dofs, n_local_cell_vector_dofs, n_local_cell_scalar_dofs - 1, n_local_cell_scalar_dofs - 1) = local_DG.block(n_local_vector_dofs + 1, n_local_vector_dofs + 1, n_local_cell_scalar_dofs - 1, n_local_cell_scalar_dofs - 1);

                // velocity-magnetic dofs
                local_LTll.block(0, local_implicit_magnetic_offset, n_local_cell_vector_dofs, n_local_cell_vector_dofs) = local_DG.block(0, local_magnetic_offset, n_local_cell_vector_dofs, n_local_cell_vector_dofs);

                // magnetic-velocity dofs
                local_LTll.block(local_implicit_magnetic_offset, 0, n_local_cell_vector_dofs, n_local_cell_vector_dofs) = local_DG.block(local_magnetic_offset, 0, n_local_cell_vector_dofs, n_local_cell_vector_dofs);

                // magnetic-magnetic dofs
                local_LTll.block(local_implicit_magnetic_offset, local_implicit_magnetic_offset, n_local_cell_vector_dofs, n_local_cell_vector_dofs) = local_DG.block(local_magnetic_offset, local_magnetic_offset, n_local_cell_vector_dofs, n_local_cell_vector_dofs);
                local_LTll.block(local_implicit_magnetic_offset, local_implicit_magnetic_offset + n_local_cell_vector_dofs, n_local_cell_vector_dofs, n_local_cell_scalar_dofs - 1) = local_DG.block(local_magnetic_offset, local_magnetic_offset + n_local_vector_dofs + 1, n_local_cell_vector_dofs, n_local_cell_scalar_dofs - 1);
                local_LTll.block(local_implicit_magnetic_offset + n_local_cell_vector_dofs, local_implicit_magnetic_offset, n_local_cell_scalar_dofs - 1, n_local_cell_vector_dofs) = local_DG.block(local_magnetic_offset + n_local_vector_dofs + 1, local_magnetic_offset, n_local_cell_scalar_dofs - 1, n_local_cell_vector_dofs);

                // magnetic pressure-pressure dofs (relaxation)
                // local_LTll.block(local_implicit_magnetic_offset + n_local_cell_vector_dofs, local_implicit_magnetic_offset + n_local_cell_vector_dofs, n_local_cell_scalar_dofs - 1, n_local_cell_scalar_dofs - 1) = local_DG.block(local_magnetic_offset + n_local_vector_dofs + 1, local_magnetic_offset + n_local_vector_dofs + 1, n_local_cell_scalar_dofs - 1, n_local_cell_scalar_dofs - 1);

                Eigen::MatrixXd local_LTlg = Eigen::MatrixXd::Zero(n_local_implicit_dofs, n_local_condensed_dofs);

                // velocity-velocity dofs
                local_LTlg.block(0, 0, n_local_cell_vector_dofs, n_local_bdry_vector_dofs) = local_DG.block(0, n_local_cell_vector_dofs, n_local_cell_vector_dofs, n_local_bdry_vector_dofs);
                local_LTlg.block(n_local_cell_vector_dofs, 0, n_local_cell_scalar_dofs - 1, n_local_bdry_vector_dofs) = local_DG.block(n_local_vector_dofs + 1, n_local_cell_vector_dofs, n_local_cell_scalar_dofs - 1, n_local_bdry_vector_dofs);

                // velocity-magnetic dofs
                local_LTlg.block(0, local_condensed_magnetic_offset, n_local_cell_vector_dofs, n_local_bdry_vector_dofs) = local_DG.block(0, local_magnetic_offset + n_local_cell_vector_dofs, n_local_cell_vector_dofs, n_local_bdry_vector_dofs);
                // local_LTlg.block(0 + n_local_cell_vector_dofs, local_condensed_magnetic_offset, n_local_cell_scalar_dofs - 1, n_local_bdry_vector_dofs) = local_DG.block(0 + n_local_vector_dofs + 1, local_magnetic_offset + n_local_cell_vector_dofs, n_local_cell_scalar_dofs - 1, n_local_bdry_vector_dofs);

                // magnetic-velocity dofs
                local_LTlg.block(local_implicit_magnetic_offset, 0, n_local_cell_vector_dofs, n_local_bdry_vector_dofs) = local_DG.block(local_magnetic_offset, 0 + n_local_cell_vector_dofs, n_local_cell_vector_dofs, n_local_bdry_vector_dofs);
                // local_LTlg.block(local_implicit_magnetic_offset + n_local_cell_vector_dofs, 0, n_local_cell_scalar_dofs - 1, n_local_bdry_vector_dofs) = local_DG.block(local_magnetic_offset + n_local_vector_dofs + 1, 0 + n_local_cell_vector_dofs, n_local_cell_scalar_dofs - 1, n_local_bdry_vector_dofs);

                // magnetic-magnetic dofs
                local_LTlg.block(local_implicit_magnetic_offset, local_condensed_magnetic_offset, n_local_cell_vector_dofs, n_local_bdry_vector_dofs) = local_DG.block(local_magnetic_offset, local_magnetic_offset + n_local_cell_vector_dofs, n_local_cell_vector_dofs, n_local_bdry_vector_dofs);
                local_LTlg.block(local_implicit_magnetic_offset + n_local_cell_vector_dofs, local_condensed_magnetic_offset, n_local_cell_scalar_dofs - 1, n_local_bdry_vector_dofs) = local_DG.block(local_magnetic_offset + n_local_vector_dofs + 1, local_magnetic_offset + n_local_cell_vector_dofs, n_local_cell_scalar_dofs - 1, n_local_bdry_vector_dofs);

                local_LTgl[iT] = Eigen::MatrixXd::Zero(n_local_condensed_dofs, n_local_implicit_dofs);

                // velocity-velocity dofs
                local_LTgl[iT].block(0, 0, n_local_bdry_vector_dofs, n_local_cell_vector_dofs) = local_DG.block(n_local_cell_vector_dofs, 0, n_local_bdry_vector_dofs, n_local_cell_vector_dofs);
                local_LTgl[iT].block(0, n_local_cell_vector_dofs, n_local_bdry_vector_dofs, n_local_cell_scalar_dofs - 1) = local_DG.block(n_local_cell_vector_dofs, n_local_vector_dofs + 1, n_local_bdry_vector_dofs, n_local_cell_scalar_dofs - 1);

                // velocity-magnetic dofs
                local_LTgl[iT].block(0, local_implicit_magnetic_offset, n_local_bdry_vector_dofs, n_local_cell_vector_dofs) = local_DG.block(0 + n_local_cell_vector_dofs, local_magnetic_offset, n_local_bdry_vector_dofs, n_local_cell_vector_dofs);
                // local_LTgl.block(0, local_implicit_magnetic_offset + n_local_cell_vector_dofs, n_local_bdry_vector_dofs, n_local_cell_scalar_dofs - 1) = local_DG.block(0 + n_local_cell_vector_dofs, local_magnetic_offset + n_local_vector_dofs + 1, n_local_bdry_vector_dofs, n_local_cell_scalar_dofs - 1);

                // magnetic-velocity dofs
                local_LTgl[iT].block(local_condensed_magnetic_offset, 0, n_local_bdry_vector_dofs, n_local_cell_vector_dofs) = local_DG.block(local_magnetic_offset + n_local_cell_vector_dofs, 0, n_local_bdry_vector_dofs, n_local_cell_vector_dofs);
                // local_LTgl.block(local_condensed_magnetic_offset, 0 + n_local_cell_vector_dofs, n_local_bdry_vector_dofs, n_local_cell_scalar_dofs - 1) = local_DG.block(local_magnetic_offset + n_local_cell_vector_dofs, 0 + n_local_vector_dofs + 1, n_local_bdry_vector_dofs, n_local_cell_scalar_dofs - 1);

                // magnetic-magnetic dofs
                local_LTgl[iT].block(local_condensed_magnetic_offset, local_implicit_magnetic_offset, n_local_bdry_vector_dofs, n_local_cell_vector_dofs) = local_DG.block(local_magnetic_offset + n_local_cell_vector_dofs, local_magnetic_offset, n_local_bdry_vector_dofs, n_local_cell_vector_dofs);
                local_LTgl[iT].block(local_condensed_magnetic_offset, local_implicit_magnetic_offset + n_local_cell_vector_dofs, n_local_bdry_vector_dofs, n_local_cell_scalar_dofs - 1) = local_DG.block(local_magnetic_offset + n_local_cell_vector_dofs, local_magnetic_offset + n_local_vector_dofs + 1, n_local_bdry_vector_dofs, n_local_cell_scalar_dofs - 1);

                local_LTgg[iT] = Eigen::MatrixXd::Zero(n_local_condensed_dofs, n_local_condensed_dofs);

                // velocity-velocity dofs
                local_LTgg[iT].block(0, 0, n_local_bdry_vector_dofs, n_local_bdry_vector_dofs) = local_DG.block(n_local_cell_vector_dofs, n_local_cell_vector_dofs, n_local_bdry_vector_dofs, n_local_bdry_vector_dofs);
                local_LTgg[iT].block(0, n_local_bdry_vector_dofs, n_local_bdry_vector_dofs, 1) = local_DG.block(n_local_cell_vector_dofs, n_local_cell_vector_dofs + n_local_bdry_vector_dofs, n_local_bdry_vector_dofs, 1);
                local_LTgg[iT].block(n_local_bdry_vector_dofs, 0, 1, n_local_bdry_vector_dofs) = local_DG.block(n_local_cell_vector_dofs + n_local_bdry_vector_dofs, n_local_cell_vector_dofs, 1, n_local_bdry_vector_dofs);

                // velocity pressure-pressure dofs (relaxation)
                // local_LTgg[iT](n_local_bdry_vector_dofs, n_local_bdry_vector_dofs) = local_DG(n_local_cell_vector_dofs + n_local_bdry_vector_dofs, n_local_cell_vector_dofs + n_local_bdry_vector_dofs);

                // velocity-magnetic dofs
                local_LTgg[iT].block(0, local_condensed_magnetic_offset, n_local_bdry_vector_dofs, n_local_bdry_vector_dofs) = local_DG.block(0 + n_local_cell_vector_dofs, local_magnetic_offset + n_local_cell_vector_dofs, n_local_bdry_vector_dofs, n_local_bdry_vector_dofs);
                // local_LTgg[iT].block(0, local_condensed_magnetic_offset + n_local_bdry_vector_dofs, n_local_bdry_vector_dofs, 1) = local_DG.block(0 + n_local_cell_vector_dofs, local_magnetic_offset + n_local_cell_vector_dofs + n_local_bdry_vector_dofs, n_local_bdry_vector_dofs, 1);
                // local_LTgg[iT].block(0 + n_local_bdry_vector_dofs, local_condensed_magnetic_offset, 1, n_local_bdry_vector_dofs) = local_DG.block(0 + n_local_cell_vector_dofs + n_local_bdry_vector_dofs, local_magnetic_offset + n_local_cell_vector_dofs, 1, n_local_bdry_vector_dofs);

                // magnetic-velocity dofs
                local_LTgg[iT].block(local_condensed_magnetic_offset, 0, n_local_bdry_vector_dofs, n_local_bdry_vector_dofs) = local_DG.block(local_magnetic_offset + n_local_cell_vector_dofs, 0 + n_local_cell_vector_dofs, n_local_bdry_vector_dofs, n_local_bdry_vector_dofs);
                // local_LTgg[iT].block(local_condensed_magnetic_offset, 0 + n_local_bdry_vector_dofs, n_local_bdry_vector_dofs, 1) = local_DG.block(local_magnetic_offset + n_local_cell_vector_dofs, 0 + n_local_cell_vector_dofs + n_local_bdry_vector_dofs, n_local_bdry_vector_dofs, 1);
                // local_LTgg[iT].block(local_condensed_magnetic_offset + n_local_bdry_vector_dofs, 0, 1, n_local_bdry_vector_dofs) = local_DG.block(local_magnetic_offset + n_local_cell_vector_dofs + n_local_bdry_vector_dofs, 0 + n_local_cell_vector_dofs, 1, n_local_bdry_vector_dofs);

                // magnetic-magnetic dofs
                local_LTgg[iT].block(local_condensed_magnetic_offset, local_condensed_magnetic_offset, n_local_bdry_vector_dofs, n_local_bdry_vector_dofs) = local_DG.block(local_magnetic_offset + n_local_cell_vector_dofs, local_magnetic_offset + n_local_cell_vector_dofs, n_local_bdry_vector_dofs, n_local_bdry_vector_dofs);
                local_LTgg[iT].block(local_condensed_magnetic_offset, local_condensed_magnetic_offset + n_local_bdry_vector_dofs, n_local_bdry_vector_dofs, 1) = local_DG.block(local_magnetic_offset + n_local_cell_vector_dofs, local_magnetic_offset + n_local_cell_vector_dofs + n_local_bdry_vector_dofs, n_local_bdry_vector_dofs, 1);
                local_LTgg[iT].block(local_condensed_magnetic_offset + n_local_bdry_vector_dofs, local_condensed_magnetic_offset, 1, n_local_bdry_vector_dofs) = local_DG.block(local_magnetic_offset + n_local_cell_vector_dofs + n_local_bdry_vector_dofs, local_magnetic_offset + n_local_cell_vector_dofs, 1, n_local_bdry_vector_dofs);

                // magnetic pressure-pressure dofs (relaxation)
                // local_LTgg[iT](local_condensed_magnetic_offset + n_local_bdry_vector_dofs, local_condensed_magnetic_offset + n_local_bdry_vector_dofs) = local_DG(local_magnetic_offset + n_local_cell_vector_dofs + n_local_bdry_vector_dofs, local_magnetic_offset + n_local_cell_vector_dofs + n_local_bdry_vector_dofs);

                Eigen::PartialPivLU<Eigen::MatrixXd> solver;
                solver.compute(local_LTll);

                Eigen::VectorXd local_rTl = Eigen::VectorXd::Zero(n_local_implicit_dofs);
                local_rTl.head(n_local_cell_vector_dofs) = local_RHS[iT].head(n_local_cell_vector_dofs);
                local_rTl.segment(local_implicit_magnetic_offset, n_local_cell_vector_dofs) = local_RHS[iT].segment(local_magnetic_offset, n_local_cell_vector_dofs);

                Eigen::VectorXd local_inv_LTll_rTl = solver.solve(local_rTl);
                local_inv_LTll_LTlg[iT] = solver.solve(local_LTlg);

                inv_LTll_rTl.segment(iT * n_local_cell_vector_dofs, n_local_cell_vector_dofs) = local_inv_LTll_rTl.segment(0, n_local_cell_vector_dofs);
                inv_LTll_rTl.segment(n_total_cell_vector_dofs + iT * (n_local_cell_scalar_dofs - 1), n_local_cell_scalar_dofs - 1) = local_inv_LTll_rTl.segment(n_local_cell_vector_dofs, n_local_cell_scalar_dofs - 1);

                inv_LTll_rTl.segment(n_implicit_velocity_dofs + iT * n_local_cell_vector_dofs, n_local_cell_vector_dofs) = local_inv_LTll_rTl.segment(local_implicit_magnetic_offset, n_local_cell_vector_dofs);
                inv_LTll_rTl.segment(n_implicit_velocity_dofs + n_total_cell_vector_dofs + iT * (n_local_cell_scalar_dofs - 1), n_local_cell_scalar_dofs - 1) = local_inv_LTll_rTl.segment(local_implicit_magnetic_offset + n_local_cell_vector_dofs, n_local_cell_scalar_dofs - 1);
            }
        };

        // Running the local constructions in parallel
        parallel_for(n_cells, build_local_mats, threading);

        Eigen::SparseMatrix<double> inv_LTll_LTlg(n_implicit_dofs, n_unknowns + 2);
        Eigen::SparseMatrix<double> LTgl(n_unknowns + 2, n_implicit_dofs);
        Eigen::SparseMatrix<double> LTgg(n_unknowns + 2, n_unknowns + 2);

        std::vector<Eigen::Triplet<double>> triplets_inv_LTll_LTlg;
        std::vector<Eigen::Triplet<double>> triplets_LTgl;
        std::vector<Eigen::Triplet<double>> triplets_LTgg;

        Eigen::VectorXd rTg = Eigen::VectorXd::Zero(n_unknowns + 2);

        for (size_t iT = 0; iT < n_cells; ++iT)
        {
            Cell *cell = mesh_ptr->cell(iT);
            size_t n_local_faces = cell->n_faces();
            size_t n_local_bdry_vector_dofs = n_local_faces * n_local_face_vector_dofs;

            size_t local_velocity_cell_offset = 0;
            size_t local_magnetic_cell_offset = n_local_cell_vector_dofs + n_local_cell_scalar_dofs - 1;
            size_t global_velocity_cell_offset = iT * n_local_cell_vector_dofs;
            size_t global_magnetic_cell_offset = n_implicit_velocity_dofs + iT * n_local_cell_vector_dofs;

            size_t local_velocity_pressure_offset = n_local_cell_vector_dofs;
            size_t local_magnetic_pressure_offset = n_local_cell_vector_dofs + n_local_cell_scalar_dofs - 1 + n_local_cell_vector_dofs;
            size_t global_velocity_pressure_offset = n_total_cell_vector_dofs + iT * (n_local_cell_scalar_dofs - 1);
            size_t global_magnetic_pressure_offset = n_implicit_velocity_dofs + n_total_cell_vector_dofs + iT * (n_local_cell_scalar_dofs - 1);

            for (size_t iTF = 0; iTF < n_local_faces; iTF++)
            {
                size_t iF = cell->face(iTF)->global_index();

                size_t local_velocity_face_iTF_offset = iTF * n_local_face_vector_dofs;
                size_t global_velocity_face_iTF_offset = iF * n_local_face_vector_dofs;
                size_t local_magnetic_face_iTF_offset = n_local_bdry_vector_dofs + 1 + iTF * n_local_face_vector_dofs;
                size_t global_magnetic_face_iTF_offset = n_velocity_unknowns + iF * n_local_face_vector_dofs;

                rTg.segment(global_velocity_face_iTF_offset, n_local_face_vector_dofs) += local_RHS[iT].segment(n_local_cell_vector_dofs + iTF * n_local_face_vector_dofs, n_local_face_vector_dofs);
                rTg.segment(global_magnetic_face_iTF_offset, n_local_face_vector_dofs) += local_RHS[iT].segment(n_local_cell_vector_dofs + n_local_bdry_vector_dofs + n_local_cell_scalar_dofs + n_local_cell_vector_dofs + iTF * n_local_face_vector_dofs, n_local_face_vector_dofs);

                for (size_t i = 0; i < n_local_face_vector_dofs; i++)
                {
                    for (size_t j = 0; j < n_local_cell_vector_dofs; j++)
                    {
                        triplets_LTgl.emplace_back(global_velocity_face_iTF_offset + i, global_velocity_cell_offset + j, local_LTgl[iT](local_velocity_face_iTF_offset + i, local_velocity_cell_offset + j));
                        triplets_LTgl.emplace_back(global_velocity_face_iTF_offset + i, global_magnetic_cell_offset + j, local_LTgl[iT](local_velocity_face_iTF_offset + i, local_magnetic_cell_offset + j));
                        triplets_LTgl.emplace_back(global_magnetic_face_iTF_offset + i, global_velocity_cell_offset + j, local_LTgl[iT](local_magnetic_face_iTF_offset + i, local_velocity_cell_offset + j));
                        triplets_LTgl.emplace_back(global_magnetic_face_iTF_offset + i, global_magnetic_cell_offset + j, local_LTgl[iT](local_magnetic_face_iTF_offset + i, local_magnetic_cell_offset + j));

                        triplets_inv_LTll_LTlg.emplace_back(global_velocity_cell_offset + j, global_velocity_face_iTF_offset + i, local_inv_LTll_LTlg[iT](local_velocity_cell_offset + j, local_velocity_face_iTF_offset + i));
                        triplets_inv_LTll_LTlg.emplace_back(global_magnetic_cell_offset + j, global_velocity_face_iTF_offset + i, local_inv_LTll_LTlg[iT](local_magnetic_cell_offset + j, local_velocity_face_iTF_offset + i));
                        triplets_inv_LTll_LTlg.emplace_back(global_velocity_cell_offset + j, global_magnetic_face_iTF_offset + i, local_inv_LTll_LTlg[iT](local_velocity_cell_offset + j, local_magnetic_face_iTF_offset + i));
                        triplets_inv_LTll_LTlg.emplace_back(global_magnetic_cell_offset + j, global_magnetic_face_iTF_offset + i, local_inv_LTll_LTlg[iT](local_magnetic_cell_offset + j, local_magnetic_face_iTF_offset + i));
                    }
                    for (size_t j = 0; j < n_local_cell_scalar_dofs - 1; j++)
                    {
                        triplets_LTgl.emplace_back(global_velocity_face_iTF_offset + i, global_velocity_pressure_offset + j, local_LTgl[iT](local_velocity_face_iTF_offset + i, local_velocity_pressure_offset + j));
                        triplets_LTgl.emplace_back(global_magnetic_face_iTF_offset + i, global_magnetic_pressure_offset + j, local_LTgl[iT](local_magnetic_face_iTF_offset + i, local_magnetic_pressure_offset + j));

                        triplets_inv_LTll_LTlg.emplace_back(global_velocity_pressure_offset + j, global_velocity_face_iTF_offset + i, local_inv_LTll_LTlg[iT](local_velocity_pressure_offset + j, local_velocity_face_iTF_offset + i));
                        triplets_inv_LTll_LTlg.emplace_back(global_magnetic_pressure_offset + j, global_magnetic_face_iTF_offset + i, local_inv_LTll_LTlg[iT](local_magnetic_pressure_offset + j, local_magnetic_face_iTF_offset + i));
                    }
                    for (size_t jTF = 0; jTF < n_local_faces; jTF++)
                    {
                        size_t jF = cell->face(jTF)->global_index();

                        size_t local_velocity_face_jTF_offset = jTF * n_local_face_vector_dofs;
                        size_t global_velocity_face_jTF_offset = jF * n_local_face_vector_dofs;
                        size_t local_magnetic_face_jTF_offset = n_local_bdry_vector_dofs + 1 + jTF * n_local_face_vector_dofs;
                        size_t global_magnetic_face_jTF_offset = n_velocity_unknowns + jF * n_local_face_vector_dofs;
                        for (size_t j = 0; j < n_local_face_vector_dofs; j++)
                        {
                            triplets_LTgg.emplace_back(global_velocity_face_iTF_offset + i, global_velocity_face_jTF_offset + j, local_LTgg[iT](local_velocity_face_iTF_offset + i, local_velocity_face_jTF_offset + j));
                            triplets_LTgg.emplace_back(global_velocity_face_iTF_offset + i, global_magnetic_face_jTF_offset + j, local_LTgg[iT](local_velocity_face_iTF_offset + i, local_magnetic_face_jTF_offset + j));
                            triplets_LTgg.emplace_back(global_magnetic_face_iTF_offset + i, global_velocity_face_jTF_offset + j, local_LTgg[iT](local_magnetic_face_iTF_offset + i, local_velocity_face_jTF_offset + j));
                            triplets_LTgg.emplace_back(global_magnetic_face_iTF_offset + i, global_magnetic_face_jTF_offset + j, local_LTgg[iT](local_magnetic_face_iTF_offset + i, local_magnetic_face_jTF_offset + j));
                        }
                    }

                    triplets_LTgg.emplace_back(global_velocity_face_iTF_offset + i, n_total_face_vector_dofs + iT, local_LTgg[iT](local_velocity_face_iTF_offset + i, n_local_bdry_vector_dofs));
                    triplets_LTgg.emplace_back(n_total_face_vector_dofs + iT, global_velocity_face_iTF_offset + i, local_LTgg[iT](n_local_bdry_vector_dofs, local_velocity_face_iTF_offset + i));

                    triplets_LTgg.emplace_back(global_magnetic_face_iTF_offset + i, n_velocity_unknowns + n_total_face_vector_dofs + iT, local_LTgg[iT](local_magnetic_face_iTF_offset + i, n_local_bdry_vector_dofs + 1 + n_local_bdry_vector_dofs));
                    triplets_LTgg.emplace_back(n_velocity_unknowns + n_total_face_vector_dofs + iT, global_magnetic_face_iTF_offset + i, local_LTgg[iT](n_local_bdry_vector_dofs + 1 + n_local_bdry_vector_dofs, local_magnetic_face_iTF_offset + i));
                }
            }

            triplets_LTgg.emplace_back(n_total_face_vector_dofs + iT, n_unknowns, std::sqrt(cell->measure()));
            triplets_LTgg.emplace_back(n_unknowns, n_total_face_vector_dofs + iT, -std::sqrt(cell->measure()));
            triplets_LTgg.emplace_back(n_velocity_unknowns + n_total_face_vector_dofs + iT, n_unknowns + 1, std::sqrt(cell->measure()));
            triplets_LTgg.emplace_back(n_unknowns + 1, n_velocity_unknowns + n_total_face_vector_dofs + iT, -std::sqrt(cell->measure()));


            // triplets_LTgg.emplace_back(n_unknowns, n_unknowns, std::sqrt(cell->measure()));
            // triplets_LTgg.emplace_back(n_unknowns + 1, n_unknowns + 1, std::sqrt(cell->measure()));
            // triplets_LTgg.emplace_back(n_unknowns, n_unknowns, 1.0);
            // triplets_LTgg.emplace_back(n_unknowns + 1, n_unknowns + 1, 1.0);

            // triplets_LTgg.emplace_back(n_total_face_vector_dofs + iT, n_unknowns, 1.0);
            // triplets_LTgg.emplace_back(n_unknowns, n_total_face_vector_dofs + iT, -1.0);
            // triplets_LTgg.emplace_back(n_velocity_unknowns + n_total_face_vector_dofs + iT, n_unknowns + 1, 1.0);
            // triplets_LTgg.emplace_back(n_unknowns + 1, n_velocity_unknowns + n_total_face_vector_dofs + iT, -1.0);

            // triplets_LTgg.emplace_back(n_total_face_vector_dofs + iT, n_unknowns, 1.0);

            triplets_LTgg.emplace_back(n_total_face_vector_dofs + iT, n_total_face_vector_dofs + iT, std::pow((cell->diam() / double(m_K + 1)), m_K + 2));
            triplets_LTgg.emplace_back(n_velocity_unknowns + n_total_face_vector_dofs + iT, n_velocity_unknowns + n_total_face_vector_dofs + iT, std::pow((cell->diam() / double(m_K + 1)), m_K + 2));

            // triplets_LTgg.emplace_back(n_total_face_vector_dofs + iT, n_total_face_vector_dofs + iT, local_LTgg[iT](n_local_bdry_vector_dofs, n_local_bdry_vector_dofs));
            // triplets_LTgg.emplace_back(n_velocity_unknowns + n_total_face_vector_dofs + iT, n_velocity_unknowns + n_total_face_vector_dofs + iT, local_LTgg[iT](n_local_bdry_vector_dofs + 1 + n_local_bdry_vector_dofs, n_local_bdry_vector_dofs + 1 + n_local_bdry_vector_dofs));
        }

        inv_LTll_LTlg.setFromTriplets(std::begin(triplets_inv_LTll_LTlg), std::end(triplets_inv_LTll_LTlg));
        LTgl.setFromTriplets(std::begin(triplets_LTgl), std::end(triplets_LTgl));
        LTgg.setFromTriplets(std::begin(triplets_LTgg), std::end(triplets_LTgg));

        for (size_t ibF = 0; ibF < mesh_ptr->n_b_faces(); ++ibF)
        {
            size_t iF = mesh_ptr->n_i_faces() + ibF;
            size_t velocity_offset = iF * n_local_face_vector_dofs;
            size_t magnetic_offset = n_velocity_unknowns + iF * n_local_face_vector_dofs;
            for (size_t i = 0; i < n_local_face_vector_dofs; ++i)
            {
                int vel_index = int(velocity_offset + i);
                int mag_index = int(magnetic_offset + i);

                if ((m_u_bc == 'D') || ((m_u_bc == 'H') && i % dim == 0))
                {
                    LTgg.row(vel_index) *= 0.0;
                    LTgg.col(vel_index) *= 0.0;
                    LTgg.coeffRef(vel_index, vel_index) = 1.0;

                    LTgl.row(vel_index) *= 0.0;
                    inv_LTll_LTlg.col(vel_index) *= 0.0;
                }
                if ((m_b_bc == 'D') || ((m_b_bc == 'H') && i % dim == 0))
                {
                    LTgg.row(mag_index) *= 0.0;
                    LTgg.col(mag_index) *= 0.0;
                    LTgg.coeffRef(mag_index, mag_index) = 1.0;

                    LTgl.row(mag_index) *= 0.0;
                    inv_LTll_LTlg.col(mag_index) *= 0.0;
                }
            }
        }
        LTgg.prune(0.0);
        inv_LTll_LTlg.prune(0.0);
        LTgl.prune(0.0);

        LTgg.row(n_unknowns) *= visc;
        LTgg.col(n_unknowns) *= visc;
        LTgg.row(n_unknowns + 1) *= diff;
        LTgg.col(n_unknowns + 1) *= diff;

        Eigen::SparseMatrix<double> condensed_SysMat = LTgg - LTgl * inv_LTll_LTlg;

        // Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> condensed_solver;
        // condensed_solver.analyzePattern(condensed_SysMat);
        // condensed_solver.factorize(condensed_SysMat);

        Eigen::PardisoLU<Eigen::SparseMatrix<double>> condensed_solver;
        condensed_solver.compute(condensed_SysMat);

        Eigen::VectorXd RHS = rTg - LTgl * inv_LTll_rTl;

        Eigen::VectorXd condensed_terms = condensed_solver.solve(RHS);
        double residual_of_linear_system = (condensed_SysMat * condensed_terms - RHS).norm();
        // std::cout << "Residual of linear system = " << (condensed_SysMat * condensed_terms - RHS).norm() << "\n";

        Eigen::VectorXd implicit_terms = inv_LTll_rTl - inv_LTll_LTlg * condensed_terms;

        Eigen::VectorXd deltaU = Eigen::VectorXd::Zero(n_total_dofs);

        deltaU.segment(0, n_total_cell_vector_dofs) = implicit_terms.segment(0, n_total_cell_vector_dofs);
        deltaU.segment(n_total_cell_vector_dofs, n_total_face_vector_dofs) = condensed_terms.segment(0, n_total_face_vector_dofs);

        deltaU.segment(n_velocity_dofs, n_total_cell_vector_dofs) = implicit_terms.segment(n_implicit_velocity_dofs, n_total_cell_vector_dofs);
        deltaU.segment(n_velocity_dofs + n_total_cell_vector_dofs, n_total_face_vector_dofs) = condensed_terms.segment(n_velocity_unknowns, n_total_face_vector_dofs);

        for (size_t iT = 0; iT < n_cells; ++iT)
        {
            deltaU(n_total_cell_vector_dofs + n_total_face_vector_dofs + iT * n_local_cell_scalar_dofs) = condensed_terms(n_total_face_vector_dofs + iT);
            deltaU(n_velocity_dofs + n_total_cell_vector_dofs + n_total_face_vector_dofs + iT * n_local_cell_scalar_dofs) = condensed_terms(n_velocity_unknowns + n_total_face_vector_dofs + iT);

            deltaU.segment(n_total_cell_vector_dofs + n_total_face_vector_dofs + iT * n_local_cell_scalar_dofs + 1, n_local_cell_scalar_dofs - 1) = implicit_terms.segment(n_total_cell_vector_dofs + iT * (n_local_cell_scalar_dofs - 1), n_local_cell_scalar_dofs - 1);
            deltaU.segment(n_velocity_dofs + n_total_cell_vector_dofs + n_total_face_vector_dofs + iT * n_local_cell_scalar_dofs + 1, n_local_cell_scalar_dofs - 1) = implicit_terms.segment(n_implicit_velocity_dofs + n_total_cell_vector_dofs + iT * (n_local_cell_scalar_dofs - 1), n_local_cell_scalar_dofs - 1);
        }

        size_t internal_vector_dofs = n_total_cell_vector_dofs + n_total_face_vector_dofs - mesh_ptr->n_b_faces() * n_local_face_vector_dofs;
        size_t bdry_vector_dofs = mesh_ptr->n_b_faces() * n_local_face_vector_dofs;

        // set boundary dofs

        if (m_u_bc == 'D')
        {
            deltaU.segment(internal_vector_dofs, bdry_vector_dofs) = Eigen::VectorXd::Zero(bdry_vector_dofs); // Homogeneous Dirchlet BCs
        }
        else
        {
            assert(m_u_bc == 'H');
            for (size_t ibF = 0; ibF < mesh_ptr->n_b_faces(); ++ibF)
            {
                size_t offset = internal_vector_dofs + ibF * n_local_face_vector_dofs;

                // set normal components to zero
                for (size_t i = 0; i < n_local_face_scalar_dofs; i++)
                {
                    deltaU(offset + dim * i) = 0.0; // homogenous bcs
                }
            }
        }

        if (m_b_bc == 'D')
        {
            deltaU.segment(n_total_cell_vector_dofs + n_total_face_vector_dofs + n_total_cell_scalar_dofs + internal_vector_dofs, bdry_vector_dofs) = Eigen::VectorXd::Zero(bdry_vector_dofs);
        }
        else
        {
            assert(m_b_bc == 'H');
            for (size_t ibF = 0; ibF < mesh_ptr->n_b_faces(); ++ibF)
            {
                size_t offset = n_total_cell_vector_dofs + n_total_face_vector_dofs + n_total_cell_scalar_dofs + internal_vector_dofs + ibF * n_local_face_vector_dofs;

                // set normal components to zero
                for (size_t i = 0; i < n_local_face_scalar_dofs; i++)
                {
                    deltaU(offset + dim * i) = 0.0; // homogenous bcs
                }
            }
        }

        // SOL.set_values(SOL.asVectorXd() + relax * deltaU);
        SOL.set_values(SOL.asVectorXd() + deltaU);
        Eigen::VectorXd residual = Eigen::VectorXd::Zero(n_total_dofs);

        std::function<void(size_t, size_t)> compute_residual = [&](size_t start, size_t end) -> void
        {
            for (size_t iT = start; iT < end; ++iT)
            {
                Cell *cell = mesh_ptr->cell(iT);
                size_t n_local_faces = cell->n_faces();
                size_t n_local_bdry_vector_dofs = n_local_faces * n_local_face_vector_dofs;
                size_t n_local_vector_dofs = n_local_cell_vector_dofs + n_local_bdry_vector_dofs;
                size_t n_local_unknowns = 2 * (n_local_vector_dofs + n_local_cell_scalar_dofs);
                size_t local_magnetic_offset = n_local_vector_dofs + n_local_cell_scalar_dofs;

                Eigen::VectorXd uT = SOL.restr(iT);
                Eigen::VectorXd vector_terms = Eigen::VectorXd::Zero(2 * n_local_vector_dofs);
                vector_terms.head(n_local_vector_dofs) = uT.head(n_local_vector_dofs);
                vector_terms.tail(n_local_vector_dofs) = uT.segment(local_magnetic_offset, n_local_vector_dofs);

                Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n_local_unknowns, n_local_unknowns);

                A.block(0, 0, n_local_vector_dofs, n_local_vector_dofs) = visc * AT[iT];
                A.block(0, n_local_vector_dofs, n_local_vector_dofs, n_local_cell_scalar_dofs) = DT[iT].transpose();
                A.block(n_local_vector_dofs, 0, n_local_cell_scalar_dofs, n_local_vector_dofs) = -DT[iT];

                A.block(local_magnetic_offset, local_magnetic_offset, n_local_vector_dofs, n_local_vector_dofs) = diff * AT[iT];
                A.block(local_magnetic_offset, local_magnetic_offset + n_local_vector_dofs, n_local_vector_dofs, n_local_cell_scalar_dofs) = DT[iT].transpose();
                A.block(local_magnetic_offset + n_local_vector_dofs, local_magnetic_offset, n_local_cell_scalar_dofs, n_local_vector_dofs) = -DT[iT];

                // Eigen::MatrixXd T = convective_term(cell, TTk[iT], vector_terms);
                Eigen::MatrixXd T = convective_term(cell, Tplus_k[iT], Tminus_k[iT], vector_terms);
                Eigen::MatrixXd BMAT = Eigen::MatrixXd::Zero(n_local_unknowns, n_local_unknowns);
                BMAT.block(0, 0, n_local_vector_dofs, n_local_vector_dofs) = T.topLeftCorner(n_local_vector_dofs, n_local_vector_dofs);
                BMAT.block(local_magnetic_offset, 0, n_local_vector_dofs, n_local_vector_dofs) = T.bottomLeftCorner(n_local_vector_dofs, n_local_vector_dofs);
                BMAT.block(0, local_magnetic_offset, n_local_vector_dofs, n_local_vector_dofs) = T.topRightCorner(n_local_vector_dofs, n_local_vector_dofs);
                BMAT.block(local_magnetic_offset, local_magnetic_offset, n_local_vector_dofs, n_local_vector_dofs) = T.bottomRightCorner(n_local_vector_dofs, n_local_vector_dofs);

                Eigen::MatrixXd GMAT = A + 0.5 * BMAT;

                Eigen::VectorXd G = GMAT * uT;
                G.segment(0, n_local_cell_vector_dofs) -= LTf[iT];
                G.segment(local_magnetic_offset, n_local_cell_vector_dofs) -= LTg[iT];

                residual.segment(iT * n_local_cell_vector_dofs, n_local_cell_vector_dofs) = G.head(n_local_cell_vector_dofs);
                // residual.segment(n_total_cell_vector_dofs + n_total_face_vector_dofs + iT * n_local_cell_scalar_dofs, n_local_cell_scalar_dofs) = G.segment(n_local_vector_dofs, n_local_cell_scalar_dofs);

                size_t global_magnetic_offset = n_total_cell_vector_dofs + n_total_face_vector_dofs + n_total_cell_scalar_dofs;

                residual.segment(global_magnetic_offset + iT * n_local_cell_vector_dofs, n_local_cell_vector_dofs) = G.segment(local_magnetic_offset, n_local_cell_vector_dofs);
                // residual.segment(global_magnetic_offset + n_total_cell_vector_dofs + n_total_face_vector_dofs + iT * n_local_cell_scalar_dofs, n_local_cell_scalar_dofs) = G.tail(n_local_cell_scalar_dofs);

                // for (size_t iTF = 0; iTF < n_local_faces; ++iTF)
                // {
                //     if (cell->face(iTF)->is_boundary())
                //     {
                //         continue;
                //     }
                //     size_t iF = cell->face(iTF)->global_index();
                //     residual.segment(n_total_cell_vector_dofs + iF * n_local_face_vector_dofs, n_local_face_vector_dofs) += G.segment(n_local_cell_vector_dofs + iTF * n_local_face_vector_dofs, n_local_face_vector_dofs);
                //     residual.segment(global_magnetic_offset + n_total_cell_vector_dofs + iF * n_local_face_vector_dofs, n_local_face_vector_dofs) += G.segment(local_magnetic_offset + n_local_cell_vector_dofs + iTF * n_local_face_vector_dofs, n_local_face_vector_dofs);
                // }
            }
        };

        // Running the local constructions in parallel
        parallel_for(n_cells, compute_residual, threading);

        residual_error = residual.norm() / exact_rhs_norm;

        ++i;

        std::cout << "     It. " << i << ":";
        std::cout << " Res = " << residual_error << ", Solver error = " << residual_of_linear_system << "\n";

        if (i > 500)
        {
            std::cout << "Failed to converge\n";
            // std::cout << "Failed to converge, best residual = " << best_residual << "\n";
            // SOL.set_values(best_guess);
            break;
        }
    } while (residual_error > tol);

    solving_timer.stop();
    std::cout << "     Solving time = " << double(solving_timer.elapsed().wall) * 1E-9 << "\n";

    return SOL;
}

#endif
