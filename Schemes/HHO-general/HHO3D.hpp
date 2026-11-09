// Helper methods and assemble and solve routines for implementing 3D Hybrid High Order schemes.

#include <basis.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <parallel_for.hpp>

#include <mesh_builder.hpp>
#include <hybridcore.hpp>
#include <elementquad.hpp>

#include <TestCase/TestCase.hpp>
#include <BoundaryConditions/BoundaryConditions.hpp>

#include <boost/timer/timer.hpp>
#include <vtu_writer.hpp>
#ifdef WITH_FHHOS4
#include <fhhos4>
#endif
#ifndef _HHO3D_HPP
#define _HHO3D_HPP

/*!	
 * @defgroup HHO3D 
 * @brief Class providing helper methods and assemble and solve routines for general 3D HHO schemes
 */

/*!
 * \addtogroup HHO3D
 * @{
 */

// ----------------------------------------------------------------------------
//                            HHO3D class definition
// ----------------------------------------------------------------------------

/** The HHO3D class contains methods applicable to a general 3D HHO model. It contains assemble and
  * solve routines usable by any scheme. The global operator and load vector are passed to the model
  * through setters and are then accessible by assemble and solve. The class also contains helper
  * methods for writing HHO schemes. 
 **/

typedef std::function<Eigen::MatrixXd(Cell *, ElementQuad &)> MatrixFType; ///< type for the global operator as a function of a Cell and an ElementQuad
typedef std::function<Eigen::VectorXd(Cell *, ElementQuad &)> VectorFType; ///< type for the load vector as a function of a Cell and an ElementQuad

class HHO3D
{
public:
    /**@brief Class constructor: initialises the model by providing a HybridCore object, and the 
      * exact solution and boundary conditions of the model.
     **/
    HHO3D(
        HybridCore &,                  ///< A reference to the HybridCore object containing mesh data
        const size_t,                  ///< Cell polynomial degree
        const size_t,                  ///< Face polynomial degree
        const bool use_threads = true, ///< Optional argument to indicate if threads should be used
        size_t doeT = 0,               ///< Optional argument to set cell quadrature degree. Default is L + K + 1
        size_t doeF = 0                ///< Optional argument to set face quadrature degree. Default is 2 * K + 1
    );

    /// A general assemble routine that calculates the statically condensed matrices required by solve
    void assemble();

    /// Solves the statically condensed system
    UVector solve();
#ifdef WITH_FHHOS4
    /// Solves the system before static condensation with fhhos4 (solve() calls it with the solver fhhos4)
    UVector solve_fhhos4();
#endif

    /// Solves the system when the model is ill posed (not yet running)
    UVector neumann_solve();

    /// Returns the energy norm of a given UVector
    double energy_norm(const UVector);

    /// Set the global operator
    void set_global_operator(const MatrixFType &);

    /// Set the load vector
    void set_load_vector(const VectorFType &);

    /// Plot the numerical and exact solutions
    void plot(
        const std::string,    ///< Plot file
        const UVector &,      ///< Numerical solution
        const FType<double> & ///< Exact solution
    );

    /// Returns the standard load vector (f, v_T)_T with no Neumann boundary conditions
    VectorFType standard_load_vector(
        const CellFType<double> & ///< Source term
    );

    /// Returns the standard load vector (f, v_T)_T
    VectorFType standard_load_vector(
        const CellFType<double> &,   ///< Source term
        const CellFType<VectorRd> &, ///< Gradient to be dotted with the normal vector on each Neumann face
        const BoundaryConditions &   ///< A reference to the BC to determine Neumann faces
    );

    /// Returns the standard load vector (f, v_T)_T
    VectorFType standard_load_vector(
        const CellFType<double> &, ///< Source term
        const FType<double> &,     ///< Function on the Neumann BC faces
        const BoundaryConditions & ///< A reference to the BC to determine Neumann faces
    );

    /// Set the Dirichlet boundary conditions
    void set_dirichlet(
        const FType<double> &, ///< Function on the Dirichlet BC faces
        const size_t           ///< Number of Dirichlet faces
    );

    /// Set the Dirichlet boundary condition to zero
    void set_dirichlet(
        const size_t ///< Number of Dirichlet faces
    );

    /// Set the linear solver: bicgstab (default) or fhhos4
    inline void set_solver_type(const std::string &solver_type)
    {
#ifndef WITH_FHHOS4
        if (solver_type == "fhhos4")
        {
            std::cerr << "HArDCore3D was built without fhhos4 (-DUSE_FHHOS4=ON)" << std::endl;
            exit(1);
        }
#endif
        m_solver_type = solver_type;
    };

    /// Return the (statically condensed) matrix system
    inline Eigen::SparseMatrix<double> get_SysMat()
    {
        return SysMat;
    };

    /// CPU time to assemble the scheme
    inline double get_assembly_time() const
    {
        return double(assembly_time) * pow(10, -9);
    };

    /// CPU time to solve the scheme
    inline double get_solving_time() const
    {
        return double(solving_time) * pow(10, -9);
    };

    /// Residual after solving the scheme
    inline double get_solving_error() const
    {
        return solving_error;
    };

private:
    HybridCore &m_hho;
    const size_t m_L;
    const size_t m_K;
    const bool m_use_threads;
    size_t m_doeT;
    size_t m_doeF;

    MatrixFType global_operator;
    VectorFType load_vector;

    const Mesh *mesh_ptr = m_hho.get_mesh();

    const size_t n_cells = mesh_ptr->n_cells();
    const size_t n_faces = mesh_ptr->n_faces();

    const size_t n_local_cell_dofs = DimPoly<Cell>(m_L);
    const size_t n_local_face_dofs = DimPoly<Face>(m_K);

    const size_t n_total_cell_dofs = n_local_cell_dofs * n_cells;
    const size_t n_total_face_dofs = n_local_face_dofs * n_faces;
    const size_t n_total_dofs = n_total_cell_dofs + n_total_face_dofs;

    double solving_error;
    size_t solving_time;
    size_t assembly_time;

    std::vector<Eigen::MatrixXd> AT;

    Eigen::VectorXd GlobRHS = Eigen::VectorXd::Zero(n_total_face_dofs);
    Eigen::VectorXd ScRHS = Eigen::VectorXd::Zero(n_total_cell_dofs);
    Eigen::VectorXd UDir;

    Eigen::SparseMatrix<double> GlobMat;
    Eigen::SparseMatrix<double> SysMat;
    Eigen::SparseMatrix<double> ScBeMat;
    std::string m_solver_type = "bicgstab";
#ifdef WITH_FHHOS4
    // With fhhos4 (no static condensation): the system before it, by blocks (T: cell unknowns, F: face unknowns, D: the
    // values of the Dirichlet faces), the matrices A_TT (block diagonal), A_TF, A_FT, A_FF, the couplings with the
    // Dirichlet values A_TD, A_FD, and the right-hand sides b_T, b_F (without the Dirichlet values: solve_fhhos4())
    Eigen::SparseMatrix<double, Eigen::RowMajor> A_TT, A_TF, A_FT, A_FF, A_TD, A_FD;
    Eigen::VectorXd b_T, b_F;
#endif
};

#endif
