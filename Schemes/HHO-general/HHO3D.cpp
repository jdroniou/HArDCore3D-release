#include "HHO3D.hpp"

HHO3D::HHO3D(HybridCore &hho, const size_t L, const size_t K, const bool use_threads, size_t doeT, size_t doeF)
    : m_hho(hho),
      m_L(L),
      m_K(K),
      m_use_threads(use_threads),
      m_doeT(doeT),
      m_doeF(doeF)
{
    GlobMat.resize(n_total_face_dofs, n_total_face_dofs);
    ScBeMat.resize(n_total_cell_dofs, n_total_face_dofs);
    AT.resize(n_cells);
    if (m_doeT == 0)
        m_doeT = m_L + m_K + 1;
    if (m_doeF == 0)
        m_doeF = 2 * m_K + 1;
}

void HHO3D::set_global_operator(const MatrixFType &AT)
{
    global_operator = [&](Cell *cell, ElementQuad &elquad) -> Eigen::MatrixXd {
        return AT(cell, elquad);
    };
}

void HHO3D::set_load_vector(const VectorFType &bT)
{
    load_vector = [&](Cell *cell, ElementQuad &elquad) -> Eigen::VectorXd {
        return bT(cell, elquad);
    };
}

void HHO3D::plot(const std::string plot_file, const UVector &sol, const FType<double> &exact_sol)
{
    // Plot if plot_file given and cells less than 5000
    if (plot_file != "" && n_cells < 5000)
    {
        VectorRd coord;

        std::vector<Vertex *> vertices = mesh_ptr->get_vertices();
        size_t n_verts = mesh_ptr->n_vertices();
        Eigen::VectorXd interpolant(n_verts);

        // Calculate interpolant of exact solution
        for (size_t i = 0; i < n_verts; i++)
        {
            coord = vertices[i]->coords();
            interpolant(i) = exact_sol(coord);
        }

        VtuWriter plotdata(mesh_ptr);

        // Calculate cell and face values on the vertices
        Eigen::VectorXd approx_sol_vertex_byF = m_hho.VertexValues(sol, "face");
        Eigen::VectorXd approx_sol_vertex_byT = m_hho.VertexValues(sol, "cell");

        // Plot solutions. T == cell, F == face, E == exact
        plotdata.write_to_vtu("T-" + plot_file + ".vtu", approx_sol_vertex_byF);
        plotdata.write_to_vtu("F-" + plot_file + ".vtu", approx_sol_vertex_byT);
        plotdata.write_to_vtu("E-" + plot_file + ".vtu", interpolant);
    }
}

VectorFType HHO3D::standard_load_vector(const CellFType<double> &source)
{
    VectorFType bT = [&](Cell *cell, ElementQuad &elquad) -> Eigen::VectorXd {
        FType<double> source_on_cell = [&](const VectorRd x) -> double { return source(x, cell); };

        const size_t n_local_faces = cell->n_faces();
        const size_t n_local_dofs = n_local_cell_dofs + n_local_faces * n_local_face_dofs;

        QuadratureRule quadT = elquad.get_quadT();
        BasisQuad<double> phiT_quadT = elquad.get_phiT_quadT();

        Eigen::VectorXd source_vec = Eigen::VectorXd::Zero(n_local_dofs);

        source_vec.head(n_local_cell_dofs) = integrate(source_on_cell, phiT_quadT, quadT, n_local_cell_dofs);
        return source_vec;
    };

    return bT;
}

VectorFType HHO3D::standard_load_vector(const CellFType<double> &source, const CellFType<VectorRd> &f, const BoundaryConditions &BC)
{
    VectorFType bT = [&](Cell *cell, ElementQuad &elquad) -> Eigen::VectorXd {
        FType<double> source_on_cell = [&](const VectorRd x) -> double { return source(x, cell); };

        const size_t n_local_faces = cell->n_faces();
        const size_t n_local_dofs = n_local_cell_dofs + n_local_faces * n_local_face_dofs;

        QuadratureRule quadT = elquad.get_quadT();
        BasisQuad<double> phiT_quadT = elquad.get_phiT_quadT();

        Eigen::VectorXd source_vec = Eigen::VectorXd::Zero(n_local_dofs);

        source_vec.head(n_local_cell_dofs) = integrate(source_on_cell, phiT_quadT, quadT, n_local_cell_dofs);

        if (cell->is_boundary())
        {
            for (size_t iTF = 0; iTF < n_local_faces; iTF++)
            {
                Face *face = cell->face(iTF);
                if (BC.type(*face) == "neu")
                {
                    const size_t offset_F = n_local_cell_dofs + iTF * n_local_face_dofs;
                    const size_t iF = face->global_index();
                    const VectorRd &nTF = cell->face_normal(iTF);

                    QuadratureRule quadF = generate_quadrature_rule(*face, 2 * m_K + 2);
                    BasisQuad<double> phiF_quadF = evaluate_quad<Function>::compute(m_hho.FaceBasis(iF), quadF);

                    FType<double> f_dot_nTF = [&](const VectorRd x) -> double {
                        return nTF.dot(f(x, cell));
                    };

                    source_vec.segment(offset_F, n_local_face_dofs) += integrate(f_dot_nTF, phiF_quadF, quadF);
                }
            }
        }
        return source_vec;
    };

    return bT;
}

VectorFType HHO3D::standard_load_vector(const CellFType<double> &source, const FType<double> &f, const BoundaryConditions &BC)
{
    VectorFType bT = [&](Cell *cell, ElementQuad &elquad) -> Eigen::VectorXd {
        FType<double> source_on_cell = [&](const VectorRd x) -> double { return source(x, cell); };

        const size_t n_local_faces = cell->n_faces();
        const size_t n_local_dofs = n_local_cell_dofs + n_local_faces * n_local_face_dofs;

        QuadratureRule quadT = elquad.get_quadT();
        BasisQuad<double> phiT_quadT = elquad.get_phiT_quadT();

        Eigen::VectorXd source_vec = Eigen::VectorXd::Zero(n_local_dofs);

        source_vec.head(n_local_cell_dofs) = integrate(source_on_cell, phiT_quadT, quadT, n_local_cell_dofs);

        if (cell->is_boundary())
        {
            for (size_t iTF = 0; iTF < n_local_faces; iTF++)
            {
                Face *face = cell->face(iTF);
                if (BC.type(*face) == "neu")
                {
                    const size_t offset_F = n_local_cell_dofs + iTF * n_local_face_dofs;
                    const size_t iF = face->global_index();

                    QuadratureRule quadF = generate_quadrature_rule(*face, 2 * m_K + 2);
                    BasisQuad<double> phiF_quadF = evaluate_quad<Function>::compute(m_hho.FaceBasis(iF), quadF);

                    source_vec.segment(offset_F, n_local_face_dofs) += integrate(f, phiF_quadF, quadF);
                }
            }
        }
        return source_vec;
    };

    return bT;
}

void HHO3D::set_dirichlet(const FType<double> &f, const size_t n_dir_faces)
{
    UDir = Eigen::VectorXd::Zero(n_local_face_dofs * n_dir_faces);
    for (size_t idF = 0; idF < n_dir_faces; idF++)
    {
        Face *face = mesh_ptr->face(n_faces - n_dir_faces + idF);
        const size_t iF = face->global_index();
        QuadratureRule quadF = generate_quadrature_rule(*face, 2 * m_K + 2);
        BasisQuad<double> phiF_quadF = evaluate_quad<Function>::compute(m_hho.FaceBasis(iF), quadF);
        UDir.segment(idF * n_local_face_dofs, n_local_face_dofs) = l2_projection<HybridCore::PolyFaceBasisType>(f, m_hho.FaceBasis(iF), quadF, phiF_quadF);
    }
}

void HHO3D::set_dirichlet(const size_t n_dir_faces)
{
    UDir = Eigen::VectorXd::Zero(n_local_face_dofs * n_dir_faces);
}

void HHO3D::assemble()
{
    boost::timer::cpu_timer timer;
    timer.start();

    std::vector<Eigen::Triplet<double>> triplets_GlobMat;
    std::vector<Eigen::Triplet<double>> triplets_ScBe;
#ifdef WITH_FHHOS4
    // With fhhos4: the triplets of A_TT, A_TF, A_FT, A_FF, A_TD, A_FD, and the local right-hand sides (one per cell);
    // the Dirichlet faces are the last ones
    std::vector<Eigen::Triplet<double>> triplets_TT, triplets_TF, triplets_FT, triplets_FF, triplets_TD, triplets_FD;
    std::vector<Eigen::VectorXd> cell_rhs(n_cells);
    const size_t n_unknowns = n_total_face_dofs - UDir.size();
    if (m_solver_type == "fhhos4")
    {
        b_T = Eigen::VectorXd::Zero(n_total_cell_dofs);
        b_F = Eigen::VectorXd::Zero(n_unknowns);
    }
#endif

    std::vector<Eigen::MatrixXd> invATT_ATF(n_cells);
    std::vector<Eigen::VectorXd> cell_source(n_cells);
    std::vector<Eigen::MatrixXd> MatF(n_cells);

    std::vector<Cell *> cells = mesh_ptr->get_cells();

    // Construct the local matrices
    std::function<void(size_t, size_t)> construct_all_local_contributions = [&](size_t start, size_t end) -> void {
        for (size_t iT = start; iT < end; iT++)
        {
            Cell *cell = cells[iT];
            const size_t n_local_faces = cell->n_faces();
            const size_t n_face_dofs = n_local_faces * n_local_face_dofs;

            ElementQuad elquad(m_hho, iT, m_doeT, m_doeF);

            const Eigen::VectorXd bT = load_vector(cell, elquad);
            AT[iT] = global_operator(cell, elquad);
#ifdef WITH_FHHOS4
            if (m_solver_type == "fhhos4")
            {
                // No static condensation (fhhos4 condenses the system): AT and bT assembled as they are, below
                cell_rhs[iT] = bT;
                continue;
            }
#endif

            const Eigen::MatrixXd ATT = AT[iT].topLeftCorner(n_local_cell_dofs, n_local_cell_dofs);
            const Eigen::MatrixXd ATF = AT[iT].topRightCorner(n_local_cell_dofs, n_face_dofs);
            const Eigen::MatrixXd AFT = AT[iT].bottomLeftCorner(n_face_dofs, n_local_cell_dofs);
            const Eigen::MatrixXd AFF = AT[iT].bottomRightCorner(n_face_dofs, n_face_dofs);

            Eigen::PartialPivLU<Eigen::MatrixXd> invATT;
            invATT.compute(ATT);

            invATT_ATF[iT] = invATT.solve(ATF);
            const Eigen::VectorXd invATT_bTcell = invATT.solve(bT.head(n_local_cell_dofs));

            MatF[iT] = AFF - AFT * invATT_ATF[iT];
            cell_source[iT] = bT.tail(n_face_dofs) - AFT * invATT_bTcell;

            ScRHS.segment(iT * n_local_cell_dofs, n_local_cell_dofs) = invATT_bTcell;
        }
    };

    // Running the local constructions in parallel
    parallel_for(n_cells, construct_all_local_contributions, m_use_threads);

    // Set the local matrices into the correct positions of the global matrices
    for (size_t iT = 0; iT < n_cells; iT++)
    {
        Cell *cell = cells[iT];
        const size_t n_local_faces = cell->n_faces();
#ifdef WITH_FHHOS4
        if (m_solver_type == "fhhos4")
        {
            // The local matrix and right-hand side by blocks, the cell unknowns first; the Dirichlet faces: their columns
            // in A_TD and A_FD, their equations dropped
            const size_t nT = n_local_cell_dofs, offset_T = iT * nT;
            b_T.segment(offset_T, nT) = cell_rhs[iT].head(nT);
            for (size_t i = 0; i < nT; i++)
                for (size_t j = 0; j < nT; j++)
                    triplets_TT.emplace_back(offset_T + i, offset_T + j, AT[iT](i, j));
            for (size_t jTF = 0; jTF < n_local_faces; jTF++)
            {
                const size_t jF = cell->face(jTF)->global_index();
                for (size_t jk = 0; jk < n_local_face_dofs; jk++)
                {
                    const size_t jLocal = nT + jTF * n_local_face_dofs + jk;
                    const size_t jGlobal = jF * n_local_face_dofs + jk;
                    const bool dirichlet = jGlobal >= n_unknowns;
                    const size_t col = dirichlet ? jGlobal - n_unknowns : jGlobal;
                    if (!dirichlet)
                    {
                        b_F(jGlobal) += cell_rhs[iT](jLocal);
                        for (size_t i = 0; i < nT; i++)
                            triplets_FT.emplace_back(jGlobal, offset_T + i, AT[iT](jLocal, i));
                    }
                    for (size_t i = 0; i < nT; i++)
                        (dirichlet ? triplets_TD : triplets_TF).emplace_back(offset_T + i, col, AT[iT](i, jLocal));
                    for (size_t iTF = 0; iTF < n_local_faces; iTF++)
                    {
                        const size_t iF = cell->face(iTF)->global_index();
                        for (size_t ik = 0; ik < n_local_face_dofs; ik++)
                        {
                            const size_t iGlobal = iF * n_local_face_dofs + ik;
                            if (iGlobal < n_unknowns)
                                (dirichlet ? triplets_FD : triplets_FF).emplace_back(iGlobal, col, AT[iT](nT + iTF * n_local_face_dofs + ik, jLocal));
                        }
                    }
                }
            }
            continue;
        }
#endif
        for (size_t iTF = 0; iTF < n_local_faces; iTF++)
        {
            const size_t iF = cell->face(iTF)->global_index();
            for (size_t ik = 0; ik < n_local_face_dofs; ik++)
            {
                const size_t iLocal = iTF * n_local_face_dofs + ik;
                const size_t iGlobal = iF * n_local_face_dofs + ik;
                GlobRHS(iGlobal) += cell_source[iT](iLocal);
                for (size_t i = 0; i < n_local_cell_dofs; i++)
                {
                    triplets_ScBe.emplace_back(iT * n_local_cell_dofs + i, iGlobal, invATT_ATF[iT](i, iLocal));
                }
                for (size_t jTF = 0; jTF < n_local_faces; jTF++)
                {
                    const size_t jF = cell->face(jTF)->global_index();
                    for (size_t jk = 0; jk < n_local_face_dofs; jk++)
                    {
                        const size_t jLocal = jTF * n_local_face_dofs + jk;
                        const size_t jGlobal = jF * n_local_face_dofs + jk;
                        triplets_GlobMat.emplace_back(iGlobal, jGlobal, MatF[iT](iLocal, jLocal));
                    }
                }
            }
        }
    }

    // Construct the global matrices
#ifdef WITH_FHHOS4
    if (m_solver_type == "fhhos4")
    {
        // The system before static condensation (the triplets of GlobMat and ScBeMat below are empty)
        A_TT.resize(n_total_cell_dofs, n_total_cell_dofs);
        A_TF.resize(n_total_cell_dofs, n_unknowns);
        A_FT.resize(n_unknowns, n_total_cell_dofs);
        A_FF.resize(n_unknowns, n_unknowns);
        A_TD.resize(n_total_cell_dofs, UDir.size());
        A_FD.resize(n_unknowns, UDir.size());
        A_TT.setFromTriplets(std::begin(triplets_TT), std::end(triplets_TT));
        A_TF.setFromTriplets(std::begin(triplets_TF), std::end(triplets_TF));
        A_FT.setFromTriplets(std::begin(triplets_FT), std::end(triplets_FT));
        A_FF.setFromTriplets(std::begin(triplets_FF), std::end(triplets_FF));
        A_TD.setFromTriplets(std::begin(triplets_TD), std::end(triplets_TD));
        A_FD.setFromTriplets(std::begin(triplets_FD), std::end(triplets_FD));
    }
#endif
    GlobMat.setFromTriplets(std::begin(triplets_GlobMat), std::end(triplets_GlobMat));
    ScBeMat.setFromTriplets(std::begin(triplets_ScBe), std::end(triplets_ScBe));

    assembly_time = timer.elapsed().wall;
}

UVector HHO3D::solve()
{
    /*
        * GlobMat = AFF - AFT * INV_ATT * ATF
        * GlobRHS = BF - AFT * INV_ATT_BT
        * 
        * ScBeMat = INV_ATT * ATF 
        * ScRHS   = INV_ATT_BT
    */

#ifdef WITH_FHHOS4
    if (m_solver_type == "fhhos4")
        return solve_fhhos4();
#endif

    boost::timer::cpu_timer timer;
    timer.start();

    const size_t n_fixed_dofs = UDir.size();
    const size_t n_unknowns = n_total_face_dofs - n_fixed_dofs;

    SysMat = GlobMat.topLeftCorner(n_unknowns, n_unknowns);

    // (D == Dirichlet, I == Internal/Neumann, T == Cell, F == Face)
    //
    //                      = BI - AIT * INV_ATT * BT        - (AID - AIT * INV_ATT * ATD * UD)
    Eigen::VectorXd B = GlobRHS.segment(0, n_unknowns) - GlobMat.topRightCorner(n_unknowns, n_fixed_dofs) * UDir;

    Eigen::VectorXd xF = Eigen::VectorXd::Zero(n_unknowns);
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(SysMat);

    // xF = INV( AII - AIT * INV_ATT * ATI ) * ( BI - AIT * INV_ATT * BT - (AID - AIT * INV_ATT * ATD * UD) )
    xF = solver.solve(B);
    std::cout << "     [solver] #iterations: " << solver.iterations() << ", estimated error: " << solver.error() << std::endl;
    solving_error = (SysMat * xF - B).norm();

    Eigen::VectorXd Xh = Eigen::VectorXd::Zero(n_total_cell_dofs + n_total_face_dofs);
    Xh.tail(n_fixed_dofs) = UDir;
    Xh.segment(n_total_cell_dofs, n_unknowns) = xF;

    //                         = INV_ATT_BT - INV_ATT * ATF * UF
    Xh.head(n_total_cell_dofs) = ScRHS - ScBeMat * Xh.tail(n_total_face_dofs);

    solving_time = timer.elapsed().wall;

    return UVector(Xh, *mesh_ptr, m_L, m_K);
}

#ifdef WITH_FHHOS4
UVector HHO3D::solve_fhhos4()
{
    boost::timer::cpu_timer timer;
    timer.start();

    const size_t n_fixed_dofs = UDir.size();
    const size_t n_unknowns = n_total_face_dofs - n_fixed_dofs;

    // The system before static condensation (see assemble()), the Dirichlet values moved to the right-hand sides
    const Eigen::VectorXd B_T = b_T - A_TD * UDir;
    const Eigen::VectorXd B_F = b_F - A_FD * UDir;

    // The interpolation of the function 1 on the bases of the cells and faces: 1/phi_0, phi_0 being their first function
    const size_t nF = n_local_face_dofs;
    Eigen::VectorXd cell_one(n_cells), face_one(n_unknowns / nF);
    for (size_t iT = 0; iT < n_cells; iT++)
        cell_one(iT) = 1. / m_hho.CellBasis(iT).function(0, mesh_ptr->cell(iT)->center_mass());
    for (size_t iF = 0; iF < n_unknowns / nF; iF++)
        face_one(iF) = 1. / m_hho.FaceBasis(iF).function(0, mesh_ptr->face(iF)->center_mass());

    // fhhos4 condenses the system and computes the cell unknowns too
    Eigen::VectorXd xT(n_total_cell_dofs), xF = Eigen::VectorXd::Zero(n_unknowns);
    fhhos4::Solver solver;
    solver.Dimension = 3;
    solver.FaceDegree = m_K;
    solver.CellDegree = m_L;
    solver.Tolerance = 1e-12;
    solver.SetupNonSym(A_TT, A_TF, A_FT, A_FF, cell_one, face_one);
    fhhos4::Result result = solver.Solve(B_T, B_F, xT, xF);
    std::cout << "     [solver] #iterations: " << result.Iterations << ", relative residual: " << result.RelativeResidual << std::endl;
    // Residual of the equations of the faces, that of the condensed system (those of the cells are solved exactly)
    solving_error = (A_FF * xF + A_FT * xT - B_F).norm();

    Eigen::VectorXd Xh(n_total_cell_dofs + n_total_face_dofs);
    Xh.head(n_total_cell_dofs) = xT;
    Xh.segment(n_total_cell_dofs, n_unknowns) = xF;
    Xh.tail(n_fixed_dofs) = UDir;

    solving_time = timer.elapsed().wall;

    return UVector(Xh, *mesh_ptr, m_L, m_K);
}
#endif

// UVector HHO3D::neumann_solve()
// {

// }

double HHO3D::energy_norm(const UVector Xh)
{
    double value = 0.0;

    for (size_t iT = 0; iT < n_cells; iT++)
    {
        Eigen::VectorXd XTF = Xh.restr(iT);
        value += XTF.transpose() * AT[iT] * XTF;
    }

    return sqrt(value);
}