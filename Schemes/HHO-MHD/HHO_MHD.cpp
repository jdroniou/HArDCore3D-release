#include "HHO_MHD.hpp"

bool program_options(int, const char **, std::string &, size_t &, size_t &, size_t &, double &, double &, char &, char &, size_t &, size_t &, std::string &, bool &, double &);

int main(const int argc, const char **argv)
{
    size_t L, K, u_id, b_id, p_id;
    char u_bc_id, b_bc_id;
    double visc, diff, tol;
    std::string mesh_name, plot_file;
    bool use_threads;

    if (!program_options(argc, argv, mesh_name, u_id, b_id, p_id, visc, diff, u_bc_id, b_bc_id, L, K, plot_file, use_threads, tol))
    {
        return 0;
    }

    std::string mesh_dir = "../../meshes/";
    std::string mesh_file = mesh_dir + mesh_name;
    std::cout << "\n[MeshBuilder] Loading " + mesh_file + "\n";
    MeshBuilder builder = MeshBuilder(mesh_file);
    std::unique_ptr<Mesh> mesh_ptr;
    try
    {
        mesh_ptr = builder.build_the_mesh();
    }
    catch (std::string msg)
    {
        std::cerr << msg;
    }

    BoundaryConditions BC("D", *mesh_ptr.get());
    BC.reorder_faces();

    // tol = 0.5 * std::pow(mesh_ptr->h_max(), K + 2);

    // Print the data
    std::cout << "\n[Scheme] Data:\n";
    std::cout << "     No. cells = " << mesh_ptr->n_cells() << ", No. internal faces = " << mesh_ptr->n_i_faces() << ", No. vertices = " << mesh_ptr->n_vertices() << "\n";
    std::cout << "     Regularity: Reg = " << mesh_ptr->regularity()[0] << ", Skew = " << mesh_ptr->regularity()[1] <<  ", Size = " << mesh_ptr->h_max() << "\n";
    std::cout << "     TestCase: Velocity = " << u_id << "; Magnetic = " << b_id  << "; Pressure = " << p_id << "; Viscosity = " << visc << "; Diffusivity = " << diff << "; Velocity BC = " << u_bc_id << "; Magnetic BC = " << b_bc_id << "\n";
    std::cout << "     Mesh = " << mesh_name << "\n";
    std::cout << "     Degrees: face = " << K << "; cell = " << L << "\n\n";
    std::cout << "     Using threads = " << (use_threads ? "true" : "false") << "\n";

    // Load hybrid core
    HybridCore hho = HybridCore(mesh_ptr.get(), std::max(K + 1, L), K, use_threads);

    // Load model
    std::cout << "\n[Scheme] Loading model\n";
    MHDModel model(hho, L, K, u_bc_id, b_bc_id);
    MHDTests mhd_tc(u_id, b_id, p_id, visc, diff);

    std::cout << "\n[Scheme] Solving\n";
    // SolutionVector sol = model.solve();

    // exit(1);
    // Eigen::VectorXd init = Eigen::VectorXd::Zero;

    // SolutionVector sol = model.solve(mhd_tc.velocity_source(), mhd_tc.magnetic_source(), visc, diff, tol, use_threads);
    SolutionVector sol = model.solve_with_static_cond(mhd_tc.velocity_source(), mhd_tc.magnetic_source(), visc, diff, tol, use_threads);
    // SolutionVector sol = model.solve(mhd_tc.velocity_source(), mhd_tc.magnetic_source(), intepolant.asVectorXd(), visc, diff, tol, use_threads);

    std::cout << "\n[Scheme] Computing Intepolant\n";
    SolutionVector intepolant = model.global_interpolant(mhd_tc.velocity(), mhd_tc.pressure(), mhd_tc.magnetic());

    std::cout << "\n[Scheme] Computing Errors\n";
    std::vector<double> errors = model.compute_errors(intepolant, sol, visc, diff);

    std::cout << "     Energy Error (velocity) = " << errors[0] << "\n";
    std::cout << "     Energy Error (magnetic) = " << errors[1] << "\n";
    std::cout << "     Energy Error (pressure) = " << errors[2] << "\n";
    std::cout << "     Energy Error (lagrange) = " << errors[3] << "\n";

    std::cout << "     L2 Error (velocity) = " << errors[4] << "\n";
    std::cout << "     L2 Error (magnetic) = " << errors[5] << "\n";

    std::ofstream out("results.txt");
    // out << "BC: " << bdry_id << '\n';
    out << "Velocity: " << u_id << '\n';
    out << "Magnetic: " << b_id << '\n';
    out << "Pressure: " << p_id << '\n';
    out << "Viscosity: " << visc << '\n';
    out << "Diffusivity: " << diff << '\n';
    out << "Mesh: " << mesh_name << "\n";
    out << "FaceDegree: " << K << "\n";
    out << "CellDegree: " << L << "\n";
    out << "VelocityEnergyError: " <<  errors[0] << "\n";
    out << "MagneticEnergyError: " <<  errors[1] << "\n";
    out << "PressureEnergyError: " <<  errors[2] << "\n";
    out << "LagrangeEnergyError: " <<  errors[3] << "\n";
    out << "VelocityL2Error: " <<  errors[4] << "\n";
    out << "MagneticL2Error: " <<  errors[5] << "\n";
    out << "hMax: " << mesh_ptr->h_max() << "\n";
    out << "NbCells: " << mesh_ptr->n_cells() << "\n";
    out << "NbInternalFaces: " << mesh_ptr->n_i_faces() << "\n";
    out << "MeshReg: " << mesh_ptr->regularity()[0] << "\n";
    // out << "DOFs: " << n_cells * n_local_cell_dofs + mesh_ptr->n_i_faces() * n_local_face_dofs << "\n";
    out.close();

    return 0;
}

bool program_options(int argc, const char **argv, std::string &mesh_name, size_t &u_id, size_t &b_id, size_t &p_id, double &visc, double &diff, char &u_bc_id, char &b_bc_id, size_t &L, size_t &K, std::string &plot_file, bool &use_threads, double &tol)
{
    namespace po = boost::program_options;

    // Program options
    po::options_description desc("Allowed options");

    desc.add_options()("help,h", "Produce help message")("mesh,m", po::value<std::string>(), "Set the mesh")("velocity,u", po::value<size_t>(), "Set the velocity")("magnetic,b", po::value<size_t>(), "Set the magnetic")("pressure,p", po::value<size_t>(), "Set the pressure")("viscosity,v", po::value<double>(), "Set the viscosity")("diffusivity,d", po::value<double>(), "Set the diffusivity")("u_bc", po::value<char>(), "Set the type of velocity boundary condition ('H' = Hodge, 'D' = Dirichlet')")("b_bc", po::value<char>(), "Set the type of magnetic boundary condition ('H' = Hodge, 'D' = Dirichlet')")("celldegree,l", po::value<size_t>(), "Set the degree of the cell polynomials")("facedegree,k", po::value<size_t>(), "Set the degree of the face polynomials")("plot", po::value<std::string>(), "Plot to file")("use_threads", po::value<bool>(), "Using multithreading")("tolerance", po::value<double>(), "The tolerance of the Newton iteration");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // Display help message
    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        return false;
    }

    // Get mesh file
    mesh_name = (vm.count("mesh") ? vm["mesh"].as<std::string>() : "Tetgen-Cube-0/RF_fmt/cube.2");

    // Get bc
    // bc_id = (vm.count("bc") ? vm["bc"].as<std::string>() : "D");

    u_id = (vm.count("velocity") ? vm["velocity"].as<size_t>() : 1);
    b_id = (vm.count("magnetic") ? vm["magnetic"].as<size_t>() : 3);
    p_id = (vm.count("pressure") ? vm["pressure"].as<size_t>() : 2);
    visc = (vm.count("viscosity") ? vm["viscosity"].as<double>() : 0.1);
    diff = (vm.count("diffusivity") ? vm["diffusivity"].as<double>() : 0.1);
    u_bc_id = (vm.count("u_bc") ? vm["u_bc"].as<char>() : 'D');
    b_bc_id = (vm.count("b_bc") ? vm["b_bc"].as<char>() : 'H');

    tol = (vm.count("tolerance") ? vm["tolerance"].as<double>() : 1E-6);

    // Get polynomial degrees
    L = (vm.count("celldegree") ? vm["celldegree"].as<size_t>() : 0);
    K = (vm.count("facedegree") ? vm["facedegree"].as<size_t>() : 0);

    // Check compatible face and cell degrees
    if ((std::abs(int(K) - int(L)) > 1) || (K < 0) || (L < 0))
    {
        std::cout << "Degrees k and l are not in acceptable range: k =" << K << ", l =" << L << "\n";
        return false;
    }

    // Get plot file
    plot_file = (vm.count("plot") ? vm["plot"].as<std::string>() : "");

    // Get use_threads
    use_threads = (vm.count("use_threads") ? vm["use_threads"].as<bool>() : true);

    return true;
}