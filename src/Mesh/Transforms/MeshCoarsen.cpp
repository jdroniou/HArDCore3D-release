#include "DirectedGraph.hpp"
#include <string>
#include <boost/timer/timer.hpp>

/*!
 * \addtogroup TransformMeshes
 * @{
 */


//const std::string mesh_root = "../../../meshes/";
const std::string mesh_root = "";

/// Main executable (MeshCoarsen) to coarsen a mesh
int main(int argc, char *argv[])
{
    if (argc == 1){
      std::cout << "Coarsen a mesh.\n Usage: mesh-coarsen <name of mesh> <number of levels of coarsening>" << std::endl;
      exit(0);
    }
    std::string mesh_file = argv[1];
    std::size_t iterations = std::stoi(argv[2]);

    std::cout.unsetf(std::ios::floatfield);
    std::cout.precision(4);

    //  Open node file

    std::ifstream inNode;
    inNode.open(mesh_root + mesh_file + ".node");
    if (!inNode)
    {
        std::cerr << "Unable to open node file\n";
        exit(1);
    }

    std::cout << "Reading node file...\n";

    //  Read mesh file
    std::string line;
    std::size_t ignore_int;
    std::size_t n_verts;
    std::vector<Eigen::VectorXd> vertices;

    bool ignore_line = true;

    while (ignore_line)
    {
        std::getline(inNode, line);
        if (!(line.substr(0, 1) == "#"))
        {
            std::string delimiter = " ";
            std::string string_n_verts = line.substr(0, line.find(delimiter));
            n_verts = std::stoi(string_n_verts);
            ignore_line = false;
        }
    }

    Eigen::VectorXd coord(3);
    std::size_t count = 0;
    double xcoord, ycoord, zcoord;

    while (count < n_verts)
    {
        inNode >> ignore_int;
        inNode >> xcoord >> ycoord >> zcoord;
        coord << xcoord, ycoord, zcoord;
        vertices.push_back(coord);
        count++;
    }
    inNode.close();

    std::cout << "Read " + std::to_string(vertices.size()) + "/" + std::to_string(n_verts) + " vertices\n";

    std::ifstream inElement;
    inElement.open(mesh_root + mesh_file + ".ele");
    if (!inElement)
    {
        std::cerr << "Unable to open element file\n";
        exit(1);
    }

    std::size_t n_cells;

    ignore_line = true;

    while (ignore_line)
    {
        std::getline(inElement, line);
        if (!(line.substr(0, 1) == "#"))
        {
            std::string delimiter = " ";
            std::string string_n_cells = line.substr(0, line.find(delimiter));
            n_cells = std::stoi(string_n_cells);
            ignore_line = false;
        }
    }

    DirectedGraph graph;

    std::size_t cell_count = 0;
    std::size_t n_cell_faces;
    std::size_t n_face_verts;

    while (cell_count < n_cells)
    {
        inElement >> ignore_int;
        inElement >> n_cell_faces;
        std::size_t face_count = 0;
        std::vector<DirectedFace> faces;
        // std::cout << n_cell_faces << std::endl;
        while (face_count < n_cell_faces)
        {
            inElement >> ignore_int;
            inElement >> n_face_verts;
            // std::cout << n_face_verts << " ";
            std::size_t vert_count = 0;
            std::vector<size_t> face_verts;
            std::size_t vert;
            while (vert_count < n_face_verts)
            {
                inElement >> vert;
                face_verts.push_back(vert);
                vert_count++;
            }
            faces.push_back(DirectedFace(face_verts));
            face_count++;
        }
        // std::cout << "\n";
        DirectedCell cell(faces, cell_count);
        graph.add_cell(cell);
        cell_count++;
    }

    std::cout << "Read " + std::to_string(graph.G.size()) + "/" + std::to_string(n_cells) + " cells\n";

    inElement.close();

    // if (!graph.test_graph())
    // {
    //     exit(1);
    // }

    std::cout << "\nCoarsening graph...\n";
    graph.randomise();
    for (size_t i = 0; i < iterations; i++)
    {
        boost::timer::cpu_timer timer;
        timer.start();
        graph.coarsen();
        double time = double(timer.elapsed().wall) * pow(10, -9);
        std::cout << std::setprecision(4) << "Iteration " + std::to_string(i + 1) + " completed in " << time << "s. " + std::to_string(graph.G.size()) + " cells remaining\n";
        graph.order();
    }
    std::cout << "Coarsening complete\n";
    assert(graph.test_graph());

    // if (!graph.test_graph())
    // {
    //     exit(1);
    // }

    // std::ofstream plotout(mesh_root + "aglomerated/" + mesh_file + ".coarse." + 			std::to_string(iterations) + ".dat");
    // graph.plotfile(&plotout, vertices);
    // plotout.close();

    // std::string partition = graph.get_partition();

    GraphNodeArray graph_node_array = graph.graph_to_array();

    // Remove unused vertices
    for (size_t iV = n_verts - 1; iV > 0; --iV)
    {
        if (!(graph_node_array.node_exists(iV)))
        {
            vertices.erase(vertices.begin() + iV);
            graph_node_array.renum_nodes(iV);
            n_verts--;
        }
    }

    std::cout << "\nWriting mesh to file\n";

    std::ofstream outNode(mesh_root + mesh_file + ".coarse." + std::to_string(iterations) + ".node");

    // std::string find = "RF_fmt";
    // std::string file_no_folders = line.substr(mesh_file.find(find) + find.length(), mesh_file.length() - (mesh_file.find(find) + find.length()));

    outNode << "# *.node file of 3D-mesh in REGN_FACE format\n";
    outNode << "# ./" + mesh_file + ".coarse." + std::to_string(iterations) + ".node\n";
    outNode << n_verts << "  3  0  0\n";
    for (size_t iV = 0; iV < n_verts; iV++)
    {
        outNode << std::setprecision(16) << iV;
        outNode << " ";
        outNode << std::setprecision(16) << vertices[iV](0);
        outNode << " ";
        outNode << std::setprecision(16) << vertices[iV](1);
        outNode << " ";
        outNode << std::setprecision(16) << vertices[iV](2);
        outNode << "\n";
    }
    outNode << "# output from MeshCoarsen.cpp";
    outNode.close();

    std::ofstream outElement(mesh_root + mesh_file + ".coarse." + std::to_string(iterations) + ".ele");
    outElement << "# *.ele file of 3D-mesh in REGN_FACE format \n";
    outElement << "# ./" + mesh_file + ".coarse." + std::to_string(iterations) + ".ele\n";
    graph_node_array.print(&outElement);
    outElement << "# output from MeshCoarsen.cpp";
    outElement.close();

    std::ofstream outPart(mesh_root + mesh_file + ".coarse." + std::to_string(iterations) + ".partition");
    std::string partition = graph.get_partition();
    outPart << partition;
    outPart.close();

    return 0;
}
