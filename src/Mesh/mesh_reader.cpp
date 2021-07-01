#include <mesh_reader.hpp>

using namespace HArDCore3D;

MeshReaderRF::MeshReaderRF(std::string file_name)
    : _file_name(file_name) {}

void MeshReaderRF::read_node_file(std::vector<std::vector<double>> &vertices)
{
    std::ifstream inNode;
    inNode.open(_file_name + ".node");
    if (!inNode)
    {
        std::cerr << "     Unable to open node file\n";
        exit(1);
    }

    // std::cout << "     Reading node file...\n";

    //  Read mesh file
    std::string line;
    std::size_t ignore_int;
    std::size_t n_verts;

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

    std::size_t count = 0;
    double xcoord, ycoord, zcoord;

    while (count < n_verts)
    {
        inNode >> ignore_int; //node number = count
        inNode >> xcoord >> ycoord >> zcoord;
        vertices.push_back(std::vector<double>{xcoord, ycoord, zcoord});
        count++;
    }
    inNode.close();
    // std::cout << "     Read " + std::to_string(vertices.size()) + "/" + std::to_string(n_verts) + " vertices\n";
}

// cells = vector of cell (std::vector<std::vector<size_t>>)
// cell = vector of face (std::vector<size_t>)
// face = vector of vert (size_t)

void MeshReaderRF::read_ele_file(std::vector<std::vector<std::vector<size_t>>> &cells) 
{
    std::ifstream inEle;
    inEle.open(_file_name + ".ele");
    if (!inEle)
    {
        std::cerr << "     Unable to open element file\n";
        exit(1);
    }

    // std::cout << "     Reading element file...\n";

    std::string line;
    std::size_t ignore_int;
    std::size_t n_cells;

    bool ignore_line = true;

    while (ignore_line)
    {
        std::getline(inEle, line);
        if (!(line.substr(0, 1) == "#"))
        {
            std::string delimiter = " ";
            std::string string_n_cells = line.substr(0, line.find(delimiter));
            n_cells = std::stoi(string_n_cells);
            ignore_line = false;
        }
    }

    std::size_t cell_count = 0;
    std::size_t n_cell_faces;
    std::size_t n_face_verts;

    while (cell_count < n_cells)
    {
        inEle >> ignore_int;
        inEle >> n_cell_faces;
        std::size_t face_count = 0;
        std::vector<std::vector<size_t>> cell;
        while (face_count < n_cell_faces)
        {
            inEle >> ignore_int;
            inEle >> n_face_verts;

            std::size_t vert_count = 0;
            std::vector<size_t> face;
            std::size_t vert;
            while (vert_count < n_face_verts)
            {
                inEle >> vert;
                face.push_back(vert);
                ++vert_count;
            }
            cell.push_back(face);
            ++face_count;
        }
        cells.push_back(cell);
        ++cell_count;
    }
    inEle.close();
    // std::cout << "     Read " + std::to_string(cells.size()) + "/" + std::to_string(n_cells) + " cells\n";
}

void MeshReaderRF::read_partition_file(std::vector<std::vector<std::size_t>> &partition)
{
    std::ifstream inPartition;
    inPartition.open(_file_name + ".partition");
    if (!inPartition)
    {
        std::cerr << "     Unable to open partition file\n";
        exit(1);
    }

    // std::cout << "     Reading partition file...\n";

    std::size_t n_parts;
    std::size_t part;
    while (inPartition >> n_parts)
    {
        std::size_t part_count = 0;
        std::vector<std::size_t> cell_partition;
        while(part_count < n_parts)
        {
            inPartition >> part;
            cell_partition.push_back(part);
            ++part_count;
        }
        partition.push_back(cell_partition);
    }
    inPartition.close();
    // std::cout << "     Read partition file.\n";
}