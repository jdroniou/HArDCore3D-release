#include <mesh_reader_fracture.hpp>

using namespace HArDCore3D;

MeshReaderFracture::MeshReaderFracture(std::string file_name, Mesh& mesh)
    : _file_name(file_name), m_mesh(mesh) {}

void MeshReaderFracture::read_frac_file(std::vector<std::vector<Vertex *>>& fracture_pointer) 
{
    std::vector<std::vector<size_t>> fracture_faces;
    std::ifstream inFrac;
    inFrac.open(_file_name + ".frac");
    if (!inFrac)
    {
        std::cerr << "     Unable to open fracture file\n";
        exit(1);
    }

    // std::cout << "     Reading element file...\n";

    std::string line;
    std::size_t ignore_int;
    std::size_t n_fracface;

    bool ignore_line = true;

    while (ignore_line)
    {
        std::getline(inFrac, line);
        if (!(line.substr(0, 1) == "#"))
        {
            std::string delimiter = " ";
            std::string string_n_fracture = line.substr(0, line.find(delimiter));
            n_fracface = std::stoi(string_n_fracture);
            ignore_line = false;
        }
    }
    if(n_fracface!=0){
        std::size_t frac_count = 0;
        std::size_t n_frac_vert;
        // std::size_t n_face_verts;

        while (frac_count < n_fracface)
        {
            inFrac >> ignore_int;
            inFrac >> n_frac_vert;
            std::size_t vert_count = 0;
            std::vector<size_t> face;
            double vert;
            while (vert_count < n_frac_vert)
            {
                inFrac >> vert;
                face.push_back(vert);
                ++vert_count;
            }
            fracture_faces.push_back(face);
            ++frac_count;
        }
    }
    inFrac.close();
    
    for(const auto& f:fracture_faces){
        std::vector<Vertex *> face_pointer_;
        for(const auto& v:f){
            face_pointer_.push_back(m_mesh.vertex(v));
        }
        fracture_pointer.push_back(face_pointer_);
    }

}

