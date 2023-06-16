#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/program_options.hpp>

/*!
 * \addtogroup TransformMeshes
 * @{
 */

/// Structure to store elements
struct Element
{
  /// Constructor
  Element(
    const size_t & _type,   ///< Type of element
    const size_t & _id,      ///< Id of the element
    const std::vector<size_t> & _vertices    ///< List of indices of vertices
    )
    : type(_type),
      id(_id),
      vertices(_vertices)
    {
      // Compute nb of vertices and faces, depending on element type
      switch (_type){
        case 4 :  // Tetrahedron
          n_vertices = 4;
          n_faces = 4;
          break;
        case 5 :  // Hexahedron
          n_vertices = 8;
          n_faces = 6;
        break;
        case 6 :  // Prism
          n_vertices = 6;
          n_faces = 5;
        break;
        case 7 :  // Pyramid
          n_vertices = 5;
          n_faces = 5;
        break;
        default :
         std::cout << "Element type " << type << " unknown" << std::endl;
         exit(1);
      }
    }

  /// Plots the element in RF format into the file provided
  void printRF(std::ofstream & file) const
  {
    // Header for the element
    file << id << " " << n_faces << std::endl;
    // Vector that lists the indices, in 'vertices', of the vertices forming the faces
    std::vector<std::vector<size_t>> faces;
    switch(type){
      case 4 : 
        faces = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
        break;
      case 5 :
        faces = {{0, 1, 2, 3}, {0, 1, 5, 4}, {0, 3, 7, 4}, {3, 2, 6, 7}, {2, 6, 5, 1}, {4, 5, 6, 7}};
        break;
      case 6 :
        faces = {{0, 1, 2}, {0, 1, 4, 3}, {0, 3, 5, 2}, {1, 2, 5, 4}, {3, 4, 5}};
        break;
      case 7 :
        faces = {{0, 1, 2, 3}, {0, 1, 4}, {0, 3, 4}, {1, 2, 4}, {2, 3, 4}};
        break;
      default :
        std::cout << "Element type " << type << " unknown" << std::endl;
        exit(1);
    }
    // Print out the faces
    for (size_t iF = 0; iF < faces.size(); iF++){
      file << iF << " " << faces[iF].size();
      for (size_t iV = 0; iV < faces[iF].size(); iV++){
        file << " " << vertices[faces[iF][iV]];
      }
      file << std::endl;
    }

  }

  const size_t type;
  const size_t id;
  const std::vector<size_t> vertices;
  size_t n_vertices;
  size_t n_faces;
};


//*** Default mesh
const std::string mesh_dir = "../../../meshes/";
std::string default_mesh = mesh_dir + "ComplexGeometries/CubeCavityWedge/MSH_fmt/hexa.1";

/// Main executable (MSHToRF) to convert a .msh mesh in RF format
int main(int argc, char *argv[])
{

  // Program options
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Converts a .msh file into RF format (two files: .node, .ele).")
    ("mesh,m", boost::program_options::value<std::string>(), "Input mesh file (complete path but without extension)")
    ("outfile,o", boost::program_options::value<std::string>(), "Output mesh file (without extension)");

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  // Display the help options
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  // Select the mesh
  std::string mesh_file = (vm.count("mesh") ? vm["mesh"].as<std::string>() : default_mesh) + ".msh";

  std::cout << "Mesh: " << mesh_file << std::endl;

  // Read .msh
  std::istringstream SS;    // To read each line as a stream
  std::ifstream mshfile;
  mshfile.open(mesh_file);
  if (!mshfile)
  {
      std::cerr << "     Unable to open .msh file (note: do not include the extension in the mesh file passed as parameter)\n";
      exit(1);
  }
  std::cout << "Opening " << mesh_file << std::endl;
  
  //------
  // Read the vertices
  //------
  // Look for nodes: ignore everything until $Nodes
  std::string line;
  while (std::getline(mshfile, line) && line.substr(0, 6) != "$Nodes")
  {
    // Just skip the line
  }

  // Nb of entities (useless), nodes, and min/max indices of nodes
  size_t ignore_int;
  size_t n_verts;
  size_t min_node;
  size_t max_node;
  std::getline(mshfile, line);
  SS.str(line);
  SS >> ignore_int >> n_verts >> min_node >> max_node;

  // Check that the vertices are labelled without any gap
  if (max_node + 1 != min_node + n_verts){
    std::cout << "Labels of vertices have gaps, stopping there." << std::endl;
    exit(1);
  }

  std::cout << "Reading vertices...";
  std::vector<std::vector<double>> vertex_coords(n_verts, std::vector<double>{0., 0., 0.});
  size_t count = 0;

  while (count < n_verts)
  {
    // Read one block
    size_t n_vert_block;
    std::getline(mshfile, line);
    SS.clear();
    SS.str(line);
    SS >> ignore_int >> ignore_int >> ignore_int >> n_vert_block;
    std::vector<size_t> idx_vertex(n_vert_block);
    
    // Read indices of vertices
    size_t idx;
    for (size_t i=0; i < n_vert_block; i++){
      std::getline(mshfile, line);
      SS.clear();
      SS.str(line);
      SS >> idx;
      idx_vertex[i] = idx - min_node;
    }
    // Read coordinates
    for (size_t i=0; i < n_vert_block; i++){
      std::getline(mshfile, line);
      SS.clear();
      SS.str(line);
      SS >> vertex_coords[idx_vertex[i]][0] >> vertex_coords[idx_vertex[i]][1] >> vertex_coords[idx_vertex[i]][2];
      count++;
    } 
  }
  // To check all vertices
//    for (size_t i=0; i<n_verts; i++){
//      std::cout << vertex_coords[i][0] << " " << vertex_coords[i][1] << " " << vertex_coords[i][2] << std::endl;
//    }
  std::cout << "...read " << count << " vertices" << std::endl;

  //------
  // Read the elements (in .msh, they can be lines, faces, cells)
  //------
  while (std::getline(mshfile, line) && line.substr(0, 9) != "$Elements")
  {
    // Just skip the line
  }

  // Nb of entities (useless), elements, and min/max indices of elements
  size_t n_eles;
  size_t min_ele;
  size_t max_ele;
  std::getline(mshfile, line);
  SS.clear();
  SS.str(line);
  SS >> ignore_int >> n_eles >> min_ele >> max_ele;

  std::cout << "Reading elements...";
  std::vector<Element> elements;  // list of elements (only cells at the moment, the rest is not useful)
  count = 0;  // counts all elements
  size_t n_cells = 0;  // counts only cells
  while (count < n_eles)
  {
    // Read one block
    size_t dim_entity;
    size_t ele_type;
    size_t n_eles_block;
    std::getline(mshfile, line);
    SS.clear();
    SS.str(line);
    SS >> dim_entity >> ignore_int >> ele_type >> n_eles_block;

    // Read the block: the number of vertices for each element depends on ele_type
    for (size_t i = 0; i < n_eles_block; i++){
      if (dim_entity != 3){
        // We actually skip this block as it's not a cell
        std::getline(mshfile, line);          
      }else{
        // Grab the vertices
        std::vector<size_t> vertices;
        std::getline(mshfile, line);
        SS.clear();
        SS.str(line);
        size_t idvertex;
        SS >> ignore_int; // We just drop the cell index, won't be using it
        while (SS >> idvertex) vertices.push_back(idvertex - min_node);   // add vertices one by one, offsetting to start at 0
        
        // Create element
        elements.emplace_back(ele_type, n_cells, vertices);

        n_cells++;
      }
      count++;
    }
  }
  // To check all elements
//    std::cout << std::endl;
//    for (size_t i=0; i<n_cells; i++){
//      std::cout << "Element " << elements[i].id << " (type " << elements[i].type << "): ";
//      for (size_t iV=0; iV < elements[i].vertices.size(); iV++){
//        std::cout << elements[i].vertices[iV] << " ";
//      }
//      std::cout << std::endl;
//    }

  std::cout << "...read " << n_cells << " cells" << std::endl;

  mshfile.close();

  //------
  // File to write RF
  //------
  const size_t extension = mesh_file.find_last_of("."); 
  const std::string original_meshfile = mesh_file.substr(0, extension).substr(mesh_file.find_last_of("/\\") + 1);  
  const std::string filename_core = (vm.count("outfile") ? vm["outfile"].as<std::string>() : original_meshfile);  
  const std::string filename_node = filename_core + ".node";
  const std::string filename_ele = filename_core + ".ele";

  // write .node file
  std::cout << "Writing .node file" << std::endl;
  std::ofstream outNode(filename_node.c_str());
  outNode << "# *.node file of 3D-mesh in REGN_FACE format" << std::endl;
  outNode << "# " << filename_node.c_str() << " created from " << original_meshfile + ".msh" << std::endl;
  outNode << n_verts << " " << "3  0  0" << std::endl;
  for (size_t iV = 0; iV < n_verts; iV++){
    outNode << std::setprecision(16) << iV;
    outNode << " ";
    outNode << std::setprecision(16) << vertex_coords[iV][0];
    outNode << " ";
    outNode << std::setprecision(16) << vertex_coords[iV][1];
    outNode << " ";
    outNode << std::setprecision(16) << vertex_coords[iV][2];
    outNode << std::endl;
  }
  outNode.close();
   
  // write .ele file
  std::cout << "Writing .ele file" << std::endl;
  std::ofstream outEle(filename_ele.c_str());
  outEle << "# *.ele file of 3D-mesh in REGN_FACE format" << std::endl;
  outEle << "# " << filename_ele.c_str() << " created from " << filename_core + ".msh" << std::endl;
  outEle << n_cells << " " << " 0" << std::endl;
  for (size_t iT = 0; iT < n_cells; iT++){
    elements[iT].printRF(outEle);
  }
  
  outEle.close();
  
  return 0;
}



