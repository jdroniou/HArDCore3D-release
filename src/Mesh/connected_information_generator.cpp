#include "connected_information_generator.hpp"
#include <set>
#include <algorithm>
#include <unordered_set>
using namespace HArDCore3D;

bool are_equal_unordered(const std::vector<size_t>& a, const std::vector<int>& b) {
    std::multiset<int> sa(a.begin(), a.end());
    std::multiset<int> sb(b.begin(), b.end());
    return sa == sb;
}

void dfs(int index, std::vector<std::set<size_t>>& v, std::vector<bool>& visited, std::unordered_set<int>& component) {
    visited[index] = true;
    
    // Add all elements from the current vector to the component set
    for (int element : v[index]) {
        component.insert(element);
    }

    // Explore the neighbors (vectors that share elements)
    for (size_t i = 0; i < v.size(); ++i) {
        if (!visited[i]) {
            // Check for common elements between v[index] and v[i]
            std::unordered_set<int> set1(v[index].begin(), v[index].end());
            std::unordered_set<int> set2(v[i].begin(), v[i].end());
            std::unordered_set<int> intersection;

            for (int e : set1) {
                if (set2.find(e) != set2.end()) {
                    intersection.insert(e);
                }
            }

            // If there is an intersection, they are connected
            if (!intersection.empty()) {
                dfs(i, v, visited, component);
            }
        }
    }
}
// ConnectedInformationGenerator::ConnectedInformationGenerator(const std::string& mesh_file, Mesh& mesh,std::vector<std::vector<Vertex *>>& fracture_pointer)
//     : mesh_file_(mesh_file), m_mesh(mesh),fracture_pointer_(fracture_pointer){}
ConnectedInformationGenerator::ConnectedInformationGenerator(const std::string& mesh_file, Mesh& mesh)
    : mesh_file_(mesh_file), m_mesh(mesh){}

void ConnectedInformationGenerator::generate() {

  MeshReaderFracture mesh_reader_fracture(mesh_file_,m_mesh);
  mesh_reader_fracture.read_frac_file(fracture_pointer);


    std::vector<std::vector<size_t>> fracture_faces;
    for(const auto& f:fracture_pointer){
      std::vector<size_t> face;
      for(const Vertex * v:f){
        face.push_back(v->global_index());
      }
      fracture_faces.push_back(face);
    }

  //printing of fracture faces in vertex format
  // std::cout<<"Fracture_Faces: {";
  // for (size_t i=0;i<fracture_pointer.size();++i){
  //   std::cout<<"{";
  //   for (size_t j=0;j<fracture_pointer[i].size();++j){
  //     std::cout<<fracture_pointer[i][j]->global_index()<<" ";
  //   }
  //   std::cout << "} ";
  // }
  // std::cout << "}\n"<< std::endl;

  //calling out global boundary vertex and storing them in Boundary_vertices
  std::vector<int> Boundary_vertices;
  for (const Vertex * node : m_mesh.get_b_vertices()) {
    int Node = node->global_index();
    // std::cout << "Vertex[" << Node << "]: {" << vertices[Node][0] << ", "
    //           << vertices[Node][1] << ", " << vertices[Node][2] << "}\n";
    Boundary_vertices.push_back(Node);
  }
  // std::cout << "Total number of Boundary vertices are: " << m_mesh.n_b_vertices() << "\n\n";


  //Calling out all fracture vertices and farcture vertices on the boundary
  std::unordered_set<int> boundary_set(Boundary_vertices.begin(), Boundary_vertices.end());
  std::set<int> unique_fracture_vertices;
  std::set<int> unique_non_boundary_fracture_vertices;  // use set to avoid duplicates
  std::vector<Vertex *>unique_fracture_vertices_pointer;
  std::vector<Vertex *>unique_non_boundary_fracture_vertices_pointer;

    //Computing Fracture vertices
    for (const auto& face : fracture_pointer) {
        for (const Vertex * vertex : face) {
                unique_fracture_vertices.insert(vertex->global_index());
        }
    }
    //Computing non-boundary Fracture vertices
    for (const auto& face : fracture_pointer) {
        for (const Vertex * vertex : face) {
            if (boundary_set.find(vertex->global_index()) == boundary_set.end()) {
                unique_non_boundary_fracture_vertices.insert(vertex->global_index());
            }
        }
    }

    for(int v:unique_fracture_vertices){
      unique_fracture_vertices_pointer.push_back(m_mesh.vertex(v));
    }
    for(int v:unique_non_boundary_fracture_vertices){
      unique_non_boundary_fracture_vertices_pointer.push_back(m_mesh.vertex(v));
    }
    
    std::vector<int> fracture_vertex(unique_fracture_vertices.begin(), unique_fracture_vertices.end());
    fracture_vertex_pointer.assign(unique_fracture_vertices_pointer.begin(), unique_fracture_vertices_pointer.end());
    // std::cout << "Fracture Vertex: {";    
    // for (const auto& v : fracture_vertex_pointer) {
    //     std::cout << v->global_index() << " ";
    // }
    // std::cout << "}\n" << std::endl;
    
    std::vector<int> fracture_non_boundary_vertex(unique_non_boundary_fracture_vertices.begin(), unique_non_boundary_fracture_vertices.end());
    fracture_non_boundary_vertex_pointer.assign(unique_non_boundary_fracture_vertices_pointer.begin(), unique_non_boundary_fracture_vertices_pointer.end());
    // std::cout << "Fracture Non-Boundary Vertex: {";
    // for (const auto& v : fracture_non_boundary_vertex_pointer) {
    //     std::cout << v->global_index() << " ";
    // }
    // std::cout << "}\n" << std::endl;



    // Storing which faces are fracture faces and which are not
    std::vector<int> fracture_face_indices;

    for (const auto& frac_face : fracture_faces) {
        // bool found = false;
        for (Face* node : m_mesh.get_faces()) {
          int faces = node->global_index();
          // const Face& face_num = *node;
          const auto vert = node->get_vertices();
          std::vector<int> vec(vert.size());
          for (size_t i = 0; i < vert.size(); ++i) {
            vec[i] = vert[i]->global_index();
          }
            if (are_equal_unordered(frac_face, vec)) {
                fracture_face_indices.push_back(faces);
                fracture_face_indices_pointer.push_back(node);
                break;
            }
        }
    }

    std::sort(fracture_face_indices.begin(), fracture_face_indices.end());

    std::vector<int> non_fracture_face_indices;
    for (size_t i = 0; i < m_mesh.n_faces(); ++i) {
        if (std::find(fracture_face_indices.begin(), fracture_face_indices.end(), i) == fracture_face_indices.end()) {
          non_fracture_face_indices.push_back(i);
          non_fracture_face_indices_pointer.push_back(m_mesh.face(i));
        }
    }


    // Print Fracture Faces and Non-Fracture Faces indices
    // std::cout << "Fracture Face Indices: {";
    // for (const auto&f : fracture_face_indices_pointer) {
    //     std::cout << f->global_index() << " ";
    // }
    // std::cout << "}\n" << std::endl;

    // std::cout << "Non-Fracture Face indices: {";
    // for (const auto&f : non_fracture_face_indices_pointer) {     
    //     std::cout << f->global_index() << " ";
    // }
    // std::cout << "}\n" << std::endl;


    // Storing which edges are fracture edges and which are not
    std::vector<int> fracture_edge_indices;
    std::vector<int> fracture_non_boundary_edge_indices;
    std::vector<int> non_fracture_edge_indices;
    std::vector<int> non_fracture_and_boundary_edge_indices;

    // bool found = false;
    for (const auto& node : m_mesh.get_edges()) {
      int edges = node->global_index();
      Edge edge_num = *node;
      auto vert = edge_num.get_vertices();
      int vert0 = vert[0]->global_index();
      int vert1 = vert[1]->global_index();
      std::vector<int> vec = {vert0, vert1};
        for (int v : vec){
            if (std::find(fracture_vertex.begin(), fracture_vertex.end(), v) == fracture_vertex.end()) {
                non_fracture_edge_indices.push_back(edges);
                break;
            }
        }
        if (std::find(non_fracture_edge_indices.begin(), non_fracture_edge_indices.end(), edges) == non_fracture_edge_indices.end()) {
          fracture_edge_indices.push_back(edges);
          for (int v : vec){
            if (std::find(fracture_non_boundary_vertex.begin(), fracture_non_boundary_vertex.end(), v) != fracture_non_boundary_vertex.end()) {
                fracture_non_boundary_edge_indices.push_back(edges);
                break;
            }
          }
        }
    }

    for (size_t i=0;i<m_mesh.n_edges();++i){
      if (std::find(fracture_non_boundary_edge_indices.begin(),fracture_non_boundary_edge_indices.end(),i)==fracture_non_boundary_edge_indices.end()){
        non_fracture_and_boundary_edge_indices.push_back(i);
      }
    }

    for(int i:fracture_edge_indices){
      fracture_edge_indices_pointer.push_back(m_mesh.edge(i));
    }
    for(int i:fracture_non_boundary_edge_indices){
      fracture_non_boundary_edge_indices_pointer.push_back(m_mesh.edge(i));
    }
    for(int i:non_fracture_edge_indices){
      non_fracture_edge_indices_pointer.push_back(m_mesh.edge(i));
    }
    for(int i:non_fracture_and_boundary_edge_indices){
      non_fracture_and_boundary_edge_indices_pointer.push_back(m_mesh.edge(i));
    }

  // Print Fracture Edges and Non-Fracture Edges indices
  // std::cout << "Fracture Edge Indices: {";
  // for (const auto& idx : fracture_edge_indices_pointer) {
  //     std::cout << idx->global_index() << " ";
  // }
  // std::cout << "}\n" << std::endl;

  // std::cout << "Non-Fracture Edge indices: {";
  // for (const auto& idx : non_fracture_edge_indices_pointer) {
  //     std::cout << idx->global_index() << " ";
  // }
  // std::cout << "}\n" << std::endl;

  // std::cout << "Non-Boundary Fracture Edge Indices: {";
  // for (const auto& idx : fracture_non_boundary_edge_indices_pointer) {
  //     std::cout << idx->global_index() << " ";
  // }
  // std::cout << "}\n" << std::endl;

  // std::cout << "Non-Fracture and Boundary Edge indices: {";
  // for (const auto& idx : non_fracture_and_boundary_edge_indices_pointer) {
  //     std::cout << idx->global_index() << " ";
  // }
  // std::cout << "}\n" << std::endl;


  vertex_first_connected.resize(m_mesh.n_vertices());
  vertex_final_connected.resize(m_mesh.n_vertices());
  edge_first_connected.resize(m_mesh.n_edges());
  edge_final_connected.resize(m_mesh.n_edges());
  face_first_connected.resize(m_mesh.n_faces());
  face_final_connected.resize(m_mesh.n_faces());
  // face_final_connected_pointer.resize(m_mesh.n_faces());
  for (const auto* cell_node :m_mesh.get_cells()){
    //looking for the connected components with respect to vertices of the mesh.
    for (const auto* v:cell_node->get_vertices()){
      int vert=v->global_index();
      if (std::find(fracture_vertex.begin(),fracture_vertex.end(),vert) == fracture_vertex.end()){
        std::set<size_t> b;
        for (size_t i=0;i<v->n_cells();++i){
          b.insert(v->get_cells()[i]->global_index());
        }
        vertex_first_connected[vert].push_back(b);
      }
      if (std::find(fracture_vertex.begin(),fracture_vertex.end(),vert) != fracture_vertex.end()){
        std::set<size_t> b;
        b.insert(cell_node->global_index());
        for(const auto* f:cell_node->get_faces()){
          int my_face=f->global_index();
          // std::cout<<my_face;
          if (std::find(non_fracture_face_indices.begin(), non_fracture_face_indices.end(), my_face) != non_fracture_face_indices.end()) {
            std::vector<size_t> my_face_vertices;
            for (const auto* mfv:f->get_vertices()){
              my_face_vertices.push_back(mfv->global_index());
            }
            if(std::find(my_face_vertices.begin(),my_face_vertices.end(),vert)!=my_face_vertices.end()){
              for(const auto* cell:f->get_cells()){
                b.insert(cell->global_index());
              }
            }
          }
        }
        vertex_first_connected[vert].push_back(b);
      }
    }
    
    //looking for connected components with respect to edges of the mesh
    for (const auto* e:cell_node->get_edges()){
      int edg=e->global_index();
      if (std::find(fracture_edge_indices.begin(),fracture_edge_indices.end(),edg) == fracture_edge_indices.end()){
        std::set<size_t> b;
        for (size_t i=0;i<e->n_cells();++i){
          b.insert(e->get_cells()[i]->global_index());
        }
        edge_first_connected[edg].push_back(b);
      }
      if (std::find(fracture_edge_indices.begin(),fracture_edge_indices.end(),edg) != fracture_edge_indices.end()){
        std::set<size_t> b;
        b.insert(cell_node->global_index());
        for(const auto* f:cell_node->get_faces()){
          int my_face=f->global_index();
          // std::cout<<my_face;
          if (std::find(non_fracture_face_indices.begin(), non_fracture_face_indices.end(), my_face) != non_fracture_face_indices.end()) {
            std::vector<size_t> my_face_edges;
            for (const auto* mfe:f->get_edges()){
              my_face_edges.push_back(mfe->global_index());
            }
            if(std::find(my_face_edges.begin(),my_face_edges.end(),edg)!=my_face_edges.end()){
              for(const auto* cell:f->get_cells()){
                b.insert(cell->global_index());
              }
            }
          }
        }
        edge_first_connected[edg].push_back(b);
      }
    }
    for (const auto* f:cell_node->get_faces()){
      int fac=f->global_index();
      std::vector<size_t> b;
      std::vector<Cell *> b1;
      for (size_t i=0;i<f->n_cells();++i){
        b.push_back(f->get_cells()[i]->global_index());
        b1.push_back(f->get_cells()[i]);
      }
      std::vector<std::vector<size_t>> c;
      std::vector<std::vector<Cell *>> c1;
      bool is_fracture = std::find(fracture_face_indices.begin(), fracture_face_indices.end(), fac) != fracture_face_indices.end();

      if (!is_fracture) {
        c.resize(1);
        c1.resize(1);
        c[0].assign(b.begin(), b.end());
        c1[0].assign(b1.begin(), b1.end());
      } 
      else {
        c.resize(2);
        c1.resize(2);
        c[0].push_back(b[0]);
        c[1].push_back(b[1]);
        c1[0].push_back(b1[0]);
        c1[1].push_back(b1[1]);
    }

    face_final_connected[fac]=c;
    face_final_connected_size[m_mesh.face(fac)]=face_final_connected[fac].size();
    face_final_connected_pointer[m_mesh.face(fac)]=c1;
    }
  }

  for(const auto* v:m_mesh.get_vertices()){
    int vertID=v->global_index();
    int n = vertex_first_connected[vertID].size();
        std::vector<bool> visited(n, false);
        std::vector<std::unordered_set<int>> connectedComponents;

        // Find all connected components
        for (int i = 0; i < n; ++i) {
            if (!visited[i]) {
                std::unordered_set<int> component;
                dfs(i, vertex_first_connected[vertID], visited, component);
                connectedComponents.push_back(component);
            }
        }

        // Store each component as a set<int> in w[v]
        for (const auto& component : connectedComponents) {
            std::set<int> sortedComponent(component.begin(), component.end());
            vertex_final_connected[vertID].push_back(sortedComponent);
        }
        std::vector<std::vector<Cell *>> c;
        vertex_final_connected_size[m_mesh.vertex(vertID)]=vertex_final_connected[vertID].size();
        for(size_t i=0;i<vertex_final_connected[vertID].size();++i){
          std::vector<Cell *> b;
          for(int j:vertex_final_connected[vertID][i]){
            b.push_back(m_mesh.cell(j));
          }
          c.push_back(b);
        }
        vertex_final_connected_pointer[m_mesh.vertex(vertID)]=c;
  }

  for(const auto* e:m_mesh.get_edges()){
    int edgeID=e->global_index();
    int n = edge_first_connected[edgeID].size();
        std::vector<bool> visited(n, false);
        std::vector<std::unordered_set<int>> connectedComponents;

        // Find all connected components
        for (int i = 0; i < n; ++i) {
            if (!visited[i]) {
                std::unordered_set<int> component;
                dfs(i, edge_first_connected[edgeID], visited, component);
                connectedComponents.push_back(component);
            }
        }

        // Store each component as a set<int> in w[v]
        for (const auto& component : connectedComponents) {
            std::set<int> sortedComponent(component.begin(), component.end());
            edge_final_connected[edgeID].push_back(sortedComponent);
        }
        std::vector<std::vector<Cell *>> c;
        edge_final_connected_size[m_mesh.edge(edgeID)]=edge_final_connected[edgeID].size();
        for(size_t i=0;i<edge_final_connected[edgeID].size();++i){
          std::vector<Cell *> b;
          for(int j:edge_final_connected[edgeID][i]){
            b.push_back(m_mesh.cell(j));
          }
          c.push_back(b);
        }
        edge_final_connected_pointer[m_mesh.edge(edgeID)]=c;
  }


  // for (const auto& v :m_mesh.get_vertices()) {
  //   std::cout << "Vertex_final_connected[" << v->global_index() << "]: ";
  //   std::cout << " { ";
  //   for (const auto& inner_vec : vertex_final_connected_pointer[v]) {
  //       std::cout << "{ ";
  //       for (const auto& val : inner_vec) {
  //           std::cout << val->global_index() << " ";
  //       }
  //       std::cout << "} ";
  //   } 
  //   std::cout << "}\n";
  //   std::cout<<"Vertex conected size: "<<vertex_final_connected_size[v]<<"\n";
  // }

  // for (const auto& e :m_mesh.get_edges()) {
  //   std::cout << "Edge_final_connected[" << e->global_index() << "]: ";
  //   std::cout << " { ";
  //   for (const auto& inner_vec : edge_final_connected_pointer[e]) {
  //       std::cout << "{ ";
  //       for (const auto& val : inner_vec) {
  //           std::cout << val->global_index() << " ";
  //       }
  //       std::cout << "} ";
  //   } 
  //   std::cout << "}\n";
  //   std::cout<<"Edge conected size: "<<edge_final_connected_size[e]<<"\n";
  // }


  // for (const auto& f :m_mesh.get_faces()) {
  //   std::cout << "Face_final_connected[" << f->global_index() << "]: ";
  //   std::cout << " { ";
  //   for (const auto& inner_vec : face_final_connected_pointer[f]) {
  //       std::cout << "{ ";
  //       for (const auto& val : inner_vec) {
  //           std::cout << val->global_index() << " ";
  //       }
  //       std::cout << "} ";
  //   } 
  //   std::cout << "}\n";
  //   std::cout<<"Face conected size: "<<face_final_connected_size[f]<<"\n";
  // }







  // for (const auto& node : m_mesh.get_vertices()) {
  //   int Node = node->global_index();
    // std::cout << "Vertex[" << Node << "]: {" << vertices[Node][0] << ", "
    //           << vertices[Node][1] << ", " << vertices[Node][2] << "}\n";
    // Vertex vert_num = *node;
    // auto cel = vert_num.get_cells();
    // std::vector<int> vec(cel.size());
    // for (int i = 0; i < cel.size(); ++i) {
    //   vec[i] = cel[i]->global_index();
    // }
    // std::cout << "Vertex_Cell[" << Node << "]: {";
    // for (size_t i = 0; i < vec.size(); ++i) {
    //   std::cout << vec[i];
    //   if (i != vec.size() - 1) std::cout << ", ";
    // }
    // std::cout << "}\n";
  // }
  
  // std::cout<<"Edge Vertex: {";
  // for (const auto& node : m_mesh.get_edges()) {
  //   int edges = node->global_index();
  //   Edge edge_num = *node;
  //   auto vert = edge_num.get_vertices();
  //   int vert0 = vert[0]->global_index();
  //   int vert1 = vert[1]->global_index();
    // std::cout << "Edge[" << edges << "]: {" << vert0 << "," << vert1 << "}\n";
    // std::cout << " {" << vert0 << "," << vert1 << "},";
  // }
  // std::cout << "}"<<std::endl;
  // std::cout << "Total number of Edges are: " << m_mesh.n_edges() << "\n\n";

  // std::cout<<"Face Vertex: {";
  // for (const auto& node : m_mesh.get_faces()) {
  //   int faces = node->global_index();
  //   Face face_num = *node;
  //   auto vert = face_num.get_vertices();
  //   std::vector<int> vec(vert.size());
  //   for (int i = 0; i < vert.size(); ++i) {
  //     vec[i] = vert[i]->global_index();
  //   }
    // std::cout << "Face[" << faces << "]: {";
    // std::cout << " {";
    // for (size_t i = 0; i < vec.size(); ++i) {
    //   std::cout << vec[i];
    //   if (i != vec.size() - 1) std::cout << ", ";
    // }
    // std::cout << "},";
  // }
  // std::cout << "}"<<std::endl;
  // std::cout << "Total number of Faces are: " << m_mesh.n_faces() << "\n\n";

  // std::cout<<"Cell Vertex: {";
  // for (const auto& node : m_mesh.get_cells()) {
  //   int cells = node->global_index();
  //   Cell cell_num = *node;
  //   auto vert = cell_num.get_vertices();
  //   std::vector<int> vec(vert.size());
  //   for (int i = 0; i < vert.size(); ++i) {
  //     vec[i] = vert[i]->global_index();
  //   }
    // std::cout << "Cell[" << cells << "]: {";
    // std::cout << " {";
    // for (size_t i = 0; i < vec.size(); ++i) {
    //   std::cout << vec[i];
    //   if (i != vec.size() - 1) std::cout << ", ";
    // }
    // std::cout << "},";
  // }
  // std::cout << "}"<<std::endl;
  // std::cout << "Total number of Cells are: " << m_mesh.n_cells() << "\n\n";

  //Storing the number of fracture faces in each cell and the corresponding fracture faces pointers
  for(const Cell*T:m_mesh.get_cells()){
    int cellID=T->global_index();
    size_t count=0;
    std::vector<Face *> face_pointers;
    for(Face*F:T->get_faces()){
      const size_t FaceID=F->global_index();
      for(const Face *fF : fracture_face_indices_pointer){
        const size_t currentface=fF->global_index();
        if(currentface==FaceID){
          count+=1;
          face_pointers.push_back(F);
          break;
        }
      }
    }
    std::pair<size_t, std::vector<Face *>> current_pair={count, face_pointers};
    Cell_fractures[m_mesh.cell(cellID)]=current_pair;
  }

  // Checking on the count of fracture faces and their respective pointers.
  // for(Cell*T:m_mesh.get_cells()){
  //   size_t count=Cell_fractures.at(T).second.size();
  //   for(size_t i=0;i<count;++i){
  //     std::cout<<T->global_index()<<" "<<Cell_fractures.at(T).second[i]<<" ";
  //   }
  //   std::cout<<"\n";
  // }
}


int ConnectedInformationGenerator::vertex_component(size_t iV, size_t iT)const
  {
    const auto& vert_ptr=m_mesh.vertex(iV);
    int connected_size_vert=vertex_final_connected_size.at(vert_ptr);
    int num=0;
    bool found = false;
    for (int i = 0; i < connected_size_vert && !found; ++i) {
      const std::vector<Cell *> &current_component = vertex_final_connected_pointer.at(vert_ptr)[i];
      for (size_t j = 0; j < current_component.size(); ++j) {
            if (current_component[j]->global_index() == iT) {
                num=i;
                found = true;
                break;
            }
        }
    }
    return num;
  }

  int ConnectedInformationGenerator::edge_component(size_t iE, size_t iT)const
  {
    const auto& edge_ptr=m_mesh.edge(iE);
    int connected_size_edge=edge_final_connected_size.at(edge_ptr);
    int num=0;
    bool found = false;
    for (int i = 0; i < connected_size_edge && !found; ++i) {
      const std::vector<Cell *> &current_component = edge_final_connected_pointer.at(edge_ptr)[i];
      for (size_t j = 0; j < current_component.size(); ++j) {
            if (current_component[j]->global_index() == iT) {
                num=i;
                found = true;
                break;
            }
        }
    }
    return num;
  }

  int ConnectedInformationGenerator::face_component(size_t iF, size_t iT)const
  {
    const auto& face_ptr=m_mesh.face(iF);
    int connected_size_face=face_final_connected_size.at(face_ptr);
    int num=0;
    bool found = false;
    for (int i = 0; i < connected_size_face && !found; ++i) {
      const std::vector<Cell *> &current_component = face_final_connected_pointer.at(face_ptr)[i];
      for (size_t j = 0; j < current_component.size(); ++j) {
            if (current_component[j]->global_index() == iT) {
                num=i;
                found = true;
                break;
            }
        }
    }
    return num;
  }

  int ConnectedInformationGenerator::findFaceIndex(const Face& F)const
  {
    const Face* ptr = &F;
    // Use std::find to find the pointer in the vector
    auto it = std::find(fracture_face_indices_pointer.begin(), fracture_face_indices_pointer.end(), ptr);

    if (it != fracture_face_indices_pointer.end()) {
        // Found, return the index
        return static_cast<int>(std::distance(fracture_face_indices_pointer.begin(), it));
    } else {
        // Not found
        return -1;
    }
  }