#ifndef CONNECTED_INFORMATION_GENERATOR_HPP
#define CONNECTED_INFORMATION_GENERATOR_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <mesh.hpp>
#include <mesh_builder.hpp>
#include "mesh_reader_fracture.hpp"
#include<set>
// #include "mesh_data.hpp"

namespace HArDCore3D
{

class ConnectedInformationGenerator {
public:
    // ConnectedInformationGenerator(const std::string& mesh_file, Mesh& mesh, std::vector<std::vector<Vertex *>>& fracture_pointer);
    ConnectedInformationGenerator(const std::string& mesh_file, Mesh& mesh);

    std::vector<std::vector<std::set<size_t>>> vertex_first_connected;
    std::vector<std::vector<std::set<int>>> vertex_final_connected;
    std::vector<std::vector<std::set<size_t>>> edge_first_connected;
    std::vector<std::vector<std::set<int>>> edge_final_connected;
    std::vector<std::vector<std::set<size_t>>> face_first_connected;
    std::vector<std::vector<std::vector<size_t>>> face_final_connected;

    std::vector<std::vector<Vertex *>> fracture_pointer;    //Storing Fracture faces in pointer format of vertex
    std::vector<Face *> fracture_face_indices_pointer;      //Storing which face indices are fracture face in pointer format
    std::vector<Face *> non_fracture_face_indices_pointer;  //Storing which face indices are non-fracture face in pointer format

    std::vector<Vertex *>fracture_vertex_pointer;           //Storing which vertex indices are fracture vertex in pointer format
    std::vector<Vertex *>fracture_non_boundary_vertex_pointer;  //Storing which vertex indices are non_boudary fracture vertex in pointer format
    std::vector<Edge *> fracture_edge_indices_pointer;      //Storing which edge indices are fracture edge in pointer format
    std::vector<Edge *> fracture_non_boundary_edge_indices_pointer; //Storing which edge indices are non_boundary fracture edge in pointer format
    std::vector<Edge *> non_fracture_edge_indices_pointer;  //Storing which edge indices are non-fracture edge in pointer format
    std::vector<Edge *> non_fracture_and_boundary_edge_indices_pointer; //Storing which edge indices are non-fracture and boundary edge in pointer format
    std::unordered_map<Vertex *, std::vector<std::vector<Cell *>>> vertex_final_connected_pointer;  //Storing connected components of vertices in pointer format
    std::unordered_map<Edge *, std::vector<std::vector<Cell *>>> edge_final_connected_pointer;  //Storing connected components of edges in pointer format
    std::unordered_map<Face *, std::vector<std::vector<Cell *>>> face_final_connected_pointer;  //Storing connected components of faces in pointer format
    std::unordered_map<Vertex *,int> vertex_final_connected_size;    //Storing the size of the connected component at each vertex
    std::unordered_map<Edge *,int> edge_final_connected_size;  //Storing the size of the connected component at each edge
    std::unordered_map<Face *,int> face_final_connected_size;  //Storing the size of the connected component at each face
    std::unordered_map<Cell *,std::pair<size_t, std::vector<Face *>>> Cell_fractures;  //Storing the number of fracture faces in each cell and the corresponding fracture faces pointers

    
    int vertex_component(size_t iV,size_t iT)const; //Storing which component does cell iT belongs in the connected compononent of vertex iV

    inline int vertex_component(const Vertex & V, const Cell & T) const
    {
        return vertex_component(V.global_index(),T.global_index());
    }

    int edge_component(size_t iE,size_t iT)const;   //Storing which component does cell iT belongs in the connected compononent of edge iE

    inline int edge_component(const Edge & E, const Cell & T) const
    {
        return edge_component(E.global_index(),T.global_index());
    }

    int face_component(size_t iF,size_t iT)const;   //Storing which component does cell iT belongs in the connected compononent of face iF

    inline int face_component(const Face & F, const Cell & T) const
    {
        return face_component(F.global_index(),T.global_index());
    }

    int findFaceIndex(const Face& F) const;     //Storing the fracture face position for face F

    // int findFaceIndex(const Face& F, const Vertex& V) const;

    void generate();

private:
    std::string mesh_file_;
    Mesh& m_mesh;
    // std::vector<std::vector<Vertex *>>& fracture_pointer_;
};
}

#endif // CONNECTED_INFORMATION_GENERATOR_HPP