#include "mesh_builder.hpp"

using namespace HArDCore3D;

// Split 2D object into simplices
Simplices<2> simplexify(std::vector<VectorRd> vertices)
{
    assert(vertices.size() >= 3);
    Simplices<2> simplices;

if(vertices.size() == 3)
    {
        Simplex<2> simplex = {vertices[0], vertices[1], vertices[2]};
        simplices.push_back(simplex);
        return simplices;
    }
    
    VectorRd center = VectorRd::Zero();
    for(auto & v : vertices)
    {
        center += v;
    }
    center /= vertices.size();

    for(size_t iV = 0; iV < vertices.size(); ++iV)
    {
        Simplex<2> simplex = {center, vertices[iV], vertices[ (iV + 1) % vertices.size() ]}; // star shaped wrt vertex average.
        simplices.push_back(simplex);
    }
    return simplices;
}

MeshBuilder::MeshBuilder() {}
MeshBuilder::MeshBuilder(const std::string mesh_file) : _mesh_file(mesh_file) {}

std::unique_ptr<Mesh> MeshBuilder::build_the_mesh()
{
    MeshReaderRF mesh_reader(_mesh_file);

    std::size_t find_coarse = _mesh_file.find(".coarse.");

    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<std::vector<size_t>>> cells;
    std::vector<std::vector<size_t>> partition;

    mesh_reader.read_node_file(vertices);
    mesh_reader.read_ele_file(cells);

    bool coarse_mesh = (find_coarse != std::string::npos);

    if (coarse_mesh) //coarse mesh
    {
        mesh_reader.read_partition_file(partition);
    }

    if (vertices.size() > 0 && cells.size() > 0)
    {
        std::unique_ptr<Mesh> mesh = std::make_unique<Mesh>(); // make a pointer to the mesh so that it outlives the builder
        std::unique_ptr<Mesh> fine_mesh = std::make_unique<Mesh>();
        if (coarse_mesh)
        {
            // build the fine mesh first
            std::string fine_mesh_file = _mesh_file.substr(0, find_coarse);
            MeshBuilder fine_mesh_builder(fine_mesh_file);
            fine_mesh = fine_mesh_builder.build_the_mesh();
            std::cout << "     Coarse ";
        }
        else
        {
            std::cout << "     ";
        }

        std::cout << "Mesh: ";

        // Create vertices
        for (auto &v : vertices)
        {
            VectorRd vert(v[0], v[1], v[2]);
            Vertex *vertex = new Vertex(mesh->n_vertices(), vert);
            mesh->add_vertex(vertex);
        }

        // Create cells
        double total_vol = 0.0;
        for (auto &c : cells)
        {
            // build edges and faces of cell
            std::vector<Face *> cell_faces;
            std::vector<Edge *> cell_edges;
            std::vector<VectorRd> cell_vertex_coords;
            std::vector<std::size_t> cell_vertex_ids;
            for (auto &f : c)
            {
                std::vector<std::size_t> vertex_ids;
                std::vector<VectorRd> vertex_coords;
                for (auto &vertID : f)
                {
                    vertex_ids.push_back(vertID);
                    VectorRd coord({vertices[vertID][0], vertices[vertID][1], vertices[vertID][2]});
                    vertex_coords.push_back(coord);
                    if (std::find(cell_vertex_ids.begin(), cell_vertex_ids.end(), vertID) == cell_vertex_ids.end())
                    {
                        cell_vertex_ids.push_back(vertID);
                        cell_vertex_coords.push_back(coord);
                    }
                }

                // make edges
                bool make_face_flag = true;
                std::vector<Edge *> face_edges;
                std::vector<std::size_t> need_to_add;
                Face *face;
                for (std::size_t i = 0; i < vertex_ids.size(); ++i)
                {
                    std::size_t plus = (i + 1) % vertex_ids.size();
                    Eigen::Vector2d e(vertex_ids[i], vertex_ids[plus]);

                    bool flag = true;

                    std::vector<Vertex *> vlist = mesh->vertex(e(0))->get_vertices();
                    for (size_t j = 0; j < vlist.size(); j++)
                    {
                        if (vlist[j]->global_index() == e(1)) // The edge exists in the mesh
                        {
                            std::vector<Edge *> elist = mesh->vertex(e(0))->get_edges();
                            Edge *edge = elist[j];
                            need_to_add.push_back(face_edges.size());
                            face_edges.push_back(edge);
                            if (std::find(cell_edges.begin(), cell_edges.end(), edge) == cell_edges.end())
                            {
                                cell_edges.push_back(edge);
                            }
                            flag = false;
                            break;
                        }
                    }

                    if (flag)
                    {
                        VectorRd coords1({vertices[e(0)][0], vertices[e(0)][1], vertices[e(0)][2]});
                        VectorRd coords2({vertices[e(1)][0], vertices[e(1)][1], vertices[e(1)][2]});

                        Simplex<1> edge_coords({coords1, coords2});
                        Edge *edge = new Edge(mesh->n_edges(), edge_coords);

                        cell_edges.push_back(edge);
                        face_edges.push_back(edge);

                        mesh->add_edge(edge);

                        Vertex *vertex1 = mesh->vertex(e(0));
                        Vertex *vertex2 = mesh->vertex(e(1));

                        edge->add_vertex(vertex1);
                        edge->add_vertex(vertex2);
                        vertex1->add_edge(edge);
                        vertex2->add_edge(edge);

                        // if vertices are not already connected, connect them
                        std::vector<Vertex *> verts_of_vert1 = vertex1->get_vertices();
                        if (std::find(verts_of_vert1.begin(), verts_of_vert1.end(), vertex2) == verts_of_vert1.end())
                        {
                            vertex1->add_vertex(vertex2);
                            vertex2->add_vertex(vertex1);
                        }

                        if (make_face_flag)
                        {
                            face = new Face(mesh->n_faces(), simplexify(vertex_coords));
                            make_face_flag = false;
                            mesh->add_face(face);
                            for (auto &vertID : vertex_ids)
                            {
                                face->add_vertex(mesh->vertex(vertID));
                                mesh->vertex(vertID)->add_face(face);
                            }
                        }

                        face->add_edge(edge);
                        edge->add_face(face);
                    }
                }

                if (make_face_flag) // face has not been made
                {
                    // check if face exists, look in faces of edges
                    bool found = false;
                    for (auto &edge_face : face_edges[0]->get_faces())
                    {
                        std::size_t count = 0;
                        for (auto &e : edge_face->get_edges())
                        {
                            if (std::find(face_edges.begin(), face_edges.end(), e) == face_edges.end())
                            {
                                break;
                            }
                            ++count;
                        }
                        if (count == edge_face->n_edges())
                        {
                            face = edge_face;
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                    {
                        // make face
                        face = new Face(mesh->n_faces(), simplexify(vertex_coords));
                        mesh->add_face(face);
                        for (auto &e : face_edges)
                        {
                            face->add_edge(e);
                            e->add_face(face);
                        }
                        for (auto &vertID : vertex_ids)
                        {
                            face->add_vertex(mesh->vertex(vertID));
                            mesh->vertex(vertID)->add_face(face);
                        }
                    }
                }
                else // need to add remaining edges to the face
                {
                    for (auto &add : need_to_add)
                    {
                        face->add_edge(face_edges[add]);
                        face_edges[add]->add_face(face);
                    }
                }
                cell_faces.push_back(face);
            }

            Simplices<3> cell_simplices;

            if (coarse_mesh)
            {
                for (auto &part : partition[mesh->n_cells()])
                {
                    for (auto &simplex : fine_mesh->cell(part)->get_simplices())
                    {
                        cell_simplices.push_back(simplex);
                    }
                }
            }
            else
            {
                // to create cell simplices, need point in cell
                // use average of vertices - guarenteed to be in cell if cell is convex
                VectorRd point_in_cell({0, 0, 0});
                for (std::size_t i = 0; i < 3; ++i)
                {
                    for (auto &vert_coord : cell_vertex_coords)
                    {
                        point_in_cell[i] += vert_coord[i];
                    }
                    point_in_cell[i] /= cell_vertex_coords.size();
                }

                for (auto &face : cell_faces)
                {
                    for (auto &face_simplex : face->get_simplices())
                    {
                        Simplex<3> cell_simplex;
                        cell_simplex[0] = point_in_cell;
                        for (std::size_t i = 0; i < 3; ++i)
                        {
                            cell_simplex[i + 1] = face_simplex[i];
                        }
                        cell_simplices.push_back(cell_simplex);
                    }
                }
            }

            Cell *cell = new Cell(mesh->n_cells(), cell_simplices);
            total_vol += cell->measure();
            mesh->add_cell(cell);

            for (auto &face : cell_faces)
            {
                cell->add_face(face);
                face->add_cell(cell);
            }

            for (auto &edge : cell_edges)
            {
                cell->add_edge(edge);
                edge->add_cell(cell);
            }

            for (auto &vertID : cell_vertex_ids)
            {
                cell->add_vertex(mesh->vertex(vertID));
                mesh->vertex(vertID)->add_cell(cell);
            }
            cell->construct_face_orientations();
        }

        // // build boundary
        build_boundary(mesh.get());

        std::cout << "added " << mesh->n_cells() << " cells; Total volume = " << total_vol << std::endl;

        // Check that all faces are flat
        std::vector<size_t> non_flat_faces;
        for (auto & F : mesh->get_faces()){
          if (!F->is_flat()){
            non_flat_faces.push_back(F->global_index());
          }
        }
        if (non_flat_faces.size()>0){
          std::cout << "     **** WARNING: there are non-flat faces: ";
          for (size_t i = 0; i < std::min(non_flat_faces.size(), size_t(5)); ++i){
            std::cout << non_flat_faces[i] << ", ";
          }
          if (non_flat_faces.size()>5){
            std::cout << "etc. ";
          }
          std::cout << "(" << non_flat_faces.size() << " faces in total)" << std::endl;
          std::cout << "  You might want to consider using FlattenFaces (in the repository) to solve that issue." << std::endl << std::endl;
        }

        return mesh;
    }
    else
    {
        throw "     Cannot build mesh. Check input file\n";
    }
    return NULL;
}

void MeshBuilder::build_boundary(Mesh *mesh)
{
    // Here we fill in the _boundary variables of the cells and vertices, and the lists of boundary
    // edges, cells and vertices
    for (auto &face : mesh->get_faces())
    {
        std::vector<Cell *> cells = face->get_cells();
        if (cells.size() == 1)
        {
            // The cell has a boundary edge, so it is a boundary cell
            cells[0]->set_boundary(true);
            face->set_boundary(true);
            // mesh->add_b_cell(cells[0]);
            mesh->add_b_face(face);

            std::vector<Vertex *> face_verts = face->get_vertices();
            std::vector<Edge *> face_edges = face->get_edges();

            for (auto &v : face->get_vertices())
            {
                std::vector<Vertex *> b_verts = mesh->get_b_vertices();
                if (std::find(b_verts.begin(), b_verts.end(), v) == b_verts.end()) // if vert not in b_verts, add to b_verts
                {
                    v->set_boundary(true);
                    mesh->add_b_vertex(v);
                }
            }
            for (auto &e : face->get_edges())
            {
                std::vector<Edge *> b_edges = mesh->get_b_edges();
                if (std::find(b_edges.begin(), b_edges.end(), e) == b_edges.end()) // if edge not in b_edges, add to b_edges
                {
                    e->set_boundary(true);
                    mesh->add_b_edge(e);
                }
            }
        }
        else
        {
            mesh->add_i_face(face);
        }
    }

    // Pass to fill in interior elements
    for (auto &cell : mesh->get_cells())
    {
        if (!(cell->is_boundary()))
        {
            mesh->add_i_cell(cell);
        }
        else
        {
            mesh->add_b_cell(cell);
        }
    }
    for (auto &edge : mesh->get_edges())
    {
        if (!(edge->is_boundary()))
        {
            mesh->add_i_edge(edge);
        }
    }
    for (auto &vertex : mesh->get_vertices())
    {
        if (!(vertex->is_boundary()))
        {
            mesh->add_i_vertex(vertex);
        }
    }
}
