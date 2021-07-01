#include "DirectedGraph.hpp"

// ----------------------------------------------------------------------------
//                            DirectedFace
// ----------------------------------------------------------------------------

DirectedFace::DirectedFace() {}
DirectedFace::~DirectedFace() {}
DirectedFace::DirectedFace(std::vector<std::size_t> _V) : V(_V) {}

void DirectedFace::add_vertex(size_t v)
{
    assert(std::find(V.begin(), V.end(), v) == V.end());
    V.push_back(v);
}

// Returns true if two faces share same vertices, and are ordered the same way.
bool operator==(const DirectedFace &lhs, const DirectedFace &rhs)
{
    std::vector<std::size_t> l_vertices = lhs.V;
    std::vector<std::size_t> r_vertices = rhs.V;

    std::sort(std::begin(l_vertices), std::end(l_vertices), std::greater<size_t>());
    std::sort(std::begin(r_vertices), std::end(r_vertices), std::greater<size_t>());

    return l_vertices == r_vertices;
}

// ----------------------------------------------------------------------------
//                            DirectedCell
// ----------------------------------------------------------------------------

DirectedCell::DirectedCell() {}
DirectedCell::DirectedCell(std::size_t cell_num) { part.push_back(cell_num); }
DirectedCell::~DirectedCell() {}
DirectedCell::DirectedCell(std::vector<DirectedFace> &_T, std::size_t cell_num) : T(_T)
{
    part.push_back(cell_num);
}

void DirectedCell::add_face(DirectedFace &f)
{
    T.push_back(f);
}

void DirectedCell::remove_face(std::size_t pos)
{
    T.erase(std::begin(T) + pos);
}

void DirectedCell::append_cell(DirectedCell &K)
{
    T.insert(std::end(T), std::begin(K.T), std::end(K.T));
    part.insert(std::end(part), std::begin(K.part), std::end(K.part));
}

bool DirectedCell::has_face(DirectedFace f)
{
    for (std::size_t it = 0; it < T.size(); it++)
    {
        if (T[it] == f)
        {
            return true;
        }
    }
    return false;
}

void DirectedCell::print()
{
    for (std::size_t it = 0; it < T.size(); it++)
    {
        for (std::size_t jt = 0; jt < T[it].V.size(); jt++)
        {
            std::cout << T[it].V[jt] << ",";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

bool DirectedCell::check_cell()
{
    bool check = true;
    for (std::size_t it = 0; it < T.size(); it++)
    {
        for (std::size_t jt = it + 1; jt < T.size(); jt++)
        {
            check = check && (!(T[it] == T[jt]));
        }
    }
    return check;
}

// ----------------------------------------------------------------------------
//                            NodeArray
// ----------------------------------------------------------------------------

CellNodeArray::CellNodeArray() {}
CellNodeArray::~CellNodeArray() {}
CellNodeArray::CellNodeArray(Array &_A) : A(_A) {}
bool CellNodeArray::node_exists(std::size_t node)
{
    for (std::size_t i = 0; i < A.size(); i++)
    {
        if (find(begin(A[i]), end(A[i]), node) != end(A[i]))
        {
            return true;
        }
    }
    return false;
}

void CellNodeArray::renum_nodes(std::size_t node)
{
    for (std::size_t i = 0; i < A.size(); i++)
    {
        for (std::size_t j = 0; j < A[i].size(); j++)
        {
            if (A[i][j] > node)
            {
                A[i][j]--;
            }
        }
    }
}

void CellNodeArray::print(std::ofstream *out, std::size_t cell_num)
{
    *out << cell_num;
    *out << " " << A.size() << "\n";
    for (std::size_t cell_i = 0; cell_i < A.size(); cell_i++)
    {
        *out << " " << cell_i;
        *out << " " << A[cell_i].size();
        for (std::size_t node_j = 0; node_j < A[cell_i].size(); node_j++)
        {
            *out << " " << A[cell_i][node_j];
        }
        *out << "\n";
    }
}

GraphNodeArray::GraphNodeArray() {}
GraphNodeArray::~GraphNodeArray() {}
GraphNodeArray::GraphNodeArray(std::vector<CellNodeArray> &_GA) : GA(_GA) {}
bool GraphNodeArray::node_exists(std::size_t node)
{
    for (std::size_t i = 0; i < GA.size(); i++)
    {
        if (GA[i].node_exists(node))
        {
            return true;
        }
    }
    return false;
}

void GraphNodeArray::renum_nodes(std::size_t node)
{
    for (std::size_t i = 0; i < GA.size(); i++)
    {
        GA[i].renum_nodes(node);
    }
}

void GraphNodeArray::print(std::ofstream *out)
{
    *out << GA.size();
    *out << " 0\n";
    for (std::size_t i = 0; i < GA.size(); i++)
    {
        GA[i].print(out, i);
    }
}

// // ----------------------------------------------------------------------------
// //                            DirectedGraph
// // ----------------------------------------------------------------------------

DirectedGraph::DirectedGraph() {}
DirectedGraph::~DirectedGraph() {}
DirectedGraph::DirectedGraph(std::vector<DirectedCell> &_G) : G(_G) {}

void DirectedGraph::add_cell(DirectedCell &T)
{
    G.push_back(T);
}

void DirectedGraph::remove_cell(std::size_t pos)
{
    G.erase(std::begin(G) + pos);
}

bool DirectedGraph::test_graph()
{
    bool test = true;
    for (std::size_t i = 0; i < G.size(); ++i)
    {
        test = test && G[i].check_cell();
    }
    return test;
}

void DirectedGraph::randomise()
{
    std::srand(std::time(nullptr));
    std::random_shuffle(begin(G), end(G));
}

void DirectedGraph::order()
{
    for (std::size_t i = 0; i < G.size(); i++)
    {
        for (std::size_t j = i + 1; j < G.size(); ++j)
        {
            if (G[i].T.size() > G[j].T.size())
            {
                std::swap(G[i], G[j]);
            }
        }
    }
}

void DirectedGraph::coarsen()
{
    for (std::size_t it = 0; it < G.size() - 1; it++)
    {
        // Merge for any number of faces shared
        // for (std::size_t jt = it + 1; jt < G.size(); jt++)
        // {

        //     std::vector<size_t> faces_shared;
        //     if (can_merge(it, jt, faces_shared))
        //     {
        //         merge_cells(it, jt, faces_shared);
        //         break;
        //     }
        // }

        // Routine to merge largest number of shared faces
        bool flag = true;
        std::vector<std::vector<size_t>> shared_faces;
        for (std::size_t jt = it + 1; jt < G.size(); jt++)
        {
            std::vector<size_t> faces_shared;
            if (can_merge(it, jt, faces_shared))
            {
                if (faces_shared.size() >= G[it].T.size() - 3)
                {
                    merge_cells(it, jt, faces_shared);
                    flag = false;
                    break;
                }
                else
                {
                    shared_faces.push_back(faces_shared);
                }
            }
            else
            {
                std::vector<size_t> null_vec;
                shared_faces.push_back(null_vec);
            }
        }
        if (flag)
        {
            std::size_t pos = 0;
            for (std::size_t jt = 1; jt < shared_faces.size(); ++jt)
            {
                if (shared_faces[jt].size() > shared_faces[pos].size())
                {
                    pos = jt;
                }
            }
            if (shared_faces[pos].size() > 0)
            {
                merge_cells(it, pos + it + 1, shared_faces[pos]);
            }
        }
    }
}

GraphNodeArray DirectedGraph::graph_to_array()
{
    std::vector<CellNodeArray> cells;
    for (auto &iT : G)
    {
        Array array(iT.T.size());
        for (size_t i = 0; i < iT.T.size(); ++i)
        {
            for (size_t j = 0; j < iT.T[i].V.size(); ++j)
            {
                array[i].push_back(iT.T[i].V[j]);
            }
        }
        cells.push_back(CellNodeArray(array));
    }

    return GraphNodeArray(cells);
}

std::string DirectedGraph::get_partition()
{
    std::string partition = "";
    for (size_t i = 0; i < G.size(); i++)
    {
        partition += std::to_string(G[i].part.size());
        partition += " ";
        for (size_t j = 0; j < G[i].part.size(); j++)
        {
            partition += std::to_string(G[i].part[j]);
            partition += " ";
        }
        partition += "\n";
    }
    return partition;
}

// void DirectedGraph::plotfile(std::ofstream *out, std::vector<Eigen::VectorXd> &vertices)
// {
//     std::vector<DirectedFace> faces;
//     for (auto &cell : G)
//     {
//         for (auto &face : cell.T)
//         {
//             if (find(std::begin(faces), std::end(faces), face.anti_face()) == std::end(faces))
//             {
//                 faces.push_back(face);
//             }
//         }
//     }
//     for (auto &face : faces)
//     {
//         *out << vertices[face.a - 1](0) << " " << vertices[face.a - 1](1) << std::endl;
//         *out << vertices[face.b - 1](0) << " " << vertices[face.b - 1](1) << std::endl;
//         *out << std::endl;
//     }
// }

bool DirectedGraph::can_merge(std::size_t pos1, std::size_t pos2, std::vector<size_t> &faces_shared)
{
    std::vector<DirectedFace> cell1 = G[pos1].T;
    std::vector<DirectedFace> cell2 = G[pos2].T;
    bool vert_flag = true;
    for (std::size_t i = 0; i < cell1.size(); i++)
    {
        for (std::size_t j = 0; j < cell2.size(); j++)
        {
            if (cell1[i] == cell2[j])
            {
                faces_shared.push_back(i);
                faces_shared.push_back(j + cell1.size());
            }
        }
    }

    bool face_flag = (faces_shared.size() > 0);

    // if(face_flag)
    // {
    //     std::cout << "found\n";
    // }
    if (face_flag)
    {
        for (std::size_t i = 0; (i < cell1.size()) && vert_flag; i++)
        {
            for (std::size_t j = 0; (j < cell2.size()) && vert_flag; j++)
            {
                if (cell1[i] == cell2[j])
                {
                    continue;
                }
                for (auto &v : cell1[i].V)
                {
                    if (std::find(std::begin(cell2[j].V), std::end(cell2[j].V), v) != std::end(cell2[j].V)) // if share vertex
                    {
                        // search for vertex in one of shared faces
                        bool flag = true;
                        for (std::size_t k = 0; k < faces_shared.size(); k += 2)
                        {
                            std::size_t pos = faces_shared[k];
                            if (std::find(std::begin(cell1[pos].V), std::end(cell1[pos].V), v) != std::end(cell1[pos].V)) // vertex in face
                            {
                                flag = false;
                                break;
                            }
                        }
                        if (flag) //vert was not found in any of shared faces
                        {
                            // std::cout << "Found shared vertex and non-shared face -- ignoring cell pairing\n";
                            vert_flag = false;
                            break;
                        }
                    }
                }
            }
        }
    }

    return vert_flag && face_flag;
}

void DirectedGraph::merge_cells(std::size_t pos1, std::size_t pos2, std::vector<size_t> &faces_shared)
{
    G[pos1].append_cell(G[pos2]);
    remove_cell(pos2);
    std::sort(std::begin(faces_shared), std::end(faces_shared), std::greater<size_t>());
    for (std::size_t iF = 0; iF < faces_shared.size(); iF++)
    {
        G[pos1].remove_face(faces_shared[iF]);
    }

    // std::vector<DirectedFace> to_remove = G[pos1].order_faces();
    // if (to_remove.size() != 0)
    // {
    //     for (size_t i = 0; i < to_remove.size(); i++)
    //     {
    //         size_t j = 0;
    //         while (j < G.size())
    //         {
    //             if (j != pos1 && G[j].has_face(to_remove[i].anti_face()))
    //             {
    //                 G[pos1].part.insert(std::end(G[pos1].part), std::begin(G[j].part), std::end(G[j].part));
    //                 remove_cell(j);
    //                 if (j < pos1)
    //                 {
    //                     pos1--;
    //                 }
    //                 j--;
    //             }
    //             j++;
    //         }
    //     }
    //     test_graph();
    //     std::cout << "Successfully connected non-simply connected cell\n";
    // }
}

// bool DirectedGraph::random_bool(double prob)
// {
//     std::srand(std::time(nullptr));
//     std::bernoulli_distribution d(prob);
//     return d(rand_engine);
// }
