#include "MeshObject.hpp"

#ifndef _MESHND_HPP
#define _MESHND_HPP


/*!
* \addtogroup Mesh
* @{
*/

namespace MeshND
{
    /// Class to describe a mesh
    template <std::size_t dimension>
    class Mesh
    {
    public:
        Mesh() {}
        ~Mesh()
        {
            for (auto &vertex : _vertices)
            {
                delete vertex;
            }
            for (auto &edge : _edges)
            {
                delete edge;
            }
            if (dimension != 2) // if dimension == 2, faces are edges (which are already deleted)
            {
                for (auto &face : _faces)
                {
                    delete face;
                }
            }
            for (auto &cell : _cells)
            {
                delete cell;
            }
        }

        void set_name(std::string name) { _mesh_name = name; } ///< set the name of the mesh
        inline std::string get_name() { return _mesh_name; }   ///< getter for the mesh name

        double h_max() const ///< max diameter of cells
        {
            double val = 0.0;
            for (auto &cell : _cells)
            {
                val = std::max(val, cell->diam());
            }
            return val;
        }
        inline std::size_t dim() const { return dimension; } ///< dimension of the mesh

        inline std::size_t n_vertices() const { return _vertices.size(); } ///< number of vertices in the mesh.
        inline std::size_t n_edges() const { return _edges.size(); }       ///< number of edges in the mesh.
        inline std::size_t n_faces() const { return _faces.size(); }       ///< number of faces in the mesh.
        inline std::size_t n_cells() const { return _cells.size(); }       ///< number of cells in the mesh.
        inline std::size_t n_elems(size_t d) const  ///< number of d-cells in the mesh
        {
          assert(d <= dimension && dimension < 4);
          if (d == dimension) {
            return n_cells();
          } else if (d == 2) {
            return n_faces(); 
          } else if (d == 1) {
            return n_edges();
          } else {
            return n_vertices();
          }
        }

        inline std::size_t n_b_vertices() const { return _b_vertices.size(); } ///< number of boundary vertices in the mesh.
        inline std::size_t n_b_edges() const { return _b_edges.size(); }       ///< number of boundary edges in the mesh.
        inline std::size_t n_b_faces() const { return _b_faces.size(); }       ///< number of boundary faces in the mesh.
        inline std::size_t n_b_cells() const { return _b_cells.size(); }       ///< number of boundary cells in the mesh.

        inline std::size_t n_i_vertices() const { return _i_vertices.size(); } ///< number of internal vertices in the mesh.
        inline std::size_t n_i_edges() const { return _i_edges.size(); }       ///< number of internal edges in the mesh.
        inline std::size_t n_i_faces() const { return _i_faces.size(); }       ///< number of internal faces in the mesh.
        inline std::size_t n_i_cells() const { return _i_cells.size(); }       ///< number of internal cells in the mesh.

        inline std::vector<Vertex<dimension> *> get_vertices() const { return _vertices; } ///< lists the vertices in the mesh.
        inline std::vector<Edge<dimension> *> get_edges() const { return _edges; }         ///< lists the edges in the mesh.
        inline std::vector<Face<dimension> *> get_faces() const { return _faces; }         ///< lists the faces in the mesh.
        inline std::vector<Cell<dimension> *> get_cells() const { return _cells; }         ///< lists the cells in the mesh.

        inline std::vector<Vertex<dimension> *> get_b_vertices() const { return _b_vertices; } ///< lists the boundary vertices in the mesh.
        inline std::vector<Edge<dimension> *> get_b_edges() const { return _b_edges; }         ///< lists the boundary edges in the mesh.
        inline std::vector<Face<dimension> *> get_b_faces() const { return _b_faces; }         ///< lists the boundary faces in the mesh.
        inline std::vector<Cell<dimension> *> get_b_cells() const { return _b_cells; }         ///< lists the boundary cells in the mesh.

        inline std::vector<Vertex<dimension> *> get_i_vertices() const { return _i_vertices; } ///< lists the internal vertices in the mesh.
        inline std::vector<Edge<dimension> *> get_i_edges() const { return _i_edges; }         ///< lists the internal edges in the mesh.
        inline std::vector<Face<dimension> *> get_i_faces() const { return _i_faces; }         ///< lists the internal faces in the mesh.
        inline std::vector<Cell<dimension> *> get_i_cells() const { return _i_cells; }         ///< lists the internal cells in the mesh.

        void add_vertex(Vertex<dimension> *vertex) ///<  adds a vertex to the mesh
        {
            assert(std::find(_vertices.begin(), _vertices.end(), vertex) == _vertices.end());
            _vertices.push_back(vertex);
        }
        void add_edge(Edge<dimension> *edge) ///<  adds a edge to the mesh
        {
            assert(std::find(_edges.begin(), _edges.end(), edge) == _edges.end());
            _edges.push_back(edge);
            if (dimension == 2)
            {
                assert(std::find(_faces.begin(), _faces.end(), reinterpret_cast<Face<dimension> *>(edge)) == _faces.end());
                _faces.push_back(reinterpret_cast<Face<dimension> *>(edge));
            }            
        }
        void add_face(Face<dimension> *face) ///<  adds a face to the mesh
        {
            assert(std::find(_faces.begin(), _faces.end(), face) == _faces.end());
            _faces.push_back(face);
            if (dimension == 2)
            {
                assert(std::find(_edges.begin(), _edges.end(), reinterpret_cast<Edge<dimension> *>(face)) == _edges.end());
                _edges.push_back(reinterpret_cast<Edge<dimension> *>(face));
            }
        }
        void add_cell(Cell<dimension> *cell) ///<  adds a cell to the mesh
        {
            assert(std::find(_cells.begin(), _cells.end(), cell) == _cells.end());
            _cells.push_back(cell);
        }

        void add_b_vertex(Vertex<dimension> *vertex) ///<  adds a boundary vertex to the mesh
        {
            assert(std::find(_b_vertices.begin(), _b_vertices.end(), vertex) == _b_vertices.end());
            _b_vertices.push_back(vertex);
        }
        void add_b_edge(Edge<dimension> *edge) ///<  adds a boundary edge to the mesh
        {
            assert(std::find(_b_edges.begin(), _b_edges.end(), edge) == _b_edges.end());
            _b_edges.push_back(edge);
            if (dimension == 2)
            {
                assert(std::find(_b_faces.begin(), _b_faces.end(), reinterpret_cast<Face<dimension> *>(edge)) == _b_faces.end());
                _b_faces.push_back(reinterpret_cast<Face<dimension> *>(edge));
            }
        }
        void add_b_face(Face<dimension> *face) ///<  adds a boundary face to the mesh
        {
            assert(std::find(_b_faces.begin(), _b_faces.end(), face) == _b_faces.end());
            _b_faces.push_back(face);
            if (dimension == 2)
            {
                assert(std::find(_b_edges.begin(), _b_edges.end(), reinterpret_cast<Edge<dimension> *>(face)) == _b_edges.end());
                _b_edges.push_back(reinterpret_cast<Edge<dimension> *>(face));
            }
        }
        void add_b_cell(Cell<dimension> *cell) ///<  adds a boundary cell to the mesh
        {
            assert(std::find(_b_cells.begin(), _b_cells.end(), cell) == _b_cells.end());
            _b_cells.push_back(cell);
        }

        void add_i_vertex(Vertex<dimension> *vertex) ///<  adds an internal vertex to the mesh
        {
            assert(std::find(_i_vertices.begin(), _i_vertices.end(), vertex) == _i_vertices.end());
            _i_vertices.push_back(vertex);
        }
        void add_i_edge(Edge<dimension> *edge) ///<  adds an internal edge to the mesh
        {
            assert(std::find(_i_edges.begin(), _i_edges.end(), edge) == _i_edges.end());
            _i_edges.push_back(edge);
            if (dimension == 2)
            {
                assert(std::find(_i_faces.begin(), _i_faces.end(), reinterpret_cast<Face<dimension> *>(edge)) == _i_faces.end());
                _i_faces.push_back(reinterpret_cast<Face<dimension> *>(edge));
            }
        }
        void add_i_face(Face<dimension> *face) ///<  adds an internal face to the mesh
        {
            assert(std::find(_i_faces.begin(), _i_faces.end(), face) == _i_faces.end());
            _i_faces.push_back(face);
            if (dimension == 2)
            {
                assert(std::find(_i_edges.begin(), _i_edges.end(), reinterpret_cast<Edge<dimension> *>(face)) == _i_edges.end());
                _i_edges.push_back(reinterpret_cast<Edge<dimension> *>(face));
            }
        }
        void add_i_cell(Cell<dimension> *cell) ///<  adds an internal cell to the mesh
        {
            assert(std::find(_i_cells.begin(), _i_cells.end(), cell) == _i_cells.end());
            _i_cells.push_back(cell);
        }

        // Note that all these assume that a MeshObject's index is equal to its position in the mesh!!
        inline Vertex<dimension> *vertex(std::size_t index) const { assert(index < _vertices.size()); return _vertices[index]; } ///<  get a constant pointer to a vertex using its global index
        inline Edge<dimension> *edge(std::size_t index) const { assert(index < _edges.size()); return _edges[index]; }        ///<  get a constant pointer to a edge using its global index
        inline Face<dimension> *face(std::size_t index) const { assert(index < _faces.size()); return _faces[index]; }        ///<  get a constant pointer to a face using its global index
        inline Cell<dimension> *cell(std::size_t index) const { assert(index < _cells.size()); return _cells[index]; }        ///<  get a constant pointer to a cell using its global index

        inline Vertex<dimension> *b_vertex(std::size_t index) const { assert(index < _b_vertices.size()); return _b_vertices[index]; } ///<  get a constant pointer to a boundary vertex using an index
        inline Edge<dimension> *b_edge(std::size_t index) const { assert(index < _b_edges.size()); return _b_edges[index]; }        ///<  get a constant pointer to boundary a edge using an index
        inline Face<dimension> *b_face(std::size_t index) const { assert(index < _b_faces.size()); return _b_faces[index]; }        ///<  get a constant pointer to boundary a face using an index
        inline Cell<dimension> *b_cell(std::size_t index) const { assert(index < _b_cells.size()); return _b_cells[index]; }        ///<  get a constant pointer to boundary a cell using an index

        inline Vertex<dimension> *i_vertex(std::size_t index) const { assert(index < _i_vertices.size()); return _i_vertices[index]; } ///<  get a constant pointer to an internal vertex using an index
        inline Edge<dimension> *i_edge(std::size_t index) const { assert(index < _i_edges.size()); return _i_edges[index]; }        ///<  get a constant pointer to an internal edge using an index
        inline Face<dimension> *i_face(std::size_t index) const { assert(index < _i_faces.size()); return _i_faces[index]; }        ///<  get a constant pointer to an internal face using an index
        inline Cell<dimension> *i_cell(std::size_t index) const { assert(index < _i_cells.size()); return _i_cells[index]; }        ///<  get a constant pointer to an internal cell using an index

        // Return the indices of the boundary of the i-th top dimensionnal cell of the MeshObject
        std::vector<size_t> get_boundary(size_t d,size_t index) const
        {
          std::vector<size_t> rv;
          static_assert(dimension < 4);
          assert(d <= dimension);
          if (d == dimension) {
            assert(index < _cells.size());
            const Cell<dimension> & T = *_cells[index];
            rv.reserve(T.n_faces());
            for (size_t j = 0; j < T.n_faces(); ++j) {
              rv.push_back(T.face(j)->global_index());
            }
          } else if (d == 2) { // Only reached for dimension = 3
            assert(index < _faces.size());
            const Face<dimension> & F = *_faces[index];
            rv.reserve(F.n_edges());
            for (size_t j = 0; j < F.n_edges(); ++j) {
              rv.push_back(F.edge(j)->global_index());
            }
          } else if (d == 1) {
            assert(index < _edges.size());
            const Edge<dimension> & E = *_edges[index];
            rv.reserve(2);
            rv.push_back(E.vertex(0)->global_index());
            rv.push_back(E.vertex(1)->global_index());
          } 
          return rv;
        }
        // Return the orientation of the j-th boundary element of the i-th d-cell 
        int boundary_orientation(size_t d,size_t i, size_t j) const 
        {
          static_assert(dimension < 4);
          assert(d <= dimension);
          if (d == dimension) {
            assert(i < _cells.size());
            return _cells[i]->induced_orientation(j);
          } else if (d == 2) { // Only reached for dimension = 3
            assert(i < _faces.size());
            return _faces[i]->induced_orientation(j);
          } else { // Special case for degenerate boundary
            assert(d == 1);
            return (j%2==0?-1:1);
          } 
        }

        std::vector<double> regularity()
        {
            /// Regularity factor =
            ///   1st component: maximum of
            ///      * diameter of cell / (measure of cell)^{1/dim}
            ///      * diameter of cell / diameter of face  [for each face of the cell]
            ///
            ///   2nd component: evaluation of max of ratio "diam of cell / radius ball inscribed in cell"

            std::vector<std::vector<double>> reg_cell(n_cells(), {0.0, 0.0});
            std::size_t count = 0;
            for (auto &T : _cells)
            {
                double hT = T->diam();
                VectorRd<dimension> xT = T->center_mass();

                reg_cell[count][0] = hT / pow(T->measure(), 1.0 / dimension);

                double rhoT = hT;
                std::vector<Face<dimension> *> faces = T->get_faces();
                for (auto &F : faces)
                {
                    double hF = F->diam();
                    VectorRd<dimension> xF = F->center_mass();
                    VectorRd<dimension> nTF = F->normal(); // sign does not matter

                    reg_cell[count][0] = std::max(reg_cell[count][0], hT / hF);

                    rhoT = std::min(rhoT, std::abs((xT - xF).dot(nTF))); // If xT is not in T, is this really a good measure?
                }
                reg_cell[count][1] = hT / rhoT;
                ++count; //could just use iterators
            }

            std::vector<double> value(2, 0.0);
            for (size_t iT = 0; iT < n_cells(); iT++)
            {
                value[0] = std::max(value[0], reg_cell[iT][0]);
                value[1] = std::max(value[1], reg_cell[iT][1]);
            }

            return value;
        }

        void renum(const char B, const std::vector<size_t> new_to_old)
        {

            switch (B)
            {
            case 'C':
            {
                std::vector<Cell<dimension> *> old_index = _cells;
                for (size_t i = 0; i < _cells.size(); i++)
                {
                    old_index[new_to_old[i]]->set_global_index(i);
                    _cells[i] = old_index[new_to_old[i]];
                }
                break;
            }
            case 'F':
            {
                std::vector<Face<dimension> *> old_index = _faces;
                for (size_t i = 0; i < _faces.size(); i++)
                {
                    old_index[new_to_old[i]]->set_global_index(i);
                    _faces[i] = old_index[new_to_old[i]];
                }
                break;
            }

            case 'E':
            {
                std::vector<Edge<dimension> *> old_index = _edges;
                for (size_t i = 0; i < _edges.size(); i++)
                {
                    old_index[new_to_old[i]]->set_global_index(i);
                    _edges[i] = old_index[new_to_old[i]];
                }
                break;
            }

            case 'V':
            {
                std::vector<Vertex<dimension> *> old_index = _vertices;
                for (size_t i = 0; i < _vertices.size(); i++)
                {
                    old_index[new_to_old[i]]->set_global_index(i);
                    _vertices[i] = old_index[new_to_old[i]];
                }
                break;
            }
            }
        }

    private:
        std::string _mesh_name;

        std::vector<Vertex<dimension> *> _vertices;
        std::vector<Edge<dimension> *> _edges;
        std::vector<Face<dimension> *> _faces;
        std::vector<Cell<dimension> *> _cells;

        // std::vector<MeshObject<dimension, 0> *> _vertices;
        // std::vector<MeshObject<dimension, 1> *> _edges;
        // std::vector<MeshObject<dimension, dimension - 1> *> _faces;
        // std::vector<MeshObject<dimension, dimension> *> _cells;

        std::vector<Vertex<dimension> *> _b_vertices;
        std::vector<Edge<dimension> *> _b_edges;
        std::vector<Face<dimension> *> _b_faces;
        std::vector<Cell<dimension> *> _b_cells;

        std::vector<Vertex<dimension> *> _i_vertices;
        std::vector<Edge<dimension> *> _i_edges;
        std::vector<Face<dimension> *> _i_faces;
        std::vector<Cell<dimension> *> _i_cells;
    };
} // namespace MeshND

#endif
