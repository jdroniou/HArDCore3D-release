#include <array>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <math.hpp>
#include <algorithm>
// #include <Eigen/FullPivLU>

#ifndef _MESHOBJECT_HPP
#define _MESHOBJECT_HPP

namespace MeshND
{
    template <size_t space_dim>
    using VectorRd = Eigen::Matrix<double, space_dim, 1>;

    template <size_t space_dim>
    using VectorZd = Eigen::Matrix<int, space_dim, 1>;

    template <size_t space_dim, size_t object_dim>
    using Simplex = std::array<VectorRd<space_dim>, object_dim + 1>;

    template <size_t space_dim, size_t object_dim>
    using Simplices = std::vector<Simplex<space_dim, object_dim>>;

    /// Method to find the center mass of an arbitrary simplex in arbitrary space
    template <size_t space_dim, size_t object_dim>
    VectorRd<space_dim> simplex_center_mass(Simplex<space_dim, object_dim> simplex);

    /// Method to find the Lebesgue measure of an arbitrary simplex in arbitrary space
    template <size_t space_dim, size_t object_dim>
    double simplex_measure(Simplex<space_dim, object_dim> simplex);

    /*!
    * @defgroup Mesh
    * @brief Class providing a template for all objects (Vertex, Edge, Face, Cell) that appear in a Mesh
    */

    /*!
    * \addtogroup Mesh
    * @{
    */

    // ----------------------------------------------------------------------------
    //                         MeshObject class definition
    // ----------------------------------------------------------------------------

    /** MeshObject is a templated class describing, in the most general sense, an object of a mesh. It takes in two template parameters:
     *  space_dim - the dimension of the space the object is embedded in, and object_dim - the dimension of the object itself.
     **/

    /// Generic class to describe a cell, face, edge or vertex
    template <size_t space_dim, size_t object_dim>
    class MeshObject 
    {
    public:
        /// Constructor for a MeshObject defined by its simplices
        MeshObject(size_t index, Simplices<space_dim, object_dim> simplices);

        /// Constructor for a Simplex
        MeshObject(size_t index, Simplex<space_dim, object_dim> simplex);

        /// Constructor for a MeshObject defined by a single coordinate (Vertex)
        MeshObject(size_t index, VectorRd<space_dim> vertex);

        /// Null constructor
        MeshObject();

        /// Destructor
        ~MeshObject();

        // inline size_t dim() const { return object_dim; }

        inline size_t global_index() const { return _index; }                                ///< Return the global index of the MeshObject
        inline double diam() const { return _diameter; }                                     ///< Return the diameter of the MeshObject
        inline VectorRd<space_dim> center_mass() const { return _center_mass; }              ///< Return the center mass of the MeshObject
        inline double measure() const { return _measure; }                                   ///< Return the Lebesgue measure of the MeshObject
        inline Simplices<space_dim, object_dim> get_simplices() const { return _simplices; } ///< Return the simplices making up the MeshObject
        inline void set_global_index(const size_t idx) { _index = idx; }                           ///< Set the global index
        inline bool is_boundary() { return _is_boundary; }                                   ///< Return true if MeshObject is a boundary object, false otherwise
        inline void set_boundary(bool val) { _is_boundary = val; }                           ///< Set the boundary value of the MeshObject

        inline std::vector<MeshObject<space_dim, 0>*> get_vertices() const { return _vertices; }       ///< Return the vertices of the MeshObject
        inline std::vector<MeshObject<space_dim, 1>*> get_edges() const { return _edges; }             ///< Return the edges of the MeshObject
        inline std::vector<MeshObject<space_dim, space_dim - 1>*> get_faces() const { return _faces; } ///< Return the faces of the MeshObject
        inline std::vector<MeshObject<space_dim, space_dim>*> get_cells() const { return _cells; }     ///< Return the cells of the MeshObject

        inline size_t n_vertices() const { return _vertices.size(); } ///< Return the number of vertices of the MeshObject
        inline size_t n_edges() const { return _edges.size(); }       ///< Return the number of edges of the MeshObject
        inline size_t n_faces() const { return _faces.size(); }       ///< Return the number of faces of the MeshObject
        inline size_t n_cells() const { return _cells.size(); }       ///< Return the number of cells of the MeshObject

        MeshObject<space_dim, 0>* vertex(const size_t i) const;           ///< Return the i-th vertex of the MeshObject
        MeshObject<space_dim, 1>* edge(const size_t i) const;             ///< Return the i-th edge of the MeshObject
        MeshObject<space_dim, space_dim - 1>* face(const size_t i) const; ///< Return the i-th face of the MeshObject
        MeshObject<space_dim, space_dim>* cell(const size_t i) const;     ///< Return the i-th cell of the MeshObject

        void add_vertex(MeshObject<space_dim, 0>* vertex);         ///< Add a vertex to the MeshObject
        void add_edge(MeshObject<space_dim, 1>* edge);             ///< Add an edge to the MeshObject
        void add_face(MeshObject<space_dim, space_dim - 1>* face); ///< Add a face to the MeshObject
        void add_cell(MeshObject<space_dim, space_dim>* cell);     ///< Add a cell to the MeshObject

        int index_vertex(const MeshObject<space_dim, 0>* vertex) const;           ///< Returns the local index of a vertex
        int index_edge(const MeshObject<space_dim, 1>* edge) const;               ///< Returns the local index of an edge
        int index_face(const MeshObject<space_dim, space_dim - 1>* face) const;   ///< Returns the local index of a face
        int index_cell(const MeshObject<space_dim, space_dim>* cell) const;       ///< Returns the local index of a cell

        VectorRd<space_dim> coords() const; ///< Return the coordinates of a Vertex

        VectorRd<space_dim> face_normal(const size_t face_index) const; ///< Return the outer normal of a Cell towards the Face located at face_index
        VectorRd<space_dim> edge_normal(const size_t edge_index) const; ///< Return the edge normal of a 2D object

        int face_orientation(const size_t face_index) const; ///< Return the orientation of a Face
        int edge_orientation(const size_t edge_index) const; ///< Return the orientation of a Edge

        VectorRd<space_dim> normal() const;     ///< Return the normal of a Face
        VectorRd<space_dim> tangent() const;    ///< Return the tangent of a Edge

        void construct_face_normals(); ///< Set the directions of the face normals of a cell

    private:
        size_t _index;
        VectorRd<space_dim> _center_mass;
        double _measure;
        double _diameter;
        bool _is_boundary;
        Simplices<space_dim, object_dim> _simplices;
        VectorRd<space_dim> _normal;       // uninitialised unless object_dim == space_dim - 1 (face)
        std::vector<int> _face_directions; // empty unless object_dim == space_dim (cell)

        std::vector<MeshObject<space_dim, 0>*> _vertices;
        std::vector<MeshObject<space_dim, 1>*> _edges;
        std::vector<MeshObject<space_dim, space_dim - 1>*> _faces;
        std::vector<MeshObject<space_dim, space_dim>*> _cells;
    };

    /// A Vertex is a MeshObject with object_dim = 0
    template <size_t space_dim>
    using Vertex = MeshObject<space_dim, 0>;

    /// An Edge is a MeshObject with object_dim = 1
    template <size_t space_dim>
    using Edge = MeshObject<space_dim, 1>;

    /// A Face is a MeshObject with object_dim = space_dim - 1
    template <size_t space_dim>
    using Face = MeshObject<space_dim, space_dim - 1>;

    /// A Cell is a MeshObject with object_dim = space_dim
    template <size_t space_dim>
    using Cell = MeshObject<space_dim, space_dim>;

    //@}

    // -----------------------------------------------------------------
    // ----------- Implementation of templated functions ---------------
    // -----------------------------------------------------------------

    template <size_t space_dim, size_t object_dim>
    VectorRd<space_dim> simplex_center_mass(Simplex<space_dim, object_dim> simplex)
    {
        VectorRd<space_dim> center_mass = Eigen::VectorXd::Zero(space_dim);

        for (auto& coord : simplex)
        {
            for (size_t i = 0; i < space_dim; ++i)
            {
                center_mass(i) += coord(i);
            }
        }
        return center_mass / (object_dim + 1);
    }

    template <size_t space_dim, size_t object_dim>
    double simplex_measure(Simplex<space_dim, object_dim> simplex)
    {
        Eigen::MatrixXd CM_mat = Eigen::MatrixXd::Zero(object_dim + 2, object_dim + 2);
        for (size_t i = 0; i < object_dim + 1; ++i)
        {
            VectorRd<space_dim> i_coords = simplex[i];
            for (size_t j = i + 1; j < object_dim + 1; ++j)
            {
                VectorRd<space_dim> j_coords = simplex[j];
                double norm = (j_coords - i_coords).norm();
                CM_mat(i, j) = norm * norm;
                CM_mat(j, i) = CM_mat(i, j);
            }
            CM_mat(i, object_dim + 1) = 1;
            CM_mat(object_dim + 1, i) = 1;
        }

        // Calculate Cayley-Menger determinant
        double det = CM_mat.determinant();
        double scaling = std::pow(Math::factorial(object_dim), 2) * std::pow(2, object_dim);

        return std::sqrt(std::abs(det / scaling));
    }

    template <size_t space_dim, size_t object_dim>
    MeshObject<space_dim, object_dim>::MeshObject(size_t index, Simplices<space_dim, object_dim> simplices)
        : _index(index), _simplices(simplices)
    {
        assert(space_dim >= object_dim);
        assert(space_dim >= 2);

        if (object_dim == 0 || object_dim == 1)
        {
            assert(simplices.size() == 1);
        }

        _is_boundary = false;

        _measure = 0.0;
        _center_mass = Eigen::VectorXd::Zero(space_dim);
        std::vector<VectorRd<space_dim>> vertex_coords;
        for (auto& simplex : simplices)
        {
            // assert(simplex.size() == _dim + 1);
            double tmp = simplex_measure<space_dim, object_dim>(simplex);
            _measure += tmp;
            _center_mass += tmp * simplex_center_mass<space_dim, object_dim>(simplex);
            for (auto& coord : simplex)
            {
                if (std::find(vertex_coords.begin(), vertex_coords.end(), coord) == vertex_coords.end())
                {
                    vertex_coords.push_back(coord); // if coord not already in vertex_coords, add it
                }
            }
        }
        _center_mass /= _measure;

        _diameter = 0.0;
        for (auto it = vertex_coords.begin(); it != vertex_coords.end(); ++it)
        {
            for (auto jt = it; jt != vertex_coords.end(); ++jt)
            {
                _diameter = std::max(_diameter, (*it - *jt).norm()); // probably more efficient methods
            }
        }

        if (object_dim == space_dim - 1) // find the normal
        {
            Simplex<space_dim, object_dim> simplex = _simplices[0];
            Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(object_dim, space_dim);
            for (size_t i = 0; i < object_dim; ++i)
            {
                for (size_t j = 0; j < space_dim; ++j)
                {
                    mat(i, j) = simplex[i + 1][j] - simplex[0][j]; // find d-1 linearly independent vectors in the face, and put them in each col of mat
                }
            }

            // all normal vectors lie in null space of mat
            Eigen::FullPivLU<Eigen::MatrixXd> lu(mat);
            _normal = lu.kernel().normalized(); // not most efficient routine
        }
    }

    template <size_t space_dim, size_t object_dim>
    MeshObject<space_dim, object_dim>::MeshObject(size_t index, Simplex<space_dim, object_dim> simplex) // simplex
        : MeshObject(index, Simplices<space_dim, object_dim>(1, simplex))
    {
    }

    template <size_t space_dim, size_t object_dim>
    MeshObject<space_dim, object_dim>::MeshObject(size_t index, VectorRd<space_dim> vertex) // vertex
        : MeshObject(index, Simplex<space_dim, object_dim>({ vertex }))
    {
        assert(object_dim == 0);
    }

    template <size_t space_dim, size_t object_dim>
    MeshObject<space_dim, object_dim>::MeshObject() {}

    template <size_t space_dim, size_t object_dim>
    MeshObject<space_dim, object_dim>::~MeshObject() {}

    template <size_t space_dim, size_t object_dim>
    void MeshObject<space_dim, object_dim>::add_vertex(MeshObject<space_dim, 0>* vertex)
    {
        assert(std::find(_vertices.begin(), _vertices.end(), vertex) == _vertices.end());
        _vertices.push_back(vertex);
    }

    template <size_t space_dim, size_t object_dim>
    void MeshObject<space_dim, object_dim>::add_edge(MeshObject<space_dim, 1>* edge)
    {
        assert(std::find(_edges.begin(), _edges.end(), edge) == _edges.end());
        _edges.push_back(edge);
        if (space_dim == 2)
        {
            assert(std::find(_faces.begin(), _faces.end(), reinterpret_cast<MeshObject<space_dim, space_dim - 1> *>(edge)) == _faces.end());
            _faces.push_back(reinterpret_cast<MeshObject<space_dim, space_dim - 1> *>(edge));
        }
    }

    template <size_t space_dim, size_t object_dim>
    void MeshObject<space_dim, object_dim>::add_face(MeshObject<space_dim, space_dim - 1>* face)
    {
        assert(std::find(_faces.begin(), _faces.end(), face) == _faces.end());
        _faces.push_back(face);
        if (space_dim == 2)
        {
            assert(std::find(_edges.begin(), _edges.end(), reinterpret_cast<MeshObject<space_dim, 1> *>(face)) == _edges.end());
            _edges.push_back(reinterpret_cast<MeshObject<space_dim, 1> *>(face));
        }
    }

    template <size_t space_dim, size_t object_dim>
    void MeshObject<space_dim, object_dim>::add_cell(MeshObject<space_dim, space_dim>* cell)
    {
        assert(std::find(_cells.begin(), _cells.end(), cell) == _cells.end());
        _cells.push_back(cell);
    }

    template <size_t space_dim, size_t object_dim>
    MeshObject<space_dim, 0>* MeshObject<space_dim, object_dim>::vertex(const size_t i) const
    {
        assert(i < _vertices.size());
        return _vertices[i];
    }

    template <size_t space_dim, size_t object_dim>
    MeshObject<space_dim, 1>* MeshObject<space_dim, object_dim>::edge(const size_t i) const
    {
        assert(i < _edges.size());
        return _edges[i];
    }

    template <size_t space_dim, size_t object_dim>
    MeshObject<space_dim, space_dim - 1>* MeshObject<space_dim, object_dim>::face(const size_t i) const
    {
        assert(i < _faces.size());
        return _faces[i];
    }

    template <size_t space_dim, size_t object_dim>
    MeshObject<space_dim, space_dim>* MeshObject<space_dim, object_dim>::cell(const size_t i) const
    {
        assert(i < _cells.size());
        return _cells[i];
    }

    template <size_t space_dim, size_t object_dim>
    void MeshObject<space_dim, object_dim>::construct_face_normals() // not very efficient - probably room for improvement
    {
        assert(object_dim == space_dim);
        for (size_t iF = 0; iF < _faces.size(); ++iF)
        {
            VectorRd<space_dim> normal = _faces[iF]->normal();
            Simplex<space_dim, object_dim - 1> face_simplex = _faces[iF]->get_simplices()[0];
            VectorRd<space_dim> center = Eigen::VectorXd::Zero(space_dim);
            double count;
            for (auto& cell_simplex : this->get_simplices())
            {
                count = 0;
                for (size_t i = 0; (i < cell_simplex.size()) && count < 2; ++i)
                {
                    if (std::find(face_simplex.begin(), face_simplex.end(), cell_simplex[i]) == face_simplex.end())
                    {
                        ++count;
                    }
                }
                if (count == 1) // only don't share one coordinate
                {
                    center = simplex_center_mass<space_dim, object_dim>(cell_simplex);
                    break;
                }
            }
            assert(count == 1);
            //    _face_directions.push_back(Math::sgn((_faces[iF]->center_mass() - _center_mass).dot(normal))); // star shaped wrt center mass
            _face_directions.push_back(Math::sgn((_faces[iF]->center_mass() - center).dot(normal)));
            assert(_face_directions[iF] != 0);
        }
    }

    template <size_t space_dim, size_t object_dim>
    VectorRd<space_dim> MeshObject<space_dim, object_dim>::face_normal(const size_t face_index) const
    {
        assert(object_dim == space_dim);
        assert(face_index < _faces.size());
        return _face_directions[face_index] * _faces[face_index]->normal();
    }

    template <size_t space_dim, size_t object_dim>
    VectorRd<space_dim> MeshObject<space_dim, object_dim>::edge_normal(const size_t edge_index) const
    {
        assert(object_dim == 2);
        if (space_dim == 2)
        {
            return this->face_normal(edge_index);
        }
        if (space_dim == 3)
        {
            VectorRd<space_dim> edge_normal = _normal.cross(_edges[edge_index]->tangent()); // unit vector as normal and tangent are orthogonal unit vectors
            return edge_normal;
        }
    }

    template <size_t space_dim, size_t object_dim>
    VectorRd<space_dim> MeshObject<space_dim, object_dim>::coords() const // only for vertices
    {
        assert(object_dim == 0);
        return this->get_simplices()[0][0];
    }

    template <size_t space_dim, size_t object_dim>
    VectorRd<space_dim> MeshObject<space_dim, object_dim>::normal() const
    {
        assert(object_dim == space_dim - 1);
        return _normal;
    }

    template <size_t space_dim, size_t object_dim>
    VectorRd<space_dim> MeshObject<space_dim, object_dim>::tangent() const
    {
        assert(object_dim == 1);
        return (this->get_vertices()[1]->coords() - this->get_vertices()[0]->coords()).normalized();
    }

    template <size_t space_dim, size_t object_dim>
    int MeshObject<space_dim, object_dim>::index_vertex(const MeshObject<space_dim, 0>* vertex) const
    {
        auto itr = std::find(_vertices.begin(), _vertices.end(), vertex);
        if (itr != _vertices.end())
        {
            return itr - _vertices.begin();
        }
        else
        {
            throw "Vertex not found";
        }
    }

    template <size_t space_dim, size_t object_dim>
    int MeshObject<space_dim, object_dim>::index_edge(const MeshObject<space_dim, 1>* edge) const
    {
        auto itr = std::find(_edges.begin(), _edges.end(), edge);
        if (itr != _edges.end())
        {
            return itr - _edges.begin();
        }
        else
        {
            throw "Edge not found";
        }
    }

    template <size_t space_dim, size_t object_dim>
    int MeshObject<space_dim, object_dim>::index_face(const MeshObject<space_dim, space_dim - 1>* face) const
    {
        auto itr = std::find(_faces.begin(), _faces.end(), face);
        if (itr != _faces.end())
        {
            return itr - _faces.begin();
        }
        else
        {
            throw "Face not found";
        }
    }

    template <size_t space_dim, size_t object_dim>
    int MeshObject<space_dim, object_dim>::index_cell(const MeshObject<space_dim, space_dim>* cell) const
    {
        auto itr = std::find(_cells.begin(), _cells.end(), cell);
        if (itr != _cells.end())
        {
            return itr - _cells.begin();
        }
        else
        {
            throw "Cell not found";
        }
    }
    template <size_t space_dim, size_t object_dim>
    int MeshObject<space_dim, object_dim>::face_orientation(const size_t face_index) const
    {
        assert(object_dim == space_dim);
        assert(face_index < _faces.size());
        return _face_directions[face_index];
    }
    template <size_t space_dim, size_t object_dim>
    int MeshObject<space_dim, object_dim>::edge_orientation(const size_t edge_index) const
    {
        assert(object_dim == 2);
        assert(edge_index < _edges.size());
        return Math::sgn((_edges[edge_index]->center_mass() - this->center_mass()).dot(this->edge_normal(edge_index))); // assuming face is star-shaped wrt face-center !!
    }


} // namespace MeshND

#endif
