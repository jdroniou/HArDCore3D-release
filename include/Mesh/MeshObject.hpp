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

        /// Returns the global index of the MeshObject
        inline size_t global_index() const { return _index; }
        /// Returns the diameter of the MeshObject
        inline double diam() const { return _diameter; }
        /// Returns the center mass of the MeshObject
        inline VectorRd<space_dim> center_mass() const { return _center_mass; }
        ///< Returns the Lebesgue measure of the MeshObject
        inline double measure() const { return _measure; }

        /// Returns the simplices making up the MeshObject
        inline Simplices<space_dim, object_dim> get_simplices() const { return _simplices; }
        inline Simplices<space_dim, object_dim>& set_simplices() { return _simplices; }
        /// Set the global index
        inline void set_global_index(const size_t idx) { _index = idx; }
        /// Returns true if MeshObject is a boundary object, false otherwise
        inline bool is_boundary() const { return _is_boundary; }
        /// Set the boundary value of the MeshObject
        inline void set_boundary(bool val) { _is_boundary = val; }
        /// Returns true if MeshObject is flat (only relevant for faces), false otherwise
        inline bool is_flat() const { return _is_flat; }

        /// Returns the vertices of the MeshObject
        inline std::vector<MeshObject<space_dim, 0>*> get_vertices() const { return _vertices; }
        /// Returns the edges of the MeshObject
        inline std::vector<MeshObject<space_dim, 1>*> get_edges() const { return _edges; }
        /// Returns the faces of the MeshObject
        inline std::vector<MeshObject<space_dim, space_dim - 1>*> get_faces() const { return _faces; }
        /// Return the cells of the MeshObject
        inline std::vector<MeshObject<space_dim, space_dim>*> get_cells() const { return _cells; }

        /// Returns the number of vertices of the MeshObject
        inline size_t n_vertices() const { return _vertices.size(); }
        /// Returns the number of edges of the MeshObject 
        inline size_t n_edges() const { return _edges.size(); }
        /// Returns the number of faces of the MeshObject
        inline size_t n_faces() const { return _faces.size(); }
        /// Returns the number of cells of the MeshObject
        inline size_t n_cells() const { return _cells.size(); }

        MeshObject<space_dim, 0>* vertex(const size_t i) const;           ///< Returns the i-th vertex of the MeshObject
        MeshObject<space_dim, 1>* edge(const size_t i) const;             ///< Returns the i-th edge of the MeshObject
        MeshObject<space_dim, space_dim - 1>* face(const size_t i) const; ///< Returns the i-th face of the MeshObject
        MeshObject<space_dim, space_dim>* cell(const size_t i) const;     ///< Returns the i-th cell of the MeshObject

        void add_vertex(MeshObject<space_dim, 0>* vertex);         ///< Add a vertex to the MeshObject
        void add_edge(MeshObject<space_dim, 1>* edge);             ///< Add an edge to the MeshObject
        void add_face(MeshObject<space_dim, space_dim - 1>* face); ///< Add a face to the MeshObject
        void add_cell(MeshObject<space_dim, space_dim>* cell);     ///< Add a cell to the MeshObject

        int index_vertex(const MeshObject<space_dim, 0>* vertex) const;           ///< Returns the local index of a vertex
        int index_edge(const MeshObject<space_dim, 1>* edge) const;               ///< Returns the local index of an edge
        int index_face(const MeshObject<space_dim, space_dim - 1>* face) const;   ///< Returns the local index of a face
        int index_cell(const MeshObject<space_dim, space_dim>* cell) const;       ///< Returns the local index of a cell

        VectorRd<space_dim> coords() const; ///< Return the coordinates of a Vertex
        void set_coords(const VectorRd<space_dim> & x); ///< Set the coordinates of a Vertex

        VectorRd<space_dim> face_normal(const size_t face_index) const; ///< Return the outer normal of a Cell towards the Face located at face_index
        VectorRd<space_dim> edge_normal(const size_t edge_index) const; ///< Return the edge normal of a 2D object

        int face_orientation(const size_t face_index) const; ///< Return the orientation of a Face
        int edge_orientation(const size_t edge_index) const; ///< Return the orientation of a Edge (multiplied by edge_normal, gives the outer normal to the face)

        int induced_orientation(const size_t index) const; ///< Return 1 if the volume form on the boundary is positively oriented, -1 else

        VectorRd<space_dim> normal() const;     ///< Return the normal of a Face
        VectorRd<space_dim> tangent() const;    ///< Return the tangent of a Edge
        Eigen::Matrix<double,space_dim,2> face_tangent() const; ///< Return the tangent space of a Face

        void construct_face_orientations(); ///< Set the directions of the face normals of a cell

    private:
        size_t _index;
        VectorRd<space_dim> _center_mass;
        double _measure;
        double _diameter;
        bool _is_boundary;
        bool _is_flat;  // To record non-planar faces
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

      _is_flat = true;
      
      if (object_dim == space_dim - 1) // find the normal to the face and check if it is flat
      {
        // We compute the normal using the most perpendicular pair of vectors created from the center to any vertex
        // Find all points in the face that are not the center of mass
        std::vector<VectorRd<space_dim>> points_in_face;
        for (size_t i = 0; i < _simplices.size(); ++i)
        {
          for (size_t j = 0; j < space_dim; ++j)
          {
            if( (std::find(points_in_face.begin(), points_in_face.end(), _simplices[i][j]) == points_in_face.end()) && ( (_simplices[i][j]-_center_mass).norm()>1e-8 * _diameter) )
            {
              points_in_face.push_back(_simplices[i][j]);
            }
          }
        }
        
        VectorRd<space_dim> perp_vec = VectorRd<space_dim>::Zero();
        double norm_perp_vec = 0.;
        for (size_t iV = 0 ; iV < points_in_face.size(); iV++){
          VectorRd<space_dim> t1 = (points_in_face[iV] - _center_mass).normalized();
          for (size_t iVp = iV+1; iVp < points_in_face.size(); iVp++){
            VectorRd<space_dim> t2 = (points_in_face[iVp] - _center_mass).normalized();
            VectorRd<space_dim> tmp = t1.cross(t2);
            double norm_tmp = tmp.norm();
            if (norm_tmp > norm_perp_vec){ 
              perp_vec = tmp;
              norm_perp_vec = norm_tmp;
            }
          }
        }
        _normal = perp_vec / norm_perp_vec;

        // Checking that the face is planar: all points_in_face should lie in the plane orthogonal to _normal at the center of mass
        for (size_t iV = 0; iV < points_in_face.size(); iV++){
          VectorRd<space_dim> v_xF = points_in_face[iV]-_center_mass;
          if (std::abs( v_xF.dot(_normal) ) > 1e-8 * v_xF.norm()){
            _is_flat = false;
          }
        }
        
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
    void MeshObject<space_dim, object_dim>::construct_face_orientations() // not very efficient - probably room for improvement
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
    void MeshObject<space_dim, object_dim>::set_coords(const VectorRd<space_dim> & x) // only for vertices
    {
        assert(object_dim == 0);
        this->set_simplices()[0][0] = x;
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
    Eigen::Matrix<double,space_dim,2> MeshObject<space_dim, object_dim>::face_tangent() const
    {
        assert(object_dim == 2);
        Eigen::Matrix<double,space_dim,2> rv;
        rv.col(0) = edge(0)->tangent();
        rv.col(1) = edge_normal(0);
        return rv;
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
    template <size_t space_dim, size_t object_dim>
    int MeshObject<space_dim, object_dim>::induced_orientation(const size_t index) const
    {
        assert((object_dim == 2 && index < _edges.size())||(object_dim == 3 && index < _faces.size()));
        if constexpr(object_dim == 2) {
          VectorRd<space_dim> const n = edge_orientation(index)*edge_normal(index), 
                                    t = edge(index)->tangent();
          Eigen::Matrix<double,space_dim,2> const ab = face_tangent();
          return Math::sgn(ab.col(0).dot(n)*ab.col(1).dot(t) - ab.col(0).dot(t)*ab.col(1).dot(n));
          // If vol_f = da^db, vol_E = dt and n is the outward unit vector
          // then da^db = (da(n)*db(t) - da(t)*db(n))dn^dt, hence i_{n}(da^db) = (da(n)*db(t) - da(t)*db(n))dn^dt
          // The induced orientation on E is then sgn(da(n)*db(t) - da(t)*db(n)) dt 
          // sgn(da(n)*db(t) - da(t)*db(n)) should be equal to da(n)*db(t) - da(t)*db(n)
        } else { // object_dim == 3
          Eigen::Matrix<double,space_dim,object_dim> nab;
          nab.rightCols(2) = face(index)->face_tangent();
          nab.leftCols(1) = face_normal(index); // This seems inconsistent with the behavior in 2D, why does this one is outer and not edge_orientation?
          return Math::sgn(nab.determinant());
          // We assume that the volume form on T is dx^dy^dz, then dx^dy^dz = det(n,a,b) dn^da^db
          // The induced orientation on F is then sign(det(n,a,b))
        }
    }

} // namespace MeshND

#endif
