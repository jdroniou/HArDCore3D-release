#include "MeshObject.hpp"
#include "MeshND.hpp"

#ifndef _MESH2D_HPP
#define _MESH2D_HPP

namespace Mesh2D
{
    using Mesh = MeshND::Mesh<2>;
    using Vertex = MeshND::Vertex<2>;
    using Edge = MeshND::Edge<2>;
    using Face = MeshND::Face<2>;
    using Cell = MeshND::Cell<2>;
    using VectorRd = MeshND::VectorRd<2>;
    using VectorZd = MeshND::VectorZd<2>;

    template <size_t object_dim>
    using Simplex = MeshND::Simplex<2, object_dim>;

    template <size_t object_dim>
    using Simplices = MeshND::Simplices<2, object_dim>;
} // namespace Mesh2D

#endif