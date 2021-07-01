#include "MeshObject.hpp"
#include "MeshND.hpp"

#ifndef _MESH3D_HPP
#define _MESH3D_HPP

namespace Mesh3D
{
    using Mesh = MeshND::Mesh<3>;
    using Vertex = MeshND::Vertex<3>;
    using Edge = MeshND::Edge<3>;
    using Face = MeshND::Face<3>;
    using Cell = MeshND::Cell<3>;
    using VectorRd = MeshND::VectorRd<3>;
    using VectorZd = MeshND::VectorZd<3>;

    template <size_t object_dim>
    using Simplex = MeshND::Simplex<3, object_dim>;

    template <size_t object_dim>
    using Simplices = MeshND::Simplices<3, object_dim>;
} // namespace Mesh3D

#endif