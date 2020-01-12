#include "quadraturerule.hpp"

#include <quad3d.hpp>
#include <quad3d_face.hpp>
#include <quad1d.hpp>

#include <edge.hpp>
#include <vertex.hpp>

namespace HArDCore3D
{

  QuadratureRule generate_quadrature_rule(
					  const Cell & T,
					  const int doe,
					  const bool force_split
					  )
  {  
    QuadRuleTetra quadCell(std::max(doe,0), true);
    QuadratureRule quad;

    size_t nfaces = T.n_faces();

    if ( (nfaces == 4) && (!force_split) ){
      // Tetrahedron
      std::vector<Vertex *> vertices = T.get_vertices();
      auto x0 = vertices[0]->coords();
      auto x1 = vertices[1]->coords();
      auto x2 = vertices[2]->coords();
      auto x3 = vertices[3]->coords();
      double xT[] = {x0.x(), x1.x(), x2.x(), x3.x()};
      double yT[] = {x0.y(), x1.y(), x2.y(), x3.y()};
      double zT[] = {x0.z(), x1.z(), x2.z(), x3.z()};

      quadCell.setup(xT, yT, zT);

      for (size_t iqn = 0; iqn < quadCell.nq(); iqn++) {
	quad.emplace_back(quadCell.xq(iqn), quadCell.yq(iqn), quadCell.zq(iqn), quadCell.wq(iqn));
      }
    } else {
      // Generic region: split into tetrahedra
      auto xC = T.center_mass();

      for (auto& face : T.get_faces()) {
	const size_t nbedges = face->n_edges();

	if (nbedges == 3){
	  // The face is a triangle, xC-face already creates a tetrahedron
	  std::vector<Vertex *> vertices = face->get_vertices();
	  auto x0 = vertices[0]->coords();
	  auto x1 = vertices[1]->coords();
	  auto x2 = vertices[2]->coords();
	  double xT[] = {xC.x(), x0.x(), x1.x(), x2.x()};
	  double yT[] = {xC.y(), x0.y(), x1.y(), x2.y()};
	  double zT[] = {xC.z(), x0.z(), x1.z(), x2.z()};

	  quadCell.setup(xT, yT, zT);

	  for (size_t iqn = 0; iqn < quadCell.nq(); iqn++) {
	    quad.emplace_back(quadCell.xq(iqn), quadCell.yq(iqn), quadCell.zq(iqn), quadCell.wq(iqn));
	  }

	} else if(nbedges == 4){
	  // The face is a rectangle, we just need to split it into two to
	  // create tetrahedras with xC
	  std::vector<Vertex *> vertices = face->get_vertices();
	  auto x0 = vertices[0]->coords();
	  for (size_t isplit=0; isplit<2; isplit++){
	    auto x1 = vertices[1+isplit]->coords();
	    auto x2 = vertices[2+isplit]->coords();
	    double xT[] = {xC.x(), x0.x(), x1.x(), x2.x()};
	    double yT[] = {xC.y(), x0.y(), x1.y(), x2.y()};
	    double zT[] = {xC.z(), x0.z(), x1.z(), x2.z()};

	    quadCell.setup(xT, yT, zT);

	    for (size_t iqn = 0; iqn < quadCell.nq(); iqn++) {
	      quad.emplace_back(quadCell.xq(iqn), quadCell.yq(iqn), quadCell.zq(iqn), quadCell.wq(iqn));
	    }
	  }

	} else {
	  // Generic face, we split it into triangles
	  auto xF = face->center_mass();

	  for (auto& edge : face->get_edges()) {

	    auto x0 = edge->vertex(0)->coords();
	    auto x1 = edge->vertex(1)->coords();

	    // Vertices of the tetrahedron from the cell center
	    // to the face center containing the current edge
	    double xT[] = {xC.x(), xF.x(), x0.x(), x1.x()};
	    double yT[] = {xC.y(), xF.y(), x0.y(), x1.y()};
	    double zT[] = {xC.z(), xF.z(), x0.z(), x1.z()};

	    quadCell.setup(xT, yT, zT);

	    for (size_t iqn = 0; iqn < quadCell.nq(); iqn++) {
	      quad.emplace_back(quadCell.xq(iqn), quadCell.yq(iqn), quadCell.zq(iqn), quadCell.wq(iqn));
	    }
	  }
	}
      }
    }

    return quad;
  }

  //------------------------------------------------------------------------------

  QuadratureRule generate_quadrature_rule(
					  const Face & F,
					  const int doe
					  )
  {
    QuadRuleTriangle quadFace(std::max(doe,0), true);
    QuadratureRule quad;

    size_t nedges = F.n_edges();

    if (nedges == 3){
      // Face is a triangle, direct rule
      std::vector<Vertex *> vertices = F.get_vertices();
      auto x0 = vertices[0]->coords();
      auto x1 = vertices[1]->coords();
      auto x2 = vertices[2]->coords();
      double xT[] = {x0.x(), x1.x(), x2.x()};
      double yT[] = {x0.y(), x1.y(), x2.y()};
      double zT[] = {x0.z(), x1.z(), x2.z()};
                  
      quadFace.setup(xT, yT, zT);

      for (size_t iqn = 0; iqn < quadFace.nq(); iqn++) {
	quad.emplace_back(quadFace.xq(iqn), quadFace.yq(iqn), quadFace.zq(iqn), quadFace.wq(iqn));
      }

    } else if (nedges == 4){
      // Face is a rectangle, we split in two
      std::vector<Vertex *> vertices = F.get_vertices();
      auto x0 = vertices[0]->coords();

      for (size_t isplit=0; isplit<2; isplit++){
	auto x1 = vertices[1+isplit]->coords();
	auto x2 = vertices[2+isplit]->coords();
	double xT[] = {x0.x(), x1.x(), x2.x()};
	double yT[] = {x0.y(), x1.y(), x2.y()};
	double zT[] = {x0.z(), x1.z(), x2.z()};

	quadFace.setup(xT, yT, zT);

	for (size_t iqn = 0; iqn < quadFace.nq(); iqn++) {
	  quad.emplace_back(quadFace.xq(iqn), quadFace.yq(iqn), quadFace.zq(iqn), quadFace.wq(iqn));
	}
      }

    } else {

      auto xF = F.center_mass();

      for (auto& edge : F.get_edges()) {

	auto x0 = edge->vertex(0)->coords();
	auto x1 = edge->vertex(1)->coords();

	// Vertices of the triangle form by the face center and the current edge
	double xT[] = {xF.x(), x0.x(), x1.x()};
	double yT[] = {xF.y(), x0.y(), x1.y()};
	double zT[] = {xF.z(), x0.z(), x1.z()};

	quadFace.setup(xT, yT, zT);

	for (size_t iqn = 0; iqn < quadFace.nq(); iqn++) {
	  quad.emplace_back(quadFace.xq(iqn), quadFace.yq(iqn), quadFace.zq(iqn), quadFace.wq(iqn));
	}
      }
    }

    return quad;
  }
 
  //------------------------------------------------------------------------------

  QuadratureRule generate_quadrature_rule(
					  const Edge & E,
					  const int doe
					  )
  {
    QuadRuleEdge quadEdge(std::max(doe,0), true);
    QuadratureRule quad;

    auto x0 = E.vertex(0)->coords();
    auto x1 = E.vertex(1)->coords();
    double xT[] = {x0.x(), x1.x()};
    double yT[] = {x0.y(), x1.y()};
    double zT[] = {x0.z(), x1.z()};
  
    quadEdge.setup(xT, yT, zT);
    for (size_t iqn = 0; iqn < quadEdge.nq(); iqn++) {
        quad.emplace_back(quadEdge.xq(iqn), quadEdge.yq(iqn), quadEdge.zq(iqn), quadEdge.wq(iqn));
    }
    return quad;
  }
}
