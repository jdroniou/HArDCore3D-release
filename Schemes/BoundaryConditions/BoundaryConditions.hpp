// Class to provide description of boundary conditions
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
*
*	This library was developed around HHO methods, although some parts of it have a more
* general purpose. If you use this code or part of it in a scientific publication, 
* please mention the following book as a reference for the underlying principles
* of HHO schemes:
*
* The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications.
* D. A. Di Pietro and J. Droniou. 2019, 516p. 
* url: https://hal.archives-ouvertes.fr/hal-02151813.
*
*/


#ifndef _BOUNDARY_CONDITIONS_HPP
#define _BOUNDARY_CONDITIONS_HPP

#include <iostream>
#include <string>
#include <mesh.hpp>

using namespace HArDCore3D;

/*!
* @defgroup TestCases
*	@brief Defines test cases (exact solutions, source terms, boundary condition, etc.)
*/

// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

// @addtogroup TestCases
//@{

/// Structure to define a plane (by a point and a normal) and check position of a mesh entity with respect to it
struct Plane {
  /// Constructor
    Plane(const VectorRd & point, const VectorRd & normal):
      point(point),
      normal(normal)
      {
        // do nothing
      };
    
    VectorRd point; ///< One point on the plane
    VectorRd normal; ///< Normal to the plane 
    
    /// Check if a point is inside the plane
    bool is_in(const VectorRd & x) const;
    /// Check if a vertex is inside the plane
    inline bool is_in(const Vertex & V) const { return is_in(V.coords()); };
    /// Check if a mesh entity is inside the plane
    template<typename MeshEntity>
    bool is_in(const MeshEntity & E) const{
      for (auto V : E.get_vertices()){
        if (! is_in(*V) ) return false;
      }
      return true;
    };

    /// Check if a point is on one side of the plane (side=+/-1 depending on the normal to the plane) 
    bool is_on_side(const VectorRd & x, const double & side) const;
    /// Check if a vertex is on one side of the plane (side=+/-1 depending on the normal to the plane) 
    inline bool is_on_side(const Vertex & V, const double & side) const { return is_on_side(V.coords(), side); };
    /// Check if a mesh entity is on one side of the plane (side=+/-1 depending on the normal to the plane) 
    template<typename MeshEntity>
    bool is_on_side(const MeshEntity & E, const double & side) const{
      for (auto V : E.get_vertices()){
        if (! is_on_side(*V, side) ) return false;
      }
      return true;
    };
    
};

/// Structure to define a flat convex hull (nodes and normal) and check if a mesh entity is inside
struct Hull {
  /// Constructor
    Hull(const std::vector<VectorRd> & nodes, const VectorRd & normal):
      nodes(nodes),
      normal(normal.normalized())
      {
        // Create center and apex
        center = VectorRd::Zero();
        for (size_t i=0; i<nodes.size(); i++){
          center += nodes[i];
        }
        center /= nodes.size();
        apex = center + normal;
      };
    
    std::vector<VectorRd> nodes; ///< Nodes defining the planar convex hull
    VectorRd normal; ///< Normal to the planar convex hull 
    VectorRd center; ///< "Center" of the hull
    VectorRd apex; ///< A point outside the plane of the convex hull (in the direction of normal)
    
    /// Check if a point is in the convex hull
    bool is_in(const VectorRd & x) const;
    /// Check if a vertex is in the convex hull
    inline bool is_in(const Vertex & V) const { return is_in(V.coords()); };
    /// Check if a mesh entity is in the convex hull
    template<typename MeshEntity>
    bool is_in(const MeshEntity & E) const{
      for (auto V : E.get_vertices()){
        if (! is_in(*V) ) return false;
      }
      return true;
    };

};


// Various definitions of planes and convex hulls used in the BCs below
const Plane plane_x_zero(VectorRd::Zero(), VectorRd(1., 0., 0.));
const Plane plane_x_quarter(VectorRd(.25,0.,0.), VectorRd(1.,0.,0.));
const Plane plane_x_half(VectorRd(.5,0.,0.), VectorRd(1.,0.,0.));
const Hull lower_bottom_corner_x_zero({VectorRd(0.,0.,0.), VectorRd(0.,.25,0.), VectorRd(0.,.25,.25), VectorRd(0.,0.,.25)}, VectorRd(1., 0., 0.));
const Hull lower_bottom_corner_x_one({VectorRd(1.,0.,0.), VectorRd(1.,.25,0.), VectorRd(1.,.25,.25), VectorRd(1.,0.,.25)}, VectorRd(1., 0., 0.));


/// The BoundaryConditions class provides definition of boundary conditions
class BoundaryConditions {

public:
  /// Initialise data
  BoundaryConditions(
    const std::string bc_id,  ///< The identifier for the boundary condition (D, N or Mx)
    Mesh& mesh          ///< reference to the mesh
  );

	/// Test the boundary condition of an face.
	const std::string type(
		                    const Face& face    ///< Face on which to check the boundary condition
	                      ) const ;
  /**< @returns "dir" or "neu" depending if the face is a Dirichlet or Neumann boundary condition, as determined by m_bc_id below. For an internal face, returns "int".
    m_bc_id = "D": all faces are Dirichlet
    m_bc_id = "N": all faces are Neumann
    m_bc_id = "Mx" (x=0,1,...): various types of mixed boundary conditions, some faces will be Dirichlet and some will be Neumman.
   **/

  /// Test the boundary condition of an edge
  const std::string type(
                         const Edge& edge    ///< Edge on which to check the boundary condition
                         ) const ;

  /// Test the boundary condition of a vertex
  const std::string type(
                         const Vertex& vertex    ///< Vertex on which to check the boundary condition
                         ) const ;


  /// Returns the number of Dirichlet faces
  inline const size_t n_dir_faces() const {
    return m_n_dir_faces;
  };

  /// Returns the number of Dirichlet edges
  inline const size_t n_dir_edges() const {
    return m_n_dir_edges;
  };

  /// Returns the number of Dirichlet vertices
  inline const size_t n_dir_vertices() const {
    return m_n_dir_vertices;
  };

  /// Returns the identifier of the BC
  inline const std::string id() { return m_bc_id; };

  /// Returns the complete name of the boundary condition
  inline const std::string name() const {
    if (m_bc_id == "D"){
      return "Dirichlet";
    }else if (m_bc_id == "N"){
      return "Neumann";
    }else if (m_bc_id == "M0"){
      return "Mixed #0 (Dirichlet at x=0, Neumann elsewhere)";
    }else if (m_bc_id == "M1"){
      return "Mixed #1 (Dirichlet at x=0 a small rectangle, Neumann elsewhere)";
    }else if (m_bc_id == "M2"){
      return "Mixed #2 (Dirichlet between planes x=.25 and x=.5, Neumann elsewhere)";
    }
    std::cout << "Unknown boundary conditions: " << m_bc_id << "\n";
    exit(1);
  };

  /// Re-order faces of the mesh to put the Dirichlet faces at the position "pos" (=end or start)
  void reorder_faces(const std::string pos = "end");

  /// Re-order edges of the mesh to put the Dirichlet edges at the position "pos" (=end or start)
  void reorder_edges(const std::string pos = "end");

  /// Re-order vertices of the mesh to put the Dirichlet vertices at the position "pos" (=end or start)
  void reorder_vertices(const std::string pos = "end");


private:
  // Parameters: id of boundary condition, reference to mesh
  const std::string m_bc_id;
  Mesh& m_mesh;

  // Number of Dirichlet faces, edges and vertices
  size_t m_n_dir_faces;
  size_t m_n_dir_edges;
  size_t m_n_dir_vertices;

};


//@}

#endif //_BOUNDARY_CONDITION_HPP
