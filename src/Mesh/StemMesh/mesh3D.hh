#ifndef _MESH3DV_HH
#define _MESH3DV_HH

#include <cassert>

#include <algorithm>
#include <string>
#include <vector>

namespace StemMesh3D {

/// The interface for working with generic 3D meshes
/** The 3D mesh class encapsulates the data structures required to store a 3D mesh and provides an interface
    for accessing data about the mesh that is necessary for implementing numerical schemes. */
class mesh_3Dv {

public:
  /// Type used for element flags
  using flag_type = int;
  /// Default flag for elements with no associated flag
  static constexpr flag_type UNSET_FLAG = -999;

private:
  // friendships and type declarations 
  friend class mesh3Dv_builder;

  static constexpr size_t DIM = 3;

  size_t nV, nE, nF, nR;
  std::string mesh_name;

  // coordinates
  std::vector<std::array<double, DIM>> V_coords;

  // bounding box & mesh size
  mutable bool bb_status;
  mutable double bb_min[DIM], bb_max[DIM];

  mutable bool ms_status;
  mutable double ms_hmax, ms_hvol, ms_havg, ms_hmin;

  // primary datasets
  std::vector<std::array<size_t, 2>> EdgeVrtx;
  std::vector<std::vector<size_t>> FaceEdge;
  std::vector<std::vector<size_t>> RegnFace;

  std::vector<std::vector<bool>> FaceEdgeOk;
  std::vector<std::vector<bool>> RegnFaceOk;

  // tranposed datasets
  std::vector<std::vector<size_t>> FaceRegn;
  std::vector<std::vector<size_t>> EdgeFace;
  std::vector<std::vector<size_t>> VrtxEdge;

  std::vector<std::vector<bool>> EdgeFaceOk;
  std::vector<std::vector<bool>> VrtxEdgeOk;

  // lists of boundary items
  std::vector<size_t> bnd_vrtx;
  std::vector<size_t> bnd_edge;
  std::vector<size_t> bnd_face;
  std::vector<size_t> bnd_regn;

  // external flags
  std::vector<flag_type> fV;
  std::vector<flag_type> fE;
  std::vector<flag_type> fF;
  std::vector<flag_type> fR;

  // aux geometrical quantities
  std::vector<std::array<double, DIM>> R_coords;    // Region barycenters
  std::vector<std::array<double, DIM>> F_coords;    // Face barycenters
  std::vector<std::array<double, DIM>> F_nor;    // Face normal vectors
  std::vector<double> R_volume;                  // Region volumes
  std::vector<double> R_diam;                  // Region diameters
  std::vector<double> F_area;                  // Face areas
  std::vector<double> F_diam;                  // Face diameters

  // sort and remove duplicates from a list
  template<typename T>
  void shrink_list(std::vector<T> &tmp_list) const;

  // private methods
  void eval_h() const;

  inline void eval_bbox() const;

public:
  mesh_3Dv() : nV(0), nE(0), nF(0), nR(0), mesh_name("mesh-3D"), bb_status(false), ms_status(false) {}

  /// Returns the number of regions (3D cells / control volumes) in the mesh
  size_t n_region() const { return nR; }

  /// Returns the number of faces in the mesh
  size_t n_face() const { return nF; }

  /// Returns the number of edges in the mesh/
  size_t n_edge() const { return nE; }

  /// Returns the number of vertices in the mesh
  size_t n_vertex() const { return nV; }

  /// Returns the number of regions on the boundary
  size_t n_bregion() const { return bnd_regn.size(); }

  /// Returns the number of faces on the boundary
  size_t n_bface() const { return bnd_face.size(); }

  /// Returns the number of edges on the boundary
  size_t n_bedge() const { return bnd_edge.size(); }

  /// Returns the number of vertices on the boundary
  size_t n_bvertex() const { return bnd_vrtx.size(); }

  /// Returns the global region number of the (ilR)'th boundary region
  size_t get_bnd_regn(
      size_t ilR             /**< The number of the boundary region (0 < ilR < n_bregion()) */
  ) const;

  /// Returns the global face number of the (ilF)'th boundary face
  size_t get_bnd_face(
      size_t ilF             /**< The number of the boundary face (0 < ilF < n_bface()) */
  ) const;

  /// Returns the global edge number of the (ilE)'th boundary edge
  size_t get_bnd_edge(
      size_t ilE             /**< The number of the boundary edge (0 < ilE < n_bedge()) */
  ) const;

  /// Returns the global vertex number of the (ilV)'th boundary vertex
  size_t get_bnd_vrtx(
      size_t ilV             /**< The number of the boundary vertex (0 < ilV < n_bvertex()) */
  ) const;

  // TOPOLOGICAL METHODS:

  /// Returns the global face number of the (ilF)'th local face of the region iR
  size_t regn_face(
      size_t iR,             /**< The global region number of the region */
      size_t ilF             /**< The local face number of the face in region iR (0 < ilF < n_regn_face(size_t iR)) */
  ) const;

  /// Returns the global edge number of the (ilE)'th local edge of the face iF
  size_t face_edge(
      size_t iF,             /**< The global face number of the face */
      size_t ilE             /**< The local edge number of the edge in the face iF (0 < ilE < n_face_edge(size_t iF)) */
  ) const;

  /// Returns the global vertex number of the (ilV)'th vertex of the edge iE
  size_t edge_vrtx(
      size_t iE,             /**< The global edge number of the edge */
      size_t ilV             /**< The local vertex number of the vertex in the edge iE (0 < ilV < n_edge_vrtx(size_t iE)) */
  ) const;

  /// Returns the global region number of the (ilR)'th region bordering the face iF
  size_t face_regn(
      size_t iF,             /**< The global face number of the face */
      size_t ilR             /**< The local region number of the region adjacent to the face iF (0 < ilR < n_face_regn(size_t iF)) */
  ) const;

  /// Returns the global face number of the (ilF)'th face bordering the edge iE
  size_t edge_face(
      size_t iE,             /**< The global edge number of the edge */
      size_t ilF             /**< The local face number of the face adjacent to the edge iE (0 < ilF < n_edge_face(size_t iE)) */
  ) const;

  /// Returns the global edge number of the (ilE)'th edge adjacent to the vertex iV
  size_t vrtx_edge(
      size_t iV,             /**< The global vertex number of the vertex */
      size_t ilE             /**< The local edge number of the edge adjacent to the vertex iV (0 < ilE < n_vrtx_edge(size_t iV)) */
  ) const;

  /// Returns the global vertex number of the (ilV)'th adjacent vertex to vertex iV.
  /** A vertex \f$V2\f$ is considered adjacent to vertex \f$V1\f$ if there is an edge with endpoints \f$V1\f$ and \f$V2\f$ */
  size_t vrtx_vrtx(
      size_t iV,             /**< The global vertex number of the vertex */
      size_t ilV             /**< The local vertex number of the adjacent vertex (0 < ilV < n_vrtx_vrtx(size_t iV)) */
  ) const;

  /// Returns the number of local faces that are adjacent to the given region
  size_t n_regn_face(
      size_t iR              /**< The global region number of the region */
  ) const;

  /// Returns the number of local edges that are adjacent to the given face
  size_t n_face_edge(
      size_t iF              /**< The global face number of the face */
  ) const;

  /// Returns the number of local vertices that are adjacent to the given edge. This should always return 2.
  size_t n_edge_vrtx(
      size_t iE              /**< The global edge number of the edge */
  ) const;

  /// Returns the number of regions that are adjacent to the given face.
  size_t n_face_regn(
      size_t iF              /**< The global face number of the face */
  ) const;

  /// Returns the number of faces that are adjacent to the given edge.
  size_t n_edge_face(
      size_t iE              /**< The global edge number of the edge */
  ) const;

  /// Returns the number of edges that are adjacent to the given vertex.
  size_t n_vrtx_edge(
      size_t iV              /**< The global vertex number of the vertex */
  ) const;
  /// Returns the number of vertices that are adjacent to the given vertex.
  /** This is the same as the number of edges that are adjacent to the given vertex. */
  size_t n_vrtx_vrtx(
      size_t iV              /**< The global vertex number of the vertex */
  ) const;

  /// Returns true if the global face corresponding to ilF is oriented correctly with respect to iR
  /** Each internal face of the mesh is adjacent to two regions. Hence, the orientation of the global
      face can only agree with one of the two local faces. This function returns true if the global
      face corresponding to ilF with respect to iR is indeed oriented correctly with respect to the given
      region. A face is said to be oriented correctly if the cross product of consecutive edges along
      the face points outward of the region. ie. the normal vector of the face points outwards. */
  bool ok_regn_face(
      size_t iR,             /**< The global region number of the region */
      size_t ilF             /**< The local face number of the adjacent face */
  ) const;
  /// Returns true if the global edge corresponding to ilE is oriented correctly with respect to iF
  /** Each internal edge of the mesh is adjacent to many faces. Hence, the orientation of the global
      edge may only agree with some of them. This function returns true if the global edge corresponding
      to ilE with respect to iF agrees with the orientation of the local edge, ie. the order of their
      vertices is the same */
  bool ok_face_edge(
      size_t iF,             /**< The global face number of the face */
      size_t ilE             /**< The local edge number of the adjacent edge */
  ) const;
  /// Returns true if the global face corresponding to ilF is oriented correctly with respect to iE
  /** Each edge of the mesh is adjacent to many faces. Hence, the orientation of the global edge
      may only agree with some of them. This function returns true if the given global edge agrees
      with the orientation of its local counterpart of the adjacent face ilF. If jF and jlE are the
      equivalent global face and local edge numbers, then then is equivalent to ok_face_edge(jF, jlE). */
  bool ok_edge_face(
      size_t iE,             /**< The global edge number of the edge */
      size_t ilF             /**< The local face number of the adjacent face */
  ) const;

  /// Returns true if the global vertex iV is the first vertex on the global edge corresponding to ilE
  bool ok_vrtx_edge(
      size_t iV,             /**< The global vertex number of the vertex */
      size_t ilE             /**< The local edge number of the adjacent edge */
  ) const;

  // ...for regions

  /// Populate the given vector with the region numbers of all regions adjacent to the given region
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent region numbers will be appended
      to the back of the vector. A region is not considered adjacent to itself. */
  void get_regn_regn(
      size_t iR,                      /**< The global region number of the region */
      std::vector<size_t> &rlist     /**< A vector to populate with adjacent region numbers */
  ) const;
  /// Populate the given vector with the face numbers of all faces adjacent to the given region
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent face numbers will be appended
      to the back of the vector. */
  void get_regn_face(
      size_t iR,                      /**< The global region number of the region */
      std::vector<size_t> &flist     /**< A vector to populate with adjacent face numbers */
  ) const;
  /// Populate the given vector with the edge numbers of all edges adjacent to the given region
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent edge numbers will be appended
      to the back of the vector. */
  void get_regn_edge(
      size_t iR,                      /**< The global region number of the region */
      std::vector<size_t> &elist     /**< A vector to populate with adjacent edge numbers */
  ) const;
  /// Populate the given vector with the vertex numbers of all vertices adjacent to the given region
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent vertex numbers will be appended
      to the back of the vector. */
  void get_regn_vrtx(
      size_t iR,                      /**< The global region number of the region */
      std::vector<size_t> &vlist     /**< A vector to populate with adjacent vertex numbers */
  ) const;

  // ...for faces

  /// Populate the given vector with the region numbers of all regions adjacent to the given face
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent region numbers will be appended
      to the back of the vector. */
  void get_face_regn(
      size_t iF,                      /**< The global face number of the face */
      std::vector<size_t> &rlist     /**< A vector to populate with adjacent region numbers */
  ) const;
  /// Populate the given vector with the face numbers of all faces adjacent to the given face
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent face numbers will be appended
      to the back of the vector. A face is not considered adjacent to itself. */
  void get_face_face(
      size_t iF,                      /**< The global face number of the face */
      std::vector<size_t> &flist     /**< A vector to populate with adjacent face numbers */
  ) const;
  /// Populate the given vector with the edge numbers of all edges adjacent to the given face
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent edge numbers will be appended
      to the back of the vector. */
  void get_face_edge(
      size_t iF,                      /**< The global face number of the face */
      std::vector<size_t> &elist     /**< A vector to populate with adjacent edge numbers */
  ) const;
  /// Populate the given vector with the vertex numbers of all vertices adjacent to the given face
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent vertex numbers will be appended
      to the back of the vector. */
  void get_face_vrtx(
      size_t iF,                      /**< The global face number of the face */
      std::vector<size_t> &vlist     /**< A vector to populate with adjacent vertex numbers */
  ) const;

  // ...for edges

  /// Populate the given vector with the region numbers of all regions adjacent to the given edge
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent region numbers will be appended
      to the back of the vector. */
  void get_edge_regn(
      size_t iE,                      /**< The global edge number of the edge */
      std::vector<size_t> &rlist     /**< A vector to populate with adjacent region numbers */
  ) const;
  /// Populate the given vector with the face numbers of all faces adjacent to the given edge
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent face numbers will be appended
      to the back of the vector. */
  void get_edge_face(
      size_t iE,                      /**< The global edge number of the edge */
      std::vector<size_t> &flist     /**< A vector to populate with adjacent face numbers */
  ) const;
  /// Populate the given vector with the edge numbers of all edges adjacent to the given edge
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent edge numbers will be appended
      to the back of the vector. An edge is not considered adjacent to itself. */
  void get_edge_edge(
      size_t iE,                      /**< The global edge number of the edge */
      std::vector<size_t> &elist     /**< A vector to populate with adjacent edge numbers */
  ) const;
  /// Populate the given vector with the vertex numbers of all vertices adjacent to the given edge
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent vertex numbers will be appended
      to the back of the vector. */
  void get_edge_vrtx(
      size_t iE,                      /**< The global edge number of the edge */
      std::vector<size_t> &vlist     /**< A vector to populate with adjacent vertex numbers */
  ) const;

  // ...for vertices

  /// Populate the given vector with the region numbers of all regions adjacent to the given vertex
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent region numbers will be appended
      to the back of the vector. */
  void get_vrtx_regn(
      size_t iV,                      /**< The global vertex number of the vertex */
      std::vector<size_t> &rlist     /**< A vector to populate with adjacent region numbers */
  ) const;
  /// Populate the given vector with the face numbers of all faces adjacent to the given vertex
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent face numbers will be appended
      to the back of the vector. */
  void get_vrtx_face(
      size_t iV,                      /**< The global vertex number of the vertex */
      std::vector<size_t> &flist     /**< A vector to populate with adjacent face numbers */
  ) const;
  /// Populate the given vector with the edge numbers of all edges adjacent to the given vertex
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent edge numbers will be appended
      to the back of the vector. */
  void get_vrtx_edge(
      size_t iV,                      /**< The global vertex number of the vertex */
      std::vector<size_t> &elist     /**< A vector to populate with adjacent edge numbers */
  ) const;
  /// Populate the given vector with the vertex numbers of all vertices adjacent to the given vertex
  /** Note that if the given vector is non-empty, it will not be cleared. The adjacent vertex numbers will be appended
      to the back of the vector. A vertex is not considered adjacent to itself. */
  void get_vrtx_vrtx(
      size_t iV,                      /**< The global vertex number of the vertex */
      std::vector<size_t> &vlist     /**< A vector to populate with adjacent vertex numbers */
  ) const;

  // Logical Methods for detecting boundary items

  /// Returns true if the vertex with the given number is on the boundary
  bool is_boundary_vrtx(
      size_t iV              /**< The global vertex number of the vertex */
  ) const;

  /// Returns true if the edge with the given number is on the boundary
  bool is_boundary_edge(
      size_t iE              /**< The global edge number of the edge */
  ) const;

  /// Returns true if the face with the given number is on the boundary
  bool is_boundary_face(
      size_t iF              /**< The global face number of the face */
  ) const;

  /// Returns true if the region with the given number is on the boundary
  bool is_boundary_regn(
      size_t iR              /**< The global region number of the region */
  ) const;

  /// Returns true if the vertex with the given number is internal (not on the boundary)
  bool is_internal_vrtx(
      size_t iV              /**< The global vertex number of the vertex */
  ) const;

  /// Returns true if the edge with the given number is internal (not on the boundary)
  bool is_internal_edge(
      size_t iE              /**< The global edge number of the edge */
  ) const;

  /// Returns true if the face with the given number is internal (not on the boundary)
  bool is_internal_face(
      size_t iF              /**< The global face number of the face */
  ) const;

  /// Returns true if the region with the given number is internal (not on the boundary)
  bool is_internal_regn(
      size_t iR              /**< The global region number of the region */
  ) const;

  /// Returns the boundary vertex number of the given vertex
  bool get_bnd_pos_vrtx(
      size_t iV              /**< The global vertex number of the vertex */
  ) const;

  /// Returns the boundary edge number of the given edge
  bool get_bnd_pos_edge(
      size_t iE              /**< The global edge number of the edge */
  ) const;

  /// Returns the boundary face number of the given face
  bool get_bnd_pos_face(
      size_t iF              /**< The global face number of the face */
  ) const;

  /// Returns the boundary region number of the given region
  bool get_bnd_pos_regn(
      size_t iR              /**< The global region number of the region */
  ) const;

  // Geometrical Methods (new)

  /// Returns the s-coordinate of the unit normal vector to the given face
  /** Note that s=0,1,2 corresponds to the x,y,z coordinates respectively. The face normal is with respect to the
      global face's edge ordering, and hence may point in or out of a given region. To ensure that the normal is
      pointing out of a given cell, check if ok_regn_face(size_t iR, size_t ilF) is false, and if so, flip the normal. */
  double get_nor(
      size_t iF,             /**< The global face number of the face */
      size_t s               /**< The coordinate index (0,1,2) */
  ) const;
  /// Returns the s-coodinate of a unit tangent vector to a given edge
  /** Note that s=0,1,2 corresponds to the x,y,z coordinates respectively. The tangent points in the direction
      given by the global edge's vertex ordering. */
  double get_tng(
      size_t iE,             /**< The global edge number of the edge */
      size_t s               /**< The coordinate index (0,1,2) */
  ) const;

  // Geometrical measures (new)

  /// Returns the 3D measure (volume) of the region with the given number
  double get_regn_measure(
      size_t iR              /**< The global region number of the region */
  ) const;

  /// Returns the diameter of the region with the given number
  double get_regn_diam(
      size_t iR              /**< The global region number of the region */
  ) const;

  /// Returns the 2D measure (area) of the face with the given number
  double get_face_measure(
      size_t iF              /**< The global face number of the face */
  ) const;

  /// Returns the diameter face with the given number
  double get_face_diam(
      size_t iF              /**< The global face number of the face */
  ) const;

  /// Returns the 1D measure (length) of the edge with the given number
  double get_edge_measure(
      size_t iF              /**< The global edge number of the edge */
  ) const;

  /// The six extremes of the bounding box of the domain
  enum class bbox_dimension {
    xmin,
    ymin,
    zmin,
    xmax,
    ymax,
    zmax
  };

  /// Compute the bounding box of the domain and return the minimum/maximum value at the given boundary
  inline double eval_bbox(
      const bbox_dimension &ret       /**< The boundary limit to evaluate */
  ) const;

  /// Returns the minimum value in the s-coordinate in the domain
  inline double min_coords(
      size_t s               /**< The coordinate index (0,1,2) */
  ) const;

  /// Returns the maximum value in the s-coordinate in the domain
  inline double max_coords(
      size_t s               /**< The coordinate index (0,1,2) */
  ) const;

  /// Compute the bounding box of the domain and assign the boundaries to the given variables
  inline void bbox(
      double &xmin,      /**< Variable to assign the minimum x-value in the domain */
      double &ymin,      /**< Variable to assign the minimum y-value in the domain */
      double &zmin,      /**< Variable to assign the minimum z-value in the domain */
      double &xmax,      /**< Variable to assign the maximum x-value in the domain */
      double &ymax,      /**< Variable to assign the maximum y-value in the domain */
      double &zmax       /**< Variable to assign the maximum z-value in the domain */
  ) const;

  // used by problems
  /// Returns the minimum x-value in the domain
  inline double xmin() const;

  /// Returns the minimum y-value in the domain
  inline double ymin() const;

  /// Returns the minimum z-value in the domain
  inline double zmin() const;

  /// Returns the maximum x-value in the domain
  inline double xmax() const;

  /// Returns the maximum y-value in the domain
  inline double ymax() const;

  /// Returns the maximum z-value in the domain
  inline double zmax() const;

  // Geometrical Methods (with some problems)

  /// Returns the k-coordinate of vertex with the given number
  /** Note that s=0,1,2 corresponds to the x,y,z coordinates respectively. */
  double coords_V(
      size_t iV,             /**< The global vertex number of the vertex */
      size_t k               /**< The coordinate index (0,1,2) */
  ) const;
  /// Returns the k-coordinate of the barycenter of the given edge
  /** Note that s=0,1,2 corresponds to the x,y,z coordinates respectively. */
  double coords_E(
      size_t iE,             /**< The global vertex number of the edge */
      size_t k               /**< The coordinate index (0,1,2) */
  ) const;
  /// Returns the k-coordinate of the barycenter of the given face
  /** Note that s=0,1,2 corresponds to the x,y,z coordinates respectively. */
  double coords_F(
      size_t iF,             /**< The global vertex number of the face */
      size_t k               /**< The coordinate index (0,1,2) */
  ) const;
  /// Returns the k-coordinate of the barycenter of the given region
  /** Note that s=0,1,2 corresponds to the x,y,z coordinates respectively. */
  double coords_R(
      size_t iR,             /**< The global vertex number of the region */
      size_t k               /**< The coordinate index (0,1,2) */
  ) const;

  /// Returns k-coordinate of the arithmetic center of the face iF
  /** The arithmetic center is the average of the edge midpoints. It is not necessarily equal
      to the barycenter / center of mass for all polygons. If you need the barycenter, use coords_F(iF, k). */
  double ari_coords_F(
      size_t iF,             /**< The global face number of the face */
      size_t k               /**< The coordinate index (0,1,2) */
  ) const;

  // various notions of mesh size

  /// Returns the maximum diameter over all regions in the mesh
  double h_max() const;

  /// Returns the minimum edge length over all edges in the mesh
  double h_min() const;

  /// Returns the maximum cube-root of the volumes of all regions in the mesh
  double h_vol() const;

  /// Returns the cube-root of the average volume of all regions in the mesh
  double h_avg() const;

  // mesh name

  /// Sets the name of the mesh to the given string
  void set_mesh_name(
      std::string _mesh_name       /**< The new name to give the mesh */
  );

  /// Returns the name of the mesh
  /** If tex_flag is true, then underscore characters will be escaped with a \\ */
  std::string get_mesh_name(
      bool tex_flag = false     /**< If true, escape underscore characters for TeX */
  ) const;

  // set external flags (useful for gmv option "explode")

  /// Set the external flag for the given vertex
  /** Each vertex may have a flag associated with it for external use. Set these as your application requires. */
  void set_fV(
      size_t iV,                 /**< The global vertex number of the vertex */
      flag_type new_fV                 /**< The value of the flag to set at the given vertex */
  );
  /// Set the external flag for the given edge
  /** Each edge may have a flag associated with it for external use. Set these as your application requires. */
  void set_fE(
      size_t iE,                 /**< The global edge number of the edge */
      flag_type new_fE                 /**< The value of the flag to set at the given edge */
  );
  /// Set the external flag for the given face
  /** Each face may have a flag associated with it for external use. Set these as your application requires. */
  void set_fF(
      size_t iF,                 /**< The global face number of the face */
      flag_type new_fF                 /**< The value of the flag to set at the given face */
  );
  /// Set the external flag for the given vertex
  /** Each region may have a flag associated with it for external use. Set these as your application requires. */
  void set_fR(
      size_t iR,                 /**< The global region number of the region */
      flag_type new_fR                 /**< The value of the flag to set at the given region */
  );

  /// Returns the external flag of the given vertex
  /** Each vertex may have a flag associated with it for external use. Set these as your application requires. */
  flag_type get_fV(
      size_t iV                  /**< The global vertex number of the vertex */
  ) const;
  /// Returns the external flag of the given edge
  /** Each edge may have a flag associated with it for external use. Set these as your application requires. */
  flag_type get_fE(
      size_t iE                  /**< The global edge number of the edge */
  ) const;
  /// Returns the external flag of the given face
  /** Each face may have a flag associated with it for external use. Set these as your application requires. */
  flag_type get_fF(
      size_t iF                  /**< The global face number of the face */
  ) const;
  /// Returns the external flag of the given region
  /** Each region may have a flag associated with it for external use. Set these as your application requires. */
  flag_type get_fR(
      size_t iR                  /**< The global region number of the region */
  ) const;
};

// External reference to static member -- required for linkage
constexpr mesh_3Dv::flag_type mesh_3Dv::UNSET_FLAG;

// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------

template<typename T>
void mesh_3Dv::shrink_list(std::vector<T> &tmp_list) const {
  sort(tmp_list.begin(), tmp_list.end());
  tmp_list.erase(std::unique(std::begin(tmp_list), std::end(tmp_list)), std::end(tmp_list));
}

// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------

void mesh_3Dv::set_mesh_name(std::string _mesh_name) { mesh_name = std::move(_mesh_name); }

std::string mesh_3Dv::get_mesh_name(bool tex_flag) const {
  std::string retval;
  for (size_t i = 0; i < mesh_name.length(); ++i) {
    if (tex_flag && mesh_name[i] == char('_')) { retval += "\\"; }
    retval += mesh_name[i];
  }
  return retval;
}

// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------

// TOPOLOGICAL METHODS:
size_t mesh_3Dv::regn_face(size_t iR, size_t ilF) const {
  assert(iR < nR);
  assert(ilF < RegnFace[iR].size());
  return RegnFace[iR][ilF];
}

size_t mesh_3Dv::face_edge(size_t iF, size_t ilE) const {
  assert(iF < nF);
  assert(ilE < FaceEdge[iF].size());
  return FaceEdge[iF][ilE];
}

size_t mesh_3Dv::edge_vrtx(size_t iE, size_t ilV) const {
  assert(iE < nE);
  assert(ilV == 0 || ilV == 1);
  return EdgeVrtx[iE][ilV];
}

size_t mesh_3Dv::face_regn(size_t iF, size_t ilR) const {
  assert(iF < nF);
  assert(ilR < FaceRegn[iF].size());
  return FaceRegn[iF][ilR];
}

size_t mesh_3Dv::edge_face(size_t iE, size_t ilF) const {
  assert(iE < nE);
  assert(ilF < EdgeFace[iE].size());
  return EdgeFace[iE][ilF];
}

size_t mesh_3Dv::vrtx_edge(size_t iV, size_t ilE) const {
  assert(iV < nV);
  assert(ilE < VrtxEdge[iV].size());
  return VrtxEdge[iV][ilE];
}

size_t mesh_3Dv::vrtx_vrtx(size_t iV, size_t ilE) const {
  assert(iV < nV);
  size_t iE = VrtxEdge[iV][ilE];
  size_t iV0 = EdgeVrtx[iE][0];
  size_t iV1 = EdgeVrtx[iE][1];
  return iV0 == iV ? iV1 : iV0;
}

size_t mesh_3Dv::n_regn_face(size_t iR) const {
  assert(iR < nR);
  return RegnFace[iR].size();
}

size_t mesh_3Dv::n_face_edge(size_t iF) const {
  assert(iF < nF);
  return FaceEdge[iF].size();
}

size_t mesh_3Dv::n_edge_vrtx(size_t iE) const {
  assert(iE < nE);
  return EdgeVrtx[iE].size();
}

// vrtx_vrtx ISO vrtx_edge
size_t mesh_3Dv::n_vrtx_vrtx(size_t iV) const {
  assert(iV < nV);
  return VrtxEdge[iV].size();
}

size_t mesh_3Dv::n_face_regn(size_t iF) const {
  assert(iF < nF);
  return FaceRegn[iF].size();
}

size_t mesh_3Dv::n_edge_face(size_t iE) const {
  assert(iE < nE);
  return EdgeFace[iE].size();
}

size_t mesh_3Dv::n_vrtx_edge(size_t iV) const {
  assert(iV < nV);
  return VrtxEdge[iV].size();
}

bool mesh_3Dv::ok_regn_face(size_t iR, size_t ilF) const {
  assert(iR < nR);
  assert(ilF < RegnFace[iR].size());
  return RegnFaceOk[iR][ilF];
}

bool mesh_3Dv::ok_face_edge(size_t iF, size_t ilE) const {
  assert(iF < nF);
  assert(ilE < FaceEdge[iF].size());
  return FaceEdgeOk[iF][ilE];
}

bool mesh_3Dv::ok_edge_face(size_t iE, size_t ilF) const {
  assert(iE < nE);
  assert(ilF < EdgeFace[iE].size());
  return EdgeFaceOk[iE][ilF];
}

bool mesh_3Dv::ok_vrtx_edge(size_t iV, size_t ilE) const {
  assert(iV < nV);
  assert(ilE < VrtxEdge[iV].size());
  return VrtxEdge[iV][ilE];
}

// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
// access method of regions

void mesh_3Dv::get_regn_regn(size_t iR, std::vector<size_t> &rlist) const {
  assert(0 <= iR && iR < nR);
  for (size_t ilF = 0; ilF < RegnFace[iR].size(); ++ilF) {
    size_t iF = RegnFace[iR][ilF];
    for (size_t ilR = 0; ilR < n_face_regn(iF); ilR++) {
      size_t iR0 = face_regn(iF, ilR);
      if (iR0 != iR) {
        rlist.push_back(iR0);
      }
    }
  }
}

void mesh_3Dv::get_regn_face(size_t iR, std::vector<size_t> &flist) const {
  assert(0 <= iR && iR < nR);
  flist.resize(RegnFace[iR].size());
  for (size_t ilF = 0; ilF < RegnFace[iR].size(); ++ilF) { flist[ilF] = RegnFace[iR][ilF]; }
}

void mesh_3Dv::get_regn_edge(size_t iR, std::vector<size_t> &elist) const {
  assert(0 <= iR && iR < nR);
  std::vector<size_t> tmp_vec;
  for (size_t ilF = 0; ilF < RegnFace[iR].size(); ++ilF) {
    size_t iF = RegnFace[iR][ilF];
    for (size_t ilE = 0; ilE < FaceEdge[iF].size(); ++ilE) {
      tmp_vec.push_back(FaceEdge[iF][ilE]);
    }
  }
  shrink_list(tmp_vec);
  elist.resize(tmp_vec.size());
  for (size_t i = 0; i < elist.size(); ++i) { elist[i] = tmp_vec[i]; }
}

void mesh_3Dv::get_regn_vrtx(size_t iR, std::vector<size_t> &vlist) const {
  assert(0 <= iR && iR < nR);
  std::vector<size_t> tmp_vec;
  for (size_t ilF = 0; ilF < RegnFace[iR].size(); ++ilF) {
    size_t iF = RegnFace[iR][ilF];
    for (size_t ilE = 0; ilE < FaceEdge[iF].size(); ++ilE) {
      size_t iE = FaceEdge[iF][ilE];
      tmp_vec.push_back(EdgeVrtx[iE][0]);
      tmp_vec.push_back(EdgeVrtx[iE][1]);
    }
  }
  shrink_list(tmp_vec);
  vlist.resize(tmp_vec.size());
  for (size_t i = 0; i < vlist.size(); ++i) { vlist[i] = tmp_vec[i]; }
}
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
// access method of faces

void mesh_3Dv::get_face_regn(size_t iF, std::vector<size_t> &rlist) const {
  assert(0 <= iF && iF < nF);
  rlist.resize(0);
  for (size_t ilR = 0; ilR < n_face_regn(iF); ilR++) {
    rlist.push_back(face_regn(iF, ilR));
  }
}

void mesh_3Dv::get_face_face(size_t iF, std::vector<size_t> &flist) const {
  assert(0 <= iF && iF < nF);
  std::vector<size_t> tmp_vec;
  for (size_t ilE = 0; ilE < FaceEdge[iF].size(); ++ilE) {
    size_t iE = FaceEdge[iF][ilE];
    for (size_t ilF = 0; ilF < EdgeFace[iE].size(); ++ilF) {
      size_t jF = EdgeFace[iE][ilF];
      if (iF != jF) { tmp_vec.push_back(jF); }
    }
  }
  shrink_list(tmp_vec);
  flist.resize(tmp_vec.size());
  for (size_t i = 0; i < flist.size(); ++i) { flist[i] = tmp_vec[i]; }
}

void mesh_3Dv::get_face_edge(size_t iF, std::vector<size_t> &elist) const {
  assert(0 <= iF && iF < nF);
  elist.resize(FaceEdge[iF].size());
  for (size_t ilE = 0; ilE < FaceEdge[iF].size(); ++ilE) { elist[ilE] = FaceEdge[iF][ilE]; }
}

void mesh_3Dv::get_face_vrtx(size_t iF, std::vector<size_t> &vlist) const {
  assert(0 <= iF && iF < nF);
  std::vector<size_t> tmp_vec;
  for (size_t ilE = 0; ilE < FaceEdge[iF].size(); ++ilE) {
    size_t iE = FaceEdge[iF][ilE];
    if (FaceEdgeOk[iF][ilE]) {
      tmp_vec.push_back(EdgeVrtx[iE][0]);
    } else {
      tmp_vec.push_back(EdgeVrtx[iE][1]);
    }
  }
  vlist.resize(tmp_vec.size());
  for (size_t i = 0; i < vlist.size(); ++i) { vlist[i] = tmp_vec[i]; }
}

// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
// access method of edges

void mesh_3Dv::get_edge_regn(size_t iE, std::vector<size_t> &rlist) const {
  assert(0 <= iE && iE < nE);
  rlist.clear();
  for (size_t ilF = 0; ilF < EdgeFace[iE].size(); ++ilF) {
    size_t iF = EdgeFace[iE][ilF];
    for (size_t ilR = 0; ilR < n_face_regn(iF); ilR++) {
      rlist.push_back(face_regn(iF, ilR));
    }
  }
  shrink_list(rlist);
}

void mesh_3Dv::get_edge_face(size_t iE, std::vector<size_t> &flist) const {
  assert(0 <= iE && iE < nE);
  flist.resize(EdgeFace[iE].size());
  for (size_t ilF = 0; ilF < EdgeFace[iE].size(); ++ilF) { flist[ilF] = EdgeFace[iE][ilF]; }
}

void mesh_3Dv::get_edge_edge(size_t iE, std::vector<size_t> &elist) const {
  assert(0 <= iE && iE < nE);
  std::vector<size_t> tmp_vec;
  for (size_t ilV = 0; ilV < EdgeVrtx[iE].size(); ++ilV) {
    size_t iV = EdgeVrtx[iE][ilV];
    for (size_t ilE = 0; ilE < VrtxEdge[iV].size(); ++ilE) {
      size_t jE = VrtxEdge[iV][ilE];
      if (iE != jE) { tmp_vec.push_back(jE); }
    }
  }
  shrink_list(tmp_vec);
  elist.resize(tmp_vec.size());
  for (size_t i = 0; i < elist.size(); ++i) { elist[i] = tmp_vec[i]; }
}

void mesh_3Dv::get_edge_vrtx(size_t iE, std::vector<size_t> &vlist) const {
  assert(0 <= iE && iE < nE);
  vlist.resize(2);
  vlist[0] = EdgeVrtx[iE][0];
  vlist[1] = EdgeVrtx[iE][1];
}

// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
// access method of vertices

void mesh_3Dv::get_vrtx_regn(size_t iV, std::vector<size_t> &rlist) const {
  assert(0 <= iV && iV < nV);
  rlist.clear();
  for (size_t ilE = 0; ilE < VrtxEdge[iV].size(); ++ilE) {
    size_t iE = VrtxEdge[iV][ilE];
    for (size_t ilF = 0; ilF < EdgeFace[iE].size(); ++ilF) {
      size_t iF = EdgeFace[iE][ilF];
      for (size_t ilR = 0; ilR < n_face_regn(iF); ilR++) {
        rlist.push_back(face_regn(iF, ilR));
      }
    }
  }
  shrink_list(rlist);
}

void mesh_3Dv::get_vrtx_face(size_t iV, std::vector<size_t> &flist) const {
  assert(0 <= iV && iV < nV);
  std::vector<size_t> tmp_vec;
  for (size_t ilE = 0; ilE < VrtxEdge[iV].size(); ++ilE) {
    size_t iE = VrtxEdge[iV][ilE];
    for (size_t ilF = 0; ilF < EdgeFace[iE].size(); ++ilF) {
      tmp_vec.push_back(EdgeFace[iE][ilF]);
    }
  }
  shrink_list(tmp_vec);
  flist.resize(tmp_vec.size());
  for (size_t i = 0; i < flist.size(); ++i) { flist[i] = tmp_vec[i]; }
}

void mesh_3Dv::get_vrtx_edge(size_t iV, std::vector<size_t> &elist) const {
  assert(0 <= iV && iV < nV);
  elist.resize(VrtxEdge[iV].size());
  for (size_t ilE = 0; ilE < VrtxEdge[iV].size(); ++ilE) { elist[ilE] = VrtxEdge[iV][ilE]; }
}

void mesh_3Dv::get_vrtx_vrtx(size_t iV, std::vector<size_t> &vlist) const {
  assert(0 <= iV && iV < nV);
  vlist.resize(VrtxEdge[iV].size());
  for (size_t ilE = 0; ilE < VrtxEdge[iV].size(); ++ilE) {
    size_t iE = VrtxEdge[iV][ilE];
    size_t iV0 = EdgeVrtx[iE][0];
    size_t iV1 = EdgeVrtx[iE][1];
    vlist[ilE] = iV0 == iV ? iV1 : iV0;
  }
}

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// (new)

double mesh_3Dv::get_regn_measure(size_t iR) const {
  assert(0 <= iR && iR < nR);
  assert(R_volume.size() == nR);
  return R_volume[iR];
}

double mesh_3Dv::get_regn_diam(size_t iR) const {
  assert(0 <= iR && iR < nR);
  assert(R_diam.size() == nR);
  return R_diam[iR];
}

double mesh_3Dv::get_face_measure(size_t iF) const {
  assert(0 <= iF && iF < nF);
  assert(F_area.size() == nF);
  return F_area[iF];
}

double mesh_3Dv::get_face_diam(size_t iF) const {
  assert(0 <= iF && iF < nF);
  assert(F_diam.size() == nF);
  return F_diam[iF];
}

double mesh_3Dv::get_edge_measure(size_t iE) const {
  assert(0 <= iE && iE < nE);
  size_t iV0 = edge_vrtx(iE, 0);
  size_t iV1 = edge_vrtx(iE, 1);
  return std::sqrt(std::pow(coords_V(iV1, 0) - coords_V(iV0, 0), 2) +
                   std::pow(coords_V(iV1, 1) - coords_V(iV0, 1), 2) +
                   std::pow(coords_V(iV1, 2) - coords_V(iV0, 2), 2));
}

double mesh_3Dv::get_nor(size_t iF, size_t s) const {
  assert(0 <= iF && iF < nF);
  assert(0 <= s && s < DIM);
  assert(F_nor.size() == nF);
  return F_nor[iF][s];
}

double mesh_3Dv::get_tng(size_t iE, size_t s) const {
  assert(0 <= iE && iE < nE);
  assert(0 <= s && s < DIM);
  size_t iV0 = edge_vrtx(iE, 0);
  size_t iV1 = edge_vrtx(iE, 1);
  double len_E = get_edge_measure(iE);
  return (coords_V(iV1, s) - coords_V(iV0, s)) / len_E;
}

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

size_t mesh_3Dv::get_bnd_regn(size_t ilR) const {
  assert(0 <= ilR && ilR < bnd_regn.size());
  return bnd_regn[ilR];
}

size_t mesh_3Dv::get_bnd_face(size_t ilF) const {
  assert(0 <= ilF && ilF < bnd_face.size());
  return bnd_face[ilF];
}

size_t mesh_3Dv::get_bnd_edge(size_t ilE) const {
  assert(0 <= ilE && ilE < bnd_edge.size());
  return bnd_edge[ilE];
}

size_t mesh_3Dv::get_bnd_vrtx(size_t ilV) const {
  assert(0 <= ilV && ilV < bnd_vrtx.size());
  return bnd_vrtx[ilV];
}

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

bool mesh_3Dv::is_boundary_vrtx(size_t iV) const {
  return std::binary_search(std::begin(bnd_vrtx), std::end(bnd_vrtx), iV);
}

bool mesh_3Dv::is_boundary_edge(size_t iE) const {
  return std::binary_search(std::begin(bnd_edge), std::end(bnd_edge), iE);
}

bool mesh_3Dv::is_boundary_face(size_t iF) const {
  return std::binary_search(std::begin(bnd_face), std::end(bnd_face), iF);
}

bool mesh_3Dv::is_boundary_regn(size_t iR) const {
  return std::binary_search(std::begin(bnd_regn), std::end(bnd_regn), iR);
}

bool mesh_3Dv::is_internal_vrtx(size_t iV) const {
  return !std::binary_search(std::begin(bnd_vrtx), std::end(bnd_vrtx), iV);
}

bool mesh_3Dv::is_internal_edge(size_t iE) const {
  return !std::binary_search(std::begin(bnd_edge), std::end(bnd_edge), iE);
}

bool mesh_3Dv::is_internal_face(size_t iF) const {
  return !std::binary_search(std::begin(bnd_face), std::end(bnd_face), iF);
}

bool mesh_3Dv::is_internal_regn(size_t iR) const {
  return !std::binary_search(std::begin(bnd_regn), std::end(bnd_regn), iR);
}

bool mesh_3Dv::get_bnd_pos_vrtx(size_t iV) const {
  auto it = std::lower_bound(std::begin(bnd_vrtx), std::end(bnd_vrtx), iV);
  assert(it != std::end(bnd_vrtx));
  return std::distance(std::begin(bnd_vrtx), it);
}

bool mesh_3Dv::get_bnd_pos_edge(size_t iE) const {
  auto it = std::lower_bound(std::begin(bnd_edge), std::end(bnd_edge), iE);
  assert(it != std::end(bnd_edge));
  return std::distance(std::begin(bnd_edge), it);
}

bool mesh_3Dv::get_bnd_pos_face(size_t iF) const {
  auto it = std::lower_bound(std::begin(bnd_face), std::end(bnd_face), iF);
  assert(it != std::end(bnd_face));
  return std::distance(std::begin(bnd_face), it);
}

bool mesh_3Dv::get_bnd_pos_regn(size_t iR) const {
  auto it = std::lower_bound(std::begin(bnd_regn), std::end(bnd_regn), iR);
  assert(it != std::end(bnd_regn));
  return std::distance(std::begin(bnd_regn), it);
}

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
double mesh_3Dv::coords_V(size_t iV, size_t k) const {
  assert(0 <= iV && iV < nV);
  assert(0 <= k && k < 3);
  return V_coords[iV][k];
}

double mesh_3Dv::coords_E(size_t iE, size_t k) const {
  assert(0 <= iE && iE < nE);
  assert(0 <= k && k < 3);
  double retval = 0.;
  for (size_t ilV = 0; ilV < EdgeVrtx[iE].size(); ++ilV) {
    size_t iV = EdgeVrtx[iE][ilV];
    retval += coords_V(iV, k);
  }
  return retval / double(EdgeVrtx[iE].size());
}

// this implementation is exact only for constants
double mesh_3Dv::ari_coords_F(size_t iF, size_t k) const {
  assert(0 <= iF && iF < nF);
  assert(0 <= k && k < 3);
  double retval = 0.;
  for (size_t ilE = 0; ilE < FaceEdge[iF].size(); ++ilE) {
    size_t iE = FaceEdge[iF][ilE];
    retval += coords_E(iE, k);
  }
  return retval / double(FaceEdge[iF].size());
}

double mesh_3Dv::coords_F(size_t iF, size_t k) const {
  assert(0 <= iF && iF < nF);
  assert(0 <= k && k < 3);
  return F_coords[iF][k];
  //return ari_coords_F( iF, k ) ; // DEBUG
}

double mesh_3Dv::coords_R(size_t iR, size_t k) const {
  assert(0 <= iR && iR < nR);
  assert(0 <= k && k < 3);
  return R_coords[iR][k];
}

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
void mesh_3Dv::eval_bbox() const {
  bb_status = true;
  for (size_t s = 0; s < DIM; ++s) {
    bb_min[s] = +1e+99;
    bb_max[s] = -1e+99;
  }
  for (size_t ilV = 0; ilV < n_bvertex(); ++ilV) {
    size_t iV = bnd_vrtx[ilV];
    for (size_t s = 0; s < DIM; ++s) {
      bb_min[s] = std::min(bb_min[s], V_coords[iV][s]);
      bb_max[s] = std::max(bb_max[s], V_coords[iV][s]);
    }
  }
}

double mesh_3Dv::eval_bbox(const bbox_dimension &ret) const {
  bb_status = true;
  for (size_t s = 0; s < DIM; ++s) {
    bb_min[s] = +1e+99;
    bb_max[s] = -1e+99;
  }
  for (size_t ilV = 0; ilV < n_bvertex(); ++ilV) {
    size_t iV = bnd_vrtx[ilV];
    for (size_t s = 0; s < DIM; ++s) {
      bb_min[s] = std::min(bb_min[s], V_coords[iV][s]);
      bb_max[s] = std::max(bb_max[s], V_coords[iV][s]);
    }
  }
  double retval;
  if (ret == bbox_dimension::xmin) { retval = bb_min[0]; }
  else if (ret == bbox_dimension::ymin) { retval = bb_min[1]; }
  else if (ret == bbox_dimension::zmin) { retval = bb_min[2]; }
  else if (ret == bbox_dimension::xmax) { retval = bb_max[0]; }
  else if (ret == bbox_dimension::ymax) { retval = bb_max[1]; }
  else if (ret == bbox_dimension::zmax) { retval = bb_max[2]; }
  else { retval = 0.0; }
  return retval;
}

void mesh_3Dv::bbox(double &xmin, double &ymin, double &zmin,
                    double &xmax, double &ymax, double &zmax) const {
  if (!bb_status) { eval_bbox(); }
  xmin = bb_min[0];
  ymin = bb_min[1];
  zmin = bb_min[2];
  xmax = bb_max[0];
  ymax = bb_max[1];
  zmax = bb_max[2];
}

inline double mesh_3Dv::min_coords(size_t s) const {
  assert(0 <= s && s < DIM);
  if (!bb_status) { eval_bbox(); }
  return bb_min[s];
}

inline double mesh_3Dv::max_coords(size_t s) const {
  assert(0 <= s && s < DIM);
  if (!bb_status) { eval_bbox(); }
  return bb_max[s];
}

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

// used by problems
inline double mesh_3Dv::xmin() const { return bb_status ? bb_min[0] : eval_bbox(bbox_dimension::xmin); }

inline double mesh_3Dv::ymin() const { return bb_status ? bb_min[1] : eval_bbox(bbox_dimension::ymin); }

inline double mesh_3Dv::zmin() const { return bb_status ? bb_min[2] : eval_bbox(bbox_dimension::zmin); }

inline double mesh_3Dv::xmax() const { return bb_status ? bb_max[0] : eval_bbox(bbox_dimension::xmax); }

inline double mesh_3Dv::ymax() const { return bb_status ? bb_max[1] : eval_bbox(bbox_dimension::ymax); }

inline double mesh_3Dv::zmax() const { return bb_status ? bb_max[2] : eval_bbox(bbox_dimension::zmax); }

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
void mesh_3Dv::eval_h() const {
  if (!ms_status) {
    ms_status = true;
    // hmin
    ms_hmin = 1.e+99;
    for (size_t iE = 0; iE < nE; ++iE) {
      size_t iV0 = EdgeVrtx[iE][0];
      size_t iV1 = EdgeVrtx[iE][1];
      double len =
          sqrt(pow(V_coords[iV0][0] - V_coords[iV1][0], 2) +
               pow(V_coords[iV0][1] - V_coords[iV1][1], 2) +
               pow(V_coords[iV0][2] - V_coords[iV1][2], 2));
      ms_hmax = std::max(len, ms_hmax);
      ms_hmin = std::min(len, ms_hmin);
    }

    // evaluate max diameter of a polyhedron (from a Jerome suggestion)
    ms_hmax = 0.;
    for (size_t iR = 0; iR < nR; ++iR) {
			ms_hmax = std::max(ms_hmax, get_regn_diam(iR));
    }

    // hvol
    ms_hvol = 0.;
    for (size_t iR = 0; iR < nR; ++iR) {
      double hvol = pow(R_volume[iR], 1. / 3.);
      ms_hvol = std::max(hvol, ms_hvol);
    }
    // havg
    double sum = 0.;
    for (size_t iR = 0; iR < nR; ++iR) {
      sum += R_volume[iR];
    }
    ms_havg = pow(sum / nR, 1. / double(DIM));
  }
}

double mesh_3Dv::h_max() const {
  if (!ms_status) { eval_h(); }
  return ms_hmax;
}

double mesh_3Dv::h_min() const {
  if (!ms_status) { eval_h(); }
  return ms_hmin;
}

double mesh_3Dv::h_vol() const {
  if (!ms_status) { eval_h(); }
  return ms_hvol;
}

double mesh_3Dv::h_avg() const {
  if (!ms_status) { eval_h(); }
  return ms_havg;
}
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

// set external flags (useful for gmv option "explode")
void mesh_3Dv::set_fV(size_t iV, flag_type new_fV) { fV[iV] = new_fV; }

void mesh_3Dv::set_fE(size_t iE, flag_type new_fE) { fE[iE] = new_fE; }

void mesh_3Dv::set_fF(size_t iF, flag_type new_fF) { fF[iF] = new_fF; }

void mesh_3Dv::set_fR(size_t iR, flag_type new_fR) { fR[iR] = new_fR; }

mesh_3Dv::flag_type mesh_3Dv::get_fV(size_t iV) const { return fV[iV]; }

mesh_3Dv::flag_type mesh_3Dv::get_fE(size_t iE) const { return fE[iE]; }

mesh_3Dv::flag_type mesh_3Dv::get_fF(size_t iF) const { return fF[iF]; }

mesh_3Dv::flag_type mesh_3Dv::get_fR(size_t iR) const { return fR[iR]; }


// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

}  // end of namespace StemMesh3D

#endif // end of _MESH3DV_HH
