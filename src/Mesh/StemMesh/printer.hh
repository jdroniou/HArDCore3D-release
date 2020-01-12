#ifndef _MESH3D_PRINTER_HH
#define _MESH3D_PRINTER_HH

#include <iostream>
#include <fstream>
#include <vector>

namespace StemMesh3D {

class mesh3Dv_printer {

private:
  mesh_3Dv &mesh;
  size_t offset;

public:
  mesh3Dv_printer(mesh_3Dv &_mesh, size_t _offset = 0) : mesh(_mesh), offset(_offset) {}

  ~mesh3Dv_printer() {}

  // print DATASETS
  void print_EdgeVrtx(std::ostream &LOGF);

  void print_FaceEdge(std::ostream &LOGF);

  void print_RegnFace(std::ostream &LOGF);

  // ---
  void print_FaceRegn(std::ostream &LOGF);

  void print_EdgeFace(std::ostream &LOGF);

  void print_VrtxEdge(std::ostream &LOGF);

  // ---
  void print_boundary_lists(std::ostream &LOGF);

  // ---
  void print_all_datasets();

  // print REGIONS
  void print_all_regions();

  void print_region(size_t iR, std::ostream &LOGF);

  // print FACES
  void print_all_faces();

  void print_face(size_t iF, std::ostream &LOGF);

  // print EDGES
  void print_all_edges();

  void print_edge(size_t iE, std::ostream &LOGF);

  // print VERTICES
  void print_all_vertices();

  void print_vertex(size_t iV, std::ostream &LOGF);
};

// --------------------------------------------------------------------------------------------
void mesh3Dv_printer::print_EdgeVrtx(std::ostream &LOGF) {
  // EdgeVrtx
  LOGF << "EdgeVrtx " << std::endl;
  LOGF << "#of vertices = " << mesh.n_edge() << std::endl;
  for (size_t iE = 0; iE < mesh.n_edge(); ++iE) {
    LOGF << "iE = " << iE + offset << ", { ";
    if (mesh.n_edge_vrtx(iE) > 0) {
      for (size_t ilV = 0; ilV < mesh.n_edge_vrtx(iE) - 1; ++ilV) {
        LOGF << mesh.edge_vrtx(iE, ilV) + offset << ", ";
      }
      LOGF << mesh.edge_vrtx(iE, mesh.n_edge_vrtx(iE) - 1) + offset;
    }
    LOGF << " } " << std::endl;
  }
  LOGF << "-------------------" << std::endl;
}

void mesh3Dv_printer::print_FaceEdge(std::ostream &LOGF) {
  // FaceEdge
  LOGF << "FaceEdge " << std::endl;
  LOGF << "#of faces = " << mesh.n_face() << std::endl;
  for (size_t iF = 0; iF < mesh.n_face(); ++iF) {
    LOGF << "iF = " << iF + offset << ", { ";
    if (mesh.n_face_edge(iF) > 0) {
      for (size_t j = 0; j < mesh.n_face_edge(iF) - 1; ++j) {
        LOGF << "(" << mesh.ok_face_edge(iF, j) << ")" << mesh.face_edge(iF, j) + offset << ", ";
      }
      size_t j = mesh.n_face_edge(iF) - 1;
      LOGF << "(" << mesh.ok_face_edge(iF, j) << ")"
           << mesh.face_edge(iF, j) + offset;
    }
    LOGF << " } " << std::endl;
  }
  LOGF << "-------------------" << std::endl;
}

void mesh3Dv_printer::print_RegnFace(std::ostream &LOGF) {
  // RegnFace
  LOGF << "RegnFace " << std::endl;
  LOGF << "#of regns = " << mesh.n_region() << std::endl;
  for (size_t iR = 0; iR < mesh.n_region(); ++iR) {
    LOGF << "iR = " << iR + offset << ", { ";
    if (mesh.n_regn_face(iR) > 0) {
      for (size_t j = 0; j < mesh.n_regn_face(iR) - 1; ++j) {
        LOGF << "(" << mesh.ok_regn_face(iR, j) << ")" << mesh.regn_face(iR, j) + offset << ", ";
      }
      size_t j = mesh.n_regn_face(iR) - 1;
      LOGF << "(" << mesh.ok_regn_face(iR, j) << ")"
           << mesh.regn_face(iR, j) + offset;
    }
    LOGF << " } " << std::endl;
  }
  LOGF << "-------------------" << std::endl;
}

// ------------------------------------------------------------------------------------------
void mesh3Dv_printer::print_VrtxEdge(std::ostream &LOGF) {
  // VrtxEdge
  LOGF << "VrtxEdge " << std::endl;
  LOGF << "#of vertices = " << mesh.n_vertex() << std::endl;
  for (size_t iV = 0; iV < mesh.n_vertex(); ++iV) {
    LOGF << "iV = " << iV + offset << ", { ";
    if (mesh.n_vrtx_edge(iV) > 0) {
      for (size_t j = 0; j < mesh.n_vrtx_edge(iV) - 1; ++j) {
        LOGF << "(" << mesh.ok_vrtx_edge(iV, j) << ")" << mesh.vrtx_edge(iV, j) + offset << ", ";
      }
      size_t j = mesh.n_vrtx_edge(iV) - 1;
      LOGF << "(" << mesh.ok_vrtx_edge(iV, j) << ")"
           << mesh.vrtx_edge(iV, j) + offset;
    }
    LOGF << " } " << std::endl;
  }
  LOGF << "-------------------" << std::endl;
}

void mesh3Dv_printer::print_EdgeFace(std::ostream &LOGF) {
  // EdgeFace
  LOGF << "EdgeFace " << std::endl;
  LOGF << "#of edges = " << mesh.n_edge() << std::endl;
  for (size_t iE = 0; iE < mesh.n_edge(); ++iE) {
    LOGF << "iE = " << iE + offset << ", { ";
    if (mesh.n_edge_face(iE) > 0) {
      for (size_t j = 0; j < mesh.n_edge_face(iE) - 1; ++j) {
        LOGF << "(" << mesh.ok_edge_face(iE, j) << ")" << mesh.edge_face(iE, j) + offset << ", ";
      }
      size_t j = mesh.n_edge_face(iE) - 1;
      LOGF << "(" << mesh.ok_edge_face(iE, j) << ")"
           << mesh.edge_face(iE, j) + offset;
    }
    LOGF << " } " << std::endl;
  }
  LOGF << "-------------------" << std::endl;
}

void mesh3Dv_printer::print_FaceRegn(std::ostream &LOGF) {
  // EdgeFace
  LOGF << "FaceRegn " << std::endl;
  LOGF << "#of faces = " << mesh.n_face() << std::endl;
  for (size_t iF = 0; iF < mesh.n_regn_face(iF); ++iF) {
    LOGF << "iF = " << iF + offset << ", { "
         << mesh.face_regn(iF, 0) + offset;
    if (mesh.n_face_regn(iF) == 2)
      LOGF << ", " << mesh.face_regn(iF, 1) + offset;
  }
  LOGF << " } " << std::endl;
  LOGF << "-------------------" << std::endl;
}

void mesh3Dv_printer::print_boundary_lists(std::ostream &LOGF) {
  // boundary lists
  LOGF << "boundary lists " << std::endl;
  LOGF << "#of boundary vertices = " << mesh.n_bvertex() << std::endl;
  for (size_t i = 0; i < mesh.n_bvertex(); ++i) {
    LOGF << "bnd_vrtx[" << i << "] = " << mesh.get_bnd_vrtx(i) + offset << std::endl;
  }
  LOGF << "-------------------" << std::endl;
  LOGF << "#of boundary edges = " << mesh.n_bedge() << std::endl;
  for (size_t i = 0; i < mesh.n_bedge(); ++i) {
    LOGF << "bnd_edge[" << i << "] = " << mesh.get_bnd_edge(i) + offset << std::endl;
  }
  LOGF << "-------------------" << std::endl;
  LOGF << "#of boundary faces = " << mesh.n_bface() << std::endl;
  for (size_t i = 0; i < mesh.n_bface(); ++i) {
    LOGF << "bnd_face[" << i << "] = " << mesh.get_bnd_face(i) + offset << std::endl;
  }
  LOGF << "-------------------" << std::endl;
  LOGF << "#of boundary regions = " << mesh.n_bregion() << std::endl;
  for (size_t i = 0; i < mesh.n_bregion(); ++i) {
    LOGF << "bnd_regn[" << i << "] = " << mesh.get_bnd_regn(i) + offset << std::endl;
  }
  LOGF << "-------------------" << std::endl;
}

// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
void mesh3Dv_printer::print_all_datasets() {
  std::ofstream LOGF("dset.log");
  print_EdgeVrtx(LOGF);
  print_FaceEdge(LOGF);
  print_RegnFace(LOGF);
  print_EdgeFace(LOGF);
  print_VrtxEdge(LOGF);
  print_FaceRegn(LOGF);
  LOGF << "-------------------" << std::endl;
  print_boundary_lists(LOGF);
  LOGF.close();
}
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
# define print_list(VEC)                                     \
  do {                                                       \
    LOGF << " { " ;               \
    if ( VEC.size()>0 ) {             \
      for ( size_t i=0 ; i<VEC.size()-1 ; ++i ) {         \
        LOGF << VEC[i] + offset << ", " ;           \
      }                   \
      LOGF << VEC.back() + offset ;             \
    }                   \
    LOGF << " }" << std::endl ;             \
  } while(0)

// --------------------------------------------------------------------------------------------
// REGIONS
// --------------------------------------------------------------------------------------------
void mesh3Dv_printer::print_all_regions() {
  std::ofstream LOGF("region.log");
  LOGF << "#of regions: nR=" << mesh.n_region() << std::endl;
  for (size_t iR = 0; iR < mesh.n_region(); ++iR) {
    print_region(iR, LOGF);
  }
  LOGF.close();
}

void mesh3Dv_printer::print_region(size_t iR, std::ostream &LOGF) {
  LOGF << "--------------------------------------------------------------" << std::endl;
  LOGF << "region: iR =" << iR + offset << std::endl;
  // faces
  size_t nRF = mesh.n_regn_face(iR);
  LOGF << "#faces: nRF=" << nRF << ", { ";
  for (size_t ilF = 0; ilF < nRF - 1; ++ilF) {
    LOGF << "("
         << mesh.ok_regn_face(iR, ilF) << ")"
         << mesh.regn_face(iR, ilF) + offset << ", ";
  }
  LOGF << "("
       << mesh.ok_regn_face(iR, nRF - 1) << ")"
       << mesh.regn_face(iR, nRF - 1) + offset << " }" << std::endl;
  // topological info
  std::vector<size_t> rlist, flist, elist, vlist;
  mesh.get_regn_regn(iR, rlist);
  mesh.get_regn_face(iR, flist);
  mesh.get_regn_edge(iR, elist);
  mesh.get_regn_vrtx(iR, vlist);
  LOGF << "---" << std::endl;
  LOGF << "#regions:  nRR=" << rlist.size();
  print_list(rlist);
  LOGF << "#faces:    nRF=" << flist.size();
  print_list(flist);
  LOGF << "#edges:    nRE=" << elist.size();
  print_list(elist);
  LOGF << "#vertices: nRV=" << vlist.size();
  print_list(vlist);
  // geometrical info
  LOGF << "---" << std::endl;
  LOGF << "#center: ( "
       << mesh.coords_R(iR, 0) << ", " << mesh.coords_R(iR, 1) << ", "
       << mesh.coords_R(iR, 2) << " )" << std::endl;
  LOGF << "#volume:     " << mesh.get_regn_measure(iR) << std::endl;
}

// --------------------------------------------------------------------------------------------
// FACES
// --------------------------------------------------------------------------------------------
void mesh3Dv_printer::print_all_faces() {
  std::ofstream LOGF("face.log");
  LOGF << "#of faces: nF=" << mesh.n_face() << std::endl;
  for (size_t iF = 0; iF < mesh.n_face(); ++iF) {
    print_face(iF, LOGF);
  }
  LOGF.close();
}

void mesh3Dv_printer::print_face(size_t iF, std::ostream &LOGF) {
  LOGF << "--------------------------------------------------------------" << std::endl;
  LOGF << "face: iF =" << iF + offset << std::endl;
  // faces
  size_t nFE = mesh.n_face_edge(iF);
  LOGF << "#edges: nFE=" << nFE << ", { ";
  for (size_t ilE = 0; ilE < nFE - 1; ++ilE) {
    LOGF << "("
         << mesh.ok_face_edge(iF, ilE) << ")"
         << mesh.face_edge(iF, ilE) + offset << ", ";
  }
  LOGF << "("
       << mesh.ok_face_edge(iF, nFE - 1) << ")"
       << mesh.face_edge(iF, nFE - 1) + offset << " }" << std::endl;
  // topological info
  std::vector<size_t> rlist, flist, elist, vlist;
  mesh.get_face_regn(iF, rlist);
  mesh.get_face_face(iF, flist);
  mesh.get_face_edge(iF, elist);
  mesh.get_face_vrtx(iF, vlist);
  LOGF << "---" << std::endl;
  LOGF << "#regions:  nFR=" << rlist.size();
  print_list(rlist);
  LOGF << "#faces:    nFF=" << flist.size();
  print_list(flist);
  LOGF << "#edges:    nFE=" << elist.size();
  print_list(elist);
  LOGF << "#vertices: nFV=" << vlist.size();
  print_list(vlist);
  // geometrical info
  LOGF << "---" << std::endl;
  LOGF << "#center: ( "
       << mesh.coords_F(iF, 0) << ", "
       << mesh.coords_F(iF, 1) << ", "
       << mesh.coords_F(iF, 2) << " )" << std::endl;
  LOGF << "#normal: ( "
       << mesh.get_nor(iF, 0) << ", "
       << mesh.get_nor(iF, 1) << ", "
       << mesh.get_nor(iF, 2) << " )" << std::endl;
  LOGF << "#area:     " << mesh.get_face_measure(iF) << std::endl;
}

// --------------------------------------------------------------------------------------------
// EDGES
// --------------------------------------------------------------------------------------------
void mesh3Dv_printer::print_all_edges() {
  std::ofstream LOGF("edge.log");
  LOGF << "#of edges: nE=" << mesh.n_edge() << std::endl;
  for (size_t iE = 0; iE < mesh.n_edge(); ++iE) {
    print_edge(iE, LOGF);
  }
  LOGF.close();
}

void mesh3Dv_printer::print_edge(size_t iE, std::ostream &LOGF) {
  LOGF << "--------------------------------------------------------------" << std::endl;
  LOGF << "edge: iE =" << iE + offset << std::endl;
  // edges
  size_t nEV = 2;
  size_t iV0 = mesh.edge_vrtx(iE, 0);
  size_t iV1 = mesh.edge_vrtx(iE, 1);
  LOGF << "#vertices: nVE=" << nEV << ", { ";
  LOGF << mesh.edge_vrtx(iE, 0) + offset << ", "
       << mesh.edge_vrtx(iE, 1) + offset << " }" << std::endl;
  // topological info
  std::vector<size_t> rlist, flist, elist, vlist;
  mesh.get_edge_regn(iE, rlist);
  mesh.get_edge_face(iE, flist);
  mesh.get_edge_edge(iE, elist);
  mesh.get_edge_vrtx(iE, vlist);
  LOGF << "---" << std::endl;
  LOGF << "#regions:  nER=" << rlist.size();
  print_list(rlist);
  LOGF << "#faces:    nEF=" << flist.size();
  print_list(flist);
  LOGF << "#edges:    nEE=" << elist.size();
  print_list(elist);
  LOGF << "#vertices: nEV=" << vlist.size();
  print_list(vlist);
  // geometrical info
  LOGF << "---" << std::endl;
  LOGF << "#center:  ( "
       << mesh.coords_E(iE, 0) << ", "
       << mesh.coords_E(iE, 1) << ", "
       << mesh.coords_E(iE, 2) << " )" << std::endl;
  LOGF << "#tangent: ( "
       << mesh.get_tng(iE, 0) << ", "
       << mesh.get_tng(iE, 1) << ", "
       << mesh.get_tng(iE, 2) << " )" << std::endl;
  LOGF << "#length:     " << mesh.get_edge_measure(iE) << std::endl;
  LOGF << "#vrtx 0:  ( "
       << mesh.coords_V(iV0, 0) << ", "
       << mesh.coords_V(iV0, 1) << ", "
       << mesh.coords_V(iV0, 2) << " )" << std::endl;
  LOGF << "#vrtx 1:  ( "
       << mesh.coords_V(iV1, 0) << ", "
       << mesh.coords_V(iV1, 1) << ", "
       << mesh.coords_V(iV1, 2) << " )" << std::endl;
}

// --------------------------------------------------------------------------------------------
// VERTICES
// --------------------------------------------------------------------------------------------
void mesh3Dv_printer::print_all_vertices() {
  std::ofstream LOGF("vertex.log");
  LOGF << "#of vertices: nV=" << mesh.n_vertex() << std::endl;
  for (size_t iV = 0; iV < mesh.n_vertex(); ++iV) {
    print_vertex(iV, LOGF);
  }
  LOGF.close();
}

void mesh3Dv_printer::print_vertex(size_t iV, std::ostream &LOGF) {
  LOGF << "--------------------------------------------------------------" << std::endl;
  LOGF << "vertex: iV =" << iV + offset << std::endl;
  // coordinates
  LOGF << "#coords:= ( "
       << mesh.coords_V(iV, 0) << ", "
       << mesh.coords_V(iV, 1) << ", "
       << mesh.coords_V(iV, 2) << ") " << std::endl;
  // topological info
  std::vector<size_t> rlist, flist, elist, vlist;
  mesh.get_vrtx_regn(iV, rlist);
  mesh.get_vrtx_face(iV, flist);
  mesh.get_vrtx_edge(iV, elist);
  mesh.get_vrtx_vrtx(iV, vlist);
  LOGF << "---" << std::endl;
  LOGF << "#regions:  nVR=" << rlist.size();
  print_list(rlist);
  LOGF << "#faces:    nVF=" << flist.size();
  print_list(flist);
  LOGF << "#edges:    nVE=" << elist.size();
  print_list(elist);
  LOGF << "#vertices: nVV=" << vlist.size();
  print_list(vlist);
}

} // end of namespace StemMesh3D

#endif // _MESH3D_PRINTER_HH
