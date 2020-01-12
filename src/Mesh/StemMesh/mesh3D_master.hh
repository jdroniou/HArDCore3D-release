#ifndef _MESH3D_MASTER_HH
#define _MESH3D_MASTER_HH

#include <algorithm>
#include <vector>
#include <iostream>

#include "mesh3D.hh"

namespace StemMesh3D {

class mesh3Dv_cmaster {

private:
  static const int DIM = 3;
  const double TOL;

private:
  mesh_3Dv &mesh;
  bool mesh_status;
  bool regn_status;
  bool face_status;
  bool edge_status;
  bool vrtx_status;
  bool ornt_status;

public:
  mesh3Dv_cmaster(mesh_3Dv &_mesh) :
      TOL(1.e-14),
      mesh(_mesh),
      mesh_status(true),
      regn_status(true),
      face_status(true),
      edge_status(true),
      vrtx_status(true),
      ornt_status(true) {}

  // check the mesh
  void check_the_mesh();

  // check the regions
  bool check_the_regions();

  bool check_regn_regn(size_t i);

  bool check_regn_face(size_t i);

  bool check_regn_edge(size_t i);

  bool check_regn_vrtx(size_t i);

  // check the faces
  bool check_the_faces();

  bool check_face_regn(size_t i);

  bool check_face_face(size_t i);

  bool check_face_edge(size_t i);

  bool check_face_vrtx(size_t i);

  // check the edges
  bool check_the_edges();

  bool check_edge_regn(size_t i);

  bool check_edge_face(size_t i);

  bool check_edge_edge(size_t i);

  bool check_edge_vrtx(size_t i);

  // check the vertices
  bool check_the_vertices();

  bool check_vrtx_regn(size_t i);

  bool check_vrtx_face(size_t i);

  bool check_vrtx_edge(size_t i);

  bool check_vrtx_vrtx(size_t i);

  // ------------------
  //void check_this_item( size_t iE ) ;
  // ------------------
  bool check_local_regions();

  bool check_local_faces();

  bool check_local_edges();

  bool check_local_vertices();

  // check orientation and local numberings
  bool check_regn_face_orientation();

  bool check_regn_face_numbering();

  bool check_face_edge_numbering();

  bool check_vrtx_vrtx_numbering();

  bool check_edge_regn_numbering();

  // check some geometric quantities
  bool check_geometric_region_factor();

  bool check_region_face_orientation();
};

void mesh3Dv_cmaster::check_the_mesh() {
  std::cout << "-------------------------------------------------------" << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
  std::cout << "BEGIN CHECKING THE MESH" << std::endl;
  if (regn_status) { mesh_status &= check_the_regions(); }
  if (face_status) { mesh_status &= check_the_faces(); }
  if (edge_status) { mesh_status &= check_the_edges(); }
  if (vrtx_status) { mesh_status &= check_the_vertices(); }
  if (mesh_status) {
    std::cout << "Mesh datasets are ok!" << std::endl << std::endl;
  } else {
    std::cout << "-->>>THERE ARE TROUBLES IN MESH DATASETS!!!" << std::endl << std::endl;
  }
  if (regn_status) { mesh_status &= check_local_regions(); }
  if (face_status) { mesh_status &= check_local_faces(); }
  if (edge_status) { mesh_status &= check_local_edges(); }
  if (vrtx_status) { mesh_status &= check_local_vertices(); }
  if (mesh_status) {
    std::cout << "Local mesh items are ok!" << std::endl << std::endl;
  } else {
    std::cout << "-->>>THERE ARE TROUBLES IN LOCAL MESH ITEMS!!!" << std::endl << std::endl;
  }
  if (ornt_status) {
    mesh_status &=
        check_regn_face_orientation() &&
        check_regn_face_numbering() &&
        check_face_edge_numbering() &&
        check_vrtx_vrtx_numbering();
  }
  std::cout << "-------------------------------------------------------" << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
  if (!check_geometric_region_factor()) {
    std::cout << "-->>>THERE ARE TROUBLES IN REGION GEOMETRIC VALUES!!!" << std::endl << std::endl;
  }
  if (!check_region_face_orientation()) {
    std::cout << "-->>>THERE ARE TROUBLES IN FACE GEOMETRIC VALUES!!!" << std::endl << std::endl;
  }
}

#define check(F0, F1)              \
  bool mesh3Dv_cmaster :: check##F0##F1( size_t i ) {      \
    bool retval(true) ;              \
    std::vector<size_t> klist ;              \
    mesh.get##F0##F1( i, klist ) ;          \
    for ( size_t il=0; il<klist.size(); ++il ) {        \
      size_t k = klist[il] ;            \
      std::vector<size_t> plist ;            \
      mesh.get##F1##F0( k, plist ) ;          \
      std::sort( plist.begin(), plist.end() ) ;        \
      retval &= std::binary_search( plist.begin(), plist.end(), i ) ;  \
    }                  \
    return retval ;              \
  }

#define check_item(I, X)            \
  if ( !X ) {              \
    std::cout << std::endl << "-->" << (#I) << " = "    << I    \
   << ", failed : " << (#X) << "  !!! " << std::endl ;    \
    local_status = false ;          \
  }

// check the regions
check(_regn, _regn);

check(_regn, _face);

check(_regn, _edge);

check(_regn, _vrtx);
// check the faces
check(_face, _regn);

check(_face, _face);

check(_face, _edge);

check(_face, _vrtx);
// check the edges
check(_edge, _regn);

check(_edge, _face);

check(_edge, _edge);

check(_edge, _vrtx);
// check the vertices
check(_vrtx, _regn);

check(_vrtx, _face);

check(_vrtx, _edge);

check(_vrtx, _vrtx);

bool mesh3Dv_cmaster::check_geometric_region_factor() {
  std::cout << "-->check geom regions  " << std::endl;
  bool local_status(true);
  size_t nR = mesh.n_region();
  double meas(0.), coords[DIM];
  for (size_t s = 0; s < DIM; ++s) { coords[s] = 0.; }
  for (size_t iR = 0; iR < nR; ++iR) {
    double meas_R = mesh.get_regn_measure(iR);
    meas += meas_R;
    for (size_t s = 0; s < DIM; ++s) {
      coords[s] += meas_R * mesh.coords_R(iR, s);
    }
  }
  for (size_t s = 0; s < DIM; ++s) { coords[s] /= meas; }
  double coords_domain[DIM];
  for (size_t s = 0; s < DIM; ++s) {
    coords_domain[s] = (mesh.min_coords(s) + mesh.max_coords(s)) / 2.;
  }
  double measure_domain = 1.;
  for (size_t s = 0; s < DIM; ++s) {
    measure_domain *= (mesh.max_coords(s) - mesh.min_coords(s));
  }
  local_status = abs(meas - measure_domain) < TOL;
  for (size_t s = 0; s < DIM; ++s) {
    local_status &= abs(coords[s] - coords_domain[s]) < TOL;
  }

  double vol = meas;
  double vol_cube = measure_domain;
  double xc = coords[0];
  double yc = coords[1];
  double zc = coords[2];
  double xc_cube = coords_domain[0];
  double yc_cube = coords_domain[1];
  double zc_cube = coords_domain[2];

  if (local_status) {

    std::cout << "--> OK (cube detected)" << std::endl;

  } else {

    std::cout << "vol = " << vol << " "
              << "vol_cube = " << vol_cube << std::endl;

    std::cout << "xc = " << xc << " "
              << "yc = " << yc << " "
              << "zc = " << zc << std::endl;

    std::cout << "xc_cube = " << xc_cube << " "
              << "yc_cube = " << yc_cube << " "
              << "zc_cube = " << zc_cube << std::endl;

    std::cout << "TOL = " << TOL << std::endl;
    std::cout << "vol-vol_cube = " << vol - vol_cube << std::endl;
    std::cout << "xc-xc_cube   = " << xc - xc_cube << std::endl;
    std::cout << "yc-yc_cube   = " << yc - yc_cube << std::endl;
    std::cout << "zc-zc_cube   = " << zc - zc_cube << std::endl;

    std::cout << std::endl << ">>>If your domain is a cube<<<<" << std::endl;
    std::cout << std::endl << ">>>THE GEOMETRIC FACTORS OF MESH REGIONS ARE NOT CONSISTENT!<<<" << std::endl
              << std::endl;
  }

  std::cout << "vol = " << vol << " "
            << "vol_cube = " << vol_cube << std::endl;

  std::cout << "xc = " << xc << " "
            << "yc = " << yc << " "
            << "zc = " << zc << std::endl;

  std::cout << "xc_cube = " << xc_cube << " "
            << "yc_cube = " << yc_cube << " "
            << "zc_cube = " << zc_cube << std::endl;

  return local_status;
}

// face orientation
bool mesh3Dv_cmaster::check_region_face_orientation() {
  std::cout << "-->check geom faces    " << std::endl;
  bool local_status(true);
  size_t nR = mesh.n_region();
  //double vol(0.), xc(0.), yc(0.), zc(0.) ;
  for (size_t iR = 0; iR < nR; ++iR) {
    size_t nRF = mesh.n_regn_face(iR);
    for (size_t ilF = 0; ilF < nRF; ++ilF) {
      size_t iF = mesh.regn_face(iR, ilF);
      double vFR[3] = {mesh.coords_F(iF, 0) - mesh.coords_R(iR, 0),
                       mesh.coords_F(iF, 1) - mesh.coords_R(iR, 1),
                       mesh.coords_F(iF, 2) - mesh.coords_R(iR, 2)};
      double nmF[3] = {mesh.get_nor(iF, 0), mesh.get_nor(iF, 1), mesh.get_nor(iF, 2)};
      double vscal = vFR[0] * nmF[0] + vFR[1] * nmF[1] + vFR[2] * nmF[2];
      double sgnRF = mesh.ok_regn_face(iR, ilF) ? +1. : -1.;
      local_status &= sgnRF * vscal > 0;
      if (sgnRF * vscal < 0) {
        for (size_t i = 0; i < 20; ++i) { std::cout << "--"; }
        std::cout << std::endl;
        std::cout << "ERROR in check_region_face_orientation\n" << std::endl;
        std::cout << "iR = " << iR << " " << "ilF = " << ilF << " "
                  << "iF = " << iF << std::endl;
        std::cout << "sgnRF = " << sgnRF << " "
                  << "vscal = " << vscal << std::endl;
      }
    }
  }
  if (local_status) {
    std::cout << "--> OK" << std::endl << std::endl;
  } else {
    std::cout << std::endl << ">>>ORIENTATION OF MESH FACES IN REGIONS IS NOT CONSISTENT!<<<" << std::endl << std::endl;
  }
  return local_status;
}

bool mesh3Dv_cmaster::check_the_regions() {
  std::cout << "-->check regions  " << std::endl;
  bool local_status(true);
  size_t nR = mesh.n_region();
  for (size_t iR = 0; iR < nR; ++iR) {
    check_item(iR, check_regn_regn(iR));
    check_item(iR, check_regn_face(iR));
    check_item(iR, check_regn_edge(iR));
    check_item(iR, check_regn_vrtx(iR));
  }
  if (local_status) {
    std::cout << "--> OK" << std::endl << std::endl;
  } else {
    std::cout << std::endl << ">>>MESH REGIONS ARE NOT CONSISTENT!<<<" << std::endl << std::endl;
  }
  return local_status;
}

bool mesh3Dv_cmaster::check_the_faces() {
  std::cout << "-->check faces    " << std::endl;
  bool local_status(true);
  size_t nF = mesh.n_face();
  for (size_t iF = 0; iF < nF; ++iF) {
    check_item(iF, check_face_regn(iF));
    check_item(iF, check_face_face(iF));
    check_item(iF, check_face_edge(iF));
    check_item(iF, check_face_vrtx(iF));
  }
  if (local_status) {
    std::cout << "--> OK" << std::endl << std::endl;
  } else {
    std::cout << std::endl << ">>>MESH FACES ARE NOT CONSISTENT!<<<" << std::endl << std::endl;
  }
  return local_status;
}

bool mesh3Dv_cmaster::check_the_edges() {
  std::cout << "-->check edges    " << std::endl;
  bool local_status(true);
  size_t nE = mesh.n_edge();
  for (size_t iE = 0; iE < nE; ++iE) {
    check_item(iE, check_edge_regn(iE));
    check_item(iE, check_edge_face(iE));
    check_item(iE, check_edge_edge(iE));
    check_item(iE, check_edge_vrtx(iE));
  }
  if (local_status) {
    std::cout << "--> OK" << std::endl << std::endl;
  } else {
    std::cout << std::endl << ">>>MESH EDGES ARE NOT CONSISTENT!<<" << std::endl;
  }
  return local_status;
}

bool mesh3Dv_cmaster::check_the_vertices() {
  std::cout << "-->check vertices " << std::endl;
  bool local_status(true);
  size_t nV = mesh.n_vertex();
  for (size_t iV = 0; iV < nV; ++iV) {
    check_item(iV, check_vrtx_regn(iV));
    check_item(iV, check_vrtx_face(iV));
    check_item(iV, check_vrtx_edge(iV));
    check_item(iV, check_vrtx_vrtx(iV));
  }
  if (local_status) {
    std::cout << "--> OK" << std::endl << std::endl;
  } else {
    std::cout << std::endl << ">>>MESH VERTICES ARE NOT CONSISTENT!<<<" << std::endl << std::endl;
  }
  return local_status;
}

bool mesh3Dv_cmaster::check_local_regions() {
  std::cout << "-->check local regions  " << std::endl;
  bool local_status(true);
  std::vector<size_t> plist;
  size_t nR = mesh.n_region();

  // check region consistency: 
  // get the list of connected/adjacent items of the region iR
  // for each item get the list of the connected regions and check
  // whether iR belongs to such lists
  for (size_t iR = 0; iR < nR; ++iR) {
    // check regions: get_regn_regn, get_regn_regn
    std::vector<size_t> regn_rlist;
    mesh.get_regn_regn(iR, regn_rlist);
    size_t nRR = regn_rlist.size();
    // check regions
    for (size_t jlR = 0; jlR < nRR; ++jlR) {
      size_t jR = regn_rlist[jlR];
      plist.resize(0);
      mesh.get_regn_regn(jR, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iR)) {
        local_status = bool(false);
        std::cout << "-->>region: " << iR << " region: " << jR << " failed!!!" << std::endl;
      }
    }
    // check faces: get_regn_face, get_face_regn
    std::vector<size_t> regn_flist;
    mesh.get_regn_face(iR, regn_flist);
    size_t nRF = regn_flist.size();
    for (size_t ilF = 0; ilF < nRF; ++ilF) {
      size_t iF = regn_flist[ilF];
      plist.resize(0);
      mesh.get_face_regn(iF, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iR)) {
        local_status = bool(false);
        std::cout << "-->>region: " << iR << " face: " << iF << " failed!!!" << std::endl;
      }
    }
    // check edges: get_regn_edge, get_edge_regn
    std::vector<size_t> regn_elist;
    mesh.get_regn_edge(iR, regn_elist);
    size_t nRE = regn_elist.size();
    for (size_t ilE = 0; ilE < nRE; ++ilE) {
      size_t iE = regn_elist[ilE];
      plist.resize(0);
      mesh.get_edge_regn(iE, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iR)) {
        local_status = bool(false);
        std::cout << "-->>region: " << iR << " edge: " << iE << " failed!!!" << std::endl;
      }
    }
    // check vertices: get_regn_vrtx, get_vrtx_regn
    std::vector<size_t> regn_vlist;
    mesh.get_regn_vrtx(iR, regn_vlist);
    size_t nRV = regn_vlist.size();
    for (size_t ilV = 0; ilV < nRV; ++ilV) {
      size_t iV = regn_vlist[ilV];
      plist.resize(0);
      mesh.get_vrtx_regn(iV, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iR)) {
        local_status = bool(false);
        std::cout << "-->>region: " << iR << " vrtx: " << iV << " failed!!!" << std::endl;
      }
    }
  }
  if (local_status) {
    std::cout << "--> OK" << std::endl << std::endl;
  } else {
    std::cout << std::endl << ">>>MESH FACES ARE NOT LOCALLY CONSISTENT!<<<" << std::endl << std::endl;
  }
  return local_status;
}

bool mesh3Dv_cmaster::check_local_faces() {
  std::cout << "-->check local faces    " << std::endl;
  bool local_status(true);
  std::vector<size_t> plist;
  size_t nF = mesh.n_face();
  for (size_t iF = 0; iF < nF; ++iF) {
    // check regions
    std::vector<size_t> face_rlist;
    mesh.get_face_regn(iF, face_rlist);
    size_t nFR = face_rlist.size();
    for (size_t ilR = 0; ilR < nFR; ++ilR) {
      size_t iR = face_rlist[ilR];
      plist.resize(0);
      mesh.get_regn_face(iR, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iF)) {
        local_status = bool(false);
        std::cout << "-->>face: " << iF << " region: " << iR << " failed!!!" << std::endl;
      }
    }
    // check faces
    std::vector<size_t> face_flist;
    mesh.get_face_face(iF, face_flist);
    size_t nFF = face_flist.size();
    for (size_t jlF = 0; jlF < nFF; ++jlF) {
      size_t jF = face_flist[jlF];
      plist.resize(0);
      mesh.get_face_face(jF, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iF)) {
        local_status = bool(false);
        std::cout << "-->>face: " << iF << " face: " << jF << " failed!!!" << std::endl;
      }
    }
    // check edges
    std::vector<size_t> face_elist;
    mesh.get_face_edge(iF, face_elist);
    size_t nFE = face_elist.size();
    for (size_t ilE = 0; ilE < nFE; ++ilE) {
      size_t iE = face_elist[ilE];
      plist.resize(0);
      mesh.get_edge_face(iE, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iF)) {
        local_status = bool(false);
        std::cout << "-->>face: " << iF << " edge: " << iE << " failed!!!" << std::endl;
      }
    }
    // check vertices
    std::vector<size_t> face_vlist;
    mesh.get_face_vrtx(iF, face_vlist);
    size_t nFV = face_vlist.size();
    for (size_t ilV = 0; ilV < nFV; ++ilV) {
      size_t iV = face_vlist[ilV];
      plist.resize(0);
      mesh.get_vrtx_face(iV, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iF)) {
        local_status = bool(false);
        std::cout << "-->>face: " << iF << " vrtx: " << iV << " failed!!!" << std::endl;
      }
    }
  }
  if (local_status) {
    std::cout << "--> OK" << std::endl << std::endl;
  } else {
    std::cout << std::endl << ">>>MESH FACES ARE NOT LOCALLY CONSISTENT!<<<" << std::endl << std::endl;
  }
  return local_status;
}

bool mesh3Dv_cmaster::check_local_edges() {
  std::cout << "-->check local edges    " << std::endl;
  bool local_status(true);
  std::vector<size_t> plist;
  size_t nE = mesh.n_edge();
  for (size_t iE = 0; iE < nE; ++iE) {
    // check regions
    std::vector<size_t> edge_rlist;
    mesh.get_edge_regn(iE, edge_rlist);
    size_t nER = edge_rlist.size();
    for (size_t ilR = 0; ilR < nER; ++ilR) {
      size_t iR = edge_rlist[ilR];
      plist.resize(0);
      mesh.get_regn_edge(iR, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iE)) {
        local_status = bool(false);
        std::cout << "-->>edge: " << iE << " region: " << iR << " failed!!!" << std::endl;
      }
    }
    // check faces
    std::vector<size_t> edge_flist;
    mesh.get_edge_face(iE, edge_flist);
    size_t nEF = edge_flist.size();
    for (size_t ilF = 0; ilF < nEF; ++ilF) {
      size_t iF = edge_flist[ilF];
      plist.resize(0);
      mesh.get_face_edge(iF, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iE)) {
        local_status = bool(false);
        std::cout << "-->>edge: " << iE << " face: " << iF << " failed!!!" << std::endl;
      }
    }
    // check edges
    std::vector<size_t> edge_elist;
    mesh.get_edge_edge(iE, edge_elist);
    size_t nEE = edge_elist.size();
    for (size_t jlE = 0; jlE < nEE; ++jlE) {
      size_t jE = edge_elist[jlE];
      plist.resize(0);
      mesh.get_edge_edge(jE, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iE)) {
        local_status = bool(false);
        std::cout << "-->>edge: " << iE << " edge: " << jE << " failed!!!" << std::endl;
      }
    }
    // check vertices
    std::vector<size_t> edge_vlist;
    mesh.get_edge_vrtx(iE, edge_vlist);
    size_t nEV = edge_vlist.size();
    for (size_t ilV = 0; ilV < nEV; ++ilV) {
      size_t iV = edge_vlist[ilV];
      plist.resize(0);
      mesh.get_vrtx_edge(iV, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iE)) {
        local_status = bool(false);
        std::cout << "-->>edge: " << iE << " vrtx: " << iV << " failed!!!" << std::endl;
      }
    }
  }
  if (local_status) {
    std::cout << "--> OK" << std::endl << std::endl;
  } else {
    std::cout << std::endl << ">>>MESH EDGES ARE NOT LOCALLY CONSISTENT!<<<" << std::endl << std::endl;
  }
  return local_status;
}

bool mesh3Dv_cmaster::check_local_vertices() {
  std::cout << "-->check local vertices " << std::endl;
  bool local_status(true);
  std::vector<size_t> plist;
  size_t nV = mesh.n_vertex();
  for (size_t iV = 0; iV < nV; ++iV) {
    // check regions
    std::vector<size_t> vrtx_rlist;
    mesh.get_vrtx_regn(iV, vrtx_rlist);
    size_t nVR = vrtx_rlist.size();
    for (size_t ilR = 0; ilR < nVR; ++ilR) {
      size_t iR = vrtx_rlist[ilR];
      plist.resize(0);
      mesh.get_regn_vrtx(iR, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iV)) {
        local_status = bool(false);
        std::cout << "-->>vertex: " << iV << " region: " << iR << " failed!!!" << std::endl;
      }
    }
    // check faces
    std::vector<size_t> vrtx_flist;
    mesh.get_vrtx_face(iV, vrtx_flist);
    size_t nVF = vrtx_flist.size();
    for (size_t ilF = 0; ilF < nVF; ++ilF) {
      size_t iF = vrtx_flist[ilF];
      plist.resize(0);
      mesh.get_face_vrtx(iF, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iV)) {
        local_status = bool(false);
        std::cout << "-->>vertex: " << iV << " face: " << iF << " failed!!!" << std::endl;
      }
    }
    // check edges
    std::vector<size_t> vrtx_elist;
    mesh.get_vrtx_edge(iV, vrtx_elist);
    size_t nVE = vrtx_elist.size();
    for (size_t ilE = 0; ilE < nVE; ++ilE) {
      size_t iE = vrtx_elist[ilE];
      plist.resize(0);
      mesh.get_edge_vrtx(iE, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iV)) {
        local_status = bool(false);
        std::cout << "-->>vertex: " << iV << " edge: " << iE << " failed!!!" << std::endl;
      }
    }
    // check vertices
    std::vector<size_t> vrtx_vlist;
    mesh.get_vrtx_vrtx(iV, vrtx_vlist);
    size_t nVV = vrtx_vlist.size();
    for (size_t jlV = 0; jlV < nVV; ++jlV) {
      size_t jV = vrtx_vlist[jlV];
      plist.resize(0);
      mesh.get_vrtx_vrtx(jV, plist);
      std::sort(plist.begin(), plist.end());
      if (!std::binary_search(plist.begin(), plist.end(), iV)) {
        local_status = bool(false);
        std::cout << "-->>vertex: " << iV << " vrtx: " << jV << " failed!!!" << std::endl;
      }
    }
  }
  if (local_status) {
    std::cout << "--> OK" << std::endl << std::endl;
  } else {
    std::cout << std::endl << ">>>MESH VERTICES ARE NOT LOCALLY CONSISTENT!<<<" << std::endl << std::endl;
  }
  return local_status;
}

// loop  om all the mesh regions
// for every region, loop on the region's faces
// for every face:
// -if INTERNAL, search the same face from the region of the opposite side and
//               check that the face has (global) opposite orientation
// -if BOUNDARY, orientation must be posistive
bool mesh3Dv_cmaster::check_regn_face_orientation() {
  std::cout << "-->check regn_face orientation " << std::endl;
  bool local_status(true);
  size_t nR = mesh.n_region();
  for (size_t iR = 0; iR < nR; ++iR) {                    // loop on mesh's regions
    size_t nRF = mesh.n_regn_face(iR);
    for (size_t ilF = 0; ilF < nRF; ++ilF) {       // get region's faces
      size_t iF = mesh.regn_face(iR, ilF);
      if (mesh.is_internal_face(iF)) {      // F is an internal face:
        size_t iRp = mesh.ok_regn_face(iR, ilF) ? // check face orientation
                     mesh.face_regn(iF, 1) :              // iF points out of iR
                     mesh.face_regn(iF, 0);              // iF points into iR
        size_t nRpF = mesh.n_regn_face(iRp);
        for (size_t ilpF = 0; ilpF < nRpF; ++ilpF) {  // orientation...
          size_t ipF = mesh.regn_face(iRp, ilpF);
          if (ipF == iF) {                       // ...got it! Check consistency!
            local_status = mesh.ok_regn_face(iR, ilF) == !mesh.ok_regn_face(iRp, ilpF);
            if (!local_status) {
              std::cout << "-->>region: " << iR << " face (size_t): " << iF << " failed!!!" << std::endl;
            }
          }
        }
      } else {                                       // F is a boundary face:
        local_status = mesh.ok_regn_face(iR, ilF);   // it MUST point outward, check it!
        if (!local_status) {
          std::cout << "-->>region: " << iR << " face (bnd): " << iF << " failed!!!" << std::endl;
        }
      }
    }
  }
  if (local_status) {
    std::cout << "--> OK" << std::endl << std::endl;
  } else {
    std::cout << std::endl << ">>>FACES ARE NOT CONSISTENTLY ORIENTED IN REGIONS!<<<" << std::endl << std::endl;
  }
  return local_status;
}

// check that:
//   for every internal region the sequence of the faces is consistent with the 
//   sequence of the adjacent regions
bool mesh3Dv_cmaster::check_regn_face_numbering() {
  // the local numbering of region faces and neighbour regions should be the same
  // if the region is internal
  std::cout << "-->check regn_face_regn numbering   " << std::endl;
  bool local_status(true);
  size_t nR = mesh.n_region();
  for (size_t iR = 0; iR < nR; ++iR) {                    // loop on mesh's regions
    if (mesh.is_internal_regn(iR)) {
      std::vector<size_t> regn_rlist;
      mesh.get_regn_regn(iR, regn_rlist);
      size_t nRF = mesh.n_regn_face(iR);
      for (size_t ilF = 0; ilF < nRF; ++ilF) {              // get region's faces
        size_t iF = mesh.regn_face(iR, ilF);
        if (mesh.is_internal_face(iF)) {             // F is an internal face:
          size_t iRp = mesh.ok_regn_face(iR, ilF) ?     // check face orientation
                       mesh.face_regn(iF, 1) :                     // iF points out of iR
                       mesh.face_regn(iF, 0);                     // iF points into iR
          local_status = iRp == regn_rlist[ilF];
          if (!local_status) {
            std::cout << "-->>region: " << iR << " face (size_t): " << iF << " failed!!!" << std::endl;
          }
        } else {                                       // an internal rfegion cannot have a a boundary face!
          std::cout << "-->>region: (size_t)" << iR << " face (bnd): " << iF << " failed!!!" << std::endl;
        }
      }
    }
  }
  if (local_status) {
    std::cout << "--> OK" << std::endl << std::endl;
  } else {
    std::cout << "\n>>>REGION FACES ARE NOT CONSISTENTLY NUMBERED WITH NEIGHBOUR REGIONS!<<<" << std::endl << std::endl;
  }
  return local_status;
}

// check that:
//   for every face the sequence of the edges is consistent with the 
//   sequence of the face vertices
bool mesh3Dv_cmaster::check_face_edge_numbering() {
  std::cout << "-->check face_edge numbering   " << std::endl;
  bool local_status(true);
  size_t nF = mesh.n_face();
  for (size_t iF = 0; iF < nF; ++iF) {                    // loop on mesh's faces
    std::vector<size_t> face_vlist;
    mesh.get_face_vrtx(iF, face_vlist);
    size_t nFE = mesh.n_face_edge(iF);
    for (size_t ilE = 0; ilE < nFE; ++ilE) {        // get face's edges
      size_t iE = mesh.face_edge(iF, ilE);       // get edge ID
      size_t iV0 = mesh.edge_vrtx(iE, 0);         // get 1st edge vertex
      size_t iV1 = mesh.edge_vrtx(iE, 1);         // get 2nd edge vertex
      if (!mesh.ok_face_edge(iF, ilE)) {      // check orientation:
        iV0 = mesh.edge_vrtx(iE, 1);           // swap the vertices if the edge
        iV1 = mesh.edge_vrtx(iE, 0);           // has opposite orientation
      }
      local_status = (iV0 == face_vlist[ilE]) && (iV1 == face_vlist[(ilE + 1) % nFE]);
      if (!local_status) {
        std::cout << "-->>face: " << iF << " edge: " << nFE << " failed!!!" << std::endl;
      }
    }
  }
  if (local_status) {
    std::cout << "--> OK" << std::endl << std::endl;
  } else {
    std::cout << std::endl << ">>>FACE EDGES AND VERTICES ARE NOT CONSISTENTLY ORIENTED IN FACES!<<<" << std::endl
              << std::endl;
  }
  return local_status;
}

bool mesh3Dv_cmaster::check_vrtx_vrtx_numbering() {
  std::cout << "-->check vrtx_vrtx numbering   " << std::endl;
  bool local_status(true);
  size_t nV = mesh.n_vertex();
  for (size_t iV = 0; iV < nV; ++iV) {                    // loop on mesh's faces

    std::vector<size_t> vrtx_elist;
    mesh.get_vrtx_edge(iV, vrtx_elist);

    std::vector<size_t> vrtx_vlist;
    mesh.get_vrtx_vrtx(iV, vrtx_vlist);
    size_t nVV = vrtx_vlist.size();

    for (size_t jlV = 0; jlV < nVV; ++jlV) {
      size_t jV = vrtx_vlist[jlV];
      size_t iE = vrtx_elist[jlV];
      local_status = jV == mesh.edge_vrtx(iE, 0) || jV == mesh.edge_vrtx(iE, 1);
      // more sophisticated check using ok_vrtx_edge()
      //local_status = mesh.ok_vrtx_edge(iV,jlV) && jV==mesh.edge_vrtx(iE,1) ;
      if (!local_status) {
        std::cout << "-->>vertex: " << iV << " edge: " << iE << " failed!!!" << std::endl;
      }
    }
  }
  if (local_status) {
    std::cout << "--> OK" << std::endl << std::endl;
  } else {
    std::cout << std::endl << ">>>CONNECTED VERTICES AND EDGES ARE NOT CONSISTENTLY ORIENTED FOR VERTICES!<<<"
              << std::endl << std::endl;
  }
  return local_status;
}

} // end of namespace StemMesh3D

#endif // end of _MESH3D_MASTER_HH
