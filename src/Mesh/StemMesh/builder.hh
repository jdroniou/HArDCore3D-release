#ifndef _MESH3D_BUILDER_HH
#define _MESH3D_BUILDER_HH

#ifndef USE_TETGEN
#define USE_TETGEN 0
#endif

#include <cassert>

#include <iomanip>

#include <tetra_quad.hpp>

#if USE_TETGEN
#include "quad_regn.hh"
#endif

#include "mesh3D.hh"

namespace StemMesh3D {

class mesh3Dv_builder {

  static constexpr size_t DIM = 3;
  static constexpr bool sorting_regn_face = true;

private:
  struct aux_struct {
  public:
    size_t iGlb;
    size_t iloc;
    bool bval;

    aux_struct(size_t _iGlb, size_t _iloc, bool _bval = false) :
        iGlb(_iGlb), iloc(_iloc), bval(_bval) {}

    ~aux_struct() {}

    bool operator<(const aux_struct &S1) const { return iGlb < S1.iGlb; }
  };

  typedef struct aux_struct regn_face_struct;
  typedef struct aux_struct face_edge_struct;

private:
  mesh_3Dv &mesh;

private:
  // build aux data structure
  void build_face_vlist(std::vector<std::vector<regn_face_struct> > &vec_regn_face,
                        std::vector<size_t> &face_vlist,
                        std::vector<size_t> &regn_flist);

  void build_face_elist(std::vector<std::vector<face_edge_struct> > &vec_face_edge,
                        std::vector<size_t> &edge_vlist,
                        std::vector<size_t> &face_vlist);

  // mesh coordinates
  void set_coords(const std::vector<double> &xV, const std::vector<double> &yV, const std::vector<double> &zV);

  void set_flags(const std::vector<mesh_3Dv::flag_type> &fV, const std::vector<mesh_3Dv::flag_type> &fE,
                 const std::vector<mesh_3Dv::flag_type> &fF, const std::vector<mesh_3Dv::flag_type> &fR);

  // build primary data structures
  void build_RegnFace(std::vector<std::vector<regn_face_struct> > vec_regn_face, size_t _nR);

  void build_FaceEdge(std::vector<std::vector<face_edge_struct> > vec_face_edge, size_t _nF);

  void build_EdgeVrtx(std::vector<size_t> &edge_vlist);

  // build transposed data structures
  void build_VrtxEdge();

  void build_EdgeFace();

  void build_FaceRegn();

  // build boundary lists
  void build_boundary_lists();

  // convenience function -- sort and remove duplicates from a list
  template<typename T>
  void shrink_list(std::vector<T> &tmp_list);

  // build geometric quantities
  void setup_geom_factors();

  void set_regn_geom_factors();

  void set_face_geom_factors();

  void new_region_intg_values(size_t iR, std::vector<size_t> &vlist, double &vol_R, double &xR, double &yR, double &zR);

  // check and reset face orientation (at face level)
  // (class check_master has a checking for orientation from inside regions)
  void fix_wrong_face_orientation();

public:
  mesh3Dv_builder(mesh_3Dv &_mesh) : mesh(_mesh) {}

  ~mesh3Dv_builder() {}

  void change_bbox(double new_xmin, double new_ymin, double new_zmin,
                   double new_xmax, double new_ymax, double new_zmax);

  void build_the_mesh(std::vector<double> &xV, std::vector<double> &yV, std::vector<double> &zV,
                      std::vector<mesh_3Dv::flag_type> &fV, std::vector<size_t> &regn_flist,
                      std::vector<mesh_3Dv::flag_type> &fR) {
    // -- set coordinates in mesh
    set_coords(xV, yV, zV);

    // -- core for building primary dataset
    const size_t nR = fR.size();
    std::vector<std::vector<regn_face_struct> > vec_regn_face(nR);
    std::vector<size_t> face_vlist;
    build_face_vlist(vec_regn_face, face_vlist, regn_flist);
    build_RegnFace(vec_regn_face, nR);

    const size_t nF = face_vlist.back();
    std::vector<std::vector<face_edge_struct> > vec_face_edge(nF);
    std::vector<size_t> edge_vlist;
    build_face_elist(vec_face_edge, edge_vlist, face_vlist);
    build_FaceEdge(vec_face_edge, nF);
    build_EdgeVrtx(edge_vlist);

    // -- build missing external flags, put flags into mesh
    size_t nE = edge_vlist.back();
    std::vector<mesh_3Dv::flag_type> fE(nE), fF(nF);
    for (size_t iE = 0; iE < nE; ++iE) { fE[iE] = mesh_3Dv::UNSET_FLAG; }
    for (size_t iF = 0; iF < nF; ++iF) { fF[iF] = mesh_3Dv::UNSET_FLAG; }
    set_flags(fV, fE, fF, fR);

    // -- tranpose datasets
    build_VrtxEdge();
    build_EdgeFace();
    build_FaceRegn();

    // -- build boundary lists
    build_boundary_lists();

    // -- compute/set last geometric quantities
    setup_geom_factors();

    // -- reset wrong face orientation
    fix_wrong_face_orientation();
  }
};

// setup coordinates
void mesh3Dv_builder::
set_coords(const std::vector<double> &xV, const std::vector<double> &yV, const std::vector<double> &zV) {
  assert(xV.size() == yV.size() && yV.size() == zV.size());
  size_t nV = xV.size();
  mesh.nV = nV;
  mesh.V_coords.resize(nV);
  for (size_t iV = 0; iV < nV; ++iV) {
    mesh.V_coords[iV][0] = xV[iV];
    mesh.V_coords[iV][1] = yV[iV];
    mesh.V_coords[iV][2] = zV[iV];
  }
}

// setup external flags
void mesh3Dv_builder::set_flags(const std::vector<mesh_3Dv::flag_type> &fV, const std::vector<mesh_3Dv::flag_type> &fE,
                                const std::vector<mesh_3Dv::flag_type> &fF,
                                const std::vector<mesh_3Dv::flag_type> &fR) {
  // --
  size_t nV = fV.size();
  mesh.fV.resize(nV);
  for (size_t iV = 0; iV < nV; ++iV) { mesh.fV[iV] = fV[iV]; }
  // --
  size_t nE = fE.size();
  mesh.fE.resize(nE);
  for (size_t iE = 0; iE < nE; ++iE) { mesh.fE[iE] = fE[iE]; }
  // --
  size_t nF = fF.size();
  mesh.fF.resize(nF);
  for (size_t iF = 0; iF < nF; ++iF) { mesh.fF[iF] = fF[iF]; }
  // --
  size_t nR = fR.size();
  mesh.fR.resize(nR);
  for (size_t iR = 0; iR < nR; ++iR) { mesh.fR[iR] = fR[iR]; }
  // --
}

// build up primary datasets
void mesh3Dv_builder::build_RegnFace(std::vector<std::vector<regn_face_struct>> vec_regn_face, size_t _nR) {
  size_t nR = _nR;
  mesh.nR = nR;
  mesh.RegnFace.assign(nR, std::vector<size_t>());
  mesh.RegnFaceOk.assign(nR, std::vector<bool>());
  for (size_t iR = 0; iR < nR; ++iR) {
    if (sorting_regn_face) {
      std::sort(std::begin(vec_regn_face[iR]), std::end(vec_regn_face[iR]));
    }
    size_t nRF = vec_regn_face[iR].size();
    mesh.RegnFace[iR].resize(nRF);
    mesh.RegnFaceOk[iR].resize(nRF);
    for (size_t il = 0; il < nRF; ++il) {
      size_t iF = vec_regn_face[iR][il].iGlb;
      size_t ilF = vec_regn_face[iR][il].iloc;
      bool bval = vec_regn_face[iR][il].bval;
      ilF = sorting_regn_face ? il : ilF;
      mesh.RegnFace[iR][ilF] = iF;
      mesh.RegnFaceOk[iR][ilF] = bval;
    }
  }
}

void mesh3Dv_builder::build_FaceEdge(std::vector<std::vector<face_edge_struct> > vec_face_edge, size_t _nF) {
  size_t nF = _nF;
  mesh.nF = nF;
  mesh.FaceEdge.assign(nF, std::vector<size_t>());
  mesh.FaceEdgeOk.assign(nF, std::vector<bool>());
  for (size_t iF = 0; iF < nF; ++iF) {
    size_t nFE = vec_face_edge[iF].size();
    mesh.FaceEdge[iF].resize(nFE);
    mesh.FaceEdgeOk[iF].resize(nFE);
    for (size_t il = 0; il < nFE; ++il) {
      size_t iE = vec_face_edge[iF][il].iGlb;
      size_t ilE = vec_face_edge[iF][il].iloc;
      bool bval = vec_face_edge[iF][il].bval;
      mesh.FaceEdge[iF][ilE] = iE;
      mesh.FaceEdgeOk[iF][ilE] = bval;
    }
  }
}

void mesh3Dv_builder::build_EdgeVrtx(std::vector<size_t> &edge_vlist) {
  size_t nE = edge_vlist.back();
  mesh.nE = nE;
  mesh.EdgeVrtx.resize(nE);
  size_t k = 0;
  for (size_t iE = 0; iE < nE; ++iE) {
    assert(edge_vlist[k] == 2);     // An edge always has two vertices
    k++;
    size_t kV0 = edge_vlist[k++];
    size_t kV1 = edge_vlist[k++];
    mesh.EdgeVrtx[iE][0] = kV0;
    mesh.EdgeVrtx[iE][1] = kV1;
  }
}

// build up transposed datasets
void mesh3Dv_builder::build_VrtxEdge() {
  const size_t nV = mesh.nV;
  const size_t nE = mesh.nE;
  std::vector<std::vector<aux_struct> > vec_aux(nV);
  for (size_t iE = 0; iE < nE; ++iE) {
    size_t nEV = 2;
    for (size_t ilV = 0; ilV < nEV; ++ilV) {
      size_t iV = mesh.EdgeVrtx[iE][ilV];
      assert(0 <= iV && iV < nV);
      vec_aux[iV].push_back(aux_struct(iE, ilV, bool(ilV == 0)));
    }
  }
  mesh.VrtxEdge.assign(nV, std::vector<size_t>());
  mesh.VrtxEdgeOk.assign(nV, std::vector<bool>());
  for (size_t iV = 0; iV < nV; ++iV) {
    size_t nVE = vec_aux[iV].size();
    mesh.VrtxEdge[iV].resize(nVE);
    mesh.VrtxEdgeOk[iV].resize(nVE);
    for (size_t ilE = 0; ilE < nVE; ++ilE) {
      size_t iE = vec_aux[iV][ilE].iGlb;
      bool bval = vec_aux[iV][ilE].bval;
      assert(0 <= iE && iE < nE);
      mesh.VrtxEdge[iV][ilE] = iE;
      mesh.VrtxEdgeOk[iV][ilE] = bval;
    }
  }
}

void mesh3Dv_builder::build_EdgeFace() {
  const size_t nE = mesh.nE;
  const size_t nF = mesh.nF;

  std::vector<std::vector<aux_struct> > vec_aux(nE);

  for (size_t iF = 0; iF < nF; ++iF) {
    size_t nFE = mesh.FaceEdge[iF].size();
    for (size_t ilE = 0; ilE < nFE; ++ilE) {
      size_t iE = mesh.FaceEdge[iF][ilE];
      assert(0 <= iE && iE < nE);
      bool bval = mesh.FaceEdgeOk[iF][ilE];
      vec_aux[iE].push_back(aux_struct(iF, ilE, bval));
    }
  }
  mesh.EdgeFace.assign(nE, std::vector<size_t>());
  mesh.EdgeFaceOk.assign(nE, std::vector<bool>());
  for (size_t iE = 0; iE < nE; ++iE) {
    size_t nEF = vec_aux[iE].size();
    mesh.EdgeFace[iE].resize(nEF);
    mesh.EdgeFaceOk[iE].resize(nEF);
    for (size_t ilF = 0; ilF < nEF; ++ilF) {
      size_t iF = vec_aux[iE][ilF].iGlb;
      bool bval = vec_aux[iE][ilF].bval;
      mesh.EdgeFace[iE][ilF] = iF;
      mesh.EdgeFaceOk[iE][ilF] = bval;
    }
  }
}

void mesh3Dv_builder::build_FaceRegn() {
  const size_t nF = mesh.nF;
  const size_t nR = mesh.nR;
  std::vector<std::vector<aux_struct> > vec_aux(nF);
  for (size_t iR = 0; iR < nR; ++iR) {
    size_t nRF = mesh.RegnFace[iR].size();
    for (size_t ilF = 0; ilF < nRF; ++ilF) {
      size_t iF = mesh.RegnFace[iR][ilF];
      bool bval = mesh.RegnFaceOk[iR][ilF];
      vec_aux[iF].push_back(aux_struct(iR, ilF, bval));
      assert(vec_aux[iF].size() == 1 || vec_aux[iF].size() == 2);
    }
  }
  mesh.FaceRegn.assign(nF, std::vector<size_t>());
  for (size_t iF = 0; iF < nF; ++iF) {
    size_t nFR = vec_aux[iF].size();
    assert(nFR == 1 || nFR == 2);
    mesh.FaceRegn[iF].resize(nFR);
    for (size_t ilR = 0; ilR < nFR; ++ilR) {
      size_t iR = vec_aux[iF][ilR].iGlb;
      //bool bval = vec_aux[iF][ilR].bval ;
      mesh.FaceRegn[iF][ilR] = iR;
    }
  }
}

template<typename T>
void mesh3Dv_builder::shrink_list(std::vector<T> &tmp_list) {
  std::sort(std::begin(tmp_list), std::end(tmp_list));
  tmp_list.erase(std::unique(std::begin(tmp_list), std::end(tmp_list)), std::end(tmp_list));
}

void mesh3Dv_builder::build_boundary_lists() {
  // introduce tmp vectors
  mesh.bnd_vrtx.clear();
  mesh.bnd_edge.clear();
  mesh.bnd_face.clear();
  mesh.bnd_regn.clear();

  // gather all boundary items
  for (size_t iF = 0; iF < mesh.nF; ++iF) {
    if (mesh.n_face_regn(iF) == 1) {
      mesh.bnd_face.push_back(iF);
      mesh.bnd_regn.push_back(mesh.FaceRegn[iF][0]);
      for (size_t ilE = 0; ilE < mesh.FaceEdge[iF].size(); ++ilE) {
        size_t iE = mesh.FaceEdge[iF][ilE];
        mesh.bnd_edge.push_back(iE);
        mesh.bnd_vrtx.push_back(mesh.EdgeVrtx[iE][0]);
        mesh.bnd_vrtx.push_back(mesh.EdgeVrtx[iE][1]);
      }
    }
  }

  // shrink
  shrink_list(mesh.bnd_vrtx);
  shrink_list(mesh.bnd_edge);
  shrink_list(mesh.bnd_face);
  shrink_list(mesh.bnd_regn);
}

// ------- REMARKS -------
//
// In this implementation, the face orientation coincides with the
// orientation of the face in the first of the two regions that this
// face belongs to.
//
// This region-based face orientation is given by the sequence of
// vertices that is input to the builder for the face in the region.
//
// The first region is determined by sorting all structures of type
// face_struct. So, the first region is the one between two with the
// smallest identifier.
//
// Implicitly, we assume that the face orientation of this region is
// correct, while in the other region the face may be correct
// (opposite orientation, reversed vertex sequence) or wrong (same
// orientation, same vertex sequence).
//
// The other face does not influence the regn-face based construction.
//
// NOTE that the oldest implementation of mesh3D_writer that outputs
// meshes in REGN_FACE format was bugged for this reason, as it did
// not distinguish between the orientation in the two regions.
//
// Nonetheless, everything worked properly because the first region
// was always written in the correct way, and the second (wrong) was
// never used to build mesh faces.

struct face_struct {
public:
  std::vector<size_t> vlist; // original vertex list
  std::vector<size_t> slist; // sorted   vertex list
  size_t iR; // region of the face iF
  size_t ilF; // local position of iF inside iR
  face_struct(std::vector<size_t> _vlist, size_t _iR, size_t _ilF) :
      vlist(std::move(_vlist)), slist(vlist), iR(_iR), ilF(_ilF) {
    sort(slist.begin(), slist.end());
  }

  ~face_struct() {}
};

bool operator<(const face_struct &F0, const face_struct &F1) {
  bool retval = F0.slist.size() < F1.slist.size();
  if (F0.slist.size() == F1.slist.size()) {
    size_t i = 0;
    while (i < F0.slist.size() && F0.slist[i] == F1.slist[i]) { ++i; }
    if (i < F0.slist.size()) {
      retval = F0.slist[i] < F1.slist[i];
    } else {
      retval = F0.iR < F1.iR;
    }
  }
  return retval;
}

bool operator==(const face_struct &F0, const face_struct &F1) {
  bool retval(false);
  if (F0.slist.size() == F1.slist.size()) {
    retval = true;
    for (size_t i = 0; i < F0.slist.size() && retval; ++i) {
      retval &= F0.slist[i] == F1.slist[i];
    }
  }
  return retval;
}

bool operator!=(const face_struct &F0, const face_struct &F1) {
  return !(F0 == F1);
}

// INPUT:  regn_flist
// OUTPUT:
// (i)  build face_vlist from regn_flist
// (ii) collect face indices for mesh::RegnFace
void mesh3Dv_builder::
build_face_vlist(std::vector<std::vector<regn_face_struct> > &vec_regn_face,   // output (construction)
                 std::vector<size_t> &face_vlist,      // output (new)
                 std::vector<size_t> &regn_flist) {   // input

  std::vector<size_t> vlist;
  std::vector<face_struct> vec_face_struct;

  // fill vec_face_struct with instances of face_struct type
  // from all the region faces (any order is OK)
  size_t nR = regn_flist.back();
  size_t k = 0;
  for (size_t iR = 0; iR < nR; ++iR) {         // loop on all the regions
    size_t nRF = regn_flist[k++];           // #of faces of region iR
    for (size_t ilF = 0; ilF < nRF; ++ilF) {   // loop on the faces of region iR
      size_t nFV = regn_flist[k++];         // #of vertices of face ilF
      vlist.resize(nFV);                 // resize vlist
      for (size_t ilV = 0; ilV < nFV; ++ilV) { // loop on the vertices of face ilF
        vlist[ilV] = regn_flist[k++];    // put face's vertices into vlist
      }
      vec_face_struct.emplace_back(vlist, iR, ilF);
    }
  }

  // build an ORDERED set of faces represented by face_struct instances
  // from all the region faces: internal faces are repeated twice
  std::sort(std::begin(vec_face_struct), std::end(vec_face_struct));

  { // make data structure face_vlist & set RegnFace
    // set the first face
    size_t iF = 0; // face index
    size_t kF = 0; // local pointer, be careful, kF==ilF is true only at this step!
    // --store the vertex list
    face_vlist.push_back(vec_face_struct[kF].vlist.size());
    for (size_t ilV = 0; ilV < vec_face_struct[kF].vlist.size(); ++ilV) {
      face_vlist.push_back(vec_face_struct[kF].vlist[ilV]);
    }
    { // --set RegnFace, face index is iF:=0
      size_t iR = vec_face_struct[kF].iR;
      size_t ilF = vec_face_struct[kF].ilF;
      vec_regn_face[iR].push_back(regn_face_struct(iF, ilF, bool(true)));
    }
    // set the remaining faces
    for (size_t il = 1; il < vec_face_struct.size(); ++il) {
      if (vec_face_struct[il] != vec_face_struct[kF]) {
        // found a new face
        iF++; // update the face counter
        kF = il; // reset  the face pointer
        // --store the vertex list
        face_vlist.push_back(vec_face_struct[kF].vlist.size());
        for (size_t ilV = 0; ilV < vec_face_struct[kF].vlist.size(); ++ilV) {
          face_vlist.push_back(vec_face_struct[kF].vlist[ilV]);
        }
      }
      // --set RegnFace, face index is iF, ilF is the other instance of the face iF
      // --take the orientation of the face iF in the first region
      size_t iR = vec_face_struct[il].iR;
      size_t ilF = vec_face_struct[il].ilF;
      bool bval = kF == il;
      vec_regn_face[iR].push_back(regn_face_struct(iF, ilF, bval));
    }
    // finally, set the number of mesh face
    face_vlist.push_back(iF + 1);
  }
}

struct edge_struct {
public:
  size_t V0, V1;
  size_t iF, ilF;
  bool bval;

  edge_struct(size_t _V0, size_t _V1, size_t _iF, size_t _ilF) : V0(_V0), V1(_V1), iF(_iF), ilF(_ilF), bval(true) {
    if (V0 > V1) {
      std::swap(V0, V1);
      bval = false;
    }
  }

  ~edge_struct() {}
};

bool operator<(const edge_struct &E0, const edge_struct &E1) {
  return
      (E0.V0 < E1.V0) ||
      (E0.V0 == E1.V0 && E0.V1 < E1.V1) ||
      (E0.V0 == E1.V0 && E0.V1 == E1.V1 && E0.iF < E1.iF);
}

bool operator==(const edge_struct &E0, const edge_struct &E1) {
  return E0.V0 == E1.V0 && E0.V1 == E1.V1; // may belong to different faces
}

bool operator!=(const edge_struct &E0, const edge_struct &E1) {
  return !(E0 == E1);
}

void mesh3Dv_builder::build_face_elist(std::vector<std::vector<face_edge_struct> > &vec_face_edge,
                                       std::vector<size_t> &edge_vlist,
                                       std::vector<size_t> &face_vlist) {

  std::vector<edge_struct> tmp_edge_struct;
  size_t nF = face_vlist.back();

  size_t kV = 0;
  for (size_t iF = 0; iF < nF; ++iF) {
    size_t nFV = face_vlist[kV++];
    for (size_t ilV = 0; ilV < nFV - 1; ++ilV) {
      size_t kV0 = face_vlist[kV];
      size_t kV1 = face_vlist[kV + 1];
      tmp_edge_struct.push_back(edge_struct(kV0, kV1, iF, ilV + 1)); // local pos. starts from 1 to store
      kV++;                                                     // the edge orientation in the face
    }
    size_t kV0 = face_vlist[kV];         // last  vertex of the face
    size_t kV1 = face_vlist[kV - (nFV - 1)]; // first vertex of the face: kV-(nFV-1)==0 !!!
    tmp_edge_struct.push_back(edge_struct(kV0, kV1, iF, nFV));
    kV++;
  }

  sort(tmp_edge_struct.begin(), tmp_edge_struct.end());

  { // set up datasets
    // first instance
    size_t kE = 0; // local pointer
    size_t iE = 0; // number of edges
    edge_vlist.push_back(2); // each edge has two vertices
    edge_vlist.push_back(tmp_edge_struct[kE].V0);
    edge_vlist.push_back(tmp_edge_struct[kE].V1);

    size_t iF = tmp_edge_struct[iE].iF;
    size_t ilF = tmp_edge_struct[iE].ilF;

    // be careful: ilF starts from 1 to keep sign (+/-) info!
    bool bval = tmp_edge_struct[iE].bval;
    vec_face_edge[iF].push_back(face_edge_struct(iE, ilF - 1, bval));

    for (size_t il = 1; il < tmp_edge_struct.size(); ++il) { // il = running on all the faces
      if (tmp_edge_struct[il] != tmp_edge_struct[kE]) {
        kE = il; // update pointer to new face
        iE++;    // update counter
        edge_vlist.push_back(2); // each edge has two vertices
        edge_vlist.push_back(tmp_edge_struct[il].V0);
        edge_vlist.push_back(tmp_edge_struct[il].V1);
      }
      size_t iF = tmp_edge_struct[il].iF;
      size_t ilF = tmp_edge_struct[il].ilF;
      bool bval = tmp_edge_struct[il].bval;
      vec_face_edge[iF].push_back(face_edge_struct(iE, ilF - 1, bval));
    }
    edge_vlist.push_back(++iE);
  }
}

// build geometric quantities
void mesh3Dv_builder::setup_geom_factors() {
  // setup of mesh arrays
  // regions:
  size_t nR = mesh.n_region();
  mesh.R_volume.resize(nR);
  mesh.R_diam.resize(nR);
  mesh.R_coords.resize(nR);
  // faces:
  size_t nF = mesh.n_face();
  mesh.F_area.resize(nF);
  mesh.F_diam.resize(nF);
  mesh.F_nor.resize(nF);
  mesh.F_coords.resize(nF);
  // ----------------------
  set_regn_geom_factors();
  set_face_geom_factors();
}

void mesh3Dv_builder::set_regn_geom_factors() {
  size_t nR = mesh.n_region();

#if USE_TETGEN
  int rule = 1 ;
  QuadRegion quad(mesh,rule) ;
#endif

  double xV[4], yV[4], zV[4];

  for (size_t iR = 0; iR < nR; ++iR) {
    // reset local vars
    double vol_R(0.), xR(0.), yR(0.), zR(0.);
    // set region iR
    std::vector<size_t> vlist;
    mesh.get_regn_vrtx(iR, vlist);
    size_t nRV = vlist.size();
		// Compute diameter
		mesh.R_diam[iR] = 0.0;
		for (size_t iVl=0; iVl < nRV; iVl++){
			size_t iV = vlist[iVl];
			for (size_t jVl=iVl+1; jVl < nRV; jVl++){
				size_t jV = vlist[jVl];
     		double x = mesh.coords_V(iV, 0) - mesh.coords_V(jV, 0);
        double y = mesh.coords_V(iV, 1) - mesh.coords_V(jV, 1);
        double z = mesh.coords_V(iV, 2) - mesh.coords_V(jV, 2);
				mesh.R_diam[iR] = std::max(mesh.R_diam[iR], sqrt(x*x+y*y+z*z));
			}
		}

    if (nRV == 4) { // R is a tetrahedron
      // get the tetrahedron vertices
      for (size_t ilV = 0; ilV < nRV; ++ilV) {
        size_t iV = vlist[ilV];
        xV[ilV] = mesh.coords_V(iV, 0);
        yV[ilV] = mesh.coords_V(iV, 1);
        zV[ilV] = mesh.coords_V(iV, 2);
        xR += xV[ilV];
        yR += yV[ilV];
        zR += zV[ilV];
      }
      // get the volume
      vol_R = HArDCore3D::tetra_volume(xV, yV, zV);
      xR /= double(nRV);
      yR /= double(nRV);
      zR /= double(nRV);
    } else { // R is NOT a tetrahedron
#if USE_TETGEN
      quad.setup(iR) ;
      quad.region_intg_values( vol_R, xR, yR, zR ) ;
      xR /= vol_R ;
      yR /= vol_R ;
      zR /= vol_R ;
#else // working only for convex cells with all convex faces
      // -----------------------------------------------------
      new_region_intg_values(iR, vlist, vol_R, xR, yR, zR);
      // -----------------------------------------------------
#endif
    }
    // set geometric values
    mesh.R_volume[iR] = vol_R;
    mesh.R_coords[iR][0] = xR;
    mesh.R_coords[iR][1] = yR;
    mesh.R_coords[iR][2] = zR;
  }
}

void mesh3Dv_builder::set_face_geom_factors() {
  size_t nF = mesh.n_face();
  double face_err = 0.;

  for (size_t iF = 0; iF < nF; ++iF) {

    std::vector<size_t> vlist;
    mesh.get_face_vrtx(iF, vlist);

    size_t nFV = vlist.size();

		// Compute diameter
		mesh.F_diam[iF] = 0.0;
		for (size_t iVl=0; iVl < nFV; iVl++){
			size_t iV = vlist[iVl];
			for (size_t jVl=iVl+1; jVl < nFV; jVl++){
				size_t jV = vlist[jVl];
     		double x = mesh.coords_V(iV, 0) - mesh.coords_V(jV, 0);
        double y = mesh.coords_V(iV, 1) - mesh.coords_V(jV, 1);
        double z = mesh.coords_V(iV, 2) - mesh.coords_V(jV, 2);
				mesh.F_diam[iF] = std::max(mesh.F_diam[iF], sqrt(x*x+y*y+z*z));
			}
		}

    double xc[3] = {0., 0., 0.}; // barycenter of the face
    double v0(0.), v1(0.), v2(0.), ar(0.);

    if (nFV == 3) { // F is a triangle

      double x[2], y[2], z[2];
      size_t iV0 = vlist[0];
      for (size_t il = 0; il < 2; ++il) {
        size_t iV = vlist[il + 1];
        x[il] = mesh.coords_V(iV, 0) - mesh.coords_V(iV0, 0);
        y[il] = mesh.coords_V(iV, 1) - mesh.coords_V(iV0, 1);
        z[il] = mesh.coords_V(iV, 2) - mesh.coords_V(iV0, 2);
      }
      v0 = y[0] * z[1] - y[1] * z[0];
      v1 = x[1] * z[0] - x[0] * z[1];
      v2 = x[0] * y[1] - x[1] * y[0];
      ar = sqrt(v0 * v0 + v1 * v1 + v2 * v2);

      // be careful that ar is twice the triangle area!
      for (size_t ilV = 0; ilV < nFV; ++ilV) {
        size_t iV = vlist[ilV];
        for (size_t s = 0; s < 3; ++s) {
          xc[s] += mesh.coords_V(iV, s) * ar / 3.;
        }
      }

    } else { // F is a planar polygon

      // get sub-triangle vertices
      double xF = mesh.ari_coords_F(iF, 0);
      double yF = mesh.ari_coords_F(iF, 1);
      double zF = mesh.ari_coords_F(iF, 2);

      for (size_t ilV = 0; ilV < nFV; ++ilV) {
        size_t iV = vlist[ilV];
        size_t iV1 = vlist[(ilV + 1) % nFV];

        double x[2] = {mesh.coords_V(iV, 0) - xF, mesh.coords_V(iV1, 0) - xF};
        double y[2] = {mesh.coords_V(iV, 1) - yF, mesh.coords_V(iV1, 1) - yF};
        double z[2] = {mesh.coords_V(iV, 2) - zF, mesh.coords_V(iV1, 2) - zF};

        double tmp_v0 = y[0] * z[1] - y[1] * z[0];
        double tmp_v1 = x[1] * z[0] - x[0] * z[1];
        double tmp_v2 = x[0] * y[1] - x[1] * y[0];
        double tmp_ar = sqrt(tmp_v0 * tmp_v0 + tmp_v1 * tmp_v1 + tmp_v2 * tmp_v2);

        v0 += tmp_v0;
        v1 += tmp_v1;
        v2 += tmp_v2;
        ar += tmp_ar;

        // be careful that: tmp_ar is twice the triangle area!
        // (the coefficient should be sub_area_triangle/6)
        for (size_t s = 0; s < 3; ++s) {
          xc[s] += (mesh.ari_coords_F(iF, s) +
                    mesh.coords_V(iV, s) + mesh.coords_V(iV1, s)) * tmp_ar / 3;
        }
      }
    }
    // --------------------------
    mesh.F_area[iF] = ar / 2.;
    // --------------------------
    mesh.F_nor[iF][0] = v0 / ar;
    mesh.F_nor[iF][1] = v1 / ar;
    mesh.F_nor[iF][2] = v2 / ar;
    // -----------------------------
    mesh.F_coords[iF][0] = xc[0] / ar;
    mesh.F_coords[iF][1] = xc[1] / ar;
    mesh.F_coords[iF][2] = xc[2] / ar;
    // -----------------------------
    face_err = std::max(face_err,
                        std::abs(std::sqrt(std::pow(v0 / ar, 2) + std::pow(v1 / ar, 2) + std::pow(v2 / ar, 2)) - 1.));
    assert(face_err < 1.e-14);
  } // end of --> for ( size_t iF=0; iF<nF; ++iF ) {...
}

void mesh3Dv_builder::change_bbox(double new_xmin, double new_ymin, double new_zmin,     // min vertex
                                  double new_xmax, double new_ymax, double new_zmax) {  // max vertex
  assert(new_xmax > new_xmin);
  assert(new_ymax > new_ymin);
  assert(new_zmax > new_zmin);
  assert(DIM == 3);

  double new_max[DIM] = {new_xmax, new_ymax, new_zmax};
  double new_min[DIM] = {new_xmin, new_ymin, new_zmin};

  double vol_scal = 1.;
  for (size_t s = 0; s < DIM; ++s) {
    vol_scal *= new_max[s] - new_min[s];
  }

  // set bounding box
  if (!mesh.bb_status) { mesh.eval_bbox(); }

  double bb_den[DIM];
  for (size_t s = 0; s < DIM; ++s) {
    bb_den[s] = mesh.bb_max[s] - mesh.bb_min[s];
  }

  // set new vertex coords
  for (size_t iV = 0; iV < mesh.nV; ++iV) {
    for (size_t s = 0; s < DIM; ++s) {
      mesh.V_coords[iV][s] = ((mesh.V_coords[iV][s] - mesh.bb_min[s]) * new_max[s] +
                              (mesh.bb_max[s] - mesh.V_coords[iV][s]) * new_min[s]) / bb_den[s];
    }
  }

  // set new region coords
  for (size_t iR = 0; iR < mesh.nR; ++iR) {
    for (size_t s = 0; s < DIM; ++s) {
      mesh.R_coords[iR][s] = ((mesh.R_coords[iR][s] - mesh.bb_min[s]) * new_max[s] +
                              (mesh.bb_max[s] - mesh.R_coords[iR][s]) * new_min[s]) / bb_den[s];
    }
    mesh.R_volume[iR] *= vol_scal;
  }

  // set face's geometric factors
  set_face_geom_factors();
}

void
mesh3Dv_builder::new_region_intg_values(size_t iR, std::vector<size_t> &vlist, double &vol_R, double &xR, double &yR,
                                        double &zR) {
  assert(0 <= iR && iR < mesh.n_region());

  // calculate an internal point (xR0,yR0,zR0)
  double xR0(0.), yR0(0.), zR0(0.);
  size_t nRV = vlist.size();
  for (size_t ilV = 0; ilV < nRV; ++ilV) {
    size_t iV = vlist[ilV];
    xR0 += mesh.coords_V(iV, 0);
    yR0 += mesh.coords_V(iV, 1);
    zR0 += mesh.coords_V(iV, 2);
  }
  xR0 /= double(nRV);
  yR0 /= double(nRV);
  zR0 /= double(nRV);

  // set flist
  std::vector<size_t> flist;
  mesh.get_regn_face(iR, flist);
  size_t nRF = flist.size();

  // calculate an internal point of F (xF,yF,zF)
  std::vector<double> xF(nRF), yF(nRF), zF(nRF);
  for (size_t ilF = 0; ilF < nRF; ++ilF) {
    size_t iF = flist[ilF];
    std::vector<size_t> F_vlist;
    mesh.get_face_vrtx(iF, F_vlist);
    size_t nFV = F_vlist.size();
    for (size_t ilV = 0; ilV < nFV; ++ilV) {
      size_t iV = F_vlist[ilV];
      xF[ilF] += mesh.coords_V(iV, 0);
      yF[ilF] += mesh.coords_V(iV, 1);
      zF[ilF] += mesh.coords_V(iV, 2);
    }
    xF[ilF] /= double(nFV);
    yF[ilF] /= double(nFV);
    zF[ilF] /= double(nFV);
  }

  // loop on internal diamonds and accumulate
  xR = yR = zR = 0.;
  double vR = 0.;
  double xV[4], yV[4], zV[4];
  for (size_t ilF = 0; ilF < nRF; ++ilF) {
    size_t iF = flist[ilF];
    std::vector<size_t> F_vlist;
    mesh.get_face_vrtx(iF, F_vlist);
    size_t nFV = F_vlist.size();
    for (size_t ilV = 0; ilV < nFV; ++ilV) {
      for (size_t s = 0; s < 2; ++s) {
        size_t iV = F_vlist[(ilV + s) % nFV];
        xV[s] = mesh.coords_V(iV, 0);
        yV[s] = mesh.coords_V(iV, 1);
        zV[s] = mesh.coords_V(iV, 2);
      }
      xV[2] = xF[ilF];
      yV[2] = yF[ilF];
      zV[2] = zF[ilF];
      // --------------
      xV[3] = xR0;
      yV[3] = yR0;
      zV[3] = zR0;
      // --------------
      double vol_tetra = HArDCore3D::tetra_volume(xV, yV, zV);
      vR += vol_tetra;
      xR += (xV[0] + xV[1] + xV[2] + xV[3]) / 4. * vol_tetra;
      yR += (yV[0] + yV[1] + yV[2] + yV[3]) / 4. * vol_tetra;
      zR += (zV[0] + zV[1] + zV[2] + zV[3]) / 4. * vol_tetra;
    }
  }

  // final setting
  vol_R = vR;
  xR /= vR;
  yR /= vR;
  zR /= vR;
}

// this method  is required when the input mesh does not necessarily satisfy the
// criterion that face orientation points for region "0" to region "1"
// (for example: the *.msh format provided by F. Hubert).
void mesh3Dv_builder::fix_wrong_face_orientation() {

  // this routine should be run after mesh construction is terminated
  size_t nF = mesh.n_face();
  for (size_t iF = 0; iF < nF; ++iF) {

    bool need_fixing = false;

    if (mesh.is_internal_face(iF)) {

      size_t iR0 = mesh.FaceRegn[iF][0];
      size_t iR1 = mesh.FaceRegn[iF][1];

      double nor_F[3], vFR_0[3], vFR_1[3];
      for (size_t s = 0; s < 3; ++s) {
        vFR_0[s] = mesh.coords_F(iF, s) - mesh.coords_R(iR0, s);
        vFR_1[s] = mesh.coords_F(iF, s) - mesh.coords_R(iR1, s);
        nor_F[s] = mesh.get_nor(iF, s);
      }

      double sgn_0 = vFR_0[0] * nor_F[0] + vFR_0[1] * nor_F[1] + vFR_0[2] * nor_F[2] > 0 ? 1. : -1;
      double sgn_1 = vFR_1[0] * nor_F[0] + vFR_1[1] * nor_F[1] + vFR_1[2] * nor_F[2] > 0 ? 1. : -1;

      assert(sgn_0 * sgn_1 < 0.);

      if (sgn_0 < 0. && sgn_1 > 0) {
        need_fixing = true;
      }

    } else if (mesh.is_boundary_face(iF)) {

      size_t iR0 = mesh.FaceRegn[iF][0];

      double nor_F[3], vFR_0[3];
      for (size_t s = 0; s < 3; ++s) {
        vFR_0[s] = mesh.coords_F(iF, s) - mesh.coords_R(iR0, s);
        nor_F[s] = mesh.get_nor(iF, s);
      }

      double sgn_0 = vFR_0[0] * nor_F[0] + vFR_0[1] * nor_F[1] + vFR_0[2] * nor_F[2] > 0 ? 1. : -1;

      if (sgn_0 < 0) {
        need_fixing = true;
      }
    }

    // reverse face orientation
    if (need_fixing) {

      std::vector<size_t> elist;
      mesh.get_face_edge(iF, elist);
      size_t nFE = elist.size();

      std::vector<bool> blist(nFE);
      for (size_t ilE = 0; ilE < nFE; ++ilE) {
        blist[ilE] = mesh.FaceEdgeOk[iF][ilE];
      }

      // reverse FaceEdge
      for (size_t ilE = 0; ilE < nFE; ++ilE) {
        size_t iE = elist[nFE - 1 - ilE];
        bool bE = blist[nFE - 1 - ilE];
        mesh.FaceEdge[iF][ilE] = iE;
        mesh.FaceEdge[iF][ilE] = !bE;
      }

      // reverse the value of EdgeFace
      for (size_t ilE = 0; ilE < nFE; ++ilE) {
        size_t iE = elist[nFE - 1 - ilE];
        size_t nEF = mesh.n_edge_face(iE);
        for (size_t ilF = 0; ilF < nEF; ++ilF) {
          if (mesh.EdgeFace[iE][ilF] == iF) {
            bool bval = mesh.EdgeFaceOk[iE][ilF];
            mesh.EdgeFaceOk[iE][ilF] = !bval;
          }
        }
      }

      for (size_t s = 0; s < 3; ++s) {
        mesh.F_nor[iF][s] *= -1.;
      }

    } // end of --> if ( need_fixing 
  } // end of --> for ( int iF=0; ...
}

} // end of namespace StemMesh3D

#endif
