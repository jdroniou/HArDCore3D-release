#ifndef _QUAD_REGION_HH
#define _QUAD_REGION_HH

#include <cassert>

#include "tetgen.h"

namespace StemMesh3D {

class QuadRegion {

private:
  mesh_3Dv &mesh;
  QuadRuleTetra quad;

  int iR;
  tetgenio in, out;
  bool status, status_tetra;

public:
  QuadRegion(mesh_3Dv &_mesh, int _rule = 0) : mesh(_mesh), status(false), quad(_rule) {}

  ~QuadRegion() {}

  void setup(int _iR);

  void setup_tetgenio();

  void region_intg_values(double &vol_R, double &xR, double &yR, double &zR);

  double region_intg(double func(double, double, double));

  // to be used with FUNCTION_OBJECT::operator()(double,double,double)
  template<class FUNCTION_OBJECT>
  double region_integral(FUNCTION_OBJECT &func_obj);

  // to be used with FUNCTION_OBJECT::operator()(int,double,double,double)
  template<class FUNCTION_OBJECT>
  double region_integral(int i, FUNCTION_OBJECT &func_obj);

  int n_tetra();

  int get_nq();

  void set_quadrule(int iT);

  void get_quadrule(int iq, double &xq, double &yq, double &zq, double &wq);
};

int QuadRegion::n_tetra() {
  assert(status);
  return out.numberoftetrahedra;
}

int QuadRegion::get_nq() { return quad.nq(); }

// safety device: 
// when set_quadrule is called for tetrahedron iT,      status_tetra is set = to true
// when get_quadrule is called for the last quadr node, status_tetra is reset = to false
void QuadRegion::set_quadrule(int iT) {
  assert(status);
  assert(0 <= iT && iT < out.numberoftetrahedra);
  status_tetra = true;

  double xV[4], yV[4], zV[4];

  // get the tetrahedron vertices
  for (int ilV = 0; ilV < 4; ++ilV) {
    int iV = out.tetrahedronlist[iT * 4 + ilV];
    xV[ilV] = out.pointlist[3 * iV + 0];
    yV[ilV] = out.pointlist[3 * iV + 1];
    zV[ilV] = out.pointlist[3 * iV + 2];
  }

  // setup quad object
  quad.setup(xV, yV, zV);
}

void QuadRegion::get_quadrule(int iq, double &xq, double &yq, double &zq, double &wq) {
  assert(status_tetra);
  quad.get_quadrule(iq, xq, yq, zq, wq);
  status_tetra = iq != quad.nq() - 1;
}

void QuadRegion::setup(int _iR) {
  iR = _iR;
  if (status) {
    in.deinitialize();
    out.deinitialize();
  }
  setup_tetgenio();
}

void QuadRegion::setup_tetgenio() {
  status = true;
  in.initialize();
  out.initialize();

  tetgenio::facet *f;
  tetgenio::polygon *p;

  vector<int> regn_vlist;
  mesh.get_regn_vrtx(iR, regn_vlist);
  int nRV = regn_vlist.size();

  vector<int> regn_flist;
  mesh.get_regn_face(iR, regn_flist);
  int nRF = regn_flist.size();

  // INSERT POINTS
  in.firstnumber = 0;
  in.numberofpoints = nRV;
  in.pointlist = new REAL[in.numberofpoints * 3];
  // node ilV
  for (int ilV = 0; ilV < nRV; ++ilV) {
    //int iV = R.vrtx_id(ilV) ;
    int iV = regn_vlist[ilV];
    in.pointlist[3 * ilV + 0] = mesh.coords_V(iV, 0);
    in.pointlist[3 * ilV + 1] = mesh.coords_V(iV, 1);
    in.pointlist[3 * ilV + 2] = mesh.coords_V(iV, 2);
  }

  // INSERT FACETS
  in.numberoffacets = nRF;
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];

  for (int ilF = 0; ilF < nRF; ++ilF) {
    int iF = regn_flist[ilF];

    vector<int> face_vlist;
    mesh.get_face_vrtx(iF, face_vlist);
    int nFV = face_vlist.size();

    // facet ilF
    f = &in.facetlist[ilF];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = nFV;
    p->vertexlist = new int[p->numberofvertices];
    for (int ilV = 0; ilV < nFV; ++ilV) {
      int jlV = -1;
      for (int k = 0; k < nRV && jlV == -1; ++k) {
        //if ( R.vrtx_id(k)==F.vrtx_id(ilV) ) { jlV=k ;}
        if (regn_vlist[k] == face_vlist[ilV]) { jlV = k; }
      }
      assert(jlV != -1);
      p->vertexlist[ilV] = jlV; // 1
    }
  }

  // CHECK
  // output the PLC to files 'barin.node' and 'barin.poly' and then exit
  if (false) {
    char fname[] = "barin";
    in.save_nodes(fname);
    in.save_poly(fname);
    exit(0);
  }

  // TETRAHEDRALIZE
  char pQ[] = "pQ";
  tetrahedralize(pQ, &in, &out);
}

// gives volume and barycenter 
// as integrals of 1,x,y,z on the region R
// BE VERY CAREFUL: this values are integrals over R
// and contains the measure |R| !!!
void QuadRegion::region_intg_values(double &vol_R, double &xR, double &yR, double &zR) {
  assert(status);

  // set the region
  //regn_3Dv R( iR, mesh ) ;

  // retvalue and aux stuff
  vol_R = 0.;
  xR = yR = zR = 0.;
  double xV[4], yV[4], zV[4];

  // apply quadrature rule to each sub-tetrahedron
  int nT = out.numberoftetrahedra;
  for (int iT = 0; iT < nT; ++iT) {

    // get the tetrahedron vertices
    for (int ilV = 0; ilV < 4; ++ilV) {
      int iV = out.tetrahedronlist[iT * 4 + ilV];
      xV[ilV] = out.pointlist[3 * iV + 0];
      yV[ilV] = out.pointlist[3 * iV + 1];
      zV[ilV] = out.pointlist[3 * iV + 2];
    }

    // update volume
    vol_R += tetra_volume(xV, yV, zV);

    // setup quad object
    quad.setup(xV, yV, zV);

    // update barycenter's coords
    int nq = quad.nq();
    double xq, yq, zq, wq;
    for (int iq = 0; iq < nq; ++iq) {
      quad.get_quadrule(iq, xq, yq, zq, wq);
      xR += wq * xq;
      yR += wq * yq;
      zR += wq * zq;
    }
  }
}

double QuadRegion::region_intg(double func(double, double, double)) {
  assert(status);

  // set the region
  //regn_3Dv R( iR, mesh ) ;

  // retvalue and aux stuff
  double val_intg = 0.;
  double xV[4], yV[4], zV[4];

  // apply quadrature rule to each sub-tetrahedron
  int nT = out.numberoftetrahedra;
  for (int iT = 0; iT < nT; ++iT) {

    // get the tetrahedron vertices
    for (int ilV = 0; ilV < 4; ++ilV) {
      int iV = out.tetrahedronlist[iT * 4 + ilV];
      xV[ilV] = out.pointlist[3 * iV + 0];
      yV[ilV] = out.pointlist[3 * iV + 1];
      zV[ilV] = out.pointlist[3 * iV + 2];
    }

    // setup quad object on tetrahedron iT
    quad.setup(xV, yV, zV);

    // update barycenter's coords
    int nq = quad.nq();
    double xq, yq, zq, wq;
    for (int iq = 0; iq < nq; ++iq) {
      quad.get_quadrule(iq, xq, yq, zq, wq);
      val_intg += wq * func(xq, yq, zq);
    }
  }

  // return integral value
  return val_intg;
}

template<class FUNCTION_OBJECT>
double QuadRegion::region_integral(FUNCTION_OBJECT &func_obj) {
  assert(status);

  // retvalue and aux stuff
  double val_intg = 0.;
  double xV[4], yV[4], zV[4];

  // apply quadrature rule to each sub-tetrahedron
  int nT = out.numberoftetrahedra;
  for (int iT = 0; iT < nT; ++iT) {

    // get the tetrahedron vertices
    for (int ilV = 0; ilV < 4; ++ilV) {
      int iV = out.tetrahedronlist[iT * 4 + ilV];
      xV[ilV] = out.pointlist[3 * iV + 0];
      yV[ilV] = out.pointlist[3 * iV + 1];
      zV[ilV] = out.pointlist[3 * iV + 2];
    }

    // setup quad object on tetrahedron iT
    quad.setup(xV, yV, zV);

    // update barycenter's coords
    int nq = quad.nq();
    double xq, yq, zq, wq;
    for (int iq = 0; iq < nq; ++iq) {
      quad.get_quadrule(iq, xq, yq, zq, wq);
      val_intg += wq * func_obj(xq, yq, zq);
    }
  }

  // return integral value
  return val_intg;
}

template<class FUNCTION_OBJECT>
double QuadRegion::region_integral(int i, FUNCTION_OBJECT &func_obj) {
  assert(status);

  // set the region
  //regn_3Dv R( iR, mesh ) ;

  // retvalue and aux stuff
  double val_intg = 0.;
  double xV[4], yV[4], zV[4];

  // apply quadrature rule to each sub-tetrahedron
  int nT = out.numberoftetrahedra;
  for (int iT = 0; iT < nT; ++iT) {

    // get the tetrahedron vertices
    for (int ilV = 0; ilV < 4; ++ilV) {
      int iV = out.tetrahedronlist[iT * 4 + ilV];
      xV[ilV] = out.pointlist[3 * iV + 0];
      yV[ilV] = out.pointlist[3 * iV + 1];
      zV[ilV] = out.pointlist[3 * iV + 2];
    }

    // setup quad object on tetrahedron iT
    quad.setup(xV, yV, zV);

    // update barycenter's coords
    int nq = quad.nq();
    double xq, yq, zq, wq;
    for (int iq = 0; iq < nq; ++iq) {
      quad.get_quadrule(iq, xq, yq, zq, wq);
      val_intg += wq * func_obj(i, xq, yq, zq);
    }
  }

  // return integral value
  return val_intg;
}

} // end of namespace StemMesh3D

#endif // end of  _QUAD_REGION_HH
