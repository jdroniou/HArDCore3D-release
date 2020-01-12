#ifndef _MESH3D_GMV_OUTPUT_HH
#define _MESH3D_GMV_OUTPUT_HH

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <string>

#include "mesh3D.hh"

namespace StemMesh3D {

class mesh3Dv_GMV_output {
protected:
  mesh_3Dv &mesh;

private:
  size_t offset;

protected:
  void write_vrtx_data(std::ostream &file_gmv);

  void write_regn_data(std::ostream &file_gmv);

  void write_materials(std::ostream &file_gmv);

protected:
  void write_vrtx_data(std::ostream &file_gmv, size_t iR, std::vector<size_t> &vrtx_idx);

  void write_regn_data(std::ostream &file_gmv, size_t iR, std::vector<size_t> &vrtx_idx);

  void write_materials(std::ostream &file_gmv, size_t iR);

public:
  mesh3Dv_GMV_output(mesh_3Dv &_mesh, size_t _offset = 1) :
      mesh(_mesh), offset(_offset) {}

  ~mesh3Dv_GMV_output() {}

  size_t get_offset() { return offset; }

  void set_offset(size_t _offset) { offset = _offset; }

  virtual void write_mesh(const std::string &_fname = "mesh3D");

  virtual void write_mesh(const std::string &_fname, size_t iR);
};

//-------------- ROUTINES FOR WRITING ALL CELLS
void mesh3Dv_GMV_output::write_mesh(const std::string &fname) {
  std::cout << "start  mesh3Dv_GMV_output::write_mesh" << std::endl;
  std::cout << "write \"" << fname << "\"" << std::endl;

  std::ofstream file_gmv(fname + ".gmv");

  // header of gmv file in ascii format
  file_gmv << "gmvinput ascii" << std::endl;

  // write vertex data file
  write_vrtx_data(file_gmv);

  // write region data file
  write_regn_data(file_gmv);

  // write materials (ONLY for GMV output)
  write_materials(file_gmv);

  // closure of gmv file in ascii format
  file_gmv << "endgmv" << std::endl;
  file_gmv.close();

  std::cout << "end of mesh3Dv_writer_REGN_FACE_format::write_mesh" << std::endl;
}

void mesh3Dv_GMV_output::write_vrtx_data(std::ostream &file_gmv) {
  // write vertex identifiers, coordinates & info
  file_gmv << "nodev " << mesh.n_vertex() << std::endl;
  for (size_t iV = 0; iV < mesh.n_vertex(); ++iV) {
    file_gmv << std::setprecision(16)
             << "     " << mesh.coords_V(iV, 0) << "   "
             << "     " << mesh.coords_V(iV, 1) << "   "
             << "     " << mesh.coords_V(iV, 2) << std::endl;
  }
  file_gmv << std::endl;
}

void mesh3Dv_GMV_output::write_regn_data(std::ostream &file_gmv) {
  // write vertex identifiers
  file_gmv << "cells " << mesh.n_region() << std::endl << std::endl;
  for (size_t iR = 0; iR < mesh.n_region(); ++iR) {
    size_t nRF = mesh.n_regn_face(iR);
    file_gmv << "general " << nRF << std::endl;
    for (size_t ilF = 0; ilF < nRF; ++ilF) {
      size_t iF = mesh.regn_face(iR, ilF);
      std::vector<size_t> face_vlist;
      mesh.get_face_vrtx(iF, face_vlist);
      size_t nFV = face_vlist.size();
      file_gmv << "  " << nFV;
    }
    file_gmv << std::endl << " ";
    for (size_t ilF = 0; ilF < nRF; ++ilF) {
      size_t iF = mesh.regn_face(iR, ilF);
      std::vector<size_t> face_vlist;
      mesh.get_face_vrtx(iF, face_vlist);
      size_t nFV = face_vlist.size();
      for (size_t ilV = 0; ilV < nFV; ++ilV) {
        size_t iV = face_vlist[ilV];
        file_gmv << " " << iV + get_offset();
      }
      file_gmv << "  ";
    }
    file_gmv << std::endl << std::endl;
  }
}

void mesh3Dv_GMV_output::write_materials(std::ostream &file_gmv) {
  // (no offset)
  size_t n_matr = 0;
  size_t nR = mesh.n_region();
  std::vector<size_t> regn_matr(nR);
  for (size_t iR = 0; iR < nR; ++iR) {
    mesh_3Dv::flag_type matr_R = mesh.get_fR(iR);
    if (matr_R > -1) {        // should be more correct to say != UNSET
      regn_matr[iR] = matr_R;
      ++n_matr;
    }
  }
  // write material field for cells (0)
  file_gmv << "material " << n_matr << " 0 " << std::endl;
  for (size_t i_matr = 0; i_matr < n_matr; ++i_matr) {
    file_gmv << "layer-" << (i_matr + offset) << std::endl;
  }
  size_t n = 0;
  for (size_t iR = 0; iR < nR; ++iR) {
    if (mesh.get_fR(iR) > -1) { // should be more correct to say != UNSET
      file_gmv << regn_matr[iR] + offset << " ";
      if (n % 40 == 0) { file_gmv << std::endl; }
    }
  }
}

//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

//-------------- ROUTINES FOR WRITING A SINGLE CELL
// origin is centered at cell barycenter coordinates
void mesh3Dv_GMV_output::write_mesh(const std::string &fname, size_t iR) {
  assert(0 <= iR && iR < mesh.n_region());
  std::cout << "start  mesh3Dv_GMV_output::write_mesh" << std::endl;
  std::cout << "write \"" << fname << "\"" << std::endl;

  // vrtx renumbering
  std::vector<size_t> vrtx_idx(mesh.n_vertex());
  for (size_t iV = 0; iV < mesh.n_vertex(); ++iV) { vrtx_idx[iV] = -1; }

  // open GMV file
  std::ofstream file_gmv(fname + ".gmv");

  // header of gmv file in ascii format
  file_gmv << "gmvinput ascii" << std::endl;

  // write vertex data file
  write_vrtx_data(file_gmv, iR, vrtx_idx);

  // write region data file
  write_regn_data(file_gmv, iR, vrtx_idx);

  // write materials (ONLY for GMV output)
  write_materials(file_gmv, iR);

  // closure of gmv file in ascii format
  file_gmv << "endgmv" << std::endl;
  file_gmv.close();

  std::cout << "end of mesh3Dv_writer_REGN_FACE_format::write_mesh" << std::endl;
}

void mesh3Dv_GMV_output::write_vrtx_data(std::ostream &file_gmv, size_t iR, std::vector<size_t> &vrtx_idx) {
  // build region iR
  std::vector<size_t> regn_vlist;
  mesh.get_regn_vrtx(iR, regn_vlist);
  size_t nRV = regn_vlist.size();

  // write vertex identifiers, coordinates & info

  file_gmv << "nodev " << nRV << std::endl;
  for (size_t ilV = 0; ilV < nRV; ++ilV) {
    size_t iV = regn_vlist[ilV];
    vrtx_idx[iV] = ilV;
    file_gmv << std::setprecision(16) << "  "
             << "   " << mesh.coords_V(iV, 0) - mesh.coords_R(iR, 0)
             << "   " << mesh.coords_V(iV, 1) - mesh.coords_R(iR, 1)
             << "   " << mesh.coords_V(iV, 2) - mesh.coords_R(iR, 2)
             << std::endl;
  }
  file_gmv << std::endl;
}

void mesh3Dv_GMV_output::write_regn_data(std::ostream &file_gmv, size_t iR, std::vector<size_t> &vrtx_idx) {
  // write vertex identifiers
  file_gmv << "cells 1" << std::endl << std::endl;
  size_t nRF = mesh.n_regn_face(iR);
  file_gmv << "general " << nRF << std::endl;
  for (size_t ilF = 0; ilF < nRF; ++ilF) {
    size_t iF = mesh.regn_face(iR, ilF);
    std::vector<size_t> face_vlist;
    mesh.get_face_vrtx(iF, face_vlist);
    size_t nFV = face_vlist.size();
    file_gmv << "  " << nFV;
  }
  file_gmv << std::endl << " ";
  for (size_t ilF = 0; ilF < nRF; ++ilF) {
    size_t iF = mesh.regn_face(iR, ilF);
    std::vector<size_t> face_vlist;
    mesh.get_face_vrtx(iF, face_vlist);
    size_t nFV = face_vlist.size();
    for (size_t ilV = 0; ilV < nFV; ++ilV) {
      size_t iV = face_vlist[ilV];
      file_gmv << " " << vrtx_idx[iV] + get_offset();
    }
    file_gmv << "  ";
  }
  file_gmv << std::endl << std::endl;
}

void mesh3Dv_GMV_output::write_materials(std::ostream &file_gmv, size_t iR) {
  // write material field for a single cell
  file_gmv << "material 1  0 " << std::endl;
  file_gmv << "layer-1" << std::endl << "1" << std::endl;
}

} // end of namespace StemMesh3D

#endif // end of _MESH3D_GMV_OUTPUT_CC
