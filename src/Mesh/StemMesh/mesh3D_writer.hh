#ifndef _MESH3D_WRITER_HH
#define _MESH3D_WRITER_HH

#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

// include class declarations
#include "mesh3D.hh"

namespace StemMesh3D {

/*

  OUTPUT DATASETS IN A SINGLE FILE USING THE MSH FORMAT (F. HUBERT)
  (USED IN THE 3-D BENCHMARK OF FVCA-VI CONFERENCE)

  // offset is set to 1 (as default value)
  //
  // The file has a very strong structure determined by the initial header lines.
  //
  // The header may be in Language::ENGLISH (default) or french (as originally was) and the
  // switching between the languages is driven by lang_flag
  //
  // The (input/ouput) is unaffected by the language choice, but must respect the 
  // file structure)
  //
  // Due to this very different structure, this class has its own implementation and
  // can be used standalone. 
  //
  // The class is also used in mesh3Dv_writer::write_mesh_MSH(...) to write a mesh in
  // MSH format using an mesh3Dv_writer object.

 */
class mesh3Dv_writer_MSH {
protected:

  enum class Language {
    FRENCH,
    ENGLISH
  };

protected:
  mesh_3Dv &mesh;
  const size_t offset;
  const Language lang_flag;

public:
  mesh3Dv_writer_MSH(mesh_3Dv &_mesh, size_t _offset = 1) :
      mesh(_mesh), offset(_offset), lang_flag(Language::ENGLISH) {}

  ~mesh3Dv_writer_MSH() {}

  // write DATASETS
  void write_header(std::ostream &OUTF, const std::string &mesh_name_str);

  void write_coords(std::ostream &OUTF);

  void write_RegnFace(std::ostream &OUTF);

  void write_RegnVrtx(std::ostream &OUTF);

  void write_FaceEdge(std::ostream &OUTF);

  void write_FaceVrtx(std::ostream &OUTF);

  void write_FaceRegn(std::ostream &OUTF);

  void write_EdgeVrtx(std::ostream &OUTF);

  // ---
  void write_mesh_MSH(const std::string &fname = "mesh3D");
};

// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
//
// general class to drive writer methods in different format
//
// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------

class mesh3Dv_writer {

  // add to MeshDefs
  enum aux_MeshDefs {
    UNDEF = -999,
    BND_CELL = -998,
    BND_FACE = -997,
    BND_NODE = -996
  };

protected:
  static const size_t DIM = 3;
  mesh_3Dv &mesh;
  size_t offset;

  void write_node_coords(const std::string &fname = "mesh3D");

  size_t n_vflag; // number of external vertex flags (=0,1)
  size_t n_fflag; // number of external face   flags (=0,1)
  size_t n_rflag; // number of external region flags (=0,1)

public:
  mesh3Dv_writer(mesh_3Dv &_mesh, size_t _offset = 0) :
      mesh(_mesh), offset(_offset), n_vflag(0), n_fflag(0), n_rflag(0) {}

  ~mesh3Dv_writer() {}

  void write_mesh_RF(const std::string &fname = "mesh3D");

  void write_mesh_Fb(const std::string &fname = "mesh3D");

  void write_mesh_MSH(const std::string &fname = "mesh3D");

  // for Sukumar output
  void write_node_coords_Suku_format(const std::string &fname);

  void write_mesh_RF_Suku_format(const std::string &fname);

  size_t get_offset() { return offset; }

  void set_offset(size_t _offset) { offset = _offset; }
};

// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
//
// methods for writing a mesh file in MSH format
//
// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------

void mesh3Dv_writer_MSH::write_header(std::ostream &OUTF, const std::string &mesh_name_str) {
  switch (lang_flag) {
    case Language::FRENCH:
      // french header
      OUTF << "Maillage cree par mesh3Dv_writer_MSH" << std::endl << std::endl << std::endl;
      OUTF << "Version" << std::endl << "1" << std::endl;
      OUTF << "Nom du maillage" << std::endl << mesh_name_str << std::endl;
      OUTF << "Infos sur le maillage" << std::endl;
      OUTF << "Nombre de sommets" << std::endl << "  " << mesh.n_vertex() << std::endl;
      OUTF << "Nombre de volumes" << std::endl << "  " << mesh.n_region() << std::endl;
      OUTF << "Nombre de faces" << std::endl << "  " << mesh.n_face() << std::endl;
      OUTF << "Nombre d\'aretes" << std::endl << "  " << mesh.n_edge() << std::endl;
      break;
    case Language::ENGLISH:
      // english header
      OUTF << "Mesh created by mesh3Dv_writer_MSH" << std::endl << std::endl << std::endl;
      OUTF << "Version" << std::endl << "1" << std::endl;
      OUTF << "Mesh name" << std::endl << mesh_name_str << std::endl;
      OUTF << "Information on the mesh" << std::endl;
      OUTF << "Number of vertices" << std::endl << "  " << mesh.n_vertex() << std::endl;
      OUTF << "Number of control volumes" << std::endl << "  " << mesh.n_region() << std::endl;
      OUTF << "Number of faces" << std::endl << "  " << mesh.n_face() << std::endl;
      OUTF << "Nomber of edges" << std::endl << "  " << mesh.n_edge() << std::endl;
      break;
  }
}

void mesh3Dv_writer_MSH::write_coords(std::ostream &OUTF) {
  // vertex coordinates
  switch (lang_flag) {
    case Language::FRENCH:
      OUTF << "Sommets " << mesh.n_vertex() << std::endl;
      break;
    case Language::ENGLISH:
      OUTF << "Vertices " << mesh.n_vertex() << std::endl;
      break;
  }
  for (size_t iV = 0; iV < mesh.n_vertex(); ++iV) {
    OUTF << std::setprecision(14)
         << std::scientific
         << "   " << mesh.coords_V(iV, 0)
         << "   " << mesh.coords_V(iV, 1)
         << "   " << mesh.coords_V(iV, 2) << std::endl;
  }
}

void mesh3Dv_writer_MSH::write_RegnFace(std::ostream &OUTF) {
  // RegnFace
  switch (lang_flag) {
    case Language::FRENCH:
      OUTF << "Volumes->faces " << mesh.n_region() << std::endl;
      break;
    case Language::ENGLISH:
      OUTF << "Volumes->faces " << mesh.n_vertex() << std::endl;
      break;
  }
  for (size_t iR = 0; iR < mesh.n_region(); ++iR) {
    std::vector<size_t> flist;
    mesh.get_regn_face(iR, flist);
    OUTF << "  " << flist.size();
    for (size_t j = 0; j < flist.size(); ++j) {
      OUTF << " " << flist[j] + offset;
    }
    OUTF << std::endl;
  }
}

void mesh3Dv_writer_MSH::write_RegnVrtx(std::ostream &OUTF) {
  // RegnVrtx
  switch (lang_flag) {
    case Language::FRENCH:
      OUTF << "Volumes->sommets " << mesh.n_region() << std::endl;
      break;
    case Language::ENGLISH:
      OUTF << "Volumes->Vertices " << mesh.n_region() << std::endl;
      break;
  }
  for (size_t iR = 0; iR < mesh.n_region(); ++iR) {
    std::vector<size_t> vlist;
    mesh.get_regn_vrtx(iR, vlist);
    OUTF << "  " << vlist.size();
    for (size_t j = 0; j < vlist.size(); ++j) {
      OUTF << " " << vlist[j] + offset;
    }
    OUTF << std::endl;
  }
}

void mesh3Dv_writer_MSH::write_FaceEdge(std::ostream &OUTF) {
  // FaceEdge
  switch (lang_flag) {
    case Language::FRENCH:
      OUTF << "Faces->Aretes " << mesh.n_face() << std::endl;
      break;
    case Language::ENGLISH:
      OUTF << "Faces->Edges " << mesh.n_face() << std::endl;
      break;
  }
  for (size_t iF = 0; iF < mesh.n_face(); ++iF) {
    std::vector<size_t> elist;
    mesh.get_face_edge(iF, elist);
    OUTF << "  " << elist.size();
    for (size_t j = 0; j < elist.size(); ++j) {
      OUTF << " " << elist[j] + offset;
    }
    OUTF << std::endl;
  }
}

void mesh3Dv_writer_MSH::write_FaceVrtx(std::ostream &OUTF) {
  // FaceVrtx
  switch (lang_flag) {
    case Language::FRENCH:
      OUTF << "Faces->Sommets " << mesh.n_face() << std::endl;
      break;
    case Language::ENGLISH:
      OUTF << "Faces->Vertices " << mesh.n_face() << std::endl;
      break;
  }
  for (size_t iF = 0; iF < mesh.n_face(); ++iF) {
    std::vector<size_t> vlist;
    mesh.get_face_vrtx(iF, vlist);
    OUTF << "  " << vlist.size();
    for (size_t j = 0; j < vlist.size(); ++j) {
      OUTF << " " << vlist[j] + offset;
    }
    OUTF << std::endl;
  }
}

void mesh3Dv_writer_MSH::write_FaceRegn(std::ostream &OUTF) {
  // FaceRegn
  switch (lang_flag) {
    case Language::FRENCH:
      OUTF << "Faces->volumes " << mesh.n_face() << std::endl;
      break;
    case Language::ENGLISH:
      OUTF << "Faces->Control volumes " << mesh.n_face() << std::endl;
      break;
  }
  for (size_t iF = 0; iF < mesh.n_face(); ++iF) {
    std::vector<size_t> rlist;
    mesh.get_face_regn(iF, rlist);
    if (rlist.size() == 1) {
      OUTF << "  " << rlist[0] + offset << "  -1" << std::endl;
    } else {
      OUTF << "  " << rlist[0] + offset
           << "  " << rlist[1] + offset << std::endl;
    }
  }
}

// ------------------------------------------------------------------------------------------
void mesh3Dv_writer_MSH::write_EdgeVrtx(std::ostream &OUTF) {
  // EdgeVrtx
  switch (lang_flag) {
    case Language::FRENCH:
      OUTF << "Aretes " << mesh.n_edge() << std::endl;
      break;
    case Language::ENGLISH:
      OUTF << "Edges " << mesh.n_edge() << std::endl;
      break;
  }
  for (size_t iE = 0; iE < mesh.n_edge(); ++iE) {
    std::vector<size_t> vlist;
    mesh.get_edge_vrtx(iE, vlist);
    assert(vlist.size() == 2);
    for (size_t j = 0; j < vlist.size(); ++j) {
      OUTF << " " << vlist[j] + offset;
    }
    OUTF << std::endl;
  }
}

// --------------------------------------------------------------------------------------------
void mesh3Dv_writer_MSH::write_mesh_MSH(const std::string &fname) {
  std::ofstream OUTF(fname + ".msh");
  write_header(OUTF, fname);
  write_coords(OUTF);
  write_RegnFace(OUTF);
  write_RegnVrtx(OUTF);
  write_FaceEdge(OUTF);
  write_FaceVrtx(OUTF);
  write_FaceRegn(OUTF);
  write_EdgeVrtx(OUTF);
  OUTF.close();
}

// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
//
// methods of class mesh3Dv_writer
//
// ----------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

void mesh3Dv_writer::write_node_coords(const std::string &fname) {
  // get number of mesh vertices
  size_t nV = mesh.n_vertex();

  // open output file for node's coordinates
  std::ofstream out_node(fname + ".node");

  // output the node file's header
  out_node << "# *.node file of 3D ATS mesh in REGN_FACE format " << std::endl;
  out_node << "# offset = " << offset << std::endl;
  out_node << "# fname  = " << fname << ".node" << std::endl;
  out_node << nV << "  3  0  " << n_vflag << std::endl;

  // output the node's coordinates
  for (size_t iV = 0; iV < nV; ++iV) {
    size_t vrtx_id = iV + offset;
    out_node << std::setw(8) << vrtx_id << "     "
             << std::setw(20) << std::setprecision(16) << mesh.coords_V(iV, 0) << "   "
             << std::setw(20) << std::setprecision(16) << mesh.coords_V(iV, 1) << "   "
             << std::setw(20) << std::setprecision(16) << mesh.coords_V(iV, 2) << "   ";
    if (n_vflag == 1) { out_node << " 1"; }
    out_node << std::endl;
  }
  out_node << "# output from mesh3Dv_writer.hh " << std::endl;
  out_node.close();
}

void mesh3Dv_writer::write_mesh_Fb(const std::string &fname) {
  std::cout << "start  mesh3Dv_writer_FACE_BASED_format::write_mesh" << std::endl;
  std::cout << "write \"" << fname << "\"" << std::endl;

  size_t nV = mesh.n_vertex();
  size_t nF = mesh.n_face();
  size_t nR = mesh.n_region();

  size_t wV = size_t(log(nV) + 2);
  //size_t wF = size_t( log(nF) + 2 ) ;
  size_t wR = size_t(log(nR) + 2);

  // write node coordinates
  write_node_coords(fname);

  // open output file for face structure
  std::ofstream out_face(fname + ".face");

  // output the faces file's header
  out_face << "# *.face file of 3D-mesh in FACE_based (Fb) format " << std::endl;
  out_face << "# offset = " << offset << std::endl;
  out_face << "# fname  = " << fname << ".face" << std::endl;
  out_face << nF << "  " << n_fflag << std::endl;

  for (size_t iF = 0; iF < nF; ++iF) {

    // get node and cell lists
    size_t face_id = iF + offset;
    std::vector<size_t> vlist, rlist;
    mesh.get_face_vrtx(iF, vlist);
    mesh.get_face_regn(iF, rlist);

    // get size (and check)
    size_t nFV = vlist.size();
    size_t nFR = rlist.size();
    assert(nFV > 3);
    assert(nFR == 1 || nFR == 2);

    // list of face's nodes
    out_face << std::setw(wV) << face_id << "   "
             << std::setw(wV) << nFV << "   ";
    for (size_t ilV = 0; ilV < nFV; ++ilV) {
      out_face << std::setw(wV) << vlist[ilV] + offset << " ";
    }

    // connected cells
    out_face << "   ";
    if (nFR == 2) {
      out_face << std::setw(wR) << rlist[0] + offset << "   ";
      out_face << std::setw(wR) << rlist[1] + offset;
    } else {
      out_face << std::setw(wR) << rlist[0] + offset << "   ";
      out_face << std::setw(wR) << BND_FACE;
    }

    // external flag
    if (n_fflag == 1) { out_face << " 1"; }

    // end of face record
    out_face << std::endl;
  }

  // close face file
  out_face.close();
  std::cout << "end of mesh3Dv_writer_FACE_BASED_format::write_mesh" << std::endl;
}

void mesh3Dv_writer::write_mesh_RF(const std::string &fname) {
  std::cout << "start  mesh3Dv_writer_REGN_FACE_format::write_mesh" << std::endl;
  std::cout << "write \"" << fname << "\"" << std::endl;

  // write node coordinates
  write_node_coords(fname);

  // set number of mesh cells
  size_t nR = mesh.n_region();

  // open output file for element's structure
  std::ofstream out_ele(fname + ".ele");

  // output the element file's header
  out_ele << "# *.ele file of 3D-mesh in REGN_FACE format " << std::endl;
  out_ele << "# offset = " << offset << std::endl;
  out_ele << "# fname  = " << fname << ".ele" << std::endl;
  out_ele << nR << "  " << n_rflag << std::endl;

  // output the element's structure
  for (size_t iR = 0; iR < mesh.n_region(); ++iR) {

    // get cell ID
    size_t cell_id = iR + offset;

    // get the list of face IDs
    std::vector<size_t> flist;
    mesh.get_regn_face(iR, flist);

    // get the list of face dirs (NOT USED)
    //std::vector<int> fdirs(flist.size()) ;
    //for ( int i=0; i<flist.size(); ++i ) { fdirs[i] = mesh.ok_regn_face(iR,i)? 1 : -1 ; }

    // set the number of the local faces
    size_t nRF = mesh.n_regn_face(iR);
    assert(nRF == flist.size());

    // write region header
    out_ele << cell_id << "  " << nRF;
    if (n_rflag == 1) { out_ele << mesh.get_fR(iR); }
    out_ele << std::endl;

    // write the face's node
    for (size_t ilF = 0; ilF < nRF; ++ilF) {

      // get the face ID
      size_t iF = flist[ilF];

      // get face's node list
      std::vector<size_t> vlist;
      mesh.get_face_vrtx(iF, vlist);

      // get the face size
      size_t nFV = vlist.size();
      out_ele << "  " << ilF + offset << "  " << nFV << "  ";

      // write the list of nodes in counterclockwise fashion
      if (mesh.ok_regn_face(iR, ilF)) {
        for (size_t ilV = 0; ilV < nFV; ++ilV) {
          size_t vrtx_ID = vlist[ilV] + offset;
          out_ele << "  " << vrtx_ID;
        }
      } else {
        for (size_t ilV = 0; ilV < nFV; ++ilV) {
          size_t vrtx_ID = vlist[nFV - 1 - ilV] + offset;
          out_ele << "  " << vrtx_ID;
        }
      }
      // end of element record
      out_ele << std::endl;
    }
  }

  out_ele << "# output from mesh3Dv_writer.hh " << std::endl;
  out_ele.close();

  std::cout << "end of mesh3Dv_writer_REGN_FACE_format::write_mesh" << std::endl;
}

void mesh3Dv_writer::write_mesh_MSH(const std::string &fname) {
  std::cout << "start  mesh3Dv_writer_MSH_format::write_mesh" << std::endl;
  std::cout << "write \"" << fname << "\"" << std::endl;

  const size_t fortran_offset = 1;
  mesh3Dv_writer_MSH mesh_writer(mesh, fortran_offset);
  mesh_writer.write_mesh_MSH(fname);

  std::cout << "end of mesh3Dv_writer_MSH_format::write_mesh" << std::endl;
}

// ---------------------------------------------------------------------

void mesh3Dv_writer::write_node_coords_Suku_format(const std::string &fname) {
  // get number of mesh vertices
  size_t nV = mesh.n_vertex();

  // open output file for node's coordinates
  std::ofstream out_node(fname + ".node");

  // output the node's coordinates
  for (size_t iV = 0; iV < nV; ++iV) {
    size_t bnd_flag = mesh.is_boundary_vrtx(iV) ? 1 : 0;
    //size_t vrtx_id = iV + offset ;
    out_node << std::setw(20) << std::setprecision(16) << mesh.coords_V(iV, 0) << "   "
             << std::setw(20) << std::setprecision(16) << mesh.coords_V(iV, 1) << "   "
             << std::setw(20) << std::setprecision(16) << mesh.coords_V(iV, 2) << "   "
             << std::setw(8) << bnd_flag;
    out_node << std::endl;
  }
  out_node.close();
}

void mesh3Dv_writer::write_mesh_RF_Suku_format(const std::string &fname) {
  std::cout << "start  mesh3Dv_writer_REGN_FACE_format::write_mesh" << std::endl;
  std::cout << "write \"" << fname << "\"" << std::endl;

  // write node coordinates
  write_node_coords_Suku_format(fname);

  // set number of mesh cells
  //size_t nR = mesh.n_region() ;

  // open output file for element's structure
  std::ofstream out_ele(fname + ".ele");

  // output the element's structure
  for (size_t iR = 0; iR < mesh.n_region(); ++iR) {

    // get cell ID
    //size_t cell_id = iR + offset ;

    // get the list of face IDs
    std::vector<size_t> flist;
    mesh.get_regn_face(iR, flist);

    // set the number of the local faces
    size_t nRF = mesh.n_regn_face(iR);
    assert(nRF == flist.size());

    // write number of faces
    out_ele << nRF << "  ";

    // write the face's node
    for (size_t ilF = 0; ilF < nRF; ++ilF) {

      // get the face ID
      size_t iF = flist[ilF];

      // get face's node list
      std::vector<size_t> vlist;
      mesh.get_face_vrtx(iF, vlist);

      // get the face size
      size_t nFV = vlist.size();
      out_ele << "  " << nFV;

      // write the list of nodes in counterclockwise fashion
      if (mesh.ok_regn_face(iR, ilF)) {
        for (size_t ilV = 0; ilV < nFV; ++ilV) {
          size_t vrtx_ID = vlist[ilV] + offset;
          out_ele << "  " << vrtx_ID;
        }
      } else {
        for (size_t ilV = 0; ilV < nFV; ++ilV) {
          size_t vrtx_ID = vlist[nFV - 1 - ilV] + offset;
          out_ele << "  " << vrtx_ID;
        }
      }
      // end of element record
      out_ele << "  ";
    }
    out_ele << std::endl;
  }
  out_ele.close();

  std::cout << "end of mesh3Dv_writer_REGN_FACE_format::write_mesh" << std::endl;
}

}  // end of namespace StemMesh3D

#endif // end of --> _MESH3D_WRITER_CC
