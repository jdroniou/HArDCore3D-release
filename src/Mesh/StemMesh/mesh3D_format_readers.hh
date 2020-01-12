#ifndef _MESH3D_READERS_HH
#define _MESH3D_READERS_HH

#include <cassert>

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "builder.hh"
#include "mesh3D.hh"

namespace StemMesh3D {

//--------------------------------------------------------------------------------------------
// auxiliary functions useful for debugging
//--------------------------------------------------------------------------------------------
void print_logfile( std::vector<double> & xV, std::vector<double> & yV,  std::vector<double> & zV, std::vector<mesh_3Dv::flag_type> & fV,
                    std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR, int offset=1 ) {

  // open log file
  std::ofstream LOGF("mesh.log") ;

  // print vertex data
  size_t nV = fV.size() ;
  LOGF << "number of vertices " << nV << std::endl ;
  for ( size_t iV=0; iV<nV; ++iV ) {
    LOGF << iV+offset << "  " << xV[iV] << "  " << yV[iV] << "  " << zV[iV] << "  " << fV[iV] << std::endl ;
  }
  LOGF << "---------------------------------" << std::endl ;

  // print regn data
  size_t nR = regn_vlist.back() ;
  LOGF << "number of regions " << nR << std::endl ;
  size_t kV = 0;
  for ( size_t iR=0; iR<nR; ++iR ) {
    size_t nRV = regn_vlist[kV++] ;
    LOGF << nRV << "\t< " ;
    for ( size_t iRV=0; iRV<nRV-1; ++iRV ) {
      LOGF << regn_vlist[kV++]+offset << ", " ;
    }
    LOGF << regn_vlist[kV++]+offset << " >  " ;
    LOGF << fR[iR] << " " << std::endl ;
  }

  // close log file
  LOGF.close() ;
}

void print_regn_face( std::vector<double> & xV, std::vector<double> & yV,  std::vector<double> & zV, std::vector<mesh_3Dv::flag_type> & fV,
                      std::vector<size_t> & regn_flist, std::vector<mesh_3Dv::flag_type> & fR, int offset=1 ) {

  // open log file
  std::ofstream LOGF("mesh.log") ;

  // print vertex data
  size_t nV = fV.size() ;
  LOGF << "number of vertices " << nV << std::endl ;
  for ( size_t iV=0; iV<nV; ++iV ) {
    LOGF << iV+offset << "  " << xV[iV] << "  " << yV[iV] << "  " << zV[iV] << "  " << fV[iV] << std::endl ;
  }
  LOGF << "---------------------------------" << std::endl ;

  // print regn data
  size_t nR = regn_flist.back() ;
  LOGF << "number of regions " << nR << std::endl ;
  size_t kV = 0;
  for ( size_t iR=0; iR<nR; ++iR ) {
    size_t nRF = regn_flist[kV++] ;
    LOGF << " iR=" << iR+offset << " -->nRF=" << nRF << std::endl ;
    for ( size_t ilF=0; ilF<nRF; ++ilF ) {
      size_t nFV = regn_flist[kV++] ;
      LOGF << "  ilF=" << ilF << " -->";
      for ( size_t ilV=0; ilV<nFV; ++ilV ) {
        LOGF << "  " << regn_flist[kV++]+offset ;
      }
      LOGF << std::endl ;
    }
  }

  // close log file
  LOGF.close() ;
}

//--------------------------------------------------------------------------------------------
// base class for 3D mesh readers to be used in public derivations
//--------------------------------------------------------------------------------------------

class mesh3D_reader {

protected:
  static std::istream & eatline(std::istream & s) {
    while ( s.get() != '\n' && s.good() ) {}
    return s ;
  }

  static std::istream & eatchar(std::istream & s) { s.get() ; return s ; }

  static std::istream & eatcomments(std::istream & s) {
    char c = s.peek() ;
    while ( ( c == '!' || c == '%' || c == '#' || c == ';' || c == '$')
            && s.good() ) { s >> eatline ; c = s.peek() ; }
    return s ;
  }

  void error_message( const std::string& s ) {
    std::cerr << "fatal error:\n"
              << "mesh_reader --> read_mesh() cannot open file " << s << std::endl ;
    exit(0) ;
  }

  virtual void read_vrtx_data( std::ifstream & input_file, std::vector<double> & xV, std::vector<double> & yV, std::vector<double> & zV, std::vector<mesh_3Dv::flag_type> & fV )=0 ;
  virtual void read_regn_data( std::ifstream & input_file, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR )=0 ;

public:
  virtual void read_the_mesh  ( std::vector<double> & xV, std::vector<double> & yV,  std::vector<double> & zV, std::vector<mesh_3Dv::flag_type> & fV,
                                std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR )=0 ;
} ;

//--------------------------------------------------------------------------------------------
// read input of MESH_3D in TETGEN format
//--------------------------------------------------------------------------------------------
class mesh3D_reader_TETGEN_format : public mesh3D_reader {

private:
  std::string file_name ;
  size_t offset ; // fortran offset=1, for example tetgen

  size_t nodelem, nodface ;

  // face/region
  void read_fixed_format_int_bnd_markers( std::ifstream & inp_file, size_t nRec, size_t nItm, std::vector<size_t> & vlist, std::vector<mesh_3Dv::flag_type> & flag_list ) ;
  void read_fixed_format_ext_bnd_markers( std::ifstream & inp_file, size_t nRec, size_t nItm, std::vector<size_t> & vlist, std::vector<mesh_3Dv::flag_type> & flag_list ) ;

protected:
  virtual void read_vrtx_data( std::ifstream & input_file, std::vector<double> & xV, std::vector<double> & yV, std::vector<double> & zV, std::vector<mesh_3Dv::flag_type> & fV ) ;
  virtual void read_regn_data( std::ifstream & input_file, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) ;

public:
  mesh3D_reader_TETGEN_format( std::string _file_name, size_t _offset=1 ) :
      file_name(std::move(_file_name)), offset(_offset), nodelem(4), nodface(3) {}
  ~mesh3D_reader_TETGEN_format() {}

  virtual void read_the_mesh( std::vector<double> & xV, std::vector<double> & yV, std::vector<double> & zV, std::vector<mesh_3Dv::flag_type> & fV,
                              std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) ;
} ;
//--------------------------------------------------------------------------------------------
void mesh3D_reader_TETGEN_format :: read_the_mesh ( std::vector<double> & xV, std::vector<double> & yV, std::vector<double> & zV, std::vector<mesh_3Dv::flag_type> & fV,
                                                    std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) {

  // read .NODE file of tetgen output
  std::string node_file_name = file_name + ".node" ;
  std::ifstream node_file( node_file_name.c_str() ) ;
  if ( node_file.good() ) {
    read_vrtx_data( node_file, xV, yV, zV, fV ) ;
  } else {
    error_message(node_file_name) ;
  }
  node_file.close() ;

  // read .ELE file of tetgen output
  std::string ele_file_name = file_name + ".ele" ;
  std::ifstream ele_file( ele_file_name.c_str() ) ;
  if ( ele_file.good() ) {
    read_regn_data( ele_file, regn_vlist, fR ) ;
  } else {
    error_message(ele_file_name) ;
  }
  ele_file.close() ;
}
//--------------------------------------------------------------------------------------------
void mesh3D_reader_TETGEN_format :: read_vrtx_data ( std::ifstream & input_file,
                                                     std::vector<double> & xV, std::vector<double> & yV, std::vector<double> & zV, std::vector<mesh_3Dv::flag_type> & fV ) {
  // declarations
  double x(0.), y(0.), z(0.) ;
  size_t id(0), nV(0), ndim(0), n_attr(0), nbmrk(0);
  mesh_3Dv::flag_type f(0) ;

  // read file header
  input_file >> eatcomments >> nV >> ndim >> n_attr >> nbmrk >> eatline ;
  assert( ndim==3 ) ;

  // resize the vertex list of mesh
  // read mesh vertices
  assert(nbmrk <= 1);
  if ( nbmrk==0 ) {
    for ( size_t iV=0 ; iV<nV ; ++iV ) {
      input_file >> eatcomments >> id >> x >> y >> z >> eatline ;
      assert( id-offset==iV ) ;
      xV.push_back( x ) ;
      yV.push_back( y ) ;
      zV.push_back( z ) ;
      fV.push_back( mesh_3Dv::UNSET_FLAG ) ;
    }
  } else if ( nbmrk==1 ) {
    for ( size_t iV=0 ; iV<nV ; ++iV ) {
      input_file >> eatcomments >> id >> x >> y >> z >> f >> eatline ;
      assert( id-offset==iV ) ;
      xV.push_back( x ) ;
      yV.push_back( y ) ;
      zV.push_back( z ) ;
      fV.push_back( f ) ;
    }
  }
}
//--------------------------------------------------------------------------------------------
void mesh3D_reader_TETGEN_format ::
read_regn_data( std::ifstream & ele_file, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) {
  // declarations
  size_t nR=0, nbmrk=0 ;

  // read file's header record
  ele_file >> eatcomments >> nR >> nodelem >> nbmrk >> eatline ;

  // read the vertex region list of mesh
  if      ( nodelem>0 && nbmrk==0 ) { read_fixed_format_int_bnd_markers( ele_file, nR, nodelem, regn_vlist, fR ) ; }
  else if ( nodelem>0 && nbmrk==1 ) { read_fixed_format_ext_bnd_markers( ele_file, nR, nodelem, regn_vlist, fR ) ; }
  else { assert(0) ; }
  regn_vlist.push_back( nR ) ;
}
//--------------------------------------------------------------------------------------------
void mesh3D_reader_TETGEN_format ::
read_fixed_format_int_bnd_markers( std::ifstream & inp_file, size_t nRec, size_t nVrt, std::vector<size_t> & vlist, std::vector<mesh_3Dv::flag_type> & flag_list ) {
  size_t id(0), kV(0) ;
  for ( size_t iR=0 ; iR<nRec ; ++iR ) {
    inp_file >> eatcomments >> id ;
    vlist.push_back( nVrt ) ;
    for ( size_t k=0; k<nVrt; ++k ) {
      inp_file >> kV ;
      vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    inp_file >> eatline ;
    flag_list.push_back( mesh_3Dv::UNSET_FLAG ) ;
  }
}
//--------------------------------------------------------------------------------------------
void mesh3D_reader_TETGEN_format ::
read_fixed_format_ext_bnd_markers( std::ifstream & inp_file, size_t nRec, size_t nVrt, std::vector<size_t> & vlist, std::vector<mesh_3Dv::flag_type> & flag_list ) {
  size_t id(0), kV(0);
  mesh_3Dv::flag_type ext_flag(mesh_3Dv::UNSET_FLAG) ;
  for ( size_t iR=0 ; iR<nRec ; ++iR ) {
    inp_file >> eatcomments >> id ;
    vlist.push_back( nVrt ) ;
    for ( size_t k=0; k<nVrt; ++k ) {
      inp_file >> kV ;
      vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    inp_file >> ext_flag >> eatline ;
    flag_list.push_back( ext_flag ) ;
  }
  vlist.push_back( nRec ) ;
}
//--------------------------------------------------------------------------------------------
class mesh3D_reader_REGN_FACE_format : public mesh3D_reader {

private:
  std::string file_name ;
  size_t offset ; // fortran offset=1, for example tetgen

  size_t nodelem, nodface ;

  // face/region
  void read_fixed_format_int_bnd_markers( std::ifstream & inp_file, size_t nRec, std::vector<size_t> & flist, std::vector<mesh_3Dv::flag_type> & flag_list ) ;
  void read_fixed_format_ext_bnd_markers( std::ifstream & inp_file, size_t nRec, std::vector<size_t> & flist, std::vector<mesh_3Dv::flag_type> & flag_list ) ;
  // face/region
  void read_fixed_format_int_bnd_markers_reversed_face( std::ifstream & inp_file, size_t nRec, std::vector<size_t> & flist, std::vector<mesh_3Dv::flag_type> & flag_list ) ;
  void read_fixed_format_ext_bnd_markers_reversed_face( std::ifstream & inp_file, size_t nRec, std::vector<size_t> & flist, std::vector<mesh_3Dv::flag_type> & flag_list ) ;

protected:
  virtual void read_vrtx_data( std::ifstream & input_file, std::vector<double> & xV, std::vector<double> & yV, std::vector<double> & zV, std::vector<mesh_3Dv::flag_type> & fV ) ;
  virtual void read_regn_data( std::ifstream & input_file, std::vector<size_t> & regn_flist, std::vector<mesh_3Dv::flag_type> & fR ) ;

public:
  mesh3D_reader_REGN_FACE_format( std::string _file_name, size_t _offset=1 ) :
      file_name(std::move(_file_name)), offset(_offset), nodelem(4), nodface(3) {}
  ~mesh3D_reader_REGN_FACE_format() {}

  virtual void read_the_mesh( std::vector<double> & xV, std::vector<double> & yV, std::vector<double> & zV, std::vector<mesh_3Dv::flag_type> & fV,
                              std::vector<size_t> & regn_flist, std::vector<mesh_3Dv::flag_type> & fR ) ;

} ;
void mesh3D_reader_REGN_FACE_format :: read_the_mesh ( std::vector<double> & xV, std::vector<double> & yV, std::vector<double> & zV, std::vector<mesh_3Dv::flag_type> & fV,
                                                       std::vector<size_t> & regn_flist, std::vector<mesh_3Dv::flag_type> & fR ) {
  // read .NODE file of tetgen output
  std::string node_file_name = file_name + ".node" ;
  std::ifstream node_file( node_file_name.c_str() ) ;
  if ( node_file.good() ) {
    read_vrtx_data( node_file, xV, yV, zV, fV ) ;
  } else {
    error_message(node_file_name) ;
  }
  node_file.close() ;

  // read .ELE file of tetgen output
  std::string ele_file_name = file_name + ".ele" ;
  std::ifstream ele_file( ele_file_name.c_str() ) ;
  if ( ele_file.good() ) {
    read_regn_data( ele_file, regn_flist, fR ) ;
  } else {
    error_message(ele_file_name) ;
  }
  ele_file.close() ;
}
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
void mesh3D_reader_REGN_FACE_format :: read_vrtx_data ( std::ifstream & input_file,
                                                        std::vector<double> & xV,
                                                        std::vector<double> & yV,
                                                        std::vector<double> & zV,
                                                        std::vector<mesh_3Dv::flag_type>    & fV ) {
  // declarations
  double x(0.), y(0.), z(0.) ;
  size_t id(0), nV(0), ndim(0), n_attr(0), nbmrk(0);
  mesh_3Dv::flag_type f(0) ;

  // read file header
  input_file >> eatcomments >> nV >> ndim >> n_attr >> nbmrk >> eatline ;
  assert( ndim==3 ) ;

  // resize the vertex list of mesh
  // read mesh vertices
  assert(nbmrk <= 1);
  if ( nbmrk==0 ) {
    for ( size_t iV=0 ; iV<nV ; ++iV ) {
      input_file >> eatcomments >> id >> x >> y >> z >> eatline ;
      assert( id-offset==iV ) ;
      xV.push_back( x ) ;
      yV.push_back( y ) ;
      zV.push_back( z ) ;
      fV.push_back( mesh_3Dv::UNSET_FLAG ) ;
    }
  } else if ( nbmrk==1 ) {
    for ( size_t iV=0 ; iV<nV ; ++iV ) {
      input_file >> eatcomments >> id >> x >> y >> z >> f >> eatline ;
      assert( id-offset==iV ) ;
      xV.push_back( x ) ;
      yV.push_back( y ) ;
      zV.push_back( z ) ;
      fV.push_back( f ) ;
    }
  }
}

// nR  nbmrk
// iR  nRF  [fR]
// ilF nFV  V0 V1 V2
void mesh3D_reader_REGN_FACE_format :: read_regn_data( std::ifstream    & ele_file,
                                                       std::vector<size_t> & regn_flist,
                                                       std::vector<mesh_3Dv::flag_type> & fR ) {
  // declarations
  size_t nR=0, nbmrk=0 ;

  // read file's header record
  ele_file >> eatcomments >> nR >> nbmrk >> eatline ;

  // read the vertex region list of mesh
  if      ( nbmrk==0 ) { read_fixed_format_int_bnd_markers( ele_file, nR, regn_flist, fR ) ; }
  else if ( nbmrk==1 ) { read_fixed_format_ext_bnd_markers( ele_file, nR, regn_flist, fR ) ; }
  else { assert(0) ; }
  regn_flist.push_back( nR ) ;
}
void mesh3D_reader_REGN_FACE_format ::
read_fixed_format_int_bnd_markers( std::ifstream & inp_file, size_t nRec, std::vector<size_t> & flist, std::vector<mesh_3Dv::flag_type> & flag_list ) {
  size_t R_ig(0), F_il(0), kV(0) ;
  size_t nRF(0), nFV(0) ;
  for ( size_t iR=0 ; iR<nRec ; ++iR ) {
    inp_file >> eatcomments >> R_ig >> nRF ;
    flist.push_back( nRF ) ;
    for ( size_t ilF=0; ilF<nRF; ++ilF ) {
      inp_file >> eatcomments >> F_il >> nFV ;
      flist.push_back( nFV ) ;
      for ( size_t ilV=0; ilV<nFV; ++ilV ) {
        inp_file >> kV ;
        assert(kV >= offset);
        flist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
      }
    }
    inp_file >> eatline ;
    flag_list.push_back( mesh_3Dv::UNSET_FLAG ) ;
  }
}
void mesh3D_reader_REGN_FACE_format ::
read_fixed_format_ext_bnd_markers( std::ifstream & inp_file, size_t nRec, std::vector<size_t> & flist, std::vector<mesh_3Dv::flag_type> & flag_list ) {
  size_t R_ig(0), F_il(0), kV(0) ;
  size_t nRF(0), nFV(0);
  mesh_3Dv::flag_type ext_flag(0) ;
  for ( size_t iR=0 ; iR<nRec ; ++iR ) {
    inp_file >> eatcomments >> R_ig >> nRF >> ext_flag ;
    flist.push_back( nRF ) ;
    for ( size_t ilF=0; ilF<nRF; ++ilF ) {
      inp_file >> eatcomments >> F_il >> nFV ;
      flist.push_back( nFV ) ;
      for ( size_t ilV=0; ilV<nFV; ++ilV ) {
        inp_file >> kV ;
        assert(kV >= offset);
        flist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
      }
    }
    flag_list.push_back( ext_flag ) ;
  }
}
//--------------------------------------------------------------------------------------------
void mesh3D_reader_REGN_FACE_format ::
read_fixed_format_int_bnd_markers_reversed_face( std::ifstream & inp_file, size_t nRec, std::vector<size_t> & flist, std::vector<mesh_3Dv::flag_type> & flag_list ) {
  size_t R_ig(0), F_il(0), kV(0) ;
  size_t nRF(0), nFV(0) ;
  std::vector<size_t> tmp ;
  for ( size_t iR=0 ; iR<nRec ; ++iR ) {
    inp_file >> eatcomments >> R_ig >> nRF ;
    flist.push_back( nRF ) ;
    for ( size_t ilF=0; ilF<nRF; ++ilF ) {
      inp_file >> eatcomments >> F_il >> nFV ;
      flist.push_back( nFV ) ;
      tmp.resize(nFV) ;
      for ( size_t ilV=0; ilV<nFV; ++ilV ) {
        inp_file >> kV ;
        assert(kV >= offset);
        tmp[ilV] = kV-offset ; // offset = 1 (FORTRAN style)
      }
      for ( size_t ilV=nFV-1; ilV>=0; --ilV ) {
        flist.push_back( tmp[ilV] ) ;
      }
    }
    inp_file >> eatline ;
    flag_list.push_back( mesh_3Dv::UNSET_FLAG ) ;
  }
}
void mesh3D_reader_REGN_FACE_format ::
read_fixed_format_ext_bnd_markers_reversed_face( std::ifstream & inp_file, size_t nRec, std::vector<size_t>& flist, std::vector<mesh_3Dv::flag_type>& flag_list ) {
  size_t R_ig(0), F_il(0), kV(0) ;
  size_t nRF(0), nFV(0);
  mesh_3Dv::flag_type ext_flag(0) ;
  std::vector<size_t> tmp ;
  for ( size_t iR=0 ; iR<nRec ; ++iR ) {
    inp_file >> eatcomments >> R_ig >> nRF >> ext_flag ;
    flist.push_back( nRF ) ;
    for ( size_t ilF=0; ilF<nRF; ++ilF ) {
      inp_file >> eatcomments >> F_il >> nFV ;
      flist.push_back( nFV ) ;
      tmp.resize(nFV) ;
      for ( size_t ilV=0; ilV<nFV; ++ilV ) {
        inp_file >> kV ;
        tmp[ilV] = kV-offset ; // offset = 1 (FORTRAN style)
      }
      for ( size_t ilV=nFV-1; ilV>=0; --ilV ) {
        flist.push_back( tmp[ilV] ) ;
      }
    }
    flag_list.push_back( ext_flag ) ;
  }
}

//--------------------------------------------------------------------------------------------
// read input of MESH_3D in MSH format
// actually, only a sub-set of the records of a *.msh file are used
//--------------------------------------------------------------------------------------------
class mesh3D_reader_MSH_format : public mesh3D_reader {

private:
  std::string file_name ;
  size_t offset ; // fortran offset=1, for example tetgen

  // mesh sizes (available in file's header)
  size_t nV, nR, nF, nE ;

  // read file's header
  void read_header( std::ifstream & input_file ) ;

  template<typename T>
  void print_vector( std::vector< std::vector<T>> & vec, size_t nvec ) ;

protected:
  virtual void read_vrtx_data( std::ifstream & input_file, std::vector<double> & xV, std::vector<double> & yV, std::vector<double> & zV, std::vector<mesh_3Dv::flag_type> & fV ) ;
  virtual void read_regn_data( std::ifstream & input_file, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) ;

public:
  mesh3D_reader_MSH_format( std::string _file_name, size_t _offset=1 ) :
      file_name(std::move(_file_name)), offset(_offset) {}
  ~mesh3D_reader_MSH_format() {}

  virtual void read_the_mesh( std::vector<double> & xV, std::vector<double> & yV, std::vector<double> & zV, std::vector<mesh_3Dv::flag_type> & fV,
                              std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) ;
} ;
//--------------------------------------------------------------------------------------------
void mesh3D_reader_MSH_format ::
read_the_mesh ( std::vector<double> & xV, std::vector<double> & yV, std::vector<double> & zV, std::vector<mesh_3Dv::flag_type> & fV,
                std::vector<size_t> & regn_flist, std::vector<mesh_3Dv::flag_type> & fR ) {

  // open *.msh file
  std::string msh_file_name = file_name + ".msh" ;
  std::ifstream msh_file( msh_file_name.c_str() ) ;
  if ( msh_file.good() ) {
    read_header( msh_file ) ;
    read_vrtx_data ( msh_file, xV, yV, zV, fV ) ;
    read_regn_data ( msh_file, regn_flist, fR ) ;
  } else {
    error_message( msh_file_name ) ;
  }
  msh_file.close() ;
}
//--------------------------------------------------------------------------------------------
void  mesh3D_reader_MSH_format ::
read_header( std::ifstream & input_file ) {
  std::string str_dummy, str_gridname ;
  int n_version ;

  input_file >> str_dummy >> eatline >> eatline >> eatline ;
  input_file >> str_dummy >> eatline >> n_version    >> eatline ;
  input_file >> str_dummy >> eatline >> str_gridname >> eatline ;
  input_file >> str_dummy >> eatline ;
  input_file >> str_dummy >> eatline >> nV >> eatline ;
  input_file >> str_dummy >> eatline >> nR >> eatline ;
  input_file >> str_dummy >> eatline >> nF >> eatline ;
  input_file >> str_dummy >> eatline >> nE >> eatline ;
}
//--------------------------------------------------------------------------------------------
void mesh3D_reader_MSH_format ::
read_vrtx_data ( std::ifstream & input_file,
                 std::vector<double> & xV, std::vector<double> & yV, std::vector<double> & zV, std::vector<mesh_3Dv::flag_type> & fV ) {
  // declarations
  double x(0.), y(0.), z(0.) ;

  // read mesh vertices
  input_file >> eatline ;
  for ( size_t iV=0 ; iV<nV ; ++iV ) {
    input_file >> x >> y >> z >> eatline ;
    xV.push_back( x ) ;
    yV.push_back( y ) ;
    zV.push_back( z ) ;
    fV.push_back( mesh_3Dv::UNSET_FLAG ) ;
  }
}
//--------------------------------------------------------------------------------------------
// re-add offset to compare with input data
template<typename T>
void mesh3D_reader_MSH_format::print_vector( std::vector< std::vector<T> > & vec, size_t nvec ) {
  for ( size_t i=0; i<nvec; ++i ) {
    size_t n = vec[i].size() ;
    std::cout << i+offset << "-> " << n << " <- " ;
    for ( size_t il=0; il<n; ++il ) {
      std::cout << vec[i][il]+offset << "  " ;
    }
    std::cout << std::endl ;
  }
}
//--------------------------------------------------------------------------------------------
void mesh3D_reader_MSH_format ::
read_regn_data( std::ifstream & input_file, std::vector<size_t> & regn_flist, std::vector<mesh_3Dv::flag_type> & fR ) {

  std::string str_dummy("") ;

  // aux STL vectors
  std::vector< std::vector<size_t> > R_flist(nR) ; // *
  //std::vector<size_t> R_vlist[nR] ;
  std::vector< std::vector<size_t> > F_vlist(nF) ; // *
  //std::vector<size_t> F_elist[nF] ;
  std::vector< std::vector<size_t> > F_rlist(nF) ;
  //std::vector<size_t> E_vlist[nE] ;

  bool log_msh = false ;
  std::ofstream LOG("./msh.log") ;

  { // Volumes->faces
    input_file >> str_dummy >> eatline ;
    if ( log_msh ) { LOG << str_dummy << std::endl ; }
    for ( size_t iR=0 ; iR<nR ; ++iR ) {
      size_t nRF;
      input_file >> nRF ;
      if ( log_msh ) { LOG << nRF << "-->" ; }
      for ( size_t ilF=0; ilF<nRF; ++ilF ) {
        size_t iF;
        input_file >> iF ;
        if ( log_msh ) { LOG << "  " << iF ; }
        assert(iF >= offset);
        R_flist[iR].push_back( iF-offset) ; // offset = 1 (FORTRAN style)
      }
      input_file >> eatline ;
      if ( log_msh ) { LOG << std::endl ; }
    }
  }

  { // Volumes->sommets
    input_file >> str_dummy >> eatline ;
    if ( log_msh ) { LOG << str_dummy << std::endl ; }
    for ( size_t iR=0 ; iR<nR ; ++iR ) {
      size_t nRV ;
      input_file >> nRV ;
      for ( size_t ilV=0; ilV<nRV; ++ilV ) {
        size_t iV;
        input_file >> iV ;
        if ( log_msh ) { LOG << "  " << iV ; }
        assert(iV >= offset);
        //R_vlist[iR].push_back( iV-offset) ; // offset = 1 (FORTRAN style)
      }
      input_file >> eatline ;
      if ( log_msh ) { LOG << std::endl ; }
    }
  }

  { // Faces->aretes
    input_file >> str_dummy >> eatline ;
    if ( log_msh ) { LOG << str_dummy << std::endl ; }
    for ( size_t iF=0; iF<nF; ++iF ) {
      size_t nFE;
      input_file >> nFE ;
      for ( size_t ilE=0; ilE<nFE; ++ilE ) {
        size_t iE;
        input_file >> iE ;
        if ( log_msh ) { LOG << "  " << iE ; }
        assert(iE >= offset);
        //F_elist[iF].push_back( iE-offset) ; // offset = 1 (FORTRAN style)
      }
      input_file >> eatline ;
      if ( log_msh ) { LOG << std::endl ; }
    }
  }

  { // Faces->sommets
    input_file >> str_dummy >> eatline ;
    if ( log_msh ) { LOG << str_dummy << std::endl ; }
    for ( size_t iF=0 ; iF<nF ; ++iF ) {
      size_t nFV;
      input_file >> nFV ;
      for ( size_t ilV=0; ilV<nFV; ++ilV ) {
        size_t iV;
        input_file >> iV ;
        if ( log_msh ) { LOG << "  " << iV ; }
        assert(iV >= offset);
        F_vlist[iF].push_back( iV-offset) ; // offset = 1 (FORTRAN style)
      }
      input_file >> eatline ;
      if ( log_msh ) { LOG << std::endl ; }
    }
  }

  { // Faces->volumes
    input_file >> str_dummy >> eatline ;
    if ( log_msh ) { LOG << str_dummy << std::endl ; }
    for ( size_t iF=0 ; iF<nF ; ++iF ) {
      for ( size_t ilR=0; ilR<2; ++ilR ) {
        size_t iR;
        input_file >> iR ;
        if ( log_msh ) { LOG << "  " << iR ; }
        assert(iR >= offset);
        F_rlist[iF].push_back( iR-offset ) ; // offset = 1 (FORTRAN style)
      }
      input_file >> eatline ;
      if ( log_msh ) { LOG << std::endl ; }
    }
  }

  { // Aretes->sommets
    input_file >> str_dummy >> eatline ;
    if ( log_msh ) { LOG << str_dummy << std::endl ; }
    for ( size_t iE=0 ; iE<nE ; ++iE ) {
      for ( size_t ilV=0; ilV<2; ++ilV ) {
        size_t iV;
        input_file >> iV ;
        if ( log_msh ) { LOG << "  " << iV ; }
        assert(iV >= offset);
        //E_vlist[iE].push_back( iV-offset) ; // offset = 1 (FORTRAN style)
      }
      input_file >> eatline ;
      if ( log_msh ) { LOG << std::endl ; }
    }
  }

  // reconstruct regn-face structure
  for ( size_t iR=0; iR<nR; ++iR ) {
    size_t nRF = R_flist[iR].size() ;
    regn_flist.push_back( nRF ) ;
    for ( size_t ilF=0; ilF<nRF; ++ilF ) {
      size_t iF = R_flist[iR][ilF] ;
      size_t  nFV = F_vlist[iF].size() ;
      regn_flist.push_back( nFV ) ;
      bool ok_regn_face = (F_rlist[iF][0] == iR);
      if ( ok_regn_face ) {
        for ( size_t ilV=0; ilV<nFV; ++ilV ) {
          size_t iV = F_vlist[iF][ilV] ;
          regn_flist.push_back( iV ) ;
        }
      } else {
        for ( size_t ilV=0; ilV<nFV; ++ilV ) {
          size_t iV = F_vlist[iF][nFV-1-ilV] ;
          regn_flist.push_back( iV ) ;
        }
      }
    }
    fR.push_back( mesh_3Dv::UNSET_FLAG ) ;
  }
  regn_flist.push_back(nR) ;
}

//--------------------------------------------------------------------------------------------
// base class for mesh readers for public derivations
//--------------------------------------------------------------------------------------------
class mesh2D_reader {

protected:
  static std::istream & eatline(std::istream & s) {
    while ( s.get() != '\n' && s.good() ) {}
    return s ;
  }

  static std::istream & eatchar(std::istream & s) { s.get() ; return s ; }

  static std::istream & eatcomments(std::istream & s) {
    char c = s.peek() ;
    while ( ( c == '!' || c == '%' || c == '#' || c == ';' || c == '$')
            && s.good() ) { s >> eatline ; c = s.peek() ; }
    return s ;
  }

  void error_message( const std::string & s ) {
    std::cerr << "fatal error:\n"
              << "mesh_reader --> read_mesh() cannot open file " << s << std::endl ;
    exit(0) ;
  }

  virtual void read_vrtx_data( std::ifstream & input_file, std::vector<double> & xV, std::vector<double> & yV, std::vector<mesh_3Dv::flag_type> & fV )=0 ;
  virtual void read_regn_data( std::ifstream & input_file, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR )=0 ;

public:
  virtual void read_mesh( std::vector<double> & xV, std::vector<double> & yV, std::vector<mesh_3Dv::flag_type> & fV, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR )=0 ;
} ;

//--------------------------------------------------------------------------------------------
// read input of MESH_2D in "General Format", which is a generalization of TRIANGLE format
//--------------------------------------------------------------------------------------------

class mesh2D_reader_GeneralFormat : public mesh2D_reader {

private:
  std::string file_name ;
  size_t offset ; // fortran offset

  size_t nodelem ;

  void read_regn_free_format_int_bnd_markers ( std::ifstream & ele_file, size_t nR, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) ;
  void read_regn_free_format_ext_bnd_markers ( std::ifstream & ele_file, size_t nR, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) ;
  void read_regn_fixed_format_int_bnd_markers( std::ifstream & ele_file, size_t nR, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) ;
  void read_regn_fixed_format_ext_bnd_markers( std::ifstream & ele_file, size_t nR, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) ;

protected:
  virtual void read_vrtx_data( std::ifstream & input_file, std::vector<double> & xV, std::vector<double> & yV, std::vector<mesh_3Dv::flag_type> & fV ) ;
  virtual void read_regn_data( std::ifstream & input_file, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) ;

public:
  mesh2D_reader_GeneralFormat( std::string _file_name, size_t _offset=1 ) : file_name(std::move(_file_name)), offset(_offset) {}
  ~mesh2D_reader_GeneralFormat() {}

  virtual void read_mesh( std::vector<double> & xV, std::vector<double> & yV, std::vector<mesh_3Dv::flag_type> & fV, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) ;
} ;
void mesh2D_reader_GeneralFormat :: read_mesh ( std::vector<double> & xV, std::vector<double> & yV, std::vector<mesh_3Dv::flag_type> & fV,
                                                std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) {
  // read NODE file of triangle output
  std::string node_file_name = file_name + ".node" ;
  std::ifstream node_file( node_file_name.c_str() ) ;
  if ( node_file.good() ) {
    read_vrtx_data( node_file, xV, yV, fV ) ;
  } else {
    error_message(node_file_name) ;
  }
  node_file.close() ;
  // read ELE file of triangle output
  std::string ele_file_name = file_name + ".ele" ;
  std::ifstream ele_file( ele_file_name.c_str() ) ;
  if ( ele_file.good() ) {
    read_regn_data( ele_file, regn_vlist, fR ) ;
  } else {
    error_message(ele_file_name) ;
  }
  ele_file.close() ;
}
void mesh2D_reader_GeneralFormat :: read_vrtx_data ( std::ifstream & input_file,
                                                     std::vector<double> & xV,
                                                     std::vector<double> & yV,
                                                     std::vector<mesh_3Dv::flag_type>    & fV ) {
  // declarations
  double x(0.), y(0.) ;
  size_t id(0), nV(0), ndim(0), n_attr(0), nbmrk(0);
  mesh_3Dv::flag_type f(0) ;

  // read file header
  input_file >> eatcomments >> nV >> ndim >> n_attr >> nbmrk >> eatline ;

  // resize the vertex list of mesh
  // read mesh vertices
  assert(nbmrk <= 1);
  if ( nbmrk==0 ) {
    for ( size_t iV=0 ; iV<nV ; ++iV ) {
      input_file >> eatcomments >> id >> x >> y >> eatline ;
      assert( id-offset==iV ) ;
      xV.push_back( x ) ;
      yV.push_back( y ) ;
      fV.push_back( mesh_3Dv::UNSET_FLAG ) ;
    }
  } else if ( nbmrk==1 ) {
    for ( size_t iV=0 ; iV<nV ; ++iV ) {
      input_file >> eatcomments >> id >> x >> y >> f >> eatline ;
      assert( id-offset==iV ) ;
      xV.push_back( x ) ;
      yV.push_back( y ) ;
      fV.push_back( f ) ;
    }
  }
}
void mesh2D_reader_GeneralFormat ::
read_regn_free_format_int_bnd_markers( std::ifstream & ele_file, size_t nR, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) {
  size_t id(0), nRV(0), kV(0) ;
  for ( size_t iR=0 ; iR<nR ; ++iR ) {
    ele_file >> eatcomments >> id >> nRV ;
    regn_vlist.push_back( nRV ) ;
    for ( size_t k=0; k<nRV; ++k ) {
      ele_file >> kV ;
      assert(kV >= offset);
      regn_vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    ele_file >> eatline ;
    fR.push_back( mesh_3Dv::UNSET_FLAG ) ;
  }
}
void mesh2D_reader_GeneralFormat ::
read_regn_free_format_ext_bnd_markers( std::ifstream & ele_file, size_t nR, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) {
  size_t id(0), nRV(0), kV(0);
  mesh_3Dv::flag_type regn_flag(mesh_3Dv::UNSET_FLAG) ;
  for ( size_t iR=0 ; iR<nR ; ++iR ) {
    ele_file >> eatcomments >> id >> nRV ;
    regn_vlist.push_back( nRV ) ;
    for ( size_t k=0; k<nRV; ++k ) {
      ele_file >> kV ;
      regn_vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    ele_file >> regn_flag >> eatline ;
    fR.push_back( regn_flag ) ;
  }
  regn_vlist.push_back( nR ) ;
}
void mesh2D_reader_GeneralFormat ::
read_regn_fixed_format_int_bnd_markers( std::ifstream & ele_file, size_t nR, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) {
  size_t id(0), kV(0) ;
  for ( size_t iR=0 ; iR<nR ; ++iR ) {
    ele_file >> eatcomments >> id ;
    regn_vlist.push_back( nodelem ) ;
    for ( size_t k=0; k<nodelem; ++k ) {
      ele_file >> kV ;
      assert(kV >= offset);
      regn_vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    ele_file >> eatline ;
    fR.push_back( mesh_3Dv::UNSET_FLAG ) ;
  }
}
void mesh2D_reader_GeneralFormat ::
read_regn_fixed_format_ext_bnd_markers( std::ifstream & ele_file, size_t nR, std::vector<size_t> & regn_vlist, std::vector<mesh_3Dv::flag_type> & fR ) {
  size_t id(0), kV(0);
  mesh_3Dv::flag_type regn_flag(mesh_3Dv::UNSET_FLAG) ;
  for ( size_t iR=0 ; iR<nR ; ++iR ) {
    ele_file >> eatcomments >> id ;
    regn_vlist.push_back( nodelem ) ;
    for ( size_t k=0; k<nodelem; ++k ) {
      ele_file >> kV ;
      assert(kV >= offset);
      regn_vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    ele_file >> regn_flag >> eatline ;
    fR.push_back( regn_flag ) ;
  }
  regn_vlist.push_back( nR ) ;
}
void mesh2D_reader_GeneralFormat :: read_regn_data( std::ifstream    & ele_file,
                                                    std::vector<size_t> & regn_vlist,
                                                    std::vector<mesh_3Dv::flag_type> & fR ) {
  // declarations
  size_t nR=0, nbmrk=0 ;

  // read file header
  ele_file >> eatcomments >> nR >> nodelem >> nbmrk >> eatline ;

  // read the vertex region list of mesh
  if      ( nodelem==0 && nbmrk==0 ) { read_regn_free_format_int_bnd_markers ( ele_file, nR, regn_vlist, fR ) ; }
  else if ( nodelem==0 && nbmrk==1 ) { read_regn_free_format_ext_bnd_markers ( ele_file, nR, regn_vlist, fR ) ; }
  else if ( nodelem >0 && nbmrk==0 ) { read_regn_fixed_format_int_bnd_markers( ele_file, nR, regn_vlist, fR ) ; }
  else if ( nodelem >0 && nbmrk==1 ) { read_regn_fixed_format_ext_bnd_markers( ele_file, nR, regn_vlist, fR ) ; }
  else { assert(0) ; }
  regn_vlist.push_back( nR ) ;
}

} // end of namesapce StemMesh3D

#endif //_MESH3D_READERS_HH
