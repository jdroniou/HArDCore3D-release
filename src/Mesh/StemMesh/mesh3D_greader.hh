#ifndef _MESH_3D_GREADER_HH
#define _MESH_3D_GREADER_HH

#include <cassert>

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "builder.hh"
#include "mesh3D.hh"

#include "mesh3D_format_readers.hh"

namespace StemMesh3D {

//--------------------------------------------------------------------------------------------
// To add a new mesh reader/builder:
// (i)  write the class for reader object
// (ii) write the corresponding method in ExtFileInput
//--------------------------------------------------------------------------------------------
// MESH READER/BUILDER
// objects of this class perform two actions:
// ------------------------------------------
// (i)  read  different mesh formats
// (ii) build the mesh into the input mesh reference
//--------------------------------------------------------------------------------------------

/// Reads mesh files from a variety of file formats and loads them into a mesh object
class ExtFileInput {
private:
  static const size_t C_offset = 0 ; // C/C++ offset
  static const size_t F_offset = 1 ; // fortran offset
  void set_regn_face_list( std::vector<size_t> & regn_flist, std::vector<size_t> & regn_vlist ) ;

public:
  ExtFileInput() {}
  ~ExtFileInput() {}

  // readers for different mesh formats, mesh builders
  /// Read a 3D mesh in the TETGEN format
  void TETGEN_format(
    mesh_3Dv & mesh,                    /**< The mesh object to load the mesh in to */
    std::string fname = "mesh",         /**< The filename of the mesh file */
    size_t offset=C_offset              /**< The indexing offset used in the file (1 = 1-based indexing etc.) */
  );

  /// Read a 3D mesh in the REGN FACE format
  void REGN_FACE_format   (
    mesh_3Dv & mesh,                    /**< The mesh object to load the mesh in to */
    std::string fname = "mesh",         /**< The filename of the mesh file */
    size_t offset=C_offset              /**< The indexing offset used in the file (1 = 1-based indexing etc.) */
  );

  /// Read a 3D mesh in the MSH format
  void MSH_format         (
    mesh_3Dv & mesh,                    /**< The mesh object to load the mesh in to */
    std::string fname = "mesh",         /**< The filename of the mesh file */
    size_t offset=F_offset              /**< The indexing offset used in the file (1 = 1-based indexing etc.) */
  );

  /// Read a general 2D mesh and build a layered 3D mesh from it
  void multi_layer_format (
    mesh_3Dv & mesh,                    /**< The mesh object to load the mesh in to */
    std::string fname = "mesh",         /**< The filename of the mesh file */
    size_t n_layer=1,                   /**< The number of layers to put in the mesh */
    size_t offset=C_offset              /**< The indexing offset used in the file (1 = 1-based indexing etc.) */
  );

} ;
// build the face list as in regn_face format for the TETGEN tetrahedra
void ExtFileInput :: set_regn_face_list( std::vector<size_t> & regn_flist, std::vector<size_t> & regn_vlist ) {

  const size_t nRF = 4 ;
  const size_t nFV = 3 ;

  std::vector<size_t> vlist(nRF) ;
  size_t perm[nRF][nFV] = { {0,2,1}, {1,2,3}, {0,1,3}, {0,3,2} } ;

  size_t nR  = regn_vlist.back() ;

  size_t k = 0;
  for ( size_t iR=0; iR<nR; ++iR ) { // loop on the regions

    size_t nRV = regn_vlist[k++] ;
    assert( nRV==nRF ) ;

    for ( size_t ilV=0; ilV<nRV; ++ilV ) {
      vlist[ilV] = regn_vlist[k++] ;
    }

    regn_flist.push_back(nRF) ;                    // number of faces of the region iR
    for ( size_t ilF=0; ilF<nRF; ++ilF ) {
      regn_flist.push_back(nFV) ;                  // number of vertices of the face ilF
      for ( size_t ip=0; ip<3; ++ip ) {
        regn_flist.push_back( vlist[ perm[ilF][ip] ] ) ;
      }
    }
  }
  regn_flist.push_back( nR ) ;
}
//--------------------------------------------------------------------------------------------
void ExtFileInput :: TETGEN_format( mesh_3Dv & mesh, std::string fname, size_t offset ) {
  std::vector<double> xV, yV, zV ;
  std::vector<mesh_3Dv::flag_type> fV, fR ;
  std::vector<size_t> regn_vlist, regn_flist ;
  std::vector<size_t> face_vlist, edge_vlist ;

  mesh3D_reader_TETGEN_format mesh_reader(std::move(fname),offset) ;
  mesh_reader.read_the_mesh( xV, yV, zV, fV, regn_vlist, fR ) ;
  //print_logfile( xV, yV, zV, fV, regn_vlist, fR, offset ) ;

  // build the list of region faces for tetrahedral meshes 
  // (to be used when data is input from TETGEN)
  set_regn_face_list( regn_flist, regn_vlist ) ;

  //  ---- start builder
  mesh3Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, zV, fV, regn_flist, fR ) ; // !!!
}
//--------------------------------------------------------------------------------------------
void ExtFileInput :: REGN_FACE_format( mesh_3Dv & mesh, std::string fname, size_t offset ) {
  std::vector<double> xV, yV, zV ;
  std::vector<mesh_3Dv::flag_type> fV, fR ;
  std::vector<size_t> regn_flist ;

  mesh3D_reader_REGN_FACE_format mesh_reader(std::move(fname),offset) ;
  mesh_reader.read_the_mesh( xV, yV, zV, fV, regn_flist, fR ) ;
  //print_regn_face( xV, yV, zV, fV, regn_flist, fR, offset ) ;

  //  ---- start builder
  mesh3Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, zV, fV, regn_flist, fR ) ; // !!!
}
//--------------------------------------------------------------------------------------------
void ExtFileInput :: MSH_format( mesh_3Dv & mesh, std::string fname, size_t offset ) {
  std::vector<double> xV, yV, zV ;
  std::vector<mesh_3Dv::flag_type> fV, fR ;
  std::vector<size_t> regn_flist ;

  mesh3D_reader_MSH_format mesh_reader(std::move(fname),offset) ;
  mesh_reader.read_the_mesh( xV, yV, zV, fV, regn_flist, fR ) ;
  //print_regn_face( xV, yV, zV, fV, regn_flist, fR, offset ) ;

  //  ---- start builder
  mesh3Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, zV, fV, regn_flist, fR ) ; // !!!
}
//--------------------------------------------------------------------------------------------
void ExtFileInput :: multi_layer_format( mesh_3Dv & mesh, std::string fname, size_t n_layer, size_t offset ) {
  std::vector<double> xV_2D, yV_2D  ;
  std::vector<mesh_3Dv::flag_type>    fV_2D, fR_2D  ;
  std::vector<size_t>    regn_vlist_2D ;

  mesh2D_reader_GeneralFormat mesh_reader(std::move(fname),offset) ;
  mesh_reader.read_mesh( xV_2D, yV_2D, fV_2D, regn_vlist_2D, fR_2D ) ;

  std::vector<double> xV, yV, zV ;
  std::vector<size_t>    regn_flist;
  std::vector<mesh_3Dv::flag_type> fV, fR ;

  size_t    nL = n_layer ;
  double dz = 1./double(nL) ;

  // randomized planes
  std::vector<double> rnd_x(nL+1), rnd_y(nL+1) ;

  /*
  bool rnd_planar_coordinate_Z = tilted_planes ;
  if ( rnd_planar_coordinate_Z ) {
    for ( size_t jL=0; jL<=nL; ++jL ) {
      //double fac = double(jL*(nL-jL)) / (double(nL*nL)/4.) ;

      double rnd_dz = 0.4*dz*(0.5-double(rand())/double(RAND_MAX)) ;
      assert(  -dz/2.<rnd_dz && rnd_dz<dz/2. ) ;

      rnd_x[jL] = jL==0||jL==nL ? 0. : rnd_dz ; //  0.45*fac*dz ;
      rnd_y[jL] = jL==0||jL==nL ? 0. : rnd_dz ; // -0.45*fac*dz ;
    }
  }
  */

  // VERTEX DATA
  size_t nV_2D = xV_2D.size() ;
  for ( size_t jL=0; jL<=nL; ++jL ) {
    for ( size_t iV=0; iV<nV_2D; ++iV ) {
      xV.push_back( xV_2D[iV] ) ;
      yV.push_back( yV_2D[iV] ) ;
      // ------------------------------------------------------------------------
      double z_vrtx = dz*double(jL) + rnd_x[jL]*xV_2D[iV] + rnd_y[jL]*yV_2D[iV] ;
      zV.push_back( z_vrtx ) ;
      // ------------------------------------------------------------------------
      fV.push_back( fV_2D[iV] ) ;
    }
  }
  assert( xV.size() == nV_2D*(nL+1) ) ;

  // REGION DATA
  // loop on region layers
  size_t nR_2D = regn_vlist_2D.back() ;
  for ( size_t iL=0; iL<nL; ++iL ) {

    size_t k = 0 ;
    for ( size_t iR=0; iR<nR_2D; ++iR ) {
      // get data of 2D region iR
      size_t nRV = regn_vlist_2D[k++] ;
      std::vector<size_t> tmp_vlist ;
      for ( size_t j=0; j<nRV; ++j ) { tmp_vlist.push_back( regn_vlist_2D[k++] ) ; }

      // set fR for the region iR + iL*nR
      fR.push_back( iL ) ;

      // number of faces = nRV+2
      regn_flist.push_back( nRV+2 ) ;

      // number of vertices of bottom face
      regn_flist.push_back( nRV ) ;
      for ( size_t j=nRV-1; j>=0; --j ) {
        regn_flist.push_back( tmp_vlist[j]+nV_2D*iL   ) ;
      }

      // number of vertices of top face
      regn_flist.push_back( nRV ) ;
      for ( size_t j=0; j<nRV; ++j ) {
        regn_flist.push_back( tmp_vlist[j]+nV_2D*(iL+1) ) ;
      }

      // loop on vertical faces
      for ( size_t j=0; j<nRV; ++j ) {
        size_t j1 = (j+1)%nRV ;
        regn_flist.push_back( 4 ) ;
        regn_flist.push_back( tmp_vlist[j ]+nV_2D*iL     ) ;
        regn_flist.push_back( tmp_vlist[j1]+nV_2D*iL     ) ;
        regn_flist.push_back( tmp_vlist[j1]+nV_2D*(iL+1) ) ;
        regn_flist.push_back( tmp_vlist[j ]+nV_2D*(iL+1) ) ;
      }

    } // end loop on regions
  } // end loop on region's layers
  regn_flist.push_back( nR_2D*nL ) ;

  // output the mesh log file
  //print_regn_face( xV, yV, zV, fV, regn_flist, fR, offset ) ;

  //  ---- start builder
  mesh3Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, zV, fV, regn_flist, fR ) ; // !!!
}

} // end of namespace StemMesh3D

#endif // end of _MESH_3D_GREADER_HH
