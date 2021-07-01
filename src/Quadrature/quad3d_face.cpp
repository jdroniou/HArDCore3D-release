
#include <quad3d_face.hpp>
#include <max_degrees_quadratures.hpp>
#include <iostream>

using namespace HArDCore3D;

/*!
*	@addtogroup Quadratures
* @{
*/

QuadRuleTriangle::QuadRuleTriangle( size_t doe, bool warn ) 
		: _npts(dunavant_order_num(std::max(int(doe),1))), 
			qsx(0),
			qsy(0),
			qsz(0),
			qswg(0),
			_xy(0),
			_w(0) {

  if (doe > MAX_DOE_FACE && warn) {std::cerr << "Warning: quadrature rule of degree " << doe  <<
                                       " requested, but the maximum available is degree " << MAX_DOE_FACE << std::endl;}
       qsx  = new double[_npts] ;
       qsy  = new double[_npts] ;
       qsz  = new double[_npts] ;
       qswg = new double[_npts] ;
	_w = new double[_npts];
        _xy = new double[_npts * 2];
    //    Input, int RULE, the index of the rule.
    //    Input, int ORDER_NUM, the order (number of points) of the rule.
    //    Output, double XY[2*ORDER_NUM], the points of the rule.
    //    Output, double W[ORDER_NUM], the weights of the rule.
    	dunavant_rule(std::max(int(doe),1), _npts, _xy, _w);
}
QuadRuleTriangle::~QuadRuleTriangle() {
  delete[] qsx  ;
  delete[] qsy  ;
  delete[] qsz  ;
  delete[] qswg ;
  delete[] _w;
  delete[] _xy;

}


double QuadRuleTriangle::area( double xV[], double yV[], double zV[] ) const {
  // set the edge vectors
  double dx[2] = { xV[1]-xV[0], xV[2]-xV[0] } ;
  double dy[2] = { yV[1]-yV[0], yV[2]-yV[0] } ;
  double dz[2] = { zV[1]-zV[0], zV[2]-zV[0] } ;
  // area
  return sqrt( pow(dy[0]*dz[1]-dy[1]*dz[0],2) +
               pow(dz[0]*dx[1]-dz[1]*dx[0],2) +
               pow(dx[0]*dy[1]-dx[1]*dy[0],2) )/2. ;
}


// --------------------------------------------------------------------------------------------


// implementation for any triangle
// Dunavant rule
// integrate polynomials of degree <= 20 EXACTLY

void QuadRuleTriangle::setup( double xV[], double yV[], double zV[] ) {
const double area_tri = area( xV, yV, zV ) ;

	for(size_t i=0; i < _npts; i++){ // from the barycentric coordinates, obtain the x,y,z coordinates, and also the quadrature weight in 3D
	qsx[i] = _xy[2*i]*xV[0]+_xy[2*i+1]*xV[1]+(1.0-_xy[2*i]-_xy[2*i+1])*xV[2];
	qsy[i] = _xy[2*i]*yV[0]+_xy[2*i+1]*yV[1]+(1.0-_xy[2*i]-_xy[2*i+1])*yV[2];
	qsz[i] = _xy[2*i]*zV[0]+_xy[2*i+1]*zV[1]+(1.0-_xy[2*i]-_xy[2*i+1])*zV[2];
	qswg[i] = _w[i] * area_tri;
	}

}

size_t QuadRuleTriangle::nq() {return _npts;}
double QuadRuleTriangle::xq( size_t i ) {  return qsx[i]  ; }
double QuadRuleTriangle::yq( size_t i ) { return qsy[i]  ; }
double QuadRuleTriangle::zq( size_t i ) {  return qsz[i]  ; }
double QuadRuleTriangle::wq( size_t i ) { return qswg[i] ; }




