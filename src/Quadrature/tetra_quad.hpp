/*
  CREDITS:
  Gianmarco Manzini
  IMATI-CNR
  via Ferrata 1
  27100 Pavia

  DATE: Nov 29, 2008

	SEPARATE hpp/cpp implementation by Jerome Droniou (jerome.droniou@monash.edu)

  PREVIOUS CREDITS: this set of routines is a C++ re-implementation of
  part of John Burkardt's software for geometric calculations (2005), which
  is freely available on Internet (see his homepage)

  COPY-LEFT: Feel free to use this software under GNU Public Licence (GPL) terms.
*/

#ifndef _TETRA_QUAD_HPP
#define _TETRA_QUAD_HPP

#include <cassert>
#include <cmath>
#include <ctime>
#include <cstdlib>

#include <iostream>
#include <iomanip>

namespace HArDCore3D {

/*!
*	@addtogroup Quadratures
* @{
*/

// -------------------------------------------------------------------------------------------
// Volume of a tetrahedron in 3D.
// input : ( x[0:3], y[0:3], z[0:3] ) coordinates of the tetrahedron vertices.
double tetra_volume(double x[4], double y[4], double z[4]);

// -------------------------------------------------------------------------------------------
// Apply a quadrature rule in a GENERIC tetrahedron
// func ( double x, double y, double z ), the function to be integrated.
// x[0:3], y[0:3], z[0:3],                the coordinates of the tetrahedron vertices.
// nq,                                    the number of quadrature nodes
// cq0[0:nq-1], cq1[0:nq-1], cq2[0:nq-1], the first three barycentric coordinates of the quad nodes
double tetra_sum(double func(double x, double y, double z),
                 double x[4], double y[4], double z[4],
                 int nq, double cq0[], double cq1[], double cq2[], double wq[]);

// -------------------------------------------------------------------------------------------
// Given the unit tetrahedron { (x,y,z) | 0<=x && 0<=y && 0<=z && x+y+z<=1 }, it sets
// - cq0, cq1, cq2, the first three barycentric coordinates of the quadrature nodes
// - wq,            the weights the weights of the quadrature nodes
// - nq,            the number of nodes
//
// References:
//
//  Hermann Engels,
//    Numerical Quadrature and Cubature,
//    Academic Press, 1980,
//
//  Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Olgierd Zienkiewicz,
//    The Finite Element Method,
//    Sixth Edition,
//    Butterworth-Heinemann, 2005,
//
size_t tetra_unit_size(size_t rule);

void tetra_unit_set(size_t rule, size_t nq, double cq0[], double cq1[], double cq2[], double wq[]);

// -------------------------------------------------------------------------------------------
// Volume of a (3D) parallelipiped
// input: ( x[0:3], y[0:3], z[0:3] ) coordinates of one corner and
//                                   its three immediate neighbors.
double parallelipiped_volume_3d(double x[4], double y[4], double z[4]);

// -------------------------------------------------------------------------------------------
// Volume of the 3D unit tetrahedron: { (x,y,z) | 0<=x && 0<=y && 0<=z && x+y+z<=1 }
double tetra_unit_volume();

// -------------------------------------------------------------------------------------------
// Apply a quadrature rule to func(x,y,z) in the UNIT tetrahedron.
// func(x,y,z)          : the function to be integrated
// nq                   : the number of nodes
// xq[iq],yq[iq],zq[iq] : the coordinates of the quad node iq=0..nq-1
// wq[iq]               : the weight      of the quad node iq=0..nq-1
double tetra_unit_sum(double func(double x, double y, double z),
                      size_t nq, double xq[], double yq[], double zq[], double wq[]);

} // end of namespace HArDCore3D

#endif // end of _TETRA_QUAD_HPP
