// Creates quadrature rule on an edge
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
*
*  This library was developed around HHO methods, although some parts of it have a more
* general purpose. If you use this code or part of it in a scientific publication, 
* please mention the following book as a reference for the underlying principles
* of HHO schemes:
*
* The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications. 
*  D. A. Di Pietro and J. Droniou. Modeling, Simulation and Applications, vol. 19. 
*  Springer International Publishing, 2020, xxxi + 525p. doi: 10.1007/978-3-030-37203-3. 
*  url: https://hal.archives-ouvertes.fr/hal-02151813.
*
*/



#ifndef QUAD1D_HPP
#define QUAD1D_HPP
#include<stddef.h>

namespace HArDCore3D {  // forward dec

class LegendreGauss;
}

/*!
*	@addtogroup Quadratures
* @{
*/


namespace HArDCore3D {

class QuadRuleEdge {

public:
    QuadRuleEdge(size_t doe, bool warn);
    ~QuadRuleEdge();

    size_t nq();
    double xq(size_t i);
    double yq(size_t i);
    double zq(size_t i);
    double wq(size_t i);
    void setup(double xV[], double yV[], double zV[]);

private:
    size_t _npts;
    LegendreGauss* _rule;
    double* _w;
    double* _xr;
    double* _yr;
    double* _zr;
    double _length;
};

class LegendreGauss {
public:
    LegendreGauss(size_t doe);
    ~LegendreGauss();

    void sub_rule_01();
    void sub_rule_02();
    void sub_rule_03();
    void sub_rule_04();
    void sub_rule_05();
    void sub_rule_06();
    void sub_rule_07();
    void sub_rule_08();
    void sub_rule_09();
    void sub_rule_10();
    void sub_rule_11();
    size_t npts();
    double wq(size_t i);
    double tq(size_t i);

private:
    size_t _doe;
    size_t _npts;
    double* _t;
    double* _w;
};
//@}
}

#endif /* QUAD1D_HPP */
