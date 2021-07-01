// Creates quadrature rule on an edge
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include <quad1d.hpp>
#include <max_degrees_quadratures.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <iomanip>

using namespace HArDCore3D;
QuadRuleEdge::QuadRuleEdge(size_t doe, bool warn):
  _rule(0),
  _w(0),
  _xr(0),
  _yr(0),
  _zr(0),
  _length(0){
  _rule = new LegendreGauss(std::max(int(doe),1));
  _npts = _rule->npts();
}

QuadRuleEdge::~QuadRuleEdge(){
    delete[] _w;
    delete[] _xr;
    delete[] _yr;
    delete[] _zr;
    delete _rule;
}

size_t QuadRuleEdge::nq() {
    return _npts;
};
double QuadRuleEdge::xq(size_t i) { return _xr[i]; }
double QuadRuleEdge::yq(size_t i) { return _yr[i]; }
double QuadRuleEdge::zq(size_t i) { return _zr[i]; }
double QuadRuleEdge::wq(size_t i) {
    return _rule->wq(i)*_length;
};
void QuadRuleEdge::setup(double x[], double y[], double z[]) {
    _xr = new double[_npts];
    _yr = new double[_npts];
    _zr = new double[_npts];
    double v[3] = {x[1] - x[0], y[1] - y[0], z[1]-z[0]};
    _length = sqrt(std::pow(v[0],2)+std::pow(v[1],2)+std::pow(v[2],2));
    for (size_t i = 0; i < _npts; i++) {
        // map the q points to the line
        _xr[i] = _rule->tq(i) * v[0] + x[0];
        _yr[i] = _rule->tq(i) * v[1] + y[0];
        _zr[i] = _rule->tq(i) * v[2] + z[0];
    }
}

LegendreGauss::LegendreGauss(size_t doe) : _doe(doe),
    _npts(std::ceil((_doe + 1)/ 2.0)),
  _t(0),
  _w(0){
    switch (_npts) {
        case 1:
            sub_rule_01();
      break;
        case 2:
            sub_rule_02();
      break;
        case 3:
            sub_rule_03();
      break;
    case 4:
            sub_rule_04();
      break;
        case 5:
            sub_rule_05();
      break;
        case 6:
            sub_rule_06();
      break;
        case 7:
            sub_rule_07();
      break;
        case 8:
            sub_rule_08();
      break;
        case 9:
            sub_rule_09();
      break;
        case 10:
            sub_rule_10();
      break;
    case 11:
      sub_rule_11();
      break;
    default:
      throw "Can't integrate edge to degree ";
    return;
    }
}

LegendreGauss::~LegendreGauss(){
  delete[] _t;
  delete[] _w;
}

void LegendreGauss::sub_rule_01() {
    _w = new double[1];
    _t = new double[1];
    _t[0] = 0.50000;
    _w[0] = 1.0;
}
void LegendreGauss::sub_rule_02() {
    _w = new double[2];
    _t = new double[2];

    _t[0] = 0.5+sqrt(3.0)/6.0;
    _t[1] = 0.5-sqrt(3.0)/6.0;
    _w[0] = 0.5;
    _w[1] = 0.5;
}
void LegendreGauss::sub_rule_03() {
    _w = new double[3];
    _t = new double[3];

    _t[0] = 0.5+ sqrt(15.0)/10.0;
    _t[1] = 0.50000;
    _t[2] = 0.5- sqrt(15.0)/10.0;
    _w[0] = 5.0/18;
    _w[1] = 8.0/18;
    _w[2] = 5.0/18;
}
void LegendreGauss::sub_rule_04() {
    _w = new double[4];
    _t = new double[4];

    _t[0] = 0.5+sqrt(525+70.0*sqrt(30))/70.0;
    _t[1] = 0.5+sqrt(525-70.0*sqrt(30))/70.0;
    _t[2] = 0.5-sqrt(525-70.0*sqrt(30))/70.0;
    _t[3] = 0.5-sqrt(525+70.0*sqrt(30))/70.0;
    _w[0] = (18-sqrt(30.0))/72.0;
    _w[1] = (18+sqrt(30.0))/72.0;
    _w[2] = (18+sqrt(30.0))/72.0;
    _w[3] = (18-sqrt(30.0))/72.0;
}
void LegendreGauss::sub_rule_05() {
    _w = new double[5];
    _t = new double[5];

    _t[0] = 0.5+0.9061798459386640/2.0;
    _t[1] = 0.5+0.5384693101056831/2.0;
    _t[2] = 0.500000;
    _t[3] = 0.5-0.5384693101056831/2.0;
    _t[4] = 0.5-0.9061798459386640/2.0;

    _w[0] = 0.2369268850561891/2.0;
    _w[1] = 0.4786286704993665/2.0;
    _w[2] = 0.5688888888888889/2.0;
    _w[3] = 0.4786286704993665/2.0;
    _w[4] = 0.2369268850561891/2.0;
}
void LegendreGauss::sub_rule_06() {
    _w = new double[6];
    _t = new double[6];

    _t[0] = 0.5+0.9324695142031521/2.0;
    _t[1] = 0.5+0.6612093864662645/2.0;
    _t[2] = 0.5+0.2386191860831969/2.0;
    _t[3] = 0.5-0.2386191860831969/2.0;
    _t[4] = 0.5-0.6612093864662645/2.0;
    _t[5] = 0.5-0.9324695142031521/2.0;

    _w[0] = 0.1713244923791704/2.0;
    _w[1] = 0.3607615730481386/2.0;
    _w[2] = 0.4679139345726910/2.0;
    _w[3] = 0.4679139345726910/2.0;
    _w[4] = 0.3607615730481386/2.0;
    _w[5] = 0.1713244923791704/2.0;
}
void LegendreGauss::sub_rule_07() {
    _w = new double[7];
    _t = new double[7];

    _t[0] = 0.5+0.9491079123427585/2.0;
    _t[1] = 0.5+0.7415311855993945/2.0;
    _t[2] = 0.5+0.4058451513773972/2.0;
    _t[3] = 0.500000;
    _t[4] = 0.5-0.4058451513773972/2.0;
    _t[5] = 0.5-0.7415311855993945/2.0;
    _t[6] = 0.5-0.9491079123427585/2.0;

    _w[0] = 0.1294849661688697/2.0;
    _w[1] = 0.2797053914892766/2.0;
    _w[2] = 0.3818300505051189/2.0;
    _w[3] = 0.4179591836734694/2.0;
    _w[4] = 0.3818300505051189/2.0;
    _w[5] = 0.2797053914892766/2.0;
    _w[6] = 0.1294849661688697/2.0;
}
void LegendreGauss::sub_rule_08() {
    _w = new double[8];
    _t = new double[8];

    _t[0] = 0.5+0.9602898564975363/2.0;
    _t[1] = 0.5+0.7966664774136267/2.0;
    _t[2] = 0.5+0.5255324099163290/2.0;
    _t[3] = 0.5+0.1834346424956498/2.0;
    _t[4] = 0.5-0.1834346424956498/2.0;
    _t[5] = 0.5-0.5255324099163290/2.0;
    _t[6] = 0.5-0.7966664774136267/2.0;
    _t[7] = 0.5-0.9602898564975363/2.0;

    _w[0] = 0.1012285362903763/2.0;
    _w[1] = 0.2223810344533745/2.0;
    _w[2] = 0.3137066458778873/2.0;
    _w[3] = 0.3626837833783620/2.0;
    _w[4] = 0.3626837833783620/2.0;
    _w[5] = 0.3137066458778873/2.0;
    _w[6] = 0.2223810344533745/2.0;
    _w[7] = 0.1012285362903763/2.0;
}
void LegendreGauss::sub_rule_09() {
    _w = new double[9];
    _t = new double[9];

    _t[0] = 0.5+0.9681602395076261/2.0;
    _t[1] = 0.5+0.8360311073266358/2.0;
    _t[2] = 0.5+0.6133714327005904/2.0;
    _t[3] = 0.5+0.3242534234038089/2.0;
    _t[4] = 0.500000;
    _t[5] = 0.5-0.3242534234038089/2.0;
    _t[6] = 0.5-0.6133714327005904/2.0;
    _t[7] = 0.5-0.8360311073266358/2.0;
    _t[8] = 0.5-0.9681602395076261/2.0;

    _w[0] = 0.0812743883615744/2.0;
    _w[1] = 0.1806481606948574/2.0;
    _w[2] = 0.2606106964029354/2.0;
    _w[3] = 0.3123470770400029/2.0;
    _w[4] = 0.3302393550012598/2.0;
    _w[5] = 0.3123470770400029/2.0;
    _w[6] = 0.2606106964029354/2.0;
    _w[7] = 0.1806481606948574/2.0;
    _w[8] = 0.0812743883615744/2.0;
}
void LegendreGauss::sub_rule_10() {
    _w = new double[10];
    _t = new double[10];

    _t[0] = 0.5+0.9739065285171717/2.0;
    _t[1] = 0.5+0.8650633666889845/2.0;
    _t[2] = 0.5+0.6794095682990244/2.0;
    _t[3] = 0.5+0.4333953941292472/2.0;
    _t[4] = 0.5+0.1488743389816312/2.0;
    _t[5] = 0.5-0.1488743389816312/2.0;
    _t[6] = 0.5-0.4333953941292472/2.0;
    _t[7] = 0.5-0.6794095682990244/2.0;
    _t[8] = 0.5-0.8650633666889845/2.0;
    _t[9] = 0.5-0.9739065285171717/2.0;

    _w[0] = 0.0666713443086881/2.0;
    _w[1] = 0.1494513491505806/2.0;
    _w[2] = 0.2190863625159820/2.0;
    _w[3] = 0.2692667193099963/2.0;
    _w[4] = 0.2955242247147529/2.0;
    _w[5] = 0.2955242247147529/2.0;
    _w[6] = 0.2692667193099963/2.0;
    _w[7] = 0.2190863625159820/2.0;
    _w[8] = 0.1494513491505806/2.0;
    _w[9] = 0.0666713443086881/2.0;
}
void LegendreGauss::sub_rule_11() {
    _w = new double[11];
    _t = new double[11];

    _t[0] = 0.5+0.9782286581460570/2.0;
    _t[1] = 0.5+0.8870625997680953/2.0;
    _t[2] = 0.5+0.7301520055740494/2.0;
    _t[3] = 0.5+0.5190961292068118/2.0;
    _t[4] = 0.5+0.2695431559523450/2.0;
    _t[5] = 0.500000000000000000000000;
    _t[6] = 0.5-0.2695431559523450/2.0;
    _t[7] = 0.5-0.5190961292068118/2.0;
    _t[8] = 0.5-0.7301520055740494/2.0;
    _t[9] = 0.5-0.8870625997680953/2.0;
  _t[10] = 0.5-0.9782286581460570/2.0;


    _w[0] = 0.0556685671161737/2.0;
    _w[1] = 0.1255803694649046/2.0;
    _w[2] = 0.1862902109277343/2.0;
    _w[3] = 0.2331937645919905/2.0;
    _w[4] = 0.2628045445102467/2.0;
    _w[5] = 0.2729250867779006/2.0;
    _w[6] = 0.2628045445102467/2.0;
    _w[7] = 0.2331937645919905/2.0;
    _w[8] = 0.1862902109277343/2.0;
    _w[9] = 0.1255803694649046/2.0;
  _w[10] = 0.0556685671161737/2.0;
}



size_t LegendreGauss::npts() { return _npts; }
double LegendreGauss::wq(size_t i) {
    if (i >= _npts) {
        throw "Trying to access quadrature point that is greater than number computed";
    }
    return _w[i];
}

double LegendreGauss::tq(size_t i) {
    if (i >= _npts) {
        throw "Trying to access quadrature point that is greater than number computed";
    }
    return _t[i];
}
