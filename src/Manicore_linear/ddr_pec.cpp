// Provide the general data needed on a mesh.
// Specifically, it provides the generic operators, the mass matrix on each cell, and the trace operators.
//
// Author: Marien Hanot (marien-lorenzo.hanot@umontpellier.fr)
//


#include "ddr_pec.hpp"

#include "parallel_for.hpp"
#include "boost/multi_array.hpp"
#include "quadraturerule.hpp"
#include "basis.hpp"

#include <iostream>

using namespace HArDCore3D;
using namespace Manicore; 

DDR_PEC::DDR_PEC(Mesh const & mesh,int r, bool use_threads, std::ostream & output) : _r(r) {
  // Initialize global data
  Initialize_exterior_module<3>::init(r);
  _dim_table.resize(3+1);
  _dim_table[0] = 0;
  for (size_t i = 1;i <= 3;++i) {
    _dim_table[i] = _dim_table[i-1] + i+1;
  }
  _list_diff.reserve(_dim_table.back());
  _list_Koszul.reserve(_dim_table.back());
  _list_diff_as_degr.reserve(_dim_table.back());
  _list_trimmed.reserve(_dim_table.back());
  _list_reduced_Koszul_m1.reserve(_dim_table.back());
  _fill_lists<1,0>();

  _nbelem[0] = mesh.n_vertices();
  _nbelem[1] = mesh.n_edges();
  _nbelem[2] = mesh.n_faces();
  _nbelem[3] = mesh.n_cells();
  _1cells.resize(_nbelem[1]);
  _2cells.resize(_nbelem[2]);
  _3cells.resize(_nbelem[3]);

  // Structure to store the computed basis_evaluations 
  std::array<std::vector<boost::multi_array<double,2>>,4> quad_evals;
  quad_evals[1].resize(_nbelem[1]);
  quad_evals[2].resize(_nbelem[2]);
  quad_evals[3].resize(_nbelem[3]);

  ///------------------------------------------------------------------------------------------------------------------------------
  // Construct cells, evaluate quad and compute the gram matrix on every elements
  std::function<void(size_t,size_t)> evaluate_quad_1 =
    [this,&quad_evals,&mesh](size_t start,size_t end)->void {
      for (size_t i1 = start; i1 < end; ++ i1) {
        const Edge &f = *mesh.edge(i1);
        // Construct cell
        _1cells[i1] = {f,_r};
        // Evaluate quad
        QuadratureRule quad_2r_f = generate_quadrature_rule(f,2*_r);
        quad_evals[1][i1].resize(boost::extents[Dimension::PolyDim(_r,1)][quad_2r_f.size()]);
        for (size_t i = 0; i < Dimension::PolyDim(_r,1);++i) {
          for (size_t iqn = 0; iqn < quad_2r_f.size(); ++iqn) {
            quad_evals[1][i1][i][iqn] = _1cells[i1]._basis.evaluate(quad_2r_f[iqn].vector(),i);
          }
        }
        // Fill cell
        _1cells[i1]._scalar_mass = compute_gram_matrix(quad_evals[1][i1],quad_2r_f);
      }
    };
  std::function<void(size_t,size_t)> evaluate_quad_2 =
    [this,&quad_evals,&mesh](size_t start,size_t end)->void {
      for (size_t i2 = start; i2 < end; ++ i2) {
        const Face &f = *mesh.face(i2);
        // Construct cell
        _2cells[i2] = {f,_r};
        // Evaluate quad
        QuadratureRule quad_2r_f = generate_quadrature_rule(f,2*_r);
        quad_evals[2][i2].resize(boost::extents[Dimension::PolyDim(_r,2)][quad_2r_f.size()]);
        for (size_t i = 0; i < Dimension::PolyDim(_r,2);++i) {
          for (size_t iqn = 0; iqn < quad_2r_f.size(); ++iqn) {
            quad_evals[2][i2][i][iqn] = _2cells[i2]._basis.evaluate(quad_2r_f[iqn].vector(),i);
          }
        }
        // Fill cell
        _2cells[i2]._scalar_mass = compute_gram_matrix(quad_evals[2][i2],quad_2r_f);
      }
    };
  std::function<void(size_t,size_t)> evaluate_quad_3 =
    [this,&quad_evals,&mesh](size_t start,size_t end)->void {
      for (size_t i3 = start; i3 < end; ++ i3) {
        const Cell &f = *mesh.cell(i3);
        // Construct cell
        _3cells[i3] = {f,_r};
        // Evaluate quad
        QuadratureRule quad_2r_f = generate_quadrature_rule(f,2*_r);
        quad_evals[3][i3].resize(boost::extents[Dimension::PolyDim(_r,3)][quad_2r_f.size()]);
        for (size_t i = 0; i < Dimension::PolyDim(_r,3);++i) {
          for (size_t iqn = 0; iqn < quad_2r_f.size(); ++iqn) {
            quad_evals[3][i3][i][iqn] = _3cells[i3]._basis.evaluate(quad_2r_f[iqn].vector(),i);
          }
        }
        // Fill cell
        _3cells[i3]._scalar_mass = compute_gram_matrix(quad_evals[3][i3],quad_2r_f);
      }
    };
  output<<"[DDR PEC] Constructing edges masses"<<std::endl;
  parallel_for(_nbelem[1],evaluate_quad_1,use_threads);
  output<<"[DDR PEC] Constructing faces masses"<<std::endl;
  parallel_for(_nbelem[2],evaluate_quad_2,use_threads);
  output<<"[DDR PEC] Constructing cells masses"<<std::endl;
  parallel_for(_nbelem[3],evaluate_quad_3,use_threads);

  ///------------------------------------------------------------------------------------------------------------------------------
  // Compute traces operator
  std::function<void(size_t,size_t)> compute_traces_1 =
    [this,&mesh](size_t start,size_t end)->void {
      for (size_t i1 = start; i1 < end; ++ i1) {
        Edge const & f = *mesh.edge(i1);
        Cell_basis<3,1> & f_cell = _1cells[i1];
        f_cell._traces.reserve(2);
        Eigen::MatrixXd trace(1,Dimension::PolyDim(_r,1));
        for (size_t i = 0; i < Dimension::PolyDim(_r,1);++i) trace(i) = f_cell._basis.evaluate(f.vertex(0)->coords(),i);
        f_cell._traces.emplace_back(trace);
        for (size_t i = 0; i < Dimension::PolyDim(_r,1);++i) trace(i) = f_cell._basis.evaluate(f.vertex(1)->coords(),i);
        f_cell._traces.emplace_back(trace);
        f_cell._ext_traces.emplace_back(std::array<Eigen::MatrixXd,1>{Eigen::Matrix<double,1,1>::Identity()});
        f_cell._ext_traces.emplace_back(std::array<Eigen::MatrixXd,1>{Eigen::Matrix<double,1,1>::Identity()});
      }
    };
  std::function<void(size_t,size_t)> compute_traces_2 =
    [this,&quad_evals,&mesh](size_t start,size_t end)->void {
      for (size_t i2 = start; i2 < end; ++ i2) {
        Face const & f = *mesh.face(i2);
        Cell_basis<3,2> & f_cell = _2cells[i2];
        f_cell._traces.reserve(f.n_edges());
        for (size_t ifp = 0; ifp < f.n_edges();++ifp) {
          // Evaluate quad on trace
          Edge const & fp = *f.edge(ifp);
          Cell_basis<3,1> const & fp_cell = _1cells[fp.global_index()];
          QuadratureRule quad_2r_fp = generate_quadrature_rule(fp,2*_r);
          boost::multi_array<double,2> tr_basis_quad(boost::extents[Dimension::PolyDim(_r,2)][quad_2r_fp.size()]);
          for (size_t i = 0; i < Dimension::PolyDim(_r,2);++i) {
            for (size_t iqn = 0; iqn < quad_2r_fp.size(); ++iqn) {
              tr_basis_quad[i][iqn] = f_cell._basis.evaluate(quad_2r_fp[iqn].vector(),i);
            }
          }
          Eigen::MatrixXd M_tr = compute_gram_matrix(quad_evals[1][fp.global_index()],tr_basis_quad,quad_2r_fp);
          
          f_cell._traces.emplace_back(fp_cell._scalar_mass.ldlt().solve(M_tr));
          // Trace on exterior algebra
          f_cell._ext_traces.emplace_back(std::array<Eigen::MatrixXd,2>{
              Compute_pullback<0,1,2>::compute(f_cell._tr*fp_cell._tr.transpose()*f_cell._scale/fp_cell._scale),
              Compute_pullback<1,1,2>::compute(f_cell._tr*fp_cell._tr.transpose()*f_cell._scale/fp_cell._scale)});
        }
      }
    };
  std::function<void(size_t,size_t)> compute_traces_3 =
    [this,&quad_evals,&mesh](size_t start,size_t end)->void {
      for (size_t i3 = start; i3 < end; ++ i3) {
        Cell const & f = *mesh.cell(i3);
        Cell_basis<3,3> & f_cell = _3cells[i3];
        f_cell._traces.reserve(f.n_faces());
        for (size_t ifp = 0; ifp < f.n_faces();++ifp) {
          // Evaluate quad on trace
          Face const & fp = *f.face(ifp);
          Cell_basis<3,2> const & fp_cell = _2cells[fp.global_index()];
          QuadratureRule quad_2r_fp = generate_quadrature_rule(fp,2*_r);
          boost::multi_array<double,2> tr_basis_quad(boost::extents[Dimension::PolyDim(_r,3)][quad_2r_fp.size()]);
          for (size_t i = 0; i < Dimension::PolyDim(_r,3);++i) {
            for (size_t iqn = 0; iqn < quad_2r_fp.size(); ++iqn) {
              tr_basis_quad[i][iqn] = f_cell._basis.evaluate(quad_2r_fp[iqn].vector(),i);
            }
          }

          Eigen::MatrixXd M_tr = compute_gram_matrix(quad_evals[2][fp.global_index()],tr_basis_quad,quad_2r_fp);
          f_cell._traces.emplace_back(fp_cell._scalar_mass.ldlt().solve(M_tr));
          // Trace on exterior algebra
          f_cell._ext_traces.emplace_back(std::array<Eigen::MatrixXd,3>{
              Compute_pullback<0,2,3>::compute(fp_cell._tr.transpose()*f_cell._scale/fp_cell._scale),
              Compute_pullback<1,2,3>::compute(fp_cell._tr.transpose()*f_cell._scale/fp_cell._scale),
              Compute_pullback<2,2,3>::compute(fp_cell._tr.transpose()*f_cell._scale/fp_cell._scale)});
        }
      }
    };
  output<<"[DDR PEC] Constructing edges traces"<<std::endl;
  parallel_for(_nbelem[1],compute_traces_1,use_threads);
  output<<"[DDR PEC] Constructing faces traces"<<std::endl;
  parallel_for(_nbelem[2],compute_traces_2,use_threads);
  output<<"[DDR PEC] Constructing cells traces"<<std::endl;
  parallel_for(_nbelem[3],compute_traces_3,use_threads);
}

template<size_t d,size_t l> void DDR_PEC::_fill_lists() {
  if constexpr (d <= 3) {
    if constexpr (l <= d) {
      _list_diff.emplace_back(Diff_full<l,d>::get(_r));
      _list_Koszul.emplace_back(Koszul_full<l,d>::get(_r));
      _list_diff_as_degr.emplace_back(Diff_full<l,d>::get_as_degr(_r));
      { // trimmed
        if constexpr(l == 0) {
          _list_trimmed.emplace_back(Eigen::MatrixXd::Identity(Dimension::PLDim(_r,0,d),Dimension::PLtrimmedDim(_r,0,d)));
        } else if (Dimension::PLtrimmedDim(_r,l,d) == 0) {
          _list_trimmed.emplace_back(Eigen::MatrixXd::Zero(Dimension::PLDim(_r,l,d),Dimension::PLtrimmedDim(_r,l,d)));
        } else {
          Eigen::MatrixXd trimmed(Dimension::PLDim(_r,l,d),Dimension::PLDim(_r,l-1,d)+Dimension::PLDim(_r-1,l+1,d));
          trimmed.block(0,0,Dimension::PLDim(_r,l,d),Dimension::PLDim(_r,l-1,d)) = Diff_full<l-1,d>::get_as_degr(_r);
          if (l < d && _r > 0) {
            trimmed.block(0,Dimension::PLDim(_r,l-1,d),Dimension::PLDim(_r,l,d),Dimension::PLDim(_r-1,l+1,d)) = Koszul_full<l+1,d>::get(_r-1);
          } else {
            assert(Dimension::PLDim(_r-1,l+1,d) == 0);
          }
          Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(trimmed);
          assert(lu_decomp.rank() == (int)Dimension::PLtrimmedDim(_r,l,d));
          _list_trimmed.emplace_back(lu_decomp.image(trimmed));
        }
      }
      if (_r > 0 && l > 0) {
        Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(Koszul_full<l,d>::get(_r-1));
        _list_reduced_Koszul_m1.emplace_back(lu_decomp.image(Koszul_full<l,d>::get(_r-1)));
      } else {
        _list_reduced_Koszul_m1.emplace_back(Eigen::Matrix<double,0,0>::Zero());
      }
      _fill_lists<d,l+1>();
    } else {
       _fill_lists<d+1,0>();
    }
  }
}

DDR_PEC::Kronecker DDR_PEC::get_mass(size_t k, size_t d, size_t i) const {
  switch(d) {
    case(1):
      return Eigen::KroneckerProduct(_1cells[i]._exterior_l2[k],_1cells[i]._scalar_mass);
    case(2):
      return Eigen::KroneckerProduct(_2cells[i]._exterior_l2[k],_2cells[i]._scalar_mass);
    case(3):
      return Eigen::KroneckerProduct(_3cells[i]._exterior_l2[k],_3cells[i]._scalar_mass);
    default:
      return Eigen::KroneckerProduct(Eigen::MatrixXd(1,1),Eigen::MatrixXd(1,1));
  }
}

DDR_PEC::Kronecker DDR_PEC::get_trace(size_t k, size_t d, size_t i, size_t j) const {
  switch(d) {
    case(1):
      return Eigen::KroneckerProduct(_1cells[i]._ext_traces[j][k],_1cells[i]._traces[j]);
    case(2):
      return Eigen::KroneckerProduct(_2cells[i]._ext_traces[j][k],_2cells[i]._traces[j]);
    case(3):
      return Eigen::KroneckerProduct(_3cells[i]._ext_traces[j][k],_3cells[i]._traces[j]);
    default:
      return Eigen::KroneckerProduct(Eigen::MatrixXd(1,1),Eigen::MatrixXd(1,1));
  }
}

double DDR_PEC::evaluate_scalar_basis(Eigen::Vector3d const &x, size_t d,size_t i_cell, int i_basis) const {
  assert(i_cell < _nbelem[d] && "Element index too high");
  assert((size_t)i_basis < Dimension::PLDim(_r,0,d) && "Basis index too high");
  switch(d) {
    case(0):
      return 1.;
    case(1):
      return _1cells[i_cell]._basis.evaluate(x,i_basis);
    case(2):
      return _2cells[i_cell]._basis.evaluate(x,i_basis);
    case(3):
      return _3cells[i_cell]._basis.evaluate(x,i_basis);
    default:
      return 0.;
  }
}

Eigen::Matrix<double,-1,-1,0,3,3> DDR_PEC::get_exterior_trace(size_t k, size_t d, size_t i) const {
  assert(d < 4 && k <= d && "Dimension or form degree too high");
  assert(i < _nbelem[d] && "Element index too high");
  if (k == 0) {
    return Eigen::Matrix<double,1,1>::Identity();
  } else {
    if (d == 3) {
      if (k==3) {
        return Eigen::Matrix<double,1,1>::Identity();
      } else {
         return Eigen::Matrix<double,3,3>::Identity();
      }
    } else if (d == 2) {
      if (k == 2) {
        return Compute_pullback<2,2,3>::compute(_2cells[i]._tr.transpose());
      } else { // d=2,k=1
        return _2cells[i]._tr;
      }
    } else if (d == 1) {
      return _1cells[i]._tr;
    } else {
      return Eigen::Matrix<double,0,3>{};
    }
  }
}

Eigen::Matrix<double,-1,-1,0,3,3> DDR_PEC::get_hodge_star(size_t k, size_t d, size_t i) const {
  assert(d < 4 && k <= d && "Dimension or form degree too high");
  assert(i < _nbelem[d] && "Element index too high");
  if (k == 0) { 
    return Eigen::Matrix<double,1,1>{1.};
  } else if (k == d) { 
    return Eigen::Matrix<double,1,1>{1.};
  } else {
    if (d == 2) { // k = 1, basis dx,dy 
      return Eigen::Matrix<double,2,2>{{0,-1.},{1.,0.}};
    } else { // d = 3, k = 1 or 2
      return Eigen::Matrix<double,3,3>{{0.,0.,1.},{0.,-1.,0.},{1.,0.,0.}};
    }
  }
}

Eigen::Matrix<double,-1,1,0,3,1> DDR_PEC::evaluate_basis(Eigen::Vector3d const &x,size_t k,size_t d,size_t i_cell,int i_basis) const {
  assert(i_cell < _nbelem[d] && "Element index too high");
  assert((size_t)i_basis < Dimension::PLDim(_r,k,d) && "Basis index too high");
  assert(k <= d && "Form degree too high for the dimension");
  auto i_div = std::div(i_basis,Dimension::PolyDim(_r,d));
  double scalar_val = evaluate_scalar_basis(x,d,i_cell,i_div.rem);
  // Quickly return for trivial cases
  if (Dimension::ExtDim(k,d)==1)  {
    return Eigen::Matrix<double,1,1>{scalar_val*std::pow(get_scaling(d,i_cell),k)};
  } else {
    Eigen::Matrix<double,-1,1,0,3,1> rv = Eigen::Matrix<double,-1,1,0,3,1>::Zero(Dimension::ExtDim(k,d));
    rv(i_div.quot) = scalar_val*std::pow(get_scaling(d,i_cell),k);
    return rv;
  }
}

Eigen::Matrix<double,-1,1,0,3,1> DDR_PEC::evaluate_basis(Eigen::Vector3d const &x,size_t k,size_t d,size_t i, Eigen::VectorXd const &b) const {
  assert(i < _nbelem[d] && "Element index too high");
  assert(b.size() == (int)Dimension::PLDim(_r,k,d) && "Wrong vector b size");
  assert(k <= d && "Form degree too high for the dimension");
  Eigen::Matrix<double,-1,1,0,3,1> rv = Eigen::Matrix<double,-1,1,0,3,1>::Zero(Dimension::ExtDim(k,d));
  // Itterate all elements
  const size_t extslides = Dimension::PLDim(_r,0,d);
  const double scaling = std::pow(get_scaling(d,i),k);
  for (size_t i_sb = 0; i_sb < Dimension::PLDim(_r,0,d); ++i_sb) { 
    double scalar_val = evaluate_scalar_basis(x,d,i,i_sb)*scaling;
    for (size_t i_ext = 0; i_ext < Dimension::ExtDim(k,d);++i_ext) {
      rv(i_ext) += scalar_val*b(i_sb+i_ext*extslides);
    }
  }
  return rv;
}

double DDR_PEC::get_scaling(size_t d,size_t i) const {
  assert(i < _nbelem[d] && "Element index too high");
  switch(d) {
    case(1):
      return _1cells[i]._scale;
    case(2):
      return _2cells[i]._scale;
    case(3):
      return _3cells[i]._scale;
    default:
      return 0.;
  }
}

