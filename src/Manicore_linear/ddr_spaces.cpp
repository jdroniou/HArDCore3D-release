// Provide the operator needed in DDR-PEC.
// Specifically, it provides the interpolator, the potential and the differential operator.
//
// Author: Marien Hanot (marien-lorenzo.hanot@umontpellier.fr)
//


#include "ddr_spaces.hpp"
#include "ddr_pec.hpp"

#include "parallel_for.hpp"
#include "quadraturerule.hpp"

using namespace HArDCore3D;
using namespace Manicore;

DDR_Spaces::DDR_Spaces(Mesh const & mesh, int r, bool use_threads, std::ostream & output) 
  : _r(r), _use_threads(use_threads), 
    _ddr(std::make_unique<DDR_PEC>(mesh,_r,_use_threads,output)),
    _ddr_po(std::make_unique<DDR_PEC>(mesh,_r+1,_use_threads,output)),
    _dofspace(std::array<GlobalDOFSpace,4>{
        GlobalDOFSpace{mesh,Dimension::PLtrimmedDim(_r,0-0,0),Dimension::PLtrimmedDim(_r,1-0,1),Dimension::PLtrimmedDim(_r,2-0,2),Dimension::PLtrimmedDim(_r,3-0,3)},
        GlobalDOFSpace{mesh,0,Dimension::PLtrimmedDim(_r,1-1,1),Dimension::PLtrimmedDim(_r,2-1,2),Dimension::PLtrimmedDim(_r,3-1,3)},
        GlobalDOFSpace{mesh,0,0,Dimension::PLtrimmedDim(_r,2-2,2),Dimension::PLtrimmedDim(_r,3-2,3)},
        GlobalDOFSpace{mesh,0,0,0,Dimension::PLtrimmedDim(_r,3-3,3)}})
  {
    // Resize ops
    _ops[0].resize(mesh.n_vertices());
    _ops[1].resize(mesh.n_edges());
    _ops[2].resize(mesh.n_faces());
    _ops[3].resize(mesh.n_cells());

    // Init P for d = k
    std::function<void(size_t,size_t,size_t)> init_P0_gen =
      [this](size_t start,size_t end,size_t k)->void {
        for (size_t i = start; i < end; ++i) {
          _ops[k][i].P[k] = Eigen::MatrixXd::Identity(Dimension::PolyDim(_r,k),Dimension::PolyDim(_r,k));
        }
      };
    // Init d for d >= k + 1
    std::function<void(size_t,size_t,size_t,size_t)> init_d_gen =
      [this](size_t start,size_t end,size_t k,size_t d)->void { // d^k_{r,f(d)}
        assert(d-k>0 && "init_d_gen called with d-k <= 0");
        for (size_t i = start; i < end; ++i) {
          // Dimension dependent setup
          auto mass_mkmo = _ddr->get_mass(d-k-1,d,i);
          auto mass_mk = _ddr->get_mass(d-k,d,i);
          size_t dofdim = _dofspace[k].dimension(d,i);
          std::vector<size_t> boundary = _dofspace[k].mesh().get_boundary(d,i);
          // RHS
          Eigen::MatrixXd RHS = Eigen::MatrixXd::Zero(Dimension::PLDim(_r,d-k-1,d),dofdim);
          // (-1)^k+1 \int_f w_f ^ du_f
          if (Dimension::PLtrimmedDim(_r,d-k,d) > 0) {
            RHS.rightCols(Dimension::PLtrimmedDim(_r,d-k,d)) = (((k+1)%2 == 0)? 1. : -1.) *_ddr->get_diff_as_degr(d-k-1,d).transpose()*mass_mk*_ddr->get_trimmed(d-k,d);
          }
          // Sum over f
          for (size_t j = 0; j < boundary.size(); ++j) {
            if (d==1) {
              RHS.col(j) = (j%2==0?-1.:1.)*_ddr->get_trace(d-k-1,d,i,j).transpose();
            } else {
              RHS += _dofspace[k].mesh().boundary_orientation(d,i,j)*
                    _ddr->get_trace(d-k-1,d,i,j).transpose()*_ddr->get_mass(d-k-1,d-1,boundary[j])*
                    _dofspace[k].extendOperator(d,i,d-1,boundary[j],_ops[d-1][boundary[j]].P[k]);
            }
          }
          _ops[d][i].full_diff[k] = mass_mkmo.ldlt().solve(RHS);
          assert(_ops[d][i].full_diff[k].rows() == (int)Dimension::PLDim(_r,d-k-1,d) && 
              _ops[d][i].full_diff[k].cols() == (int)((d == 1)? _dofspace[k].dimensionEdge(i) :
                                                 (d == 2)? _dofspace[k].dimensionFace(i) :
                                                           _dofspace[k].dimensionCell(i)) && "Wrong dimension for full_diff operator");

        }
      };
    // Init P for d >= k + 1
    std::function<void(size_t,size_t,size_t,size_t)> init_P_gen =
      [this](size_t start,size_t end,size_t k,size_t d)->void { // d^k_{r,f(d)}
        assert(d-k>0 && "init_d_gen called with d-k <= 0");
        for (size_t i = start; i < end; ++i) {
          // Dimension dependent setup
          auto mass_mk = _ddr->get_mass(d-k,d,i);
          auto mass_rpo_dmkmo = _ddr_po->get_mass(d-k-1,d,i);
          size_t dofdim = _dofspace[k].dimension(d,i);
          std::vector<size_t> boundary = _dofspace[k].mesh().get_boundary(d,i);
          auto as_degrpo = Eigen::KroneckerProduct(
                Eigen::MatrixXd::Identity(Dimension::ExtDim(d-k-1,d),Dimension::ExtDim(d-k-1,d)),
                Eigen::MatrixXd::Identity(Dimension::PolyDim(_r+1,d),Dimension::PolyDim(_r,d)));
          const size_t dim_uf = Dimension::kPLDim(_r,d-k,d);
          const size_t dim_vf = Dimension::kPLDim(_r-1,d-k+1,d);
          assert(Dimension::PLDim(_r,d-k,d) == dim_uf+dim_vf && "Wrong size for test functions in P");
          // LHS
          Eigen::MatrixXd LHS(Dimension::PLDim(_r,d-k,d),Dimension::PLDim(_r,d-k,d));
          LHS.topRows(dim_uf) = (_ddr_po->get_diff(d-k-1,d)*_ddr_po->get_reduced_Koszul_m1(d-k,d)).transpose();
          if (dim_vf > 0) {
            LHS.bottomRows(dim_vf) = _ddr->get_reduced_Koszul_m1(d-k+1,d).transpose();
          }
          LHS = LHS*(((k+1)%2==0)? 1.: -1.)* mass_mk;

          // RHS
          Eigen::MatrixXd RHS = Eigen::MatrixXd::Zero(Dimension::PLDim(_r,d-k,d),dofdim);
          RHS.topRows(dim_uf) = _ddr_po->get_reduced_Koszul_m1(d-k,d).transpose()
                               *mass_rpo_dmkmo*as_degrpo*_ops[d][i].full_diff[k];
          if (dim_vf > 0) {
            RHS.bottomRightCorner(dim_vf,Dimension::PLtrimmedDim(_r,d-k,d)) = (((k+1)%2==0)? 1.: -1.) * 
                                  _ddr->get_reduced_Koszul_m1(d-k+1,d).transpose()*
                                  mass_mk*_ddr->get_trimmed(d-k,d);
          }
          for (size_t j = 0; j < boundary.size(); ++j) {
            auto as_degrpo_tr = Eigen::KroneckerProduct(
                Eigen::MatrixXd::Identity(Dimension::ExtDim(d-k-1,d-1),Dimension::ExtDim(d-k-1,d-1)),
                Eigen::MatrixXd::Identity(Dimension::PolyDim(_r+1,d-1),Dimension::PolyDim(_r,d-1)));
            if (d==1) {
              RHS.block(0,j,dim_uf,1) -= (j%2==0?-1.:1.)
                                          *_ddr_po->get_reduced_Koszul_m1(d-k,d).transpose()
                                          *_ddr_po->get_trace(d-k-1,d,i,j).transpose();
            } else {
                RHS.topRows(dim_uf) -= _dofspace[k].mesh().boundary_orientation(d,i,j)
                 *_ddr_po->get_reduced_Koszul_m1(d-k,d).transpose()
                 *_ddr_po->get_trace(d-k-1,d,i,j).transpose()*_ddr_po->get_mass(d-k-1,d-1,boundary[j])
                 *as_degrpo_tr*_dofspace[k].extendOperator(d,i,d-1,boundary[j],_ops[d-1][boundary[j]].P[k]);
            }
          }
          _ops[d][i].P[k] = LHS.partialPivLu().solve(RHS);
        }
      };

    std::function<void(size_t,size_t,size_t,size_t)> init_projd_gen =
      [this](size_t start,size_t end,size_t k,size_t d)->void { // d^k_{r,f(d)}
        for (size_t i = start; i < end; ++i) {
          // Dimension dependent setup
          assert(d-k > 0 && "Projector called on wrong dimension/forms");
          auto mass_mkmo = _ddr->get_mass(d-k-1,d,i);
          // LHS
          Eigen::MatrixXd trimmed = _ddr->get_trimmed(d-k-1,d);
          Eigen::MatrixXd LHS = trimmed.transpose()*mass_mkmo*trimmed;
          // RHS
          Eigen::MatrixXd RHS = trimmed.transpose()*mass_mkmo*_ops[d][i].full_diff[k];
          _ops[d][i].diff[k] = LHS.ldlt().solve(RHS);
        }
      };

    output<<"[DDR Spaces] Initializing potentials for d=k"<<std::endl;
    for (size_t i = 0; i <= 3; ++i) {
      parallel_for(_ops[i].size(),std::bind(init_P0_gen,std::placeholders::_1,std::placeholders::_2,i),_use_threads);
    }
    output<<"[DDR Spaces] Initializing diff for d=k+1"<<std::endl;
    for (size_t i = 1; i <= 3; ++i) {
      parallel_for(_ops[i].size(),std::bind(init_d_gen,std::placeholders::_1,std::placeholders::_2,i-1,i),_use_threads);
    }
    output<<"[DDR Spaces] Initializing potentials for d=k+1"<<std::endl;
    for (size_t i = 1; i <= 3; ++i) {
      parallel_for(_ops[i].size(),std::bind(init_P_gen,std::placeholders::_1,std::placeholders::_2,i-1,i),_use_threads);
    }
    output<<"[DDR Spaces] Initializing diff for d=k+2"<<std::endl;
    for (size_t i = 2; i <= 3; ++i) {
      parallel_for(_ops[i].size(),std::bind(init_d_gen,std::placeholders::_1,std::placeholders::_2,i-2,i),_use_threads);
    }
    output<<"[DDR Spaces] Initializing potentials for d=k+2"<<std::endl;
    for (size_t i = 2; i <= 3; ++i) {
      parallel_for(_ops[i].size(),std::bind(init_P_gen,std::placeholders::_1,std::placeholders::_2,i-2,i),_use_threads);
    }
    output<<"[DDR Spaces] Initializing diff for d=k+3"<<std::endl;
    for (size_t i = 3; i <= 3; ++i) {
      parallel_for(_ops[i].size(),std::bind(init_d_gen,std::placeholders::_1,std::placeholders::_2,i-3,i),_use_threads);
    }
    output<<"[DDR Spaces] Initializing potentials for d=k+3"<<std::endl;
    for (size_t i = 3; i <= 3; ++i) {
      parallel_for(_ops[i].size(),std::bind(init_P_gen,std::placeholders::_1,std::placeholders::_2,i-3,i),_use_threads);
    }
    // Init diff
    output<<"[DDR Spaces] Initializing ul diff"<<std::endl;
    for (size_t i = 0; i <= 3; ++i) {
      for (size_t j = 0; j < i; ++j) {
        parallel_for(_ops[i].size(),std::bind(init_projd_gen,std::placeholders::_1,std::placeholders::_2,j,i),_use_threads);
      }
    }
  }

Eigen::VectorXd DDR_Spaces::interpolate(const DDR_function_type & func, size_t k) const {
  Eigen::VectorXd qh(_dofspace[k].dimension());
  GlobalDOFSpace const & dofsp = _dofspace[k];
  // General l2_projection
  auto l2_proj = [this,func,k](QuadratureRule const & quad,Eigen::MatrixXd const & mass,size_t d,size_t i_cell)->Eigen::VectorXd {
    if (Dimension::PLtrimmedDim(_r,d-k,d) == 0) return Eigen::VectorXd::Zero(0);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(mass.cols());
    const bool is_scalar = (func.func.index() == 0);
    for (size_t iqn = 0; iqn < quad.size(); iqn++) {
      Eigen::Matrix<double,-1,1,0,3,1> fv;
      if (is_scalar) {
        fv = Eigen::Matrix<double,1,1>{std::get<DDR_function_type::scalar>(func.func)(quad[iqn].vector())};
      } else {
        fv = std::get<DDR_function_type::vector>(func.func)(quad[iqn].vector());
      }
      fv = _ddr->get_hodge_star(k,d,i_cell)*_ddr->get_exterior_trace(k,d,i_cell)*fv;
      for (int i = 0; i < mass.cols(); ++i) {
        b(i) += quad[iqn].w*fv.dot(_ddr->evaluate_basis(quad[iqn].vector(),d-k,d,i_cell,i));
      }
    }
    Eigen::MatrixXd trimmed = _ddr->get_trimmed(d-k,d);

    return (trimmed.transpose()*mass*trimmed).ldlt().solve(trimmed.transpose()*b);
  };

  // Interpolate at vertices
  auto interpolate_vertices = [&qh,dofsp,func](size_t start,size_t end)->void{
    assert(dofsp.numLocalDofsVertex() == 1 && "Interpolation at vertices only done for scalar field");
    for (size_t i = start; i < end; ++i) {
      const Vertex & V = *dofsp.mesh().vertex(i);
      qh[dofsp.globalOffset(V)] = std::get<DDR_function_type::scalar>(func.func)(V.coords());
    }
  };
  // Interpolate at cells
  auto interpolate_cells = [this,&qh,dofsp,func,k,l2_proj](size_t d, size_t start,size_t end)->void{
    assert(1 <= d && d <= 3);
    for (size_t i_cell = start; i_cell < end; ++i_cell) {
      const int dqr = (func.dqr[1] >= 0)? func.dqr[d] : 2*_r + 3;
      QuadratureRule quad_dqr_f = (d == 1) ? generate_quadrature_rule(*dofsp.mesh().edge(i_cell),dqr) :
                                  (d == 2) ? generate_quadrature_rule(*dofsp.mesh().face(i_cell),dqr) :
                                             generate_quadrature_rule(*dofsp.mesh().cell(i_cell),dqr);
      Eigen::MatrixXd mass = _ddr->get_mass(d-k,d,i_cell);
      qh.segment(dofsp.globalOffset(d,i_cell),Dimension::PLtrimmedDim(_r,d-k,d)) = l2_proj(quad_dqr_f,mass,d,i_cell);
    }
  };

  switch(k) { // fallthrough
    case(0):
      parallel_for(dofsp.mesh().n_elems(0),interpolate_vertices,_use_threads);
    case(1):
      parallel_for(dofsp.mesh().n_elems(1),std::bind(interpolate_cells,1,std::placeholders::_1,std::placeholders::_2),_use_threads);
    case(2):
      parallel_for(dofsp.mesh().n_elems(2),std::bind(interpolate_cells,2,std::placeholders::_1,std::placeholders::_2),_use_threads);
    case(3):
      parallel_for(dofsp.mesh().n_elems(3),std::bind(interpolate_cells,3,std::placeholders::_1,std::placeholders::_2),_use_threads);
    default:
      ;
  }
  return qh;
}

Eigen::MatrixXd DDR_Spaces::compose_diff(size_t k,size_t d,size_t i) const {
  assert(d < 4 && k < d && i < _ops[d].size() && "Access of diff out of range");
  assert(_dofspace[k+1].numLocalDofs(0) == 0 && "Expected no dofs on vertices");
  switch(d) {
    case(1):
      return _ops[1][i].diff[k];
    case(2):
      {
        const Face &F = *_dofspace[k].mesh().face(i);
        Eigen::MatrixXd diff(_dofspace[k+1].dimensionFace(i),_dofspace[k].dimensionFace(i));
        if (_dofspace[k+1].numLocalDofs(1) > 0) {
          for (size_t j = 0; j < F.n_edges(); ++j) {
            diff.middleRows(_dofspace[k+1].localOffset(F,*F.edge(j)),_dofspace[k+1].numLocalDofs(1)) = 
                _dofspace[k].extendOperator(F,*F.edge(j),_ops[1][F.edge(j)->global_index()].diff[k]);
          }
        }
        if (_dofspace[k+1].numLocalDofs(2) > 0) {
          diff.bottomRows(_dofspace[k+1].numLocalDofs(2)) = _ops[2][i].diff[k];
        }
        return diff;
      }
    case(3):
      {
        const Cell &T = *_dofspace[k].mesh().cell(i);
        Eigen::MatrixXd diff(_dofspace[k+1].dimensionCell(i),_dofspace[k].dimensionCell(i));
        if (_dofspace[k+1].numLocalDofs(1) > 0) {
          for (size_t j = 0; j < T.n_edges(); ++j) {
            diff.middleRows(_dofspace[k+1].localOffset(T,*T.edge(j)),_dofspace[k+1].numLocalDofs(1)) = 
                _dofspace[k].extendOperator(T,*T.edge(j),_ops[1][T.edge(j)->global_index()].diff[k]);
          }
        }
        if (_dofspace[k+1].numLocalDofs(2) > 0) {
          for (size_t j = 0; j < T.n_faces(); ++j) {
            diff.middleRows(_dofspace[k+1].localOffset(T,*T.face(j)),_dofspace[k+1].numLocalDofs(2)) = 
                _dofspace[k].extendOperator(T,*T.face(j),_ops[2][T.face(j)->global_index()].diff[k]);
          }
        }
        if (_dofspace[k+1].numLocalDofs(3) > 0) {
          diff.bottomRows(_dofspace[k+1].numLocalDofs(3)) = _ops[3][i].diff[k];
        }
        return diff;
      }
    default:
      return Eigen::Matrix<double,0,0>{};
  }
}

const Eigen::MatrixXd & DDR_Spaces::full_diff(size_t k, size_t d, size_t i) const {
  assert(d < 4 && k < d && i < _ops[d].size() && "Access of potential out of range");
  return _ops[d][i].full_diff[k];
}

const Eigen::MatrixXd & DDR_Spaces::potential(size_t k, size_t d, size_t i) const {
  assert(d < 4 && k <= d && i < _ops[d].size() && "Access of potential out of range");
  return _ops[d][i].P[k];
}

Eigen::Matrix<double,-1,1,0,3,1> DDR_Spaces::evaluate_basis(Eigen::Vector3d const &x,size_t k,size_t d,size_t i,Eigen::VectorXd const &b) const {
  return _ddr->evaluate_basis(x,k,d,i,b);
}

