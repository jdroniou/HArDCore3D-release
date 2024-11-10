#ifndef EXTERIOR_CELL_HPP
#define EXTERIOR_CELL_HPP

#include "exterior_evaluation.hpp"
#include "MeshObject.hpp"

namespace HArDCore3D {

  template<size_t space_dim,size_t object_dim>
  class Cell_basis {
    public:
      Cell_basis(const MeshND::MeshObject<space_dim,object_dim> & f, int r);
      Cell_basis() : _r(-1) {;}

      friend class DDR_PEC;
    private:
      template<size_t l> void _cmp_pb_recurse();
      static Eigen::Matrix<double, object_dim, space_dim> _init_trace(const MeshND::MeshObject<space_dim,object_dim> & f);

      int _r;
      double _scale; // 1/h_f, scaling from global chart to element's  
      Eigen::Matrix<double, object_dim, space_dim> _tr;
      Manicore::Monomial_scalar_basis_linear_ID<space_dim,object_dim> _basis;
      std::array<Eigen::MatrixXd,object_dim+1> _exterior_l2;
      // Must be initialized by a friend class
      Eigen::MatrixXd _scalar_mass;
      std::vector<Eigen::MatrixXd> _traces; // Maps the scalar basis on the cell to the one on the boundary
      std::vector<std::array<Eigen::MatrixXd,object_dim>> _ext_traces;

      // the mass is given by kronecker(_exterior_l2[k],_scalar_mass);
      // the mapping to a trace is given by Kronecker(_ext_traces[j][k],_trace[j])
  };

}

#endif

