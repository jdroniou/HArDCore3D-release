#ifndef DDR_PEC_HPP
#define DDR_PEC_HPP

#include "exterior_cell.hpp"
#include "mesh.hpp"

namespace HArDCore3D {
  
  class DDR_PEC {
    public:
      DDR_PEC(Mesh const & mesh,int r,bool use_threads = true, std::ostream & output = std::cout);

      typedef Eigen::KroneckerProduct<Eigen::MatrixXd,Eigen::MatrixXd> Kronecker;
      // Return the mass matrix for the k-forms on the i-th d-cell 
      Kronecker get_mass(size_t k, size_t d, size_t i) const;
      // Return the trace for the k-forms on the i-th d-cell onto its j-th (d-1)-neighbour
      Kronecker get_trace(size_t k, size_t d, size_t i, size_t j) const;

      // Getter for the generic operators matrices
      const Eigen::MatrixXd & get_diff(size_t l, size_t d) const {return _list_diff[_cmp_ind(l,d)];}
      const Eigen::MatrixXd & get_Koszul(size_t l, size_t d) const {return _list_Koszul[_cmp_ind(l,d)];}
      const Eigen::MatrixXd & get_diff_as_degr(size_t l, size_t d) const {return _list_diff_as_degr[_cmp_ind(l,d)];}
      const Eigen::MatrixXd & get_trimmed(size_t l, size_t d) const {return _list_trimmed[_cmp_ind(l,d)];}
      const Eigen::MatrixXd & get_reduced_Koszul_m1(size_t l, size_t d) const {return _list_reduced_Koszul_m1[_cmp_ind(l,d)];}

      // Trace from the global space to the element
      Eigen::Matrix<double,-1,-1,0,3,3> get_exterior_trace(size_t k, size_t d, size_t i) const;
      Eigen::Matrix<double,-1,-1,0,3,3> get_hodge_star(size_t k, size_t d, size_t i) const;

      //Gives the evaluation into the global basis (before the scaling).
      double evaluate_scalar_basis(Eigen::Vector3d const &x, size_t d,size_t i_cell, int i_basis) const;
      Eigen::Matrix<double,-1,1,0,3,1> evaluate_basis(Eigen::Vector3d const &x,size_t k,size_t d,size_t i_cell,int i_basis) const;
      Eigen::Matrix<double,-1,1,0,3,1> evaluate_basis(Eigen::Vector3d const &x,size_t k,size_t d,size_t i, Eigen::VectorXd const &b) const;

      double get_scaling(size_t d,size_t i) const;

    private:
      template<size_t d,size_t l> void _fill_lists();
      inline int _cmp_ind(size_t l, size_t d) const {return _dim_table[d-1]+l;}

      int _r;
      std::array<size_t,4> _nbelem;
      std::vector<Cell_basis<3,1>> _1cells;
      std::vector<Cell_basis<3,2>> _2cells;
      std::vector<Cell_basis<3,3>> _3cells;
      std::vector<size_t> _dim_table;
      std::vector<Eigen::MatrixXd> _list_diff, // Image of dPL(r,l,d) inside PL(r-1,l+1,d)
                                   _list_Koszul, // Image of kPL(r,l,d) inside PL(r+1,l-1,d)
                                   _list_diff_as_degr, // Image of dPL(r,l,d) inside PL(r,l+1,d)
                                   _list_trimmed, // Image of PLtrimmed(r,l,d)
                                   _list_reduced_Koszul_m1; // Image of kPL(r-1,l,d)

  };

}
#endif

