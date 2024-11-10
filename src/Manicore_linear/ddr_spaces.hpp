#ifndef DDR_SPACES_HPP
#define DDR_SPACES_HPP

#include "mesh.hpp"
#include "globaldofspace.hpp"

#include <memory>
#include <variant>

namespace HArDCore3D {

  class DDR_PEC;
  class DDR_Spaces {
    public:
      DDR_Spaces(Mesh const & mesh, int r, bool use_threads = true, std::ostream & output = std::cout);

      struct DDR_function_type {
        typedef std::function<double(const Eigen::Vector3d &)> scalar;
        typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> vector;
        std::variant<scalar,vector> func;
        int dqr[4] = {-1,-1,-1,-1}; /// Quadrature degree on {vertex (unused),edge,face,cell}. Default to 2*r + 3 if negative 
      };
      /// Return the interpolate of the given function as a k-form
      Eigen::VectorXd interpolate(const DDR_function_type &,size_t k) const;

      const Eigen::MatrixXd & full_diff(size_t k, size_t d,size_t i) const; // Return \star d in PL(r,d-k-1,d)
      Eigen::MatrixXd compose_diff(size_t k, size_t d,size_t i) const; // Return \star \ul{d}^k in PLtrimmed(r,d-k-1,d) on the cell (including its boundary)
      const Eigen::MatrixXd & potential(size_t k, size_t d,size_t i) const; // Return \star P^k in PL(r,d-k,d)

      inline GlobalDOFSpace const & dofspace(size_t k) const {
        return _dofspace[k];
      }
      inline int degree() const {return _r;}

      // Evaluate PL(r,k,d) in the cell i, at x 
      Eigen::Matrix<double,-1,1,0,3,1> evaluate_basis(Eigen::Vector3d const &x,size_t k,size_t d,size_t i, Eigen::VectorXd const &b) const;


    private:
      // More consistent notation to access, but generate many useless uninitialized entries
      struct DDR_Operators { // one for each form degree k
        std::array<Eigen::MatrixXd,4> full_diff; // equal to \star d, PL(r,d-k-1,d) valued
        std::array<Eigen::MatrixXd,4> diff; // equal to \star \ul{d}, PLtrimmed(r,d-k-1,d) valued
        std::array<Eigen::MatrixXd,4> P; // equal to \star P, PL(r,d-k,d) valued
      };

      int _r;
      bool _use_threads;
      std::unique_ptr<DDR_PEC> _ddr, _ddr_po;
      std::array<GlobalDOFSpace,3+1> _dofspace; // one for all form degree
      std::array<std::vector<DDR_Operators>,4> _ops; // one for each cell dimension
  };
}
#endif

