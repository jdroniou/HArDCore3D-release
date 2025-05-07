#ifndef DDR_SPACES_HPP
#define DDR_SPACES_HPP

#include "mesh.hpp"
#include "globaldofspace.hpp"
#include "ddr_pec.hpp"

#include <memory>
#include <variant>

namespace HArDCore3D {

  /*!
   *  \addtogroup Manicore
   * @{
   */

  class DDR_Spaces 
  {
    public:

      /// Constructor
      DDR_Spaces(Mesh const & mesh, int r, bool use_threads = true, std::ostream & output = std::cout);

      struct DDR_function_type {
        typedef std::function<double(const Eigen::Vector3d &)> scalar;
        typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &)> vector;
        std::variant<scalar,vector> func;
        int dqr[4] = {-1,-1,-1,-1}; // Quadrature degree on {vertex (unused),edge,face,cell}. Default to 2*r + 3 if negative 
      };

      /// Return the mesh
      const Mesh & mesh() const
      {
        return _mesh;
      }

      /// Return the polynomial degree
      inline int degree() const {return _r;}

      /// Returns the interpolate of the given function as a k-form
      Eigen::VectorXd interpolate(const DDR_function_type &, size_t k) const;

      /// Returns \star d in PL(r,d-k-1,d), the Hodge star of the full discrete exterior derivative of a discrete k-form on the mesh element of dimension d with index i_cell 
      const Eigen::MatrixXd & full_diff(size_t k, size_t d, size_t i_cell) const; 
      /// Returns \star \ul{d}^k in PLtrimmed(r,d-k-1,d) on the cell (including its boundary), the Hodge star of the discrete exterior derivative of a discrete k-form on the mesh element of dimension d with index i_cell (discrete subspace of X^{k+1})
      Eigen::MatrixXd compose_diff(size_t k, size_t d, size_t i_cell) const; 
      /// Returns \star P^k in PL(r,d-k,d), the Hodge star of the potential of a discrete k-form on the mesh element of dimension d with index i_cell
      const Eigen::MatrixXd &  potential(size_t k, size_t d, size_t i_cell) const; 

      /// Returns the dofspace associated to discrete k-forms 
      inline GlobalDOFSpace const & dofspace(size_t k) const {
        return _dofspace[k];
      }

      /// Evaluate at point x, the sum of the coefficients b times the basis of PL(r,k,d) in the cell with index i_cell 
      Eigen::Matrix<double,-1,1,0,3,1> evaluate_basis(Eigen::Vector3d const &x, size_t k, size_t d, size_t i_cell, Eigen::VectorXd const &b) const;
      
      /// Evaluate at point x, the i th basis function of PL(r,k,d) in the cell with index i_cell
      Eigen::Matrix<double,-1,1,0,3,1> evaluate_basis(Eigen::Vector3d const &x, size_t k, size_t d, size_t i_cell, int i_basis) const;
    
      /// Returns the mass matrix for the k-forms on the mesh element of dimension d with index i_cell
      DDR_PEC::Kronecker get_mass(size_t k, size_t d, size_t i_cell) const;

      /// Gets the continuous Hodge star on k forms on the mesh element of dimension d with index i
      Eigen::Matrix<double,-1,-1,0,3,3> get_hodge_star(size_t k, size_t d, size_t i_cell) const;

      /// Takes a k form on the polynomial k basis to a d-k form on the d-k basis (deals with scaling) on the mesh element of dimension d with index i_cell
      Eigen::MatrixXd discrete_hodge_star(size_t k, size_t d, size_t i_cell) const;

      // Takes the hodge star of a d-k form, i.e. a k form on the k basis to the original d-k form on the d-k basis (deals with scaling and sign change) on the mesh element of dimension d with index i_cell
      /// e.g. discrete_inv_hodge_star(d-k, d, i_cell) * potential(k, d, i_cell) is the potential P^k in PL(r,k,d)
      Eigen::MatrixXd discrete_inv_hodge_star(size_t k, size_t d, size_t i_cell) const;

      /// Trace of a k form from the global space to the mesh element of dimension d with index i_cell
      Eigen::Matrix<double,-1,-1,0,3,3> get_exterior_trace(size_t k, size_t d, size_t i_cell) const;

      /// Computes the L2 product for discrete k-forms on the highest dimensional cell with index i_cell
      Eigen::MatrixXd computeL2Product(size_t k, size_t i_cell, const double penalty_factor = 0.1) const;

      /// Computes the averaged vertex values of a discrete k-form u using the potential
      std::vector<VectorRd> computeVertexValues(size_t k, const Eigen::VectorXd & u) const;

    private:
      // More consistent notation to access, but generate many useless uninitialized entries
      struct DDR_Operators { // one for each form degree k
        std::array<Eigen::MatrixXd,4> full_diff; // equal to \star d, PL(r,d-k-1,d) valued
        std::array<Eigen::MatrixXd,4> diff; // equal to \star \ul{d}, PLtrimmed(r,d-k-1,d) valued
        std::array<Eigen::MatrixXd,4> P; // equal to \star P, PL(r,d-k,d) valued
      };

      // Pointer to the mesh
      const Mesh & _mesh;

      int _r;
      bool _use_threads;
      std::unique_ptr<DDR_PEC> _ddr, _ddr_po;
      std::array<GlobalDOFSpace,3+1> _dofspace; // one for all form degree 0 <= k <= n=3 dimensions
      std::array<std::vector<DDR_Operators>,4> _ops; // one for each cell dimension 0 <= d <= n=3 dimensions
  };
}
#endif

