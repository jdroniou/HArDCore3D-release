#ifndef LIE_ALGEBRA
#define LIE_ALGEBRA

#include <Eigen/Dense>
#include <vector>

namespace HArDCore3D
{
  
  /*!
   *	\addtogroup LADDRCore
   * @{
   */

  /// Lie algebra class: mass matrix, structure constants and Lie bracket
  class LieAlgebra /// Assuming orthonormal basis
  {
  public:

    typedef Eigen::MatrixXcd LieAlgValue;
    typedef std::function<double(LieAlgValue &, LieAlgValue &)> LieAlgProduct; 

    /// Constructors
    LieAlgebra();
    LieAlgebra(std::vector<LieAlgValue> & basis, LieAlgProduct & product); /// Assuming orthonormal basis

    /// Returns the dimension of the Lie algebra
    inline const size_t dimension() const
    {
      return m_basis.size();
    }

    /// Returns the basis of the Lie algebra
    inline const std::vector<LieAlgValue> & basis() const
    {
      return m_basis;
    }

    /// Returns the Gram matrix of the Lie algebra
    inline const Eigen::MatrixXd & massMatrix() const
    {
      return m_mass_matrix;
    }

    /// Computes the Lie bracket of two Lie algebra elements
    inline const LieAlgValue lieBracket(const LieAlgValue & A, const LieAlgValue & B) const
    {
      return A*B - B*A;
    }
    /// Computes the structure constants of the Lie algebra
    inline const std::vector<Eigen::MatrixXd> & structureConst() const
    {
      return m_strucConst;
    }

  private:

    /// Computes the mass matrix of the Lie algebra 
    Eigen::MatrixXd _compute_mass_matrix();
    /// Computes the structure constants of the Lie algebra
    std::vector<Eigen::MatrixXd> _compute_structure_constants();

    std::vector<LieAlgValue> m_basis;
    LieAlgProduct m_product;
    Eigen::MatrixXd m_mass_matrix;
    std::vector<Eigen::MatrixXd> m_strucConst;
  };

}

#endif
