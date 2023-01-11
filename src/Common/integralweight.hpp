// Structure to provide weight functions for integrals: provides value, and polynomial degree (cell by cell)

/*
 *
 *      This library was developed around HHO methods, although some parts of it have a more
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

#ifndef INTEGRALWEIGHT_HPP
#define INTEGRALWEIGHT_HPP

#include <mesh.hpp>

namespace HArDCore3D
{

  /*!
   *	\addtogroup Common
   * @{
   */
   
    typedef std::function<double (const Cell &, const Eigen::Vector3d &)> IntegralWeightValueType;
    typedef std::function<size_t (const Cell &)> IntegralWeightDegreeType;
   
    /// Structure for weights (scalar, at the moment) in integral.
    /** Each weight is represented as a piecewise function defined cell-by-cell, together with its local polynomial degree to determine the offset for quadrature rules (degree=0 means that the weight is constant and some calculations are easier) */
    struct IntegralWeight
    {
        /// Generic constructor
        IntegralWeight(
              const IntegralWeightValueType _value, ///< Value of weight
              const IntegralWeightDegreeType _deg ///< Local degree of weight
              )
        : value(_value),
          deg(_deg)
        {
          // Do nothing
        }

       /// Constructor for constant weights
       IntegralWeight(double val)
        : IntegralWeight( [val](const Cell &T, const Eigen::Vector3d &x)->double {return val;}, [](const Cell &T)->size_t {return 0;} )
       {
        // Do nothing
       }

       /// Constructor when the dependency on the cell T is not explicit in the value (and degree is constant)
       IntegralWeight(const std::function<double(const VectorRd &)> & val, const size_t & deg)
        : IntegralWeight( [val](const Cell &T, const Eigen::VectorXd &x)->double {return val(x);}, [deg](const Cell &T)->size_t {return deg;} )
       {
        // Do nothing
       }

       /// Constructor when the dependency on the cell T is not explicit in the value, and degree is not provided (it is assumed to be 0)
       IntegralWeight(const std::function<double(const VectorRd &)> & val)
        : IntegralWeight(val, 0)
       {
        // Do nothing
       }
        
      IntegralWeightValueType value;
      IntegralWeightDegreeType deg;
    };

    /// Operator to multiply an IntegralWeight by a number
    inline IntegralWeight operator* (double const & r, IntegralWeight const & weight)
    {
      IntegralWeightValueType r_times_value 
          = [weight, r](const Cell & T, const Eigen::Vector3d & x)->double { return r * weight.value(T,x);};
      IntegralWeightDegreeType deg = weight.deg;
      
      return IntegralWeight(r_times_value, deg);
    }

    /// Operator to add an IntegralWeight to another one
    inline IntegralWeight operator+ (IntegralWeight const & weight1, IntegralWeight const & weight2)
    {
      IntegralWeightValueType plus_value 
          = [weight1, weight2](const Cell & T, const Eigen::Vector3d & x)->double 
                    { return weight1.value(T,x) +  weight2.value(T,x);};
      IntegralWeightDegreeType max_deg 
          = [weight1, weight2](const Cell & T)->size_t
                    { return std::max(weight1.deg(T), weight2.deg(T)); };
      
      return IntegralWeight(plus_value, max_deg);
    }

  //@}

} // end of namespace HArDCore3D

#endif
