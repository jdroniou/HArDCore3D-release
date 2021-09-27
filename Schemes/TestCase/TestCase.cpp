// Class to Provides various test cases (diffusion, exact solution, and their derivatives)
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "TestCase.hpp"

using namespace HArDCore3D;

// ----------------------------------------------------------------------------
//                          Implementation
// ----------------------------------------------------------------------------

// Class
TestCase::TestCase(std::vector<int> iTC)
    : m_iTC(iTC),
      m_sol(create_solution(iTC[0])),
      m_diff(create_diffusion(iTC[1]))

{
  if (m_iTC.size() == 2)
  {
    m_iTC.push_back(0);
    m_iTC.push_back(0);
  }
  validate();
  
  // Creation of advection, reaction
  m_advec = create_advection(m_iTC[2]);
  m_reac = create_reaction(m_iTC[3]);

}

/////////////////// SOLUTION ///////////////////////////////:

// Solution
Solution TestCase::create_solution(const size_t n)
{
  FType<double> u = [&](const VectorRd x) -> double {
      return 0;
  };
  CellFType<VectorRd> G = [&](const VectorRd x, const Cell *cell) -> VectorRd {
      return VectorRd::Zero();
  };
  CellFType<MatrixRd> H = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
      return MatrixRd::Zero();
  };

  switch(n)
  {
    /// m_iTC[0]=1: \f$u(x,y,z)=sin(\pi x)  sin(\pi y)  sin (\pi z)\f$
  case 1:
    u = [this](const VectorRd p) -> double {
      return sin(pi * p.x()) * sin(pi * p.y()) * sin(pi * p.z());
    };
    G = [this](const VectorRd p, const Cell *T) -> VectorRd {
      VectorRd val = VectorRd::Zero();
      val(0) = pi * cos(pi * p.x()) * sin(pi * p.y()) * sin(pi * p.z());
      val(1) = pi * sin(pi * p.x()) * cos(pi * p.y()) * sin(pi * p.z());
      val(2) = pi * sin(pi * p.x()) * sin(pi * p.y()) * cos(pi * p.z());
      return val;
    };
    H = [this](const VectorRd p, const Cell *T) -> MatrixRd {
      MatrixRd val = MatrixRd::Zero();
      double x = p.x();
      double y = p.y();
      double z = p.z();
      val.row(0) << -pi * pi * sin(pi * x) * sin(pi * y) * sin(pi * z), pi * pi * cos(pi * x) * cos(pi * y) * sin(pi * z), pi * pi * cos(pi * x) * sin(pi * y) * cos(pi * z);
      val.row(1) << pi * pi * cos(pi * x) * cos(pi * y) * sin(pi * z), -pi * pi * sin(pi * x) * sin(pi * y) * sin(pi * z), pi * pi * sin(pi * x) * cos(pi * y) * cos(pi * z);
      val.row(2) << pi * pi * cos(pi * x) * sin(pi * y) * cos(pi * z), pi * pi * sin(pi * x) * cos(pi * y) * cos(pi * z), -pi * pi * sin(pi * x) * sin(pi * y) * sin(pi * z);
      return val;
    };
    break;

    /// m_iTC[0]=2: \f$u(x,y,z)=cos(\pi x)  cos(\pi y)  cos(\pi z)\f$
  case 2:
    u = [this](const VectorRd p) -> double {
      return cos(pi * p.x()) * cos(pi * p.y()) * cos(pi * p.z());
    };
    G = [this](const VectorRd p, const Cell *T) -> VectorRd {
      VectorRd val = VectorRd::Zero();
      val(0) = -pi * sin(pi * p.x()) * cos(pi * p.y()) * cos(pi * p.z());
      val(1) = -pi * cos(pi * p.x()) * sin(pi * p.y()) * cos(pi * p.z());
      val(2) = -pi * cos(pi * p.x()) * cos(pi * p.y()) * sin(pi * p.z());
      return val;
    };
    H = [this](const VectorRd p, const Cell *T) -> MatrixRd {
      MatrixRd val = MatrixRd::Zero();
      double x = p.x();
      double y = p.y();
      double z = p.z();
      val.row(0) << -pi * pi * cos(pi * x) * cos(pi * y) * cos(pi * z), pi * pi * sin(pi * x) * sin(pi * y) * cos(pi * z), pi * pi * sin(pi * x) * cos(pi * y) * sin(pi * z);
      val.row(1) << pi * pi * sin(pi * x) * sin(pi * y) * cos(pi * z), -pi * pi * cos(pi * x) * cos(pi * y) * cos(pi * z), pi * pi * cos(pi * x) * sin(pi * y) * sin(pi * z);
      val.row(2) << pi * pi * sin(pi * x) * cos(pi * y) * sin(pi * z), pi * pi * cos(pi * x) * sin(pi * y) * sin(pi * z), -pi * pi * cos(pi * x) * cos(pi * y) * cos(pi * z);
      return val;
    };
    break;

    /// m_iTC[0]=3: \f$u(x,y,z)= x\f$
  case 3:
    u = [&](const VectorRd p) -> double {
      return p.x();
    };
    G = [&](const VectorRd p, const Cell *T) -> VectorRd {
      return VectorRd(1.0, 0.0, 0.0);
    };
    break;

    /// m_iTC[0]=4: \f$u(x,y,z)= y\f$
  case 4:
    u = [&](const VectorRd p) -> double {
      return p.y();
    };
    G = [&](const VectorRd p, const Cell *T) -> VectorRd {
      return VectorRd(0.0, 1.0, 0.0);
    };
    break;

    /// m_iTC[0]=5: \f$u(x,y,z)= z\f$
  case 5:
    u = [&](const VectorRd p) -> double {
      return p.z();
    };
    G = [&](const VectorRd p, const Cell *T) -> VectorRd {
      return VectorRd(0.0, 0.0, 1.0);
    };
    break;

    /// m_iTC[0]=6: \f$u(x,y,z)= x^2+y^2+z^2\f$
  case 6:
    u = [&](const VectorRd p) -> double {
      return pow(p.x(), 2) + pow(p.y(), 2) + pow(p.z(), 2);
    };
    G = [&](const VectorRd p, const Cell *T) -> VectorRd {
      return VectorRd(2 * p.x(), 2 * p.y(), 2 * p.z());
    };
    H = [&](const VectorRd p, const Cell *T) -> MatrixRd {
      MatrixRd val = MatrixRd::Zero();
      val.row(0) << 2, 0, 0;
      val.row(1) << 0, 2, 0;
      val.row(2) << 0, 0, 2;
      return val;
    };
    break;

  default:
    break;
  }
  
  return Solution(u, G, H);
}

//////////////////////////// DIFFUSION /////////////////////////////

Diffusion TestCase::create_diffusion(const size_t n)
{
  CellFType<MatrixRd> K = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
      return MatrixRd::Identity();
  };
  CellFType<VectorRd> divK = [&](const VectorRd x, const Cell *cell) -> VectorRd {
      return VectorRd::Zero();
  };
  size_t deg = 0;

  switch(n)
  {
    /// m_iTC[1]=1: Diff = Id
  case 1:
    break;
    
    /// m_iTC[1]=2: Diff = \f$\left[\begin{array}{ccc}y^2+z^2+1 & -xy & -xz\\ -xy & x^2+y^2+1 & -yz\\ -xz & -yz  & x^2+y^2+1\end{array}\right]\f$
  case 2:
    K = [&](const VectorRd p, const Cell *T) -> MatrixRd {
      MatrixRd val = MatrixRd::Zero();
      val.row(0) << pow(p.y(), 2) + pow(p.z(), 2) + 1, -p.x() * p.y(), -p.x() * p.z();
      val.row(1) << -p.x() * p.y(), pow(p.x(), 2) + pow(p.z(), 2) + 1, -p.y() * p.z();
      val.row(2) << -p.x() * p.z(), -p.y() * p.z(), pow(p.x(), 2) + pow(p.y(), 2) + 1;
      return val;
    };
    divK = [&](const VectorRd p, const Cell *T) -> VectorRd {
      return VectorRd(-2 * p.x(), -2 * p.y(), -2 * p.z());
    };
    deg = 2;
    break;
    
  default:
    break;
  }
  
  return Diffusion(K, divK, deg);
}

//////////////////////////// ADVECTION /////////////////////////////

// Advection term
Advection TestCase::create_advection(const size_t n)
{
  // Default is zero
  bool is_zero = true;
  bool is_div_zero = true;
  bool is_div_const = true;
  CellFType<VectorRd> A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
    VectorRd a(0, 0, 0);
    return a;
  };
    CellFType<double> divA = [&](const VectorRd x, const Cell *cell) -> double {
    return 0;
  };

  switch(n)
  {
  case 1:
    is_zero = false;  
    A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
      return VectorRd(1., 1., 1.);
    };
    break;

  case 2:
    is_zero = false;
    A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
      return x;
    };
    
    is_div_zero = false;
    divA = [&](const VectorRd x, const Cell *cell) -> double {
      return 3;
    };
    break;

  case 3:
    is_zero = false;
    A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
      return VectorRd(x(1) * x(2), x(0) * x(2), x(0) * x(1));
    };
    break;

  case 4:
    is_zero = false;
    A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
      double b = pow(x(0), 2) + pow(x(1), 2) + pow(x(2), 2);
      return VectorRd(b, b, b);
    };
    
    is_div_zero = false;
    is_div_const = false;
    divA = [&](const VectorRd x, const Cell *cell) -> double {
      return 2 * x(0) + 2 * x(1) + 2 * x(2);
    };    
    break;

  case 5:
    is_zero = false;
    A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
      return VectorRd(pow(eps, -1), 0, 0);
    };
    break;

  case 6:
    is_zero = false;
    A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
      return VectorRd(exp(x(0)), exp(x(1)), exp(x(2)));
    };

    is_div_zero = false;
    is_div_const = false;
    divA = [&](const VectorRd x, const Cell *cell) -> double {
      return exp(x(0)) + exp(x(1)) + exp(x(2));
    };
    break;

  default:
    break;
  }
  
  return Advection(A, divA, is_zero, is_div_zero, is_div_const);
}

//////////////////////////// REACTION /////////////////////////////

Reaction TestCase::create_reaction(const size_t n)
{
  // Default must be zero
  bool is_zero = true;
  bool is_const = true;
  CellFType<double> R = [&](const VectorRd x, const Cell *cell) -> double {
      return 0.;
  };

  switch (n)
  {
  case 1:
    is_zero = false;
    R = [&](const VectorRd x, const Cell *cell) -> double {
      return 1;
    };
    break;

  case 2:
    is_zero = false;
    is_const = false;
    R = [&](const VectorRd x, const Cell *cell) -> double {
      if (cell->center_mass().y() < 0.5)
      {
        return x(0) + eps;
      }
      else
      {
        return x(1) + eps;
      }
    };
    break;

  case 3:
    is_zero = false;
    is_const = false;
    R = [&](const VectorRd x, const Cell *cell) -> double {
      return x(0) + x(1) + x(2) + eps;
    };
    break;

  case 4:
    is_zero = false;
    is_const = false;
    R = [&](const VectorRd x, const Cell *cell) -> double {
      return exp(x(0) * x(1) * x(2));
    };
    break;

  case 5:
  {
    is_zero = false;
    is_const = false;
    R = [&](const VectorRd x, const Cell *cell) -> double {
      return sin(pi * x(0)) * sin(pi * x(1)) * sin(pi * x(2)) + eps;
    };
    break;
  }

  case 6:
    is_zero = false;
    is_const = false;
    R = [&](const VectorRd x, const Cell *cell) -> double {
      double xcoord = cell->center_mass().x();
      double ycoord = cell->center_mass().y();
      double zcoord = cell->center_mass().z();
      if ((xcoord < 0.5 && ycoord < 0.5 && zcoord < 0.5))
      {
        return 1;
      }
      else
      {
        return x(0) + x(1) + x(2);
      }
    };
    break;

  case 7:
    is_zero = false;
    R = [&](const VectorRd x, const Cell *cell) -> double {
      return eps;
    };
    break;

  case 8:
    is_zero = false;
    R = [&](const VectorRd x, const Cell *cell) -> double {
      return pow(eps, -1);
    };
    break;

  default:
    break;
  }

  return Reaction(R, is_zero, is_const);
}

///////////////////////////// SOURCE TERM ///////////////////////////

// Diffusion source term
CellFType<double> TestCase::diff_source()
{
  CellFType<double> f = [&](const VectorRd x, const Cell *cell) -> double {
      MatrixRd AHu = m_diff.value(x, cell) * m_sol.hessian(x, cell);
      return -AHu.trace() - m_diff.divergence(x, cell).dot(m_sol.gradient(x, cell));
  };

  return f;
}

// Diffusion-Advection-Reaction source term
CellFType<double> TestCase::diff_advec_reac_source()
{
  CellFType<double> f = [&](const VectorRd x, const Cell *cell) -> double {
      return diff_source()(x, cell) 
          + m_advec.divergence(x, cell) * m_sol.value(x) + m_advec.value(x, cell).dot(m_sol.gradient(x, cell)) 
          + m_reac.value(x, cell) * m_sol.value(x);
  };

  return f;
}

///////////////////////////// VALIDATION ////////////////////////////

void TestCase::validate()
{
  if (m_iTC[0] > 6 || m_iTC[0] < 1 || m_iTC[1] > 2 || m_iTC[1] < 1 || m_iTC[2] > 6 || m_iTC[2] < 0 || m_iTC[3] > 8 || m_iTC[3] < 0 || (m_iTC[1] == 3 && m_iTC[0] != 2))
  {
    std::cout << "Incorrect choice of test cases: iTC = " << m_iTC[0] << ", " << m_iTC[1] << ", " << m_iTC[2] << ", " << m_iTC[3] << "\n";
    exit(EXIT_FAILURE);
  }
}
