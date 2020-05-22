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
      _deg_diff(0)
{
  if (m_iTC.size() == 2)
  {
    m_iTC.push_back(0);
    m_iTC.push_back(0);
  }
  validate();
  if (m_iTC[1] == 2)
  {
    _deg_diff = 2;
  }
}

/////////////////// SOLUTION ///////////////////////////////:

// Solution
FType<double> TestCase::sol()
{
  FType<double> u = [&](const VectorRd p) -> double { return 0; };
  switch(m_iTC[0])
  {
    /// m_iTC[0]=1: \f$u(x,y,z)=sin(\pi x)  sin(\pi y)  sin (\pi z)\f$
  case 1:
    u = [this](const VectorRd p) -> double {
      return sin(pi * p.x()) * sin(pi * p.y()) * sin(pi * p.z());
    };
    break;
    /// m_iTC[0]=2: \f$u(x,y,z)=cos(\pi x)  cos(\pi y)  cos(\pi z)\f$
  case 2:
    u = [this](const VectorRd p) -> double {
      return cos(pi * p.x()) * cos(pi * p.y()) * cos(pi * p.z());
    };
    break;
    /// m_iTC[0]=3: \f$u(x,y,z)= x\f$
  case 3:
    u = [&](const VectorRd p) -> double {
      return p.x();
    };
    break;
    /// m_iTC[0]=4: \f$u(x,y,z)= y\f$
  case 4:
    u = [&](const VectorRd p) -> double {
      return p.y();
    };
    break;
    /// m_iTC[0]=5: \f$u(x,y,z)= z\f$
  case 5:
    u = [&](const VectorRd p) -> double {
      return p.z();
    };
    break;
    /// m_iTC[0]=6: \f$u(x,y,z)= x^2+y^2+z^2\f$
  case 6:
    u = [&](const VectorRd p) -> double {
      return pow(p.x(), 2) + pow(p.y(), 2) + pow(p.z(), 2);
    };
    break;

  default:
    break;
  }
  return u;
}

// Gradient of the solution
CellFType<VectorRd> TestCase::grad_sol()
{
  CellFType<VectorRd> G = [&](const VectorRd p, const Cell *T) -> VectorRd { return VectorRd::Zero(); };
  switch(m_iTC[0])
  {
  case 1:
    G = [this](const VectorRd p, const Cell *T) -> VectorRd {
      VectorRd val = VectorRd::Zero();
      val(0) = pi * cos(pi * p.x()) * sin(pi * p.y()) * sin(pi * p.z());
      val(1) = pi * sin(pi * p.x()) * cos(pi * p.y()) * sin(pi * p.z());
      val(2) = pi * sin(pi * p.x()) * sin(pi * p.y()) * cos(pi * p.z());
      return val;
    };
    break;

  case 2:
    G = [this](const VectorRd p, const Cell *T) -> VectorRd {
      VectorRd val = VectorRd::Zero();
      val(0) = -pi * sin(pi * p.x()) * cos(pi * p.y()) * cos(pi * p.z());
      val(1) = -pi * cos(pi * p.x()) * sin(pi * p.y()) * cos(pi * p.z());
      val(2) = -pi * cos(pi * p.x()) * cos(pi * p.y()) * sin(pi * p.z());
      return val;
    };
    break;

  case 3:
    G = [&](const VectorRd p, const Cell *T) -> VectorRd {
      return VectorRd(1.0, 0.0, 0.0);
    };
    break;

  case 4:
    G = [&](const VectorRd p, const Cell *T) -> VectorRd {
      return VectorRd(0.0, 1.0, 0.0);
    };
    break;

  case 5:
    G = [&](const VectorRd p, const Cell *T) -> VectorRd {
      return VectorRd(0.0, 0.0, 1.0);
    };
    break;

  case 6:
    G = [&](const VectorRd p, const Cell *T) -> VectorRd {
      return VectorRd(2 * p.x(), 2 * p.y(), 2 * p.z());
    };
    break;

  default:
    break;
  }
  return G;
}

// Hessian of the solution
CellFType<MatrixRd> TestCase::hess_sol()
{
  CellFType<MatrixRd> H = [&](const VectorRd p, const Cell *T) -> MatrixRd { return MatrixRd::Zero(); };
  switch(m_iTC[0])
  {

  case 1:
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

  case 2:
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

  case 3:
    break;
  case 4:
    break;
  case 5:
    break;

  case 6:
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
  return H;
}

//////////////////////////// DIFFUSION /////////////////////////////

// Diffusion matrix
CellFType<MatrixRd> TestCase::diff()
{
  CellFType<MatrixRd> K = [&](const VectorRd p, const Cell *T) -> MatrixRd { return MatrixRd::Identity(); };
  switch(m_iTC[1])
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
    break;
  default:
    break;
  }
  return K;
}

// Divergence by row of the diffusion matrix
CellFType<VectorRd> TestCase::div_diff()
{
  CellFType<VectorRd> divK = [&](const VectorRd p, const Cell *T) -> VectorRd { return VectorRd::Zero(); };
  switch(m_iTC[1])
  {
  case 1:
    break;
  case 2:
    divK = [&](const VectorRd p, const Cell *T) -> VectorRd {
      return VectorRd(-2 * p.x(), -2 * p.y(), -2 * p.z());
    };
    break;
  default:
    break;
  }
  return divK;
}

CellFType<VectorRd> TestCase::advec()
{
  // Default must be zero
  CellFType<VectorRd> A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
    VectorRd a(0, 0, 0);
    return a;
  };
  switch (m_iTC[2])
  {
  case 1:
    A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
      VectorRd a(1, 1, 1);
      return a;
    };
    break;
  case 2:
    A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
      return x;
    };
    break;
  case 3:
    A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
      VectorRd a(x(1) * x(2), x(0) * x(2), x(0) * x(1));
      return a;
    };
    break;
  case 4:
    A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
      double b = pow(x(0), 2) + pow(x(1), 2) + pow(x(2), 2);
      VectorRd a(b, b, b);
      return a;
    };
    break;
  case 5:
    A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
      VectorRd a(pow(eps, -1), 0, 0);
      return a;
    };
    break;
  case 6:
    A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
      VectorRd a(exp(x(0)), exp(x(1)), exp(x(2)));
      return a;
    };
    break;
  default:
    break;
  }
  return A;
}

CellFType<double> TestCase::div_advec()
{
  CellFType<double> divA = [&](const VectorRd x, const Cell *cell) -> double {
    return 0;
  };
  switch (m_iTC[2])
  {
  case 1:
    break;
  case 2:
    divA = [&](const VectorRd x, const Cell *cell) -> double {
      return 3;
    };
    div_advec_zero = false;
    break;
  case 3:
    break;
  case 4:
    divA = [&](const VectorRd x, const Cell *cell) -> double {
      return 2 * x(0) + 2 * x(1) + 2 * x(2);
    };
    div_advec_zero = false;
    div_advec_const = false;
    break;
  case 5:
    break;
  case 6:
    divA = [&](const VectorRd x, const Cell *cell) -> double {
      return exp(x(0)) + exp(x(1)) + exp(x(2));
    };
    div_advec_zero = false;
    div_advec_const = false;
    break;
  default:
    break;
  }
  return divA;
}

CellFType<double> TestCase::reac()
{
  // Default must be zero
  CellFType<double> R = [&](const VectorRd x, const Cell *cell) -> double {
    return 0;
  };
  switch (m_iTC[3])
  {
  case 1:
    R = [&](const VectorRd x, const Cell *cell) -> double {
      return 1;
    };
    break;
  case 2:
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
    reac_const = false;
    break;
  case 3:
    R = [&](const VectorRd x, const Cell *cell) -> double {
      return x(0) + x(1) + x(2) + eps;
    };
    reac_const = false;
    break;
  case 4:
    R = [&](const VectorRd x, const Cell *cell) -> double {
      return exp(x(0) * x(1) * x(2));
    };
    reac_const = false;
    break;
  case 5:
  {
    R = [&](const VectorRd x, const Cell *cell) -> double {
      return sin(pi * x(0)) * sin(pi * x(1)) * sin(pi * x(2)) + eps;
    };
    reac_const = false;
    break;
  }
  case 6:
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
    reac_const = false;
    break;
  case 7:
    R = [&](const VectorRd x, const Cell *cell) -> double {
      return eps;
    };
    break;
  case 8:
    R = [&](const VectorRd x, const Cell *cell) -> double {
      return pow(eps, -1);
    };
    break;
  default:
    break;
  }
  return R;
}

///////////////////////////// SOURCE TERM ///////////////////////////

// Source term
CellFType<double> TestCase::diff_source()
{

  return [this](const VectorRd p, const Cell *T) -> double {
    Eigen::Matrix3d AHu = diff()(p, T) * hess_sol()(p, T);
    return -AHu.trace() - div_diff()(p, T).dot(grad_sol()(p, T));
  };
}

CellFType<double> TestCase::diff_advec_reac_source()
{
  CellFType<double> f = [&](const VectorRd x, const Cell *cell) -> double {
    return diff_source()(x, cell) + div_advec()(x, cell) * sol()(x) + advec()(x, cell).dot(grad_sol()(x, cell)) + reac()(x, cell) * sol()(x);
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
