#include "MHDTests.hpp"

using namespace HArDCore3D;

// ----------------------------------------------------------------------------
//                          Implementation
// ----------------------------------------------------------------------------

MHDTests::MHDTests(
    size_t u_id,
    size_t b_id,
    size_t p_id,
    double visc,
    double diff) : m_u_id(u_id), m_b_id(b_id), m_p_id(p_id), m_visc(visc), m_diff(diff)
{
    validate();
}

FType<VectorRd> MHDTests::velocity_source()
{
    FType<VectorRd> f = [&](const VectorRd x) -> VectorRd
    {
        return -m_visc * laplace_velocity()(x) + u_dot_grad_u()(x) - b_dot_grad_b()(x) + grad_pressure()(x);
        // return -m_visc * laplace_velocity()(x) + u_dot_grad_u()(x) - b_dot_grad_b()(x);
        // return -m_visc * laplace_velocity()(x) + u_dot_grad_u()(x);
        // return -m_visc * laplace_velocity()(x);
        // return -m_visc * laplace_velocity()(x) + u_dot_grad_u()(x) + grad_pressure()(x);
        // return -m_visc * laplace_velocity()(x) + grad_pressure()(x);
    };
    return f;
}

FType<VectorRd> MHDTests::magnetic_source()
{
    FType<VectorRd> f = [&](const VectorRd x) -> VectorRd
    {
        return -m_diff * laplace_magnetic()(x) + u_dot_grad_b()(x) - b_dot_grad_u()(x);
        // return -m_diff * laplace_magnetic()(x);
    };
    return f;
}

FType<double> MHDTests::pressure()
{
    FType<double> p;
    switch (m_p_id)
    {
    case 0:
        p = [&](const VectorRd x) -> double
        {
            return 0.0;
        };
        break;
    case 1:
        p = [&](const VectorRd x) -> double
        {
            return x(0) + x(1) + x(2) - 3.0 / 2.0;
        };
        break;
    case 2:
        p = [&](const VectorRd x) -> double
        {
            return std::sin(2.0 * PI * x(0)) * std::sin(2.0 * PI * x(1)) * std::sin(2.0 * PI * x(2));
        };
        break;
    }
    return p;
}

FType<VectorRd> MHDTests::grad_pressure()
{
    FType<VectorRd> Dp;
    switch (m_p_id)
    {
    case 0:
        Dp = [&](const VectorRd x) -> VectorRd
        {
            return VectorRd(0.0, 0.0, 0.0);
        };
        break;
    case 1:
        Dp = [&](const VectorRd x) -> VectorRd
        {
            return VectorRd(1.0, 1.0, 1.0);
        };
        break;
    case 2:
        Dp = [&](const VectorRd x) -> VectorRd
        {
            return 2.0 * PI * VectorRd(std::cos(2.0 * PI * x(0)) * std::sin(2.0 * PI * x(1)) * std::sin(2.0 * PI * x(2)), std::sin(2.0 * PI * x(0)) * std::cos(2.0 * PI * x(1)) * std::sin(2.0 * PI * x(2)), std::sin(2.0 * PI * x(0)) * std::sin(2.0 * PI * x(1)) * std::cos(2.0 * PI * x(2)));
        };
        break;
    }
    return Dp;
}

FType<VectorRd> MHDTests::velocity()
{
    FType<VectorRd> u;
    switch (m_u_id)
    {
    case 0:
        u = [&](const VectorRd x) -> VectorRd
        {
            return VectorRd(0.0, 0.0, 0.0);
        };
        break;
    case 1:
        u = [&](const VectorRd x) -> VectorRd
        {
            return std::sin(PI * x(0)) * std::sin(PI * x(1)) * std::sin(PI * x(2)) * VectorRd(std::sin(PI * x(0)) * std::sin(PI * (x(1) - x(2))), std::sin(PI * x(1)) * std::sin(PI * (x(2) - x(0))), std::sin(PI * x(2)) * std::sin(PI * (x(0) - x(1))));
        };
        break;
    case 2:
        u = [&](const VectorRd x) -> VectorRd
        {
            return std::pow(std::sin(PI * x(0)) * std::sin(PI * x(1)) * std::sin(PI * x(2)), 2) * VectorRd(std::sin(PI * x(0)) * std::sin(PI * (x(1) - x(2))), std::sin(PI * x(1)) * std::sin(PI * (x(2) - x(0))), std::sin(PI * x(2)) * std::sin(PI * (x(0) - x(1))));
        };
        break;
    case 3:
        u = [&](const VectorRd x) -> VectorRd
        {
            return VectorRd(x(0) * (1 - x(0)), x(1) * (1 - x(1)), x(2) * (1 - x(2))); // not zero div!! Just for debugging
        };
        break;
    }
    return u;
}

FType<MatrixRd> MHDTests::grad_velocity()
{
    FType<MatrixRd> Du;
    switch (m_u_id)
    {
    case 0:
        Du = [&](const VectorRd x) -> MatrixRd
        {
            return MatrixRd::Zero();
        };
        break;
    case 1:
        Du = [&](const VectorRd x) -> MatrixRd
        {
            double d1u1 = std::sin(2.0 * PI * x(0)) * std::sin(PI * (x(1) - x(2)));
            double d2u1 = std::sin(PI * x(0)) * std::sin(PI * (2.0 * x(1) - x(2)));
            double d3u1 = std::sin(PI * x(0)) * std::sin(PI * (x(1) - 2.0 * x(2)));
            double d1u2 = std::sin(PI * x(1)) * std::sin(PI * (x(2) - 2.0 * x(0)));
            double d2u2 = std::sin(2.0 * PI * x(1)) * std::sin(PI * (x(2) - x(0)));
            double d3u2 = std::sin(PI * x(1)) * std::sin(PI * (2.0 * x(2) - x(0)));
            double d1u3 = std::sin(PI * x(2)) * std::sin(PI * (2.0 * x(0) - x(1)));
            double d2u3 = std::sin(PI * x(2)) * std::sin(PI * (x(0) - 2.0 * x(1)));
            double d3u3 = std::sin(2.0 * PI * x(2)) * std::sin(PI * (x(0) - x(1)));

            MatrixRd m = MatrixRd::Zero();
            m.col(0) = PI * std::sin(PI * x(1)) * std::sin(PI * x(2)) * VectorRd(d1u1, d1u2, d1u3);
            m.col(1) = PI * std::sin(PI * x(0)) * std::sin(PI * x(2)) * VectorRd(d2u1, d2u2, d2u3);
            m.col(2) = PI * std::sin(PI * x(0)) * std::sin(PI * x(1)) * VectorRd(d3u1, d3u2, d3u3);
            return m;
        };
        break;
    case 2:
        Du = [&](const VectorRd x) -> MatrixRd
        {
            double d1u1 = 1.5 * std::sin(2.0 * PI * x(0)) * std::sin(PI * (x(1) - x(2)));
            double d2u1 = -0.5 * std::sin(PI * x(0)) * (-3.0 * std::sin(PI * (2 * x(1) - x(2))) + std::sin(x(2)));
            double d3u1 = 0.5 * std::sin(PI * x(0)) * (-3.0 * std::sin(PI * (2 * x(2) - x(1))) + std::sin(x(1)));
            double d1u2 = 0.5 * std::sin(PI * x(1)) * (-3.0 * std::sin(PI * (2 * x(0) - x(2))) + std::sin(x(2)));
            double d2u2 = 1.5 * std::sin(2.0 * PI * x(1)) * std::sin(PI * (x(2) - x(0)));
            double d3u2 = -0.5 * std::sin(PI * x(1)) * (-3.0 * std::sin(PI * (2.0 * x(2) - x(0))) + std::sin(PI * x(0)));
            double d1u3 = -0.5 * std::sin(PI * x(2)) * (-3.0 * std::sin(PI * (2 * x(0) - x(1))) + std::sin(x(1)));
            double d2u3 = 0.5 * std::sin(PI * x(0)) * (-3.0 * std::sin(PI * (2 * x(1) - x(0))) + std::sin(x(0)));
            double d3u3 = 1.5 * std::sin(2.0 * PI * x(2)) * std::sin(PI * (x(0) - x(1)));

            MatrixRd m = MatrixRd::Zero();
            m.col(0) = PI * std::sin(PI * x(0)) * std::pow(std::sin(PI * x(1)) * std::sin(PI * x(2)), 2) * VectorRd(d1u1, d1u2, d1u3);
            m.col(1) = PI * std::sin(PI * x(1)) * std::pow(std::sin(PI * x(0)) * std::sin(PI * x(2)), 2) * VectorRd(d2u1, d2u2, d2u3);
            m.col(2) = PI * std::sin(PI * x(2)) * std::pow(std::sin(PI * x(0)) * std::sin(PI * x(1)), 2) * VectorRd(d3u1, d3u2, d3u3);
            return m;
        };
        break;
    case 3:
        Du = [&](const VectorRd x) -> MatrixRd
        {
            double d1u1 = 1.0 - 2.0 * x(0);
            double d2u1 = 0.0;
            double d3u1 = 0.0;
            double d1u2 = 0.0;
            double d2u2 = 1.0 - 2.0 * x(1);
            double d3u2 = 0.0;
            double d1u3 = 0.0;
            double d2u3 = 0.0;
            double d3u3 = 1.0 - 2.0 * x(2);

            MatrixRd m = MatrixRd::Zero();
            m.col(0) = VectorRd(d1u1, d1u2, d1u3);
            m.col(1) = VectorRd(d2u1, d2u2, d2u3);
            m.col(2) = VectorRd(d3u1, d3u2, d3u3);
            return m;
        };
        break;
    }
    return Du;
}

FType<VectorRd> MHDTests::magnetic()
{
    FType<VectorRd> b;
    switch (m_b_id)
    {
    case 0:
        b = [&](const VectorRd x) -> VectorRd
        {
            return VectorRd::Zero();
        };
        break;
    case 1:
        b = [&](const VectorRd x) -> VectorRd
        {
            return std::sin(PI * x(0)) * std::sin(PI * x(1)) * std::sin(PI * x(2)) * VectorRd(std::sin(PI * x(0)) * std::sin(PI * (x(1) - x(2))), std::sin(PI * x(1)) * std::sin(PI * (x(2) - x(0))), std::sin(PI * x(2)) * std::sin(PI * (x(0) - x(1))));
        };
        break;
    case 2:
        b = [&](const VectorRd x) -> VectorRd
        {
            return std::pow(std::sin(PI * x(0)) * std::sin(PI * x(1)) * std::sin(PI * x(2)), 2) * VectorRd(std::sin(PI * x(0)) * std::sin(PI * (x(1) - x(2))), std::sin(PI * x(1)) * std::sin(PI * (x(2) - x(0))), std::sin(PI * x(2)) * std::sin(PI * (x(0) - x(1))));
        };
        break;
    case 3:
        b = [&](const VectorRd x) -> VectorRd
        {
            double b1 = -0.5 * std::sin(PI * x(0)) * std::cos(PI * x(1)) * std::cos(PI * x(2));
            double b2 = std::cos(PI * x(0)) * std::sin(PI * x(1)) * std::cos(PI * x(2));
            double b3 = -0.5 * std::cos(PI * x(0)) * std::cos(PI * x(1)) * std::sin(PI * x(2));
            return VectorRd(b1, b2, b3);
        };
        break;
    }
    return b;
}

FType<MatrixRd> MHDTests::grad_magnetic()
{
    FType<MatrixRd> Db;
    switch (m_b_id)
    {
    case 0:
        Db = [&](const VectorRd x) -> MatrixRd
        {
            return MatrixRd::Zero();
        };
        break;
    case 1:
        Db = [&](const VectorRd x) -> MatrixRd
        {
            double d1b1 = std::sin(2.0 * PI * x(0)) * std::sin(PI * (x(1) - x(2)));
            double d2b1 = std::sin(PI * x(0)) * std::sin(PI * (2.0 * x(1) - x(2)));
            double d3b1 = std::sin(PI * x(0)) * std::sin(PI * (x(1) - 2.0 * x(2)));
            double d1b2 = std::sin(PI * x(1)) * std::sin(PI * (x(2) - 2.0 * x(0)));
            double d2b2 = std::sin(2.0 * PI * x(1)) * std::sin(PI * (x(2) - x(0)));
            double d3b2 = std::sin(PI * x(1)) * std::sin(PI * (2.0 * x(2) - x(0)));
            double d1b3 = std::sin(PI * x(2)) * std::sin(PI * (2.0 * x(0) - x(1)));
            double d2b3 = std::sin(PI * x(2)) * std::sin(PI * (x(0) - 2.0 * x(1)));
            double d3b3 = std::sin(2.0 * PI * x(2)) * std::sin(PI * (x(0) - x(1)));

            MatrixRd m = MatrixRd::Zero();
            m.col(0) = PI * std::sin(PI * x(1)) * std::sin(PI * x(2)) * VectorRd(d1b1, d1b2, d1b3);
            m.col(1) = PI * std::sin(PI * x(0)) * std::sin(PI * x(2)) * VectorRd(d2b1, d2b2, d2b3);
            m.col(2) = PI * std::sin(PI * x(0)) * std::sin(PI * x(1)) * VectorRd(d3b1, d3b2, d3b3);
            return m;
        };
        break;
    case 2:
        Db = [&](const VectorRd x) -> MatrixRd
        {
            double d1b1 = 1.5 * std::sin(2.0 * PI * x(0)) * std::sin(PI * (x(1) - x(2)));
            double d2b1 = -0.5 * std::sin(PI * x(0)) * (-3.0 * std::sin(PI * (2 * x(1) - x(2))) + std::sin(x(2)));
            double d3b1 = 0.5 * std::sin(PI * x(0)) * (-3.0 * std::sin(PI * (2 * x(2) - x(1))) + std::sin(x(1)));
            double d1b2 = 0.5 * std::sin(PI * x(1)) * (-3.0 * std::sin(PI * (2 * x(0) - x(2))) + std::sin(x(2)));
            double d2b2 = 1.5 * std::sin(2.0 * PI * x(1)) * std::sin(PI * (x(2) - x(0)));
            double d3b2 = -0.5 * std::sin(PI * x(1)) * (-3.0 * std::sin(PI * (2.0 * x(2) - x(0))) + std::sin(PI * x(0)));
            double d1b3 = -0.5 * std::sin(PI * x(2)) * (-3.0 * std::sin(PI * (2 * x(0) - x(1))) + std::sin(x(1)));
            double d2b3 = 0.5 * std::sin(PI * x(0)) * (-3.0 * std::sin(PI * (2 * x(1) - x(0))) + std::sin(x(0)));
            double d3b3 = 1.5 * std::sin(2.0 * PI * x(2)) * std::sin(PI * (x(0) - x(1)));

            MatrixRd m = MatrixRd::Zero();
            m.col(0) = PI * std::sin(PI * x(0)) * std::pow(std::sin(PI * x(1)) * std::sin(PI * x(2)), 2) * VectorRd(d1b1, d1b2, d1b3);
            m.col(1) = PI * std::sin(PI * x(1)) * std::pow(std::sin(PI * x(0)) * std::sin(PI * x(2)), 2) * VectorRd(d2b1, d2b2, d2b3);
            m.col(2) = PI * std::sin(PI * x(2)) * std::pow(std::sin(PI * x(0)) * std::sin(PI * x(1)), 2) * VectorRd(d3b1, d3b2, d3b3);
            return m;
        };
        break;
    case 3:
        Db = [&](const VectorRd x) -> MatrixRd
        {
            double d1b1 = -0.5 * PI * std::cos(PI * x(0)) * std::cos(PI * x(1)) * std::cos(PI * x(2));
            double d2b1 = 0.5 * PI * std::sin(PI * x(0)) * std::sin(PI * x(1)) * std::cos(PI * x(2));
            double d3b1 = 0.5 * PI * std::sin(PI * x(0)) * std::cos(PI * x(1)) * std::sin(PI * x(2));

            double d1b2 = -PI * std::sin(PI * x(0)) * std::sin(PI * x(1)) * std::cos(PI * x(2));
            double d2b2 = PI * std::cos(PI * x(0)) * std::cos(PI * x(1)) * std::cos(PI * x(2));
            double d3b2 = -PI * std::cos(PI * x(0)) * std::sin(PI * x(1)) * std::sin(PI * x(2));

            double d1b3 = 0.5 * PI * std::sin(PI * x(0)) * std::cos(PI * x(1)) * std::sin(PI * x(2));
            double d2b3 = 0.5 * PI * std::cos(PI * x(0)) * std::sin(PI * x(1)) * std::sin(PI * x(2));
            double d3b3 = -0.5 * PI * std::cos(PI * x(0)) * std::cos(PI * x(1)) * std::cos(PI * x(2));

            MatrixRd m = MatrixRd::Zero();
            m.col(0) = VectorRd(d1b1, d1b2, d1b3);
            m.col(1) = VectorRd(d2b1, d2b2, d2b3);
            m.col(2) = VectorRd(d3b1, d3b2, d3b3);
            return m;
        };
        break;
    }
    return Db;
}

FType<VectorRd> MHDTests::u_dot_grad_u()
{
    if (m_u_id == 3)
    {
        FType<VectorRd> uDu_temp = [&](const VectorRd x) -> VectorRd
        {
            return VectorRd(x(0) * (1 - x(0)) * (1 - 2.0 * x(0)), x(1) * (1 - x(1)) * (1 - 2.0 * x(1)), x(2) * (1 - x(2)) * (1 - 2.0 * x(2)));
        };
    }

    FType<VectorRd> uDu = [&](const VectorRd x) -> VectorRd
    {
        MatrixRd grad_u = this->grad_velocity()(x);

        double d1u1 = grad_u(0, 0);
        double d2u1 = grad_u(0, 1);
        double d3u1 = grad_u(0, 2);

        double d1u2 = grad_u(1, 0);
        double d2u2 = grad_u(1, 1);
        double d3u2 = grad_u(1, 2);

        double d1u3 = grad_u(2, 0);
        double d2u3 = grad_u(2, 1);
        double d3u3 = grad_u(2, 2);

        double u1 = this->velocity()(x)(0);
        double u2 = this->velocity()(x)(1);
        double u3 = this->velocity()(x)(2);

        return VectorRd(u1 * d1u1 + u2 * d2u1 + u3 * d3u1, u1 * d1u2 + u2 * d2u2 + u3 * d3u2, u1 * d1u3 + u2 * d2u3 + u3 * d3u3);
    };

    return uDu;
}

FType<VectorRd> MHDTests::u_dot_grad_b()
{
    FType<VectorRd> uDb = [&](const VectorRd x) -> VectorRd
    {
        MatrixRd grad_b = this->grad_magnetic()(x);

        double d1b1 = grad_b(0, 0);
        double d2b1 = grad_b(0, 1);
        double d3b1 = grad_b(0, 2);

        double d1b2 = grad_b(1, 0);
        double d2b2 = grad_b(1, 1);
        double d3b2 = grad_b(1, 2);

        double d1b3 = grad_b(2, 0);
        double d2b3 = grad_b(2, 1);
        double d3b3 = grad_b(2, 2);

        double u1 = this->velocity()(x)(0);
        double u2 = this->velocity()(x)(1);
        double u3 = this->velocity()(x)(2);

        return VectorRd(u1 * d1b1 + u2 * d2b1 + u3 * d3b1, u1 * d1b2 + u2 * d2b2 + u3 * d3b2, u1 * d1b3 + u2 * d2b3 + u3 * d3b3);
    };

    return uDb;
}

FType<VectorRd> MHDTests::b_dot_grad_u()
{
    FType<VectorRd> bDu = [&](const VectorRd x) -> VectorRd
    {
        MatrixRd grad_u = this->grad_velocity()(x);

        double d1u1 = grad_u(0, 0);
        double d2u1 = grad_u(0, 1);
        double d3u1 = grad_u(0, 2);

        double d1u2 = grad_u(1, 0);
        double d2u2 = grad_u(1, 1);
        double d3u2 = grad_u(1, 2);

        double d1u3 = grad_u(2, 0);
        double d2u3 = grad_u(2, 1);
        double d3u3 = grad_u(2, 2);

        double b1 = this->magnetic()(x)(0);
        double b2 = this->magnetic()(x)(1);
        double b3 = this->magnetic()(x)(2);

        return VectorRd(b1 * d1u1 + b2 * d2u1 + b3 * d3u1, b1 * d1u2 + b2 * d2u2 + b3 * d3u2, b1 * d1u3 + b2 * d2u3 + b3 * d3u3);
    };

    return bDu;
}

FType<VectorRd> MHDTests::b_dot_grad_b()
{
    FType<VectorRd> bDb = [&](const VectorRd x) -> VectorRd
    {
        MatrixRd grad_b = this->grad_magnetic()(x);

        double d1b1 = grad_b(0, 0);
        double d2b1 = grad_b(0, 1);
        double d3b1 = grad_b(0, 2);

        double d1b2 = grad_b(1, 0);
        double d2b2 = grad_b(1, 1);
        double d3b2 = grad_b(1, 2);

        double d1b3 = grad_b(2, 0);
        double d2b3 = grad_b(2, 1);
        double d3b3 = grad_b(2, 2);

        double b1 = this->magnetic()(x)(0);
        double b2 = this->magnetic()(x)(1);
        double b3 = this->magnetic()(x)(2);

        return VectorRd(b1 * d1b1 + b2 * d2b1 + b3 * d3b1, b1 * d1b2 + b2 * d2b2 + b3 * d3b2, b1 * d1b3 + b2 * d2b3 + b3 * d3b3);
    };

    return bDb;
}

FType<VectorRd> MHDTests::laplace_velocity()
{
    FType<VectorRd> Lu;
    switch (m_u_id)
    {
    case 0:
        Lu = [&](const VectorRd x) -> VectorRd
        {
            return VectorRd::Zero();
        };
        break;
    case 1:
        Lu = [&](const VectorRd x) -> VectorRd
        {
            double l1 = (2.0 * std::cos(2.0 * PI * x(0)) - 1.0) * std::sin(PI * x(1)) * std::sin(PI * (x(1) - x(2))) * std::sin(PI * x(2)) - 0.5 * std::sin(PI * x(0)) * std::sin(PI * x(0)) * std::sin(2.0 * PI * (x(1) - x(2)));
            double l2 = (2.0 * std::cos(2.0 * PI * x(1)) - 1.0) * std::sin(PI * x(0)) * std::sin(PI * (x(2) - x(0))) * std::sin(PI * x(2)) - 0.5 * std::sin(PI * x(1)) * std::sin(PI * x(1)) * std::sin(2.0 * PI * (x(2) - x(0)));
            double l3 = (2.0 * std::cos(2.0 * PI * x(2)) - 1.0) * std::sin(PI * x(0)) * std::sin(PI * (x(0) - x(1))) * std::sin(PI * x(1)) - 0.5 * std::sin(PI * x(2)) * std::sin(PI * x(2)) * std::sin(2.0 * PI * (x(0) - x(1)));
            return 2.0 * PI * PI * VectorRd(l1, l2, l3);
        };
        break;
    case 2:
        Lu = [&](const VectorRd x) -> VectorRd
        {
            // double SinPix = std::sin(PI * x(0));
            // double SinPiy = std::sin(PI * x(1));
            // double SinPiz = std::sin(PI * x(2));

            double SinTwoPix = std::sin(2.0 * PI * x(0));
            double SinTwoPiy = std::sin(2.0 * PI * x(1));
            double SinTwoPiz = std::sin(2.0 * PI * x(2));

            double SinPixSquared = std::pow(std::sin(PI * x(0)), 2);
            double SinPiySquared = std::pow(std::sin(PI * x(1)), 2);
            double SinPizSquared = std::pow(std::sin(PI * x(2)), 2);

            double CosPixSquared = std::pow(std::cos(PI * x(0)), 2);
            double CosPiySquared = std::pow(std::cos(PI * x(1)), 2);
            double CosPizSquared = std::pow(std::cos(PI * x(2)), 2);

            double l1 = SinPixSquared * (2.0 * SinPiySquared * CosPizSquared + 2.0 * CosPiySquared * SinPizSquared - 13.0 * SinPiySquared * SinPizSquared - SinTwoPiy * SinTwoPiz) + 6.0 * CosPixSquared * SinPiySquared * SinPizSquared;

            double l2 = SinPiySquared * (2.0 * SinPixSquared * CosPizSquared + 2.0 * CosPixSquared * SinPizSquared - 13.0 * SinPixSquared * SinPizSquared - SinTwoPix * SinTwoPiz) + 6.0 * CosPiySquared * SinPixSquared * SinPizSquared;

            double l3 = SinPizSquared * (2.0 * SinPixSquared * CosPiySquared + 2.0 * CosPixSquared * SinPiySquared - 13.0 * SinPixSquared * SinPiySquared - SinTwoPix * SinTwoPiy) + 6.0 * CosPizSquared * SinPixSquared * SinPizSquared;

            l1 *= PI * PI * std::sin(PI * x(0)) * std::sin(PI * (x(1) - x(2)));

            l2 *= PI * PI * std::sin(PI * x(1)) * std::sin(PI * (x(2) - x(0)));

            l3 *= PI * PI * std::sin(PI * x(2)) * std::sin(PI * (x(0) - x(1)));

            return VectorRd(l1, l2, l3);
        };
        break;
    case 3:
        Lu = [&](const VectorRd x) -> VectorRd
        {
            return -2.0 * VectorRd(1.0, 1.0, 1.0);
        };
        break;
    }
    return Lu;
}

FType<VectorRd> MHDTests::laplace_magnetic()
{
    FType<VectorRd> Lb;
    switch (m_b_id)
    {
    case 0:
        Lb = [&](const VectorRd x) -> VectorRd
        {
            return VectorRd::Zero();
        };
        break;
    case 1:
        Lb = [&](const VectorRd x) -> VectorRd
        {
            double l1 = (2.0 * std::cos(2.0 * PI * x(0)) - 1.0) * std::sin(PI * x(1)) * std::sin(PI * (x(1) - x(2))) * std::sin(PI * x(2)) - 0.5 * std::sin(PI * x(0)) * std::sin(PI * x(0)) * std::sin(2.0 * PI * (x(1) - x(2)));
            double l2 = (2.0 * std::cos(2.0 * PI * x(1)) - 1.0) * std::sin(PI * x(0)) * std::sin(PI * (x(2) - x(0))) * std::sin(PI * x(2)) - 0.5 * std::sin(PI * x(1)) * std::sin(PI * x(1)) * std::sin(2.0 * PI * (x(2) - x(0)));
            double l3 = (2.0 * std::cos(2.0 * PI * x(2)) - 1.0) * std::sin(PI * x(0)) * std::sin(PI * (x(0) - x(1))) * std::sin(PI * x(1)) - 0.5 * std::sin(PI * x(2)) * std::sin(PI * x(2)) * std::sin(2.0 * PI * (x(0) - x(1)));
            return 2.0 * PI * PI * VectorRd(l1, l2, l3);
        };
        break;
    case 2:
        Lb = [&](const VectorRd x) -> VectorRd
        {
            // double SinPix = std::sin(PI * x(0));
            // double SinPiy = std::sin(PI * x(1));
            // double SinPiz = std::sin(PI * x(2));

            double SinTwoPix = std::sin(2.0 * PI * x(0));
            double SinTwoPiy = std::sin(2.0 * PI * x(1));
            double SinTwoPiz = std::sin(2.0 * PI * x(2));

            double SinPixSquared = std::pow(std::sin(PI * x(0)), 2);
            double SinPiySquared = std::pow(std::sin(PI * x(1)), 2);
            double SinPizSquared = std::pow(std::sin(PI * x(2)), 2);

            double CosPixSquared = std::pow(std::cos(PI * x(0)), 2);
            double CosPiySquared = std::pow(std::cos(PI * x(1)), 2);
            double CosPizSquared = std::pow(std::cos(PI * x(2)), 2);

            double l1 = SinPixSquared * (2.0 * SinPiySquared * CosPizSquared + 2.0 * CosPiySquared * SinPizSquared - 13.0 * SinPiySquared * SinPizSquared - SinTwoPiy * SinTwoPiz) + 6.0 * CosPixSquared * SinPiySquared * SinPizSquared;

            double l2 = SinPiySquared * (2.0 * SinPixSquared * CosPizSquared + 2.0 * CosPixSquared * SinPizSquared - 13.0 * SinPixSquared * SinPizSquared - SinTwoPix * SinTwoPiz) + 6.0 * CosPiySquared * SinPixSquared * SinPizSquared;

            double l3 = SinPizSquared * (2.0 * SinPixSquared * CosPiySquared + 2.0 * CosPixSquared * SinPiySquared - 13.0 * SinPixSquared * SinPiySquared - SinTwoPix * SinTwoPiy) + 6.0 * CosPizSquared * SinPixSquared * SinPizSquared;

            l1 *= PI * PI * std::sin(PI * x(0)) * std::sin(PI * (x(1) - x(2)));

            l2 *= PI * PI * std::sin(PI * x(1)) * std::sin(PI * (x(2) - x(0)));

            l3 *= PI * PI * std::sin(PI * x(2)) * std::sin(PI * (x(0) - x(1)));

            return VectorRd(l1, l2, l3);
        };
        break;
    case 3:
        Lb = [&](const VectorRd x) -> VectorRd
        {
            double l1 = 1.5 * PI * PI * std::sin(PI * x(0)) * std::cos(PI * x(1)) * std::cos(PI * x(2));
            double l2 = -3.0 * PI * PI * std::cos(PI * x(0)) * std::sin(PI * x(1)) * std::cos(PI * x(2));
            double l3 = 1.5 * PI * PI * std::cos(PI * x(0)) * std::cos(PI * x(1)) * std::sin(PI * x(2));

            return VectorRd(l1, l2, l3);
        };
        break;
    }
    return Lb;
}

void MHDTests::validate()
{
    assert(m_p_id == 0 || m_p_id == 1 || m_p_id == 2);
    assert(m_u_id == 0 || m_u_id == 1 || m_u_id == 2 || m_u_id == 3);
    assert(m_b_id == 0 || m_b_id == 1 || m_b_id == 2 || m_b_id == 3);
    // assert(m_b_id == 0 || m_b_id == 2); // assuming homogeneous Hodge conditions - which 1 does not support
    assert(m_visc > 0);
    assert(m_diff > 0);
}