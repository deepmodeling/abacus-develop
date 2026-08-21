#ifndef ABACUS_D3_TYPES_H
#define ABACUS_D3_TYPES_H

#include <array>
#include <vector>

namespace vdw
{
namespace d3
{

struct Vec3
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vec3() = default;
    Vec3(double x_in, double y_in, double z_in) : x(x_in), y(y_in), z(z_in) {}

    Vec3& operator+=(const Vec3& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    Vec3& operator-=(const Vec3& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }
};

inline Vec3 operator+(Vec3 lhs, const Vec3& rhs)
{
    lhs += rhs;
    return lhs;
}

inline Vec3 operator-(Vec3 lhs, const Vec3& rhs)
{
    lhs -= rhs;
    return lhs;
}

inline Vec3 operator-(const Vec3& value)
{
    return Vec3(-value.x, -value.y, -value.z);
}

inline Vec3 operator*(const Vec3& value, double factor)
{
    return Vec3(value.x * factor, value.y * factor, value.z * factor);
}

inline Vec3 operator*(double factor, const Vec3& value)
{
    return value * factor;
}

inline Vec3 operator/(const Vec3& value, double divisor)
{
    return value * (1.0 / divisor);
}

struct Matrix3
{
    double value[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
};

struct Structure
{
    std::vector<int> atomic_numbers;
    std::vector<Vec3> positions;
    std::array<Vec3, 3> lattice;
    std::array<bool, 3> periodic = {{false, false, false}};
};

enum class Damping
{
    Zero,
    Rational
};

struct Parameters
{
    Damping damping = Damping::Rational;
    double s6 = 1.0;
    double s8 = 1.0;
    double s9 = 0.0;
    double rs6 = 1.0;
    double rs8 = 1.0;
    double a1 = 0.4;
    double a2 = 5.0;
    double alp = 14.0;
};

struct Cutoffs
{
    double disp2 = 60.0;
    double disp3 = 40.0;
    double cn = 40.0;
    double width2 = 0.0;
    double width3 = 0.0;
};

struct Result
{
    double energy = 0.0;
    std::vector<Vec3> gradient;
    Matrix3 virial;
};

} // namespace d3
} // namespace vdw

#endif // ABACUS_D3_TYPES_H
