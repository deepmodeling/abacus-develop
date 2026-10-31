#include "snap_descriptor.h"

#include "source_base/constants.h"

#include <cmath>
#include <stdexcept>

namespace ModuleMD
{
namespace SNAP
{

namespace
{
struct Dual
{
    double value;
    double derivative;

    Dual(const double value_in = 0.0, const double derivative_in = 0.0)
        : value(value_in), derivative(derivative_in)
    {
    }
};

Dual operator+(const Dual& left, const Dual& right)
{
    return Dual(left.value + right.value, left.derivative + right.derivative);
}

Dual operator-(const Dual& left, const Dual& right)
{
    return Dual(left.value - right.value, left.derivative - right.derivative);
}

Dual operator-(const Dual& value)
{
    return Dual(-value.value, -value.derivative);
}

Dual operator*(const Dual& left, const Dual& right)
{
    return Dual(left.value * right.value,
                left.derivative * right.value + left.value * right.derivative);
}

Dual operator/(const Dual& left, const Dual& right)
{
    const double inverse = 1.0 / right.value;
    return Dual(left.value * inverse,
                (left.derivative * right.value - left.value * right.derivative) * inverse * inverse);
}

Dual operator+(const Dual& left, const double right)
{
    return left + Dual(right);
}

Dual operator-(const Dual& left, const double right)
{
    return left - Dual(right);
}

Dual operator*(const Dual& left, const double right)
{
    return left * Dual(right);
}

Dual operator*(const double left, const Dual& right)
{
    return Dual(left) * right;
}

Dual operator/(const Dual& left, const double right)
{
    return left / Dual(right);
}

Dual operator/(const double left, const Dual& right)
{
    return Dual(left) / right;
}

Dual& operator*=(Dual& left, const Dual& right)
{
    left = left * right;
    return left;
}

Dual& operator+=(Dual& left, const Dual& right)
{
    left = left + right;
    return left;
}

Dual dual_sqrt(const Dual& value)
{
    const double root = std::sqrt(value.value);
    return Dual(root, value.derivative / (2.0 * root));
}

Dual dual_cos(const Dual& value)
{
    return Dual(std::cos(value.value), -std::sin(value.value) * value.derivative);
}

Dual dual_tan(const Dual& value)
{
    const double cosine = std::cos(value.value);
    return Dual(std::tan(value.value), value.derivative / (cosine * cosine));
}
}

UWithDerivative DescriptorEvaluator::compute_u_with_derivative(const Neighbor& neighbor, const int direction) const
{
    if (direction < 0 || direction > 2)
    {
        throw std::invalid_argument("SNAP derivative direction must be 0, 1, or 2.");
    }
    const Dual x(neighbor.x, direction == 0 ? 1.0 : 0.0);
    const Dual y(neighbor.y, direction == 1 ? 1.0 : 0.0);
    const Dual z(neighbor.z, direction == 2 ? 1.0 : 0.0);
    const Dual rsq = x * x + y * y + z * z;
    const Dual distance = dual_sqrt(rsq);
    const double cutoff_radius = neighbor.cutoff_radius > 0.0 ? neighbor.cutoff_radius : parameters_.rcut;

    UWithDerivative result;
    result.values.assign(u_size_, std::complex<double>(0.0, 0.0));
    result.derivatives.assign(u_size_, std::complex<double>(0.0, 0.0));
    if (distance.value <= 1.0e-14 || distance.value >= cutoff_radius)
    {
        return result;
    }

    Dual scale;
    if (!parameters_.switch_flag)
    {
        scale = Dual(1.0);
    }
    else if (distance.value <= parameters_.rmin)
    {
        scale = Dual(1.0);
    }
    else
    {
        const Dual angle = (distance - parameters_.rmin) * ModuleBase::PI
                           / (cutoff_radius - parameters_.rmin);
        scale = 0.5 * (dual_cos(angle) + 1.0);
    }
    scale = scale * neighbor.weight;

    const Dual theta0 = (distance - parameters_.rmin) * parameters_.rfac0 * ModuleBase::PI
                        / (cutoff_radius - parameters_.rmin);
    const Dual z0 = distance / dual_tan(theta0);
    const Dual r0inv = 1.0 / dual_sqrt(rsq + z0 * z0);
    const Dual ar = r0inv * z0;
    const Dual ai = -r0inv * z;
    const Dual br = r0inv * y;
    const Dual bi = -r0inv * x;
    std::vector<std::complex<Dual> > raw(u_size_, std::complex<Dual>(Dual(0.0), Dual(0.0)));
    raw[0] = std::complex<Dual>(Dual(1.0), Dual(0.0));
    const int dimension = parameters_.twojmax + 1;
    for (int j = 1; j <= parameters_.twojmax; ++j)
    {
        for (int mb = 0; 2 * mb <= j; ++mb)
        {
            raw[u_index(j, 0, mb)] = std::complex<Dual>(Dual(0.0), Dual(0.0));
            for (int ma = 0; ma < j; ++ma)
            {
                const std::size_t current = u_index(j, ma, mb);
                const std::size_t previous = u_index(j - 1, ma, mb);
                const double first_root = rootpq_[static_cast<std::size_t>((j - ma) * dimension + (j - mb))];
                const double second_root = rootpq_[static_cast<std::size_t>((ma + 1) * dimension + (j - mb))];
                raw[current] += Dual(first_root) * std::complex<Dual>(ar, -ai) * raw[previous];
                raw[current + 1] = Dual(-second_root) * std::complex<Dual>(br, -bi) * raw[previous];
            }
        }
        for (int mb = 0; 2 * mb <= j; ++mb)
        {
            for (int ma = 0; ma <= j; ++ma)
            {
                const std::size_t left = u_index(j, ma, mb);
                const std::size_t right = u_index(j, j - ma, j - mb);
                const double sign = ((ma - mb) % 2 == 0) ? 1.0 : -1.0;
                raw[right] = Dual(sign) * std::conj(raw[left]);
            }
        }
    }
    for (std::size_t i = 0; i < raw.size(); ++i)
    {
        const std::complex<Dual> value = scale * raw[i];
        result.values[i] = std::complex<double>(value.real().value, value.imag().value);
        result.derivatives[i] = std::complex<double>(value.real().derivative, value.imag().derivative);
    }
    return result;
}

} // namespace SNAP
} // namespace ModuleMD
