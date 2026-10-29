#include "snap_descriptor.h"

#include "source_base/constants.h"

#include <algorithm>
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

DescriptorEvaluator::DescriptorEvaluator(const ModelParameters& parameters)
    : parameters_(parameters), cg_table_(parameters.twojmax), u_size_(0)
{
    std::string error;
    if (!parameters_.valid(error))
    {
        throw std::invalid_argument(error);
    }

    u_block_.resize(static_cast<std::size_t>(parameters_.twojmax + 1), 0);
    for (int j = 0; j <= parameters_.twojmax; ++j)
    {
        u_block_[static_cast<std::size_t>(j)] = u_size_;
        u_size_ += static_cast<std::size_t>((j + 1) * (j + 1));
    }
    compute_rootpq();

    for (int j1 = 0; j1 <= parameters_.twojmax; ++j1)
    {
        for (int j2 = 0; j2 <= j1; ++j2)
        {
            for (int j = j1 - j2; j <= std::min(parameters_.twojmax, j1 + j2); j += 2)
            {
                if (j >= j1)
                {
                    BispectrumIndex index;
                    index.j1 = j1;
                    index.j2 = j2;
                    index.j = j;
                    bispectrum_indices_.push_back(index);
                }
            }
        }
    }
}

void DescriptorEvaluator::compute_rootpq()
{
    const int dimension = parameters_.twojmax + 1;
    rootpq_.assign(static_cast<std::size_t>(dimension * dimension), 0.0);
    for (int p = 1; p <= parameters_.twojmax; ++p)
    {
        for (int q = 1; q <= parameters_.twojmax; ++q)
        {
            rootpq_[static_cast<std::size_t>(p * dimension + q)] = std::sqrt(static_cast<double>(p) / q);
        }
    }
}

std::size_t DescriptorEvaluator::u_index(const int j, const int ma, const int mb) const
{
    return u_block_[static_cast<std::size_t>(j)] + static_cast<std::size_t>(mb * (j + 1) + ma);
}

double DescriptorEvaluator::cutoff(const double distance, const double cutoff_radius) const
{
    if (!parameters_.switch_flag)
    {
        return distance < cutoff_radius ? 1.0 : 0.0;
    }
    if (distance <= parameters_.rmin)
    {
        return 1.0;
    }
    if (distance >= cutoff_radius)
    {
        return 0.0;
    }
    const double angle = ModuleBase::PI * (distance - parameters_.rmin)
                         / (cutoff_radius - parameters_.rmin);
    return 0.5 * (std::cos(angle) + 1.0);
}

std::vector<std::complex<double> > DescriptorEvaluator::compute_u(const Neighbor& neighbor) const
{
    const double rsq = neighbor.x * neighbor.x + neighbor.y * neighbor.y + neighbor.z * neighbor.z;
    const double distance = std::sqrt(rsq);
    const double cutoff_radius = neighbor.cutoff_radius > 0.0 ? neighbor.cutoff_radius : parameters_.rcut;
    std::vector<std::complex<double> > result(u_size_, std::complex<double>(0.0, 0.0));
    if (distance <= 1.0e-14 || cutoff(distance, cutoff_radius) == 0.0)
    {
        return result;
    }

    const double theta0 = (distance - parameters_.rmin) * parameters_.rfac0 * ModuleBase::PI
                          / (cutoff_radius - parameters_.rmin);
    const double z0 = distance / std::tan(theta0);
    const double r0inv = 1.0 / std::sqrt(rsq + z0 * z0);
    const double ar = r0inv * z0;
    const double ai = -r0inv * neighbor.z;
    const double br = r0inv * neighbor.y;
    const double bi = -r0inv * neighbor.x;
    const int dimension = parameters_.twojmax + 1;

    result[0] = std::complex<double>(1.0, 0.0);
    for (int j = 1; j <= parameters_.twojmax; ++j)
    {
        for (int mb = 0; 2 * mb <= j; ++mb)
        {
            result[u_index(j, 0, mb)] = std::complex<double>(0.0, 0.0);
            for (int ma = 0; ma < j; ++ma)
            {
                const std::size_t current = u_index(j, ma, mb);
                const std::size_t previous = u_index(j - 1, ma, mb);
                const double first_root = rootpq_[static_cast<std::size_t>((j - ma) * dimension + (j - mb))];
                const double second_root = rootpq_[static_cast<std::size_t>((ma + 1) * dimension + (j - mb))];
                result[current] += first_root * std::complex<double>(ar, -ai) * result[previous];
                result[current + 1] = -second_root * std::complex<double>(br, -bi) * result[previous];
            }
        }

        for (int mb = 0; 2 * mb <= j; ++mb)
        {
            for (int ma = 0; ma <= j; ++ma)
            {
                const std::size_t left = u_index(j, ma, mb);
                const std::size_t right = u_index(j, j - ma, j - mb);
                const double sign = ((ma - mb) % 2 == 0) ? 1.0 : -1.0;
                result[right] = sign * std::conj(result[left]);
            }
        }
    }
    return result;
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

std::vector<std::complex<double> > DescriptorEvaluator::accumulate_u(const std::vector<Neighbor>& neighbors) const
{
    return accumulate_u(neighbors, 0);
}

std::vector<std::complex<double> > DescriptorEvaluator::accumulate_u(const std::vector<Neighbor>& neighbors,
                                                                      const std::size_t central_element) const
{
    const std::size_t channels = channel_count();
    const std::size_t self_channel = parameters_.chem_flag ? central_element : 0;
    if (self_channel >= channels)
    {
        throw std::out_of_range("SNAP central element is outside the descriptor channels.");
    }
    std::vector<std::complex<double> > result(channels * u_size_, std::complex<double>(0.0, 0.0));
    for (std::size_t channel = 0; channel < channels; ++channel)
    {
        if (channel != self_channel && !parameters_.wselfall_flag)
        {
            continue;
        }
        for (int j = 0; j <= parameters_.twojmax; ++j)
        {
            for (int ma = 0; ma <= j; ++ma)
            {
                result[channel * u_size_ + u_index(j, ma, ma)] = parameters_.wself;
            }
        }
    }

    for (std::size_t index = 0; index < neighbors.size(); ++index)
    {
        const Neighbor& neighbor = neighbors[index];
        const std::size_t channel = parameters_.chem_flag ? neighbor.element_index : 0;
        if (channel >= channels)
        {
            throw std::out_of_range("SNAP neighbor element is outside the descriptor channels.");
        }
        const double distance = std::sqrt(neighbor.x * neighbor.x + neighbor.y * neighbor.y + neighbor.z * neighbor.z);
        const double cutoff_radius = neighbor.cutoff_radius > 0.0 ? neighbor.cutoff_radius : parameters_.rcut;
        const double scale = cutoff(distance, cutoff_radius) * neighbor.weight;
        if (scale == 0.0)
        {
            continue;
        }
        const std::vector<std::complex<double> > one_neighbor = compute_u(neighbor);
        const std::size_t offset = channel * u_size_;
        for (std::size_t i = 0; i < one_neighbor.size(); ++i)
        {
            result[offset + i] += scale * one_neighbor[i];
        }
    }
    return result;
}

std::vector<double> DescriptorEvaluator::compute_bispectrum(const std::vector<std::complex<double> >& u) const
{
    const std::size_t channels = channel_count();
    if (u.size() != channels * u_size_)
    {
        throw std::invalid_argument("SNAP u array has an unexpected size.");
    }
    const std::size_t bispectrum_size = bispectrum_indices_.size();
    std::vector<double> result(channels * channels * channels * bispectrum_size, 0.0);
    const auto compute_one = [this, &u](const std::size_t first_channel,
                                         const std::size_t second_channel,
                                         const std::size_t third_channel,
                                         const std::size_t bispectrum_index) {
        const BispectrumIndex& index = bispectrum_indices_[bispectrum_index];
        const std::size_t first_offset = first_channel * u_size_;
        const std::size_t second_offset = second_channel * u_size_;
        const std::size_t third_offset = third_channel * u_size_;
        const auto compute_z = [this, &u, &index, first_offset, second_offset](const int ma, const int mb) {
        std::complex<double> z(0.0, 0.0);
        const int m = 2 * ma - index.j;
        const int mp = 2 * mb - index.j;
        for (int ma1 = 0; ma1 <= index.j1; ++ma1)
        {
            const int m1 = 2 * ma1 - index.j1;
            for (int ma2 = 0; ma2 <= index.j2; ++ma2)
            {
                const int m2 = 2 * ma2 - index.j2;
                const double c1 = cg_table_.coefficient(index.j1, m1, index.j2, m2, index.j, m);
                if (c1 == 0.0)
                {
                    continue;
                }
                for (int mb1 = 0; mb1 <= index.j1; ++mb1)
                {
                    const int mp1 = 2 * mb1 - index.j1;
                    for (int mb2 = 0; mb2 <= index.j2; ++mb2)
                    {
                        const int mp2 = 2 * mb2 - index.j2;
                        const double c2 = cg_table_.coefficient(index.j1, mp1, index.j2, mp2, index.j, mp);
                        if (c2 != 0.0)
                        {
                            z += c1 * c2
                                 * u[first_offset + u_index(index.j1, ma1, mb1)]
                                 * u[second_offset + u_index(index.j2, ma2, mb2)];
                        }
                    }
                }
            }
        }
        return z;
    };
        std::complex<double> value(0.0, 0.0);
        for (int mb = 0; 2 * mb < index.j; ++mb)
        {
            for (int ma = 0; ma <= index.j; ++ma)
            {
                value += std::conj(u[third_offset + u_index(index.j, ma, mb)]) * compute_z(ma, mb);
            }
        }
        if (index.j % 2 == 0)
        {
            const int mb = index.j / 2;
            for (int ma = 0; ma < mb; ++ma)
            {
                value += std::conj(u[third_offset + u_index(index.j, ma, mb)]) * compute_z(ma, mb);
            }
            value += 0.5 * std::conj(u[third_offset + u_index(index.j, mb, mb)]) * compute_z(mb, mb);
        }
        return 2.0 * value.real() / (parameters_.bnorm_flag ? index.j + 1.0 : 1.0);
    };
    for (std::size_t first_channel = 0; first_channel < channels; ++first_channel)
    {
        for (std::size_t second_channel = 0; second_channel < channels; ++second_channel)
        {
            for (std::size_t third_channel = 0; third_channel < channels; ++third_channel)
            {
                const std::size_t triple = (first_channel * channels + second_channel) * channels + third_channel;
                for (std::size_t ib = 0; ib < bispectrum_size; ++ib)
                {
                    result[triple * bispectrum_size + ib]
                        = compute_one(first_channel, second_channel, third_channel, ib);
                }
            }
        }
    }
    return result;
}

std::vector<double> DescriptorEvaluator::compute_bispectrum_derivative(
    const std::vector<std::complex<double> >& u,
    const std::vector<std::complex<double> >& du) const
{
    const std::size_t channels = channel_count();
    if (u.size() != channels * u_size_ || du.size() != channels * u_size_)
    {
        throw std::invalid_argument("SNAP u derivative array has an unexpected size.");
    }
    const std::size_t bispectrum_size = bispectrum_indices_.size();
    std::vector<double> result(channels * channels * channels * bispectrum_size, 0.0);
    const auto compute_one = [this, &u, &du](const std::size_t first_channel,
                                               const std::size_t second_channel,
                                               const std::size_t third_channel,
                                               const std::size_t bispectrum_index) {
        const BispectrumIndex& index = bispectrum_indices_[bispectrum_index];
        const std::size_t first_offset = first_channel * u_size_;
        const std::size_t second_offset = second_channel * u_size_;
        const std::size_t third_offset = third_channel * u_size_;
        const auto compute_z = [this, &u, &du, &index, first_offset, second_offset](const int ma,
                                                                                       const int mb,
                                                                                       std::complex<double>& z,
                                                                                       std::complex<double>& dz) {
        z = std::complex<double>(0.0, 0.0);
        dz = std::complex<double>(0.0, 0.0);
        const int m = 2 * ma - index.j;
        const int mp = 2 * mb - index.j;
        for (int ma1 = 0; ma1 <= index.j1; ++ma1)
        {
            const int m1 = 2 * ma1 - index.j1;
            for (int ma2 = 0; ma2 <= index.j2; ++ma2)
            {
                const int m2 = 2 * ma2 - index.j2;
                const double c1 = cg_table_.coefficient(index.j1, m1, index.j2, m2, index.j, m);
                if (c1 == 0.0)
                {
                    continue;
                }
                for (int mb1 = 0; mb1 <= index.j1; ++mb1)
                {
                    const int mp1 = 2 * mb1 - index.j1;
                    for (int mb2 = 0; mb2 <= index.j2; ++mb2)
                    {
                        const int mp2 = 2 * mb2 - index.j2;
                        const double c2 = cg_table_.coefficient(index.j1, mp1, index.j2, mp2, index.j, mp);
                        if (c2 == 0.0)
                        {
                            continue;
                        }
                        const double coefficient = c1 * c2;
                        const std::complex<double> u1 = u[first_offset + u_index(index.j1, ma1, mb1)];
                        const std::complex<double> u2 = u[second_offset + u_index(index.j2, ma2, mb2)];
                        const std::complex<double> du1 = du[first_offset + u_index(index.j1, ma1, mb1)];
                        const std::complex<double> du2 = du[second_offset + u_index(index.j2, ma2, mb2)];
                        z += coefficient * u1 * u2;
                        dz += coefficient * (du1 * u2 + u1 * du2);
                    }
                }
            }
        }
        };
        std::complex<double> derivative(0.0, 0.0);
        for (int mb = 0; 2 * mb < index.j; ++mb)
        {
            for (int ma = 0; ma <= index.j; ++ma)
            {
                std::complex<double> z(0.0, 0.0);
                std::complex<double> dz(0.0, 0.0);
                compute_z(ma, mb, z, dz);
                const std::complex<double> du3 = du[third_offset + u_index(index.j, ma, mb)];
                derivative += std::conj(du3) * z
                              + std::conj(u[third_offset + u_index(index.j, ma, mb)]) * dz;
            }
        }
        if (index.j % 2 == 0)
        {
            const int mb = index.j / 2;
            for (int ma = 0; ma < mb; ++ma)
            {
                std::complex<double> z(0.0, 0.0);
                std::complex<double> dz(0.0, 0.0);
                compute_z(ma, mb, z, dz);
                const std::complex<double> du3 = du[third_offset + u_index(index.j, ma, mb)];
                derivative += std::conj(du3) * z
                              + std::conj(u[third_offset + u_index(index.j, ma, mb)]) * dz;
            }
            std::complex<double> z(0.0, 0.0);
            std::complex<double> dz(0.0, 0.0);
            compute_z(mb, mb, z, dz);
            const std::complex<double> du3 = du[third_offset + u_index(index.j, mb, mb)];
            derivative += 0.5 * (std::conj(du3) * z
                                 + std::conj(u[third_offset + u_index(index.j, mb, mb)]) * dz);
        }
        return 2.0 * derivative.real() / (parameters_.bnorm_flag ? index.j + 1.0 : 1.0);
    };
    for (std::size_t first_channel = 0; first_channel < channels; ++first_channel)
    {
        for (std::size_t second_channel = 0; second_channel < channels; ++second_channel)
        {
            for (std::size_t third_channel = 0; third_channel < channels; ++third_channel)
            {
                const std::size_t triple = (first_channel * channels + second_channel) * channels + third_channel;
                for (std::size_t ib = 0; ib < bispectrum_size; ++ib)
                {
                    result[triple * bispectrum_size + ib]
                        = compute_one(first_channel, second_channel, third_channel, ib);
                }
            }
        }
    }
    return result;
}

} // namespace SNAP
} // namespace ModuleMD
