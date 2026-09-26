#include "xc_rvv10.h"

#include "source_base/constants.h"

#include <cmath>
#include <stdexcept>
#include <vector>

// Independent implementation of the published rVV10 equations:
// Sabatini, Gorni, de Gironcoli, PRB 87, 041108(R) (2013),
// DOI: 10.1103/PhysRevB.87.041108. No QE routines are incorporated.
// See docs/developers_guide/rvv10_core.md for units and cutoff semantics.
namespace Rvv10
{
LocalModel::LocalModel(double b, double c)
{
    if (!std::isfinite(b) || b <= 0 || !std::isfinite(c) || c < 0)
    {
        throw std::invalid_argument("rVV10 requires finite b > 0 and C >= 0");
    }
    const double pi = ModuleBase::PI;
    k_prefactor_ = 3 * pi * b / std::pow(9 * pi, 1.0 / 6.0);
    gap_prefactor_ = 2 * std::sqrt(c);
    const double b_power = std::pow(b, -1.5);
    weight_prefactor_ = b_power / (3 * std::pow(pi, 1.25));
    beta_ = std::pow(3.0, 0.75) * b_power / 16;
    if (!std::isfinite(k_prefactor_) || k_prefactor_ <= 0 || !std::isfinite(weight_prefactor_) || weight_prefactor_ <= 0
        || !std::isfinite(beta_) || beta_ <= 0)
    {
        throw std::overflow_error("rVV10 parameter scale is not representable in double precision");
    }
}

double LocalModel::beta() const
{
    return beta_;
}

LocalValues LocalModel::evaluate(double n, double sigma) const
{
    if (!std::isfinite(n) || !std::isfinite(sigma) || sigma < 0)
    {
        throw std::invalid_argument("rVV10 requires finite density and finite sigma >= 0");
    }
    LocalValues v = {q_cut, 0.0, 0.0, 0.0, 0.0};
    if (n <= density_cutoff)
    {
        return v;
    }

    v.weight = weight_prefactor_ * std::pow(n, 0.75);
    v.dweight_dn = 0.75 * v.weight / n;
    if (!std::isfinite(v.weight) || !std::isfinite(v.dweight_dn))
    {
        throw std::overflow_error("rVV10 density weight is not representable");
    }

    // omega_Ry^2 = 16*pi*n/3 + 4*C*sigma^2/n^4; q = omega_Ry/kappa_Ry.
    // Work with the two components of q to avoid explicitly squaring sigma.
    const double plasma_q = std::sqrt(16 * ModuleBase::PI / 3) * std::cbrt(n) / k_prefactor_;
    if (plasma_q >= 2 * q_cut)
    {
        return v;
    }
    const double kappa = k_prefactor_ * std::pow(n, 1.0 / 6.0);
    if (!std::isfinite(kappa) || kappa <= 0)
    {
        throw std::overflow_error("rVV10 kappa is not representable");
    }
    const double gradient_q = (sigma == 0 || gap_prefactor_ == 0)
                                  ? 0.0
                                  : (gap_prefactor_ / kappa) * (sigma / n) / n;
    const double q = std::hypot(plasma_q, gradient_q);
    if (std::isnan(q))
    {
        throw std::overflow_error("rVV10 local q is not representable");
    }
    // At q/q_cut >= 2 the 12-term exponent exceeds 700: q0 rounds to q_cut.
    // This also handles +infinity from an enormous but finite sigma safely.
    if (q >= 2 * q_cut)
    {
        return v;
    }

    const double x = q / q_cut;
    double power = 1.0;
    double exponent = 0.0;
    double slope = 0.0;
    for (int m = 1; m <= 12; ++m)
    {
        slope += power;
        power *= x;
        exponent += power / m;
    }
    v.q0 = -q_cut * std::expm1(-exponent);
    if (v.q0 <= q_min)
        v.q0 = q_min;
    const double saturation_derivative = slope * std::exp(-exponent);
    const double fraction = gradient_q / q;
    // q_n = (q/n)*(1/3 - (5/2)*omega_gap^2/omega^2), at fixed sigma.
    v.dq0_dn = saturation_derivative * (q / n) * (1.0 / 3.0 - 2.5 * fraction * fraction);
    v.dq0_dsigma = sigma == 0 ? 0.0
                              : saturation_derivative * fraction * (gap_prefactor_ / kappa) / n / n;
    if (!std::isfinite(v.dq0_dn) || !std::isfinite(v.dq0_dsigma))
    {
        throw std::overflow_error("rVV10 q derivatives are not representable");
    }
    return v;
}

const std::array<double, channel_count>& q_mesh()
{
    // Numerical mesh used by the pinned QE reference; not an algorithm copy.
    static const std::array<double, channel_count> mesh = {{q_min,
                                                            3.e-4,
                                                            5.893850845618885e-4,
                                                            1.008103720396345e-3,
                                                            1.613958359589310e-3,
                                                            2.490584839564653e-3,
                                                            3.758997979748929e-3,
                                                            5.594297198907115e-3,
                                                            8.249838297569416e-3,
                                                            1.209220822453922e-2,
                                                            1.765183095571029e-2,
                                                            2.569619042667097e-2,
                                                            3.733577865542191e-2,
                                                            5.417739477463518e-2,
                                                            7.854595729872216e-2,
                                                            0.113805449932145,
                                                            0.164823306218807,
                                                            0.238642339497217,
                                                            0.345452975434964,
                                                            q_cut}};
    return mesh;
}

double real_kernel(double q1, double q2, double distance)
{
    if (!std::isfinite(q1) || !std::isfinite(q2) || !std::isfinite(distance) || q1 <= 0 || q2 <= 0 || distance < 0)
    {
        throw std::invalid_argument("rVV10 kernel requires finite q1,q2 > 0 and distance >= 0");
    }
    // Form sqrt(q)*R first.  Computing R^2 can overflow even when qR^2 is
    // representable (the kernel is routinely evaluated at very large R).
    const double qr1 = std::sqrt(q1) * distance;
    const double qr2 = std::sqrt(q2) * distance;
    const double a = 1 + qr1 * qr1;
    const double b = 1 + qr2 * qr2;
    return ((-24.0 / a) / b) / (a + b);
}

namespace
{
constexpr int radial_intervals = 1024;
constexpr double radial_limit = 100.0;
const double reciprocal_step = 2 * ModuleBase::PI / radial_limit;
} // namespace

KernelTable::KernelTable() : spline_(radial_intervals + 1, 0.0, reciprocal_step)
{
    // Fixed numerical quadrature, independent of b,C and the simulation cell.
    // Precompute radial transform weights once for all 210 independent pairs.
    const int size = radial_intervals + 1;
    const double dr = radial_limit / radial_intervals;
    std::vector<double> transform(size * size, 0.0);
    for (int k = 0; k < size; ++k)
    {
        const double g = k * reciprocal_step;
        for (int j = 1; j < size; ++j)
        {
            const double r = j * dr;
            const double weight = j == radial_intervals ? 0.5 : 1.0;
            transform[k * size + j] = 4 * ModuleBase::PI * dr * weight * (k == 0 ? r * r : r * std::sin(g * r) / g);
        }
    }
    typedef ModuleBase::CubicSpline Spline;
    const Spline::BoundaryCondition natural(Spline::BoundaryType::second_deriv, 0.0);
    spline_.reserve(pair_count);
    std::vector<double> radial(size);
    std::vector<double> reciprocal(size);
    int pair = 0;
    for (int a = 0; a < channel_count; ++a)
    {
        for (int b = 0; b <= a; ++b)
        {
            for (int j = 0; j < size; ++j)
            {
                radial[j] = real_kernel(q_mesh()[a], q_mesh()[b], j * dr);
            }
            for (int k = 0; k < size; ++k)
            {
                double sum = 0.0;
                for (int j = 1; j < size; ++j)
                {
                    sum += transform[k * size + j] * radial[j];
                }
                reciprocal[k] = sum;
            }
            indices_[pair] = pair;
            ++pair;
            spline_.add(reciprocal.data(), natural, natural);
        }
    }
}

double KernelTable::maximum_g() const
{
    return radial_intervals * reciprocal_step;
}

KernelValues KernelTable::evaluate(double g) const
{
    KernelValues result;
    interpolate(g, result.value.data(), result.derivative.data());
    return result;
}

std::array<double, channel_count * channel_count> KernelTable::values(double g) const
{
    std::array<double, channel_count * channel_count> result;
    interpolate(g, result.data(), nullptr);
    return result;
}

void KernelTable::interpolate(double g, double* value, double* derivative) const
{
    if (!std::isfinite(g) || g < 0 || g >= maximum_g())
    {
        throw std::out_of_range("rVV10 |G| is outside the kernel table; extrapolation is unsupported");
    }
    std::array<double, pair_count> values;
    std::array<double, pair_count> derivatives;
    spline_.multi_eval(pair_count, indices_.data(), g, values.data(), derivative ? derivatives.data() : nullptr, nullptr);
    int pair = 0;
    for (int a = 0; a < channel_count; ++a)
    {
        for (int b = 0; b <= a; ++b)
        {
            value[a * channel_count + b] = value[b * channel_count + a] = values[pair];
            if (derivative)
                derivative[a * channel_count + b] = derivative[b * channel_count + a] = derivatives[pair];
            ++pair;
        }
    }
}

SplineBasis::SplineBasis() : spline_(channel_count, q_mesh().data())
{
    typedef ModuleBase::CubicSpline Spline;
    const Spline::BoundaryCondition natural(Spline::BoundaryType::second_deriv, 0.0);
    spline_.reserve(channel_count);
    for (int a = 0; a < channel_count; ++a)
    {
        std::array<double, channel_count> cardinal = {};
        cardinal[a] = 1.0;
        spline_.add(cardinal.data(), natural, natural);
    }
}

Channel SplineBasis::channel(const LocalValues& local, int index) const
{
    if (index < 0 || index >= channel_count)
        throw std::out_of_range("rVV10 channel index is outside its q mesh");
    if (local.weight == 0)
        return {0.0, 0.0, 0.0};
    if (!std::isfinite(local.q0) || local.q0 < q_min || local.q0 > q_cut)
        throw std::out_of_range("rVV10 spline argument is outside its q mesh");
    double p;
    double dp;
    // A channel-major FFT loop needs one spline, not all 20 at every pass.
    // Keep the existing ABACUS evaluator and avoid a full real-space cache.
    spline_.multi_eval(1, &index, local.q0, &p, &dp, nullptr);
    // Multiply q derivatives by p' before the weight so a zero derivative in a
    // saturated or density-cutoff point cannot form an unused inf*0 product.
    return {local.weight * p,
            local.dweight_dn * p + (local.dq0_dn * dp) * local.weight,
            (local.dq0_dsigma * dp) * local.weight};
}
} // namespace Rvv10
