#include "rpa_abfs_preorthogonalization.h"

#include "source_base/constants.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace ModuleRI
{
namespace
{

constexpr double kGridTolerance = 1.0e-12;
constexpr double kOrthogonalityTolerance = 1.0e-20;
constexpr int kMaximumProjectionPasses = 4;
constexpr double kIdentityTolerance = 1.0e-10;

int expanded_size(const RpaAbfsOrbitalSet& abfs)
{
    int result = 0;
    for (const auto& type: abfs)
    {
        for (std::size_t l = 0; l != type.size(); ++l)
        {
            result += static_cast<int>(type[l].size() * (2 * l + 1));
        }
    }
    return result;
}

std::vector<double> simpson_weights(const Numerical_Orbital_Lm& reference)
{
    const int nr = reference.getNr();
    if (nr < 3 || nr % 2 == 0)
    {
        throw std::invalid_argument(
            "RPA ABFS Coulomb preorthogonalization requires an odd radial mesh with at least three points");
    }

    std::vector<double> weights(nr, 0.0);
    for (int ir = 0; ir != nr; ++ir)
    {
        const double coefficient = (ir == 0 || ir == nr - 1) ? 1.0 : (ir % 2 == 0 ? 2.0 : 4.0);
        weights[ir] = coefficient * reference.getRab(ir) / 3.0;
    }
    return weights;
}

void validate_channel(const std::vector<Numerical_Orbital_Lm>& channel)
{
    if (channel.empty())
    {
        return;
    }
    const Numerical_Orbital_Lm& reference = channel.front();
    for (std::size_t i = 0; i != channel.size(); ++i)
    {
        const Numerical_Orbital_Lm& orbital = channel[i];
        if (orbital.getType() != reference.getType() || orbital.getL() != reference.getL()
            || orbital.getNr() != reference.getNr())
        {
            throw std::invalid_argument(
                "RPA ABFS Coulomb preorthogonalization requires one consistent (type,l) radial channel");
        }
        for (int ir = 0; ir != reference.getNr(); ++ir)
        {
            if (std::abs(orbital.getRadial(ir) - reference.getRadial(ir))
                    > kGridTolerance * std::max(1.0, std::abs(reference.getRadial(ir)))
                || std::abs(orbital.getRab(ir) - reference.getRab(ir))
                       > kGridTolerance * std::max(1.0, std::abs(reference.getRab(ir))))
            {
                throw std::invalid_argument(
                    "RPA ABFS Coulomb preorthogonalization requires identical radial grids within a channel");
            }
        }
    }
}

double metric_product(const std::vector<double>& lhs,
                      const std::vector<std::vector<double>>& metric,
                      const std::vector<double>& rhs)
{
    double result = 0.0;
    for (std::size_t i = 0; i != lhs.size(); ++i)
    {
        for (std::size_t j = 0; j != rhs.size(); ++j)
        {
            result += lhs[i] * metric[i][j] * rhs[j];
        }
    }
    return result;
}

std::vector<double> metric_overlaps(const std::vector<std::vector<double>>& accepted,
                                    const std::vector<std::vector<double>>& metric,
                                    const std::vector<double>& candidate)
{
    std::vector<double> overlaps;
    overlaps.reserve(accepted.size());
    for (const auto& vector: accepted)
    {
        overlaps.push_back(metric_product(vector, metric, candidate));
    }
    return overlaps;
}

double squared_norm(const std::vector<double>& values)
{
    return std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
}

Numerical_Orbital_Lm rebuild_orbital(const std::vector<Numerical_Orbital_Lm>& original,
                                     const std::vector<double>& coefficients,
                                     const int output_index,
                                     const bool force_flag)
{
    const Numerical_Orbital_Lm& reference = original.front();
    std::vector<double> psi(reference.getNr(), 0.0);
    for (std::size_t input_index = 0; input_index != original.size(); ++input_index)
    {
        for (int ir = 0; ir != reference.getNr(); ++ir)
        {
            psi[ir] += coefficients[input_index] * original[input_index].getPsi(ir);
        }
    }

    Numerical_Orbital_Lm output;
    output.set_orbital_info(reference.getLabel(),
                            reference.getType(),
                            reference.getL(),
                            output_index,
                            reference.getNr(),
                            reference.getRab(),
                            reference.getRadial(),
                            Numerical_Orbital_Lm::Psi_Type::Psi,
                            psi.data(),
                            reference.getNk(),
                            reference.getDk(),
                            reference.getDruniform(),
                            false,
                            true,
                            force_flag);
    return output;
}

RpaAbfsPreorthChannelReport preorthogonalize_channel(std::vector<Numerical_Orbital_Lm>& channel,
                                                     const double threshold,
                                                     const bool force_flag)
{
    validate_channel(channel);
    const int atom_type = channel.front().getType();
    const int angular_momentum = channel.front().getL();
    const auto metric = compute_onsite_coulomb_metric(channel);
    std::vector<std::vector<double>> accepted;
    std::vector<int> rejected;
    double minimum_residual = std::numeric_limits<double>::infinity();

    for (std::size_t input_index = 0; input_index != channel.size(); ++input_index)
    {
        std::vector<double> candidate(channel.size(), 0.0);
        candidate[input_index] = 1.0;
        const double input_norm = metric_product(candidate, metric, candidate);
        if (!std::isfinite(input_norm) || input_norm <= 0.0)
        {
            throw std::invalid_argument(
                "RPA ABFS Coulomb preorthogonalization found a non-positive input norm for type "
                + std::to_string(atom_type) + ", l " + std::to_string(angular_momentum) + ", radial "
                + std::to_string(input_index));
        }
        for (double& value: candidate)
        {
            value /= std::sqrt(input_norm);
        }

        std::vector<double> overlaps = metric_overlaps(accepted, metric, candidate);
        double residual_norm = 1.0 - squared_norm(overlaps);
        if (!std::isfinite(residual_norm) || residual_norm < -threshold)
        {
            throw std::runtime_error("RPA ABFS Coulomb preorthogonalization produced an invalid residual norm");
        }
        residual_norm = std::max(0.0, residual_norm);
        minimum_residual = std::min(minimum_residual, residual_norm);
        if (residual_norm <= threshold)
        {
            rejected.push_back(static_cast<int>(input_index));
            continue;
        }

        bool orthogonal = accepted.empty();
        for (int pass = 0; pass != kMaximumProjectionPasses && !orthogonal; ++pass)
        {
            for (std::size_t i = 0; i != accepted.size(); ++i)
            {
                for (std::size_t j = 0; j != candidate.size(); ++j)
                {
                    candidate[j] -= overlaps[i] * accepted[i][j];
                }
            }
            const double norm = metric_product(candidate, metric, candidate);
            if (!std::isfinite(norm) || norm <= 0.0)
            {
                throw std::invalid_argument("RPA ABFS Coulomb preorthogonalization lost a positive residual norm");
            }
            for (double& value: candidate)
            {
                value /= std::sqrt(norm);
            }
            overlaps = metric_overlaps(accepted, metric, candidate);
            orthogonal = squared_norm(overlaps) <= kOrthogonalityTolerance;
        }
        if (!orthogonal)
        {
            throw std::runtime_error("RPA ABFS Coulomb preorthogonalization did not reach the orthogonality tolerance");
        }
        accepted.push_back(std::move(candidate));
    }

    if (accepted.empty())
    {
        throw std::invalid_argument("RPA ABFS Coulomb preorthogonalization rejected an entire radial channel");
    }

    const std::vector<Numerical_Orbital_Lm> original = channel;
    channel.clear();
    channel.reserve(accepted.size());
    for (std::size_t output_index = 0; output_index != accepted.size(); ++output_index)
    {
        channel.push_back(
            rebuild_orbital(original, accepted[output_index], static_cast<int>(output_index), force_flag));
    }

    const auto transformed_metric = compute_onsite_coulomb_metric(channel);
    double maximum_identity_error = 0.0;
    for (std::size_t i = 0; i != transformed_metric.size(); ++i)
    {
        for (std::size_t j = 0; j != transformed_metric.size(); ++j)
        {
            const double target = i == j ? 1.0 : 0.0;
            maximum_identity_error = std::max(maximum_identity_error, std::abs(transformed_metric[i][j] - target));
        }
    }
    if (!std::isfinite(maximum_identity_error) || maximum_identity_error > kIdentityTolerance)
    {
        throw std::runtime_error("RPA ABFS Coulomb preorthogonalization failed its on-site identity check");
    }

    RpaAbfsPreorthChannelReport report;
    report.atom_type = atom_type;
    report.angular_momentum = angular_momentum;
    report.input_count = static_cast<int>(original.size());
    report.output_count = static_cast<int>(channel.size());
    report.rejected_indices = std::move(rejected);
    report.minimum_residual_norm = minimum_residual;
    report.maximum_identity_error = maximum_identity_error;
    return report;
}

} // namespace

std::vector<std::vector<double>> compute_onsite_coulomb_metric(const std::vector<Numerical_Orbital_Lm>& channel)
{
    validate_channel(channel);
    if (channel.empty())
    {
        return {};
    }

    const Numerical_Orbital_Lm& reference = channel.front();
    const int nr = reference.getNr();
    const int l = reference.getL();
    const std::vector<double> weights = simpson_weights(reference);
    std::vector<std::vector<double>> metric(channel.size(), std::vector<double>(channel.size(), 0.0));

    for (std::size_t b = 0; b != channel.size(); ++b)
    {
        std::vector<double> prefix(nr, 0.0);
        std::vector<double> suffix(nr, 0.0);
        for (int ir = 1; ir != nr; ++ir)
        {
            const double r = reference.getRadial(ir);
            if (!(r > 0.0) || !std::isfinite(r))
            {
                throw std::invalid_argument(
                    "RPA ABFS Coulomb preorthogonalization requires positive finite radii after the origin");
            }
            prefix[ir] = prefix[ir - 1] + weights[ir] * std::pow(r, l + 2) * channel[b].getPsi(ir);
        }
        for (int ir = nr - 2; ir >= 1; --ir)
        {
            const double r_next = reference.getRadial(ir + 1);
            suffix[ir] = suffix[ir + 1] + weights[ir + 1] * std::pow(r_next, 1 - l) * channel[b].getPsi(ir + 1);
        }

        for (std::size_t a = 0; a != channel.size(); ++a)
        {
            double value = 0.0;
            for (int ir = 1; ir != nr; ++ir)
            {
                const double r = reference.getRadial(ir);
                const double potential = std::pow(r, -l - 1) * prefix[ir] + std::pow(r, l) * suffix[ir];
                value += weights[ir] * r * r * channel[a].getPsi(ir) * potential;
            }
            metric[a][b] = ModuleBase::FOUR_PI * value / (2 * l + 1);
        }
    }

    for (std::size_t i = 0; i != metric.size(); ++i)
    {
        for (std::size_t j = i + 1; j != metric.size(); ++j)
        {
            const double symmetric = 0.5 * (metric[i][j] + metric[j][i]);
            metric[i][j] = symmetric;
            metric[j][i] = symmetric;
        }
    }
    return metric;
}

RpaAbfsPreorthReport apply_rpa_abfs_preorthogonalization(RpaAbfsOrbitalSet& abfs,
                                                         const std::string& mode,
                                                         const double threshold,
                                                         const bool force_flag)
{
    if (mode != "none" && mode != "onsite_coulomb")
    {
        throw std::invalid_argument("rpa_abfs_preorth must be none or onsite_coulomb");
    }
    if (!std::isfinite(threshold) || threshold <= 0.0 || threshold >= 1.0)
    {
        throw std::invalid_argument("rpa_abfs_preorth_threshold must be finite and strictly inside (0,1)");
    }

    RpaAbfsPreorthReport report;
    report.mode = mode;
    report.threshold = threshold;
    report.input_expanded_size = expanded_size(abfs);
    report.output_expanded_size = report.input_expanded_size;
    if (mode == "none")
    {
        return report;
    }

    for (auto& type: abfs)
    {
        for (auto& channel: type)
        {
            if (!channel.empty())
            {
                report.channels.push_back(preorthogonalize_channel(channel, threshold, force_flag));
            }
        }
    }
    report.output_expanded_size = expanded_size(abfs);
    return report;
}

RpaAbfsPreorthReport finalize_rpa_abfs_from_input(RpaAbfsOrbitalSet& abfs,
                                                  const Input_para& input,
                                                  const bool force_flag)
{
    return apply_rpa_abfs_preorthogonalization(abfs,
                                               input.rpa_abfs_preorth,
                                               input.rpa_abfs_preorth_threshold,
                                               force_flag);
}

std::string format_rpa_abfs_preorth_report(const RpaAbfsPreorthReport& report)
{
    std::ostringstream output;
    output << std::setprecision(16) << "RPA_ABFS_PREORTH mode=" << report.mode << " threshold=" << report.threshold
           << " input_expanded_size=" << report.input_expanded_size
           << " output_expanded_size=" << report.output_expanded_size << '\n';
    for (const auto& channel: report.channels)
    {
        output << "RPA_ABFS_PREORTH_CHANNEL type=" << channel.atom_type << " l=" << channel.angular_momentum
               << " input_count=" << channel.input_count << " output_count=" << channel.output_count
               << " rejected_indices=";
        if (channel.rejected_indices.empty())
        {
            output << "none";
        }
        else
        {
            for (std::size_t i = 0; i != channel.rejected_indices.size(); ++i)
            {
                output << (i == 0 ? "" : ",") << channel.rejected_indices[i];
            }
        }
        output << " minimum_residual_norm=" << channel.minimum_residual_norm
               << " maximum_identity_error=" << channel.maximum_identity_error << '\n';
    }
    return output.str();
}

} // namespace ModuleRI
