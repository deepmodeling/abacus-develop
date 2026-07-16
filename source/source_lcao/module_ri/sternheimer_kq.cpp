#include "sternheimer_kq.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace ModuleRI
{
namespace
{

void validate_finite_kpoint(const SternheimerReducedKPoint& kpoint, const char* name)
{
    for (const double coordinate: kpoint)
    {
        if (!std::isfinite(coordinate))
        {
            throw std::invalid_argument(std::string("Sternheimer ") + name + " coordinates must be finite.");
        }
    }
}

} // namespace

SternheimerFoldedKPoint fold_sternheimer_kpoint(const SternheimerReducedKPoint& kpoint)
{
    validate_finite_kpoint(kpoint, "k-point");

    SternheimerFoldedKPoint result{};
    for (std::size_t direction = 0; direction != kpoint.size(); ++direction)
    {
        const double shift = std::floor(kpoint[direction] + 0.5);
        if (shift < static_cast<double>(std::numeric_limits<int>::min())
            || shift > static_cast<double>(std::numeric_limits<int>::max()))
        {
            throw std::invalid_argument("Sternheimer reciprocal-lattice shift exceeds integer range.");
        }
        result.kpoint[direction] = kpoint[direction] - shift;
        result.reciprocal_shift[direction] = static_cast<int>(shift);
    }
    return result;
}

double periodic_sternheimer_kpoint_distance(const SternheimerReducedKPoint& lhs, const SternheimerReducedKPoint& rhs)
{
    validate_finite_kpoint(lhs, "left k-point");
    validate_finite_kpoint(rhs, "right k-point");

    double distance = 0.0;
    for (std::size_t direction = 0; direction != lhs.size(); ++direction)
    {
        const double difference = lhs[direction] - rhs[direction];
        distance = std::max(distance, std::abs(difference - std::round(difference)));
    }
    return distance;
}

std::complex<double> sternheimer_bloch_phase(const SternheimerReducedKPoint& kpoint,
                                             const std::array<int, 3>& lattice_translation)
{
    validate_finite_kpoint(kpoint, "Bloch k-point");

    double reduced_phase = 0.0;
    for (std::size_t direction = 0; direction != kpoint.size(); ++direction)
    {
        reduced_phase += kpoint[direction] * lattice_translation[direction];
    }
    const double argument = 2.0 * std::acos(-1.0) * reduced_phase;
    return {std::cos(argument), std::sin(argument)};
}

std::vector<SternheimerKQPair> build_sternheimer_kq_map(const std::vector<SternheimerReducedKPoint>& kpoints,
                                                        const SternheimerReducedKPoint& qpoint,
                                                        const double tolerance)
{
    if (kpoints.empty())
    {
        throw std::invalid_argument("Sternheimer k+q mapping requires a nonempty k-point mesh.");
    }
    if (!std::isfinite(tolerance) || tolerance <= 0.0)
    {
        throw std::invalid_argument("Sternheimer k+q mapping tolerance must be positive and finite.");
    }
    validate_finite_kpoint(qpoint, "q-point");
    for (const auto& kpoint: kpoints)
    {
        validate_finite_kpoint(kpoint, "mesh k-point");
    }

    std::vector<SternheimerKQPair> mapping;
    mapping.reserve(kpoints.size());
    for (std::size_t source_index = 0; source_index != kpoints.size(); ++source_index)
    {
        SternheimerReducedKPoint k_plus_q{};
        for (std::size_t direction = 0; direction != k_plus_q.size(); ++direction)
        {
            k_plus_q[direction] = kpoints[source_index][direction] + qpoint[direction];
        }
        const SternheimerFoldedKPoint folded = fold_sternheimer_kpoint(k_plus_q);

        int target_index = -1;
        for (std::size_t candidate_index = 0; candidate_index != kpoints.size(); ++candidate_index)
        {
            if (periodic_sternheimer_kpoint_distance(folded.kpoint, kpoints[candidate_index]) <= tolerance)
            {
                if (target_index >= 0)
                {
                    throw std::invalid_argument("Sternheimer k+q mapping found duplicate periodic target k-points.");
                }
                target_index = static_cast<int>(candidate_index);
            }
        }
        if (target_index < 0)
        {
            throw std::invalid_argument("Sternheimer q-point is not commensurate with the supplied k-point mesh.");
        }

        const SternheimerReducedKPoint& selected_target = kpoints[static_cast<std::size_t>(target_index)];
        std::array<int, 3> reciprocal_shift{};
        for (std::size_t direction = 0; direction != k_plus_q.size(); ++direction)
        {
            const double shift = k_plus_q[direction] - selected_target[direction];
            const double rounded_shift = std::round(shift);
            if (std::abs(shift - rounded_shift) > tolerance
                || rounded_shift < static_cast<double>(std::numeric_limits<int>::min())
                || rounded_shift > static_cast<double>(std::numeric_limits<int>::max()))
            {
                throw std::runtime_error(
                    "Sternheimer k+q mapping could not express the selected target by a reciprocal-lattice shift.");
            }
            reciprocal_shift[direction] = static_cast<int>(rounded_shift);
        }
        mapping.push_back({static_cast<int>(source_index), target_index, selected_target, reciprocal_shift});
    }
    return mapping;
}

} // namespace ModuleRI
