#ifndef STERNHEIMER_WEAK_Q_UNIT_H
#define STERNHEIMER_WEAK_Q_UNIT_H

#include "source_lcao/module_ri/sternheimer_weak_matrix_audit.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace ModuleRI
{

// Pure one-based state/frequency selection. Columns and shard indices are zero-based.
struct SternheimerWeakQUnit
{
    struct Raw
    {
        const char* source_k = nullptr;
        const char* band_begin = nullptr;
        const char* band_end = nullptr;
        const char* frequency_begin = nullptr;
        const char* frequency_end = nullptr;
        const char* shard_index = nullptr;
        const char* shard_count = nullptr;
    };

    int source_k;
    int band_begin;
    int band_end;
    int frequency_begin;
    int frequency_end;
    SternheimerWeakAuditShard shard;

    static SternheimerWeakQUnit parse(const Raw& raw, const std::vector<int>& occupied_counts,
                                     const int nfreq, const int channels)
    {
        if (occupied_counts.empty() || nfreq <= 0
            || occupied_counts.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())
            || std::any_of(occupied_counts.begin(), occupied_counts.end(), [](int n) { return n <= 0; }))
            throw std::invalid_argument("Weak q unit requires full occupied counts and positive frequency count.");
        const int source = integer(raw.source_k, 1, "WEAK_Q_SOURCE_K");
        if (source > static_cast<int>(occupied_counts.size()))
            throw std::invalid_argument("WEAK_Q_SOURCE_K is outside the full k mesh.");
        const int bands = occupied_counts[source - 1];
        SternheimerWeakQUnit unit{source,
            integer(raw.band_begin, 1, "WEAK_Q_BAND_BEGIN"),
            integer(raw.band_end, bands, "WEAK_Q_BAND_END"),
            integer(raw.frequency_begin, 1, "WEAK_Q_FREQ_BEGIN"),
            integer(raw.frequency_end, nfreq, "WEAK_Q_FREQ_END"),
            SternheimerWeakAuditShard::parse(raw.shard_index, raw.shard_count, channels)};
        if (unit.band_begin > unit.band_end || unit.band_end > bands
            || unit.frequency_begin > unit.frequency_end || unit.frequency_end > nfreq)
            throw std::invalid_argument("Weak q unit band/frequency ranges must be inclusive and within the full input.");
        unit.expected_equations();
        return unit;
    }

    std::vector<int> columns() const { return shard.columns(); }
    bool all_source_bands(int count) const { return band_begin == 1 && band_end == count; }
    bool all_frequencies(int count) const { return frequency_begin == 1 && frequency_end == count; }

    std::uint64_t expected_equations() const
    {
        std::uint64_t result = 2;
        for (const auto count : {std::uint64_t(band_end - band_begin + 1),
                                 std::uint64_t(frequency_end - frequency_begin + 1),
                                 std::uint64_t(columns().size())})
        {
            if (count == 0 || result > std::numeric_limits<std::uint64_t>::max() / count)
                throw std::overflow_error("Weak q unit equation count overflow.");
            result *= count;
        }
        return result;
    }

    std::string stem(int iq) const
    {
        if (iq <= 0) throw std::invalid_argument("Weak q unit requires a positive q index.");
        return "STERNHEIMER_WEAK_Q_iq_" + std::to_string(iq) + "_k_" + std::to_string(source_k)
            + "_b_" + std::to_string(band_begin) + "_" + std::to_string(band_end)
            + "_f_" + std::to_string(frequency_begin) + "_" + std::to_string(frequency_end)
            + "_shard_" + std::to_string(shard.index) + "_of_" + std::to_string(shard.count);
    }

  private:
    static int integer(const char* raw, int fallback, const char* name)
    {
        if (raw == nullptr) return fallback;
        const auto invalid = [name]() {
            return std::invalid_argument(std::string(name) + " must be a positive decimal integer.");
        };
        if (*raw == '\0') throw invalid();
        int result = 0;
        for (const char* p = raw; *p; ++p)
        {
            if (*p < '0' || *p > '9' || result > (std::numeric_limits<int>::max() - (*p - '0')) / 10)
                throw invalid();
            result = 10 * result + (*p - '0');
        }
        if (result <= 0) throw invalid();
        return result;
    }
};

inline void validate_sternheimer_weak_q_fold(const std::array<double, 3>& source,
                                            const std::array<double, 3>& q,
                                            const std::array<double, 3>& target,
                                            const std::array<int, 3>& reciprocal_shift)
{
    for (int a = 0; a < 3; ++a)
    {
        const double shift = source[a] + q[a] - target[a];
        if (!std::isfinite(source[a]) || !std::isfinite(q[a]) || !std::isfinite(target[a])
            || !std::isfinite(shift) || std::abs(shift - reciprocal_shift[a]) > 1e-10)
            throw std::invalid_argument("Weak q unit requires canonical target and an integer k+q fold matching the plan.");
    }
}

inline double sternheimer_weak_q_weight(double kweight, double occupation)
{
    if (!std::isfinite(kweight) || kweight <= 0 || !std::isfinite(occupation)
        || std::abs(occupation - 1.0) > 1e-10 || !std::isfinite(kweight * occupation))
        throw std::invalid_argument("Weak q unit requires positive kweight and integer insulating occupation.");
    return kweight * occupation; // ABACUS kweight already includes spin degeneracy.
}

inline double sternheimer_weak_q_omega_ry(double omega_ha)
{
    if (!std::isfinite(omega_ha) || omega_ha <= 0 || !std::isfinite(2.0 * omega_ha))
        throw std::invalid_argument("Weak q unit requires finite positive imaginary frequency in Ha.");
    return 2.0 * omega_ha;
}

inline double sternheimer_weak_q_left_vertex_factor() { return 0.5; }

// Compact radial interpolation extends through the zero node AFTER the last nonzero node.
// Only exact zeros are ignored; this is not a density-tail truncation tolerance.
inline double sternheimer_weak_q_radial_support(const std::vector<double>& radii,
                                                const std::vector<double>& values)
{
    if (radii.empty() || radii.size() != values.size())
        throw std::invalid_argument("Weak q unit invalid radial density support.");
    std::size_t last = 0;
    bool nonzero = false;
    for (std::size_t i = 0; i < radii.size(); ++i)
    {
        if (!std::isfinite(radii[i]) || radii[i] < 0 || !std::isfinite(values[i])
            || (i > 0 && radii[i] <= radii[i - 1]))
            throw std::invalid_argument("Weak q unit radial support requires finite ordered radial data.");
        if (values[i] != 0) { nonzero = true; last = i; }
    }
    return nonzero ? radii[std::min(last + 1, radii.size() - 1)] : 0.0;
}

inline bool sternheimer_weak_q_z_support_contained(double center, double radius, double height)
{
    if (!std::isfinite(center) || !std::isfinite(radius) || radius < 0 || !std::isfinite(height) || height <= 0)
        throw std::invalid_argument("Weak q unit invalid slab support geometry.");
    const double margin = 1e-12 * height;
    return center - radius > margin && center + radius < height - margin;
}

inline void validate_sternheimer_weak_q_residual(bool converged, double relative, double absolute)
{
    if (!converged || !std::isfinite(relative) || relative < 0 || relative > 1e-8
        || !std::isfinite(absolute) || absolute < 0)
        throw std::runtime_error("Weak q unit original augmented residual did not pass 1e-8.");
}

// Each column already contains the Ha left-vertex contraction of ONE signed solve.
// Accumulate into column-major owned-column storage, with no Hermitian completion.
inline void accumulate_sternheimer_weak_q_column(std::vector<std::complex<double>>& partial,
                                                const std::vector<std::complex<double>>& column,
                                                int channels, std::size_t owned, double weight)
{
    if (channels <= 0 || column.size() != static_cast<std::size_t>(channels)
        || partial.size() % static_cast<std::size_t>(channels) != 0
        || owned >= partial.size() / static_cast<std::size_t>(channels)
        || !std::isfinite(weight) || weight <= 0)
        throw std::invalid_argument("Weak q unit invalid partial-column layout or weight.");
    const auto finite = [](std::complex<double> value) {
        return std::isfinite(value.real()) && std::isfinite(value.imag());
    };
    for (int row = 0; row < channels; ++row)
    {
        const auto index = static_cast<std::size_t>(row) + static_cast<std::size_t>(channels) * owned;
        if (!finite(column[row]) || !finite(partial[index]) || !finite(partial[index] + weight * column[row]))
            throw std::invalid_argument("Weak q unit nonfinite weighted response.");
    }
    for (int row = 0; row < channels; ++row)
        partial[static_cast<std::size_t>(row) + static_cast<std::size_t>(channels) * owned] += weight * column[row];
}

} // namespace ModuleRI
#endif
