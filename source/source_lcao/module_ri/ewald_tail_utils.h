#ifndef EWALD_TAIL_UTILS_H
#define EWALD_TAIL_UTILS_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace EwaldVqDetail
{
using TailCell = std::array<int, 3>;

enum class TailMode
{
    Off,
    Warn,
    Enlarge,
    Strict
};

enum class TailDecision
{
    Accept,
    Warn,
    Expand,
    Fail
};

enum class TailBlockCoverage
{
    Common,
    AnalyticCancellation,
    Incomplete
};

struct TailBlockPlan
{
    TailBlockCoverage coverage = TailBlockCoverage::Incomplete;
    bool use_bare_multipole = false;
    bool use_gaussian_multipole = false;
};

struct TailKey
{
    int atom_i = 0;
    int atom_j = 0;
    TailCell cell{0, 0, 0};
};

inline bool operator==(const TailKey& lhs, const TailKey& rhs)
{
    return lhs.atom_i == rhs.atom_i && lhs.atom_j == rhs.atom_j && lhs.cell == rhs.cell;
}

inline bool operator<(const TailKey& lhs, const TailKey& rhs)
{
    return std::tie(lhs.atom_i, lhs.atom_j, lhs.cell)
           < std::tie(rhs.atom_i, rhs.atom_j, rhs.cell);
}

struct TailTolerances
{
    double abs_tol = 1e-12;
    double rel_tol = 1e-10;
    double sum_rel_tol = 1e-8;
};

struct TailStats
{
    std::size_t production_blocks = 0;
    std::size_t shell_blocks = 0;
    double reference_norm = 0.0;
    double max_norm = 0.0;
    double sum_norm = 0.0;
    double hermitian_residual = 0.0;
    TailKey worst_key{};
    double worst_bare_norm = 0.0;
    double worst_gaussian_norm = 0.0;
    double worst_distance = 0.0;
    bool coverage_complete = true;
};

struct TailSample
{
    TailKey key{};
    double bare_norm = 0.0;
    double gaussian_norm = 0.0;
    double difference_norm = 0.0;
    double distance = 0.0;
};

inline TailMode parse_tail_mode(const std::string& mode)
{
    if (mode == "off") return TailMode::Off;
    if (mode == "warn") return TailMode::Warn;
    if (mode == "enlarge") return TailMode::Enlarge;
    if (mode == "strict") return TailMode::Strict;
    throw std::invalid_argument("Unknown Ewald tail-check mode: " + mode);
}

inline const char* tail_mode_name(const TailMode mode)
{
    switch (mode)
    {
        case TailMode::Off: return "off";
        case TailMode::Warn: return "warn";
        case TailMode::Enlarge: return "enlarge";
        case TailMode::Strict: return "strict";
    }
    return "unknown";
}

inline const char* tail_decision_name(const TailDecision decision)
{
    switch (decision)
    {
        case TailDecision::Accept: return "accepted";
        case TailDecision::Warn: return "warning";
        case TailDecision::Expand: return "expanded";
        case TailDecision::Fail: return "failed";
    }
    return "unknown";
}

inline bool uses_adaptive_realspace_range(const TailMode mode)
{
    return mode == TailMode::Enlarge || mode == TailMode::Strict;
}

inline double adaptive_exact_rmesh_times(const double requested,
                                         const std::vector<double>& lcaos_cutoff,
                                         const std::vector<double>& abfs_cutoff)
{
    if (requested <= 0.0)
    {
        throw std::invalid_argument("Ewald radial-mesh factor must be positive");
    }
    if (lcaos_cutoff.size() != abfs_cutoff.size())
    {
        throw std::invalid_argument("LCAO and ABFS cutoff arrays must have the same size");
    }

    // A factor of one is also the minimum exact Gaussian pair range.  For the
    // bare pair, require c*r_I^ABFS + r_J^LCAO to reach the ABFS nonoverlap
    // boundary r_I^ABFS + r_J^ABFS for every ordered type pair.
    double exact = std::max(requested, 1.0);
    for (std::size_t type_i = 0; type_i != abfs_cutoff.size(); ++type_i)
    {
        if (abfs_cutoff[type_i] <= 0.0) continue;
        for (std::size_t type_j = 0; type_j != abfs_cutoff.size(); ++type_j)
        {
            const double required = (abfs_cutoff[type_i] + abfs_cutoff[type_j]
                                     - lcaos_cutoff[type_j])
                                    / abfs_cutoff[type_i];
            exact = std::max(exact, required);
        }
    }
    return exact;
}

inline void validate_odd_period(const TailCell& period)
{
    for (const int value: period)
    {
        if (value <= 0 || value % 2 == 0)
        {
            throw std::invalid_argument("Ewald real-space periods must be positive odd integers");
        }
    }
}

inline TailCell expand_odd_period(const TailCell& period,
                                  const std::array<bool, 3>& periodic,
                                  const int guard_cells)
{
    validate_odd_period(period);
    if (guard_cells < 0)
    {
        throw std::invalid_argument("Ewald tail guard cells must be non-negative");
    }
    TailCell expanded = period;
    for (std::size_t dim = 0; dim != expanded.size(); ++dim)
    {
        if (periodic[dim]) expanded[dim] += 2 * guard_cells;
    }
    return expanded;
}

inline bool cell_in_period(const TailCell& cell, const TailCell& period)
{
    validate_odd_period(period);
    for (std::size_t dim = 0; dim != cell.size(); ++dim)
    {
        if (cell[dim] < -period[dim] / 2 || cell[dim] > period[dim] / 2) return false;
    }
    return true;
}

inline std::vector<TailCell> make_shell_cells(const TailCell& production_period,
                                              const TailCell& probe_period)
{
    validate_odd_period(production_period);
    validate_odd_period(probe_period);
    for (std::size_t dim = 0; dim != production_period.size(); ++dim)
    {
        if (probe_period[dim] < production_period[dim])
        {
            throw std::invalid_argument("Ewald probe period cannot be smaller than production period");
        }
    }

    std::vector<TailCell> shell;
    for (int x = -probe_period[0] / 2; x <= probe_period[0] / 2; ++x)
    {
        for (int y = -probe_period[1] / 2; y <= probe_period[1] / 2; ++y)
        {
            for (int z = -probe_period[2] / 2; z <= probe_period[2] / 2; ++z)
            {
                const TailCell cell{x, y, z};
                if (!cell_in_period(cell, production_period)) shell.push_back(cell);
            }
        }
    }
    return shell;
}

inline std::vector<TailKey> make_shell_keys(const std::vector<int>& atoms,
                                            const TailCell& production_period,
                                            const TailCell& probe_period)
{
    const auto cells = make_shell_cells(production_period, probe_period);
    std::vector<TailKey> keys;
    keys.reserve(atoms.size() * atoms.size() * cells.size());
    for (const int atom_i: atoms)
    {
        for (const int atom_j: atoms)
        {
            for (const TailCell& cell: cells)
            {
                keys.push_back(TailKey{atom_i, atom_j, cell});
            }
        }
    }
    return keys;
}

inline TailKey hermitian_partner(const TailKey& key)
{
    return TailKey{key.atom_j,
                   key.atom_i,
                   TailCell{-key.cell[0], -key.cell[1], -key.cell[2]}};
}

inline TailKey canonical_hermitian_key(const TailKey& key)
{
    const TailKey partner = hermitian_partner(key);
    return partner < key ? partner : key;
}

inline std::uint64_t stable_tail_key_hash(const TailKey& key)
{
    const TailKey canonical = canonical_hermitian_key(key);
    std::uint64_t hash = 1469598103934665603ULL;
    const auto mix = [&hash](const std::int64_t value)
    {
        const std::uint64_t encoded = static_cast<std::uint64_t>(value);
        for (int byte = 0; byte != 8; ++byte)
        {
            hash ^= (encoded >> (8 * byte)) & 0xffULL;
            hash *= 1099511628211ULL;
        }
    };
    mix(canonical.atom_i);
    mix(canonical.atom_j);
    for (const int value: canonical.cell) mix(value);
    return hash;
}

inline int shell_owner(const TailKey& key, const int mpi_size)
{
    if (mpi_size <= 0) throw std::invalid_argument("MPI size must be positive");
    return static_cast<int>(stable_tail_key_hash(key) % static_cast<std::uint64_t>(mpi_size));
}

inline double tail_reference_scale(const TailStats& stats)
{
    return std::max(stats.reference_norm, std::numeric_limits<double>::min());
}

inline double relative_tail_max(const TailStats& stats)
{
    return stats.max_norm / tail_reference_scale(stats);
}

inline double relative_tail_sum(const TailStats& stats)
{
    return stats.sum_norm / tail_reference_scale(stats);
}

inline bool tail_passes(const TailStats& stats, const TailTolerances& tolerances)
{
    if (!stats.coverage_complete) return false;
    const double reference = tail_reference_scale(stats);
    return stats.max_norm <= tolerances.abs_tol + tolerances.rel_tol * reference
           && stats.sum_norm <= tolerances.sum_rel_tol * reference
           && stats.hermitian_residual <= tolerances.abs_tol + tolerances.rel_tol * reference;
}

inline TailDecision decide_tail(const TailStats& stats,
                                const TailTolerances& tolerances,
                                const TailMode mode,
                                const int expansions,
                                const int max_expansions)
{
    if (mode == TailMode::Off) return TailDecision::Accept;
    const double hermitian_limit
        = tolerances.abs_tol + tolerances.rel_tol * tail_reference_scale(stats);
    if (stats.hermitian_residual > hermitian_limit)
    {
        return TailDecision::Fail;
    }
    if (!stats.coverage_complete)
    {
        if (mode == TailMode::Warn || expansions >= max_expansions)
        {
            return TailDecision::Fail;
        }
        return TailDecision::Expand;
    }
    if (tail_passes(stats, tolerances)) return TailDecision::Accept;
    if (mode == TailMode::Warn) return TailDecision::Warn;
    if (expansions >= max_expansions)
    {
        return TailDecision::Fail;
    }
    return TailDecision::Expand;
}

inline TailBlockCoverage classify_tail_block(const bool has_bare,
                                              const bool has_gaussian,
                                              const bool far_multipoles_equivalent)
{
    if (has_bare && has_gaussian) return TailBlockCoverage::Common;
    if (far_multipoles_equivalent) return TailBlockCoverage::AnalyticCancellation;
    return TailBlockCoverage::Incomplete;
}

inline TailBlockPlan plan_tail_block(const bool has_bare,
                                     const bool has_gaussian,
                                     const double distance,
                                     const double bare_support,
                                     const double gaussian_support)
{
    const bool use_bare_multipole = distance >= bare_support;
    const bool use_gaussian_multipole = distance >= gaussian_support;
    if ((!has_bare && !use_bare_multipole) || (!has_gaussian && !use_gaussian_multipole))
    {
        return TailBlockPlan{};
    }
    if (use_bare_multipole && use_gaussian_multipole)
    {
        return TailBlockPlan{TailBlockCoverage::AnalyticCancellation, false, false};
    }
    return TailBlockPlan{TailBlockCoverage::Common,
                         use_bare_multipole,
                         use_gaussian_multipole};
}

inline void accumulate_tail_sample(TailStats& stats, const TailSample& sample)
{
    ++stats.shell_blocks;
    stats.sum_norm += sample.difference_norm;
    if (sample.difference_norm > stats.max_norm)
    {
        stats.max_norm = sample.difference_norm;
        stats.worst_key = sample.key;
        stats.worst_bare_norm = sample.bare_norm;
        stats.worst_gaussian_norm = sample.gaussian_norm;
        stats.worst_distance = sample.distance;
    }
}

inline std::map<int, std::vector<std::pair<int, TailCell>>>
group_tail_keys_by_atom(const std::vector<TailKey>& keys)
{
    std::map<int, std::vector<std::pair<int, TailCell>>> grouped;
    for (const TailKey& key: keys)
    {
        grouped[key.atom_i].emplace_back(key.atom_j, key.cell);
    }
    return grouped;
}

inline bool far_multipoles_equivalent(const double distance,
                                      const double bare_support_radius,
                                      const double gaussian_support_radius)
{
    return distance >= bare_support_radius && distance >= gaussian_support_radius;
}
} // namespace EwaldVqDetail

#endif
