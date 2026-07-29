#ifndef STERNHEIMER_ABACUS_ST_SMOKE_H
#define STERNHEIMER_ABACUS_ST_SMOKE_H

#include "source_lcao/module_ri/sternheimer_abfs_perturbation.h"
#include "source_lcao/module_ri/sternheimer_abacus_fd_adapter.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

class UnitCell;
class LCAO_Orbitals;

namespace ModulePW
{
class PW_Basis;
}

namespace elecstate
{
class ElecState;
class Potential;
} // namespace elecstate

namespace ModuleRI
{

struct SternheimerLCAOOccupiedKPoint
{
    int local_k_index = -1;
    int global_k_index = -1;
    int zero_order_k_index = -1;
    int symmetry_spatial_isym = 0;
    bool symmetry_time_reversal = false;
    int spin_index = -1;
    SternheimerReducedKPoint kpoint{0.0, 0.0, 0.0};
    double kweight = 0.0;
    std::vector<double> eigenvalues;
    std::vector<double> occupations;
    std::vector<std::vector<std::complex<double>>> coefficients;
    std::vector<double> unoccupied_eigenvalues;
    std::vector<std::vector<std::complex<double>>> unoccupied_coefficients;
};

inline SternheimerLCAOOccupiedKPoint make_sternheimer_full_kpoint_record(
    const SternheimerLCAOOccupiedKPoint& ibz_record,
    const int full_k_index,
    const SternheimerReducedKPoint& kpoint,
    const double full_kweight,
    std::vector<std::vector<std::complex<double>>> occupied_coefficients,
    std::vector<std::vector<std::complex<double>>> unoccupied_coefficients = {})
{
    if (full_k_index < 0 || !std::isfinite(full_kweight) || full_kweight <= 0.0
        || occupied_coefficients.size() != ibz_record.coefficients.size()
        || (!ibz_record.unoccupied_coefficients.empty()
            && unoccupied_coefficients.size() != ibz_record.unoccupied_coefficients.size()))
    {
        throw std::invalid_argument("Invalid full-grid Sternheimer LCAO record data.");
    }
    SternheimerLCAOOccupiedKPoint record = ibz_record;
    record.local_k_index = full_k_index;
    record.global_k_index = full_k_index;
    record.kpoint = kpoint;
    record.kweight = full_kweight;
    record.coefficients = std::move(occupied_coefficients);
    record.unoccupied_coefficients = std::move(unoccupied_coefficients);
    return record;
}

inline std::vector<int> sternheimer_canonical_q_indices_one_based(
    const std::vector<SternheimerLCAOOccupiedKPoint>& records)
{
    if (records.empty())
    {
        throw std::invalid_argument("Cannot identify canonical Sternheimer q points on an empty grid.");
    }
    int zero_order_count = 0;
    for (const auto& record: records)
    {
        if (record.global_k_index < 0 || record.zero_order_k_index < 0)
        {
            throw std::invalid_argument("Sternheimer q-point metadata contain a negative index.");
        }
        zero_order_count = std::max(zero_order_count, record.zero_order_k_index + 1);
    }
    std::vector<int> representatives(static_cast<std::size_t>(zero_order_count), -1);
    for (const auto& record: records)
    {
        if (record.symmetry_spatial_isym != 0 || record.symmetry_time_reversal)
        {
            continue;
        }
        int& representative = representatives[static_cast<std::size_t>(record.zero_order_k_index)];
        if (representative >= 0)
        {
            throw std::invalid_argument("A Sternheimer q star has multiple canonical full-q points.");
        }
        representative = record.global_k_index + 1;
    }
    if (std::any_of(representatives.begin(), representatives.end(), [](const int index) { return index <= 0; }))
    {
        throw std::invalid_argument("A Sternheimer q star has no canonical full-q point.");
    }
    return representatives;
}

inline bool sternheimer_lcao_sos_diagnostic_enabled()
{
    const char* raw = std::getenv("ABACUS_STERNHEIMER_LCAO_SOS_DIAG");
    if (raw == nullptr)
    {
        return false;
    }
    std::string value(raw);
    std::transform(value.begin(), value.end(), value.begin(), [](const unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return !(value.empty() || value == "0" || value == "false" || value == "off" || value == "no");
}

struct SternheimerPeriodicResponsePlan
{
    int iq = 1;
    SternheimerReducedKPoint qpoint{0.0, 0.0, 0.0};
    std::vector<int> record_index_by_global_k;
    std::vector<SternheimerKQPair> kq_pairs;
    double kweight_sum = 0.0;
};

inline int sternheimer_periodic_band_count(const int available_bands, const int requested_bands)
{
    if (available_bands <= 0)
    {
        throw std::invalid_argument("Periodic Sternheimer requires at least one occupied band.");
    }
    return requested_bands > 0 ? std::min(available_bands, requested_bands) : available_bands;
}

inline bool sternheimer_write_periodic_v1(const bool use_supercell_translation_sum,
                                           const bool bands_are_truncated)
{
    return !use_supercell_translation_sum && !bands_are_truncated;
}

inline void validate_sternheimer_periodic_output_mode(const bool write_periodic_v1,
                                                       const bool write_partial_kresolved)
{
    if (!write_periodic_v1 && write_partial_kresolved)
    {
        throw std::invalid_argument(
            "Diagnostic-only periodic Sternheimer output is incompatible with symmetry or k-resolved partial v1.");
    }
}

struct SternheimerFixedQKOrbit
{
    int representative_ik_full = -1;
    std::vector<int> members;
};

struct SternheimerFixedQKRoute
{
    int iq = 0;
    int representative_ik_full = -1;
    int member_ik_full = -1;
    int spatial_isym = -1;
    bool time_reversal = false;
    std::array<int, 3> fold_G{0, 0, 0};
};

struct SternheimerQStarPermutation
{
    int spatial_isym = -1;
    bool time_reversal = false;
    std::vector<int> mapped_index_by_full_q;
    std::vector<std::array<int, 3>> fold_G_by_full_q;
};

struct SternheimerQStarRoute
{
    int representative_iq = 0;
    int member_iq = 0;
    int spatial_isym = -1;
    bool time_reversal = false;
    std::array<int, 3> fold_G{0, 0, 0};
};

inline std::string format_sternheimer_fixed_q_routes(
    const std::vector<SternheimerFixedQKRoute>& routes)
{
    if (routes.empty())
    {
        throw std::invalid_argument("Sternheimer fixed-q route manifest cannot be empty.");
    }
    std::vector<SternheimerFixedQKRoute> ordered = routes;
    std::sort(ordered.begin(), ordered.end(), [](const auto& lhs, const auto& rhs) {
        return std::tie(lhs.iq, lhs.member_ik_full, lhs.representative_ik_full)
               < std::tie(rhs.iq, rhs.member_ik_full, rhs.representative_ik_full);
    });
    std::ostringstream output;
    output << "version 1\n";
    output << "# iq representative_ik member_ik spatial_isym time_reversal fold_Gx fold_Gy fold_Gz\n";
    for (std::size_t index = 0; index != ordered.size(); ++index)
    {
        const auto& route = ordered[index];
        if (route.iq <= 0 || route.representative_ik_full < 0 || route.member_ik_full < 0
            || route.spatial_isym < 0)
        {
            throw std::invalid_argument("Invalid Sternheimer fixed-q route record.");
        }
        if (index > 0
            && std::tie(route.iq, route.member_ik_full)
                   == std::tie(ordered[index - 1].iq, ordered[index - 1].member_ik_full))
        {
            throw std::invalid_argument("Duplicate Sternheimer fixed-q route member.");
        }
        output << route.iq << ' ' << route.representative_ik_full << ' '
               << route.member_ik_full << ' ' << route.spatial_isym << ' '
               << static_cast<int>(route.time_reversal) << ' ' << route.fold_G[0] << ' '
               << route.fold_G[1] << ' ' << route.fold_G[2] << '\n';
    }
    return output.str();
}

inline std::vector<SternheimerFixedQKOrbit> build_sternheimer_fixed_q_k_orbits_from_permutations(
    const int full_kpoint_count,
    const std::vector<std::vector<int>>& little_group_permutations)
{
    if (full_kpoint_count <= 0 || little_group_permutations.empty())
    {
        throw std::invalid_argument("Fixed-q Sternheimer k orbits require a nonempty full k grid and little group.");
    }
    for (const auto& permutation: little_group_permutations)
    {
        if (permutation.size() != static_cast<std::size_t>(full_kpoint_count))
        {
            throw std::invalid_argument("Fixed-q Sternheimer little-group permutation has an invalid size.");
        }
        std::vector<bool> seen(static_cast<std::size_t>(full_kpoint_count), false);
        for (const int mapped_index: permutation)
        {
            if (mapped_index < 0 || mapped_index >= full_kpoint_count
                || seen[static_cast<std::size_t>(mapped_index)])
            {
                throw std::invalid_argument("Fixed-q Sternheimer little-group operation is not a permutation.");
            }
            seen[static_cast<std::size_t>(mapped_index)] = true;
        }
    }

    std::vector<bool> assigned(static_cast<std::size_t>(full_kpoint_count), false);
    std::vector<SternheimerFixedQKOrbit> orbits;
    for (int seed = 0; seed != full_kpoint_count; ++seed)
    {
        if (assigned[static_cast<std::size_t>(seed)])
        {
            continue;
        }
        std::vector<bool> in_orbit(static_cast<std::size_t>(full_kpoint_count), false);
        std::vector<int> pending = {seed};
        in_orbit[static_cast<std::size_t>(seed)] = true;
        for (std::size_t pending_index = 0; pending_index != pending.size(); ++pending_index)
        {
            const int member = pending[pending_index];
            for (const auto& permutation: little_group_permutations)
            {
                const int mapped = permutation[static_cast<std::size_t>(member)];
                if (!in_orbit[static_cast<std::size_t>(mapped)])
                {
                    in_orbit[static_cast<std::size_t>(mapped)] = true;
                    pending.push_back(mapped);
                }
            }
        }
        std::sort(pending.begin(), pending.end());
        SternheimerFixedQKOrbit orbit;
        orbit.representative_ik_full = pending.front();
        orbit.members = std::move(pending);
        for (const int member: orbit.members)
        {
            if (assigned[static_cast<std::size_t>(member)])
            {
                throw std::invalid_argument("Fixed-q Sternheimer k orbits overlap.");
            }
            assigned[static_cast<std::size_t>(member)] = true;
        }
        orbits.push_back(std::move(orbit));
    }
    if (std::any_of(assigned.begin(), assigned.end(), [](const bool value) { return !value; }))
    {
        throw std::invalid_argument("Fixed-q Sternheimer k orbits do not cover the full k grid.");
    }
    return orbits;
}

inline std::vector<SternheimerQStarRoute> build_sternheimer_qstar_routes_from_permutations(
    const int full_qpoint_count,
    const std::vector<SternheimerQStarPermutation>& permutations)
{
    if (full_qpoint_count <= 0 || permutations.empty())
    {
        throw std::invalid_argument(
            "Discrete Sternheimer q-star routes require a nonempty full q grid and symmetry group.");
    }
    std::vector<std::vector<int>> index_permutations;
    index_permutations.reserve(permutations.size());
    for (const auto& permutation: permutations)
    {
        if (permutation.spatial_isym < 0
            || permutation.mapped_index_by_full_q.size()
                   != static_cast<std::size_t>(full_qpoint_count)
            || permutation.fold_G_by_full_q.size()
                   != static_cast<std::size_t>(full_qpoint_count))
        {
            throw std::invalid_argument("Invalid discrete Sternheimer q-star permutation.");
        }
        index_permutations.push_back(permutation.mapped_index_by_full_q);
    }

    const auto orbits = build_sternheimer_fixed_q_k_orbits_from_permutations(
        full_qpoint_count, index_permutations);
    std::vector<SternheimerQStarRoute> routes;
    routes.reserve(static_cast<std::size_t>(full_qpoint_count));
    for (const auto& orbit: orbits)
    {
        for (const int member: orbit.members)
        {
            const auto inverse = std::find_if(
                permutations.begin(), permutations.end(), [&](const auto& permutation) {
                    return permutation.mapped_index_by_full_q[static_cast<std::size_t>(member)]
                           == orbit.representative_ik_full;
                });
            if (inverse == permutations.end())
            {
                throw std::invalid_argument(
                    "A discrete Sternheimer q-star member has no inverse route to its representative.");
            }
            routes.push_back({orbit.representative_ik_full + 1,
                              member + 1,
                              inverse->spatial_isym,
                              inverse->time_reversal,
                              inverse->fold_G_by_full_q[static_cast<std::size_t>(member)]});
        }
    }
    std::sort(routes.begin(), routes.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.member_iq < rhs.member_iq;
    });
    return routes;
}

inline std::string format_sternheimer_qstar_routes(
    const std::vector<SternheimerQStarRoute>& routes)
{
    if (routes.empty())
    {
        throw std::invalid_argument("Sternheimer q-star route manifest cannot be empty.");
    }
    std::vector<SternheimerQStarRoute> ordered = routes;
    std::sort(ordered.begin(), ordered.end(), [](const auto& lhs, const auto& rhs) {
        return std::tie(lhs.member_iq, lhs.representative_iq)
               < std::tie(rhs.member_iq, rhs.representative_iq);
    });
    std::ostringstream output;
    output << "version 1\n";
    output << "# representative_iq member_iq spatial_isym time_reversal fold_Gx fold_Gy fold_Gz\n";
    for (std::size_t index = 0; index != ordered.size(); ++index)
    {
        const auto& route = ordered[index];
        if (route.representative_iq <= 0 || route.member_iq <= 0 || route.spatial_isym < 0)
        {
            throw std::invalid_argument("Invalid Sternheimer q-star route record.");
        }
        if (index > 0 && route.member_iq == ordered[index - 1].member_iq)
        {
            throw std::invalid_argument("Duplicate Sternheimer q-star route member.");
        }
        output << route.representative_iq << ' ' << route.member_iq << ' '
               << route.spatial_isym << ' ' << static_cast<int>(route.time_reversal) << ' '
               << route.fold_G[0] << ' ' << route.fold_G[1] << ' ' << route.fold_G[2] << '\n';
    }
    return output.str();
}

struct SternheimerPartialResponseRecord
{
    int iq = 0;
    int ik_full = -1;
    int ifrequency = 0;
    std::string filename;
    std::vector<std::complex<double>> matrix;
};

inline std::string sternheimer_partial_response_filename(const int iq,
                                                         const int ik_full,
                                                         const int ifrequency)
{
    if (iq <= 0 || ik_full < 0 || ifrequency <= 0)
    {
        throw std::invalid_argument("Invalid Sternheimer partial-response index.");
    }
    std::ostringstream filename;
    filename << "v1_sternheimer_chi0_iq_" << iq << "_ik_" << ik_full << "_ifreq_"
             << ifrequency << ".dat";
    return filename.str();
}

inline SternheimerPartialResponseRecord make_sternheimer_partial_response_record(
    const int iq,
    const int ik_full,
    const int ifrequency,
    const std::vector<std::complex<double>>& branch_matrix,
    const int num_channels)
{
    if (num_channels <= 0
        || branch_matrix.size()
               != static_cast<std::size_t>(num_channels) * static_cast<std::size_t>(num_channels))
    {
        throw std::invalid_argument("Invalid Sternheimer partial-response matrix dimensions.");
    }

    SternheimerPartialResponseRecord record;
    record.iq = iq;
    record.ik_full = ik_full;
    record.ifrequency = ifrequency;
    record.filename = sternheimer_partial_response_filename(iq, ik_full, ifrequency);
    record.matrix.assign(branch_matrix.size(), std::complex<double>(0.0, 0.0));
    for (int row = 0; row != num_channels; ++row)
    {
        for (int column = 0; column != num_channels; ++column)
        {
            const std::size_t index
                = static_cast<std::size_t>(row) * static_cast<std::size_t>(num_channels)
                  + static_cast<std::size_t>(column);
            const std::size_t transpose
                = static_cast<std::size_t>(column) * static_cast<std::size_t>(num_channels)
                  + static_cast<std::size_t>(row);
            record.matrix[index] = branch_matrix[index] + std::conj(branch_matrix[transpose]);
        }
    }
    return record;
}

inline std::string format_sternheimer_partial_manifest(
    const std::vector<SternheimerPartialResponseRecord>& records)
{
    if (records.empty())
    {
        throw std::invalid_argument("Sternheimer partial-response manifest cannot be empty.");
    }
    std::vector<SternheimerPartialResponseRecord> ordered = records;
    std::sort(ordered.begin(), ordered.end(), [](const auto& lhs, const auto& rhs) {
        return std::tie(lhs.iq, lhs.ik_full, lhs.ifrequency)
               < std::tie(rhs.iq, rhs.ik_full, rhs.ifrequency);
    });

    std::ostringstream manifest;
    manifest << "# iq ik_full ifreq response_file\n";
    for (std::size_t index = 0; index != ordered.size(); ++index)
    {
        const auto& record = ordered[index];
        if (record.iq <= 0 || record.ik_full < 0 || record.ifrequency <= 0
            || record.filename.empty()
            || record.filename.find_first_of(" \t\r\n") != std::string::npos)
        {
            throw std::invalid_argument("Invalid Sternheimer partial-response manifest record.");
        }
        if (index > 0
            && std::tie(record.iq, record.ik_full, record.ifrequency)
                   == std::tie(ordered[index - 1].iq,
                               ordered[index - 1].ik_full,
                               ordered[index - 1].ifrequency))
        {
            throw std::invalid_argument("Duplicate Sternheimer partial-response manifest key.");
        }
        manifest << record.iq << ' ' << record.ik_full << ' ' << record.ifrequency << ' '
                 << record.filename << '\n';
    }
    return manifest.str();
}

inline std::string format_sternheimer_full_kpoint_manifest(
    const std::vector<SternheimerLCAOOccupiedKPoint>& records)
{
    if (records.empty())
    {
        throw std::invalid_argument("Sternheimer full-k-point manifest cannot be empty.");
    }
    std::vector<const SternheimerLCAOOccupiedKPoint*> ordered(records.size(), nullptr);
    for (const auto& record: records)
    {
        if (record.global_k_index < 0
            || record.global_k_index >= static_cast<int>(records.size())
            || !std::isfinite(record.kpoint[0]) || !std::isfinite(record.kpoint[1])
            || !std::isfinite(record.kpoint[2])
            || ordered[static_cast<std::size_t>(record.global_k_index)] != nullptr)
        {
            throw std::invalid_argument("Invalid Sternheimer full-k-point manifest record.");
        }
        ordered[static_cast<std::size_t>(record.global_k_index)] = &record;
    }

    std::ostringstream manifest;
    manifest << std::setprecision(17) << "# ik_full kx ky kz\n";
    for (std::size_t ik = 0; ik != ordered.size(); ++ik)
    {
        if (ordered[ik] == nullptr)
        {
            throw std::invalid_argument("Sternheimer full-k-point manifest indices are not contiguous.");
        }
        const auto& kpoint = ordered[ik]->kpoint;
        manifest << ik << ' ' << kpoint[0] << ' ' << kpoint[1] << ' ' << kpoint[2] << '\n';
    }
    return manifest.str();
}

inline std::vector<SternheimerABFBlochGridChannel> limit_sternheimer_abf_channels_per_atom(
    const std::vector<SternheimerABFBlochGridChannel>& channels,
    const int max_channels_per_atom)
{
    if (max_channels_per_atom <= 0)
    {
        return channels;
    }

    int max_atom_index = -1;
    for (const SternheimerABFBlochGridChannel& channel: channels)
    {
        if (channel.atom_index < 0)
        {
            throw std::invalid_argument("A Sternheimer ABFS channel has an invalid atom index.");
        }
        max_atom_index = std::max(max_atom_index, channel.atom_index);
    }

    std::vector<int> selected_per_atom(static_cast<std::size_t>(max_atom_index + 1), 0);
    std::vector<SternheimerABFBlochGridChannel> limited;
    limited.reserve(channels.size());
    for (const SternheimerABFBlochGridChannel& channel: channels)
    {
        int& atom_count = selected_per_atom[static_cast<std::size_t>(channel.atom_index)];
        if (atom_count >= max_channels_per_atom)
        {
            continue;
        }
        SternheimerABFBlochGridChannel selected = channel;
        selected.channel_index = static_cast<int>(limited.size());
        selected.atom_local_index = atom_count;
        ++atom_count;
        limited.push_back(std::move(selected));
    }
    return limited;
}

inline SternheimerPeriodicResponsePlan build_sternheimer_periodic_response_plan(
    const std::vector<SternheimerLCAOOccupiedKPoint>& records,
    const int q_index,
    const bool single_gamma_supercell_translation = false)
{
    if (records.empty())
    {
        throw std::invalid_argument("Sternheimer response plan requires occupied k-point records.");
    }
    if (q_index < 0 || q_index > static_cast<int>(records.size()))
    {
        throw std::invalid_argument("sternheimer_q_index is outside the full k-point mesh.");
    }

    SternheimerPeriodicResponsePlan plan;
    plan.record_index_by_global_k.assign(records.size(), -1);
    std::vector<SternheimerReducedKPoint> kpoints(records.size());
    for (std::size_t record_index = 0; record_index != records.size(); ++record_index)
    {
        const SternheimerLCAOOccupiedKPoint& record = records[record_index];
        if (record.global_k_index < 0 || record.global_k_index >= static_cast<int>(records.size()))
        {
            throw std::invalid_argument("Sternheimer response plan found an invalid global k-point index.");
        }
        int& mapped_record = plan.record_index_by_global_k[static_cast<std::size_t>(record.global_k_index)];
        if (mapped_record >= 0)
        {
            throw std::invalid_argument("Sternheimer response plan found a duplicate global k-point index.");
        }
        mapped_record = static_cast<int>(record_index);
        kpoints[static_cast<std::size_t>(record.global_k_index)] = record.kpoint;
        plan.kweight_sum += record.kweight;
    }

    constexpr double tolerance = 1.0e-10;
    if (q_index == 0 || single_gamma_supercell_translation)
    {
        if (single_gamma_supercell_translation && q_index != 1)
        {
            throw std::invalid_argument(
                "A single-Gamma supercell translation response must use output q index 1.");
        }
        if (records.size() != 1
            || std::any_of(records.front().kpoint.begin(),
                           records.front().kpoint.end(),
                           [](const double coordinate) { return std::abs(coordinate) > tolerance; }))
        {
            throw std::invalid_argument(
                "sternheimer_q_index=0 is reserved for the single-k Gamma compatibility path.");
        }
    }
    else
    {
        if (records.size() <= 1)
        {
            throw std::invalid_argument("A nonzero Sternheimer q point requires more than one k point.");
        }
        for (const SternheimerLCAOOccupiedKPoint& record: records)
        {
            if (record.spin_index != 0)
            {
                throw std::invalid_argument("The first solid Sternheimer driver supports only nspin=1.");
            }
            for (const double occupation: record.occupations)
            {
                if (std::abs(occupation - 1.0) > tolerance)
                {
                    throw std::invalid_argument(
                        "The first nspin=1 solid Sternheimer driver requires fully occupied insulating bands.");
                }
            }
        }
        plan.iq = q_index;
        const int q_record_index
            = plan.record_index_by_global_k[static_cast<std::size_t>(q_index - 1)];
        plan.qpoint = records[static_cast<std::size_t>(q_record_index)].kpoint;
    }
    plan.kq_pairs = build_sternheimer_kq_map(kpoints, plan.qpoint, tolerance);
    return plan;
}

inline void validate_sternheimer_periodic_kmesh(const std::array<int, 3>& kmesh,
                                                const int global_kpoint_count)
{
    if (global_kpoint_count <= 0
        || std::any_of(kmesh.begin(), kmesh.end(), [](const int dimension) { return dimension <= 0; }))
    {
        throw std::invalid_argument("Periodic Sternheimer requires positive Monkhorst-Pack dimensions.");
    }
    const long long mesh_size = static_cast<long long>(kmesh[0]) * kmesh[1] * kmesh[2];
    if (mesh_size != global_kpoint_count)
    {
        throw std::invalid_argument(
            "Periodic Sternheimer Monkhorst-Pack dimensions do not match the full k-point count.");
    }
}

inline double sternheimer_periodic_gamma_inverse_k2(const SternheimerReducedKPoint& qpoint,
                                                     const std::string& singularity_correction,
                                                     const double massidda_chi)
{
    constexpr double tolerance = 1.0e-10;
    const bool gamma = std::all_of(qpoint.begin(), qpoint.end(), [](const double coordinate) {
        return std::abs(coordinate) <= tolerance;
    });
    if (!gamma)
    {
        return 0.0;
    }
    if (singularity_correction != "massidda")
    {
        throw std::invalid_argument("Periodic Sternheimer q=0 requires Massidda singularity correction.");
    }
    if (!std::isfinite(massidda_chi) || massidda_chi <= 0.0)
    {
        throw std::invalid_argument("Periodic Sternheimer q=0 requires a positive finite Massidda value.");
    }
    return massidda_chi;
}

inline int sternheimer_kpoint_owner_group(const int global_kpoint_index,
                                          const int global_kpoint_count,
                                          const int kpoint_groups)
{
    if (global_kpoint_count <= 0 || kpoint_groups <= 0 || kpoint_groups > global_kpoint_count
        || global_kpoint_index < 0 || global_kpoint_index >= global_kpoint_count)
    {
        throw std::invalid_argument("Invalid Sternheimer k-point partition dimensions.");
    }

    const int base_count = global_kpoint_count / kpoint_groups;
    const int extra_groups = global_kpoint_count % kpoint_groups;
    const int enlarged_span = extra_groups * (base_count + 1);
    if (global_kpoint_index < enlarged_span)
    {
        return global_kpoint_index / (base_count + 1);
    }
    return extra_groups + (global_kpoint_index - enlarged_span) / base_count;
}

struct SternheimerNestedMPIAssignment
{
    int kpoint_group = 0;
    int frequency_slot = 0;
    int owner_rank = 0;
};

inline SternheimerNestedMPIAssignment sternheimer_nested_mpi_assignment(
    const int global_kpoint_index,
    const int global_kpoint_count,
    const int ifrequency_zero_based,
    const int frequency_count,
    const int kpoint_groups,
    const int mpi_ranks,
    const int frequency_rank_shift = 0)
{
    if (frequency_count <= 0 || ifrequency_zero_based < 0
        || ifrequency_zero_based >= frequency_count)
    {
        throw std::invalid_argument("Invalid Sternheimer nested-MPI frequency dimensions.");
    }
    if (mpi_ranks != kpoint_groups * frequency_count)
    {
        throw std::invalid_argument(
            "Nested Sternheimer MPI requires NPROC=k-point-groups*frequency-count.");
    }

    const int kpoint_group = sternheimer_kpoint_owner_group(global_kpoint_index,
                                                             global_kpoint_count,
                                                             kpoint_groups);
    int normalized_shift = frequency_rank_shift % frequency_count;
    if (normalized_shift < 0)
    {
        normalized_shift += frequency_count;
    }
    const int frequency_slot = (ifrequency_zero_based + normalized_shift) % frequency_count;
    return {kpoint_group,
            frequency_slot,
            kpoint_group * frequency_count + frequency_slot};
}

inline std::vector<std::size_t> sternheimer_owned_kq_pair_indices(const SternheimerPeriodicResponsePlan& plan,
                                                                  const int kpoint_group,
                                                                  const int kpoint_groups)
{
    if (kpoint_group < 0 || kpoint_group >= kpoint_groups)
    {
        throw std::invalid_argument("Invalid Sternheimer k-point group index.");
    }
    const int global_kpoint_count = static_cast<int>(plan.record_index_by_global_k.size());
    std::vector<std::size_t> owned;
    for (std::size_t pair_index = 0; pair_index != plan.kq_pairs.size(); ++pair_index)
    {
        if (sternheimer_kpoint_owner_group(plan.kq_pairs[pair_index].source_index,
                                           global_kpoint_count,
                                           kpoint_groups)
            == kpoint_group)
        {
            owned.push_back(pair_index);
        }
    }
    return owned;
}

inline void validate_sternheimer_lcao_occupied_kpoints(
    const std::vector<SternheimerLCAOOccupiedKPoint>& records,
    const int local_kpoint_count,
    const int global_kpoint_count,
    const int spin_channel_count,
    const int basis_size,
    const int zero_order_kpoint_count = -1)
{
    if (local_kpoint_count <= 0 || global_kpoint_count <= 0 || spin_channel_count <= 0 || basis_size <= 0)
    {
        throw std::invalid_argument("Sternheimer LCAO k-point dimensions must be positive.");
    }
    if (records.size() != static_cast<std::size_t>(global_kpoint_count)
        || local_kpoint_count != global_kpoint_count)
    {
        throw std::invalid_argument(
            "Sternheimer LCAO occupied k-point records are incomplete; the first solid implementation requires "
            "KPAR=1.");
    }

    std::vector<bool> seen_local(static_cast<std::size_t>(local_kpoint_count), false);
    std::vector<bool> seen_global(static_cast<std::size_t>(global_kpoint_count), false);
    const int zero_order_count
        = zero_order_kpoint_count > 0 ? zero_order_kpoint_count : local_kpoint_count;
    for (const SternheimerLCAOOccupiedKPoint& record: records)
    {
        if (record.local_k_index < 0 || record.local_k_index >= local_kpoint_count)
        {
            throw std::invalid_argument("Sternheimer LCAO local k-point index is out of range.");
        }
        if (record.global_k_index < 0 || record.global_k_index >= global_kpoint_count)
        {
            throw std::invalid_argument("Sternheimer LCAO global k-point index is out of range.");
        }
        if (record.zero_order_k_index < 0 || record.zero_order_k_index >= zero_order_count)
        {
            throw std::invalid_argument("Sternheimer LCAO zero-order k-point index is out of range.");
        }
        if (seen_local[static_cast<std::size_t>(record.local_k_index)])
        {
            throw std::invalid_argument("Sternheimer LCAO local k-point index is duplicated.");
        }
        if (seen_global[static_cast<std::size_t>(record.global_k_index)])
        {
            throw std::invalid_argument("Sternheimer LCAO global k-point index is duplicated.");
        }
        seen_local[static_cast<std::size_t>(record.local_k_index)] = true;
        seen_global[static_cast<std::size_t>(record.global_k_index)] = true;

        if (record.spin_index < 0 || record.spin_index >= spin_channel_count)
        {
            throw std::invalid_argument("Sternheimer LCAO occupied spin index is out of range.");
        }
        for (const double coordinate: record.kpoint)
        {
            if (!std::isfinite(coordinate))
            {
                throw std::invalid_argument("Sternheimer LCAO reduced k-point coordinate is not finite.");
            }
        }
        if (!std::isfinite(record.kweight) || record.kweight <= 0.0)
        {
            throw std::invalid_argument("Sternheimer LCAO k-point weight must be finite and positive.");
        }
        if (record.coefficients.empty()
            || record.eigenvalues.size() != record.coefficients.size()
            || record.occupations.size() != record.coefficients.size())
        {
            throw std::invalid_argument("Sternheimer LCAO occupied k-point band data are inconsistent.");
        }
        for (std::size_t ib = 0; ib != record.coefficients.size(); ++ib)
        {
            if (!std::isfinite(record.eigenvalues[ib])
                || !std::isfinite(record.occupations[ib])
                || record.occupations[ib] <= 0.0)
            {
                throw std::invalid_argument("Sternheimer LCAO occupied eigenvalue or occupation is invalid.");
            }
            const auto& band_coefficients = record.coefficients[ib];
            if (band_coefficients.size() != static_cast<std::size_t>(basis_size))
            {
                throw std::invalid_argument("Sternheimer LCAO coefficient basis size is inconsistent.");
            }
            for (const std::complex<double>& coefficient: band_coefficients)
            {
                if (!std::isfinite(coefficient.real()) || !std::isfinite(coefficient.imag()))
                {
                    throw std::invalid_argument("Sternheimer LCAO coefficient is not finite.");
                }
            }
        }
        if (record.unoccupied_eigenvalues.size() != record.unoccupied_coefficients.size())
        {
            throw std::invalid_argument("Sternheimer LCAO unoccupied k-point band data are inconsistent.");
        }
        for (std::size_t ib = 0; ib != record.unoccupied_coefficients.size(); ++ib)
        {
            if (!std::isfinite(record.unoccupied_eigenvalues[ib]))
            {
                throw std::invalid_argument("Sternheimer LCAO unoccupied eigenvalue is invalid.");
            }
            const auto& band_coefficients = record.unoccupied_coefficients[ib];
            if (band_coefficients.size() != static_cast<std::size_t>(basis_size))
            {
                throw std::invalid_argument("Sternheimer LCAO unoccupied coefficient basis size is inconsistent.");
            }
            for (const std::complex<double>& coefficient: band_coefficients)
            {
                if (!std::isfinite(coefficient.real()) || !std::isfinite(coefficient.imag()))
                {
                    throw std::invalid_argument("Sternheimer LCAO unoccupied coefficient is not finite.");
                }
            }
        }
    }
}

inline void validate_sternheimer_full_lcao_occupied_kpoints(
    const std::vector<SternheimerLCAOOccupiedKPoint>& records,
    const int zero_order_kpoint_count,
    const int spin_channel_count,
    const int basis_size)
{
    if (zero_order_kpoint_count <= 0)
    {
        throw std::invalid_argument("Sternheimer LCAO zero-order k-point count must be positive.");
    }
    const int full_kpoint_count = static_cast<int>(records.size());
    validate_sternheimer_lcao_occupied_kpoints(records,
                                               full_kpoint_count,
                                               full_kpoint_count,
                                               spin_channel_count,
                                               basis_size,
                                               zero_order_kpoint_count);
}

inline int sternheimer_lcao_total_occupied_bands(
    const std::vector<SternheimerLCAOOccupiedKPoint>& records)
{
    int count = 0;
    for (const SternheimerLCAOOccupiedKPoint& record: records)
    {
        count += static_cast<int>(record.coefficients.size());
    }
    return count;
}

inline double sternheimer_lcao_weighted_occupation(const SternheimerLCAOOccupiedKPoint& record,
                                                   const int band_index)
{
    if (band_index < 0 || band_index >= static_cast<int>(record.occupations.size()))
    {
        throw std::out_of_range("Sternheimer LCAO occupied band index is out of range.");
    }
    return record.kweight * record.occupations[static_cast<std::size_t>(band_index)];
}

inline std::vector<const SternheimerLCAOOccupiedKPoint*> select_sternheimer_gamma_spin_records(
    const std::vector<SternheimerLCAOOccupiedKPoint>& records,
    const int spin_channel_count,
    const double tolerance = 1.0e-12)
{
    if (spin_channel_count <= 0 || records.size() != static_cast<std::size_t>(spin_channel_count))
    {
        throw std::invalid_argument(
            "Gamma Sternheimer response requires exactly one LCAO record per spin channel.");
    }
    std::vector<const SternheimerLCAOOccupiedKPoint*> selected(
        static_cast<std::size_t>(spin_channel_count), nullptr);
    for (const SternheimerLCAOOccupiedKPoint& record: records)
    {
        if (record.spin_index < 0 || record.spin_index >= spin_channel_count)
        {
            throw std::invalid_argument("Gamma Sternheimer response spin index is out of range.");
        }
        if (std::any_of(record.kpoint.begin(), record.kpoint.end(), [tolerance](const double coordinate) {
                return std::abs(coordinate) > tolerance;
            }))
        {
            throw std::invalid_argument("Gamma Sternheimer response received a non-Gamma LCAO record.");
        }
        const std::size_t spin_index = static_cast<std::size_t>(record.spin_index);
        if (selected[spin_index] != nullptr)
        {
            throw std::invalid_argument("Gamma Sternheimer response has duplicate spin records.");
        }
        selected[spin_index] = &record;
    }
    if (std::any_of(selected.begin(), selected.end(), [](const auto* record) { return record == nullptr; }))
    {
        throw std::invalid_argument("Gamma Sternheimer response is missing a spin record.");
    }
    return selected;
}

struct SternheimerABACUSSTChannelResult
{
    int band_index = -1;
    int channel_index = -1;
    int atom_index = -1;
    int angular_momentum = 0;
    int radial_index = 0;
    int magnetic_index = 0;
    double fd_eigenvalue = 0.0;
    double occupation = 0.0;
    double rhs_norm = 0.0;
    double projected_rhs_norm = 0.0;
    bool solver_converged = false;
    int solver_iterations = 0;
    double solver_relative_residual = 0.0;
    double equation_residual_norm = 0.0;
    std::complex<double> polarizability = {0.0, 0.0};
};

struct SternheimerABACUSSTSmokeResult
{
    SternheimerABACUSFDGridData grid_data;
    double omega = 0.0;
    double pca_threshold = 0.0;
    double ccp_rmesh_times = 0.0;
    std::string perturbation_source;
    int num_available_channels = 0;
    std::vector<SternheimerABACUSSTChannelResult> channels;
};

inline std::string format_sternheimer_abacus_st_report(const SternheimerABACUSSTSmokeResult& result)
{
    std::ostringstream out;
    const int grid_size = result.grid_data.grid.nx * result.grid_data.grid.ny * result.grid_data.grid.nz;
    out << std::setprecision(16);
    out << "# ABACUS Sternheimer FD linear-response smoke test\n";
    out << "grid " << result.grid_data.grid.nx << ' ' << result.grid_data.grid.ny << ' '
        << result.grid_data.grid.nz << " size " << grid_size << " dV " << result.grid_data.volume_element << '\n';
    out << "omega_Ry " << result.omega << '\n';
    out << "pca_threshold " << result.pca_threshold << '\n';
    out << "ccp_rmesh_times " << result.ccp_rmesh_times << '\n';
    out << "perturbation_source " << result.perturbation_source << '\n';
    out << "available_channels " << result.num_available_channels << '\n';
    out << "band channel atom l radial m fd_eigenvalue_Ry occupation rhs_norm projected_rhs_norm "
           "solver_converged solver_iterations solver_relative_residual equation_residual_norm "
           "polarizability_real polarizability_imag\n";
    for (const SternheimerABACUSSTChannelResult& channel: result.channels)
    {
        out << channel.band_index << ' ' << channel.channel_index << ' ' << channel.atom_index << ' '
            << channel.angular_momentum << ' ' << channel.radial_index << ' ' << channel.magnetic_index << ' '
            << channel.fd_eigenvalue << ' ' << channel.occupation << ' ' << channel.rhs_norm << ' '
            << channel.projected_rhs_norm << ' ' << (channel.solver_converged ? "yes" : "no") << ' '
            << channel.solver_iterations << ' ' << channel.solver_relative_residual << ' '
            << channel.equation_residual_norm << ' ' << channel.polarizability.real() << ' '
            << channel.polarizability.imag() << '\n';
    }
    return out.str();
}

bool sternheimer_abacus_st_smoke_enabled();

inline bool sternheimer_abfs_diag_only_enabled()
{
    const char* raw = std::getenv("ABACUS_STERNHEIMER_FD_ST_ABFS_DIAG_ONLY");
    if (raw == nullptr)
    {
        return false;
    }
    std::string value(raw);
    std::transform(value.begin(), value.end(), value.begin(), [](const unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return !(value.empty() || value == "0" || value == "false" || value == "off" || value == "no");
}

inline bool sternheimer_grid_coulomb_diagnostic_enabled(const int num_channels)
{
    if (num_channels <= 32)
    {
        return true;
    }
    const char* raw = std::getenv("ABACUS_STERNHEIMER_GRID_COULOMB_DIAG");
    if (raw == nullptr)
    {
        return false;
    }
    std::string value(raw);
    std::transform(value.begin(), value.end(), value.begin(), [](const unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return !(value.empty() || value == "0" || value == "false" || value == "off" || value == "no");
}

void run_sternheimer_abacus_st_smoke(const elecstate::Potential& potential,
                                     const ModulePW::PW_Basis& pw_basis,
                                     const UnitCell& ucell,
                                     const elecstate::ElecState& elec_state,
                                     const std::string& output_dir);

void run_sternheimer_abacus_chi0_output(const elecstate::Potential& potential,
                                        const ModulePW::PW_Basis& pw_basis,
                                        const UnitCell& ucell,
                                        const elecstate::ElecState& elec_state,
                                        const std::string& output_dir);

void run_sternheimer_abacus_lcao_chi0_output(
    const elecstate::Potential& potential,
    const ModulePW::PW_Basis& pw_basis,
    const UnitCell& ucell,
    const elecstate::ElecState& elec_state,
    const LCAO_Orbitals& orbitals,
    const std::vector<SternheimerLCAOOccupiedKPoint>& occupied_kpoints,
    const std::array<int, 3>& kmesh,
    const std::string& output_dir);

} // namespace ModuleRI

#endif
