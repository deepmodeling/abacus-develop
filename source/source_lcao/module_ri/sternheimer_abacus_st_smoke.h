#ifndef STERNHEIMER_ABACUS_ST_SMOKE_H
#define STERNHEIMER_ABACUS_ST_SMOKE_H

#include "source_lcao/module_ri/sternheimer_abacus_fd_adapter.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
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
    int spin_index = -1;
    SternheimerReducedKPoint kpoint{0.0, 0.0, 0.0};
    double kweight = 0.0;
    std::vector<double> eigenvalues;
    std::vector<double> occupations;
    std::vector<std::vector<std::complex<double>>> coefficients;
};

inline void validate_sternheimer_lcao_occupied_kpoints(
    const std::vector<SternheimerLCAOOccupiedKPoint>& records,
    const int local_kpoint_count,
    const int global_kpoint_count,
    const int spin_channel_count,
    const int basis_size)
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
    }
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
    const std::string& output_dir);

} // namespace ModuleRI

#endif
