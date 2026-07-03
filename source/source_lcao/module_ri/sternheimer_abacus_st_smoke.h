#ifndef STERNHEIMER_ABACUS_ST_SMOKE_H
#define STERNHEIMER_ABACUS_ST_SMOKE_H

#include "source_lcao/module_ri/sternheimer_abacus_fd_adapter.h"

#include <complex>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

class UnitCell;

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

} // namespace ModuleRI

#endif
