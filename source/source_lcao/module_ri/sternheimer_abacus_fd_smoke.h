#ifndef STERNHEIMER_ABACUS_FD_SMOKE_H
#define STERNHEIMER_ABACUS_FD_SMOKE_H

#include "source_lcao/module_ri/sternheimer_abacus_dft_zero_order.h"

#include <iomanip>
#include <sstream>
#include <string>

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

inline std::string format_sternheimer_fd_zero_order_report(const SternheimerABACUSDFTZeroOrderResult& result)
{
    std::ostringstream out;
    const int grid_size = result.grid_data.grid.nx * result.grid_data.grid.ny * result.grid_data.grid.nz;
    out << std::setprecision(16);
    out << "# ABACUS Sternheimer FD zero-order smoke test\n";
    out << "grid " << result.grid_data.grid.nx << ' ' << result.grid_data.grid.ny << ' '
        << result.grid_data.grid.nz << " size " << grid_size << " dV "
        << result.grid_data.volume_element << '\n';
    out << "spacing " << result.grid_data.grid.hx << ' ' << result.grid_data.grid.hy << ' '
        << result.grid_data.grid.hz << '\n';
    out << "band fd_eigenvalue_Ry dft_eigenvalue_Ry fd_minus_dft_Ry occupation fd_residual_norm "
           "eigenvalue_within_tolerance\n";

    for (const SternheimerDFTZeroOrderBandComparison& band: result.comparison.bands)
    {
        const int band_index = band.band_index;
        const double occupation = band_index >= 0 && band_index < static_cast<int>(result.dft_occupations.size())
                                      ? result.dft_occupations[band_index]
                                      : 0.0;
        out << band.band_index << ' ' << band.fd_eigenvalue << ' ' << band.dft_eigenvalue << ' '
            << band.fd_minus_dft << ' ' << occupation << ' ' << band.fd_residual_norm << ' '
            << (band.eigenvalue_within_tolerance ? "yes" : "no") << '\n';
    }

    out << "max_abs_fd_minus_dft_Ry " << result.comparison.max_abs_fd_minus_dft << '\n';
    out << "all_eigenvalues_within_tolerance "
        << (result.comparison.all_eigenvalues_within_tolerance ? "yes" : "no") << '\n';
    return out.str();
}

bool sternheimer_fd_zero_order_smoke_enabled();

void run_sternheimer_fd_zero_order_smoke(const elecstate::Potential& potential,
                                         const ModulePW::PW_Basis& pw_basis,
                                         const UnitCell& ucell,
                                         const elecstate::ElecState& elec_state,
                                         const std::string& output_dir);

} // namespace ModuleRI

#endif
