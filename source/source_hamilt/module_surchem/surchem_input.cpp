#include "surchem_input.h"
#include "surchem.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_cell/klist.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_quit.h"

#include <algorithm>
#include <cmath>

namespace ModuleSurchem
{
SurchemParameters make_parameters(const Input_para& inp,
                                  const UnitCell& ucell,
                                  const double electron_count,
                                  const bool use_uspp,
                                  const int pool_process_count)
{
    SurchemParameters parameters;
    parameters.eb_k = inp.eb_k;
    parameters.tau = inp.tau;
    parameters.sigma_k = inp.sigma_k;
    parameters.nc_k = inp.nc_k;
    parameters.use_sccs = inp.imp_sol == 2;
    parameters.use_legacy_solvent = inp.imp_sol == 1;
    if (inp.assume_isolated == "pcc_0d" || inp.assume_isolated == "pcc_2d")
    {
        parameters.pcc_boundary = ModuleSccs::parse_boundary(inp.assume_isolated);
    }
    if (!parameters.use_sccs && parameters.pcc_boundary == ModuleSccs::Boundary::Periodic)
    {
        return parameters;
    }
    if (use_uspp)
    {
        ModuleBase::WARNING_QUIT("surchem", "SCCS/PCC currently supports only norm-conserving pseudopotentials");
    }
    parameters.pool_process_count = pool_process_count;
    parameters.debug = inp.sccs_debug;
    parameters.expected_electron_count = electron_count;
    for (int atom_type = 0; atom_type < ucell.ntype; ++atom_type)
    {
        parameters.expected_ionic_charge
            += ucell.atoms[atom_type].ncpp.zv * ucell.atoms[atom_type].na;
    }
    if (!parameters.use_sccs)
    {
        if (parameters.pcc_boundary == ModuleSccs::Boundary::Pcc2d
            && std::abs(parameters.expected_ionic_charge
                        - parameters.expected_electron_count)
                   > parameters.normalization_tolerance)
        {
            ModuleBase::WARNING(
                "surchem",
                "charged pcc_2d slab: absolute energies at different y cell lengths are not directly comparable");
        }
        return parameters;
    }
    const ModuleSccs::Preset preset = ModuleSccs::parse_preset(inp.sccs_preset);
    if (preset == ModuleSccs::Preset::Custom)
    {
        parameters.sccs_config.cavity.density_min = inp.sccs_rho_min;
        parameters.sccs_config.cavity.density_max = inp.sccs_rho_max;
        parameters.sccs_config.cavity.epsilon_bulk = inp.sccs_epsilon;
        parameters.sccs_config.surface_tension
            = ModuleSccs::dyn_per_cm_to_hartree_per_bohr2(inp.sccs_gamma);
        parameters.sccs_config.pressure
            = ModuleSccs::gpa_to_hartree_per_bohr3(inp.sccs_pressure);
    }
    else if (preset == ModuleSccs::Preset::Vacuum)
    {
        parameters.sccs_config = ModuleSccs::vacuum_preset();
    }
    else
    {
        parameters.sccs_config = ModuleSccs::water_preset(preset);
    }
    parameters.sccs_config.boundary = parameters.pcc_boundary;
    parameters.sccs_config.max_iterations = inp.sccs_maxiter;
    parameters.sccs_config.mixing_method = inp.sccs_mixing_type;
    parameters.sccs_config.mixing_history = inp.sccs_mixing_ndim;
    parameters.sccs_config.mixing = inp.sccs_mixing;
    parameters.sccs_config.adaptive_mixing = inp.sccs_mixing_adaptive;
    parameters.sccs_config.mixing_min = inp.sccs_mixing_min;
    parameters.sccs_config.mixing_max = inp.sccs_mixing_max;
    parameters.sccs_config.tolerance_rms = inp.sccs_tol_rms;
    parameters.sccs_config.tolerance_max = inp.sccs_tol_max;
    parameters.sccs_config.surface_regularization = inp.sccs_surface_eta;
    parameters.start_drho = inp.sccs_start_drho;
    parameters.start_nmax = inp.sccs_start_nmax;
    ModuleSccs::validate_config(parameters.sccs_config);

    const double net_charge
        = parameters.expected_ionic_charge - parameters.expected_electron_count;
    if (parameters.sccs_config.boundary == ModuleSccs::Boundary::Pcc2d
        && std::abs(net_charge) > parameters.normalization_tolerance)
    {
        ModuleBase::WARNING(
            "surchem",
            "charged SCCS pcc_2d slab: the open-boundary field energy grows "
            "linearly with the cell length, so absolute total energies at "
            "different y cell lengths are not directly comparable");
    }
    return parameters;
}

void validate_kpoints(const SurchemParameters& parameters,
                           const K_Vectors& kv)
{
    const ModuleSccs::Boundary boundary
        = parameters.use_sccs ? parameters.sccs_config.boundary : parameters.pcc_boundary;
    if (boundary != ModuleSccs::Boundary::Pcc2d)
    {
        return;
    }

    double maximum_open_direction_k = 0.0;
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        maximum_open_direction_k
            = std::max(maximum_open_direction_k, std::abs(kv.kvec_d[ik].y));
    }
    Parallel_Reduce::reduce_max(maximum_open_direction_k);
    if (maximum_open_direction_k > 1.0e-12)
    {
        ModuleBase::WARNING_QUIT(
            "surchem",
            "pcc_2d requires Gamma-only sampling along the second lattice direction");
    }
}

} // namespace ModuleSurchem
