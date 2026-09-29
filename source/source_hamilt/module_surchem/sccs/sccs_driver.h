#ifndef SCCS_DRIVER_H
#define SCCS_DRIVER_H

#include "sccs_charge.h"
#include "sccs_functional.h"
#include "sccs_nonel.h"
#include "sccs_parameters.h"
#include "../pcc/pcc_2d.h"
#include "sccs_response.h"

#include <cstdint>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{

struct SccsState
{
    // Unshifted sqrt-CG solution used for the next warm start.
    std::vector<double> potential;
    int local_grid_size = 0;
    int global_grid_size = 0;
    int nx = 0;
    int ny = 0;
    int nz = 0;
    int local_plane_count = 0;
    int local_plane_start = 0;
    std::uint64_t grid_position_signature = 0;
    ModulePcc::Boundary boundary = ModulePcc::Boundary::Periodic;
    double tpiba = 0.0;
    double volume_element = 0.0;
    ModulePcc::PccGeometry pcc_geometry;
    ModulePcc::Pcc2dGeometry pcc_2d_geometry;
    CavityParameters cavity;
    bool valid = false;
};

struct SccsResult
{
    bool reused_fixed_sources = false;
    CoulombTransformCounts forward_transforms;
    ChargeDensity charge;
    SccsResponse response;
    ElectrostaticFunctionalResult electrostatic;
    NonElectrostaticResult non_electrostatic;
    std::vector<double> electron_potential_hartree;
    // PCC0D moments about the PCC origin; zero for other boundaries.
    ModulePcc::MultipoleMoments solute_moments;
    ModulePcc::MultipoleMoments polarization_moments;
    ModulePcc::MultipoleMoments screened_moments;
    // PCC2D moments along the open direction; zero for other boundaries.
    ModulePcc::Pcc2dMoments solute_moments_2d;
    ModulePcc::Pcc2dMoments polarization_moments_2d;
    ModulePcc::Pcc2dMoments screened_moments_2d;
    double smooth_vacuum_pcc_energy = 0.0;
    ModulePcc::MultipoleMoments point_solute_moments;
    ModulePcc::Pcc2dMoments point_solute_moments_2d;
    double vacuum_pcc_energy = 0.0;
};

// cavity_core_density: ENVIRON 'full' core-electron Gaussians added to the
// electron density that defines the cavity (config.core_electrons), else empty.
// The reciprocal unit tpiba and the volume element omega/nxyz come from basis.
SccsResult evaluate_pw_sccs(
    const std::vector<double>& electron_density,
    const std::vector<double>& ionic_density,
    const std::vector<double>& cavity_core_density,
    double expected_electron_count,
    double expected_ionic_charge,
    double normalization_tolerance,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    const SccsConfig& config,
    const ModulePcc::PccGeometry& pcc_geometry,
    const ModulePcc::Pcc2dGeometry& pcc_2d_geometry,
    const ModulePW::PW_Basis& basis,
    const ModuleSurchem::ChargeReduction& reduction,
    SccsState& state);

} // namespace ModuleSccs

#endif
