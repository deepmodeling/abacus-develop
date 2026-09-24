#ifndef SCCS_DRIVER_H
#define SCCS_DRIVER_H

#include "sccs_charge.h"
#include "sccs_functional.h"
#include "sccs_nonel.h"
#include "sccs_parameters.h"
#include "sccs_pcc_2d.h"
#include "sccs_periodic.h"

#include <cstdint>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{

struct SccsState
{
    std::vector<double> polarization_charge;
    int local_grid_size = 0;
    int global_grid_size = 0;
    int nx = 0;
    int ny = 0;
    int nz = 0;
    int local_plane_count = 0;
    int local_plane_start = 0;
    std::uint64_t grid_position_signature = 0;
    Boundary boundary = Boundary::Periodic;
    double tpiba = 0.0;
    double volume_element = 0.0;
    ModuleBase::Vector3<double> origin;
    PccGeometry pcc_geometry;
    Pcc2dGeometry pcc_2d_geometry;
    CavityParameters cavity;
    bool valid = false;

    void reset();
};

struct SccsResult
{
    ChargeDensity charge;
    PeriodicSccsResult response;
    ElectrostaticField vacuum_field;
    ElectrostaticFunctionalResult electrostatic;
    NonElectrostaticResult non_electrostatic;
    std::vector<double> electron_potential_hartree;
    MultipoleMoments solute_moments;
    MultipoleMoments polarization_moments;
    MultipoleMoments screened_moments;
    Pcc2dMoments solute_moments_2d;
    Pcc2dMoments polarization_moments_2d;
    Pcc2dMoments screened_moments_2d;
    double smooth_vacuum_pcc_energy = 0.0;
    MultipoleMoments point_solute_moments;
    Pcc2dMoments point_solute_moments_2d;
    double vacuum_pcc_energy = 0.0;
    double ionic_shape_pcc_energy = 0.0;
};

SccsResult evaluate_pw_sccs(
    const std::vector<double>& electron_density,
    const std::vector<double>& ionic_density,
    double expected_electron_count,
    double expected_ionic_charge,
    double normalization_tolerance,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    const ModuleBase::Vector3<double>& origin,
    const SccsConfig& config,
    const PccGeometry& pcc_geometry,
    const Pcc2dGeometry& pcc_2d_geometry,
    const ModulePW::PW_Basis& basis,
    double tpiba,
    double volume_element,
    const ChargeReduction& charge_reduction,
    const PolarizationReduction& polarization_reduction,
    SccsState& state);

} // namespace ModuleSccs

#endif
