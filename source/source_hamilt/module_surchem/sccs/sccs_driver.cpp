#include "sccs_driver.h"

#include "../pcc/sccs_pcc_2d_coulomb.h"
#include "../pcc/sccs_pcc_coulomb.h"
#include "sccs_pw_coulomb.h"
#include "sccs_pw_nonel.h"

#include "source_basis/module_pw/pw_basis.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace ModuleSccs
{
namespace
{

const double relative_polarization_charge_tolerance = 1.0e-4;

bool same_cavity(const CavityParameters& left, const CavityParameters& right)
{
    return left.density_min == right.density_min
           && left.density_max == right.density_max
           && left.epsilon_bulk == right.epsilon_bulk;
}

bool same_vector(const ModuleBase::Vector3<double>& left,
                 const ModuleBase::Vector3<double>& right)
{
    return left.x == right.x && left.y == right.y && left.z == right.z;
}

void hash_double(std::uint64_t& hash, const double value)
{
    unsigned char bytes[sizeof(double)];
    std::memcpy(bytes, &value, sizeof(double));
    for (std::size_t index = 0; index < sizeof(double); ++index)
    {
        hash ^= static_cast<std::uint64_t>(bytes[index]);
        hash *= static_cast<std::uint64_t>(1099511628211ULL);
    }
}

std::uint64_t grid_position_signature(
    const std::vector<ModuleBase::Vector3<double>>& positions)
{
    std::uint64_t hash = static_cast<std::uint64_t>(1469598103934665603ULL);
    for (std::size_t index = 0; index < positions.size(); ++index)
    {
        hash_double(hash, positions[index].x);
        hash_double(hash, positions[index].y);
        hash_double(hash, positions[index].z);
    }
    return hash;
}

bool same_state_signature(const SccsState& state,
                          const Boundary boundary,
                          const PccGeometry& pcc_geometry,
                          const Pcc2dGeometry& pcc_2d_geometry,
                          const CavityParameters& cavity,
                          const ModulePW::PW_Basis& basis,
                          const double tpiba,
                          const double volume_element,
                          const ModuleBase::Vector3<double>& origin,
                          const std::uint64_t position_signature)
{
    return state.valid && state.boundary == boundary
           && state.local_grid_size == basis.nrxx
           && state.global_grid_size == basis.nxyz && state.nx == basis.nx
           && state.ny == basis.ny && state.nz == basis.nz
           && state.local_plane_count == basis.nplane
           && state.local_plane_start == basis.startz_current
           && state.grid_position_signature == position_signature
           && state.tpiba == tpiba && state.volume_element == volume_element
           && same_vector(state.origin, origin)
           && state.pcc_geometry.parameters.cube_length
                  == pcc_geometry.parameters.cube_length
           && state.pcc_geometry.parameters.madelung
                  == pcc_geometry.parameters.madelung
           && same_vector(state.pcc_geometry.origin, pcc_geometry.origin)
           && same_vector(state.pcc_geometry.axis_a, pcc_geometry.axis_a)
           && same_vector(state.pcc_geometry.axis_b, pcc_geometry.axis_b)
           && same_vector(state.pcc_geometry.axis_c, pcc_geometry.axis_c)
           && state.pcc_2d_geometry.parameters.periodic_area
                  == pcc_2d_geometry.parameters.periodic_area
           && state.pcc_2d_geometry.parameters.cell_length_y
                  == pcc_2d_geometry.parameters.cell_length_y
           && state.pcc_2d_geometry.origin_y == pcc_2d_geometry.origin_y
           && same_cavity(state.cavity, cavity)
           && state.polarization_charge.size() == static_cast<std::size_t>(basis.nrxx);
}

std::vector<double> initial_polarization(const SccsState& state,
                                         const Boundary boundary,
                                         const PccGeometry& pcc_geometry,
                                         const Pcc2dGeometry& pcc_2d_geometry,
                                         const CavityParameters& cavity,
                                         const ModulePW::PW_Basis& basis,
                                         const double tpiba,
                                         const double volume_element,
                                         const ModuleBase::Vector3<double>& origin,
                                         const std::uint64_t position_signature)
{
    if (!same_state_signature(state,
                              boundary,
                              pcc_geometry,
                              pcc_2d_geometry,
                              cavity,
                              basis,
                              tpiba,
                              volume_element,
                              origin,
                              position_signature))
    {
        return std::vector<double>();
    }
    return state.polarization_charge;
}

MultipoleMoments add_moments(const MultipoleMoments& left, const MultipoleMoments& right)
{
    MultipoleMoments result;
    result.charge = left.charge + right.charge;
    result.dipole.x = left.dipole.x + right.dipole.x;
    result.dipole.y = left.dipole.y + right.dipole.y;
    result.dipole.z = left.dipole.z + right.dipole.z;
    result.quadrupole_trace = left.quadrupole_trace + right.quadrupole_trace;
    return result;
}

Pcc2dMoments add_moments_2d(const Pcc2dMoments& left, const Pcc2dMoments& right)
{
    Pcc2dMoments result;
    result.charge = left.charge + right.charge;
    result.dipole_y = left.dipole_y + right.dipole_y;
    result.quadrupole_yy = left.quadrupole_yy + right.quadrupole_yy;
    return result;
}

} // namespace

void SccsState::reset()
{
    polarization_charge.clear();
    local_grid_size = 0;
    global_grid_size = 0;
    nx = 0;
    ny = 0;
    nz = 0;
    local_plane_count = 0;
    local_plane_start = 0;
    grid_position_signature = 0;
    boundary = Boundary::Periodic;
    tpiba = 0.0;
    volume_element = 0.0;
    origin = ModuleBase::Vector3<double>();
    pcc_geometry = PccGeometry();
    pcc_2d_geometry = Pcc2dGeometry();
    cavity = CavityParameters();
    valid = false;
}

SccsResult evaluate_pw_sccs(
    const std::vector<double>& electron_density,
    const std::vector<double>& ionic_density,
    const double expected_electron_count,
    const double expected_ionic_charge,
    const double normalization_tolerance,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    const ModuleBase::Vector3<double>& origin,
    const SccsConfig& config,
    const PccGeometry& pcc_geometry,
    const Pcc2dGeometry& pcc_2d_geometry,
    const ModulePW::PW_Basis& basis,
    const double tpiba,
    const double volume_element,
    const ChargeReduction& charge_reduction,
    const PolarizationReduction& polarization_reduction,
    SccsState& state)
{
    validate_config(config);
    if (positions.size() != static_cast<std::size_t>(basis.nrxx))
    {
        throw std::invalid_argument("SCCS grid positions must match the local PW grid");
    }

    SccsResult result;
    result.charge = assemble_charge_density(electron_density,
                                            ionic_density,
                                            volume_element,
                                            expected_electron_count,
                                            expected_ionic_charge,
                                            normalization_tolerance,
                                            charge_reduction);
    std::unique_ptr<CoulombOperator> coulomb;
    if (config.boundary == Boundary::Pcc0d)
    {
        validate_pcc_geometry(pcc_geometry);
        coulomb.reset(new PccCoulombOperator(basis,
                                             tpiba,
                                             positions,
                                             volume_element,
                                             pcc_geometry,
                                             charge_reduction));
    }
    else if (config.boundary == Boundary::Pcc2d)
    {
        validate_pcc_2d_parameters(pcc_2d_geometry.parameters);
        coulomb.reset(new Pcc2dCoulombOperator(basis,
                                               tpiba,
                                               positions,
                                               volume_element,
                                               pcc_2d_geometry,
                                               charge_reduction));
    }
    else
    {
        coulomb.reset(new PeriodicCoulombOperator(basis, tpiba));
    }

    PolarizationSolverParameters solver_parameters;
    solver_parameters.max_iterations = config.max_iterations;
    solver_parameters.mixing_method = config.mixing_method;
    solver_parameters.mixing_history = config.mixing_history;
    solver_parameters.mixing = config.mixing;
    solver_parameters.adaptive_mixing = config.adaptive_mixing;
    solver_parameters.mixing_min = config.mixing_min;
    solver_parameters.mixing_max = config.mixing_max;
    solver_parameters.tolerance_rms = config.tolerance_rms;
    solver_parameters.tolerance_max = config.tolerance_max;
    const std::uint64_t position_signature = grid_position_signature(positions);
    const std::vector<double> initial = initial_polarization(state,
                                                             config.boundary,
                                                             pcc_geometry,
                                                             pcc_2d_geometry,
                                                             config.cavity,
                                                             basis,
                                                             tpiba,
                                                             volume_element,
                                                             origin,
                                                             position_signature);
    result.response = solve_sccs_response(result.charge.electron,
                                          result.charge.solute,
                                          config.cavity,
                                          solver_parameters,
                                          initial,
                                          basis,
                                          tpiba,
                                          *coulomb,
                                          polarization_reduction);
    if (result.response.polarization.status != PolarizationStatus::Converged)
    {
        throw std::runtime_error("SCCS polarization iteration did not converge");
    }

    coulomb->apply(result.charge.solute, result.vacuum_field);
    result.electrostatic = evaluate_electrostatic_functional(result.charge.solute,
                                                              result.response.polarization.field,
                                                              result.vacuum_field,
                                                              result.response.depsilon_drho,
                                                              volume_element,
                                                              charge_reduction);

    NonElectrostaticParameters non_electrostatic_parameters;
    non_electrostatic_parameters.surface_tension = config.surface_tension;
    non_electrostatic_parameters.pressure = config.pressure;
    non_electrostatic_parameters.surface_regularization = config.surface_regularization;
    result.non_electrostatic = evaluate_pw_non_electrostatic(basis,
                                                             tpiba,
                                                             volume_element,
                                                             non_electrostatic_parameters,
                                                             result.response.solute,
                                                             result.response.dsolute_drho,
                                                             charge_reduction);

    result.electron_potential_hartree.resize(electron_density.size());
    for (std::size_t index = 0; index < electron_density.size(); ++index)
    {
        result.electron_potential_hartree[index]
            = result.electrostatic.electron_potential[index]
              + result.non_electrostatic.density_potential[index];
    }

    if (config.boundary == Boundary::Pcc0d)
    {
        result.solute_moments
            = reduced_pcc_density_moments(result.charge.solute,
                                          positions,
                                          volume_element,
                                          pcc_geometry,
                                          charge_reduction);
        result.polarization_moments
            = reduced_pcc_density_moments(
                result.response.polarization.polarization_charge,
                positions,
                volume_element,
                pcc_geometry,
                charge_reduction);
    }
    else
    {
        result.solute_moments = reduced_density_moments(result.charge.solute,
                                                        positions,
                                                        volume_element,
                                                        origin,
                                                        charge_reduction);
        result.polarization_moments
            = reduced_density_moments(result.response.polarization.polarization_charge,
                                      positions,
                                      volume_element,
                                      origin,
                                      charge_reduction);
    }
    result.screened_moments = add_moments(result.solute_moments, result.polarization_moments);
    if (config.boundary == Boundary::Pcc0d)
    {
        result.smooth_vacuum_pcc_energy
            = pcc_self_energy(result.solute_moments, pcc_geometry.parameters);
    }
    else if (config.boundary == Boundary::Pcc2d)
    {
        result.solute_moments_2d
            = reduced_pcc_2d_density_moments(result.charge.solute,
                                             positions,
                                             volume_element,
                                             pcc_2d_geometry,
                                             charge_reduction);
        result.polarization_moments_2d
            = reduced_pcc_2d_density_moments(result.response.polarization.polarization_charge,
                                             positions,
                                             volume_element,
                                             pcc_2d_geometry,
                                             charge_reduction);
        result.screened_moments_2d = add_moments_2d(result.solute_moments_2d,
                                                    result.polarization_moments_2d);
        const double polarization_charge_tolerance
            = std::max(normalization_tolerance,
                       std::max(config.tolerance_max * volume_element
                                    * static_cast<double>(basis.nxyz),
                                relative_polarization_charge_tolerance
                                    * std::max(1.0,
                                               std::max(expected_electron_count,
                                                        expected_ionic_charge))));
        const double expected_polarization_charge
            = -(1.0 - 1.0 / config.cavity.epsilon_bulk)
              * result.solute_moments_2d.charge;
        if (std::abs(result.polarization_moments_2d.charge
                     - expected_polarization_charge)
            > polarization_charge_tolerance)
        {
            std::ostringstream message;
            message << "SCCS pcc_2d polarization charge "
                    << result.polarization_moments_2d.charge
                    << " differs from expected " << expected_polarization_charge
                    << " by more than tolerance " << polarization_charge_tolerance;
            throw std::runtime_error(message.str());
        }
        result.smooth_vacuum_pcc_energy
            = pcc_2d_self_energy(result.solute_moments_2d,
                                 pcc_2d_geometry.parameters);
    }

    state.polarization_charge = result.response.polarization.polarization_charge;
    state.local_grid_size = basis.nrxx;
    state.global_grid_size = basis.nxyz;
    state.nx = basis.nx;
    state.ny = basis.ny;
    state.nz = basis.nz;
    state.local_plane_count = basis.nplane;
    state.local_plane_start = basis.startz_current;
    state.grid_position_signature = position_signature;
    state.boundary = config.boundary;
    state.tpiba = tpiba;
    state.volume_element = volume_element;
    state.origin = origin;
    state.pcc_geometry = pcc_geometry;
    state.pcc_2d_geometry = pcc_2d_geometry;
    state.cavity = config.cavity;
    state.valid = true;
    return result;
}

} // namespace ModuleSccs
