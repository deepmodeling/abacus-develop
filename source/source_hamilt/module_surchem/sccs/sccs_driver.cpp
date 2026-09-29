#include "sccs_driver.h"

#include "../pcc/pcc_moments.h"
#include "sccs_pcc_2d_coulomb.h"
#include "sccs_pcc_coulomb.h"
#include "sccs_pw_coulomb.h"
#include "sccs_pw_nonel.h"

#include "source_base/timer.h"
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
           && left.epsilon_bulk == right.epsilon_bulk
           && left.lowpass_p1 == right.lowpass_p1
           && left.lowpass_p2 == right.lowpass_p2;
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
                          const ModulePcc::Boundary boundary,
                          const ModulePcc::PccGeometry& pcc_geometry,
                          const ModulePcc::Pcc2dGeometry& pcc_2d_geometry,
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
           && state.potential.size() == static_cast<std::size_t>(basis.nrxx);
}

ModulePcc::MultipoleMoments add_moments(const ModulePcc::MultipoleMoments& left, const ModulePcc::MultipoleMoments& right)
{
    ModulePcc::MultipoleMoments result;
    result.charge = left.charge + right.charge;
    result.dipole.x = left.dipole.x + right.dipole.x;
    result.dipole.y = left.dipole.y + right.dipole.y;
    result.dipole.z = left.dipole.z + right.dipole.z;
    result.quadrupole_trace = left.quadrupole_trace + right.quadrupole_trace;
    return result;
}

ModulePcc::Pcc2dMoments add_moments_2d(const ModulePcc::Pcc2dMoments& left, const ModulePcc::Pcc2dMoments& right)
{
    ModulePcc::Pcc2dMoments result;
    result.charge = left.charge + right.charge;
    result.dipole_y = left.dipole_y + right.dipole_y;
    result.quadrupole_yy = left.quadrupole_yy + right.quadrupole_yy;
    return result;
}

} // namespace

void SccsState::reset()
{
    potential.clear();
    local_grid_size = 0;
    global_grid_size = 0;
    nx = 0;
    ny = 0;
    nz = 0;
    local_plane_count = 0;
    local_plane_start = 0;
    grid_position_signature = 0;
    boundary = ModulePcc::Boundary::Periodic;
    tpiba = 0.0;
    volume_element = 0.0;
    origin = ModuleBase::Vector3<double>();
    pcc_geometry = ModulePcc::PccGeometry();
    pcc_2d_geometry = ModulePcc::Pcc2dGeometry();
    cavity = CavityParameters();
    valid = false;
}

// Assemble q = rho_ion - n, solve the ENVIRON-style sqrt-CG response, then
// combine electrostatic and cavity terms. Energies/potentials here are in Ha;
// the surchem adapter adds the point-ion vacuum PCC and converts to Ry.
SccsResult evaluate_pw_sccs(
    const std::vector<double>& electron_density,
    const std::vector<double>& ionic_density,
    const std::vector<double>& cavity_core_density,
    const double expected_electron_count,
    const double expected_ionic_charge,
    const double normalization_tolerance,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    const ModuleBase::Vector3<double>& origin,
    const SccsConfig& config,
    const ModulePcc::PccGeometry& pcc_geometry,
    const ModulePcc::Pcc2dGeometry& pcc_2d_geometry,
    const ModulePW::PW_Basis& basis,
    const double tpiba,
    const double volume_element,
    const ModuleSurchem::ChargeReduction& charge_reduction,
    const PolarizationReduction& polarization_reduction,
    SccsState& state)
{
    ModuleBase::timer::start("ModuleSccs", "evaluate_pw_sccs");
    validate_config(config);
    if (positions.size() != static_cast<std::size_t>(basis.nrxx))
    {
        throw std::invalid_argument("SCCS grid positions must match the local PW grid");
    }
    const std::size_t expected_core_size
        = config.core_electrons ? static_cast<std::size_t>(basis.nrxx) : 0;
    if (cavity_core_density.size() != expected_core_size)
    {
        throw std::invalid_argument("SCCS core-electron density must be given exactly when core_electrons is set");
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
    if (config.boundary == ModulePcc::Boundary::Pcc0d)
    {
        ModulePcc::validate_pcc_geometry(pcc_geometry);
        coulomb.reset(new PccCoulombOperator(basis,
                                             tpiba,
                                             positions,
                                             volume_element,
                                             pcc_geometry,
                                             charge_reduction));
    }
    else if (config.boundary == ModulePcc::Boundary::Pcc2d)
    {
        ModulePcc::validate_pcc_2d_parameters(pcc_2d_geometry.parameters);
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
    solver_parameters.tolerance_rms = config.tolerance_rms;
    solver_parameters.tolerance_max = config.tolerance_max;
    solver_parameters.check_fixed_point = config.check_fixed_point;
    const std::uint64_t position_signature = grid_position_signature(positions);
    // A changed grid, cavity or PCC origin invalidates the warm-start potential.
    // Borrow the cache until the solve succeeds; state is updated only below.
    const bool reuse_state = same_state_signature(state,
                                                   config.boundary,
                                                   pcc_geometry,
                                                   pcc_2d_geometry,
                                                   config.cavity,
                                                   basis,
                                                   tpiba,
                                                   volume_element,
                                                   origin,
                                                   position_signature);
    const std::vector<double> empty_initial;
    const std::vector<double>& initial_potential
        = reuse_state ? state.potential : empty_initial;
    // Every boundary uses the ENVIRON sqrt-CG. With PCC the preconditioner's
    // Poisson solve includes the analytic open-boundary term, so a charged
    // solute keeps its screening charge and the physical potential gauge.
    // The cavity follows the electrons, plus the core electrons in ENVIRON
    // 'full' mode; the solute charge is unchanged.
    std::vector<double> cavity_density = result.charge.electron;
    for (std::size_t index = 0; index < cavity_core_density.size(); ++index)
    {
        cavity_density[index] += cavity_core_density[index];
    }
    result.response = solve_chain_sccs_response(cavity_density,
                                                result.charge.solute,
                                                config.cavity,
                                                solver_parameters,
                                                initial_potential,
                                                basis,
                                                tpiba,
                                                *coulomb,
                                                polarization_reduction);
    if (result.response.polarization.status != PolarizationStatus::Converged)
    {
        throw std::runtime_error("SCCS polarization iteration did not converge");
    }
    result.forward_transforms = coulomb->transform_counts();

    std::vector<double> vacuum_potential;
    coulomb->apply_potential(result.charge.solute, vacuum_potential);
    result.electrostatic = evaluate_electrostatic_functional(result.charge.solute,
                                                              result.response.polarization.field.potential,
                                                              vacuum_potential,
                                                              result.response.cavity_potential,
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for (std::size_t index = 0; index < electron_density.size(); ++index)
    {
        result.electron_potential_hartree[index]
            = result.electrostatic.electron_potential[index]
              + result.non_electrostatic.density_potential[index];
    }

    if (config.boundary == ModulePcc::Boundary::Pcc0d)
    {
        result.solute_moments
            = ModulePcc::reduced_pcc_density_moments(result.charge.solute,
                                          positions,
                                          volume_element,
                                          pcc_geometry,
                                          charge_reduction);
        result.polarization_moments
            = ModulePcc::reduced_pcc_density_moments(
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
    if (config.boundary == ModulePcc::Boundary::Pcc0d)
    {
        result.smooth_vacuum_pcc_energy
            = ModulePcc::pcc_self_energy(result.solute_moments, pcc_geometry.parameters);
    }
    else if (config.boundary == ModulePcc::Boundary::Pcc2d)
    {
        result.solute_moments_2d
            = ModulePcc::reduced_pcc_2d_density_moments(result.charge.solute,
                                             positions,
                                             volume_element,
                                             pcc_2d_geometry,
                                             charge_reduction);
        result.polarization_moments_2d
            = ModulePcc::reduced_pcc_2d_density_moments(result.response.polarization.polarization_charge,
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
        // Gauss's law for the solution's far field. The dielectric_of_potential
        // density integral (polarization_moments_2d) is only a diagnostic here:
        // its finite-grid error can exceed this tolerance for sharp cavities.
        const double far_field_charge = result.response.far_field_polarization_charge;
        if (std::abs(far_field_charge - expected_polarization_charge)
            > polarization_charge_tolerance)
        {
            std::ostringstream message;
            message << "SCCS pcc_2d far-field polarization charge "
                    << far_field_charge
                    << " differs from expected " << expected_polarization_charge
                    << " by more than tolerance " << polarization_charge_tolerance;
            throw std::runtime_error(message.str());
        }
        result.smooth_vacuum_pcc_energy
            = ModulePcc::pcc_2d_self_energy(result.solute_moments_2d,
                                 pcc_2d_geometry.parameters);
    }

    state.potential = result.response.restart_potential;
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
    ModuleBase::timer::end("ModuleSccs", "evaluate_pw_sccs");
    return result;
}

} // namespace ModuleSccs
