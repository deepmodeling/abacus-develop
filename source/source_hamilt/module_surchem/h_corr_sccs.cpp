#include "surchem.h"

#include "sccs_pcc_2d_coulomb.h"
#include "sccs_pcc_coulomb.h"
#include "sccs_pw_charge.h"
#include "sccs_pw_reduction.h"

#include "source_base/timer.h"
#include "source_base/timer_wrapper.h"
#include "source_base/tool_title.h"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{

ModuleSccs::MultipoleMoments add_moments(const ModuleSccs::MultipoleMoments& left,
                                         const ModuleSccs::MultipoleMoments& right)
{
    ModuleSccs::MultipoleMoments result;
    result.charge = left.charge + right.charge;
    result.dipole.x = left.dipole.x + right.dipole.x;
    result.dipole.y = left.dipole.y + right.dipole.y;
    result.dipole.z = left.dipole.z + right.dipole.z;
    result.quadrupole_trace = left.quadrupole_trace + right.quadrupole_trace;
    return result;
}

ModuleSccs::Pcc2dMoments add_moments_2d(const ModuleSccs::Pcc2dMoments& left,
                                       const ModuleSccs::Pcc2dMoments& right)
{
    ModuleSccs::Pcc2dMoments result;
    result.charge = left.charge + right.charge;
    result.dipole_y = left.dipole_y + right.dipole_y;
    result.quadrupole_yy = left.quadrupole_yy + right.quadrupole_yy;
    return result;
}

double ionic_system_center_y(const UnitCell& cell, const double cell_length_y)
{
    std::vector<double> positions_y;
    std::vector<double> masses;
    positions_y.reserve(cell.nat);
    masses.reserve(cell.nat);
    for (int atom_type = 0; atom_type < cell.ntype; ++atom_type)
    {
        for (int atom = 0; atom < cell.atoms[atom_type].na; ++atom)
        {
            positions_y.push_back(cell.atoms[atom_type].tau[atom].y * cell.lat0);
            masses.push_back(cell.atoms[atom_type].mass);
        }
    }
    return ModuleSccs::pcc_2d_system_center_y(positions_y,
                                               masses,
                                               cell_length_y);
}

ModuleBase::Vector3<double> ionic_system_center(
    const UnitCell& cell,
    const ModuleSccs::PccGeometry& geometry)
{
    std::vector<ModuleBase::Vector3<double>> positions;
    std::vector<double> masses;
    positions.reserve(cell.nat);
    masses.reserve(cell.nat);
    for (int atom_type = 0; atom_type < cell.ntype; ++atom_type)
    {
        for (int atom = 0; atom < cell.atoms[atom_type].na; ++atom)
        {
            positions.push_back(cell.atoms[atom_type].tau[atom] * cell.lat0);
            masses.push_back(cell.atoms[atom_type].mass);
        }
    }
    return ModuleSccs::pcc_system_center(positions, masses, geometry);
}

std::vector<ModuleSccs::PointCharge> ionic_point_charges(const UnitCell& cell)
{
    std::vector<ModuleSccs::PointCharge> charges;
    charges.reserve(cell.nat);
    for (int atom_type = 0; atom_type < cell.ntype; ++atom_type)
    {
        for (int atom = 0; atom < cell.atoms[atom_type].na; ++atom)
        {
            ModuleSccs::PointCharge point;
            point.charge = cell.atoms[atom_type].ncpp.zv;
            point.position = cell.atoms[atom_type].tau[atom] * cell.lat0;
            charges.push_back(point);
        }
    }
    return charges;
}

} // namespace

void surchem::v_correction_sccs(const UnitCell& cell,
                                const ModulePW::PW_Basis& rho_basis,
                                const int nspin,
                                const double* const* rho,
                                const double* vlocal,
                                ModuleBase::matrix& v)
{
    ModuleBase::TITLE("surchem", "v_correction_sccs");
    ModuleBase::timer::start("surchem", "v_correction_sccs");
    const ModuleBase::TimePoint start_time = ModuleBase::get_time();
    if (!this->parameters_set_ || !this->parameters_.use_sccs)
    {
        throw std::logic_error("SCCS correction requires an initialized SCCS configuration");
    }
    if (nspin != 1 && nspin != 2)
    {
        throw std::invalid_argument("SCCS currently supports nspin=1 or nspin=2");
    }

    std::vector<std::vector<double>> spin_density(nspin);
    for (int spin = 0; spin < nspin; ++spin)
    {
        spin_density[spin].assign(rho[spin], rho[spin] + rho_basis.nrxx);
    }
    const std::vector<double> electron_density
        = ModuleSccs::sum_electron_density(spin_density, nspin);
    const std::vector<double> local_potential(vlocal, vlocal + rho_basis.nrxx);
    const std::vector<double> ionic_density
        = ModuleSccs::ionic_charge_from_local_potential(local_potential,
                                                        this->parameters_.expected_ionic_charge,
                                                        cell.omega,
                                                        cell.tpiba,
                                                        rho_basis);
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSccs::pw_grid_positions(rho_basis, cell.latvec, cell.lat0);
    ModuleBase::Vector3<double> origin
        = ModuleSccs::cell_center(cell.latvec, cell.lat0);
    ModuleSccs::PccGeometry pcc_geometry;
    ModuleSccs::Pcc2dGeometry pcc_2d_geometry;
    if (this->parameters_.sccs_config.boundary == ModuleSccs::Boundary::Pcc0d)
    {
        pcc_geometry = ModuleSccs::pcc_geometry(cell.latvec,
                                                cell.lat0,
                                                1.0e-10);
        pcc_geometry.origin = ionic_system_center(cell, pcc_geometry);
        origin = pcc_geometry.origin;
    }
    else if (this->parameters_.sccs_config.boundary == ModuleSccs::Boundary::Pcc2d)
    {
        pcc_2d_geometry = ModuleSccs::pcc_2d_geometry(cell.latvec,
                                                      cell.lat0,
                                                      1.0e-10);
        pcc_2d_geometry.origin_y
            = ionic_system_center_y(cell,
                                    pcc_2d_geometry.parameters.cell_length_y);
    }

    const double volume_element = cell.omega / static_cast<double>(rho_basis.nxyz);
    const ModuleSccs::PoolChargeReduction charge_reduction;
    const ModuleSccs::PoolPolarizationReduction polarization_reduction(
        this->parameters_.pool_process_count);
    this->sccs_result_
        = ModuleSccs::evaluate_pw_sccs(electron_density,
                                       ionic_density,
                                       this->parameters_.expected_electron_count,
                                       this->parameters_.expected_ionic_charge,
                                       this->parameters_.normalization_tolerance,
                                       positions,
                                       origin,
                                       this->parameters_.sccs_config,
                                       pcc_geometry,
                                       pcc_2d_geometry,
                                       rho_basis,
                                       cell.tpiba,
                                       volume_element,
                                       charge_reduction,
                                       polarization_reduction,
                                       this->sccs_state_);

    if (this->parameters_.sccs_config.boundary == ModuleSccs::Boundary::Pcc0d)
    {
        std::vector<double> electronic_charge(electron_density.size());
        for (std::size_t index = 0; index < electron_density.size(); ++index)
        {
            electronic_charge[index] = -electron_density[index];
        }
        const ModuleSccs::MultipoleMoments electronic_moments
            = ModuleSccs::reduced_pcc_density_moments(electronic_charge,
                                                      positions,
                                                      volume_element,
                                                      pcc_geometry,
                                                      charge_reduction);
        const ModuleSccs::MultipoleMoments ionic_moments
            = ModuleSccs::point_charge_moments(ionic_point_charges(cell),
                                               pcc_geometry);
        this->sccs_result_.point_solute_moments = add_moments(ionic_moments,
                                                              electronic_moments);
        if (std::abs(this->sccs_result_.point_solute_moments.charge
                     - this->sccs_result_.charge.net_charge)
            > this->parameters_.normalization_tolerance)
        {
            throw std::runtime_error("SCCS point-ion and smooth-source net charges disagree");
        }
        this->sccs_result_.vacuum_pcc_energy
            = ModuleSccs::pcc_self_energy(this->sccs_result_.point_solute_moments,
                                          pcc_geometry.parameters);
        for (std::size_t index = 0; index < positions.size(); ++index)
        {
            this->sccs_result_.electron_potential_hartree[index]
                -= ModuleSccs::pcc_potential(this->sccs_result_.point_solute_moments,
                                              ModuleSccs::pcc_relative_position(
                                                  positions[index],
                                                  pcc_geometry),
                                              pcc_geometry.parameters);
        }
    }
    else if (this->parameters_.sccs_config.boundary == ModuleSccs::Boundary::Pcc2d)
    {
        std::vector<double> electronic_charge(electron_density.size());
        for (std::size_t index = 0; index < electron_density.size(); ++index)
        {
            electronic_charge[index] = -electron_density[index];
        }
        const ModuleSccs::Pcc2dMoments electronic_moments
            = ModuleSccs::reduced_pcc_2d_density_moments(electronic_charge,
                                                         positions,
                                                         volume_element,
                                                         pcc_2d_geometry,
                                                         charge_reduction);
        const ModuleSccs::Pcc2dMoments ionic_moments
            = ModuleSccs::pcc_2d_point_charge_moments(ionic_point_charges(cell),
                                                      pcc_2d_geometry);
        const ModuleSccs::Pcc2dMoments smooth_ionic_moments
            = ModuleSccs::reduced_pcc_2d_density_moments(
                this->sccs_result_.charge.ionic,
                positions,
                volume_element,
                pcc_2d_geometry,
                charge_reduction);
        this->sccs_result_.point_solute_moments_2d = add_moments_2d(ionic_moments,
                                                                   electronic_moments);
        if (std::abs(this->sccs_result_.point_solute_moments_2d.charge
                     - this->sccs_result_.charge.net_charge)
            > this->parameters_.normalization_tolerance)
        {
            throw std::runtime_error("SCCS point-ion and smooth-source net charges disagree");
        }
        this->sccs_result_.vacuum_pcc_energy
            = ModuleSccs::pcc_2d_self_energy(this->sccs_result_.point_solute_moments_2d,
                                             pcc_2d_geometry.parameters);
        this->sccs_result_.ionic_shape_pcc_energy
            = ModuleSccs::pcc_2d_ionic_shape_energy(
                this->sccs_result_.polarization_moments_2d.charge,
                smooth_ionic_moments,
                ionic_moments,
                pcc_2d_geometry.parameters);
        for (std::size_t index = 0; index < positions.size(); ++index)
        {
            this->sccs_result_.electron_potential_hartree[index]
                -= ModuleSccs::pcc_2d_potential(this->sccs_result_.point_solute_moments_2d,
                                                 ModuleSccs::pcc_2d_relative_y(
                                                     positions[index].y,
                                                     pcc_2d_geometry),
                                                 pcc_2d_geometry.parameters);
        }
    }

    if (v.nr != nspin || v.nc != rho_basis.nrxx)
    {
        v.create(nspin, rho_basis.nrxx);
    }
    ModuleBase::GlobalFunc::ZEROS(v.c, nspin * rho_basis.nrxx);
    for (int spin = 0; spin < nspin; ++spin)
    {
        for (int ir = 0; ir < rho_basis.nrxx; ++ir)
        {
            v(spin, ir) = 2.0 * this->sccs_result_.electron_potential_hartree[ir];
        }
    }

    surchem::Ael = 2.0 * (this->sccs_result_.electrostatic.reaction_energy
                          + this->sccs_result_.vacuum_pcc_energy
                          + this->sccs_result_.ionic_shape_pcc_energy);
    surchem::Acav = 2.0 * (this->sccs_result_.non_electrostatic.surface_energy
                           + this->sccs_result_.non_electrostatic.volume_energy);
    this->sccs_elapsed_seconds_
        = ModuleBase::get_duration(start_time, ModuleBase::get_time());
    ModuleBase::timer::end("surchem", "v_correction_sccs");
}

void surchem::v_correction_pcc(const UnitCell& cell,
                               const ModulePW::PW_Basis& rho_basis,
                               const int nspin,
                               const double* const* rho,
                               ModuleBase::matrix& v)
{
    if (!this->uses_pcc() || this->uses_sccs())
    {
        throw std::logic_error("standalone PCC requires PCC without SCCS");
    }
    if (nspin != 1 && nspin != 2)
    {
        throw std::invalid_argument("PCC supports nspin=1 or nspin=2");
    }
    std::vector<double> electronic_charge(rho_basis.nrxx, 0.0);
    for (int spin = 0; spin < nspin; ++spin)
    {
        for (int ir = 0; ir < rho_basis.nrxx; ++ir)
        {
            electronic_charge[ir] -= rho[spin][ir];
        }
    }
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSccs::pw_grid_positions(rho_basis, cell.latvec, cell.lat0);
    const double volume_element = cell.omega / static_cast<double>(rho_basis.nxyz);
    const ModuleSccs::PoolChargeReduction reduction;
    const std::vector<ModuleSccs::PointCharge> ions = ionic_point_charges(cell);
    if (v.nr != nspin || v.nc != rho_basis.nrxx)
    {
        v.create(nspin, rho_basis.nrxx);
    }
    ModuleBase::GlobalFunc::ZEROS(v.c, nspin * rho_basis.nrxx);
    this->pcc_result_valid_ = false;
    double energy_hartree = 0.0;
    if (this->parameters_.pcc_boundary == ModuleSccs::Boundary::Pcc0d)
    {
        this->pcc_geometry_ = ModuleSccs::pcc_geometry(cell.latvec, cell.lat0, 1.0e-10);
        this->pcc_geometry_.origin = ionic_system_center(cell, this->pcc_geometry_);
        const ModuleSccs::MultipoleMoments electronic_moments
            = ModuleSccs::reduced_pcc_density_moments(electronic_charge,
                                                      positions,
                                                      volume_element,
                                                      this->pcc_geometry_,
                                                      reduction);
        const ModuleSccs::MultipoleMoments ionic_moments
            = ModuleSccs::point_charge_moments(ions, this->pcc_geometry_);
        this->pcc_moments_ = add_moments(ionic_moments, electronic_moments);
        if (std::abs(this->pcc_moments_.charge
                     - (this->parameters_.expected_ionic_charge
                        - this->parameters_.expected_electron_count))
            > this->parameters_.normalization_tolerance)
        {
            throw std::runtime_error("standalone PCC charge does not match the requested electron count");
        }
        energy_hartree = ModuleSccs::pcc_self_energy(this->pcc_moments_,
                                                     this->pcc_geometry_.parameters);
        for (int ir = 0; ir < rho_basis.nrxx; ++ir)
        {
            const double electron_potential
                = -2.0 * ModuleSccs::pcc_potential(
                    this->pcc_moments_,
                    ModuleSccs::pcc_relative_position(positions[ir], this->pcc_geometry_),
                    this->pcc_geometry_.parameters);
            for (int spin = 0; spin < nspin; ++spin)
            {
                v(spin, ir) = electron_potential;
            }
        }
    }
    else
    {
        this->pcc_2d_geometry_ = ModuleSccs::pcc_2d_geometry(cell.latvec, cell.lat0, 1.0e-10);
        this->pcc_2d_geometry_.origin_y
            = ionic_system_center_y(cell, this->pcc_2d_geometry_.parameters.cell_length_y);
        const ModuleSccs::Pcc2dMoments electronic_moments
            = ModuleSccs::reduced_pcc_2d_density_moments(electronic_charge,
                                                         positions,
                                                         volume_element,
                                                         this->pcc_2d_geometry_,
                                                         reduction);
        const ModuleSccs::Pcc2dMoments ionic_moments
            = ModuleSccs::pcc_2d_point_charge_moments(ions, this->pcc_2d_geometry_);
        this->pcc_2d_moments_ = add_moments_2d(ionic_moments, electronic_moments);
        if (std::abs(this->pcc_2d_moments_.charge
                     - (this->parameters_.expected_ionic_charge
                        - this->parameters_.expected_electron_count))
            > this->parameters_.normalization_tolerance)
        {
            throw std::runtime_error("standalone PCC charge does not match the requested electron count");
        }
        energy_hartree = ModuleSccs::pcc_2d_self_energy(this->pcc_2d_moments_,
                                                        this->pcc_2d_geometry_.parameters);
        for (int ir = 0; ir < rho_basis.nrxx; ++ir)
        {
            const double electron_potential
                = -2.0 * ModuleSccs::pcc_2d_potential(
                    this->pcc_2d_moments_,
                    ModuleSccs::pcc_2d_relative_y(positions[ir].y, this->pcc_2d_geometry_),
                    this->pcc_2d_geometry_.parameters);
            for (int spin = 0; spin < nspin; ++spin)
            {
                v(spin, ir) = electron_potential;
            }
        }
    }
    surchem::Ael = 2.0 * energy_hartree;
    surchem::Acav = 0.0;
    this->pcc_result_valid_ = true;
}

void surchem::cal_force_pcc(const UnitCell& cell, ModuleBase::matrix& force) const
{
    if (!this->pcc_result_valid_)
    {
        throw std::logic_error("PCC force requires a current PCC potential");
    }
    int atom_index = 0;
    for (int atom_type = 0; atom_type < cell.ntype; ++atom_type)
    {
        for (int atom = 0; atom < cell.atoms[atom_type].na; ++atom)
        {
            ModuleSccs::PointCharge ion;
            ion.charge = cell.atoms[atom_type].ncpp.zv;
            ion.position = cell.atoms[atom_type].tau[atom] * cell.lat0;
            const ModuleBase::Vector3<double> correction
                = this->parameters_.pcc_boundary == ModuleSccs::Boundary::Pcc0d
                      ? ModuleSccs::pcc_point_charge_force(this->pcc_moments_,
                                                           ion,
                                                           this->pcc_geometry_)
                      : ModuleSccs::pcc_2d_point_charge_force(this->pcc_2d_moments_,
                                                              ion,
                                                              this->pcc_2d_geometry_);
            force(atom_index, 0) += 2.0 * correction.x;
            force(atom_index, 1) += 2.0 * correction.y;
            force(atom_index, 2) += 2.0 * correction.z;
            ++atom_index;
        }
    }
}
