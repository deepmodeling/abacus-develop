#include "../surchem.h"
#include "sccs_pcc_coulomb.h"
#include "sccs_pcc_2d_coulomb.h"
#include "../sccs/sccs_pw_charge.h"
#include "../sccs/sccs_pw_reduction.h"
#include "source_base/timer_wrapper.h"

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
            const double y = cell.atoms[atom_type].tau[atom].y * cell.lat0;
            positions_y.push_back(y);
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
            const ModuleBase::Vector3<double> position = cell.atoms[atom_type].tau[atom] * cell.lat0;
            positions.push_back(position);
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

void surchem::v_correction_pcc(const UnitCell& cell,
                               const ModulePW::PW_Basis& rho_basis,
                               const int nspin,
                               const double* const* rho,
                               ModuleBase::matrix& v)
{
    const ModuleBase::TimePoint start_time = ModuleBase::get_time();
    if (!this->uses_pcc())
    {
        throw std::logic_error("PCC correction requires a PCC boundary");
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
    const int potential_size = nspin * rho_basis.nrxx;
    ModuleBase::GlobalFunc::ZEROS(v.c, potential_size);
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
        const double expected_charge = this->parameters_.expected_ionic_charge
                                       - this->parameters_.expected_electron_count;
        const double charge_error = this->pcc_moments_.charge - expected_charge;
        if (std::abs(charge_error) > this->parameters_.normalization_tolerance)
        {
            throw std::runtime_error("standalone PCC charge does not match the requested electron count");
        }
        energy_hartree = ModuleSccs::pcc_self_energy(this->pcc_moments_,
                                                     this->pcc_geometry_.parameters);
        for (int ir = 0; ir < rho_basis.nrxx; ++ir)
        {
            const ModuleBase::Vector3<double> relative_position
                = ModuleSccs::pcc_relative_position(positions[ir], this->pcc_geometry_);
            const double electron_potential
                = -2.0 * ModuleSccs::pcc_potential(
                    this->pcc_moments_,
                    relative_position,
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
        this->pcc_ionic_moments_2d_ = ionic_moments;
        this->pcc_2d_moments_ = add_moments_2d(ionic_moments, electronic_moments);
        const double expected_charge = this->parameters_.expected_ionic_charge
                                       - this->parameters_.expected_electron_count;
        const double charge_error = this->pcc_2d_moments_.charge - expected_charge;
        if (std::abs(charge_error) > this->parameters_.normalization_tolerance)
        {
            throw std::runtime_error("standalone PCC charge does not match the requested electron count");
        }
        energy_hartree = ModuleSccs::pcc_2d_self_energy(this->pcc_2d_moments_,
                                                        this->pcc_2d_geometry_.parameters);
        for (int ir = 0; ir < rho_basis.nrxx; ++ir)
        {
            const double relative_y
                = ModuleSccs::pcc_2d_relative_y(positions[ir].y, this->pcc_2d_geometry_);
            const double electron_potential
                = -2.0 * ModuleSccs::pcc_2d_potential(
                    this->pcc_2d_moments_,
                    relative_y,
                    this->pcc_2d_geometry_.parameters);
            for (int spin = 0; spin < nspin; ++spin)
            {
                v(spin, ir) = electron_potential;
            }
        }
    }
    this->pcc_energy_rydberg_ = 2.0 * energy_hartree;
    surchem::Ael = this->pcc_energy_rydberg_;
    surchem::Acav = 0.0;
    this->pcc_result_valid_ = true;
    const ModuleBase::TimePoint end_time = ModuleBase::get_time();
    this->pcc_elapsed_seconds_ = ModuleBase::get_duration(start_time, end_time);
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
