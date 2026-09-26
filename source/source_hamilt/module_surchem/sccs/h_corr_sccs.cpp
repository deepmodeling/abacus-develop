#include "../surchem.h"

#include "../pcc/sccs_pcc_2d_coulomb.h"
#include "../pcc/sccs_pcc_coulomb.h"
#include "sccs_pw_charge.h"
#include "sccs_pw_reduction.h"

#include "source_base/timer.h"
#include "source_base/timer_wrapper.h"
#include "source_base/tool_title.h"

#include <cmath>
#include <stdexcept>
#include <vector>

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
    ModuleBase::matrix pcc_potential;
    double vacuum_pcc_energy = 0.0;
    if (this->uses_pcc())
    {
        this->v_correction_pcc(cell, rho_basis, nspin, rho, pcc_potential);
        vacuum_pcc_energy = 0.5 * this->pcc_energy_rydberg_;
        if (this->parameters_.pcc_boundary == ModuleSccs::Boundary::Pcc0d)
        {
            origin = this->pcc_geometry_.origin;
        }
    }
    const ModuleSccs::PccGeometry& pcc_geometry = this->pcc_geometry_;
    const ModuleSccs::Pcc2dGeometry& pcc_2d_geometry = this->pcc_2d_geometry_;

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

    if (this->uses_pcc())
    {
        this->sccs_result_.point_solute_moments = this->pcc_moments_;
        this->sccs_result_.point_solute_moments_2d = this->pcc_2d_moments_;
        this->sccs_result_.vacuum_pcc_energy = vacuum_pcc_energy;
        for (int ir = 0; ir < rho_basis.nrxx; ++ir)
        {
            this->sccs_result_.electron_potential_hartree[ir] += 0.5 * pcc_potential(0, ir);
        }
        if (this->parameters_.pcc_boundary == ModuleSccs::Boundary::Pcc2d)
        {
            const ModuleSccs::Pcc2dMoments smooth_ionic_moments
                = ModuleSccs::reduced_pcc_2d_density_moments(this->sccs_result_.charge.ionic,
                                                            positions,
                                                            volume_element,
                                                            pcc_2d_geometry,
                                                            charge_reduction);
            this->sccs_result_.ionic_shape_pcc_energy
                = ModuleSccs::pcc_2d_ionic_shape_energy(this->sccs_result_.polarization_moments_2d.charge,
                                                       smooth_ionic_moments,
                                                       this->pcc_ionic_moments_2d_,
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
