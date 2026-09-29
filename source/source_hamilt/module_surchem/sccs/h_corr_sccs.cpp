#include "../surchem.h"

#include "../common/charge_reduction.h"
#include "../common/pw_grid.h"
#include "sccs_gaussian_ion.h"

#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

// Adapt the total spin density and local-pseudopotential ionic source to SCCS.
// The reaction energy uses smooth ions; standalone PCC supplies the point-ion
// vacuum correction once. Convert the final Ha quantities to ABACUS Ry units.
void surchem::v_correction_sccs(const UnitCell& cell,
                                const ModulePW::PW_Basis& rho_basis,
                                const int nspin,
                                const double* const* rho,
                                const double* vlocal,
                                ModuleBase::matrix& v)
{
    ModuleBase::TITLE("surchem", "v_correction_sccs");
    ModuleBase::timer::start("surchem", "v_correction_sccs");
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
        const double* spin_density_end = rho[spin] + rho_basis.nrxx;
        spin_density[spin].assign(rho[spin], spin_density_end);
    }
    const std::vector<double> electron_density
        = ModuleSccs::sum_electron_density(spin_density, nspin);
    FixedSourceCache& cache = this->fixed_source_cache_;
    const double* local_potential_end = vlocal + rho_basis.nrxx;
    const bool same_lattice
        = cache.valid && cache.lattice_vectors.e11 == cell.latvec.e11
          && cache.lattice_vectors.e12 == cell.latvec.e12
          && cache.lattice_vectors.e13 == cell.latvec.e13
          && cache.lattice_vectors.e21 == cell.latvec.e21
          && cache.lattice_vectors.e22 == cell.latvec.e22
          && cache.lattice_vectors.e23 == cell.latvec.e23
          && cache.lattice_vectors.e31 == cell.latvec.e31
          && cache.lattice_vectors.e32 == cell.latvec.e32
          && cache.lattice_vectors.e33 == cell.latvec.e33;
    const bool reuse_fixed_sources
        = cache.valid && cache.basis == &rho_basis
          && cache.nx == rho_basis.nx && cache.ny == rho_basis.ny
          && cache.nz == rho_basis.nz && cache.nrxx == rho_basis.nrxx
          && cache.nplane == rho_basis.nplane
          && cache.startz == rho_basis.startz_current
          && cache.lattice_constant == cell.lat0
          && cache.cell_volume == cell.omega && cache.tpiba == cell.tpiba
          && cache.ionic_charge == this->parameters_.expected_ionic_charge
          && same_lattice
          && cache.local_potential.size() == static_cast<std::size_t>(rho_basis.nrxx)
          && std::equal(vlocal, local_potential_end, cache.local_potential.begin());
    if (!reuse_fixed_sources)
    {
        cache.valid = false;
        cache.local_potential.assign(vlocal, local_potential_end);
        cache.ionic_density
            = ModuleSccs::gaussian_ionic_density(cell, rho_basis, ModuleSccs::gaussian_ion_spread);
        cache.core_density.clear();
        if (this->parameters_.sccs_config.core_electrons)
        {
            cache.core_density = ModuleSccs::gaussian_core_density(
                cell, rho_basis, this->parameters_.sccs_config.core_spread);
        }
        cache.positions = ModuleSurchem::pw_grid_positions(rho_basis, cell.latvec, cell.lat0);
        cache.lattice_vectors = cell.latvec;
        cache.basis = &rho_basis;
        cache.lattice_constant = cell.lat0;
        cache.cell_volume = cell.omega;
        cache.tpiba = cell.tpiba;
        cache.ionic_charge = this->parameters_.expected_ionic_charge;
        cache.nx = rho_basis.nx;
        cache.ny = rho_basis.ny;
        cache.nz = rho_basis.nz;
        cache.nrxx = rho_basis.nrxx;
        cache.nplane = rho_basis.nplane;
        cache.startz = rho_basis.startz_current;
        cache.valid = true;
    }
    const std::vector<double>& ionic_density = cache.ionic_density;
    const std::vector<double>& core_density = cache.core_density;
    const std::vector<ModuleBase::Vector3<double>>& positions = cache.positions;
    ModuleBase::matrix pcc_potential;
    double vacuum_pcc_energy = 0.0;
    if (this->uses_pcc())
    {
        this->v_correction_pcc(cell, rho_basis, nspin, rho, pcc_potential);
        vacuum_pcc_energy = 0.5 * this->pcc_energy_rydberg_;
    }
    else
    {
        surchem::Epcc = 0.0;
    }
    const ModulePcc::PccGeometry& pcc_geometry = this->pcc_geometry_;
    const ModulePcc::Pcc2dGeometry& pcc_2d_geometry = this->pcc_2d_geometry_;

    const double volume_element = cell.omega / static_cast<double>(rho_basis.nxyz);
    const ModuleSurchem::PoolChargeReduction reduction(this->parameters_.pool_process_count);
    this->sccs_result_
        = ModuleSccs::evaluate_pw_sccs(electron_density,
                                       ionic_density,
                                       core_density,
                                       this->parameters_.expected_electron_count,
                                       this->parameters_.expected_ionic_charge,
                                       this->parameters_.normalization_tolerance,
                                       positions,
                                       this->parameters_.sccs_config,
                                       pcc_geometry,
                                       pcc_2d_geometry,
                                       rho_basis,
                                       cell.tpiba,
                                       volume_element,
                                       reduction,
                                       this->sccs_state_);
    this->sccs_result_.reused_fixed_sources = reuse_fixed_sources;

    if (this->uses_pcc())
    {
        this->sccs_result_.point_solute_moments = this->pcc_moments_;
        this->sccs_result_.point_solute_moments_2d = this->pcc_2d_moments_;
        this->sccs_result_.vacuum_pcc_energy = vacuum_pcc_energy;
        for (int ir = 0; ir < rho_basis.nrxx; ++ir)
        {
            this->sccs_result_.electron_potential_hartree[ir] += 0.5 * pcc_potential(0, ir);
        }

    }

    if (v.nr != nspin || v.nc != rho_basis.nrxx)
    {
        v.create(nspin, rho_basis.nrxx);
    }
    for (int spin = 0; spin < nspin; ++spin)
    {
        for (int ir = 0; ir < rho_basis.nrxx; ++ir)
        {
            v(spin, ir) = 2.0 * this->sccs_result_.electron_potential_hartree[ir];
        }
    }

    // The vacuum PCC energy stays in surchem::Epcc, set by v_correction_pcc.
    surchem::Ael = 2.0 * this->sccs_result_.electrostatic.reaction_energy;
    surchem::Acav = 2.0 * (this->sccs_result_.non_electrostatic.surface_energy
                           + this->sccs_result_.non_electrostatic.volume_energy);
    ModuleBase::timer::end("surchem", "v_correction_sccs");
}

void surchem::v_correction_solvent(const UnitCell& cell,
                                   const ModulePW::PW_Basis& rho_basis,
                                   const int nspin,
                                   const double* const* rho,
                                   const double* vlocal,
                                   ModuleBase::matrix& v)
{
    try
    {
        if (this->sccs_is_active())
        {
            this->v_correction_sccs(cell, rho_basis, nspin, rho, vlocal, v);
        }
        else if (this->uses_pcc())
        {
            this->v_correction_pcc(cell, rho_basis, nspin, rho, v);
        }
        else
        {
            // Delayed SCCS without PCC contributes nothing before activation.
            if (v.nr != nspin || v.nc != rho_basis.nrxx)
            {
                v.create(nspin, rho_basis.nrxx);
            }
            ModuleBase::GlobalFunc::ZEROS(v.c, nspin * rho_basis.nrxx);
            surchem::Ael = 0.0;
            surchem::Acav = 0.0;
            surchem::Epcc = 0.0;
        }
    }
    catch (const std::exception& error)
    {
        ModuleBase::WARNING_QUIT("surchem::v_correction", error.what());
    }
}
