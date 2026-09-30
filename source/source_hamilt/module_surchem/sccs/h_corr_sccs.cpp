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

namespace
{

bool same_lattice(const ModuleBase::Matrix3& left, const ModuleBase::Matrix3& right)
{
    return left.e11 == right.e11 && left.e12 == right.e12 && left.e13 == right.e13
           && left.e21 == right.e21 && left.e22 == right.e22 && left.e23 == right.e23
           && left.e31 == right.e31 && left.e32 == right.e32 && left.e33 == right.e33;
}

} // namespace

bool surchem::FixedSourceCache::matches(const UnitCell& cell,
                                        const ModulePW::PW_Basis& rho_basis,
                                        const double* vlocal,
                                        const double valence_charge) const
{
    if (!this->valid || this->basis != &rho_basis
        || this->local_potential.size() != static_cast<std::size_t>(rho_basis.nrxx))
    {
        return false;
    }
    const bool same_grid = this->nx == rho_basis.nx && this->ny == rho_basis.ny
                           && this->nz == rho_basis.nz && this->nrxx == rho_basis.nrxx
                           && this->nplane == rho_basis.nplane
                           && this->startz == rho_basis.startz_current;
    const bool same_cell = this->lattice_constant == cell.lat0 && this->cell_volume == cell.omega
                           && this->tpiba == cell.tpiba
                           && same_lattice(this->lattice_vectors, cell.latvec);
    const double* vlocal_end = vlocal + rho_basis.nrxx;
    return same_grid && same_cell && this->ionic_charge == valence_charge
           && std::equal(vlocal, vlocal_end, this->local_potential.begin());
}

bool surchem::update_fixed_sources(const UnitCell& cell,
                                   const ModulePW::PW_Basis& rho_basis,
                                   const double* vlocal)
{
    FixedSourceCache& cache = this->fixed_source_cache_;
    const double valence_charge = this->parameters_.expected_ionic_charge;
    if (cache.matches(cell, rho_basis, vlocal, valence_charge))
    {
        return true;
    }
    const double* vlocal_end = vlocal + rho_basis.nrxx;
    cache.valid = false;
    cache.local_potential.assign(vlocal, vlocal_end);
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
    cache.ionic_charge = valence_charge;
    cache.nx = rho_basis.nx;
    cache.ny = rho_basis.ny;
    cache.nz = rho_basis.nz;
    cache.nrxx = rho_basis.nrxx;
    cache.nplane = rho_basis.nplane;
    cache.startz = rho_basis.startz_current;
    cache.valid = true;
    return false;
}

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
    const bool reuse_fixed_sources = this->update_fixed_sources(cell, rho_basis, vlocal);
    const FixedSourceCache& cache = this->fixed_source_cache_;
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

    // The electrostatic output potential adds the reaction potential of the
    // solute charge to the vacuum PCC term; the cavity derivatives are not
    // electrostatic and stay out of it.
    const std::vector<double>& reaction_potential = this->sccs_result_.electrostatic.reaction_potential;
    this->electrostatic_correction_ry_.resize(rho_basis.nrxx);
    for (int ir = 0; ir < rho_basis.nrxx; ++ir)
    {
        const double pcc_part = this->uses_pcc() ? pcc_potential(0, ir) : 0.0;
        this->electrostatic_correction_ry_[ir] = pcc_part - 2.0 * reaction_potential[ir];
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
            this->electrostatic_correction_ry_.assign(rho_basis.nrxx, 0.0);
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
