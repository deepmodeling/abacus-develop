#include "source_basis/module_pw/pw_basis.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_lcao/module_ri/sternheimer_abacus_fd_adapter.h"
#include "source_lcao/module_ri/sternheimer_abacus_fd_nonlocal.h"

#include <stdexcept>
#include <utility>

namespace ModuleRI
{

SternheimerABACUSFDGridData make_sternheimer_fd_grid(const ModulePW::PW_Basis& pw_basis,
                                                     const double orthogonality_tolerance)
{
    return make_sternheimer_fd_grid_from_lattice(pw_basis.nx,
                                                 pw_basis.ny,
                                                 pw_basis.nz,
                                                 pw_basis.nrxx,
                                                 pw_basis.lat0,
                                                 pw_basis.latvec,
                                                 orthogonality_tolerance);
}

std::vector<double> copy_sternheimer_local_potential(const elecstate::Potential& potential,
                                                     const ModulePW::PW_Basis& pw_basis,
                                                     const int spin)
{
    const SternheimerABACUSFDGridData grid_data = make_sternheimer_fd_grid(pw_basis);
    const int grid_size = grid_data.grid.size();
    const ModuleBase::matrix& veff = potential.get_eff_v();

    if (spin < 0 || spin >= veff.nr)
    {
        throw std::invalid_argument("Sternheimer ABACUS FD potential spin index is out of range.");
    }
    if (veff.nc != grid_size)
    {
        throw std::invalid_argument("Sternheimer ABACUS FD potential size does not match the real-space grid.");
    }

    const double* veff_spin = potential.get_eff_v(spin);
    if (veff_spin == nullptr)
    {
        throw std::invalid_argument("Sternheimer ABACUS FD potential is not allocated.");
    }
    return std::vector<double>(veff_spin, veff_spin + grid_size);
}

SternheimerFDHamiltonian make_sternheimer_fd_hamiltonian(const elecstate::Potential& potential,
                                                         const ModulePW::PW_Basis& pw_basis,
                                                         const int spin,
                                                         const double kinetic_prefactor)
{
    const SternheimerABACUSFDGridData grid_data = make_sternheimer_fd_grid(pw_basis);
    return make_sternheimer_fd_hamiltonian_from_local_potential(grid_data,
                                                                copy_sternheimer_local_potential(potential,
                                                                                                  pw_basis,
                                                                                                  spin),
                                                                kinetic_prefactor);
}

SternheimerFDHamiltonian make_sternheimer_fd_hamiltonian(const elecstate::Potential& potential,
                                                         const ModulePW::PW_Basis& pw_basis,
                                                         const UnitCell& ucell,
                                                         const int spin,
                                                         const double kinetic_prefactor)
{
    const SternheimerABACUSFDGridData grid_data = make_sternheimer_fd_grid(pw_basis);
    auto nonlocal_projector
        = make_sternheimer_fd_nonlocal_projector_from_unitcell(ucell, grid_data.grid, grid_data.volume_element);
    return make_sternheimer_fd_hamiltonian_from_local_potential(grid_data,
                                                                copy_sternheimer_local_potential(potential,
                                                                                                  pw_basis,
                                                                                                  spin),
                                                                kinetic_prefactor,
                                                                std::move(nonlocal_projector));
}

} // namespace ModuleRI
