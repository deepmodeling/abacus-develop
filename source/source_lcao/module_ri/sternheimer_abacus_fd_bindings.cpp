#include "source_basis/module_pw/pw_basis.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_lcao/module_ri/sternheimer_abacus_fd_adapter.h"
#include "source_lcao/module_ri/sternheimer_abacus_fd_nonlocal.h"

#include <stdexcept>
#include <utility>

#ifdef __MPI
#include <mpi.h>
#endif

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

SternheimerABACUSFDGridData make_sternheimer_fd_full_grid(const ModulePW::PW_Basis& pw_basis,
                                                          const double orthogonality_tolerance)
{
    return make_sternheimer_fd_grid_from_lattice(pw_basis.nx,
                                                 pw_basis.ny,
                                                 pw_basis.nz,
                                                 pw_basis.nxyz,
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

std::vector<double> copy_sternheimer_full_local_potential(const elecstate::Potential& potential,
                                                          const ModulePW::PW_Basis& pw_basis,
                                                          const int spin)
{
    const ModuleBase::matrix& veff = potential.get_eff_v();
    std::vector<double> local_potential;
    if (pw_basis.nrxx > 0)
    {
        if (spin < 0 || spin >= veff.nr)
        {
            throw std::invalid_argument("Sternheimer ABACUS FD full-grid potential spin index is out of range.");
        }
        if (veff.nc != pw_basis.nrxx)
        {
            throw std::invalid_argument("Sternheimer ABACUS FD full-grid potential size does not match local nrxx.");
        }

        const double* veff_spin = potential.get_eff_v(spin);
        if (veff_spin == nullptr)
        {
            throw std::invalid_argument("Sternheimer ABACUS FD full-grid potential is not allocated.");
        }
        local_potential.assign(veff_spin, veff_spin + pw_basis.nrxx);
    }
    std::vector<double> full_potential = embed_sternheimer_local_z_slab(local_potential,
                                                                        pw_basis.nxy,
                                                                        pw_basis.nz,
                                                                        pw_basis.nplane,
                                                                        pw_basis.startz_current);

#ifdef __MPI
    if (pw_basis.poolnproc > 1)
    {
        MPI_Allreduce(MPI_IN_PLACE,
                      full_potential.data(),
                      pw_basis.nxyz,
                      MPI_DOUBLE,
                      MPI_SUM,
                      pw_basis.pool_world);
    }
#endif

    return full_potential;
}

std::vector<double> copy_sternheimer_fixed_local_potential(const elecstate::Potential& potential,
                                                           const ModulePW::PW_Basis& pw_basis)
{
    const SternheimerABACUSFDGridData grid_data = make_sternheimer_fd_grid(pw_basis);
    if (grid_data.grid.size() != pw_basis.nrxx || pw_basis.nrxx <= 0)
    {
        throw std::invalid_argument(
            "Sternheimer ABACUS FD fixed potential size does not match the real-space grid.");
    }
    const double* fixed = potential.get_fixed_v();
    if (fixed == nullptr)
    {
        throw std::invalid_argument("Sternheimer ABACUS FD fixed potential is not allocated.");
    }
    return std::vector<double>(fixed, fixed + pw_basis.nrxx);
}

std::vector<double> copy_sternheimer_full_fixed_local_potential(
    const elecstate::Potential& potential,
    const ModulePW::PW_Basis& pw_basis)
{
    std::vector<double> local_fixed;
    if (pw_basis.nrxx > 0)
    {
        const double* fixed = potential.get_fixed_v();
        if (fixed == nullptr)
        {
            throw std::invalid_argument("Sternheimer ABACUS FD fixed potential is not allocated.");
        }
        local_fixed.assign(fixed, fixed + pw_basis.nrxx);
    }
    std::vector<double> full_fixed = embed_sternheimer_local_z_slab(local_fixed,
                                                                    pw_basis.nxy,
                                                                    pw_basis.nz,
                                                                    pw_basis.nplane,
                                                                    pw_basis.startz_current);
#ifdef __MPI
    if (pw_basis.poolnproc > 1)
    {
        MPI_Allreduce(MPI_IN_PLACE,
                      full_fixed.data(),
                      pw_basis.nxyz,
                      MPI_DOUBLE,
                      MPI_SUM,
                      pw_basis.pool_world);
    }
#endif
    return full_fixed;
}

SternheimerFDHamiltonian make_sternheimer_fd_hamiltonian(const elecstate::Potential& potential,
                                                         const ModulePW::PW_Basis& pw_basis,
                                                         const int spin,
                                                         const double kinetic_prefactor,
                                                         const int finite_difference_order)
{
    const SternheimerABACUSFDGridData grid_data = make_sternheimer_fd_grid(pw_basis);
    return make_sternheimer_fd_hamiltonian_from_local_potential(grid_data,
                                                                copy_sternheimer_local_potential(potential,
                                                                                                  pw_basis,
                                                                                                  spin),
                                                                kinetic_prefactor,
                                                                nullptr,
                                                                finite_difference_order);
}

SternheimerFDHamiltonian make_sternheimer_fd_hamiltonian(const elecstate::Potential& potential,
                                                         const ModulePW::PW_Basis& pw_basis,
                                                         const UnitCell& ucell,
                                                         const int spin,
                                                         const double kinetic_prefactor,
                                                         const SternheimerReducedKPoint& kpoint,
                                                         const int finite_difference_order)
{
    SternheimerABACUSFDGridData grid_data = make_sternheimer_fd_grid(pw_basis);
    grid_data.grid.kpoint = kpoint;
    auto nonlocal_projector
        = make_sternheimer_fd_nonlocal_projector_from_unitcell(ucell, grid_data.grid, grid_data.volume_element);
    return make_sternheimer_fd_hamiltonian_from_local_potential(grid_data,
                                                                copy_sternheimer_local_potential(potential,
                                                                                                  pw_basis,
                                                                                                  spin),
                                                                kinetic_prefactor,
                                                                std::move(nonlocal_projector),
                                                                finite_difference_order);
}

SternheimerFDHamiltonian make_sternheimer_fd_full_hamiltonian(const elecstate::Potential& potential,
                                                              const ModulePW::PW_Basis& pw_basis,
                                                              const UnitCell& ucell,
                                                              const int spin,
                                                              const double kinetic_prefactor,
                                                              const SternheimerReducedKPoint& kpoint,
                                                              const int finite_difference_order)
{
    SternheimerABACUSFDGridData grid_data = make_sternheimer_fd_full_grid(pw_basis);
    grid_data.grid.kpoint = kpoint;
    auto nonlocal_projector
        = make_sternheimer_fd_nonlocal_projector_from_unitcell(ucell, grid_data.grid, grid_data.volume_element);
    return make_sternheimer_fd_hamiltonian_from_local_potential(
        grid_data,
        copy_sternheimer_full_local_potential(potential, pw_basis, spin),
        kinetic_prefactor,
        std::move(nonlocal_projector),
        finite_difference_order);
}

} // namespace ModuleRI
