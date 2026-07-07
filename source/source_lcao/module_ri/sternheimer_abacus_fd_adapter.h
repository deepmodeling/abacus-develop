#ifndef STERNHEIMER_ABACUS_FD_ADAPTER_H
#define STERNHEIMER_ABACUS_FD_ADAPTER_H

#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"

#include <memory>
#include <vector>

class UnitCell;

namespace ModuleBase
{
class Matrix3;
}

namespace ModulePW
{
class PW_Basis;
}

namespace elecstate
{
class Potential;
}

namespace ModuleRI
{

struct SternheimerABACUSFDGridData
{
    SternheimerFDHamiltonian::Grid grid;
    double volume_element = 0.0;
};

SternheimerABACUSFDGridData make_sternheimer_fd_grid_from_lattice(int nx,
                                                                  int ny,
                                                                  int nz,
                                                                  int nrxx,
                                                                  double lat0,
                                                                  const ModuleBase::Matrix3& latvec,
                                                                  double orthogonality_tolerance = 1.0e-12);

SternheimerABACUSFDGridData make_sternheimer_fd_grid(const ModulePW::PW_Basis& pw_basis,
                                                     double orthogonality_tolerance = 1.0e-12);

SternheimerABACUSFDGridData make_sternheimer_fd_full_grid(const ModulePW::PW_Basis& pw_basis,
                                                          double orthogonality_tolerance = 1.0e-12);

std::vector<double> copy_sternheimer_local_potential(const elecstate::Potential& potential,
                                                     const ModulePW::PW_Basis& pw_basis,
                                                     int spin);

std::vector<double> copy_sternheimer_full_local_potential(const elecstate::Potential& potential,
                                                          const ModulePW::PW_Basis& pw_basis,
                                                          int spin);

SternheimerFDHamiltonian make_sternheimer_fd_hamiltonian_from_local_potential(
    const SternheimerABACUSFDGridData& grid_data,
    std::vector<double> local_potential,
    double kinetic_prefactor = 1.0,
    std::shared_ptr<const SternheimerFDNonlocalProjector> nonlocal_projector = nullptr);

SternheimerFDHamiltonian make_sternheimer_fd_hamiltonian(const elecstate::Potential& potential,
                                                         const ModulePW::PW_Basis& pw_basis,
                                                         int spin,
                                                         double kinetic_prefactor = 1.0);

SternheimerFDHamiltonian make_sternheimer_fd_hamiltonian(const elecstate::Potential& potential,
                                                         const ModulePW::PW_Basis& pw_basis,
                                                         const UnitCell& ucell,
                                                         int spin,
                                                         double kinetic_prefactor = 1.0);

SternheimerFDHamiltonian make_sternheimer_fd_full_hamiltonian(const elecstate::Potential& potential,
                                                              const ModulePW::PW_Basis& pw_basis,
                                                              const UnitCell& ucell,
                                                              int spin,
                                                              double kinetic_prefactor = 1.0);

} // namespace ModuleRI

#endif
