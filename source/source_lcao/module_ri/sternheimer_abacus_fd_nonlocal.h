#ifndef STERNHEIMER_ABACUS_FD_NONLOCAL_H
#define STERNHEIMER_ABACUS_FD_NONLOCAL_H

#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"
#include "source_lcao/module_ri/sternheimer_fd_nonlocal_projector.h"

#include <memory>

class UnitCell;

namespace ModuleRI
{

std::shared_ptr<SternheimerFDNonlocalProjector>
make_sternheimer_fd_nonlocal_projector_from_unitcell(const UnitCell& ucell,
                                                     const SternheimerFDHamiltonian::Grid& grid,
                                                     double volume_element);

} // namespace ModuleRI

#endif
