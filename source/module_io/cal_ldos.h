#ifndef CAL_LDOS_H
#define CAL_LDOS_H

#include "module_elecstate/elecstate_pw.h"

namespace ModuleIO
{

void cal_ldos(const elecstate::ElecStatePW<std::complex<double>>* pelec,
              const psi::Psi<std::complex<double>>& psi,
              const Parallel_Grid& pgrid,
              const UnitCell& ucell);

} // namespace ModuleIO

#endif // CAL_LDOS_H