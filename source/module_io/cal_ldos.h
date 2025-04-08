#ifndef CAL_LDOS_H
#define CAL_LDOS_H

#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/elecstate_pw.h"

namespace ModuleIO
{
template <typename T>
class Cal_ldos
{
  public:
    Cal_ldos(){};
    ~Cal_ldos(){};

    static void cal_ldos_pw(const elecstate::ElecStatePW<std::complex<double>>* pelec,
                            const psi::Psi<std::complex<double>>& psi,
                            const Parallel_Grid& pgrid,
                            const UnitCell& ucell);

    static void cal_ldos_lcao(const elecstate::ElecStateLCAO<T>* pelec,
                              const psi::Psi<T>& psi,
                              const Parallel_Grid& pgrid,
                              const UnitCell& ucell);
}; // namespace Cal_ldos
} // namespace ModuleIO

#endif // CAL_LDOS_H