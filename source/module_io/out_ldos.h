#ifndef OUT_LDOS_H
#define OUT_LDOS_H

#include "module_elecstate/elecstate.h"
#include "module_io/cube_io.h"

namespace ModuleIO
{
template <typename T>
inline void out_ldos(const elecstate::ElecState* pelec,
                     const psi::Psi<T>& psi,
                     const Parallel_Grid& pgrid,
                     const UnitCell& ucell)
{
    std::vector<double> ldos(pelec->charge->nrxx);

    pelec->cal_ldos(psi, ldos);

    std::stringstream fn;
    fn << PARAM.globalv.global_out_dir << "LDOS_" << PARAM.inp.stm_bias << "eV"
       << ".cube";

    ModuleIO::write_vdata_palgrid(pgrid, ldos.data(), 0, PARAM.inp.nspin, 0, fn.str(), 0, &ucell, 11, 0);
}

} // namespace ModuleIO

#endif // OUT_LDOS_H
