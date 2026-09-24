#ifndef CAL_DM_H
#define CAL_DM_H

#include "source_estate/module_dm/dm_from_psi.h"
#include "source_base/timer.h"
#include "source_base/matrix.h"
#include "source_base/complexmatrix.h"

namespace elecstate
{

// for gamma_only(double case) and multi-k(complex<double> case)
// Thin adapter over module_dm::dmk_from_psi for callers that own the DM blocks as
// ModuleBase::matrix / ModuleBase::ComplexMatrix storage (e.g. DeePKS bandgap terms).
inline void cal_dm(const Parallel_Orbitals* ParaV, const ModuleBase::matrix& wg, const psi::Psi<double>& wfc, std::vector<ModuleBase::matrix>& dm)
{
    ModuleBase::TITLE("elecstate", "cal_dm");
    ModuleBase::timer::start("elecstate","cal_dm");

    for (int ik = 0; ik < wfc.get_nk(); ++ik)
    {
        dm[ik].create(ParaV->ncol, ParaV->nrow);
        module_dm::dmk_from_psi(ParaV, wg, ik, wfc, dm[ik].c);
    }
    ModuleBase::timer::end("elecstate","cal_dm");
}

inline void cal_dm(const Parallel_Orbitals* ParaV, const ModuleBase::matrix& wg, const psi::Psi<std::complex<double>>& wfc, std::vector<ModuleBase::ComplexMatrix>& dm)
{
    ModuleBase::TITLE("elecstate", "cal_dm");
    ModuleBase::timer::start("elecstate","cal_dm");

    for (int ik = 0; ik < wfc.get_nk(); ++ik)
    {
        dm[ik].create(ParaV->ncol, ParaV->nrow);
        module_dm::dmk_from_psi(ParaV, wg, ik, wfc, dm[ik].c);
    }

    ModuleBase::timer::end("elecstate","cal_dm");
}

}//namespace elecstate

#endif
