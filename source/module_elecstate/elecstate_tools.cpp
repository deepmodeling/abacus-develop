#include "elecstate_tools.h"

namespace elecstate
{
    void calEBand(const ModuleBase::matrix& ekb,const ModuleBase::matrix& wg,fenergy& f_en)
    {
        ModuleBase::TITLE("ElecState", "calEBand");
        // calculate ebands using wg and ekb
        double eband = 0.0;
    #ifdef _OPENMP
    #pragma omp parallel for collapse(2) reduction(+ : eband)
    #endif
        for (int ik = 0; ik < ekb.nr; ++ik)
        {
            for (int ibnd = 0; ibnd < ekb.nc; ibnd++)
            {
                eband += ekb(ik, ibnd) * wg(ik, ibnd);
            }
        }
        f_en.eband = eband;

    #ifdef __MPI
        const int npool = GlobalV::KPAR * PARAM.inp.bndpar;
        Parallel_Reduce::reduce_double_allpool(npool, GlobalV::NPROC_IN_POOL, f_en.eband);
    #endif
        return;
    }
}