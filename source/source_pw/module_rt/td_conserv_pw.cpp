#include "source_pw/module_rt/td_conserv_pw.h"

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"

namespace pw
{
TDConserv cal_td_conserv(double local_npsi,
                         const Charge& charge,
                         double omega,
                         int pool_rank,
                         double orth_before,
                         double orth_after,
                         bool full_gram)
{
    ModuleBase::timer::start("pw", "cal_td_conserv");
    TDConserv result;
    // Orbital norms are replicated within each k-point pool.
    if (pool_rank == 0)
    {
        result.npsi = local_npsi;
    }
    Parallel_Reduce::reduce_all(result.npsi);

    const int density_spins = charge.nspin == 2 ? 2 : 1;
    for (int spin = 0; spin < density_spins; ++spin)
    {
        for (int ir = 0; ir < charge.nrxx; ++ir)
        {
            result.nrho += charge.rho[spin][ir];
        }
    }
    result.nrho *= omega / charge.rhopw->nxyz;
    Parallel_Reduce::reduce_pool(result.nrho);

    double maxima[] = {orth_before, orth_after};
    if (full_gram)
    {
        Parallel_Reduce::reduce_max(maxima, 2);
    }
    result.orth_before = maxima[0];
    result.orth_after = maxima[1];
    ModuleBase::timer::end("pw", "cal_td_conserv");
    return result;
}
} // namespace pw
