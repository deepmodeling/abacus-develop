#include "chg_parallel.h"

#ifdef __MPI

#include <cassert>

#include "charge.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_comm.h"
#include "source_base/timer.h"

namespace module_charge
{

void reduce_diff_pools(double* array_rho, const Charge& chr, const int kpar,
                       const bool all_ks_run, const int bndpar)
{
    ModuleBase::TITLE("Charge", "reduce_diff_pools");
    ModuleBase::timer::start("Charge", "reduce_diff_pools");
    // A rank may own zero real-space grid points (nrxx == 0); in that case
    // the buffer is legitimately null and the MPI calls below use count 0.
    // Only a null buffer with a non-zero nrxx is a genuine bug.
    assert(array_rho != nullptr || chr.nrxx == 0);
    assert(kpar >= 1);
    assert(bndpar >= 1);
    if (kpar > 1)
    {
        assert(chr.pgrid != nullptr);
        chr.pgrid->reduce_across_pools(array_rho);
    }
    if (all_ks_run && bndpar > 1)
    {
        // nrxx may be 0 on ranks with empty grid partitions; MPI_Allreduce
        // with count 0 is valid and ignores the buffer.
        MPI_Allreduce(MPI_IN_PLACE, array_rho, chr.nrxx, MPI_DOUBLE, MPI_SUM, BP_WORLD);
    }
    ModuleBase::timer::end("Charge", "reduce_diff_pools");
}

void rho_mpi(Charge& chr, const int kpar, const bool all_ks_run,
             const int bndpar, const int nspin)
{
    ModuleBase::TITLE("Charge", "rho_mpi");
    assert(kpar >= 1);
    assert(bndpar >= 1);
    assert(nspin > 0);
    if (kpar * bndpar <= 1)
    {
        return;
    }
    ModuleBase::timer::start("Charge", "rho_mpi");

    assert(chr.rho != nullptr);
    for (int is = 0; is < nspin; ++is)
    {
        reduce_diff_pools(chr.rho[is], chr, kpar, all_ks_run, bndpar);
        if (chr.kin_r != nullptr)
        {
            reduce_diff_pools(chr.kin_r[is], chr, kpar, all_ks_run, bndpar);
        }
    }

    ModuleBase::timer::end("Charge", "rho_mpi");
    return;
}

void kin_r_mpi(Charge& chr, const int kpar, const bool all_ks_run,
               const int bndpar, const int nspin)
{
    ModuleBase::TITLE("Charge", "kin_r_mpi");
    assert(kpar >= 1);
    assert(bndpar >= 1);
    assert(nspin > 0);
    if (kpar * bndpar <= 1)
    {
        return;
    }
    ModuleBase::timer::start("Charge", "kin_r_mpi");

    if (chr.kin_r != nullptr)
    {
        for (int is = 0; is < nspin; ++is)
        {
            reduce_diff_pools(chr.kin_r[is], chr, kpar, all_ks_run, bndpar);
        }
    }

    ModuleBase::timer::end("Charge", "kin_r_mpi");
    return;
}

} // namespace module_charge

#endif
