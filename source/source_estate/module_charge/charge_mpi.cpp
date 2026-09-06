#include "charge.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_comm.h"
#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"
#ifdef __MPI
void Charge::reduce_diff_pools(double* array_rho) const
{
    ModuleBase::TITLE("Charge", "reduce_diff_pools");
    ModuleBase::timer::start("Charge", "reduce_diff_pools");
    if (GlobalV::KPAR > 1)
    {
        assert(this->pgrid != nullptr);
        this->pgrid->reduce_across_pools(array_rho);
    }
    if (PARAM.globalv.all_ks_run && PARAM.inp.bndpar > 1)
    {
        MPI_Allreduce(MPI_IN_PLACE, array_rho, this->nrxx, MPI_DOUBLE, MPI_SUM, BP_WORLD);
    }
    ModuleBase::timer::end("Charge", "reduce_diff_pools");
}

void Charge::rho_mpi()
{
    ModuleBase::TITLE("Charge", "rho_mpi");
    if (GlobalV::KPAR * PARAM.inp.bndpar <= 1)
    {
        return;
    }
    ModuleBase::timer::start("Charge", "rho_mpi");

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        reduce_diff_pools(this->rho[is]);
        if (XC_Functional::get_ked_flag() || PARAM.inp.out_elf[0] > 0)
        {
            reduce_diff_pools(this->kin_r[is]);
        }
    }

    ModuleBase::timer::end("Charge", "rho_mpi");
    return;
}

void Charge::kin_r_mpi()
{
    ModuleBase::TITLE("Charge", "kin_r_mpi");
    if (GlobalV::KPAR * PARAM.inp.bndpar <= 1)
    {
        return;
    }
    ModuleBase::timer::start("Charge", "kin_r_mpi");

    if (XC_Functional::get_ked_flag() || PARAM.inp.out_elf[0] > 0)
    {
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            reduce_diff_pools(this->kin_r[is]);
        }
    }

    ModuleBase::timer::end("Charge", "kin_r_mpi");
    return;
}
#endif
