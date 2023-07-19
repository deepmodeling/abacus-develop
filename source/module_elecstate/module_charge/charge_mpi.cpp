#include "charge.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_elecstate/elecstate_getters.h"

#ifdef __MPI
void Charge::init_chgmpi(const int& nbz, const int& bz)
{
    if (GlobalV::NPROC_IN_STOGROUP % GlobalV::KPAR == 0)
    {
        this->use_intel_pool = true;
    }
    else
    {
        this->use_intel_pool = false;
        delete[] rec;
        rec = new int[GlobalV::NPROC_IN_POOL];
        delete[] dis;
        dis = new int[GlobalV::NPROC_IN_POOL];

        const int ncxy = this->rhopw->nx * this->rhopw->ny;
        for (int ip = 0; ip < GlobalV::NPROC_IN_POOL; ip++)
        {
            rec[ip] = this->rhopw->numz[ip] * ncxy;
            dis[ip] = this->rhopw->startz[ip] * ncxy;
        }
    }
}

void Charge::reduce_diff_pools(double* array_rho) const
{
    if (this->use_intel_pool)
    {
        MPI_Allreduce(MPI_IN_PLACE, array_rho, this->nrxx, MPI_DOUBLE, MPI_SUM, INTER_POOL);
    }
    else
    {
        double* array_tmp = new double[this->rhopw->nxyz];
        double* array_tot = new double[this->rhopw->nxyz];
        double* array_tot_aux = new double[this->rhopw->nxyz];
        //==================================
        // Collect the rho in each pool
        //==================================
        for (int ir = 0; ir < this->rhopw->nrxx; ++ir)
        {
            array_tmp[ir] = array_rho[ir] / GlobalV::NPROC_IN_POOL;
        }
        MPI_Allgatherv(array_tmp, this->rhopw->nrxx, MPI_DOUBLE, array_tot, rec, dis, MPI_DOUBLE, POOL_WORLD);

        const int ncxy = this->rhopw->nx * this->rhopw->ny;
        for (int ip = 0; ip < GlobalV::NPROC_IN_POOL; ++ip)
        {
            for (int ir = 0; ir < ncxy; ++ir)
            {
                for (int iz = 0; iz < this->rhopw->numz[ip]; ++iz)
                {
                    // -------------------------------------------------
                    // very carefully with the order of charge density.
                    // the data (ir,iz) is now in processor 'ip'.
                    // different POOL has different ordering.
                    // we want to collect them in each processor
                    // in a unit format,
                    // and then reduce among all POOLS to yield
                    // the correct charge density.
                    // we know the division of 'z' is indipendent
                    // in each processor, so the 'unit format'
                    // must have no relationship with 'z' divide method.
                    // -------------------------------------------------
                    // rot_tot_aux : suitable among all pools.
                    // (1) the data save along z direction.
                    // (2) and each element number of group 'z data'
                    // is 'this->rhopw->nz'
                    // (3) however, the data rearrange is occured
                    // between [ this->rhopw->startz[ip], this->rhopw->startz[ip]+this->rhopw->numz[ip] )
                    // (4) this->rhopw->startz[ip] + iz yields correct z coordiante.
                    // -------------------------------------------------
                    // rot_tot: suitable for local pool.
                    // (1) the data save along z direction, only
                    // in a small distance.
                    // (2) however, the number of z in each processor
                    // 'ip' is this->rhopw->numz[ip]
                    // (3) the index of data increases with the ip,
                    // so the data on large 'ip' processor must
                    // have large 'start position', which we label
                    // this->rhopw->startz[ip] * ncxy.
                    // -------------------------------------------------
                    array_tot_aux[this->rhopw->nz * ir + this->rhopw->startz[ip] + iz]
                        = array_tot[this->rhopw->numz[ip] * ir + this->rhopw->startz[ip] * ncxy + iz];
                }
            }
        }

        //==================================
        // Reduce all the rho in each cpu
        //==================================
        if (GlobalV::ESOLVER_TYPE == "sdft") // qinarui add it temporarily.
            MPI_Allreduce(array_tot_aux, array_tot, this->rhopw->nxyz, MPI_DOUBLE, MPI_SUM, STO_WORLD);
        else
            MPI_Allreduce(array_tot_aux, array_tot, this->rhopw->nxyz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        //=====================================
        // Change the order of rho in each cpu
        //=====================================
        for (int ir = 0; ir < ncxy; ir++)
        {
            for (int iz = 0; iz < this->rhopw->numz[GlobalV::RANK_IN_POOL]; iz++)
            {
                array_rho[this->rhopw->numz[GlobalV::RANK_IN_POOL] * ir + iz]
                    = array_tot[this->rhopw->nz * ir + this->rhopw->startz_current + iz];
            }
        }
        delete[] array_tot_aux;
        delete[] array_tot;
        delete[] array_tmp;
    }
}

void Charge::rho_mpi()
{
    ModuleBase::TITLE("Charge", "rho_mpi");
    if (GlobalV::NPROC == 1)
        return;
    if (GlobalV::ESOLVER_TYPE == "sdft" && GlobalV::NPROC_IN_STOGROUP == 1)
        return; // qinarui add it temporarily.
    ModuleBase::timer::tick("Charge", "rho_mpi");

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        reduce_diff_pools(this->rho[is]);
        if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
        {
            reduce_diff_pools(this->kin_r[is]);
        }
    }

    ModuleBase::timer::tick("Charge", "rho_mpi");
    return;
}
#endif