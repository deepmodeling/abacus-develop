#include "charge.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_elecstate/elecstate_getters.h"

#ifdef __MPI
void Charge::rho_mpi(const int& nbz, const int& bz)
{
    ModuleBase::TITLE("Charge", "rho_mpi");
    if (GlobalV::NPROC == 1)
        return;
    if (GlobalV::ESOLVER_TYPE == "sdft" && GlobalV::NPROC_IN_STOGROUP == 1)
        return; // qinarui add it temporarily.
    ModuleBase::timer::tick("Charge", "rho_mpi");
    int ir; // counters on real space mesh point.
    int iz; // counters on z direction of fft grid.
    int ip; // counters on processors

    //=========================================
    // There are two steps to do before getting
    // the final charge:
    // (1) sum up the plane rhos in each pool.
    // (2) sum up all the rhos from all pools.
    //=========================================

    //=================================================
    // Searching in all planes, and calculated each
    // plane belong to which processor.
    // Count number of planes for each cpu in this pool
    // num_z: how many planes on processor 'ip'
    //=================================================
    int* num_z = new int[GlobalV::NPROC_IN_POOL];
    ModuleBase::GlobalFunc::ZEROS(num_z, GlobalV::NPROC_IN_POOL);
    for (iz = 0; iz < nbz; iz++)
    {
        ip = iz % GlobalV::NPROC_IN_POOL;
        num_z[ip]++;
    }

    // mohan update 2011-04-26
    for (int ip = 0; ip < GlobalV::NPROC_IN_POOL; ip++)
    {
        num_z[ip] *= bz;
    }

    //=======================================
    // Find current number of planes (nz)
    // start_z: start position of z in
    // processor ip.
    //=======================================
    int* start_z = new int[GlobalV::NPROC_IN_POOL];
    ModuleBase::GlobalFunc::ZEROS(start_z, GlobalV::NPROC_IN_POOL);
    for (ip = 1; ip < GlobalV::NPROC_IN_POOL; ip++)
    {
        start_z[ip] = start_z[ip - 1] + num_z[ip - 1];
    }

    //====================================================
    // Find "number of data" in each processor in each pool
    //====================================================
    int* rec = new int[GlobalV::NPROC_IN_POOL];
    ModuleBase::GlobalFunc::ZEROS(rec, GlobalV::NPROC_IN_POOL);
    const int ncxy = this->rhopw->nx * this->rhopw->ny;
    for (ip = 0; ip < GlobalV::NPROC_IN_POOL; ip++)
    {
        rec[ip] = num_z[ip] * ncxy;
    }

    //======================================================
    // Find current "index of data" in each cpu in this pool
    // also, we mean start position of data.
    //======================================================
    int* dis = new int[GlobalV::NPROC_IN_POOL];
    ModuleBase::GlobalFunc::ZEROS(dis, GlobalV::NPROC_IN_POOL);
    for (ip = 1; ip < GlobalV::NPROC_IN_POOL; ip++)
    {
        dis[ip] = dis[ip - 1] + rec[ip - 1];
    }

    //==========================================
    // Collection of rho in each pool
    // ( according to different k distribution,
    // so the rho in each pool is different
    //==========================================
    double* rho_tmp = new double[this->rhopw->nrxx];
    double* rho_tot = new double[this->rhopw->nxyz];
    double* rho_tot_aux = new double[this->rhopw->nxyz];
    ModuleBase::GlobalFunc::ZEROS(rho_tot_aux, this->rhopw->nxyz);

    double* tau_tmp;
    double* tau_tot;
    double* tau_tot_aux;

    if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
    {
        tau_tmp = new double[this->rhopw->nrxx];
        tau_tot = new double[this->rhopw->nxyz];
        tau_tot_aux = new double[this->rhopw->nxyz];
        ModuleBase::GlobalFunc::ZEROS(tau_tot_aux, this->rhopw->nxyz);
    }

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(rho_tot, this->rhopw->nxyz);
        if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
            ModuleBase::GlobalFunc::ZEROS(tau_tot, this->rhopw->nxyz);

        for (ir = 0; ir < this->rhopw->nrxx; ir++)
        {
            rho_tmp[ir] = this->rho[is][ir] / static_cast<double>(GlobalV::NPROC_IN_POOL);
            if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
            {
                tau_tmp[ir] = this->kin_r[is][ir] / static_cast<double>(GlobalV::NPROC_IN_POOL);
            }
        }

        MPI_Allgatherv(rho_tmp, this->rhopw->nrxx, MPI_DOUBLE, rho_tot, rec, dis, MPI_DOUBLE, POOL_WORLD);
        if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
        {
            MPI_Allgatherv(tau_tmp, this->rhopw->nrxx, MPI_DOUBLE, tau_tot, rec, dis, MPI_DOUBLE, POOL_WORLD);
        }
        //=================================================================
        // Change the order of rho_tot in each pool , make them consistent
        // this is the most complicated part !!
        //=================================================================
        ModuleBase::GlobalFunc::ZEROS(rho_tot_aux, this->rhopw->nxyz);
        if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
        {
            ModuleBase::GlobalFunc::ZEROS(tau_tot_aux, this->rhopw->nxyz);
        }

        for (ip = 0; ip < GlobalV::NPROC_IN_POOL; ip++)
        {
            for (ir = 0; ir < ncxy; ir++)
            {
                for (iz = 0; iz < num_z[ip]; iz++)
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
                    // between [ start_z[ip], start_z[ip]+num_z[ip] )
                    // (4) start_z[ip] + iz yields correct z coordiante.
                    // -------------------------------------------------
                    // rot_tot: suitable for local pool.
                    // (1) the data save along z direction, only
                    // in a small distance.
                    // (2) however, the number of z in each processor
                    // 'ip' is num_z[ip]
                    // (3) the index of data increases with the ip,
                    // so the data on large 'ip' processor must
                    // have large 'start position', which we label
                    // start_z[ip] * ncxy.
                    // -------------------------------------------------
                    rho_tot_aux[this->rhopw->nz * ir + start_z[ip] + iz]
                        = rho_tot[num_z[ip] * ir + start_z[ip] * ncxy + iz];
                    if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
                    {
                        tau_tot_aux[this->rhopw->nz * ir + start_z[ip] + iz]
                            = tau_tot[num_z[ip] * ir + start_z[ip] * ncxy + iz];
                    }
                }
            }
        }
        //==================================
        // Reduce all the rho in each cpu
        //==================================
        if (GlobalV::ESOLVER_TYPE == "sdft") // qinarui add it temporarily.
        {
            MPI_Allreduce(rho_tot_aux, rho_tot, this->rhopw->nxyz, MPI_DOUBLE, MPI_SUM, STO_WORLD);
            if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
            {
                MPI_Allreduce(tau_tot_aux, tau_tot, this->rhopw->nxyz, MPI_DOUBLE, MPI_SUM, STO_WORLD);
            }
        }
        else
            MPI_Allreduce(rho_tot_aux, rho_tot, this->rhopw->nxyz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
        {
            MPI_Allreduce(tau_tot_aux, tau_tot, this->rhopw->nxyz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }

        //=====================================
        // Change the order of rho in each cpu
        //=====================================
        for (ir = 0; ir < ncxy; ir++)
        {
            for (iz = 0; iz < num_z[GlobalV::RANK_IN_POOL]; iz++)
            {
                this->rho[is][num_z[GlobalV::RANK_IN_POOL] * ir + iz]
                    = rho_tot[this->rhopw->nz * ir + this->rhopw->startz_current + iz];
                if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
                {
                    this->kin_r[is][num_z[GlobalV::RANK_IN_POOL] * ir + iz]
                        = tau_tot[this->rhopw->nz * ir + this->rhopw->startz_current + iz];
                }
            }
        }
    }
    delete[] rho_tot_aux;
    delete[] rho_tot;
    delete[] rho_tmp;

    if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
    {
        delete[] tau_tot_aux;
        delete[] tau_tot;
        delete[] tau_tmp;
    }
    delete[] rec;
    delete[] dis;

    delete[] num_z;
    delete[] start_z;
    ModuleBase::timer::tick("Charge", "rho_mpi");
    return;
}
#endif