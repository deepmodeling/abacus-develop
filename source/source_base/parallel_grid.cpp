#include "parallel_grid.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"

#ifdef __MPI
#include "source_base/parallel_comm.h" // use POOL_WORLD

#include <mpi.h>
#endif

Parallel_Grid::Parallel_Grid()
{
}

Parallel_Grid::~Parallel_Grid()
{
}

void Parallel_Grid::init(const int& ncx_in,
                         const int& ncy_in,
                         const int& ncz_in,
                         const int& nczp_in,
                         const int& nrxx_in,
                         const int& nbz_in,
                         const int& bz_in,
                         const int nprocgroup)
{

    ModuleBase::TITLE("Parallel_Grid", "init");

    this->ncx = ncx_in;
    this->ncy = ncy_in;
    this->ncz = ncz_in;
    this->nczp = nczp_in;
    this->nrxx = nrxx_in;
    this->nbz = nbz_in;
    this->bz = bz_in;

    if (nczp < 0)
    {
        GlobalV::ofs_warning << " nczp = " << nczp << std::endl;
        ModuleBase::WARNING_QUIT("Parallel_Grid::init", "nczp<0");
    }

    assert(ncx > 0);
    assert(ncy > 0);
    assert(ncz > 0);

    this->ncxy = ncx * ncy;
    this->ncxyz = ncxy * ncz;

#ifndef __MPI
    return;
#endif

    // enable to call this function again liuyu 2023-03-10
    if (!this->numz.empty())
    {
        this->nproc_in_pool.clear();
        this->numz.clear();
        this->startz.clear();
        this->whichpro.clear();
        this->whichpro_loc.clear();
    }

    // (2)
    assert(this->numz.empty());
    assert(GlobalV::KPAR > 0);

    this->nproc_in_pool.resize(GlobalV::KPAR);

    const int remain_pro = nprocgroup % GlobalV::KPAR;
    for (int i = 0; i < GlobalV::KPAR; i++)
    {
        nproc_in_pool[i] = nprocgroup / GlobalV::KPAR;
        if (i < remain_pro)
        {
            this->nproc_in_pool[i]++;
        }
    }

    this->numz.resize(GlobalV::KPAR);
    this->startz.resize(GlobalV::KPAR);
    this->whichpro.resize(GlobalV::KPAR);
    this->whichpro_loc.resize(GlobalV::KPAR);

    for (int ip = 0; ip < GlobalV::KPAR; ip++)
    {
        const int nproc = nproc_in_pool[ip];
        this->numz[ip].assign(nproc, 0);
        this->startz[ip].assign(nproc, 0);
        this->whichpro[ip].assign(this->ncz, 0);
        this->whichpro_loc[ip].assign(this->ncz, 0);
    }

    this->z_distribution();

    return;
}

void Parallel_Grid::z_distribution()
{
    assert(!this->numz.empty());

    int* startp = new int[GlobalV::KPAR];
    startp[0] = 0;
    for (int ip = 0; ip < GlobalV::KPAR; ip++)
    {
        //		GlobalV::ofs_running << "\n now POOL=" << ip;
        const int nproc = nproc_in_pool[ip];

        if (ip > 0)
        {
            startp[ip] = startp[ip - 1] + nproc_in_pool[ip - 1];
        }

        // (1) how many z on each 'proc' in each 'pool'
        for (int iz = 0; iz < nbz; iz++)
        {
            const int proc = iz % nproc;
            numz[ip][proc] += bz;
        }

        //		for(int proc=0; proc<nproc; proc++)
        //		{
        //			GlobalV::ofs_running << "\n proc=" << proc << " numz=" << numz[ip][proc];
        //		}

        // (2) start position of z in each 'proc' in each 'pool'
        startz[ip][0] = 0;
        for (int proc = 1; proc < nproc; proc++)
        {
            startz[ip][proc] = startz[ip][proc - 1] + numz[ip][proc - 1];
        }

        //		for(int proc=0; proc<nproc; proc++)
        //		{
        //			GlobalV::ofs_running << "\n proc=" << proc << " startz=" << startz[ip][proc];
        //		}

        // (3) each z belongs to which 'proc' ( global index )
        for (int iz = 0; iz < ncz; iz++)
        {
            for (int proc = 0; proc < nproc; proc++)
            {
                if (iz >= startz[ip][nproc - 1])
                {
                    whichpro[ip][iz] = startp[ip] + nproc - 1;
                    whichpro_loc[ip][iz] = nproc - 1;
                    break;
                }
                else if (iz >= startz[ip][proc] && iz < startz[ip][proc + 1])
                {
                    whichpro[ip][iz] = startp[ip] + proc;
                    whichpro_loc[ip][iz] = proc;
                    break;
                }
            }
        }

        //		for(int iz=0; iz<ncz; iz++)
        //		{
        //			GlobalV::ofs_running << "\n iz=" << iz << " whichpro=" << whichpro[ip][iz];
        //		}
    }

    delete[] startp;
    return;
}

void Parallel_Grid::reduce_across_pools(double* data) const
{
#ifdef __MPI
    if (GlobalV::KPAR <= 1)
    {
        return;
    }

    assert(data != nullptr);
    if (KP_WORLD != MPI_COMM_NULL)
    {
        // Equal-sized pools give corresponding ranks identical z-slab layouts,
        // so their local buffers can be summed directly without redistribution.
        MPI_Allreduce(MPI_IN_PLACE, data, this->nrxx, MPI_DOUBLE, MPI_SUM, KP_WORLD);
        return;
    }

    // Uneven pool sizes have no KP_WORLD and may assign different z-slabs to
    // corresponding ranks. Validate the local distribution before rebuilding
    // a common global layout for the cross-pool reduction.
    assert(!this->numz.empty());
    assert(GlobalV::MY_POOL >= 0 && GlobalV::MY_POOL < static_cast<int>(this->numz.size()));
    assert(GlobalV::RANK_IN_POOL >= 0 && GlobalV::RANK_IN_POOL < static_cast<int>(this->numz[GlobalV::MY_POOL].size()));
    assert(this->nczp == this->numz[GlobalV::MY_POOL][GlobalV::RANK_IN_POOL]);
    assert(this->nrxx == this->ncxy * this->nczp);

    const int pool_size = this->nproc_in_pool[GlobalV::MY_POOL];
    std::vector<int> receive_counts(pool_size);
    std::vector<int> displacements(pool_size);
    for (int ip = 0; ip < pool_size; ++ip)
    {
        receive_counts[ip] = this->numz[GlobalV::MY_POOL][ip] * this->ncxy;
        displacements[ip] = this->startz[GlobalV::MY_POOL][ip] * this->ncxy;
    }

    std::vector<double> local_data(this->nrxx);
    // The allgather below replicates one complete pool grid on every rank in
    // that pool. INT_BGROUP then sums all of those replicas, so divide each
    // local slab by the pool size to make each pool contribute exactly once.
    const double pool_normalization = 1.0 / static_cast<double>(pool_size);
    for (int ir = 0; ir < this->nrxx; ++ir)
    {
        local_data[ir] = data[ir] * pool_normalization;
    }

    std::vector<double> pool_data(this->ncxyz);
    // Collect the rank-local [xy][local_z] slabs into rank-contiguous blocks.
    MPI_Allgatherv(local_data.data(),
                   this->nrxx,
                   MPI_DOUBLE,
                   pool_data.data(),
                   receive_counts.data(),
                   displacements.data(),
                   MPI_DOUBLE,
                   POOL_WORLD);

    std::vector<double> global_layout(this->ncxyz);
    // Convert the rank-contiguous allgather result to the canonical
    // [xy][global_z] order required for element-wise reduction across pools.
    for (int ip = 0; ip < pool_size; ++ip)
    {
        const int local_nz = this->numz[GlobalV::MY_POOL][ip];
        const int global_z_start = this->startz[GlobalV::MY_POOL][ip];
        const int gathered_start = global_z_start * this->ncxy;
        for (int ixy = 0; ixy < this->ncxy; ++ixy)
        {
            for (int iz = 0; iz < local_nz; ++iz)
            {
                global_layout[ixy * this->ncz + global_z_start + iz] = pool_data[gathered_start + ixy * local_nz + iz];
            }
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, global_layout.data(), this->ncxyz, MPI_DOUBLE, MPI_SUM, INT_BGROUP);

    // Return only the z-slab owned by this rank under its pool's distribution.
    const int local_z_start = this->startz[GlobalV::MY_POOL][GlobalV::RANK_IN_POOL];
    for (int ixy = 0; ixy < this->ncxy; ++ixy)
    {
        for (int iz = 0; iz < this->nczp; ++iz)
        {
            data[ixy * this->nczp + iz] = global_layout[ixy * this->ncz + local_z_start + iz];
        }
    }
#else
    (void)data;
#endif
}

#ifdef __MPI
void Parallel_Grid::bcast(const double* const data_global, double* data_local, const int& rank, const bool is_sdft) const
{
    std::vector<double> zpiece(ncxy);
    for (int iz = 0; iz < this->ncz; ++iz)
    {
        ModuleBase::GlobalFunc::ZEROS(zpiece.data(), ncxy);
        if (rank == 0)
        {
            for (int ix = 0; ix < ncx; ix++)
            {
                for (int iy = 0; iy < ncy; iy++)
                {
                    const int ixy = ix * ncy + iy;
                    zpiece[ixy] = data_global[ixy * ncz + iz];
                }
            }
        }
        if (is_sdft)
        {
            this->zpiece_to_stogroup(zpiece.data(), iz, data_local);
        }
        else
        {
            this->zpiece_to_all(zpiece.data(), iz, data_local);
        }
    }
}

void Parallel_Grid::zpiece_to_all(double* zpiece, const int& iz, double* rho) const
{
    assert(!this->numz.empty());
    // ModuleBase::TITLE("Parallel_Grid","zpiece_to_all");
    MPI_Status ierror;

    const int znow = iz - this->startz[GlobalV::MY_POOL][GlobalV::RANK_IN_POOL];
    const int proc = this->whichpro[GlobalV::MY_POOL][iz];

    if (GlobalV::MY_POOL == 0)
    {
        // case 1: the first part of rho in processor 0.
        // and send zpeice to to other pools.
        if (proc == 0 && GlobalV::MY_RANK == 0)
        {
            for (int ir = 0; ir < ncxy; ir++)
            {
                rho[ir * nczp + znow] = zpiece[ir];
            }
            for (int ipool = 1; ipool < GlobalV::KPAR; ipool++)
            {
                MPI_Send(zpiece, ncxy, MPI_DOUBLE, this->whichpro[ipool][iz], iz, MPI_COMM_WORLD);
            }
        }

        // case 2: processor n (n!=0) receive rho from processor 0.
        // and the receive tag is iz.
        else if (proc == GlobalV::RANK_IN_POOL)
        {
            MPI_Recv(zpiece, ncxy, MPI_DOUBLE, 0, iz, MPI_COMM_WORLD, &ierror);
            for (int ir = 0; ir < ncxy; ir++)
            {
                rho[ir * nczp + znow] = zpiece[ir];
            }
        }

        // case 2: > first part rho: processor 0 send the rho
        // to all pools. The tag is iz, because processor may
        // send more than once, and the only tag to distinguish
        // them is iz.
        else if (GlobalV::RANK_IN_POOL == 0)
        {
            for (int ipool = 0; ipool < GlobalV::KPAR; ipool++)
            {
                MPI_Send(zpiece, ncxy, MPI_DOUBLE, this->whichpro[ipool][iz], iz, MPI_COMM_WORLD);
            }
        }
    } // GlobalV::MY_POOL == 0
    else
    {
        // GlobalV::ofs_running << "\n Receive charge density iz=" << iz << std::endl;
        //  the processors in other pools always receive rho from
        //  processor 0. the tag is 'iz'
        if (proc == GlobalV::MY_RANK)
        {
            MPI_Recv(zpiece, ncxy, MPI_DOUBLE, 0, iz, MPI_COMM_WORLD, &ierror);
            for (int ir = 0; ir < ncxy; ir++)
            {
                rho[ir * nczp + znow] = zpiece[ir];
            }
        }
    }

    // GlobalV::ofs_running << "\n iz = " << iz << " Done.";
    return;
}
#endif

#ifdef __MPI
void Parallel_Grid::zpiece_to_stogroup(double* zpiece, const int& iz, double* rho) const
{
    assert(!this->numz.empty());
    // TITLE("Parallel_Grid","zpiece_to_all");
    MPI_Status ierror;

    const int znow = iz - this->startz[GlobalV::MY_POOL][GlobalV::RANK_IN_POOL];
    const int proc = this->whichpro[GlobalV::MY_POOL][iz];

    if (GlobalV::MY_POOL == 0)
    {
        // case 1: the first part of rho in processor 0.
        // and send zpeice to to other pools.
        if (proc == 0 && GlobalV::RANK_IN_BPGROUP == 0)
        {
            for (int ir = 0; ir < ncxy; ir++)
            {
                rho[ir * nczp + znow] = zpiece[ir];
            }
            for (int ipool = 1; ipool < GlobalV::KPAR; ipool++)
            {
                MPI_Send(zpiece, ncxy, MPI_DOUBLE, this->whichpro[ipool][iz], iz, INT_BGROUP);
            }
        }

        // case 2: processor n (n!=0) receive rho from processor 0.
        // and the receive tag is iz.
        else if (proc == GlobalV::RANK_IN_POOL)
        {
            MPI_Recv(zpiece, ncxy, MPI_DOUBLE, 0, iz, INT_BGROUP, &ierror);
            for (int ir = 0; ir < ncxy; ir++)
            {
                rho[ir * nczp + znow] = zpiece[ir];
            }
        }

        // case 2: > first part rho: processor 0 send the rho
        // to all pools. The tag is iz, because processor may
        // send more than once, and the only tag to distinguish
        // them is iz.
        else if (GlobalV::RANK_IN_POOL == 0)
        {
            for (int ipool = 0; ipool < GlobalV::KPAR; ipool++)
            {
                MPI_Send(zpiece, ncxy, MPI_DOUBLE, this->whichpro[ipool][iz], iz, INT_BGROUP);
            }
        }
    } // MY_POOL == 0
    else
    {
        // ofs_running << "\n Receive charge density iz=" << iz << endl;
        //  the processors in other pools always receive rho from
        //  processor 0. the tag is 'iz'
        if (proc == GlobalV::RANK_IN_BPGROUP)
        {
            MPI_Recv(zpiece, ncxy, MPI_DOUBLE, 0, iz, INT_BGROUP, &ierror);
            for (int ir = 0; ir < ncxy; ir++)
            {
                rho[ir * nczp + znow] = zpiece[ir];
            }
        }
    }

    // ofs_running << "\n iz = " << iz << " Done.";
    return;
}
void Parallel_Grid::reduce(double* rhotot, const double* const rhoin, const bool reduce_all_pool) const
{
    // ModuleBase::TITLE("Parallel_Grid","reduce");

    // if not the first pool, wait here until processpr 0
    // send the Barrier command.
    if (!reduce_all_pool && GlobalV::MY_POOL != 0)
    {
        return;
    }

    double* zpiece = new double[this->ncxy];

    for (int iz = 0; iz < this->ncz; iz++)
    {
        const int znow = iz - this->startz[GlobalV::MY_POOL][GlobalV::RANK_IN_POOL];
        const int proc = this->whichpro[GlobalV::MY_POOL][iz];
        const int proc_loc = this->whichpro_loc[GlobalV::MY_POOL][iz]; // Obtain the local processor index in the pool
        ModuleBase::GlobalFunc::ZEROS(zpiece, this->ncxy);
        int tag = iz;
        MPI_Status ierror;

        // Local processor 0 collects data from all other processors in the pool
        // proc = proc_loc if GlobalV::MY_POOL == 0
        if (proc_loc == GlobalV::RANK_IN_POOL)
        {
            for (int ir = 0; ir < ncxy; ir++)
            {
                zpiece[ir] = rhoin[ir * this->nczp + znow];
            }
            // Send data to the root of the pool
            if (GlobalV::RANK_IN_POOL != 0)
            {
                MPI_Send(zpiece, ncxy, MPI_DOUBLE, 0, tag, POOL_WORLD);
            }
        }

        // The root of the pool receives data from other processors
        if (GlobalV::RANK_IN_POOL == 0 && proc_loc != GlobalV::RANK_IN_POOL)
        {
            MPI_Recv(zpiece, ncxy, MPI_DOUBLE, proc_loc, tag, POOL_WORLD, &ierror);
        }

        if (GlobalV::RANK_IN_POOL == 0)
        {
            for (int ixy = 0; ixy < this->ncxy; ++ixy)
            {
                rhotot[ixy * ncz + iz] = zpiece[ixy];
            }
        }
    }

    delete[] zpiece;

    return;
}
#endif
