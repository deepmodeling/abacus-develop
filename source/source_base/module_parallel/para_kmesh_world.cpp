#include "para_kmesh_world.h"

#include <algorithm>
#include <cassert>

namespace Parallel
{

ParaKmeshWorld::ParaKmeshWorld(int nkstot, int nspin)
    : ParaWorld("kmesh"), kpar_(1), my_pool_(0), rank_in_pool_(0),
      nproc_(1), nspin_(nspin), nkstot_(nkstot)
{
    distribute_kpoints();
    nks_local_ = nkstot_;
    startk_global_ = 0;
}

#ifdef __MPI
ParaKmeshWorld::ParaKmeshWorld(const MPI_Comm& comm, int kpar, int my_pool, int nproc, int nkstot, int nspin)
    : ParaWorld("kmesh", comm), kpar_(kpar), my_pool_(my_pool),
      rank_in_pool_(rank()), nproc_(nproc), nspin_(nspin), nkstot_(nkstot)
{
    distribute_kpoints();
    nks_local_ = nks_pool_[my_pool_];
    startk_global_ = startk_pool_[my_pool_];
}
#endif

void ParaKmeshWorld::distribute_kpoints()
{
    // k-points per pool (evenly divided, remainder to front)
    nks_pool_.resize(kpar_, 0);
    const int nks_ave = nkstot_ / kpar_;
    const int nks_rem = nkstot_ % kpar_;
    for (int i = 0; i < kpar_; ++i)
    {
        nks_pool_[i] = nks_ave + (i < nks_rem ? 1 : 0);
    }

    // global start index per pool
    startk_pool_.resize(kpar_, 0);
    for (int i = 1; i < kpar_; ++i)
    {
        startk_pool_[i] = startk_pool_[i - 1] + nks_pool_[i - 1];
    }

    // pool index per k-point
    whichpool_.resize(nkstot_, 0);
    for (int p = 0; p < kpar_; ++p)
    {
        for (int ik = 0; ik < nks_pool_[p]; ++ik)
        {
            whichpool_[startk_pool_[p] + ik] = p;
        }
    }

    // first world rank per pool
    startpro_pool_.resize(kpar_, 0);
    const int nproc_ave = nproc_ / kpar_;
    const int nproc_rem = nproc_ % kpar_;
    for (int i = 1; i < kpar_; ++i)
    {
        startpro_pool_[i] = startpro_pool_[i - 1] + nproc_ave + (i - 1 < nproc_rem ? 1 : 0);
    }
}

int ParaKmeshWorld::nks_pool(int pool) const
{
    assert(pool >= 0 && pool < kpar_);
    return nks_pool_[pool];
}

int ParaKmeshWorld::startk_pool(int pool) const
{
    assert(pool >= 0 && pool < kpar_);
    return startk_pool_[pool];
}

int ParaKmeshWorld::which_pool(int ik_global) const
{
    assert(ik_global >= 0 && ik_global < nkstot_);
    return whichpool_[ik_global];
}

int ParaKmeshWorld::startpro_pool(int pool) const
{
    assert(pool >= 0 && pool < kpar_);
    return startpro_pool_[pool];
}

int ParaKmeshWorld::max_nks_pool() const
{
    return *std::max_element(nks_pool_.begin(), nks_pool_.end());
}

void ParaKmeshWorld::pool_collection(double& value, const double* wk, int ik) const
{
#ifdef __MPI
    const int ik_local = ik - startk_pool_[my_pool_];
    const int pool = whichpool_[ik];

    if (rank_in_pool_ == 0)
    {
        if (my_pool_ == 0)
        {
            if (pool == 0)
            {
                value = wk[ik_local];
            }
            else
            {
                MPI_Status status;
                MPI_Recv(&value, 1, MPI_DOUBLE, startpro_pool_[pool], ik, MPI_COMM_WORLD, &status);
            }
        }
        else
        {
            if (my_pool_ == pool)
            {
                MPI_Send(&wk[ik_local], 1, MPI_DOUBLE, 0, ik, MPI_COMM_WORLD);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
#else
    value = wk[ik];
#endif
}

template <typename T>
void ParaKmeshWorld::pool_collection(T* value, const T* w, int dim, int ik) const
{
#ifdef __MPI
    const int ik_local = ik - startk_pool_[my_pool_];
    const int begin = ik_local * dim;
    const T* src = &w[begin];

    // nspin==2 restricts to pool 0 (legacy behavior from Parallel_Kpoints)
    const int pool = (nspin_ == 2) ? 0 : whichpool_[ik];

    if (rank_in_pool_ == 0)
    {
        if (my_pool_ == 0)
        {
            if (pool == 0)
            {
                std::copy(src, src + dim, value);
            }
            else
            {
                MPI_Status status;
                MPI_Recv(value, dim * sizeof(T), MPI_BYTE, startpro_pool_[pool], ik * 2, MPI_COMM_WORLD, &status);
            }
        }
        else
        {
            if (my_pool_ == pool)
            {
                MPI_Send(src, dim * sizeof(T), MPI_BYTE, 0, ik * 2, MPI_COMM_WORLD);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
#else
    const int begin = ik * dim;
    std::copy(&w[begin], &w[begin] + dim, value);
#endif
}

void ParaKmeshWorld::gather_kvec(const std::vector<double>& vec_local, std::vector<double>& vec_global) const
{
#ifdef __MPI
    int world_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    const bool is_pool_root = (world_rank == startpro_pool_[my_pool_]);

    vec_global.resize(nkstot_ * 3, 0.0);
    if (is_pool_root)
    {
        for (int i = 0; i < nks_local_; ++i)
        {
            const int gk = startk_global_ + i;
            vec_global[gk * 3 + 0] = vec_local[i * 3 + 0];
            vec_global[gk * 3 + 1] = vec_local[i * 3 + 1];
            vec_global[gk * 3 + 2] = vec_local[i * 3 + 2];
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, vec_global.data(), nkstot_ * 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    vec_global = vec_local;
#endif
}

// explicit instantiation
template void ParaKmeshWorld::pool_collection<double>(double*, const double*, int, int) const;
template void ParaKmeshWorld::pool_collection<std::complex<double>>(std::complex<double>*, const std::complex<double>*, int, int) const;

} // namespace Parallel
