#include "para_rgrid_world.h"

#include <cassert>

namespace Parallel
{

ParaRgridWorld::ParaRgridWorld(int ncx, int ncy, int ncz)
    : ParaWorld("rgrid"), ncx_(ncx), ncy_(ncy), ncz_(ncz)
{
    assert(ncx > 0 && ncy > 0 && ncz > 0);
    distribute_z();
    nczp_ = numz_[0];
}

#ifdef __MPI
ParaRgridWorld::ParaRgridWorld(const MPI_Comm& comm, int ncx, int ncy, int ncz)
    : ParaWorld("rgrid", comm), ncx_(ncx), ncy_(ncy), ncz_(ncz)
{
    assert(ncx > 0 && ncy > 0 && ncz > 0);
    distribute_z();
    nczp_ = numz_[rank()];
}
#endif

void ParaRgridWorld::distribute_z()
{
    const int np = size();
    numz_.resize(np, 0);
    startz_.resize(np, 0);
    whichpro_.resize(ncz_, 0);

    // Evenly distribute z-planes, remainder to front processes
    const int base = ncz_ / np;
    const int rem = ncz_ % np;
    int acc = 0;
    for (int p = 0; p < np; ++p)
    {
        numz_[p] = base + (p < rem ? 1 : 0);
        startz_[p] = acc;
        acc += numz_[p];
    }

    // Build owner table
    for (int p = 0; p < np; ++p)
    {
        for (int iz = 0; iz < numz_[p]; ++iz)
        {
            whichpro_[startz_[p] + iz] = p;
        }
    }
}

int ParaRgridWorld::numz(int p) const
{
    assert(p >= 0 && p < static_cast<int>(numz_.size()));
    return numz_[p];
}

int ParaRgridWorld::startz(int p) const
{
    assert(p >= 0 && p < static_cast<int>(startz_.size()));
    return startz_[p];
}

int ParaRgridWorld::whichpro(int iz) const
{
    assert(iz >= 0 && iz < ncz_);
    return whichpro_[iz];
}

// ===== Cross-domain operations =====

void ParaRgridWorld::reduce_across_pools(double* data, const ParaWorld& kmesh_world) const
{
#ifdef __MPI
    if (!kmesh_world.valid()) return;
    if (kmesh_world.size() <= 1) return;

    assert(data != nullptr);

    // Equal-sized pools: corresponding ranks have identical z-slab layouts,
    // so local buffers can be summed directly without redistribution.
    MPI_Allreduce(MPI_IN_PLACE, data, nrxx(), MPI_DOUBLE, MPI_SUM, kmesh_world.comm());
#else
    (void)data;
    (void)kmesh_world;
#endif
}

void ParaRgridWorld::bcast_data(const double* data_global, double* data_local,
                                 const ParaWorld& comm_world, int root) const
{
    // Serial or single-process: just copy local slab
    if (!comm_world.valid() || comm_world.size() == 1)
    {
        const int ncxy = ncx_ * ncy_;
        const int z_start = startz_[rank()];
        for (int ixy = 0; ixy < ncxy; ++ixy)
        {
            for (int iz = 0; iz < nczp_; ++iz)
            {
                data_local[ixy * nczp_ + iz] = data_global[ixy * ncz_ + z_start + iz];
            }
        }
        return;
    }

#ifdef __MPI
    // Broadcast z-plane by z-plane
    std::vector<double> zpiece(ncx_ * ncy_);
    for (int iz = 0; iz < ncz_; ++iz)
    {
        if (comm_world.rank() == root)
        {
            for (int ix = 0; ix < ncx_; ++ix)
            {
                for (int iy = 0; iy < ncy_; ++iy)
                {
                    zpiece[ix * ncy_ + iy] = data_global[(ix * ncy_ + iy) * ncz_ + iz];
                }
            }
        }
        MPI_Bcast(zpiece.data(), ncx_ * ncy_, MPI_DOUBLE, root, comm_world.comm());

        // Store z-plane if this process owns it
        const int znow = iz - startz_[comm_world.rank()];
        if (znow >= 0 && znow < nczp_)
        {
            for (int ixy = 0; ixy < ncx_ * ncy_; ++ixy)
            {
                data_local[ixy * nczp_ + znow] = zpiece[ixy];
            }
        }
    }
#else
    (void)data_global;
    (void)data_local;
    (void)root;
#endif
}

void ParaRgridWorld::reduce_data(double* rhotot, const double* rhoin,
                                  const ParaWorld& comm_world) const
{
    // Serial: just copy local slab to global grid
    if (!comm_world.valid() || comm_world.size() == 1)
    {
        const int ncxy = ncx_ * ncy_;
        const int z_start = startz_[comm_world.rank()];
        for (int ixy = 0; ixy < ncxy; ++ixy)
        {
            for (int iz = 0; iz < nczp_; ++iz)
            {
                rhotot[ixy * ncz_ + z_start + iz] = rhoin[ixy * nczp_ + iz];
            }
        }
        return;
    }

#ifdef __MPI
    // Gather local z-slabs from all processes
    const int np = comm_world.size();
    std::vector<int> local_z_counts(np);
    std::vector<int> receive_counts(np);
    std::vector<int> displacements(np, 0);

    int my_nczp = nczp_;
    MPI_Allgather(&my_nczp, 1, MPI_INT, local_z_counts.data(), 1, MPI_INT, comm_world.comm());

    int total_z = 0;
    for (int p = 0; p < np; ++p)
    {
        receive_counts[p] = local_z_counts[p] * ncx_ * ncy_;
        if (p > 0)
        {
            displacements[p] = displacements[p - 1] + receive_counts[p - 1];
        }
        total_z += local_z_counts[p];
    }
    assert(total_z == ncz_);

    std::vector<double> gathered;
    if (comm_world.rank() == 0)
    {
        gathered.resize(ncx_ * ncy_ * ncz_);
    }

    MPI_Gatherv(rhoin, nrxx(), MPI_DOUBLE,
                gathered.data(), receive_counts.data(), displacements.data(),
                MPI_DOUBLE, 0, comm_world.comm());

    if (comm_world.rank() == 0)
    {
        // Convert from rank-contiguous [xy][local_z] to canonical [xy][global_z]
        int global_z_start = 0;
        for (int p = 0; p < np; ++p)
        {
            const int local_nz = local_z_counts[p];
            for (int ixy = 0; ixy < ncx_ * ncy_; ++ixy)
            {
                for (int iz = 0; iz < local_nz; ++iz)
                {
                    rhotot[ixy * ncz_ + global_z_start + iz]
                        = gathered[displacements[p] + ixy * local_nz + iz];
                }
            }
            global_z_start += local_nz;
        }
    }
#else
    (void)rhotot;
    (void)rhoin;
#endif
}

} // namespace Parallel
