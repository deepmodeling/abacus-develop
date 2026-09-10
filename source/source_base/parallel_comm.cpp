#if defined __MPI

#include "mpi.h"
#include "parallel_global.h"

// Two-level pool terminology used across the parallel layer:
//   - k-pool: a group of processes that share one subset of k-points. The
//     processes are split into KPAR k-pools first (divide_pools); this split
//     is independent of bndpar. MY_POOL / KP_WORLD refer to this level.
//   - band-pool: a sub-group of one k-pool, created afterwards by dividing
//     the k-pool into BNDPAR band groups. NPROC_IN_POOL / RANK_IN_POOL /
//     POOL_WORLD refer to this level, i.e. the term "pool" in those globals
//     means the (k-pool, band-group) cell, NOT the k-pool itself.
MPI_Comm POOL_WORLD;   // one band-pool (k-pool x band-group cell): plane waves are distributed, k-points and the band window are shared.
MPI_Comm KP_WORLD;     // links k-pools: only k-points differ; same rank_in_pool position in every k-pool. Valid ONLY when k-pools are equally sized (NPROC % KPAR == 0), otherwise MPI_COMM_NULL.
MPI_Comm BP_WORLD;     // links band groups inside one k-pool: only the band window differs; k-points and plane-wave slab are the same. One communicator per rank position.
MPI_Comm INT_BGROUP;   // same band-group index across all k-pools (plus the plane-wave ranks of that band group): k-points differ, the band window is the same. Always valid, also for uneven k-pools.
MPI_Comm GRID_WORLD; // mohan add 2012-01-13
MPI_Comm DIAG_WORLD; // mohan add 2012-01-13

MPICommGroup::MPICommGroup(MPI_Comm parent_comm)
    : parent_comm(parent_comm)
{
    MPI_Comm_size(parent_comm, &this->gsize);
    MPI_Comm_rank(parent_comm, &this->grank);
}

MPICommGroup::~MPICommGroup()
{
    if (group_comm != MPI_COMM_NULL)
    {
        MPI_Comm_free(&group_comm);
    }
    if (inter_comm != MPI_COMM_NULL)
    {
        MPI_Comm_free(&inter_comm);
    }
}

void MPICommGroup::divide_group_comm(const int& ngroup, const bool assert_even)
{
    this->ngroups = ngroup;
    Parallel_Global::divide_mpi_groups(this->gsize,
                                       ngroup,
                                       this->grank,
                                       this->nprocs_in_group,
                                       this->my_group,
                                       this->rank_in_group,
                                       assert_even);

    MPI_Comm_split(parent_comm, my_group, rank_in_group, &group_comm);
    if(this->gsize % ngroup == 0)
    {
        this->is_even = true;
    }

    if (this->is_even)
    {
        MPI_Comm_split(parent_comm, my_inter, rank_in_inter, &inter_comm);
    }
}

#endif