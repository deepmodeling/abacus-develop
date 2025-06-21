#if defined __MPI

#include "mpi.h"
#include "parallel_global.h"

MPI_Comm POOL_WORLD; //groups for different plane waves. In this group, only plane waves are different. K-points and bands are the same.
MPI_Comm KP_WORLD;   // groups for differnt k. In this group, only k-points are different. Bands and plane waves are the same.
MPI_Comm BP_WORLD;   // groups for differnt bands. In this group, only bands are different. K-points and plane waves are the same.
MPI_Comm INT_BGROUP; // internal comm groups for same bands. In this group, only bands are the same. K-points and plane waves are different.
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