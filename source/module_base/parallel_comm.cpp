#if defined __MPI

#include "mpi.h"
#include "parallel_global.h"

MPI_Comm POOL_WORLD;
MPI_Comm INTER_POOL; // communicator among different pools
MPI_Comm STO_WORLD;
MPI_Comm PARAPW_WORLD;
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