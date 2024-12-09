#ifndef PARALLEL_COMM_H
#define PARALLEL_COMM_H

#ifdef __MPI
#include "mpi.h"
extern MPI_Comm POOL_WORLD;
extern MPI_Comm INTER_POOL; // communicator among different pools
extern MPI_Comm STO_WORLD;
extern MPI_Comm PARAPW_WORLD;
extern MPI_Comm GRID_WORLD; // mohan add 2012-01-13
extern MPI_Comm DIAG_WORLD; // mohan add 2012-01-13


class MPICommGroup
{
public:
    MPICommGroup(MPI_Comm parent_comm);
    ~MPICommGroup();
    void divide_group_comm(const int& ngroup, const bool assert_even = true);
public:
  bool is_even = false; ///< whether the group is even

  MPI_Comm parent_comm = MPI_COMM_NULL; ///< parent communicator
  int gsize = 0;                        ///< size of parent communicator
  int grank = 0;                        ///< rank of parent communicator

  MPI_Comm group_comm = MPI_COMM_NULL; ///< group communicator
  int ngroups = 0;                     ///< number of groups
  int nprocs_in_group = 0;             ///<  number of processes in the group
  int my_group = 0;                    ///< the group index
  int rank_in_group = 0;               ///< the rank in the group

  MPI_Comm inter_comm = MPI_COMM_NULL; ///< inter communicator
  bool has_inter_comm = false;         ///< whether has inter communicator
  int& nprocs_in_inter = ngroups;      ///< number of processes in the inter communicator
  int& my_inter = rank_in_group;  ///< the rank in the inter communicator
  int& rank_in_inter = my_group;  ///< the inter group index
};

#endif

#endif // PARALLEL_COMM_H