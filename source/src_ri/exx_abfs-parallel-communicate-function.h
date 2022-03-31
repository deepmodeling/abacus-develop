#ifdef __MPI
#ifndef EXX_ABFS_PARALLEL_COMMUNICATE_FUNCTION_H
#define EXX_ABFS_PARALLEL_COMMUNICATE_FUNCTION_H

#include "exx_abfs-parallel.h"
#include "module_orbital/parallel_orbitals.h"
#include <vector>
#include <utility>
#ifdef __MPI
#include "mpi.h"
#endif
using namespace std;

class Exx_Abfs::Parallel::Communicate::Function
{
public:
#ifdef __MPI
	static std::vector<std::pair<std::vector<bool>,std::vector<bool>>> get_atom_in_2D_list(const MPI_Comm &mpi_comm, const Parallel_Orbitals &pv);
#endif
};

#endif
#endif