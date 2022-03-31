#ifndef EXX_ABFS_UTIL_H
#define EXX_ABFS_UTIL_H

#include "exx_abfs.h"
#ifdef __MPI
#include "mpi.h"
#endif

#include <vector>
#include <string>

#ifdef __MPI
class Exx_Abfs::Util
{
public:
	static void bcast( std::vector<std::string> &v, const int rank_src, const MPI_Comm &mpi_comm );
};
#endif

#endif