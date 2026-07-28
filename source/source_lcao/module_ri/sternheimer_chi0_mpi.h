#ifndef STERNHEIMER_CHI0_MPI_H
#define STERNHEIMER_CHI0_MPI_H

#include <complex>
#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

namespace module_ri
{
namespace sternheimer_chi0
{

#ifdef __MPI
void reduce_branch_to_root(std::vector<std::complex<double>>& branch, int root, MPI_Comm communicator);
#endif

} // namespace sternheimer_chi0
} // namespace module_ri

#endif
