#ifndef PARALLEL_COMMON_H
#define PARALLEL_COMMON_H

#ifdef __MPI
#include "mpi.h"
#endif

#include <complex>
#include <string>

namespace Parallel_Common
{

//(1) bcast array
void bcast_complex_double(std::complex<double>* object, const int n);
void bcast_string(std::string* object, const int n);
void bcast_double(double* object, const int n);
void bcast_int(int* object, const int n);
void bcast_char(char* object, const int n);

//(2) bcast single
void bcast_complex_double(std::complex<double>& object);
void bcast_string(std::string& object);
void bcast_double(double& object);
void bcast_int(int& object);
void bcast_bool(bool& object);

#ifdef __MPI
//(3) communicator- and root-aware counterparts for the sub-communicator
//    stages (RI/RPA) that run outside MPI_COMM_WORLD and outside root 0
void bcast_int(int* object, const int n, const MPI_Comm mpi_comm, const int root);
void bcast_double(double* object, const int n, const MPI_Comm mpi_comm, const int root);
/// @brief Rank of the caller on an arbitrary communicator.
int my_rank(const MPI_Comm mpi_comm);
/// @brief Number of ranks on an arbitrary communicator.
int nproc(const MPI_Comm mpi_comm);
/// @brief Barrier on an arbitrary communicator.
void barrier(const MPI_Comm mpi_comm);
#endif

} // namespace Parallel_Common

#endif
