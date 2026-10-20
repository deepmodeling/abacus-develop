#ifndef PARALLEL_REDUCE_H
#define PARALLEL_REDUCE_H

#ifdef __MPI
#include <mpi.h>
#endif

#include <cassert>
#include <complex>

namespace Parallel_Reduce
{

#ifdef __MPI
template <typename T>
struct MPI_Type;

template <>
struct MPI_Type<int>
{
    static const MPI_Datatype value;
};

template <>
struct MPI_Type<double>
{
    static const MPI_Datatype value;
};

template <>
struct MPI_Type<float>
{
    static const MPI_Datatype value;
};

template <>
struct MPI_Type<std::complex<double>>
{
    static const MPI_Datatype value;
};

template <>
struct MPI_Type<std::complex<float>>
{
    static const MPI_Datatype value;
};

template <>
struct MPI_Type<long long>
{
    static const MPI_Datatype value;
};

template <>
struct MPI_Type<unsigned long long>
{
    static const MPI_Datatype value;
};
#endif

/// reduce in all process
template <typename T>
void reduce_all(T& object);
template <typename T>
void reduce_all(T* object, const int n);
template <typename T>
void reduce_pool(T& object);
template <typename T>
void reduce_pool(T* object, const int n);
template <typename T>
void reduce_min(T& v);
template <typename T>
void reduce_max(T& v);
/** @brief Reduce an array of maxima across all ranks. */
void reduce_max(double* values, const int count);
template <typename T>
void reduce_min_pool(const int& nproc_in_pool, T& v);
template <typename T>
void reduce_max_pool(const int& nproc_in_pool, T& v);

// reduce double only in this pool
// (each pool contain different k points)
void reduce_double_grid(double* object, const int n);
void reduce_double_diag(double* object, const int n);

void reduce_double_allpool(const int& npool, const int& nproc_in_pool, double& object);
void reduce_double_allpool(const int& npool, const int& nproc_in_pool, double* object, const int n);

void gather_int_all(int& v, int* all);

#ifdef __MPI
/// Communicator-aware counterparts of the MPI_COMM_WORLD-bound helpers above.
/// The RI and RPA stages run on caller-provided sub-communicators, which those
/// helpers cannot express; route every reduction on such a communicator
/// through these overloads so the raw MPI calls stay centralized here.
template <typename T>
void reduce_all(const MPI_Comm mpi_comm, T& object);
template <typename T>
void reduce_all(const MPI_Comm mpi_comm, T* object, const int n);
template <typename T>
void reduce_max(const MPI_Comm mpi_comm, T& object);
template <typename T>
void reduce_min(const MPI_Comm mpi_comm, T& object);
/// @brief Largest value together with the rank that supplied it (MPI_MAXLOC).
void reduce_max_loc(const MPI_Comm mpi_comm, double& value, int& rank);
/// @brief Allgather one int per rank on an arbitrary communicator.
void gather_int_all(const MPI_Comm mpi_comm, int& v, int* all);
#endif

} // namespace Parallel_Reduce

#endif
