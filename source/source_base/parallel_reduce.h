#ifndef PARALLEL_REDUCE_H
#define PARALLEL_REDUCE_H

#ifdef __MPI
#include <mpi.h>
#else
// The communicator-aware overloads below name MPI_Comm. A serial build has no
// such type, so reuse the stand-in that Parallel_Common declares rather than
// defining a second one that a caller could not interchange.
#include "parallel_common.h"
#endif

#include <cassert>
#include <complex>

namespace Parallel_Reduce
{
#ifndef __MPI
using MPI_Comm = Parallel_Common::MPI_Comm;
#endif

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

/// Communicator-aware counterparts of the MPI_COMM_WORLD-bound helpers above.
/// The RI and RPA stages run on caller-provided sub-communicators, which those
/// helpers cannot express; route every reduction on such a communicator
/// through these overloads so the raw MPI calls stay centralized here.
///
/// These declarations are visible in serial builds too: the implementations
/// guard their bodies instead of the whole function, so a caller does not have
/// to reproduce the #ifdef at every use site.
/// The scalar overloads are plain functions rather than templates: with
/// OpenMPI, `MPI_Comm` is a pointer, so a template `reduce_all(const MPI_Comm,
/// T&)` is ambiguous against the pre-existing `reduce_all(T*, const int)` for
/// an `int` payload. The set of supported scalar types is closed and matches
/// the instantiations in the translation unit, so naming them explicitly costs
/// nothing and keeps every call unambiguous.
void reduce_all(const MPI_Comm mpi_comm, int& object);
void reduce_all(const MPI_Comm mpi_comm, double& object);
void reduce_all(const MPI_Comm mpi_comm, unsigned long long& object);
template <typename T>
void reduce_all(const MPI_Comm mpi_comm, T* object, const int n);
/// Same reasoning as reduce_all: the supported scalar set is closed.
void reduce_max(const MPI_Comm mpi_comm, int& object);
void reduce_max(const MPI_Comm mpi_comm, double& object);
void reduce_max(const MPI_Comm mpi_comm, unsigned long long& object);
void reduce_min(const MPI_Comm mpi_comm, int& object);
void reduce_min(const MPI_Comm mpi_comm, double& object);
void reduce_min(const MPI_Comm mpi_comm, unsigned long long& object);
/// @brief Largest value together with the rank that supplied it (MPI_MAXLOC).
/// "loc" is MPI's contraction of "location"; the name is spelled out here so a
/// reader who does not know MPI_MAXLOC can still tell what the call does.
void reduce_max_with_rank(const MPI_Comm mpi_comm, double& value, int& rank);
/// @brief Allgather one int per rank on an arbitrary communicator.
void gather_int_all(const MPI_Comm mpi_comm, int& v, int* all);
/// @brief Allgatherv with per-rank counts and displacements on an arbitrary communicator.
void gatherv_int_all(const MPI_Comm mpi_comm,
                     const int* send_buffer,
                     const int send_count,
                     int* receive_buffer,
                     const int* receive_counts,
                     const int* receive_displacements);

} // namespace Parallel_Reduce

#endif
