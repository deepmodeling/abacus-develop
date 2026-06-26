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

bool check_if_equal(double& v); // mohan add 2009-11-11

} // namespace Parallel_Reduce

#endif
