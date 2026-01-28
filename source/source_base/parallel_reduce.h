#ifndef PARALLEL_REDUCE_H
#define PARALLEL_REDUCE_H

#ifdef __MPI
#include <mpi.h>
#endif

#include <cassert>
#include <complex>

// using std::complex;

namespace Parallel_Reduce
{
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

void reduce_int_diag(int& object); // mohan add 2012-01-12

void reduce_int_grid(int* object, const int n); // mohan add 2012-01-12

// reduce double only in this pool
// (each pool contain different k points)
void reduce_double_grid(double* object, const int n);
void reduce_double_diag(double* object, const int n);

void reduce_double_allpool(const int& npool, const int& nproc_in_pool, double& object);
void reduce_double_allpool(const int& npool, const int& nproc_in_pool, double* object, const int n);

void gather_int_all(int& v, int* all);

bool check_if_equal(double& v); // mohan add 2009-11-11

template <class T, class TI>
inline void ZEROS(std::complex<T>* u, const TI n)
{
    assert(n >= 0);
    for (TI i = 0; i < n; i++)
    {
        u[i] = std::complex<T>(0.0, 0.0);
    }
    return;
}

template <class T, class TI>
inline void ZEROS(T* u, const TI n)
{
    assert(n >= 0);
    for (TI i = 0; i < n; i++)
    {
        u[i] = 0;
    }
}
} // namespace Parallel_Reduce

#endif
