#include "parallel_reduce.h"

#include "parallel_comm.h"

#include <vector>

template <>
void Parallel_Reduce::reduce_all<int>(int& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<long long>(long long& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

void Parallel_Reduce::reduce_int_diag(int& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_INT, MPI_SUM, DIAG_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<double>(double& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<float>(float& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<int>(int* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<long long>(long long* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

void Parallel_Reduce::reduce_int_grid(int* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_SUM, GRID_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<double>(double* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

void Parallel_Reduce::reduce_double_grid(double* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, GRID_WORLD);
#endif
    return;
}

void Parallel_Reduce::reduce_double_diag(double* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, DIAG_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<float>(float& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_FLOAT, MPI_SUM, POOL_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<double>(double& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_SUM, POOL_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<int>(int* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_SUM, POOL_WORLD);
#endif
}

template <>
void Parallel_Reduce::reduce_pool<double>(double* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, POOL_WORLD);
#endif
    return;
}

// (1) the value is same in each pool.
// (2) we need to reduce the value from different pool.
void Parallel_Reduce::reduce_double_allpool(const int& npool, const int& nproc_in_pool, double& object)
{
    if (npool == 1) 
    {
        return;
    }
#ifdef __MPI
    double swap = object / nproc_in_pool;
    MPI_Allreduce(&swap, &object, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}

// (1) the value is same in each pool.
// (2) we need to reduce the value from different pool.
void Parallel_Reduce::reduce_double_allpool(const int& npool, const int& nproc_in_pool, double* object, const int n)
{
    if (npool == 1) 
    {
        return;
    }
#ifdef __MPI
    std::vector<double> swap(n, 0.0);
    for (int i = 0; i < n; i++)
    {
        swap[i] = object[i] / nproc_in_pool;
    }
    MPI_Allreduce(swap.data(), object, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}

template <>
void Parallel_Reduce::reduce_all<std::complex<double>>(std::complex<double>& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

// LiuXh add 2019-07-16
template <>
void Parallel_Reduce::reduce_all<std::complex<double>>(std::complex<double>* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}


template <>
void Parallel_Reduce::reduce_all<std::complex<float>>(std::complex<float>& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_C_FLOAT_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

// LiuXh add 2019-07-16
template <>
void Parallel_Reduce::reduce_all<std::complex<float>>(std::complex<float>* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_C_FLOAT_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<std::complex<double>>(std::complex<double>& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<std::complex<float>>(std::complex<float>* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_C_FLOAT_COMPLEX, MPI_SUM, POOL_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<std::complex<double>>(std::complex<double>* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
#endif
    return;
}

void Parallel_Reduce::gather_int_all(int& v, int* all)
{
#ifdef __MPI
    assert(all != nullptr);
    MPI_Allgather(&v, 1, MPI_INT, all, 1, MPI_INT, MPI_COMM_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_min<int>(int& v)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif
}

template <>
void Parallel_Reduce::reduce_min<float>(float& v)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#endif
}

template <>
void Parallel_Reduce::reduce_min<double>(double& v)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
}

template <>
void Parallel_Reduce::reduce_max<float>(float& v)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
#endif
}

template <>
void Parallel_Reduce::reduce_max<double>(double& v)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}

#ifdef __MPI
template <>
void Parallel_Reduce::reduce_sum<int>(int& object, MPI_Comm comm)
{
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_INT, MPI_SUM, comm);
}

template <>
void Parallel_Reduce::reduce_sum<float>(float& object, MPI_Comm comm)
{
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_FLOAT, MPI_SUM, comm);
}

template <>
void Parallel_Reduce::reduce_sum<double>(double& object, MPI_Comm comm)
{
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_SUM, comm);
}

template <>
void Parallel_Reduce::reduce_sum<int>(int* object, const int n, MPI_Comm comm)
{
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_SUM, comm);
}

template <>
void Parallel_Reduce::reduce_sum<float>(float* object, const int n, MPI_Comm comm)
{
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_FLOAT, MPI_SUM, comm);
}

template <>
void Parallel_Reduce::reduce_sum<double>(double* object, const int n, MPI_Comm comm)
{
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, comm);
}

template <>
void Parallel_Reduce::reduce_max<float>(float& v, MPI_Comm comm)
{
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_FLOAT, MPI_MAX, comm);
}

template <>
void Parallel_Reduce::reduce_max<double>(double& v, MPI_Comm comm)
{
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, comm);
}

void Parallel_Reduce::reduce_bool_or(bool& v, MPI_Comm comm)
{
    int reduced = v ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &reduced, 1, MPI_INT, MPI_LOR, comm);
    v = (reduced != 0);
}

void Parallel_Reduce::reduce_bool_and(bool& v, MPI_Comm comm)
{
    int reduced = v ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &reduced, 1, MPI_INT, MPI_LAND, comm);
    v = (reduced != 0);
}
#endif

template <>
void Parallel_Reduce::reduce_max_pool<double>(const int& nproc_in_pool, double& v)
{
#ifdef __MPI
    if (nproc_in_pool == 1) 
    {
        return;
    }
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, POOL_WORLD);
#endif
}
template <>
void Parallel_Reduce::reduce_min_pool<double>(const int& nproc_in_pool, double& v)
{
#ifdef __MPI
    if (nproc_in_pool == 1) 
    {
        return;
    }
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MIN, POOL_WORLD);
#endif
}
