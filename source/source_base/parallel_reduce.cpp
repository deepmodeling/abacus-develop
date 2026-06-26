#include "parallel_reduce.h"

#include "parallel_comm.h"

#include <vector>

#ifdef __MPI
const MPI_Datatype Parallel_Reduce::MPI_Type<int>::value = MPI_INT;
const MPI_Datatype Parallel_Reduce::MPI_Type<double>::value = MPI_DOUBLE;
const MPI_Datatype Parallel_Reduce::MPI_Type<float>::value = MPI_FLOAT;
const MPI_Datatype Parallel_Reduce::MPI_Type<std::complex<double>>::value = MPI_DOUBLE_COMPLEX;
const MPI_Datatype Parallel_Reduce::MPI_Type<std::complex<float>>::value = MPI_C_FLOAT_COMPLEX;
const MPI_Datatype Parallel_Reduce::MPI_Type<long long>::value = MPI_LONG_LONG;
#endif

template <typename T>
void Parallel_Reduce::reduce_all(T& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_Type<T>::value, MPI_SUM, MPI_COMM_WORLD);
#endif
}

template <typename T>
void Parallel_Reduce::reduce_all(T* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_Type<T>::value, MPI_SUM, MPI_COMM_WORLD);
#endif
}

template <typename T>
void Parallel_Reduce::reduce_pool(T& object)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_Type<T>::value, MPI_SUM, POOL_WORLD);
#endif
}

template <typename T>
void Parallel_Reduce::reduce_pool(T* object, const int n)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_Type<T>::value, MPI_SUM, POOL_WORLD);
#endif
}

template <typename T>
void Parallel_Reduce::reduce_min(T& v)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_Type<T>::value, MPI_MIN, MPI_COMM_WORLD);
#endif
}

template <typename T>
void Parallel_Reduce::reduce_max(T& v)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_Type<T>::value, MPI_MAX, MPI_COMM_WORLD);
#endif
}

template <typename T>
void Parallel_Reduce::reduce_min_pool(const int& nproc_in_pool, T& v)
{
    if (nproc_in_pool == 1)
    {
        return;
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_Type<T>::value, MPI_MIN, POOL_WORLD);
#endif
}

template <typename T>
void Parallel_Reduce::reduce_max_pool(const int& nproc_in_pool, T& v)
{
    if (nproc_in_pool == 1)
    {
        return;
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_Type<T>::value, MPI_MAX, POOL_WORLD);
#endif
}

template void Parallel_Reduce::reduce_all<int>(int&);
template void Parallel_Reduce::reduce_all<double>(double&);
template void Parallel_Reduce::reduce_all<float>(float&);
template void Parallel_Reduce::reduce_all<std::complex<double>>(std::complex<double>&);
template void Parallel_Reduce::reduce_all<std::complex<float>>(std::complex<float>&);
template void Parallel_Reduce::reduce_all<long long>(long long&);

template void Parallel_Reduce::reduce_all<int>(int*, const int);
template void Parallel_Reduce::reduce_all<double>(double*, const int);
template void Parallel_Reduce::reduce_all<float>(float*, const int);
template void Parallel_Reduce::reduce_all<std::complex<double>>(std::complex<double>*, const int);
template void Parallel_Reduce::reduce_all<std::complex<float>>(std::complex<float>*, const int);
template void Parallel_Reduce::reduce_all<long long>(long long*, const int);

template void Parallel_Reduce::reduce_pool<float>(float&);
template void Parallel_Reduce::reduce_pool<double>(double&);
template void Parallel_Reduce::reduce_pool<std::complex<double>>(std::complex<double>&);

template void Parallel_Reduce::reduce_pool<int>(int*, const int);
template void Parallel_Reduce::reduce_pool<double>(double*, const int);
template void Parallel_Reduce::reduce_pool<std::complex<float>>(std::complex<float>*, const int);
template void Parallel_Reduce::reduce_pool<std::complex<double>>(std::complex<double>*, const int);

template void Parallel_Reduce::reduce_min<int>(int&);
template void Parallel_Reduce::reduce_min<float>(float&);
template void Parallel_Reduce::reduce_min<double>(double&);

template void Parallel_Reduce::reduce_max<float>(float&);
template void Parallel_Reduce::reduce_max<double>(double&);

template void Parallel_Reduce::reduce_min_pool<double>(const int&, double&);
template void Parallel_Reduce::reduce_max_pool<double>(const int&, double&);

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

void Parallel_Reduce::gather_int_all(int& v, int* all)
{
#ifdef __MPI
    assert(all != nullptr);
    MPI_Allgather(&v, 1, MPI_INT, all, 1, MPI_INT, MPI_COMM_WORLD);
#endif
    return;
}