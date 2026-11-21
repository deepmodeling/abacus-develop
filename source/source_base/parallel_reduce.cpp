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
const MPI_Datatype Parallel_Reduce::MPI_Type<unsigned long long>::value = MPI_UNSIGNED_LONG_LONG;
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

void Parallel_Reduce::reduce_max(double* values, const int count)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, values, count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
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

template void Parallel_Reduce::reduce_max<int>(int&);
template void Parallel_Reduce::reduce_max<float>(float&);
template void Parallel_Reduce::reduce_max<double>(double&);

template void Parallel_Reduce::reduce_min_pool<double>(const int&, double&);
template void Parallel_Reduce::reduce_max_pool<double>(const int&, double&);

#ifdef __MPI
template <typename T>
void Parallel_Reduce::reduce_all(const MPI_Comm mpi_comm, T& object)
{
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_Type<T>::value, MPI_SUM, mpi_comm);
}

template <typename T>
void Parallel_Reduce::reduce_all(const MPI_Comm mpi_comm, T* object, const int n)
{
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_Type<T>::value, MPI_SUM, mpi_comm);
}

template <typename T>
void Parallel_Reduce::reduce_max(const MPI_Comm mpi_comm, T& object)
{
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_Type<T>::value, MPI_MAX, mpi_comm);
}

template <typename T>
void Parallel_Reduce::reduce_min(const MPI_Comm mpi_comm, T& object)
{
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_Type<T>::value, MPI_MIN, mpi_comm);
}

void Parallel_Reduce::reduce_max_loc(const MPI_Comm mpi_comm, double& value, int& rank)
{
    struct MaxLoc
    {
        double value;
        int rank;
    };
    MaxLoc local{value, rank};
    MaxLoc global{};
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, mpi_comm);
    value = global.value;
    rank = global.rank;
}

void Parallel_Reduce::gather_int_all(const MPI_Comm mpi_comm, int& v, int* all)
{
    assert(all != nullptr);
    MPI_Allgather(&v, 1, MPI_INT, all, 1, MPI_INT, mpi_comm);
}

void Parallel_Reduce::gatherv_int_all(const MPI_Comm mpi_comm,
                                      const int* send_buffer,
                                      const int send_count,
                                      int* receive_buffer,
                                      const int* receive_counts,
                                      const int* receive_displacements)
{
    assert(receive_buffer != nullptr);
    MPI_Allgatherv(send_buffer, send_count, MPI_INT, receive_buffer, receive_counts, receive_displacements,
                   MPI_INT, mpi_comm);
}

template void Parallel_Reduce::reduce_all<int>(const MPI_Comm, int&);
template void Parallel_Reduce::reduce_all<double>(const MPI_Comm, double&);
template void Parallel_Reduce::reduce_all<unsigned long long>(const MPI_Comm, unsigned long long&);

template void Parallel_Reduce::reduce_all<int>(const MPI_Comm, int*, const int);
template void Parallel_Reduce::reduce_all<double>(const MPI_Comm, double*, const int);
template void Parallel_Reduce::reduce_all<unsigned long long>(const MPI_Comm, unsigned long long*, const int);

template void Parallel_Reduce::reduce_max<int>(const MPI_Comm, int&);
template void Parallel_Reduce::reduce_max<double>(const MPI_Comm, double&);
template void Parallel_Reduce::reduce_max<unsigned long long>(const MPI_Comm, unsigned long long&);

template void Parallel_Reduce::reduce_min<int>(const MPI_Comm, int&);
template void Parallel_Reduce::reduce_min<double>(const MPI_Comm, double&);
template void Parallel_Reduce::reduce_min<unsigned long long>(const MPI_Comm, unsigned long long&);
#endif


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
