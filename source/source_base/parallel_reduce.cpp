#include "parallel_reduce.h"
#include "parallel_comm.h"
#include <vector>
#include <cassert>

// Helper to safely check MPI status for WORLD comm
// optimize: use a global static variable to avoid the overhead of function-static thread-safety checks
static bool is_mpi_active = false;
inline static bool is_mpi_initialized() {
#ifdef __MPI
    if (is_mpi_active) return true;
    int flag = 0;
    MPI_Initialized(&flag);
    if (flag) is_mpi_active = true;
    return flag != 0;
#else
    return false;
#endif
}

template <>
void Parallel_Reduce::reduce_all<int>(int& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<long long>(long long& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

void Parallel_Reduce::reduce_int_diag(int& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_INT, MPI_SUM, DIAG_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<double>(double& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<float>(float& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<int>(int* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<long long>(long long* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

void Parallel_Reduce::reduce_int_grid(int* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_SUM, GRID_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<double>(double* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

void Parallel_Reduce::reduce_double_grid(double* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, GRID_WORLD);
#endif
    return;
}

void Parallel_Reduce::reduce_double_diag(double* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, DIAG_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<float>(float& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_FLOAT, MPI_SUM, POOL_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<double>(double& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_SUM, POOL_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<double>(double* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, POOL_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<int>(int& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_INT, MPI_SUM, POOL_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<int>(int* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_SUM, POOL_WORLD);
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
    if (!is_mpi_initialized()) return;
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
    if (!is_mpi_initialized()) return;
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
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

// LiuXh add 2019-07-16
template <>
void Parallel_Reduce::reduce_all<std::complex<double>>(std::complex<double>* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}


template <>
void Parallel_Reduce::reduce_all<std::complex<float>>(std::complex<float>& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_C_FLOAT_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<std::complex<double>>(std::complex<double>& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<std::complex<float>>(std::complex<float>* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_C_FLOAT_COMPLEX, MPI_SUM, POOL_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_pool<std::complex<double>>(std::complex<double>* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_all<std::complex<float>>(std::complex<float>* object, int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_C_FLOAT_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}

template <>
void Parallel_Reduce::reduce_max_all<double>(double& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}

template<>
void Parallel_Reduce::reduce_max_pool<double>(double& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_MAX, POOL_WORLD);
#endif
}

template<>
void Parallel_Reduce::reduce_max_pool<int>(int& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_INT, MPI_MAX, POOL_WORLD);
#endif
}

template<>
void Parallel_Reduce::reduce_max_pool<int>(int* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_MAX, POOL_WORLD);
#endif
}

template<>
void Parallel_Reduce::reduce_min_pool<double>(double& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_MIN, POOL_WORLD);
#endif
}

template<>
void Parallel_Reduce::reduce_bgroup<double>(double* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, INT_BGROUP);
#endif
}

template<>
void Parallel_Reduce::reduce_kp<int>(int* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    if (KP_WORLD != MPI_COMM_NULL)
    {
        MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_SUM, KP_WORLD);
    }
#endif
}

template<>
void Parallel_Reduce::reduce_bp<int>(int* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_SUM, BP_WORLD);
#endif
}

template<>
void Parallel_Reduce::reduce_bgroup<int>(int* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_SUM, INT_BGROUP);
#endif
}

template <>
void Parallel_Reduce::reduce_min_all<int>(int& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif
}

template <>
void Parallel_Reduce::reduce_min_all<double>(double& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
}

void Parallel_Reduce::gather_int_all(int& v, int* all)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allgather(&v, 1, MPI_INT, all, 1, MPI_INT, MPI_COMM_WORLD);
#else
    all[0] = v;
#endif
}

template<>
void Parallel_Reduce::reduce_kp<double>(double* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    if (KP_WORLD != MPI_COMM_NULL)
    {
        MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, KP_WORLD);
    }
#endif
}

template<>
void Parallel_Reduce::reduce_bp<double>(double* object, const int n)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, BP_WORLD);
#endif
}

template<>
void Parallel_Reduce::reduce_or_bp<bool>(bool& object)
{
#ifdef __MPI
    if (!is_mpi_initialized()) return;
    int v_int = object ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &v_int, 1, MPI_INT, MPI_LOR, BP_WORLD);
    object = (v_int != 0);
#endif
}
