#ifdef __MPI
#include "mpi.h"
#endif
#include "depend_mock.h"
#include <complex>

namespace GlobalV
{ 
    std::ofstream ofs_running;
}
#ifdef __MPI
MPI_Comm POOL_WORLD;
namespace Parallel_Reduce
{
    template<typename T> void reduce_all(T& object);
    template<typename T> void reduce_all(T* object, const int n);
    template<typename T> void reduce_pool(T& object);
    template<typename T> void reduce_pool(T* object, const int n);

    template<>
    void reduce_all<int>(int& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); }
    template<>
    void reduce_all<long long>(long long& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD); }
    template<>
    void reduce_all<double>(double& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); }
    template<>
    void reduce_all<float>(float& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); }
    template<>
    void reduce_all<std::complex<double>>(std::complex<double>& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD); }
    template<>
    void reduce_all<std::complex<float>>(std::complex<float>& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_C_FLOAT_COMPLEX, MPI_SUM, MPI_COMM_WORLD); }

    template<>
    void reduce_all<int>(int* object, const int n) { MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD); }
    template<>
    void reduce_all<long long>(long long* object, const int n) { MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD); }
    template<>
    void reduce_all<double>(double* object, const int n) { MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); }
    template<>
    void reduce_all<std::complex<double>>(std::complex<double>* object, const int n) { MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD); }
    template<>
    void reduce_all<std::complex<float>>(std::complex<float>* object, const int n) { MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_C_FLOAT_COMPLEX, MPI_SUM, MPI_COMM_WORLD); }

    template<>
    void reduce_pool<float>(float& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_FLOAT, MPI_SUM, POOL_WORLD); }
    template<>
    void reduce_pool<double>(double& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_SUM, POOL_WORLD); }
    template<>
    void reduce_pool<std::complex<double>>(std::complex<double>& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD); }

    template<>
    void reduce_pool<int>(int* object, const int n) { MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_SUM, POOL_WORLD); }
    template<>
    void reduce_pool<double>(double* object, const int n) { MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, POOL_WORLD); }

    void reduce_max_pool(int* object, const int n) { MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_INT, MPI_MAX, POOL_WORLD); }
    void reduce_min_pool(double& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_MIN, POOL_WORLD); }

    // Other stubs can remain as is if not used or if they don't break logic
    void reduce_or_all(bool& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD); }
    
    template <typename T>
    void reduce_max_all(T& object);
    template<> void reduce_max_all<double>(double& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); }
    template<> void reduce_max_all<float>(float& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD); }
    template<> void reduce_max_all<int>(int& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); }

    template <typename T>
    void reduce_min_all(T& object);
    template<> void reduce_min_all<double>(double& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); }
    template<> void reduce_min_all<float>(float& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD); }
    template<> void reduce_min_all<int>(int& object) { MPI_Allreduce(MPI_IN_PLACE, &object, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD); }

    void reduce_or_bp(bool& object) { return; };

    void reduce_double_bgroup(double& object) { return; };
    void reduce_double_bgroup(double* object, const int n) { return; };

    void reduce_double_bp(double& object) { return; };
    void reduce_double_bp(double* object, const int n) { return; };

    void reduce_double_kp(double* object, const int n) { return; };

    void reduce_double_allpool(const int& npool, const int& nproc_in_pool, double& object) { return; };
    void reduce_double_allpool(const int& npool, const int& nproc_in_pool, double* object, const int n) { return; };

    void gather_min_int_all(const int& nproc, int& v) { return; };
    void gather_max_double_all(const int& nproc, double& v) { return; };
    void gather_min_double_all(const int& nproc, double& v) { return; };
    void gather_max_double_pool(const int& nproc_in_pool, double& v) { return; };
    void gather_min_double_pool(const int& nproc_in_pool, double& v) { return; };
    void gather_int_all(int& v, int* all) { return; };
}
#endif