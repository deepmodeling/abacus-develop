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
    template<typename T> void reduce_all(T& object) { return; };
    template<typename T> void reduce_all(T* object, const int n) { return; };
    template<typename T> void reduce_pool(T& object) { return; };
    template<typename T> void reduce_pool(T* object, const int n) { return; };

    template<>
    void reduce_all<int>(int& object) { return; };
    template<>
    void reduce_all<long long>(long long& object) { return; };
    template<>
    void reduce_all<double>(double& object) { return; };
    template<>
    void reduce_all<float>(float& object) { return; };
    template<>
    void reduce_all<std::complex<double>>(std::complex<double>& object) { return; };
    template<>
    void reduce_all<std::complex<float>>(std::complex<float>& object) { return; };

    template<>
    void reduce_all<int>(int* object, const int n) { return; };
    template<>
    void reduce_all<long long>(long long* object, const int n) { return; };
    template<>
    void reduce_all<double>(double* object, const int n) { return; };
    template<>
    void reduce_all<std::complex<double>>(std::complex<double>* object, const int n) { return; };
    template<>
    void reduce_all<std::complex<float>>(std::complex<float>* object, const int n) { return; };

    template<>
    void reduce_pool<float>(float& object) { return; };
    template<>
    void reduce_pool<double>(double& object) { return; };
    template<>
    void reduce_pool<std::complex<double>>(std::complex<double>& object) { return; };

    template<>
    void reduce_pool<int>(int* object, const int n) { return; };
    template<>
    void reduce_pool<double>(double* object, const int n) { return; };
    template<>
    void reduce_pool<std::complex<float>>(std::complex<float>* object, const int n) { return; };
    template<>
    void reduce_pool<std::complex<double>>(std::complex<double>* object, const int n) { return; };

    void reduce_or_all(bool& object) { return; };
    
    template <typename T>
    void reduce_max_all(T& object) { return; };
    template<> void reduce_max_all<double>(double& object) { return; };
    template<> void reduce_max_all<float>(float& object) { return; };
    template<> void reduce_max_all<int>(int& object) { return; };

    template <typename T>
    void reduce_min_all(T& object) { return; };
    template<> void reduce_min_all<double>(double& object) { return; };
    template<> void reduce_min_all<float>(float& object) { return; };
    template<> void reduce_min_all<int>(int& object) { return; };

    void reduce_max_pool(int* object, const int n) { return; };
    void reduce_min_pool(double& object) { return; };

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