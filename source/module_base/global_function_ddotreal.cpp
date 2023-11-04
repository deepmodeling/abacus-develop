#include "blas_connector.h"
#include "global_function.h"
#include "module_base/parallel_reduce.h"

namespace ModuleBase
{
namespace GlobalFunc
{

template double ddot_real(const int& dim,
                          const std::complex<double>* psi_L,
                          const std::complex<double>* psi_R,
                          const bool reduce);
template float ddot_real(const int& dim,
                         const std::complex<float>* psi_L,
                         const std::complex<float>* psi_R,
                         const bool reduce);

template <typename T>
T ddot_real(const int& dim, const std::complex<T>* psi_L, const std::complex<T>* psi_R, const bool reduce)
{
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // qianrui modify 2021-3-14
    // Note that  ddot_(2*dim,a,1,b,1) = REAL( zdotc_(dim,a,1,b,1) )
    int dim2 = 2 * dim;
    T *pL, *pR;
    pL = (T*)psi_L;
    pR = (T*)psi_R;
    T result = BlasConnector::dot(dim2, pL, 1, pR, 1);
    if (reduce)
        Parallel_Reduce::reduce_pool(result);
    return result;
    //======================================================================
    /*std::complex<double> result(0,0);
    for (int i=0;i<dim;i++)
    {
        result += conj( psi_L[i] ) * psi_R[i];
    }
    Parallel_Reduce::reduce_complex_double_pool( result );
    return result.real();*/
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
}

} // namespace GlobalFunc
} // namespace ModuleBase