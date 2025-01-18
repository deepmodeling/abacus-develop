#ifndef STO_CHE_H
#define STO_CHE_H
#include "module_base/math_chebyshev.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_base/module_container/ATen/kernels/blas.h"

template <typename REAL, typename Device = base_device::DEVICE_CPU>
class StoChe
{
  public:
    StoChe(const int& nche, const int& method, const REAL& emax_sto, const REAL& emin_sto);
    ~StoChe();

  public:
    int nche = 0;               ///< order of Chebyshev expansion
    REAL* spolyv = nullptr;     ///< [Device] coefficients of Chebyshev expansion
    REAL* spolyv_cpu = nullptr; ///< [CPU] coefficients of Chebyshev expansion
    int method_sto = 0;         ///< method for the stochastic calculation

    // Chebyshev expansion
    // It stores the plan of FFTW and should be initialized at the beginning of the calculation
    ModuleBase::Chebyshev<REAL, Device>* p_che = nullptr;

    REAL emax_sto = 0.0; ///< maximum energy for normalization
    REAL emin_sto = 0.0; ///< minimum energy for normalization

  private:
    Device* ctx = {};
    using resmem_var_op = base_device::memory::resize_memory_op<REAL, Device>;
    using syncmem_var_d2h_op = base_device::memory::synchronize_memory_op<REAL, base_device::DEVICE_CPU, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<REAL, Device>;
};

/**
 * @brief calculate v^T*M*v
 *
 * @param v v
 * @param M M
 * @param n the dimension of v
 * @return REAL
 */
template <typename REAL, typename Device>
REAL vTMv(const REAL* v, const REAL* M, const int n)
{
    Device* ctx = {};
    base_device::DEVICE_CPU* cpu_ctx = {};
    using ct_Device = typename container::PsiToContainer<Device>::type;
    const char normal = 'N';
    const REAL one = 1;
    const int inc = 1;
    const REAL zero = 0;
    REAL* y = nullptr;
    base_device::memory::resize_memory_op<REAL, Device>()(y, n);
    hsolver::gemv_op<REAL, Device>()(ctx, normal, n, n, &one, M, n, v, inc, &zero, y, inc);
    REAL result = 0;
    REAL* dot_device = nullptr;
    base_device::memory::resize_memory_op<REAL, Device>()(dot_device, 1);
    container::kernels::blas_dot<REAL, ct_Device>()(n, y, 1, v, 1, dot_device);
    base_device::memory::synchronize_memory_op<REAL, base_device::DEVICE_CPU, Device>()(&result,
                                                                                        dot_device,
                                                                                        1);
    base_device::memory::delete_memory_op<REAL, Device>()(y);
    base_device::memory::delete_memory_op<REAL, Device>()(dot_device);
    return result;
}

#endif