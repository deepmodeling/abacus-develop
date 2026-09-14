#ifndef STO_CHE_H
#define STO_CHE_H
#include "source_base/math_chebyshev.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_container/ATen/kernels/blas.h"

#include <memory>
#include <vector>

template <typename REAL, typename Device = base_device::DEVICE_CPU>
class StoChe
{
  private:
#ifdef __DSP
    using resmem_var_op = base_device::memory::resize_memory_op_mt<REAL, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op_mt<REAL, Device>;
#else
    using resmem_var_op = base_device::memory::resize_memory_op<REAL, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<REAL, Device>;
#endif
    using syncmem_var_d2h_op = base_device::memory::synchronize_memory_op<REAL, base_device::DEVICE_CPU, Device>;

    /// @brief deleter for device-side spolyv buffer
    struct delmem_var_deleter
    {
        void operator()(REAL* p) const
        {
            if (p)
            {
                delmem_var_op()(p);
            }
        }
    };

  public:
    StoChe() = default;
    ~StoChe() = default;

    StoChe(const StoChe&) = delete;
    StoChe& operator=(const StoChe&) = delete;
    StoChe(StoChe&&) noexcept = default;
    StoChe& operator=(StoChe&&) noexcept = default;

    /// @brief (Re)allocate Chebyshev expansion buffers. Safe to call multiple times.
    void init(const int& nche, const int& method, const REAL& emax_sto, const REAL& emin_sto);

    int nche = 0;                                      ///< order of Chebyshev expansion
    std::unique_ptr<REAL, delmem_var_deleter> spolyv;  ///< [Device] coefficients of Chebyshev expansion
    std::vector<REAL> spolyv_cpu;                       ///< [CPU] coefficients of Chebyshev expansion (method==1 only)
    int method_sto = 0;                                 ///< method for the stochastic calculation

    /// Chebyshev expansion. Stores the plan of FFTW and should be initialized at the beginning of the calculation
    std::unique_ptr<ModuleBase::Chebyshev<REAL, Device>> p_che;

    REAL emax_sto = 0.0; ///< maximum energy for normalization
    REAL emin_sto = 0.0; ///< minimum energy for normalization

  private:
    Device* ctx = {};
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
    ModuleBase::gemv_op<REAL, Device>()(normal, n, n, &one, M, n, v, inc, &zero, y, inc);
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
