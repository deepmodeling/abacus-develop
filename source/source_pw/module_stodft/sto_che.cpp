#include "sto_che.h"
#include "source_base/module_device/device.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_container/ATen/kernels/blas.h"

template <typename REAL, typename Device>
void StoChe<REAL, Device>::init(const int& nche, const int& method, const REAL& emax_sto, const REAL& emin_sto)
{
    // release old resources first (safe for default-constructed state too)
    this->p_che.reset();
    this->spolyv.reset();
    this->spolyv_cpu.clear();

    this->nche = nche;
    this->method_sto = method;
    this->p_che = std::make_unique<ModuleBase::Chebyshev<REAL, Device>>(nche);
    if (method == 1)
    {
        REAL* spolyv_ptr = nullptr;
        resmem_var_op()(spolyv_ptr, nche);
        this->spolyv.reset(spolyv_ptr);
        this->spolyv_cpu.resize(nche);
    }
    else
    {
        REAL* spolyv_ptr = nullptr;
        resmem_var_op()(spolyv_ptr, nche * nche);
        this->spolyv.reset(spolyv_ptr);
    }

    this->emax_sto = emax_sto;
    this->emin_sto = emin_sto;
}

template class StoChe<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class StoChe<double, base_device::DEVICE_GPU>;
#endif
