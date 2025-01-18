#include "sto_che.h"
#include "module_base/blas_connector.h"
#include "module_base/module_device/device.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_base/module_container/ATen/kernels/blas.h"

template <typename REAL, typename Device>
StoChe<REAL, Device>::~StoChe()
{
    delete p_che;
    delete[] spolyv_cpu;
    delmem_var_op()(spolyv);
}

template <typename REAL, typename Device>
StoChe<REAL, Device>::StoChe(const int& nche, const int& method, const REAL& emax_sto, const REAL& emin_sto)
{
    this->nche = nche;
    this->method_sto = method;
    p_che = new ModuleBase::Chebyshev<REAL, Device>(nche);
    if (method == 1)
    {
        resmem_var_op()(spolyv, nche);
        spolyv_cpu = new REAL[nche];
    }
    else
    {
        resmem_var_op()(spolyv, nche * nche);
    }

    this->emax_sto = emax_sto;
    this->emin_sto = emin_sto;
}

template class StoChe<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class StoChe<double, base_device::DEVICE_GPU>;
#endif