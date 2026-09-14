// Gradient of the wavefunction in reciprocal space
// (XC_Functional::grad_wfc), used by stress_func_mgga.cpp.

#include "xc_functional.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <ATen/core/tensor_types.h>
#include <source_hamilt/module_xc/kernels/xc_functional_op.h>

template <typename T, typename Device, typename Real>
void XC_Functional::grad_wfc(
    const int ik,
    const Real tpiba,
    const ModulePW::PW_Basis_K* wfc_basis,
    const T* rhog,
    T* grad)
{
    using ct_Device = typename ct::PsiToContainer<Device>::type;
    const int npw_k = wfc_basis->npwk[ik];

    auto porter = std::move(ct::Tensor(
        ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, {wfc_basis->nmaxgr}));
    auto gcar = ct::TensorMap(
        &wfc_basis->gcar[0][0], ct::DataType::DT_DOUBLE, ct::DeviceType::CpuDevice, {wfc_basis->nks * wfc_basis->npwk_max, 3}).to_device<ct_Device>();
    auto kvec_c = ct::TensorMap(
        &wfc_basis->kvec_c[0][0],ct::DataType::DT_DOUBLE, ct::DeviceType::CpuDevice, {wfc_basis->nks, 3}).to_device<ct_Device>();

    auto xc_functional_grad_wfc_solver
        = hamilt::xc_functional_grad_wfc_op<T, Device>();

    for(int ipol=0; ipol<3; ipol++)
    {
        xc_functional_grad_wfc_solver(
            ik, ipol, npw_k, wfc_basis->npwk_max,
            tpiba,
            gcar.template data<Real>(),
            kvec_c.template data<Real>(),
            rhog, porter.data<T>());

        // bring the gdr from G --> R
        Device * ctx = nullptr;
        wfc_basis->recip_to_real(ctx, porter.data<T>(), porter.data<T>(), ik);

        xc_functional_grad_wfc_solver(
            ipol, wfc_basis->nrxx,
            porter.data<T>(), grad);
    }
}

template void XC_Functional::grad_wfc<std::complex<double>, base_device::DEVICE_CPU, double>(
    const int ik,
    const double tpiba,
    const ModulePW::PW_Basis_K* wfc_basis,
    const std::complex<double>* rhog,
    std::complex<double>* grad);
#if __CUDA || __ROCM
template void XC_Functional::grad_wfc<std::complex<double>, base_device::DEVICE_GPU, double>(
    const int ik,
    const double tpiba,
    const ModulePW::PW_Basis_K* wfc_basis,
    const std::complex<double>* rhog,
    std::complex<double>* grad);
#endif
