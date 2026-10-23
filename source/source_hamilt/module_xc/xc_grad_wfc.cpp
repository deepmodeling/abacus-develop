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
    // gcar/kvec_c are stored as Vector3<double> on the host. Build
    // precision-matched host buffers and upload them, so single precision
    // never reinterprets a double tensor and does not depend on the
    // device-only PW_Basis_K caches being populated in the stress path.
    const int gcar_size = wfc_basis->nks * wfc_basis->npwk_max;
    std::vector<Real> gcar_host(3 * gcar_size);
    for (int i = 0; i < gcar_size; ++i)
    {
        gcar_host[3 * i] = static_cast<Real>(wfc_basis->gcar[i].x);
        gcar_host[3 * i + 1] = static_cast<Real>(wfc_basis->gcar[i].y);
        gcar_host[3 * i + 2] = static_cast<Real>(wfc_basis->gcar[i].z);
    }
    std::vector<Real> kvec_host(3 * wfc_basis->nks);
    for (int i = 0; i < wfc_basis->nks; ++i)
    {
        kvec_host[3 * i] = static_cast<Real>(wfc_basis->kvec_c[i].x);
        kvec_host[3 * i + 1] = static_cast<Real>(wfc_basis->kvec_c[i].y);
        kvec_host[3 * i + 2] = static_cast<Real>(wfc_basis->kvec_c[i].z);
    }
    auto gcar = ct::TensorMap(gcar_host.data(),
                              ct::DataTypeToEnum<Real>::value,
                              ct::DeviceType::CpuDevice,
                              {gcar_size, 3})
                    .to_device<ct_Device>();
    auto kvec_c = ct::TensorMap(kvec_host.data(),
                                ct::DataTypeToEnum<Real>::value,
                                ct::DeviceType::CpuDevice,
                                {wfc_basis->nks, 3})
                      .to_device<ct_Device>();

    grad_wfc<T, Device, Real>(ik, tpiba, wfc_basis, gcar.template data<Real>(), kvec_c.template data<Real>(), rhog, grad);
}

template <typename T, typename Device, typename Real>
void XC_Functional::grad_wfc(
    const int ik,
    const Real tpiba,
    const ModulePW::PW_Basis_K* wfc_basis,
    const Real* gcar,
    const Real* kvec_c,
    const T* rhog,
    T* grad)
{
    using ct_Device = typename ct::PsiToContainer<Device>::type;
    const int npw_k = wfc_basis->npwk[ik];

    auto porter = std::move(ct::Tensor(
        ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, {wfc_basis->nmaxgr}));

    auto xc_functional_grad_wfc_solver
        = hamilt::xc_functional_grad_wfc_op<T, Device>();

    for(int ipol=0; ipol<3; ipol++)
    {
        xc_functional_grad_wfc_solver(
            ik, ipol, npw_k, wfc_basis->npwk_max,
            tpiba,
            gcar,
            kvec_c,
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
template void XC_Functional::grad_wfc<std::complex<double>, base_device::DEVICE_CPU, double>(
    const int ik,
    const double tpiba,
    const ModulePW::PW_Basis_K* wfc_basis,
    const double* gcar,
    const double* kvec_c,
    const std::complex<double>* rhog,
    std::complex<double>* grad);
#if __CUDA || __ROCM
template void XC_Functional::grad_wfc<std::complex<double>, base_device::DEVICE_GPU, double>(
    const int ik,
    const double tpiba,
    const ModulePW::PW_Basis_K* wfc_basis,
    const std::complex<double>* rhog,
    std::complex<double>* grad);
template void XC_Functional::grad_wfc<std::complex<double>, base_device::DEVICE_GPU, double>(
    const int ik,
    const double tpiba,
    const ModulePW::PW_Basis_K* wfc_basis,
    const double* gcar,
    const double* kvec_c,
    const std::complex<double>* rhog,
    std::complex<double>* grad);
template void XC_Functional::grad_wfc<std::complex<float>, base_device::DEVICE_GPU, float>(
    const int ik,
    const float tpiba,
    const ModulePW::PW_Basis_K* wfc_basis,
    const std::complex<float>* rhog,
    std::complex<float>* grad);
template void XC_Functional::grad_wfc<std::complex<float>, base_device::DEVICE_GPU, float>(
    const int ik,
    const float tpiba,
    const ModulePW::PW_Basis_K* wfc_basis,
    const float* gcar,
    const float* kvec_c,
    const std::complex<float>* rhog,
    std::complex<float>* grad);
#endif
