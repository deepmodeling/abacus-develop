#include "module_base/timer.h"
#include "module_basis/module_pw/kernels/pw_op.h"
#include "pw_basis_k.h"
#include "pw_gatherscatter.h"

#include <cassert>
#include <complex>
#if defined (__DSP)
namespace ModulePW
{
template <typename FPTYPE>
void PW_Basis_K::real2recip_dsp(const std::complex<FPTYPE>* in,
                                std::complex<FPTYPE>* out,
                                const int ik,
                                const bool add,
                                const FPTYPE factor) const
{
    const base_device::DEVICE_CPU* ctx;
    const base_device::DEVICE_GPU* gpux;
    assert(this->gamma_only == false);
    auto* auxr = this->fft_bundle.get_auxr_3d_data<double>();

    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];
    // copy the in into the auxr with complex<double>
    memcpy(auxr, in, this->nrxx * 2 * 8);

    // 3d fft
    this->fft_bundle.resource_handler(1);
    this->fft_bundle.fft3D_forward(auxr, 
                                   auxr);
    this->fft_bundle.resource_handler(0);
    // copy the result from the auxr to the out ,while consider the add
    set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_CPU>()(npw_k,
                                                                   this->nxyz,
                                                                   add,
                                                                   factor,
                                                                   this->ig2ixyz_k_cpu.data() + startig,
                                                                   auxr,
                                                                   out);
}
template <typename FPTYPE>
void PW_Basis_K::recip2real_dsp(const std::complex<FPTYPE>* in,
                                std::complex<FPTYPE>* out,
                                const int ik,
                                const bool add,
                                const FPTYPE factor) const
{
    assert(this->gamma_only == false);
    const base_device::DEVICE_CPU* ctx;
    const base_device::DEVICE_GPU* gpux;
    // memset the auxr of 0 in the auxr,here the len of the auxr is nxyz
    auto* auxr = this->fft_bundle.get_auxr_3d_data<double>();
    memset(auxr, 0, this->nxyz * 2 * 8);

    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];
    // copy the mapping form the type of stick to the 3dfft
    set_3d_fft_box_op<double, base_device::DEVICE_CPU>()(npw_k, this->ig2ixyz_k_cpu.data() + startig, in, auxr);
    // use 3d fft backward
    this->fft_bundle.resource_handler(1);
    this->fft_bundle.fft3D_backward(auxr, auxr);
    this->fft_bundle.resource_handler(0);
    if (add)
    {
        const int one = 1;
        const std::complex<double> factor1 = std::complex<double>(factor, 0);
        zaxpy_(&nrxx, &factor1, auxr, &one, out, &one);
    }
    else
    {
        memcpy(out, auxr, nrxx * 2 * 8);
    }
}
template <>
void PW_Basis_K::convolution(const base_device::DEVICE_CPU* ctx,
                             const int ik,
                             const int size,
                             const std::complex<float>* input,
                             const float* input1,
                             std::complex<float>* output,
                             const bool add,
                             const float factor) const
{
}

template <>
void PW_Basis_K::convolution(const base_device::DEVICE_CPU* ctx,
                             const int ik,
                             const int size,
                             const std::complex<double>* input,
                             const double* input1,
                             std::complex<double>* output,
                             const bool add,
                             const double factor) const
{
    ModuleBase::timer::tick(this->classname, "convolution");

    assert(this->gamma_only == false);
    const base_device::DEVICE_GPU* gpux;
    // memset the auxr of 0 in the auxr,here the len of the auxr is nxyz
    auto* auxr = this->fft_bundle.get_auxr_3d_data<double>();
    memset(auxr, 0, this->nxyz * 2 * 8);
    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];

    // copy the mapping form the type of stick to the 3dfft
    set_3d_fft_box_op<double, base_device::DEVICE_CPU>()(npw_k, this->ig2ixyz_k_cpu.data() + startig, input, auxr);

    // use 3d fft backward
    this->fft_bundle.fft3D_backward(auxr, auxr);

    for (int ir = 0; ir < size; ir++)
    {
        auxr[ir] *= input1[ir];
    }

    // 3d fft
    this->fft_bundle.fft3D_forward(auxr, auxr);
    // copy the result from the auxr to the out ,while consider the add
    set_real_to_recip_output_op<double, base_device::DEVICE_CPU>()(npw_k,
                                                                   this->nxyz,
                                                                   add,
                                                                   factor,
                                                                   this->ig2ixyz_k_cpu.data() + startig,
                                                                   auxr,
                                                                   output);
    ModuleBase::timer::tick(this->classname, "convolution");
}

// template void PW_Basis_K::real2recip_dsp<float>(const std::complex<float>* in,
//                                             std::complex<float>* out,
//                                             const int ik,
//                                             const bool add,
//                                             const float factor) const; // in:(nplane,nx*ny)  ; out(nz, ns)
// template void PW_Basis_K::recip2real_dsp<float>(const std::complex<float>* in,
//                                             std::complex<float>* out,
//                                             const int ik,
//                                             const bool add,
//                                             const float factor) const; // in:(nz, ns)  ; out(nplane,nx*ny)

template void PW_Basis_K::real2recip_dsp<double>(const std::complex<double>* in,
                                                 std::complex<double>* out,
                                                 const int ik,
                                                 const bool add,
                                                 const double factor) const; // in:(nplane,nx*ny)  ; out(nz, ns)
template void PW_Basis_K::recip2real_dsp<double>(const std::complex<double>* in,
                                                 std::complex<double>* out,
                                                 const int ik,
                                                 const bool add,
                                                 const double factor) const;
} // namespace ModulePW
#endif
