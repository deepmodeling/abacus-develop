#include "pw_basis_k.h"
#include "source_base/timer.h"
#include "source_basis/module_pw/kernels/pw_op.h"
#include "source_pw/module_pwdft/kernels/veff_op.h"
#include "pw_gatherscatter.h"
namespace ModulePW
{
template <typename FPTYPE>
void PW_Basis_K::convolution_cpu(const int ik,
                             const int size,
                             const std::complex<FPTYPE>* input,
                             const FPTYPE* input1,
                             std::complex<FPTYPE>* output,
                             const bool add,
                             const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "convolution");
    assert(this->gamma_only == false);
    // ModuleBase::GlobalFunc::ZEROS(fft_bundle.get_auxg_data<double>(), this->nst * this->nz);
    // memset the auxr of 0 in the auxr,here the len of the auxr is nxyz
    auto* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
    auto* auxr=this->fft_bundle.get_auxr_data<FPTYPE>();

    memset(auxg, 0, this->nst * this->nz * sizeof(FPTYPE)*2);
    const int startig = ik * this->npwk_max;
    const int npwk = this->npwk[ik];

    // copy the mapping form the type of stick to the 3dfft
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
    #endif
    for (int igl = 0; igl < npwk; ++igl)
    {
        auxg[this->igl2isz_k[igl + startig]] = input[igl];
    }

    // use 3d fft backward
    this->fft_bundle.fftzbac(auxg, auxg);

    this->gathers_scatterp(auxg, auxr);

    this->fft_bundle.fftxybac(auxr, auxr);
    for (int ir = 0; ir < size; ir++)
    {
        auxr[ir] *= input1[ir];
    }

    // 3d fft
    this->fft_bundle.fftxyfor(auxr, auxr);

    this->gatherp_scatters(auxr, auxg);

    this->fft_bundle.fftzfor(auxg, auxg);
    // copy the result from the auxr to the out ,while consider the add
        FPTYPE tmpfac = factor / FPTYPE(this->nxyz);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
        for (int igl = 0; igl < npwk; ++igl)
        {
            output[igl] += tmpfac * auxg[this->igl2isz_k[igl + startig]];
        }
    ModuleBase::timer::tick(this->classname, "convolution");
}

template <typename FPTYPE>
void PW_Basis_K::convolution_cpu(const int ik,
                             const int size,
                             const int max_npw,
                             const std::complex<FPTYPE>* input,
                             std::complex<FPTYPE>* tmp,
                             std::complex<FPTYPE>* input1,
                             std::complex<FPTYPE>* output,
                             const bool add,
                             const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "convolution");
    assert(this->gamma_only == false);
    base_device::DEVICE_CPU* cpu_ctx ;
    memset(tmp, 0, 2*size * 2 * sizeof(FPTYPE));
    const int startig = ik * this->npwk_max;
    const int startig2 = ik * this->npwk_max;
    const int npwk = this->npwk[ik];
    auto *augr = tmp;
    auto *augr1 = &tmp[size];
    auto *augx = this->fft_bundle.get_auxg_data<FPTYPE>();
    auto *augx1 = this->fft_bundle.get_auxr_data<FPTYPE>();
    memset(augx, 0, this->nst * this->nz * sizeof(FPTYPE)*2);
    memset(augx1, 0, this->nst * this->nz * sizeof(FPTYPE)*2);
    for (int igl = 0; igl < npwk; ++igl)
    {
        const int idx=this->igl2isz_k[igl + startig];
        augr[idx] = input[igl];
        augr1[idx] = input[igl+max_npw];
    }
    // use 3d fft backward
    this->fft_bundle.fftzbac(augr, augr);
    this->fft_bundle.fftzbac(augr1, augr1);

    this->gathers_scatterp(augr, augx);
    this->gathers_scatterp(augr1, augx1);

    this->fft_bundle.fftxybac(augx, augx);
    this->fft_bundle.fftxybac(augx1, augx1);
    hamilt::veff_pw_op<FPTYPE,base_device::DEVICE_CPU>()(cpu_ctx, size, augx, augx1, input1);
    this->fft_bundle.fftxyfor(augx, augx);
    this->fft_bundle.fftxyfor(augx1, augx1);

    this->gatherp_scatters(augx, augr);
    this->gatherp_scatters(augx1,augr1);
    
    this->fft_bundle.fftzfor(augr, augr);
    this->fft_bundle.fftzfor(augr1,augr1);

    FPTYPE tmpfac = factor / FPTYPE(this->nxyz);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
    for (int igl = 0; igl < npwk; ++igl)
    {
        const int idx=this->igl2isz_k[igl + startig];
        output[igl] += tmpfac * augr[idx];
        output[igl+max_npw] += tmpfac * augr1[idx];
    }
    ModuleBase::timer::tick(this->classname, "convolution");
    }

#if (defined(__CUDA) || defined(__ROCM))
template <typename FPTYPE>
void PW_Basis_K::convolution_gpu(const int ik,
                             const int size,
                             const std::complex<FPTYPE>* input,
                             const FPTYPE* input1,
                             std::complex<FPTYPE>* output,
                             const bool add,
                             const FPTYPE factor) const
{       
    ModuleBase::timer::tick(this->classname, "convolution");

    assert(this->gamma_only == false);
    const base_device::DEVICE_GPU* gpux;
    // memset the auxr of 0 in the auxr,here the len of the auxr is nxyz
    
    base_device::memory::set_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(
        this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
        0,
        this->nxyz);
    auto* auxr = this->fft_bundle.get_auxr_3d_data<FPTYPE>();
    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];

    // copy the mapping form the type of stick to the 3dfft
    set_3d_fft_box_op<FPTYPE, base_device::DEVICE_GPU>()(npw_k, this->ig2ixyz_k + startig, input, auxr);

    // use 3d fft backward
    this->fft_bundle.fft3D_backward(auxr, auxr);
    
    hamilt::veff_pw_op<FPTYPE,base_device::DEVICE_GPU>()(gpux,size,auxr,input1);

    // 3d fft
    this->fft_bundle.fft3D_forward(auxr, auxr);
    // copy the result from the auxr to the out ,while consider the add
    set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>()(npw_k,
                                                                   this->nxyz,
                                                                   add,
                                                                   factor,
                                                                   this->ig2ixyz_k + startig,
                                                                   auxr,
                                                                   output);
    ModuleBase::timer::tick(this->classname, "convolution");
}
template <typename FPTYPE>
void PW_Basis_K::convolution_gpu(const int ik,
                             const int size,
                             const int max_npw,
                             const std::complex<FPTYPE>* input,
                             std::complex<FPTYPE>* tmp,
                             std::complex<FPTYPE>* input1,
                             std::complex<FPTYPE>* output,   
                             const bool add,
                             const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "convolution");
    assert(this->gamma_only == false);
    base_device::DEVICE_GPU* gpu_ctx;
    // memset the auxr of 0 in the auxr,here the len of the auxr is nxyz
    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];
    auto *auxg = this->fft_bundle.get_auxr_3d_data<FPTYPE>();
    auto *auxg1 = &this->fft_bundle.get_auxr_3d_data<FPTYPE>()[size];
    base_device::memory::set_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(
        auxg,
        0,
        this->nxyz);
    base_device::memory::set_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(
        auxg1,
        0,
        this->nxyz);
    set_3d_fft_box_op<FPTYPE, base_device::DEVICE_GPU>()(npw_k,
                                                         this->ig2ixyz_k + startig,
                                                         input,
                                                         auxg);
    set_3d_fft_box_op<FPTYPE, base_device::DEVICE_GPU>()(npw_k,
                                                         this->ig2ixyz_k + startig,
                                                         input + max_npw,
                                                         auxg1);
    this->fft_bundle.fft3D_backward(auxg, auxg);
    this->fft_bundle.fft3D_backward(auxg1, auxg1);

    hamilt::veff_pw_op<FPTYPE,base_device::DEVICE_GPU>()(gpu_ctx,size,auxg,auxg1,input1);

    this->fft_bundle.fft3D_forward(auxg, auxg);
    this->fft_bundle.fft3D_forward(auxg1, auxg1);
    // use 3d fft backward
    set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>()(npw_k,
                                                                   this->nxyz,
                                                                   true,
                                                                   factor,
                                                                   this->ig2ixyz_k + startig,
                                                                   auxg,
                                                                   output);
    set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>()(npw_k,
                                                                   this->nxyz,
                                                                   true,
                                                                   factor,
                                                                   this->ig2ixyz_k + startig,
                                                                   auxg1,
                                                                   output + max_npw);
    ModuleBase::timer::tick(this->classname, "convolution");
}

template void PW_Basis_K::convolution_gpu<float>(const int ik,
                                                 const int size,
                                                 const std::complex<float>* input,
                                                 const float* input1,
                                                 std::complex<float>* output,
                                                 const bool add,
                                                 const float factor) const;
template void PW_Basis_K::convolution_gpu<double>(const int ik,
                                                  const int size,
                                                  const std::complex<double>* input,
                                                  const double* input1,
                                                  std::complex<double>* output,
                                                  const bool add,
                                                  const double factor) const;
template void PW_Basis_K::convolution_gpu<float>(const int ik,
                                                 const int size,
                                                 const int max_npw,
                                                 const std::complex<float>* input,
                                                 std::complex<float>* tmp,
                                                 std::complex<float>* input1,
                                                 std::complex<float>* output,
                                                 const bool add,    
                                                 const float factor) const;
template void PW_Basis_K::convolution_gpu<double>(const int ik,
                                                   const int size,
                                                   const int max_npw,
                                                   const std::complex<double>* input,
                                                   std::complex<double>* tmp,
                                                   std::complex<double>* input1,
                                                   std::complex<double>* output,
                                                   const bool add,
                                                   const double factor) const;


#endif

template void PW_Basis_K::convolution_cpu<float>(const int ik,
                                                    const int size,
                                                    const int max_npw,
                                                    const std::complex<float>* input,
                                                    std::complex<float>* tmp,
                                                    std::complex<float>* input1,
                                                    std::complex<float>* output,
                                                    const bool add,
                                                    const float factor) const;
template void PW_Basis_K::convolution_cpu<double>(const int ik,
                                                   const int size,
                                                   const int max_npw,
                                                   const std::complex<double>* input,
                                                   std::complex<double>* tmp,
                                                   std::complex<double>* input1,
                                                   std::complex<double>* output,
                                                   const bool add,
                                                   const double factor) const;

template void PW_Basis_K::convolution_cpu<float>(const int ik,
                                                 const int size,
                                                 const std::complex<float>* input,
                                                 const float* input1,
                                                 std::complex<float>* output,
                                                 const bool add,
                                                 const float factor) const;
template void PW_Basis_K::convolution_cpu<double>(const int ik,
                                                  const int size,
                                                  const std::complex<double>* input,
                                                  const double* input1,
                                                  std::complex<double>* output,
                                                  const bool add,
                                                  const double factor) const;
}

