#include "source_base/timer.h"
#include "source_basis/module_pw/kernels/pw_op.h"
#include "pw_basis_k.h"
#include "pw_gatherscatter.h"
#include "source_pw/module_pwdft/kernels/veff_op.h"
#include <cassert>
#include <complex>

namespace ModulePW
{

/**
 * @brief transform real space to reciprocal space
 * @details real wave function f(k,r):
 *          f(k,r)=1/V*\sum_{g} c(k,g)*exp(i(g+k)*r) \equiv exp(ikr)f'(k.r)
 *          c(k,g)=\int dr*f(k,r)*exp(-i(g+k)*r)
 *          However, we use f'(k,r)!!! :
 *          f'(k,r)=1/V*\sum_{g} c(k,g)*exp(ig*r)
 *          c(k,g)=\int dr*f'(k,r)*exp(-ig*r)
 *
 *          This function tranform f'(r) to c(k,g).
 * @param in: (nplane,ny,nx), std::complex<double> data
 * @param out: (nz, ns),  std::complex<double> data
 */
template <typename FPTYPE>
void PW_Basis_K::real2recip(const std::complex<FPTYPE>* in,
                            std::complex<FPTYPE>* out,
                            const int ik,
                            const bool add,
                            const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "real2recip");

    assert(this->gamma_only == false);
    auto* auxr = this->fft_bundle.get_auxr_data<FPTYPE>();
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
    for (int ir = 0; ir < this->nrxx; ++ir)
    {
        auxr[ir] = in[ir];
    }
    this->fft_bundle.fftxyfor(fft_bundle.get_auxr_data<FPTYPE>(), fft_bundle.get_auxr_data<FPTYPE>());

    this->gatherp_scatters(this->fft_bundle.get_auxr_data<FPTYPE>(), this->fft_bundle.get_auxg_data<FPTYPE>());

    this->fft_bundle.fftzfor(fft_bundle.get_auxg_data<FPTYPE>(), fft_bundle.get_auxg_data<FPTYPE>());

    const int startig = ik * this->npwk_max;
    const int npwk = this->npwk[ik];
    auto* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
    if (add)
    {
        FPTYPE tmpfac = factor / FPTYPE(this->nxyz);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
        for (int igl = 0; igl < npwk; ++igl)
        {
            out[igl] += tmpfac * auxg[this->igl2isz_k[igl + startig]];
        }
    }
    else
    {
        FPTYPE tmpfac = 1.0 / FPTYPE(this->nxyz);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
        for (int igl = 0; igl < npwk; ++igl)
        {
            out[igl] = tmpfac * auxg[this->igl2isz_k[igl + startig]];
        }
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
}

/**
 * @brief transform real space to reciprocal space
 * @details real wave function f(k,r):
 *          f(k,r)=1/V*\sum_{g} c(k,g)*exp(i(g+k)*r) \equiv exp(ikr)f'(k.r)
 *          c(k,g)=\int dr*f(k,r)*exp(-i(g+k)*r)
 *          However, we use f'(k,r)!!! :
 *          f'(k,r)=1/V*\sum_{g} c(k,g)*exp(ig*r)
 *          c(k,g)=\int dr*f'(k,r)*exp(-ig*r)
 *
 *          This function tranform f'(r) to c(k,g).
 * @param in: (nplane,ny,nx), double data
 * @param out: (nz, ns),  std::complex<double> data
 */
template <typename FPTYPE>
void PW_Basis_K::real2recip(const FPTYPE* in,
                            std::complex<FPTYPE>* out,
                            const int ik,
                            const bool add,
                            const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "real2recip");
    assert(this->gamma_only == true);
    // for(int ir = 0 ; ir < this->nrxx ; ++ir)
    // {
    //     this->fft_bundle.get_rspace_data<FPTYPE>()[ir] = in[ir];
    // }
    // r2c in place
    const int npy = this->ny * this->nplane;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 4096 / sizeof(FPTYPE))
#endif
    for (int ix = 0; ix < this->nx; ++ix)
    {
        for (int ipy = 0; ipy < npy; ++ipy)
        {
            this->fft_bundle.get_rspace_data<FPTYPE>()[ix * npy + ipy] = in[ix * npy + ipy];
        }
    }

    this->fft_bundle.fftxyr2c(fft_bundle.get_rspace_data<FPTYPE>(), fft_bundle.get_auxr_data<FPTYPE>());

    this->gatherp_scatters(this->fft_bundle.get_auxr_data<FPTYPE>(), this->fft_bundle.get_auxg_data<FPTYPE>());

    this->fft_bundle.fftzfor(fft_bundle.get_auxg_data<FPTYPE>(), fft_bundle.get_auxg_data<FPTYPE>());

    const int startig = ik * this->npwk_max;
    const int npwk = this->npwk[ik];
    auto* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
    if (add)
    {
        FPTYPE tmpfac = factor / FPTYPE(this->nxyz);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
        for (int igl = 0; igl < npwk; ++igl)
        {
            out[igl] += tmpfac * auxg[this->igl2isz_k[igl + startig]];
        }
    }
    else
    {
        FPTYPE tmpfac = 1.0 / FPTYPE(this->nxyz);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
        for (int igl = 0; igl < npwk; ++igl)
        {
            out[igl] = tmpfac * auxg[this->igl2isz_k[igl + startig]];
        }
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
    return;
}

/**
 * @brief transform reciprocal space to real space
 * @details real wave function f(k,r):
 *          f(k,r)=1/V*\sum_{g} c(k,g)*exp(i(g+k)*r) \equiv exp(ikr)f'(k.r)
 *          c(k,g)=\int dr*f(k,r)*exp(-i(g+k)*r)
 *          However, we use f'(k,r)!!! :
 *          f'(k,r)=1/V*\sum_{g} c(k,g)*exp(ig*r)
 *          c(k,g)=\int dr*f'(k,r)*exp(-ig*r)
 *
 *          This function tranform c(k,g) to f'(r).
 * @param in: (nz,ns), std::complex<double>
 * @param out: (nplane, ny, nx), std::complex<double>
 */
template <typename FPTYPE>
void PW_Basis_K::recip2real(const std::complex<FPTYPE>* in,
                            std::complex<FPTYPE>* out,
                            const int ik,
                            const bool add,
                            const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "recip2real");
    assert(this->gamma_only == false);
    ModuleBase::GlobalFunc::ZEROS(fft_bundle.get_auxg_data<FPTYPE>(), this->nst * this->nz);

    const int startig = ik * this->npwk_max;
    const int npwk = this->npwk[ik];
    auto* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
    for (int igl = 0; igl < npwk; ++igl)
    {
        auxg[this->igl2isz_k[igl + startig]] = in[igl];
    }
    this->fft_bundle.fftzbac(fft_bundle.get_auxg_data<FPTYPE>(), fft_bundle.get_auxg_data<FPTYPE>());

    this->gathers_scatterp(this->fft_bundle.get_auxg_data<FPTYPE>(), this->fft_bundle.get_auxr_data<FPTYPE>());

    this->fft_bundle.fftxybac(fft_bundle.get_auxr_data<FPTYPE>(), fft_bundle.get_auxr_data<FPTYPE>());
    auto* auxr = this->fft_bundle.get_auxr_data<FPTYPE>();
    if (add)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            out[ir] += factor * auxr[ir];
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            out[ir] = auxr[ir];
        }
    }
    ModuleBase::timer::tick(this->classname, "recip2real");
}

/**
 * @brief transform reciprocal space to real space
 * @details real wave function f(k,r):
 *          f(k,r)=1/V*\sum_{g} c(k,g)*exp(i(g+k)*r) \equiv exp(ikr)f'(k.r)
 *          c(k,g)=\int dr*f(k,r)*exp(-i(g+k)*r)
 *          However, we use f'(k,r)!!! :
 *          f'(k,r)=1/V*\sum_{g} c(k,g)*exp(ig*r)
 *          c(k,g)=\int dr*f'(k,r)*exp(-ig*r)
 *
 *          This function tranform c(k,g) to f'(r).
 * @param in: (nz,ns), std::complex<double>
 * @param out: (nplane, ny, nx), double
 */
template <typename FPTYPE>
void PW_Basis_K::recip2real(const std::complex<FPTYPE>* in,
                            FPTYPE* out,
                            const int ik,
                            const bool add,
                            const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "recip2real");
    assert(this->gamma_only == true);
    ModuleBase::GlobalFunc::ZEROS(fft_bundle.get_auxg_data<FPTYPE>(), this->nst * this->nz);

    const int startig = ik * this->npwk_max;
    const int npwk = this->npwk[ik];
    auto* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
    for (int igl = 0; igl < npwk; ++igl)
    {
        auxg[this->igl2isz_k[igl + startig]] = in[igl];
    }
    this->fft_bundle.fftzbac(fft_bundle.get_auxg_data<FPTYPE>(), fft_bundle.get_auxg_data<FPTYPE>());

    this->gathers_scatterp(this->fft_bundle.get_auxg_data<FPTYPE>(), this->fft_bundle.get_auxr_data<FPTYPE>());

    this->fft_bundle.fftxyc2r(fft_bundle.get_auxr_data<FPTYPE>(), fft_bundle.get_rspace_data<FPTYPE>());

    // for(int ir = 0 ; ir < this->nrxx ; ++ir)
    // {
    //     out[ir] = this->fft_bundle.get_rspace_data<FPTYPE>()[ir] / this->nxyz;
    // }

    // r2c in place
    const int npy = this->ny * this->nplane;
    auto* rspace = this->fft_bundle.get_rspace_data<FPTYPE>();
    if (add)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 4096 / sizeof(FPTYPE))
#endif
        for (int ix = 0; ix < this->nx; ++ix)
        {
            for (int ipy = 0; ipy < npy; ++ipy)
            {
                out[ix * npy + ipy] += factor * rspace[ix * npy + ipy];
            }
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 4096 / sizeof(FPTYPE))
#endif
        for (int ix = 0; ix < this->nx; ++ix)
        {
            for (int ipy = 0; ipy < npy; ++ipy)
            {
                out[ix * npy + ipy] = rspace[ix * npy + ipy];
            }
        }
    }
    ModuleBase::timer::tick(this->classname, "recip2real");
}

template <>
void PW_Basis_K::real_to_recip(const base_device::DEVICE_CPU* /*dev*/,
                               const std::complex<float>* in,
                               std::complex<float>* out,
                               const int ik,
                               const bool add,
                               const float factor) const
{
    this->real2recip(in, out, ik, add, factor);
}
template <>
void PW_Basis_K::real_to_recip(const base_device::DEVICE_CPU* /*dev*/,
                               const std::complex<double>* in,
                               std::complex<double>* out,
                               const int ik,
                               const bool add,
                               const double factor) const
{
    #if defined(__DSP)
        this->real2recip_dsp(in,out,ik,add,factor);
    #else
        this->real2recip(in, out, ik, add, factor);
    #endif
}

template <>
void PW_Basis_K::recip_to_real(const base_device::DEVICE_CPU* /*dev*/,
                               const std::complex<float>* in,
                               std::complex<float>* out,
                               const int ik,
                               const bool add,
                               const float factor) const
{
    this->recip2real(in, out, ik, add, factor);
}
template <>
void PW_Basis_K::recip_to_real(const base_device::DEVICE_CPU* /*dev*/,
                               const std::complex<double>* in,
                               std::complex<double>* out,
                               const int ik,
                               const bool add,
                               const double factor) const
{
    #if defined(__DSP)
        this->recip2real_dsp(in,out,ik,add,factor);
    #else
        this->recip2real(in, out, ik, add, factor);
    #endif
}
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
        // const int idx = this->igl2isz_k[igl + startig];
        augr[this->igl2isz_k[igl + startig]] = input[igl];
    }
    for (int igl =0 ; igl < npwk ; ++igl)
    {
        augr1[this->igl2isz_k[igl + startig]] = input[igl+max_npw];
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
            // const int idx = this->igl2isz_k[igl + startig];
            output[igl] += tmpfac * augr[this->igl2isz_k[igl + startig]];
        }
    
        for (int igl =0 ; igl < npwk ; ++igl)
        {
            output[igl+max_npw] += tmpfac * augr1[this->igl2isz_k[igl + startig]];
        }
    ModuleBase::timer::tick(this->classname, "convolution");
    }


#if (defined(__CUDA) || defined(__ROCM))
template <>
void PW_Basis_K::real_to_recip(const base_device::DEVICE_GPU* ctx,
                               const std::complex<float>* in,
                               std::complex<float>* out,
                               const int ik,
                               const bool add,
                               const float factor) const
{
    ModuleBase::timer::tick(this->classname, "real_to_recip gpu");
    assert(this->gamma_only == false);
    assert(this->poolnproc == 1);

    base_device::memory::synchronize_memory_op<std::complex<float>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(
        this->fft_bundle.get_auxr_3d_data<float>(),
        in,
        this->nrxx);

    this->fft_bundle.fft3D_forward(this->fft_bundle.get_auxr_3d_data<float>(), this->fft_bundle.get_auxr_3d_data<float>());

    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];
    set_real_to_recip_output_op<float, base_device::DEVICE_GPU>()(npw_k,
                                                                  this->nxyz,
                                                                  add,
                                                                  factor,
                                                                  this->ig2ixyz_k + startig,
                                                                  this->fft_bundle.get_auxr_3d_data<float>(),
                                                                  out);
    ModuleBase::timer::tick(this->classname, "real_to_recip gpu");
}
template <>
void PW_Basis_K::real_to_recip(const base_device::DEVICE_GPU* ctx,
                               const std::complex<double>* in,
                               std::complex<double>* out,
                               const int ik,
                               const bool add,
                               const double factor) const
{
    ModuleBase::timer::tick(this->classname, "real_to_recip gpu");
    assert(this->gamma_only == false);
    assert(this->poolnproc == 1);

    base_device::memory::synchronize_memory_op<std::complex<double>,
                                               base_device::DEVICE_GPU,
                                               base_device::DEVICE_GPU>()(this->fft_bundle.get_auxr_3d_data<double>(),
                                                                          in,
                                                                          this->nrxx);

    this->fft_bundle.fft3D_forward(this->fft_bundle.get_auxr_3d_data<double>(), this->fft_bundle.get_auxr_3d_data<double>());

    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];
    set_real_to_recip_output_op<double, base_device::DEVICE_GPU>()(npw_k,
                                                                   this->nxyz,
                                                                   add,
                                                                   factor,
                                                                   this->ig2ixyz_k + startig,
                                                                   this->fft_bundle.get_auxr_3d_data<double>(),
                                                                   out);
    ModuleBase::timer::tick(this->classname, "real_to_recip gpu");
}

template <>
void PW_Basis_K::recip_to_real(const base_device::DEVICE_GPU* ctx,
                               const std::complex<float>* in,
                               std::complex<float>* out,
                               const int ik,
                               const bool add,
                               const float factor) const
{
    ModuleBase::timer::tick(this->classname, "recip_to_real gpu");
    assert(this->gamma_only == false);
    assert(this->poolnproc == 1);
    // ModuleBase::GlobalFunc::ZEROS(fft_bundle.get_auxr_3d_data<float>(), this->nxyz);
    base_device::memory::set_memory_op<std::complex<float>, base_device::DEVICE_GPU>()(
        this->fft_bundle.get_auxr_3d_data<float>(),
        0,
        this->nxyz);

    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];

    set_3d_fft_box_op<float, base_device::DEVICE_GPU>()(npw_k,
                                                        this->ig2ixyz_k + startig,
                                                        in,
                                                        this->fft_bundle.get_auxr_3d_data<float>());
    this->fft_bundle.fft3D_backward(this->fft_bundle.get_auxr_3d_data<float>(), this->fft_bundle.get_auxr_3d_data<float>());

    set_recip_to_real_output_op<float, base_device::DEVICE_GPU>()(this->nrxx,
                                                                  add,
                                                                  factor,
                                                                  this->fft_bundle.get_auxr_3d_data<float>(),
                                                                  out);

    ModuleBase::timer::tick(this->classname, "recip_to_real gpu");
}
template <>
void PW_Basis_K::recip_to_real(const base_device::DEVICE_GPU* ctx,
                               const std::complex<double>* in,
                               std::complex<double>* out,
                               const int ik,
                               const bool add,
                               const double factor) const
{
    ModuleBase::timer::tick(this->classname, "recip_to_real gpu");
    assert(this->gamma_only == false);
    assert(this->poolnproc == 1);
    // ModuleBase::GlobalFunc::ZEROS(fft_bundle.get_auxr_3d_data<double>(), this->nxyz);
    base_device::memory::set_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(
        this->fft_bundle.get_auxr_3d_data<double>(),
        0,
        this->nxyz);

    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];

    set_3d_fft_box_op<double, base_device::DEVICE_GPU>()(npw_k,
                                                         this->ig2ixyz_k + startig,
                                                         in,
                                                         this->fft_bundle.get_auxr_3d_data<double>());
    this->fft_bundle.fft3D_backward(this->fft_bundle.get_auxr_3d_data<double>(), this->fft_bundle.get_auxr_3d_data<double>());

    set_recip_to_real_output_op<double, base_device::DEVICE_GPU>()(this->nrxx,
                                                                   add,
                                                                   factor,
                                                                   this->fft_bundle.get_auxr_3d_data<double>(),
                                                                   out);

    ModuleBase::timer::tick(this->classname, "recip_to_real gpu");
}

template <typename FPTYPE>
void PW_Basis_K::real2recip_gpu(const std::complex<FPTYPE>* in,
                               std::complex<FPTYPE>* out,
                               const int ik,
                               const bool add,
                               const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "real_to_recip gpu");
    assert(this->gamma_only == false);
    assert(this->poolnproc == 1);

    base_device::memory::synchronize_memory_op<std::complex<FPTYPE>,
                                               base_device::DEVICE_GPU,
                                               base_device::DEVICE_GPU>()(this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                                                          in,
                                                                          this->nrxx);

    this->fft_bundle.fft3D_forward(this->fft_bundle.get_auxr_3d_data<FPTYPE>(), this->fft_bundle.get_auxr_3d_data<FPTYPE>());

    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];
    set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>()(npw_k,
                                                                   this->nxyz,
                                                                   add,
                                                                   factor,
                                                                   this->ig2ixyz_k + startig,
                                                                   this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                                                   out);
    ModuleBase::timer::tick(this->classname, "real_to_recip gpu");
}
template <typename FPTYPE>
void PW_Basis_K::recip2real_gpu(const std::complex<FPTYPE>* in,
                               std::complex<FPTYPE>* out,
                               const int ik,
                               const bool add,
                               const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "recip_to_real gpu");
    assert(this->gamma_only == false);
    assert(this->poolnproc == 1);
    // ModuleBase::GlobalFunc::ZEROS(fft_bundle.get_auxr_3d_data<FPTYPE>(), this->nxyz);
    base_device::memory::set_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(
        this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
        0,
        this->nxyz);

    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];

    set_3d_fft_box_op<FPTYPE, base_device::DEVICE_GPU>()(npw_k,
                                                         this->ig2ixyz_k + startig,
                                                         in,
                                                         this->fft_bundle.get_auxr_3d_data<FPTYPE>());
    this->fft_bundle.fft3D_backward(this->fft_bundle.get_auxr_3d_data<FPTYPE>(), this->fft_bundle.get_auxr_3d_data<FPTYPE>());

    set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>()(this->nrxx,
                                                                   add,
                                                                   factor,
                                                                   this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                                                   out);

    ModuleBase::timer::tick(this->classname, "recip_to_real gpu");
}

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

template void PW_Basis_K::real2recip_gpu<float>(const std::complex<float>*,
                                                std::complex<float>*,
                                                const int,
                                                const bool,
                                                const float) const;

template void PW_Basis_K::real2recip_gpu<double>(const std::complex<double>*,
                                                 std::complex<double>*,
                                                 const int,
                                                 const bool,
                                                 const double) const;

template void PW_Basis_K::recip2real_gpu<float>(const std::complex<float>*,
                                                std::complex<float>*,
                                                const int,
                                                const bool,
                                                const float) const;

template void PW_Basis_K::recip2real_gpu<double>(const std::complex<double>*,
                                                 std::complex<double>*,
                                                 const int,
                                                 const bool,
                                                 const double) const;
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

template void PW_Basis_K::real2recip<float>(const float* in,
                                            std::complex<float>* out,
                                            const int ik,
                                            const bool add,
                                            const float factor) const; // in:(nplane,nx*ny)  ; out(nz, ns)
template void PW_Basis_K::real2recip<float>(const std::complex<float>* in,
                                            std::complex<float>* out,
                                            const int ik,
                                            const bool add,
                                            const float factor) const; // in:(nplane,nx*ny)  ; out(nz, ns)
template void PW_Basis_K::recip2real<float>(const std::complex<float>* in,
                                            float* out,
                                            const int ik,
                                            const bool add,
                                            const float factor) const; // in:(nz, ns)  ; out(nplane,nx*ny)
template void PW_Basis_K::recip2real<float>(const std::complex<float>* in,
                                            std::complex<float>* out,
                                            const int ik,
                                            const bool add,
                                            const float factor) const; // in:(nz, ns)  ; out(nplane,nx*ny)

template void PW_Basis_K::real2recip<double>(const double* in,
                                             std::complex<double>* out,
                                             const int ik,
                                             const bool add,
                                             const double factor) const; // in:(nplane,nx*ny)  ; out(nz, ns)
template void PW_Basis_K::real2recip<double>(const std::complex<double>* in,
                                             std::complex<double>* out,
                                             const int ik,
                                             const bool add,
                                             const double factor) const; // in:(nplane,nx*ny)  ; out(nz, ns)
template void PW_Basis_K::recip2real<double>(const std::complex<double>* in,
                                             double* out,
                                             const int ik,
                                             const bool add,
                                             const double factor) const; // in:(nz, ns)  ; out(nplane,nx*ny)
template void PW_Basis_K::recip2real<double>(const std::complex<double>* in,
                                             std::complex<double>* out,
                                             const int ik,
                                             const bool add,
                                             const double factor) const; // in:(nz, ns)  ; out(nplane,nx*ny)
} // namespace ModulePW
