#include "module_base/timer.h"
#include "module_basis/module_pw/kernels/pw_op.h"
#include "pw_basis_k.h"
#include "pw_gatherscatter.h"

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
 * @param in: (nplane,ny,nx), complex<double> data
 * @param out: (nz, ns),  complex<double> data
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
 * @param out: (nz, ns),  complex<double> data
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
 * @param in: (nz,ns), complex<double>
 * @param out: (nplane, ny, nx), complex<double>
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
 * @param in: (nz,ns), complex<double>
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
    this->real2recip(in, out, ik, add, factor);
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
    this->recip2real(in, out, ik, add, factor);
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
        ctx,
        ctx,
        this->fft_bundle.get_auxr_3d_data<float>(),
        in,
        this->nrxx);

    this->fft_bundle.fft3D_forward(ctx, this->fft_bundle.get_auxr_3d_data<float>(), this->fft_bundle.get_auxr_3d_data<float>());

    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];
    set_real_to_recip_output_op<float, base_device::DEVICE_GPU>()(ctx,
                                                                  npw_k,
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
                                               base_device::DEVICE_GPU>()(ctx,
                                                                          ctx,
                                                                          this->fft_bundle.get_auxr_3d_data<double>(),
                                                                          in,
                                                                          this->nrxx);

    this->fft_bundle.fft3D_forward(ctx, this->fft_bundle.get_auxr_3d_data<double>(), this->fft_bundle.get_auxr_3d_data<double>());

    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];
    set_real_to_recip_output_op<double, base_device::DEVICE_GPU>()(ctx,
                                                                   npw_k,
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
        ctx,
        this->fft_bundle.get_auxr_3d_data<float>(),
        0,
        this->nxyz);

    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];

    set_3d_fft_box_op<float, base_device::DEVICE_GPU>()(ctx,
                                                        npw_k,
                                                        this->ig2ixyz_k + startig,
                                                        in,
                                                        this->fft_bundle.get_auxr_3d_data<float>());
    this->fft_bundle.fft3D_backward(ctx, this->fft_bundle.get_auxr_3d_data<float>(), this->fft_bundle.get_auxr_3d_data<float>());

    set_recip_to_real_output_op<float, base_device::DEVICE_GPU>()(ctx,
                                                                  this->nrxx,
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
        ctx,
        this->fft_bundle.get_auxr_3d_data<double>(),
        0,
        this->nxyz);

    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];

    set_3d_fft_box_op<double, base_device::DEVICE_GPU>()(ctx,
                                                         npw_k,
                                                         this->ig2ixyz_k + startig,
                                                         in,
                                                         this->fft_bundle.get_auxr_3d_data<double>());
    this->fft_bundle.fft3D_backward(ctx, this->fft_bundle.get_auxr_3d_data<double>(), this->fft_bundle.get_auxr_3d_data<double>());

    set_recip_to_real_output_op<double, base_device::DEVICE_GPU>()(ctx,
                                                                   this->nrxx,
                                                                   add,
                                                                   factor,
                                                                   this->fft_bundle.get_auxr_3d_data<double>(),
                                                                   out);

    ModuleBase::timer::tick(this->classname, "recip_to_real gpu");
}
#endif

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
