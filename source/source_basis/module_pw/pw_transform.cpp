#include "source_base/global_function.h"
#include "source_base/timer.h"
#include "source_basis/module_pw/kernels/pw_op.h"
#include "source_base/module_fft/fft_bundle.h"
#include "pw_basis.h"
#include "pw_gatherscatter.h"

#include <algorithm>
#include <cassert>
#include <complex>

namespace ModulePW
{
namespace
{
constexpr int pw_transform_cache_block = 1024;

inline int block_end(const int begin, const int size)
{
    return std::min(begin + pw_transform_cache_block, size);
}
} // namespace

//     const base_device::DEVICE_CPU* PW_Basis::get_default_device_ctx() {
//         static const base_device::DEVICE_CPU* default_device_cpu;
//     return default_device_cpu;
// }
/**
 * @brief transform real space to reciprocal space
 * @details c(g)=\int dr*f(r)*exp(-ig*r)
 *          Here we calculate c(g)
 * @param in: (nplane,ny,nx), std::complex<double> data
 * @param out: (nz, ns),  std::complex<double> data
 */
template <typename FPTYPE>
void PW_Basis::real2recip(const std::complex<FPTYPE>* in,
                          std::complex<FPTYPE>* out,
                          const bool add,
                          const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "real2recip");

    assert(this->gamma_only == false);
    const int nrxx_ = this->nrxx;
    const int npw_ = this->npw;
    const int nxyz_ = this->nxyz;
    const int* ig2isz_ = this->ig2isz;
    const std::complex<FPTYPE>* in_ = in;
    std::complex<FPTYPE>* auxr = this->fft_bundle.get_auxr_data<FPTYPE>();
    std::complex<FPTYPE>* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
    ModuleBase::timer::start(this->classname, "real2recip_copy_r");
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ib = 0; ib < nrxx_; ib += pw_transform_cache_block)
    {
        const int iend = block_end(ib, nrxx_);
#ifdef _OPENMP
#pragma omp simd
#endif
        for (int ir = ib; ir < iend; ++ir)
        {
            auxr[ir] = in_[ir];
        }
    }
    ModuleBase::timer::end(this->classname, "real2recip_copy_r");
    this->fft_bundle.fftxyfor(auxr, auxr);

    this->gatherp_scatters(auxr, auxg);

    this->fft_bundle.fftzfor(auxg, auxg);

    ModuleBase::timer::start(this->classname, "real2recip_copy_g");
    if (add)
    {
        FPTYPE tmpfac = factor / FPTYPE(nxyz_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ib = 0; ib < npw_; ib += pw_transform_cache_block)
        {
            const int iend = block_end(ib, npw_);
#ifdef _OPENMP
#pragma omp simd
#endif
            for (int ig = ib; ig < iend; ++ig)
            {
                out[ig] += tmpfac * auxg[ig2isz_[ig]];
            }
        }
    }
    else
    {
        FPTYPE tmpfac = 1.0 / FPTYPE(nxyz_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ib = 0; ib < npw_; ib += pw_transform_cache_block)
        {
            const int iend = block_end(ib, npw_);
#ifdef _OPENMP
#pragma omp simd
#endif
            for (int ig = ib; ig < iend; ++ig)
            {
                out[ig] = tmpfac * auxg[ig2isz_[ig]];
            }
        }
    }
    ModuleBase::timer::end(this->classname, "real2recip_copy_g");
    ModuleBase::timer::end(this->classname, "real2recip");
}

/**
 * @brief transform real space to reciprocal space
 * @details c(g)=\int dr*f(r)*exp(-ig*r)
 *          Here we calculate c(g)
 * @param in: (nplane,ny,nx), double data
 * @param out: (nz, ns),  std::complex<double> data
 */
template <typename FPTYPE>
void PW_Basis::real2recip(const FPTYPE* in, std::complex<FPTYPE>* out, const bool add, const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "real2recip");
    const int nrxx_ = this->nrxx;
    const int npw_ = this->npw;
    const int nxyz_ = this->nxyz;
    const int* ig2isz_ = this->ig2isz;
    const int nx_ = this->nx;
    const int ny_ = this->ny;
    const int nplane_ = this->nplane;
    const FPTYPE* in_ = in;
    std::complex<FPTYPE>* auxr = this->fft_bundle.get_auxr_data<FPTYPE>();
    std::complex<FPTYPE>* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
    FPTYPE* rspace = this->fft_bundle.get_rspace_data<FPTYPE>();
    ModuleBase::timer::start(this->classname, "real2recip_copy_r");
    if (this->gamma_only)
    {
        const int npy = ny_ * nplane_;
        const int nreal = nx_ * npy;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ib = 0; ib < nreal; ib += pw_transform_cache_block)
        {
            const int iend = block_end(ib, nreal);
#ifdef _OPENMP
#pragma omp simd
#endif
            for (int ir = ib; ir < iend; ++ir)
            {
                rspace[ir] = in_[ir];
            }
        }

        ModuleBase::timer::end(this->classname, "real2recip_copy_r");
        this->fft_bundle.fftxyr2c(rspace, auxr);
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ib = 0; ib < nrxx_; ib += pw_transform_cache_block)
        {
            const int iend = block_end(ib, nrxx_);
#ifdef _OPENMP
#pragma omp simd
#endif
            for (int ir = ib; ir < iend; ++ir)
            {
                auxr[ir] = std::complex<FPTYPE>(in_[ir], 0);
            }
        }
        ModuleBase::timer::end(this->classname, "real2recip_copy_r");
        this->fft_bundle.fftxyfor(auxr, auxr);
    }
    this->gatherp_scatters(auxr, auxg);

    this->fft_bundle.fftzfor(auxg, auxg);

    ModuleBase::timer::start(this->classname, "real2recip_copy_g");
    if (add)
    {
        FPTYPE tmpfac = factor / FPTYPE(nxyz_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ib = 0; ib < npw_; ib += pw_transform_cache_block)
        {
            const int iend = block_end(ib, npw_);
#ifdef _OPENMP
#pragma omp simd
#endif
            for (int ig = ib; ig < iend; ++ig)
            {
                out[ig] += tmpfac * auxg[ig2isz_[ig]];
            }
        }
    }
    else
    {
        FPTYPE tmpfac = 1.0 / FPTYPE(nxyz_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ib = 0; ib < npw_; ib += pw_transform_cache_block)
        {
            const int iend = block_end(ib, npw_);
#ifdef _OPENMP
#pragma omp simd
#endif
            for (int ig = ib; ig < iend; ++ig)
            {
                out[ig] = tmpfac * auxg[ig2isz_[ig]];
            }
        }
    }
    ModuleBase::timer::end(this->classname, "real2recip_copy_g");
    ModuleBase::timer::end(this->classname, "real2recip");
}

/**
 * @brief transform reciprocal space to real space
 * @details f(r)=1/V * \sum_{g} c(g)*exp(ig*r)
 *          Here we calculate f(r)
 * @param in: (nz,ns), std::complex<double>
 * @param out: (nplane, ny, nx), std::complex<double>
 */
template <typename FPTYPE>
void PW_Basis::recip2real(const std::complex<FPTYPE>* in,
                          std::complex<FPTYPE>* out,
                          const bool add,
                          const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "recip2real");
    assert(this->gamma_only == false);
    const int nst_ = this->nst;
    const int nz_ = this->nz;
    const int npw_ = this->npw;
    const int nrxx_ = this->nrxx;
    const int* ig2isz_ = this->ig2isz;
    const int nstnz_ = nst_ * nz_;
    std::complex<FPTYPE>* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
    std::complex<FPTYPE>* auxr = this->fft_bundle.get_auxr_data<FPTYPE>();
    ModuleBase::timer::start(this->classname, "recip2real_copy_g");
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ib = 0; ib < nstnz_; ib += pw_transform_cache_block)
    {
        const int iend = block_end(ib, nstnz_);
#ifdef _OPENMP
#pragma omp simd
#endif
        for (int i = ib; i < iend; ++i)
        {
            auxg[i] = std::complex<FPTYPE>(0, 0);
        }
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ib = 0; ib < npw_; ib += pw_transform_cache_block)
    {
        const int iend = block_end(ib, npw_);
#ifdef _OPENMP
#pragma omp simd
#endif
        for (int ig = ib; ig < iend; ++ig)
        {
            auxg[ig2isz_[ig]] = in[ig];
        }
    }
    ModuleBase::timer::end(this->classname, "recip2real_copy_g");
    this->fft_bundle.fftzbac(auxg, auxg);

    this->gathers_scatterp(auxg, auxr);

    this->fft_bundle.fftxybac(auxr, auxr);

    ModuleBase::timer::start(this->classname, "recip2real_copy_r");
    if (add)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ib = 0; ib < nrxx_; ib += pw_transform_cache_block)
        {
            const int iend = block_end(ib, nrxx_);
#ifdef _OPENMP
#pragma omp simd
#endif
            for (int ir = ib; ir < iend; ++ir)
            {
                out[ir] += factor * auxr[ir];
            }
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ib = 0; ib < nrxx_; ib += pw_transform_cache_block)
        {
            const int iend = block_end(ib, nrxx_);
#ifdef _OPENMP
#pragma omp simd
#endif
            for (int ir = ib; ir < iend; ++ir)
            {
                out[ir] = auxr[ir];
            }
        }
    }
    ModuleBase::timer::end(this->classname, "recip2real_copy_r");
    ModuleBase::timer::end(this->classname, "recip2real");
}

/**
 * @brief transform reciprocal space to real space
 * @details f(r)=1/V * \sum_{g} c(g)*exp(ig*r)
 *          Here we calculate f(r)
 * @param in: (nz,ns), std::complex<double>
 * @param out: (nplane, ny, nx), double
 */
template <typename FPTYPE>
void PW_Basis::recip2real(const std::complex<FPTYPE>* in, FPTYPE* out, const bool add, const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "recip2real");
    const int nst_ = this->nst;
    const int nz_ = this->nz;
    const int npw_ = this->npw;
    const int nrxx_ = this->nrxx;
    const int nx_ = this->nx;
    const int ny_ = this->ny;
    const int nplane_ = this->nplane;
    const int* ig2isz_ = this->ig2isz;
    const int nstnz_ = nst_ * nz_;
    std::complex<FPTYPE>* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
    std::complex<FPTYPE>* auxr = this->fft_bundle.get_auxr_data<FPTYPE>();
    FPTYPE* rspace = this->fft_bundle.get_rspace_data<FPTYPE>();
    ModuleBase::timer::start(this->classname, "recip2real_copy_g");
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ib = 0; ib < nstnz_; ib += pw_transform_cache_block)
    {
        const int iend = block_end(ib, nstnz_);
#ifdef _OPENMP
#pragma omp simd
#endif
        for (int i = ib; i < iend; ++i)
        {
            auxg[i] = std::complex<FPTYPE>(0, 0);
        }
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ib = 0; ib < npw_; ib += pw_transform_cache_block)
    {
        const int iend = block_end(ib, npw_);
#ifdef _OPENMP
#pragma omp simd
#endif
        for (int ig = ib; ig < iend; ++ig)
        {
            auxg[ig2isz_[ig]] = in[ig];
        }
    }
    ModuleBase::timer::end(this->classname, "recip2real_copy_g");
    this->fft_bundle.fftzbac(auxg, auxg);

    this->gathers_scatterp(auxg, auxr);

    if (this->gamma_only)
    {
        this->fft_bundle.fftxyc2r(auxr, rspace);

        const int npy = ny_ * nplane_;
        const int nreal = nx_ * npy;

        ModuleBase::timer::start(this->classname, "recip2real_copy_r");
        if (add)
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int ib = 0; ib < nreal; ib += pw_transform_cache_block)
            {
                const int iend = block_end(ib, nreal);
#ifdef _OPENMP
#pragma omp simd
#endif
                for (int ir = ib; ir < iend; ++ir)
                {
                    out[ir] += factor * rspace[ir];
                }
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int ib = 0; ib < nreal; ib += pw_transform_cache_block)
            {
                const int iend = block_end(ib, nreal);
#ifdef _OPENMP
#pragma omp simd
#endif
                for (int ir = ib; ir < iend; ++ir)
                {
                    out[ir] = rspace[ir];
                }
            }
        }
        ModuleBase::timer::end(this->classname, "recip2real_copy_r");
    }
    else
    {
        this->fft_bundle.fftxybac(auxr, auxr);
        ModuleBase::timer::start(this->classname, "recip2real_copy_r");
        if (add)
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int ib = 0; ib < nrxx_; ib += pw_transform_cache_block)
            {
                const int iend = block_end(ib, nrxx_);
#ifdef _OPENMP
#pragma omp simd
#endif
                for (int ir = ib; ir < iend; ++ir)
                {
                    out[ir] += factor * auxr[ir].real();
                }
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int ib = 0; ib < nrxx_; ib += pw_transform_cache_block)
            {
                const int iend = block_end(ib, nrxx_);
#ifdef _OPENMP
#pragma omp simd
#endif
                for (int ir = ib; ir < iend; ++ir)
                {
                    out[ir] = auxr[ir].real();
                }
            }
        }
        ModuleBase::timer::end(this->classname, "recip2real_copy_r");
    }
    ModuleBase::timer::end(this->classname, "recip2real");
}
template void PW_Basis::real2recip<float>(const float* in,
                                          std::complex<float>* out,
                                          const bool add,
                                          const float factor) const; // in:(nplane,nx*ny)  ; out(nz, ns)
template void PW_Basis::real2recip<float>(const std::complex<float>* in,
                                          std::complex<float>* out,
                                          const bool add,
                                          const float factor) const; // in:(nplane,nx*ny)  ; out(nz, ns)
template void PW_Basis::recip2real<float>(const std::complex<float>* in,
                                          float* out,
                                          const bool add,
                                          const float factor) const; // in:(nz, ns)  ; out(nplane,nx*ny)
template void PW_Basis::recip2real<float>(const std::complex<float>* in,
                                          std::complex<float>* out,
                                          const bool add,
                                          const float factor) const;

template void PW_Basis::real2recip<double>(const double* in,
                                           std::complex<double>* out,
                                           const bool add,
                                           const double factor) const; // in:(nplane,nx*ny)  ; out(nz, ns)
template void PW_Basis::real2recip<double>(const std::complex<double>* in,
                                           std::complex<double>* out,
                                           const bool add,
                                           const double factor) const; // in:(nplane,nx*ny)  ; out(nz, ns)
template void PW_Basis::recip2real<double>(const std::complex<double>* in,
                                           double* out,
                                           const bool add,
                                           const double factor) const; // in:(nz, ns)  ; out(nplane,nx*ny)
template void PW_Basis::recip2real<double>(const std::complex<double>* in,
                                           std::complex<double>* out,
                                           const bool add,
                                           const double factor) const;
} // namespace ModulePW
