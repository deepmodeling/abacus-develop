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
constexpr int pw_transform_cache_block = 128;

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
 * @brief Transform real-space data to reciprocal-space plane-wave coefficients (complex input).
 * @details
 * Performs the forward 3D FFT with MPI parallel transposition for non-gamma-only calculations.
 * Computes: c(g) = (1/N) * sum_r f(r) * exp(-i g·r)
 *
 * This is the #1 hotspot function in ABACUS, accounting for 22-30% of total SCF runtime.
 * The implementation uses a 2D domain decomposition strategy:
 * 1. Copy input real-space data to FFT buffer (z-slab distributed, O(nrxx))
 * 2. In-place 2D FFT on each xy-plane (fftxyfor, independent per process)
 * 3. MPI_Alltoallv transposition: xy-planes → z-sticks (gatherp_scatters)
 *    - Communication volume: O(nst_per * nz * sizeof(complex))
 * 4. In-place 1D FFT along each z-stick (fftzfor)
 * 5. Extract plane-wave coefficients: out[ig] = auxg[ig2isz[ig]] / nxyz
 *
 * @tparam FPTYPE Floating-point precision (float or double)
 * @param in  Input real-space array, shape (nplane, ny, nx) in z-slab distribution
 *            Each MPI process holds nplane xy-planes, for nrxx = nplane*nx*ny local elements
 * @param out Output reciprocal-space array, shape (npw) — plane-wave coefficients
 *            Only stores coefficients for G-vectors on this process (stick distribution)
 * @param add If true, add scaled result to existing out[]; if false, overwrite out[]
 * @param factor Scaling factor applied when add=true: out[ig] += factor * c(g)
 * @note The 1/nxyz normalization is always applied regardless of add/factor
 * @note For gamma-only calculations, use the real-input overload (r2c FFT path)
 * @see recip2real() for the inverse transform, gatherp_scatters() for MPI communication
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
 * @brief Transform real-valued real-space data to reciprocal-space (gamma-only or non-gamma).
 * @details
 * Two code paths depending on gamma_only flag:
 * - gamma_only=true:  Uses r2c FFT (fftxyr2c). Only half the FFT grid is stored (fftnx = nx/2+1),
 *   exploiting Hermitian symmetry to save ~50% memory and computation.
 * - gamma_only=false: Converts real input to complex, then follows the same 3D FFT path as
 *   the complex-input overload (fftxyfor → gatherp_scatters → fftzfor).
 *
 * @tparam FPTYPE Floating-point precision (float or double)
 * @param in  Input real-space array (real-valued), shape (nplane, ny, nx) in z-slab distribution
 * @param out Output reciprocal-space plane-wave coefficients (complex)
 * @param add  If true, accumulate scaled result into out[]; if false, overwrite
 * @param factor Scaling factor for add mode
 * @see real2recip(const std::complex<FPTYPE>*, ...) for the complex-input variant
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
 * @brief Transform reciprocal-space plane-wave coefficients to real-space data (complex output).
 * @details
 * Performs the inverse 3D FFT — the reverse of real2recip():
 * f(r) = sum_g c(g) * exp(i g·r)
 *
 * Algorithm (reverse of real2recip):
 * 1. Zero-fill the FFT stick buffer (nst*nz elements), then scatter: auxg[ig2isz[ig]] = in[ig]
 * 2. Backward 1D FFT along each z-stick (fftzbac)
 * 3. MPI_Alltoallv transposition: sticks → xy-planes (gathers_scatterp, reverse direction)
 * 4. Backward 2D FFT on each xy-plane (fftxybac)
 * 5. Copy/extract real-space result: out[ir] = auxr[ir]
 *
 * @tparam FPTYPE Floating-point precision (float or double)
 * @param in  Input reciprocal-space array, shape (npw) — plane-wave coefficients in stick distribution
 * @param out Output real-space array, shape (nplane, ny, nx) in z-slab distribution
 * @param add If true, add scaled result to existing out[]; if false, overwrite
 * @param factor Scaling factor for add mode: out[ir] += factor * f(r)
 * @note No 1/nxyz normalization factor is applied (unlike real2recip)
 * @see real2recip() for the forward transform, gathers_scatterp() for MPI communication
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
 * @brief Transform reciprocal-space to real-valued real-space (gamma-only or non-gamma).
 * @details
 * Two code paths:
 * - gamma_only=true:  Uses c2r FFT (fftxyc2r) to exploit Hermitian symmetry. After backward 1D FFT
 *   and MPI transposition, applies c2r FFT producing real-valued output directly.
 * - gamma_only=false: Follows the standard complex path (fftzbac → gathers_scatterp → fftxybac),
 *   then extracts the real part of the complex result.
 *
 * @tparam FPTYPE Floating-point precision (float or double)
 * @param in  Input reciprocal-space plane-wave coefficients (complex)
 * @param out Output real-space array (real-valued)
 * @param add If true, accumulate scaled result into out[]; if false, overwrite
 * @param factor Scaling factor for add mode
 * @see recip2real(const std::complex<FPTYPE>*, std::complex<FPTYPE>*, ...) for complex output
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
