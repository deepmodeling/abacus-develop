#include "source_base/global_function.h"
#include "source_base/timer.h"
#include "source_basis/module_pw/kernels/pw_op.h"
#include "source_base/module_fft/fft_bundle.h"
#include "pw_basis.h"
#include "pw_gatherscatter.h"

#include <cassert>
#include <complex>

namespace ModulePW
{
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ir = 0; ir < nrxx_; ++ir)
    {
        this->fft_bundle.get_auxr_data<FPTYPE>()[ir] = in[ir];
    }
    this->fft_bundle.fftxyfor(fft_bundle.get_auxr_data<FPTYPE>(), fft_bundle.get_auxr_data<FPTYPE>());

    this->gatherp_scatters(this->fft_bundle.get_auxr_data<FPTYPE>(), this->fft_bundle.get_auxg_data<FPTYPE>());

    this->fft_bundle.fftzfor(fft_bundle.get_auxg_data<FPTYPE>(), fft_bundle.get_auxg_data<FPTYPE>());

    if (add)
    {
        FPTYPE tmpfac = factor / FPTYPE(nxyz_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ig = 0; ig < npw_; ++ig)
        {
            out[ig] += tmpfac * this->fft_bundle.get_auxg_data<FPTYPE>()[ig2isz_[ig]];
        }
    }
    else
    {
        FPTYPE tmpfac = 1.0 / FPTYPE(nxyz_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ig = 0; ig < npw_; ++ig)
        {
            out[ig] = tmpfac * this->fft_bundle.get_auxg_data<FPTYPE>()[ig2isz_[ig]];
        }
    }
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
    if (this->gamma_only)
    {
        const int npy = ny_ * nplane_;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
        for (int ix = 0; ix < nx_; ++ix)
        {
            for (int ipy = 0; ipy < npy; ++ipy)
            {
                this->fft_bundle.get_rspace_data<FPTYPE>()[ix * npy + ipy] = in[ix * npy + ipy];
            }
        }

        this->fft_bundle.fftxyr2c(fft_bundle.get_rspace_data<FPTYPE>(), fft_bundle.get_auxr_data<FPTYPE>());
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0; ir < nrxx_; ++ir)
        {
            this->fft_bundle.get_auxr_data<FPTYPE>()[ir] = std::complex<FPTYPE>(in[ir], 0);
        }
        this->fft_bundle.fftxyfor(fft_bundle.get_auxr_data<FPTYPE>(), fft_bundle.get_auxr_data<FPTYPE>());
    }
    this->gatherp_scatters(this->fft_bundle.get_auxr_data<FPTYPE>(), this->fft_bundle.get_auxg_data<FPTYPE>());

    this->fft_bundle.fftzfor(fft_bundle.get_auxg_data<FPTYPE>(), fft_bundle.get_auxg_data<FPTYPE>());

    if (add)
    {
        FPTYPE tmpfac = factor / FPTYPE(nxyz_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ig = 0; ig < npw_; ++ig)
        {
            out[ig] += tmpfac * this->fft_bundle.get_auxg_data<FPTYPE>()[ig2isz_[ig]];
        }
    }
    else
    {
        FPTYPE tmpfac = 1.0 / FPTYPE(nxyz_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ig = 0; ig < npw_; ++ig)
        {
            out[ig] = tmpfac * this->fft_bundle.get_auxg_data<FPTYPE>()[ig2isz_[ig]];
        }
    }
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < nst_ * nz_; ++i)
    {
        fft_bundle.get_auxg_data<FPTYPE>()[i] = std::complex<FPTYPE>(0, 0);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ig = 0; ig < npw_; ++ig)
    {
        this->fft_bundle.get_auxg_data<FPTYPE>()[ig2isz_[ig]] = in[ig];
    }
    this->fft_bundle.fftzbac(fft_bundle.get_auxg_data<FPTYPE>(), fft_bundle.get_auxg_data<FPTYPE>());

    this->gathers_scatterp(this->fft_bundle.get_auxg_data<FPTYPE>(), this->fft_bundle.get_auxr_data<FPTYPE>());

    this->fft_bundle.fftxybac(fft_bundle.get_auxr_data<FPTYPE>(), this->fft_bundle.get_auxr_data<FPTYPE>());

    if (add)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0; ir < nrxx_; ++ir)
        {
            out[ir] += factor * this->fft_bundle.get_auxr_data<FPTYPE>()[ir];
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0; ir < nrxx_; ++ir)
        {
            out[ir] = this->fft_bundle.get_auxr_data<FPTYPE>()[ir];
        }
    }
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < nst_ * nz_; ++i)
    {
        fft_bundle.get_auxg_data<FPTYPE>()[i] = std::complex<FPTYPE>(0, 0);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ig = 0; ig < npw_; ++ig)
    {
        this->fft_bundle.get_auxg_data<FPTYPE>()[ig2isz_[ig]] = in[ig];
    }
    this->fft_bundle.fftzbac(fft_bundle.get_auxg_data<FPTYPE>(), fft_bundle.get_auxg_data<FPTYPE>());

    this->gathers_scatterp(this->fft_bundle.get_auxg_data<FPTYPE>(), this->fft_bundle.get_auxr_data<FPTYPE>());

    if (this->gamma_only)
    {
        this->fft_bundle.fftxyc2r(fft_bundle.get_auxr_data<FPTYPE>(), fft_bundle.get_rspace_data<FPTYPE>());

        const int npy = ny_ * nplane_;

        if (add)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
            for (int ix = 0; ix < nx_; ++ix)
            {
                for (int ipy = 0; ipy < npy; ++ipy)
                {
                    out[ix * npy + ipy] += factor * this->fft_bundle.get_rspace_data<FPTYPE>()[ix * npy + ipy];
                }
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
            for (int ix = 0; ix < nx_; ++ix)
            {
                for (int ipy = 0; ipy < npy; ++ipy)
                {
                    out[ix * npy + ipy] = this->fft_bundle.get_rspace_data<FPTYPE>()[ix * npy + ipy];
                }
            }
        }
    }
    else
    {
        this->fft_bundle.fftxybac(fft_bundle.get_auxr_data<FPTYPE>(), fft_bundle.get_auxr_data<FPTYPE>());
        if (add)
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int ir = 0; ir < nrxx_; ++ir)
            {
                out[ir] += factor * this->fft_bundle.get_auxr_data<FPTYPE>()[ir].real();
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int ir = 0; ir < nrxx_; ++ir)
            {
                out[ir] = this->fft_bundle.get_auxr_data<FPTYPE>()[ir].real();
            }
        }
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