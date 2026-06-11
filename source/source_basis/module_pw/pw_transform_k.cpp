#include "source_base/timer.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_pw/kernels/pw_op.h"
#include "pw_basis_k.h"
#include "pw_gatherscatter.h"

#include <cassert>
#include <complex>
#include <vector>

#if defined(__ROCM) && !defined(__CUDA)
#error "PW_Basis_K GPU symmetry remap and batch transforms are not implemented for ROCm in this merge."
#endif

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
    ModuleBase::timer::start(this->classname, "real2recip");

    assert(this->gamma_only == false);
    auto* auxr = this->fft_bundle.get_auxr_data<FPTYPE>();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
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
#pragma omp parallel for schedule(static)
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
#pragma omp parallel for schedule(static)
#endif
        for (int igl = 0; igl < npwk; ++igl)
        {
            out[igl] = tmpfac * auxg[this->igl2isz_k[igl + startig]];
        }
    }
    ModuleBase::timer::end(this->classname, "real2recip");
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
    ModuleBase::timer::start(this->classname, "real2recip");
    assert(this->gamma_only == true);
    // for(int ir = 0 ; ir < this->nrxx ; ++ir)
    // {
    //     this->fft_bundle.get_rspace_data<FPTYPE>()[ir] = in[ir];
    // }
    // r2c in place
    const int npy = this->ny * this->nplane;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
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
#pragma omp parallel for schedule(static)
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
#pragma omp parallel for schedule(static)
#endif
        for (int igl = 0; igl < npwk; ++igl)
        {
            out[igl] = tmpfac * auxg[this->igl2isz_k[igl + startig]];
        }
    }
    ModuleBase::timer::end(this->classname, "real2recip");
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
    ModuleBase::timer::start(this->classname, "recip2real");
    assert(this->gamma_only == false);
    ModuleBase::GlobalFunc::ZEROS(fft_bundle.get_auxg_data<FPTYPE>(), this->nst * this->nz);

    const int startig = ik * this->npwk_max;
    const int npwk = this->npwk[ik];
    auto* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
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
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            out[ir] += factor * auxr[ir];
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            out[ir] = auxr[ir];
        }
    }
    ModuleBase::timer::end(this->classname, "recip2real");
}

template <typename FPTYPE>
void PW_Basis_K::recip2real_remapped(const std::complex<FPTYPE>* in,
                                     std::complex<FPTYPE>* out,
                                     const int npw_full,
                                     const int* rep_igl,
                                     const int* fft_isz,
                                     const std::complex<double>* phase,
                                     const bool add,
                                     const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "recip2real_remapped");
    assert(this->gamma_only == false);
#if defined(__ROCM) && !defined(__CUDA)
    if (this->device == "gpu")
    {
        ModuleBase::WARNING_QUIT("PW_Basis_K::recip2real_remapped",
                                 "PW symmetry remap transforms are not implemented for ROCm yet");
    }
#endif
#if defined(__CUDA)
    if (this->device == "gpu")
    {
        if (this->poolnproc != 1)
        {
            ModuleBase::WARNING_QUIT("PW_Basis_K::recip2real_remapped",
                                     "GPU symmetry remap transforms require poolnproc == 1");
        }
        std::complex<FPTYPE>* phase_device = nullptr;
        int* rep_igl_device = nullptr;
        int* fft_ixyz_device = nullptr;
        std::vector<std::complex<FPTYPE>> phase_fptype(npw_full);
        std::vector<int> fft_ixyz(npw_full);
        for (int ig = 0; ig < npw_full; ++ig)
        {
            phase_fptype[ig] = std::complex<FPTYPE>(static_cast<FPTYPE>(phase[ig].real()),
                                                    static_cast<FPTYPE>(phase[ig].imag()));
            const int isz = fft_isz[ig];
            const int iz = isz % this->nz;
            const int is = isz / this->nz;
            const int ixy = this->is2fftixy[is];
            const int iy = ixy % this->ny;
            const int ix = ixy / this->ny;
            fft_ixyz[ig] = iz + iy * this->nz + ix * this->ny * this->nz;
        }
        base_device::memory::resize_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(phase_device, npw_full);
        base_device::memory::resize_memory_op<int, base_device::DEVICE_GPU>()(rep_igl_device, npw_full);
        base_device::memory::resize_memory_op<int, base_device::DEVICE_GPU>()(fft_ixyz_device, npw_full);
        base_device::memory::synchronize_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            phase_device,
            phase_fptype.data(),
            npw_full);
        base_device::memory::synchronize_memory_op<int, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            rep_igl_device,
            rep_igl,
            npw_full);
        base_device::memory::synchronize_memory_op<int, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            fft_ixyz_device,
            fft_ixyz.data(),
            npw_full);

        base_device::memory::set_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(
            this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
            0,
            this->nxyz);
        set_3d_fft_box_remapped_op<FPTYPE, base_device::DEVICE_GPU>()(
            npw_full,
            rep_igl_device,
            fft_ixyz_device,
            phase_device,
            false,
            in,
            this->fft_bundle.get_auxr_3d_data<FPTYPE>());
        this->fft_bundle.fft3D_backward(this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                        this->fft_bundle.get_auxr_3d_data<FPTYPE>());
        set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>()(this->nrxx,
                                                                       add,
                                                                       factor,
                                                                       this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                                                       out);
        base_device::memory::delete_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(phase_device);
        base_device::memory::delete_memory_op<int, base_device::DEVICE_GPU>()(rep_igl_device);
        base_device::memory::delete_memory_op<int, base_device::DEVICE_GPU>()(fft_ixyz_device);
        ModuleBase::timer::end(this->classname, "recip2real_remapped");
        return;
    }
#endif
    ModuleBase::GlobalFunc::ZEROS(fft_bundle.get_auxg_data<FPTYPE>(), this->nst * this->nz);

    auto* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int igl_full = 0; igl_full < npw_full; ++igl_full)
    {
        auxg[fft_isz[igl_full]] = in[rep_igl[igl_full]]
                                  * std::complex<FPTYPE>(static_cast<FPTYPE>(phase[igl_full].real()),
                                                         static_cast<FPTYPE>(phase[igl_full].imag()));
    }
    this->fft_bundle.fftzbac(fft_bundle.get_auxg_data<FPTYPE>(), fft_bundle.get_auxg_data<FPTYPE>());

    this->gathers_scatterp(this->fft_bundle.get_auxg_data<FPTYPE>(), this->fft_bundle.get_auxr_data<FPTYPE>());

    this->fft_bundle.fftxybac(fft_bundle.get_auxr_data<FPTYPE>(), fft_bundle.get_auxr_data<FPTYPE>());
    auto* auxr = this->fft_bundle.get_auxr_data<FPTYPE>();
    if (add)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            out[ir] += factor * auxr[ir];
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            out[ir] = auxr[ir];
        }
    }
    ModuleBase::timer::end(this->classname, "recip2real_remapped");
}

template <typename FPTYPE>
void PW_Basis_K::real2recip_remapped(const std::complex<FPTYPE>* in,
                                     std::complex<FPTYPE>* out,
                                     const int npw_full,
                                     const int* rep_igl,
                                     const int* fft_isz,
                                     const std::complex<double>* phase,
                                     const bool add,
                                     const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "real2recip_remapped");
    assert(this->gamma_only == false);
#if defined(__ROCM) && !defined(__CUDA)
    if (this->device == "gpu")
    {
        ModuleBase::WARNING_QUIT("PW_Basis_K::real2recip_remapped",
                                 "PW symmetry remap transforms are not implemented for ROCm yet");
    }
#endif
#if defined(__CUDA)
    if (this->device == "gpu")
    {
        if (this->poolnproc != 1)
        {
            ModuleBase::WARNING_QUIT("PW_Basis_K::real2recip_remapped",
                                     "GPU symmetry remap transforms require poolnproc == 1");
        }
        std::complex<FPTYPE>* phase_device = nullptr;
        int* rep_igl_device = nullptr;
        int* fft_ixyz_device = nullptr;
        std::vector<std::complex<FPTYPE>> phase_fptype(npw_full);
        std::vector<int> fft_ixyz(npw_full);
        for (int ig = 0; ig < npw_full; ++ig)
        {
            phase_fptype[ig] = std::complex<FPTYPE>(static_cast<FPTYPE>(phase[ig].real()),
                                                    static_cast<FPTYPE>(phase[ig].imag()));
            const int isz = fft_isz[ig];
            const int iz = isz % this->nz;
            const int is = isz / this->nz;
            const int ixy = this->is2fftixy[is];
            const int iy = ixy % this->ny;
            const int ix = ixy / this->ny;
            fft_ixyz[ig] = iz + iy * this->nz + ix * this->ny * this->nz;
        }
        base_device::memory::resize_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(phase_device, npw_full);
        base_device::memory::resize_memory_op<int, base_device::DEVICE_GPU>()(rep_igl_device, npw_full);
        base_device::memory::resize_memory_op<int, base_device::DEVICE_GPU>()(fft_ixyz_device, npw_full);
        base_device::memory::synchronize_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            phase_device,
            phase_fptype.data(),
            npw_full);
        base_device::memory::synchronize_memory_op<int, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            rep_igl_device,
            rep_igl,
            npw_full);
        base_device::memory::synchronize_memory_op<int, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            fft_ixyz_device,
            fft_ixyz.data(),
            npw_full);

        base_device::memory::synchronize_memory_op<std::complex<FPTYPE>,
                                                   base_device::DEVICE_GPU,
                                                   base_device::DEVICE_GPU>()(
            this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
            in,
            this->nrxx);
        this->fft_bundle.fft3D_forward(this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                       this->fft_bundle.get_auxr_3d_data<FPTYPE>());
        set_real_to_recip_remapped_output_op<FPTYPE, base_device::DEVICE_GPU>()(
            npw_full,
            this->nxyz,
            add,
            factor,
            rep_igl_device,
            fft_ixyz_device,
            phase_device,
            false,
            this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
            out);
        base_device::memory::delete_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(phase_device);
        base_device::memory::delete_memory_op<int, base_device::DEVICE_GPU>()(rep_igl_device);
        base_device::memory::delete_memory_op<int, base_device::DEVICE_GPU>()(fft_ixyz_device);
        ModuleBase::timer::end(this->classname, "real2recip_remapped");
        return;
    }
#endif

    auto* auxr = this->fft_bundle.get_auxr_data<FPTYPE>();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ir = 0; ir < this->nrxx; ++ir)
    {
        auxr[ir] = in[ir];
    }
    this->fft_bundle.fftxyfor(fft_bundle.get_auxr_data<FPTYPE>(), fft_bundle.get_auxr_data<FPTYPE>());

    this->gatherp_scatters(this->fft_bundle.get_auxr_data<FPTYPE>(), this->fft_bundle.get_auxg_data<FPTYPE>());

    this->fft_bundle.fftzfor(fft_bundle.get_auxg_data<FPTYPE>(), fft_bundle.get_auxg_data<FPTYPE>());

    if (!add)
    {
        for (int ig = 0; ig < npw_full; ++ig)
        {
            out[rep_igl[ig]] = 0.0;
        }
    }

    const FPTYPE tmpfac = (add ? factor : FPTYPE(1.0)) / FPTYPE(this->nxyz);
    auto* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
    for (int igl_full = 0; igl_full < npw_full; ++igl_full)
    {
        const std::complex<FPTYPE> phase_conj(static_cast<FPTYPE>(phase[igl_full].real()),
                                              -static_cast<FPTYPE>(phase[igl_full].imag()));
        const std::complex<FPTYPE> contrib = tmpfac * auxg[fft_isz[igl_full]] * phase_conj;
        out[rep_igl[igl_full]] += contrib;
    }
    ModuleBase::timer::end(this->classname, "real2recip_remapped");
}

template <typename FPTYPE>
void PW_Basis_K::recip2real_remapped_conjugate(const std::complex<FPTYPE>* in,
                                               std::complex<FPTYPE>* out,
                                               const int npw_full,
                                               const int* rep_igl,
                                               const int* fft_isz,
                                               const std::complex<double>* phase,
                                               const bool add,
                                               const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "recip2real_remapped_conj");
    assert(this->gamma_only == false);
#if defined(__ROCM) && !defined(__CUDA)
    if (this->device == "gpu")
    {
        ModuleBase::WARNING_QUIT("PW_Basis_K::recip2real_remapped_conjugate",
                                 "PW symmetry remap transforms are not implemented for ROCm yet");
    }
#endif
#if defined(__CUDA)
    if (this->device == "gpu")
    {
        if (this->poolnproc != 1)
        {
            ModuleBase::WARNING_QUIT("PW_Basis_K::recip2real_remapped_conjugate",
                                     "GPU symmetry remap transforms require poolnproc == 1");
        }
        std::complex<FPTYPE>* phase_device = nullptr;
        int* rep_igl_device = nullptr;
        int* fft_ixyz_device = nullptr;
        std::vector<std::complex<FPTYPE>> phase_fptype(npw_full);
        std::vector<int> fft_ixyz(npw_full);
        for (int ig = 0; ig < npw_full; ++ig)
        {
            phase_fptype[ig] = std::complex<FPTYPE>(static_cast<FPTYPE>(phase[ig].real()),
                                                    static_cast<FPTYPE>(phase[ig].imag()));
            const int isz = fft_isz[ig];
            const int iz = isz % this->nz;
            const int is = isz / this->nz;
            const int ixy = this->is2fftixy[is];
            const int iy = ixy % this->ny;
            const int ix = ixy / this->ny;
            fft_ixyz[ig] = iz + iy * this->nz + ix * this->ny * this->nz;
        }
        base_device::memory::resize_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(phase_device, npw_full);
        base_device::memory::resize_memory_op<int, base_device::DEVICE_GPU>()(rep_igl_device, npw_full);
        base_device::memory::resize_memory_op<int, base_device::DEVICE_GPU>()(fft_ixyz_device, npw_full);
        base_device::memory::synchronize_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            phase_device,
            phase_fptype.data(),
            npw_full);
        base_device::memory::synchronize_memory_op<int, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            rep_igl_device,
            rep_igl,
            npw_full);
        base_device::memory::synchronize_memory_op<int, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            fft_ixyz_device,
            fft_ixyz.data(),
            npw_full);

        base_device::memory::set_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(
            this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
            0,
            this->nxyz);
        set_3d_fft_box_remapped_op<FPTYPE, base_device::DEVICE_GPU>()(
            npw_full,
            rep_igl_device,
            fft_ixyz_device,
            phase_device,
            true,
            in,
            this->fft_bundle.get_auxr_3d_data<FPTYPE>());
        this->fft_bundle.fft3D_backward(this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                        this->fft_bundle.get_auxr_3d_data<FPTYPE>());
        set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>()(this->nrxx,
                                                                       add,
                                                                       factor,
                                                                       this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                                                       out);
        base_device::memory::delete_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(phase_device);
        base_device::memory::delete_memory_op<int, base_device::DEVICE_GPU>()(rep_igl_device);
        base_device::memory::delete_memory_op<int, base_device::DEVICE_GPU>()(fft_ixyz_device);
        ModuleBase::timer::end(this->classname, "recip2real_remapped_conj");
        return;
    }
#endif
    ModuleBase::GlobalFunc::ZEROS(fft_bundle.get_auxg_data<FPTYPE>(), this->nst * this->nz);

    auto* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int igl_full = 0; igl_full < npw_full; ++igl_full)
    {
        const std::complex<FPTYPE> phase_conj(static_cast<FPTYPE>(phase[igl_full].real()),
                                              -static_cast<FPTYPE>(phase[igl_full].imag()));
        auxg[fft_isz[igl_full]] = std::conj(in[rep_igl[igl_full]]) * phase_conj;
    }
    this->fft_bundle.fftzbac(fft_bundle.get_auxg_data<FPTYPE>(), fft_bundle.get_auxg_data<FPTYPE>());

    this->gathers_scatterp(this->fft_bundle.get_auxg_data<FPTYPE>(), this->fft_bundle.get_auxr_data<FPTYPE>());

    this->fft_bundle.fftxybac(fft_bundle.get_auxr_data<FPTYPE>(), fft_bundle.get_auxr_data<FPTYPE>());
    auto* auxr = this->fft_bundle.get_auxr_data<FPTYPE>();
    if (add)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            out[ir] += factor * auxr[ir];
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            out[ir] = auxr[ir];
        }
    }
    ModuleBase::timer::end(this->classname, "recip2real_remapped_conj");
}

template <typename FPTYPE, typename Device>
void PW_Basis_K::recip2real_remapped_batch(const Device* ctx,
                                           const std::complex<FPTYPE>* in_batch,
                                           std::complex<FPTYPE>* out_batch,
                                           const int npw_full,
                                           const int* rep_igl,
                                           const int* fft_isz,
                                           const int* fft_ixyz,
                                           const std::complex<double>* phase,
                                           const std::complex<FPTYPE>* phase_device,
                                           int batch_count,
                                           const bool conjugate,
                                           const bool add,
                                           const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "recip2real_remapped_batch");
    assert(this->gamma_only == false);
#if defined(__ROCM) && !defined(__CUDA)
    if (this->device == "gpu")
    {
        ModuleBase::WARNING_QUIT("PW_Basis_K::recip2real_remapped_batch",
                                 "PW symmetry remap batch transforms are not implemented for ROCm yet");
    }
#endif
#if defined(__CUDA)
    if (this->device == "gpu")
    {
        if (this->poolnproc != 1)
        {
            ModuleBase::WARNING_QUIT("PW_Basis_K::recip2real_remapped_batch",
                                     "GPU symmetry remap batch transforms require poolnproc == 1");
        }
        if (!this->fft_bundle.is_batch_fft_available<FPTYPE>()
            || batch_count > this->fft_bundle.get_batch_size<FPTYPE>())
        {
            ModuleBase::WARNING_QUIT("PW_Basis_K::recip2real_remapped_batch",
                                     "GPU symmetry remap batch transform requires available batch FFT");
        }
        if (fft_ixyz == nullptr || phase_device == nullptr)
        {
            ModuleBase::WARNING_QUIT("PW_Basis_K::recip2real_remapped_batch",
                                     "GPU symmetry remap batch transform requires cached device remap data");
        }

        std::complex<FPTYPE>* batch_in = this->fft_bundle.get_batch_input_buffer<FPTYPE>();
        std::complex<FPTYPE>* batch_out = this->fft_bundle.get_batch_output_buffer<FPTYPE>();
        const int nxyz = this->nxyz;
        base_device::memory::set_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(
            batch_in,
            0,
            static_cast<std::size_t>(batch_count) * nxyz);
        set_3d_fft_box_remapped_op<FPTYPE, base_device::DEVICE_GPU>().operator_batch(
            npw_full,
            nxyz,
            this->npwk_max,
            rep_igl,
            fft_ixyz,
            phase_device,
            conjugate,
            batch_count,
            in_batch,
            batch_in);
        this->fft_bundle.fft3D_backward_batch(batch_in, batch_out, batch_count);
        set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>().operator_batch(
            this->nrxx,
            nxyz,
            batch_count,
            add,
            factor,
            batch_out,
            out_batch);
        ModuleBase::timer::end(this->classname, "recip2real_remapped_batch");
        return;
    }
#endif
    for (int ib = 0; ib < batch_count; ++ib)
    {
        const std::complex<FPTYPE>* in = in_batch + static_cast<std::size_t>(ib) * this->npwk_max;
        std::complex<FPTYPE>* out = out_batch + static_cast<std::size_t>(ib) * this->nrxx;
        if (conjugate)
        {
            this->recip2real_remapped_conjugate(in, out, npw_full, rep_igl, fft_isz, phase, add, factor);
        }
        else
        {
            this->recip2real_remapped(in, out, npw_full, rep_igl, fft_isz, phase, add, factor);
        }
    }
    ModuleBase::timer::end(this->classname, "recip2real_remapped_batch");
}

template <typename FPTYPE>
void PW_Basis_K::real2recip_remapped_conjugate(const std::complex<FPTYPE>* in,
                                               std::complex<FPTYPE>* out,
                                               const int npw_full,
                                               const int* rep_igl,
                                               const int* fft_isz,
                                               const std::complex<double>* phase,
                                               const bool add,
                                               const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "real2recip_remapped_conj");
    assert(this->gamma_only == false);
#if defined(__ROCM) && !defined(__CUDA)
    if (this->device == "gpu")
    {
        ModuleBase::WARNING_QUIT("PW_Basis_K::real2recip_remapped_conjugate",
                                 "PW symmetry remap transforms are not implemented for ROCm yet");
    }
#endif
#if defined(__CUDA)
    if (this->device == "gpu")
    {
        if (this->poolnproc != 1)
        {
            ModuleBase::WARNING_QUIT("PW_Basis_K::real2recip_remapped_conjugate",
                                     "GPU symmetry remap transforms require poolnproc == 1");
        }
        std::complex<FPTYPE>* phase_device = nullptr;
        int* rep_igl_device = nullptr;
        int* fft_ixyz_device = nullptr;
        std::vector<std::complex<FPTYPE>> phase_fptype(npw_full);
        std::vector<int> fft_ixyz(npw_full);
        for (int ig = 0; ig < npw_full; ++ig)
        {
            phase_fptype[ig] = std::complex<FPTYPE>(static_cast<FPTYPE>(phase[ig].real()),
                                                    static_cast<FPTYPE>(phase[ig].imag()));
            const int isz = fft_isz[ig];
            const int iz = isz % this->nz;
            const int is = isz / this->nz;
            const int ixy = this->is2fftixy[is];
            const int iy = ixy % this->ny;
            const int ix = ixy / this->ny;
            fft_ixyz[ig] = iz + iy * this->nz + ix * this->ny * this->nz;
        }
        base_device::memory::resize_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(phase_device, npw_full);
        base_device::memory::resize_memory_op<int, base_device::DEVICE_GPU>()(rep_igl_device, npw_full);
        base_device::memory::resize_memory_op<int, base_device::DEVICE_GPU>()(fft_ixyz_device, npw_full);
        base_device::memory::synchronize_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            phase_device,
            phase_fptype.data(),
            npw_full);
        base_device::memory::synchronize_memory_op<int, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            rep_igl_device,
            rep_igl,
            npw_full);
        base_device::memory::synchronize_memory_op<int, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            fft_ixyz_device,
            fft_ixyz.data(),
            npw_full);

        base_device::memory::synchronize_memory_op<std::complex<FPTYPE>,
                                                   base_device::DEVICE_GPU,
                                                   base_device::DEVICE_GPU>()(
            this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
            in,
            this->nrxx);
        this->fft_bundle.fft3D_forward(this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                       this->fft_bundle.get_auxr_3d_data<FPTYPE>());
        set_real_to_recip_remapped_output_op<FPTYPE, base_device::DEVICE_GPU>()(
            npw_full,
            this->nxyz,
            add,
            factor,
            rep_igl_device,
            fft_ixyz_device,
            phase_device,
            true,
            this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
            out);
        base_device::memory::delete_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(phase_device);
        base_device::memory::delete_memory_op<int, base_device::DEVICE_GPU>()(rep_igl_device);
        base_device::memory::delete_memory_op<int, base_device::DEVICE_GPU>()(fft_ixyz_device);
        ModuleBase::timer::end(this->classname, "real2recip_remapped_conj");
        return;
    }
#endif

    auto* auxr = this->fft_bundle.get_auxr_data<FPTYPE>();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ir = 0; ir < this->nrxx; ++ir)
    {
        auxr[ir] = in[ir];
    }
    this->fft_bundle.fftxyfor(fft_bundle.get_auxr_data<FPTYPE>(), fft_bundle.get_auxr_data<FPTYPE>());

    this->gatherp_scatters(this->fft_bundle.get_auxr_data<FPTYPE>(), this->fft_bundle.get_auxg_data<FPTYPE>());

    this->fft_bundle.fftzfor(fft_bundle.get_auxg_data<FPTYPE>(), fft_bundle.get_auxg_data<FPTYPE>());

    if (!add)
    {
        for (int ig = 0; ig < npw_full; ++ig)
        {
            out[rep_igl[ig]] = 0.0;
        }
    }

    const FPTYPE tmpfac = (add ? factor : FPTYPE(1.0)) / FPTYPE(this->nxyz);
    auto* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
    for (int igl_full = 0; igl_full < npw_full; ++igl_full)
    {
        const std::complex<FPTYPE> phase_value(static_cast<FPTYPE>(phase[igl_full].real()),
                                               static_cast<FPTYPE>(phase[igl_full].imag()));
        const std::complex<FPTYPE> contrib = tmpfac * std::conj(auxg[fft_isz[igl_full]]) * phase_value;
        out[rep_igl[igl_full]] += contrib;
    }
    ModuleBase::timer::end(this->classname, "real2recip_remapped_conj");
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
    ModuleBase::timer::start(this->classname, "recip2real");
    assert(this->gamma_only == true);
    ModuleBase::GlobalFunc::ZEROS(fft_bundle.get_auxg_data<FPTYPE>(), this->nst * this->nz);

    const int startig = ik * this->npwk_max;
    const int npwk = this->npwk[ik];
    auto* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
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
#pragma omp parallel for collapse(2) schedule(static)
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
#pragma omp parallel for collapse(2) schedule(static)
#endif
        for (int ix = 0; ix < this->nx; ++ix)
        {
            for (int ipy = 0; ipy < npy; ++ipy)
            {
                out[ix * npy + ipy] = rspace[ix * npy + ipy];
            }
        }
    }
    ModuleBase::timer::end(this->classname, "recip2real");
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

#if (defined(__CUDA) || defined(__ROCM))
template <>
void PW_Basis_K::real_to_recip(const base_device::DEVICE_GPU* ctx,
                               const std::complex<float>* in,
                               std::complex<float>* out,
                               const int ik,
                               const bool add,
                               const float factor) const
{
    ModuleBase::timer::start(this->classname, "real_to_recip gpu");
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
    ModuleBase::timer::end(this->classname, "real_to_recip gpu");
}
template <>
void PW_Basis_K::real_to_recip(const base_device::DEVICE_GPU* ctx,
                               const std::complex<double>* in,
                               std::complex<double>* out,
                               const int ik,
                               const bool add,
                               const double factor) const
{
    ModuleBase::timer::start(this->classname, "real_to_recip gpu");
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
    ModuleBase::timer::end(this->classname, "real_to_recip gpu");
}

template <>
void PW_Basis_K::recip_to_real(const base_device::DEVICE_GPU* ctx,
                               const std::complex<float>* in,
                               std::complex<float>* out,
                               const int ik,
                               const bool add,
                               const float factor) const
{
    ModuleBase::timer::start(this->classname, "recip_to_real gpu");
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

    ModuleBase::timer::end(this->classname, "recip_to_real gpu");
}
template <>
void PW_Basis_K::recip_to_real(const base_device::DEVICE_GPU* ctx,
                               const std::complex<double>* in,
                               std::complex<double>* out,
                               const int ik,
                               const bool add,
                               const double factor) const
{
    ModuleBase::timer::start(this->classname, "recip_to_real gpu");
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

    ModuleBase::timer::end(this->classname, "recip_to_real gpu");
}

template <typename FPTYPE>
void PW_Basis_K::real2recip_gpu(const std::complex<FPTYPE>* in,
                               std::complex<FPTYPE>* out,
                               const int ik,
                               const bool add,
                               const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "real_to_recip gpu");
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
    ModuleBase::timer::end(this->classname, "real_to_recip gpu");
}
template <typename FPTYPE>
void PW_Basis_K::recip2real_gpu(const std::complex<FPTYPE>* in,
                               std::complex<FPTYPE>* out,
                               const int ik,
                               const bool add,
                               const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "recip_to_real gpu");
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

    ModuleBase::timer::end(this->classname, "recip_to_real gpu");
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
template void PW_Basis_K::recip2real_remapped<float>(const std::complex<float>* in,
                                                     std::complex<float>* out,
                                                     const int npw_full,
                                                     const int* rep_igl,
                                                     const int* fft_isz,
                                                     const std::complex<double>* phase,
                                                     const bool add,
                                                     const float factor) const;
template void PW_Basis_K::real2recip_remapped<float>(const std::complex<float>* in,
                                                     std::complex<float>* out,
                                                     const int npw_full,
                                                     const int* rep_igl,
                                                     const int* fft_isz,
                                                     const std::complex<double>* phase,
                                                     const bool add,
                                                     const float factor) const;
template void PW_Basis_K::recip2real_remapped_conjugate<float>(const std::complex<float>* in,
                                                               std::complex<float>* out,
                                                               const int npw_full,
                                                               const int* rep_igl,
                                                               const int* fft_isz,
                                                               const std::complex<double>* phase,
                                                               const bool add,
                                                               const float factor) const;
template void PW_Basis_K::real2recip_remapped_conjugate<float>(const std::complex<float>* in,
                                                               std::complex<float>* out,
                                                               const int npw_full,
                                                               const int* rep_igl,
                                                               const int* fft_isz,
                                                               const std::complex<double>* phase,
                                                               const bool add,
                                                               const float factor) const;

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
template void PW_Basis_K::recip2real_remapped<double>(const std::complex<double>* in,
                                                      std::complex<double>* out,
                                                      const int npw_full,
                                                      const int* rep_igl,
                                                      const int* fft_isz,
                                                      const std::complex<double>* phase,
                                                      const bool add,
                                                      const double factor) const;
template void PW_Basis_K::real2recip_remapped<double>(const std::complex<double>* in,
                                                      std::complex<double>* out,
                                                      const int npw_full,
                                                      const int* rep_igl,
                                                      const int* fft_isz,
                                                      const std::complex<double>* phase,
                                                      const bool add,
                                                      const double factor) const;
template void PW_Basis_K::recip2real_remapped_conjugate<double>(const std::complex<double>* in,
                                                                std::complex<double>* out,
                                                                const int npw_full,
                                                                const int* rep_igl,
                                                                const int* fft_isz,
                                                                const std::complex<double>* phase,
                                                                const bool add,
                                                                const double factor) const;
template void PW_Basis_K::real2recip_remapped_conjugate<double>(const std::complex<double>* in,
                                                                std::complex<double>* out,
                                                                const int npw_full,
                                                                const int* rep_igl,
                                                                const int* fft_isz,
                                                                const std::complex<double>* phase,
                                                                const bool add,
                                                                const double factor) const;
// ============================================================================
// Batch Transform Implementations
// ============================================================================

#if defined(__CUDA) || defined(__ROCM)

template <typename FPTYPE, typename Device>
void PW_Basis_K::real_to_recip_batch(const Device* ctx,
                                     const std::complex<FPTYPE>* in_batch,
                                     std::complex<FPTYPE>* out_batch,
                                     const int ik,
                                     int batch_count,
                                     const bool add,
                                     const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "real_to_recip_batch gpu");

    // Check if batch FFT is available
    if (!this->fft_bundle.is_batch_fft_available<FPTYPE>() ||
        batch_count > this->fft_bundle.get_batch_size<FPTYPE>())
    {
        // Fallback to sequential transforms using GPU-specific template overload
        for (int ib = 0; ib < batch_count; ++ib)
        {
            this->real_to_recip<std::complex<FPTYPE>, Device>(
                               in_batch + ib * this->nrxx,
                               out_batch + ib * this->npwk_max,
                               ik,
                               add,
                               factor);
        }
        ModuleBase::timer::end(this->classname, "real_to_recip_batch gpu");
        return;
    }

    // Batch FFT path
    assert(this->gamma_only == false);
    assert(this->poolnproc == 1);

    // Get batch output buffer from FFT bundle
    std::complex<FPTYPE>* batch_out = this->fft_bundle.get_batch_output_buffer<FPTYPE>();
    const int nxyz = this->nxyz;

    // Perform batch forward FFT directly on input (out-of-place FFT doesn't modify input)
    // const_cast is safe because cuFFT out-of-place C2C transform only reads from input
    this->fft_bundle.fft3D_forward_batch(const_cast<std::complex<FPTYPE>*>(in_batch), batch_out, batch_count);

    // Extract results using batch operator (all batch elements share same ik)
    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];

#if defined(__CUDA)
    // Call operator_batch once
    set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>().operator_batch(
        npw_k,
        nxyz,
        this->npwk_max,     // out_stride: matches output buffer layout
        batch_count,
        add,
        factor,
        this->ig2ixyz_k + startig,
        batch_out,
        out_batch);
#endif

    ModuleBase::timer::end(this->classname, "real_to_recip_batch gpu");
}

template <typename FPTYPE, typename Device>
void PW_Basis_K::recip_to_real_batch(const Device* ctx,
                                     const std::complex<FPTYPE>* in_batch,
                                     std::complex<FPTYPE>* out_batch,
                                     const int ik,
                                     int batch_count,
                                     const bool add,
                                     const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "recip_to_real_batch gpu");

    // Check if batch FFT is available
    if (!this->fft_bundle.is_batch_fft_available<FPTYPE>() ||
        batch_count > this->fft_bundle.get_batch_size<FPTYPE>())
    {
        // Fallback to sequential transforms using GPU-specific template overload
        for (int ib = 0; ib < batch_count; ++ib)
        {
            this->recip_to_real<std::complex<FPTYPE>, Device>(
                               in_batch + ib * this->npwk_max,
                               out_batch + ib * this->nrxx,
                               ik,
                               add,
                               factor);
        }
        ModuleBase::timer::end(this->classname, "recip_to_real_batch gpu");
        return;
    }

    // Batch FFT path
    assert(this->gamma_only == false);
    assert(this->poolnproc == 1);

    // Get batch buffers from FFT bundle
    std::complex<FPTYPE>* batch_in = this->fft_bundle.get_batch_input_buffer<FPTYPE>();
    std::complex<FPTYPE>* batch_out = this->fft_bundle.get_batch_output_buffer<FPTYPE>();
    const int nxyz = this->nxyz;

    // Zero batch input buffer
    base_device::memory::set_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(
        batch_in,
        0,
        batch_count * nxyz);

    // Populate FFT input for all batch elements (all share same ik)
    const int startig = ik * this->npwk_max;
    const int npw_k = this->npwk[ik];

    // Use batched kernel with separate count (npw_k) and stride (npwk_max)
    set_3d_fft_box_op<FPTYPE, base_device::DEVICE_GPU>().operator_batch(
        npw_k,              // count: actual plane waves for this k-point
        nxyz,
        this->npwk_max,     // in_stride: matches input data layout
        this->ig2ixyz_k + startig,
        batch_count,
        in_batch,
        batch_in);

    // Perform batch backward FFT
    this->fft_bundle.fft3D_backward_batch(batch_in, batch_out, batch_count);

#if defined(__CUDA)
    // Extract results using the scalar factor shared by all batch elements.
    set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>().operator_batch(
        this->nrxx,
        nxyz,
        batch_count,
        add,
        factor,
        batch_out,
        out_batch);
#endif

    ModuleBase::timer::end(this->classname, "recip_to_real_batch gpu");
}

// Template instantiations for batch transforms
template void PW_Basis_K::real_to_recip_batch<float, base_device::DEVICE_GPU>(
    const base_device::DEVICE_GPU*,
    const std::complex<float>*,
    std::complex<float>*,
    const int,
    int,
    const bool,
    const float) const;

template void PW_Basis_K::real_to_recip_batch<double, base_device::DEVICE_GPU>(
    const base_device::DEVICE_GPU*,
    const std::complex<double>*,
    std::complex<double>*,
    const int,
    int,
    const bool,
    const double) const;

template void PW_Basis_K::recip_to_real_batch<float, base_device::DEVICE_GPU>(
    const base_device::DEVICE_GPU*,
    const std::complex<float>*,
    std::complex<float>*,
    const int,
    int,
    const bool,
    const float) const;

template void PW_Basis_K::recip_to_real_batch<double, base_device::DEVICE_GPU>(
    const base_device::DEVICE_GPU*,
    const std::complex<double>*,
    std::complex<double>*,
    const int,
    int,
    const bool,
    const double) const;

template void PW_Basis_K::recip2real_remapped_batch<float, base_device::DEVICE_GPU>(
    const base_device::DEVICE_GPU*,
    const std::complex<float>*,
    std::complex<float>*,
    const int,
    const int*,
    const int*,
    const int*,
    const std::complex<double>*,
    const std::complex<float>*,
    int,
    const bool,
    const bool,
    const float) const;

template void PW_Basis_K::recip2real_remapped_batch<double, base_device::DEVICE_GPU>(
    const base_device::DEVICE_GPU*,
    const std::complex<double>*,
    std::complex<double>*,
    const int,
    const int*,
    const int*,
    const int*,
    const std::complex<double>*,
    const std::complex<double>*,
    int,
    const bool,
    const bool,
    const double) const;

#else
// CPU fallback (not implemented - just use sequential)
template <typename FPTYPE, typename Device>
void PW_Basis_K::real_to_recip_batch(const Device* ctx,
                                     const std::complex<FPTYPE>* in_batch,
                                     std::complex<FPTYPE>* out_batch,
                                     const int ik,
                                     int batch_count,
                                     const bool add,
                                     const FPTYPE factor) const
{
    // Fallback to sequential transforms on CPU
    for (int ib = 0; ib < batch_count; ++ib)
    {
        this->real_to_recip(ctx,
                           in_batch + ib * this->nrxx,
                           out_batch + ib * this->npwk_max,
                           ik,
                           add,
                           factor);
    }
}

template <typename FPTYPE, typename Device>
void PW_Basis_K::recip_to_real_batch(const Device* ctx,
                                     const std::complex<FPTYPE>* in_batch,
                                     std::complex<FPTYPE>* out_batch,
                                     const int ik,
                                     int batch_count,
                                     const bool add,
                                     const FPTYPE factor) const
{
    // Fallback to sequential transforms on CPU
    for (int ib = 0; ib < batch_count; ++ib)
    {
        this->recip_to_real(ctx,
                           in_batch + ib * this->npwk_max,
                           out_batch + ib * this->nrxx,
                           ik,
                           add,
                           factor);
    }
}
#endif

// Explicit template instantiations for CPU batch transforms
template void PW_Basis_K::real_to_recip_batch<float, base_device::DEVICE_CPU>(
    const base_device::DEVICE_CPU* ctx,
    const std::complex<float>* in_batch,
    std::complex<float>* out_batch,
    const int ik,
    int batch_count,
    const bool add,
    const float factor) const;

template void PW_Basis_K::real_to_recip_batch<double, base_device::DEVICE_CPU>(
    const base_device::DEVICE_CPU* ctx,
    const std::complex<double>* in_batch,
    std::complex<double>* out_batch,
    const int ik,
    int batch_count,
    const bool add,
    const double factor) const;

template void PW_Basis_K::recip_to_real_batch<float, base_device::DEVICE_CPU>(
    const base_device::DEVICE_CPU* ctx,
    const std::complex<float>* in_batch,
    std::complex<float>* out_batch,
    const int ik,
    int batch_count,
    const bool add,
    const float factor) const;

template void PW_Basis_K::recip_to_real_batch<double, base_device::DEVICE_CPU>(
    const base_device::DEVICE_CPU* ctx,
    const std::complex<double>* in_batch,
    std::complex<double>* out_batch,
    const int ik,
    int batch_count,
    const bool add,
    const double factor) const;

template void PW_Basis_K::recip2real_remapped_batch<float, base_device::DEVICE_CPU>(
    const base_device::DEVICE_CPU*,
    const std::complex<float>*,
    std::complex<float>*,
    const int,
    const int*,
    const int*,
    const int*,
    const std::complex<double>*,
    const std::complex<float>*,
    int,
    const bool,
    const bool,
    const float) const;

template void PW_Basis_K::recip2real_remapped_batch<double, base_device::DEVICE_CPU>(
    const base_device::DEVICE_CPU*,
    const std::complex<double>*,
    std::complex<double>*,
    const int,
    const int*,
    const int*,
    const int*,
    const std::complex<double>*,
    const std::complex<double>*,
    int,
    const bool,
    const bool,
    const double) const;

} // namespace ModulePW
