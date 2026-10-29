#include "module_fft/fft_bundle.h"
#include <complex>
#include "pw_basis.h"
#include "kernels/pw_op.h"
#include <algorithm>
#include <type_traits>
#include <cassert>
#include "module_base/global_function.h"
#include "module_base/timer.h"
#include "pw_gatherscatter.h"

namespace ModulePW {

#if defined(__CUDA)
template <typename Input, typename Real>
void PW_Basis::real2recip_gpu(const Input* in, std::complex<Real>* out, bool add, Real factor) const
{
    using Complex = std::complex<Real>;
    int scratch_size = this->npw;
    if constexpr (std::is_same<Input, Real>::value)
    {
        scratch_size = std::max(scratch_size, (this->nrxx + 1) / 2);
    }
    Complex* scratch = nullptr;
    base_device::memory::resize_memory_op<Complex, base_device::DEVICE_GPU>()(gpu_ctx, scratch, scratch_size);
    auto* grid = this->fft_bundle.get_auxr_3d_data<Real>();
    if constexpr (std::is_same<Input, Complex>::value)
    {
        base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            gpu_ctx, cpu_ctx, grid, in, this->nrxx);
    }
    else
    {
        // A complex allocation also stores the contiguous real input before its FFT coefficients are gathered.
        auto* real_input = reinterpret_cast<Real*>(scratch);
        base_device::memory::synchronize_memory_op<Real, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            gpu_ctx, cpu_ctx, real_input, in, this->nrxx);
        pack_real_grid_gpu(this->nrxx, real_input, grid);
    }
    this->fft_bundle.fft3D_forward(gpu_ctx, grid, grid);
    if (add)
    {
        base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            gpu_ctx, cpu_ctx, scratch, out, this->npw);
    }
    set_real_to_recip_output_op<Real, base_device::DEVICE_GPU>()(
        gpu_ctx, this->npw, this->nxyz, add, factor, this->d_ig2fftbox, grid, scratch);
    base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(
        cpu_ctx, gpu_ctx, out, scratch, this->npw);
    base_device::memory::delete_memory_op<Complex, base_device::DEVICE_GPU>()(gpu_ctx, scratch);
}

template <typename Real, typename Output>
void PW_Basis::recip2real_gpu(const std::complex<Real>* in, Output* out, bool add, Real factor) const
{
    using Complex = std::complex<Real>;
    int scratch_size = this->npw;
    if constexpr (std::is_same<Output, Real>::value)
    {
        scratch_size = std::max(scratch_size, (this->nrxx + 1) / 2);
    }
    else if (add)
    {
        scratch_size = std::max(scratch_size, this->nrxx);
    }
    Complex* scratch = nullptr;
    base_device::memory::resize_memory_op<Complex, base_device::DEVICE_GPU>()(gpu_ctx, scratch, scratch_size);
    base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
        gpu_ctx, cpu_ctx, scratch, in, this->npw);
    auto* grid = this->fft_bundle.get_auxr_3d_data<Real>();
    base_device::memory::set_memory_op<Complex, base_device::DEVICE_GPU>()(gpu_ctx, grid, 0, this->nxyz);
    set_3d_fft_box_op<Real, base_device::DEVICE_GPU>()(gpu_ctx, this->npw, this->d_ig2fftbox, scratch, grid);
    this->fft_bundle.fft3D_backward(gpu_ctx, grid, grid);
    auto* result = reinterpret_cast<Output*>(scratch);
    if (add)
    {
        base_device::memory::synchronize_memory_op<Output, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            gpu_ctx, cpu_ctx, result, out, this->nrxx);
    }
    if constexpr (std::is_same<Output, Real>::value)
    {
        extract_real_grid_gpu(this->nrxx, add, factor, grid, result);
    }
    else if (add)
    {
        set_recip_to_real_output_op<Real, base_device::DEVICE_GPU>()(
            gpu_ctx, this->nrxx, true, factor, grid, result);
    }
    else
    {
        result = grid;
    }
    base_device::memory::synchronize_memory_op<Output, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(
        cpu_ctx, gpu_ctx, out, result, this->nrxx);
    base_device::memory::delete_memory_op<Complex, base_device::DEVICE_GPU>()(gpu_ctx, scratch);
}
#endif


#if defined(__CUDA) || defined(__ROCM)
template <typename FPTYPE>
void PW_Basis::fft3d_gpu(bool inverse) const
{
    using Complex = std::complex<FPTYPE>;
    auto* host_grid = this->fft_bundle.get_auxr_data<FPTYPE>();
    auto* device_grid = this->fft_bundle.get_auxr_3d_data<FPTYPE>();
    base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
        gpu_ctx, cpu_ctx, device_grid, host_grid, this->nxyz);
    if (inverse)
    {
        this->fft_bundle.fft3D_backward<FPTYPE>(gpu_ctx, device_grid, device_grid);
    }
    else
    {
        this->fft_bundle.fft3D_forward<FPTYPE>(gpu_ctx, device_grid, device_grid);
    }
    base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(
        cpu_ctx, gpu_ctx, host_grid, device_grid, this->nxyz);
}
#endif

template <typename FPTYPE>
void PW_Basis::fft_forward() const
{
    auto* auxr = this->fft_bundle.get_auxr_data<FPTYPE>();
    auto* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == "gpu" && this->poolnproc == 1 && !this->gamma_only)
    {
        // A complete complex grid can be transformed before gathering its sticks.
        this->fft3d_gpu<FPTYPE>(false);
        this->gatherp_scatters(auxr, auxg);
        return;
    }
#endif
    if (this->gamma_only)
    {
        this->fft_bundle.fftxyr2c(this->fft_bundle.get_rspace_data<FPTYPE>(), auxr);
    }
    else
    {
        this->fft_bundle.fftxyfor(auxr, auxr);
    }
    this->gatherp_scatters(auxr, auxg);
    this->fft_bundle.fftzfor(auxg, auxg);
}

template <typename FPTYPE>
void PW_Basis::fft_backward() const
{
    auto* auxr = this->fft_bundle.get_auxr_data<FPTYPE>();
    auto* auxg = this->fft_bundle.get_auxg_data<FPTYPE>();
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == "gpu" && this->poolnproc == 1 && !this->gamma_only)
    {
        // Assemble the reciprocal grid before the full inverse transform.
        this->gathers_scatterp(auxg, auxr);
        this->fft3d_gpu<FPTYPE>(true);
        return;
    }
#endif
    this->fft_bundle.fftzbac(auxg, auxg);
    this->gathers_scatterp(auxg, auxr);
    if (this->gamma_only)
    {
        this->fft_bundle.fftxyc2r(auxr, this->fft_bundle.get_rspace_data<FPTYPE>());
    }
    else
    {
        this->fft_bundle.fftxybac(auxr, auxr);
    }
}

/**
 * @brief transform real space to reciprocal space
 * @details c(g)=\int dr*f(r)*exp(-ig*r)
 *          Here we calculate c(g)
 * @param in: (nplane,ny,nx), complex<double> data
 * @param out: (nz, ns),  complex<double> data
 */
template <typename FPTYPE>
void PW_Basis::real2recip(const std::complex<FPTYPE>* in,
                          std::complex<FPTYPE>* out,
                          const bool add,
                          const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "real2recip");
#if defined(__CUDA)
    if (this->d_ig2fftbox != nullptr)
    {
        this->real2recip_gpu(in, out, add, factor);
        ModuleBase::timer::tick(this->classname, "real2recip");
        return;
    }
#endif

    assert(this->gamma_only == false);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->fft_bundle.get_auxr_data<FPTYPE>()[ir] = in[ir];
    }
    this->fft_forward<FPTYPE>();

    if(add)
    {
        FPTYPE tmpfac = factor / FPTYPE(this->nxyz);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] += tmpfac * this->fft_bundle.get_auxg_data<FPTYPE>()[this->ig2isz[ig]];
        }
    }
    else
    {
        FPTYPE tmpfac = 1.0 / FPTYPE(this->nxyz);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] = tmpfac * this->fft_bundle.get_auxg_data<FPTYPE>()[this->ig2isz[ig]];
        }
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
}

/**
 * @brief transform real space to reciprocal space
 * @details c(g)=\int dr*f(r)*exp(-ig*r)
 *          Here we calculate c(g)
 * @param in: (nplane,ny,nx), double data
 * @param out: (nz, ns),  complex<double> data
 */
template <typename FPTYPE>
void PW_Basis::real2recip(const FPTYPE* in, std::complex<FPTYPE>* out, const bool add, const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "real2recip");
#if defined(__CUDA)
    if (this->d_ig2fftbox != nullptr)
    {
        this->real2recip_gpu(in, out, add, factor);
        ModuleBase::timer::tick(this->classname, "real2recip");
        return;
    }
#endif
    if(this->gamma_only)
    {
        const int npy = this->ny * this->nplane;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ix = 0 ; ix < this->nx ; ++ix)
        {
            for(int ipy = 0 ; ipy < npy ; ++ipy)
            {
                this->fft_bundle.get_rspace_data<FPTYPE>()[ix*npy + ipy] = in[ix*npy + ipy];
            }
        }

    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ir = 0 ; ir < this->nrxx ; ++ir)
        {
            this->fft_bundle.get_auxr_data<FPTYPE>()[ir] = std::complex<FPTYPE>(in[ir],0);
        }
    }
    this->fft_forward<FPTYPE>();

    if(add)
    {
        FPTYPE tmpfac = factor / FPTYPE(this->nxyz);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] += tmpfac * this->fft_bundle.get_auxg_data<FPTYPE>()[this->ig2isz[ig]];
        }
    }
    else
    {
        FPTYPE tmpfac = 1.0 / FPTYPE(this->nxyz);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] = tmpfac * this->fft_bundle.get_auxg_data<FPTYPE>()[this->ig2isz[ig]];
        }
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
}

/**
 * @brief transform reciprocal space to real space
 * @details f(r)=1/V * \sum_{g} c(g)*exp(ig*r)
 *          Here we calculate f(r)
 * @param in: (nz,ns), complex<double>
 * @param out: (nplane, ny, nx), complex<double>
 */
template <typename FPTYPE>
void PW_Basis::recip2real(const std::complex<FPTYPE>* in,
                          std::complex<FPTYPE>* out,
                          const bool add,
                          const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "recip2real");
#if defined(__CUDA)
    if (this->d_ig2fftbox != nullptr)
    {
        this->recip2real_gpu(in, out, add, factor);
        ModuleBase::timer::tick(this->classname, "recip2real");
        return;
    }
#endif
    assert(this->gamma_only == false);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
    for(int i = 0 ; i < this->nst * this->nz ; ++i)
    {
        fft_bundle.get_auxg_data<FPTYPE>()[i] = std::complex<FPTYPE>(0, 0);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        this->fft_bundle.get_auxg_data<FPTYPE>()[this->ig2isz[ig]] = in[ig];
    }
    this->fft_backward<FPTYPE>();
    
    if(add)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ir = 0 ; ir < this->nrxx ; ++ir)
        {
            out[ir] += factor * this->fft_bundle.get_auxr_data<FPTYPE>()[ir];
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for(int ir = 0 ; ir < this->nrxx ; ++ir)
        {
            out[ir] = this->fft_bundle.get_auxr_data<FPTYPE>()[ir];
        }
    }
    ModuleBase::timer::tick(this->classname, "recip2real");
}

/**
 * @brief transform reciprocal space to real space
 * @details f(r)=1/V * \sum_{g} c(g)*exp(ig*r)
 *          Here we calculate f(r)
 * @param in: (nz,ns), complex<double>
 * @param out: (nplane, ny, nx), double
 */
template <typename FPTYPE>
void PW_Basis::recip2real(const std::complex<FPTYPE>* in, FPTYPE* out, const bool add, const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "recip2real");
#if defined(__CUDA)
    if (this->d_ig2fftbox != nullptr)
    {
        this->recip2real_gpu(in, out, add, factor);
        ModuleBase::timer::tick(this->classname, "recip2real");
        return;
    }
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
    for(int i = 0 ; i < this->nst * this->nz ; ++i)
    {
        fft_bundle.get_auxg_data<FPTYPE>()[i] = std::complex<double>(0, 0);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        this->fft_bundle.get_auxg_data<FPTYPE>()[this->ig2isz[ig]] = in[ig];
    }
    this->fft_backward<FPTYPE>();

    if(this->gamma_only)
    {
        // r2c in place
        const int npy = this->ny * this->nplane;

        if(add)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 4096/sizeof(FPTYPE))
#endif
            for(int ix = 0 ; ix < this->nx ; ++ix)
            {
                for(int ipy = 0 ; ipy < npy ; ++ipy)
                {
                    out[ix*npy + ipy] += factor * this->fft_bundle.get_rspace_data<FPTYPE>()[ix*npy + ipy];
                }
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 4096/sizeof(FPTYPE))
#endif
            for(int ix = 0 ; ix < this->nx ; ++ix)
            {
                for(int ipy = 0 ; ipy < npy ; ++ipy)
                {
                    out[ix*npy + ipy] = this->fft_bundle.get_rspace_data<FPTYPE>()[ix*npy + ipy];
                }
            }
        }
    }
    else
    {
        if(add)
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
            for(int ir = 0 ; ir < this->nrxx ; ++ir)
            {
                out[ir] += factor * this->fft_bundle.get_auxr_data<FPTYPE>()[ir].real();
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
            for(int ir = 0 ; ir < this->nrxx ; ++ir)
            {
                out[ir] = this->fft_bundle.get_auxr_data<FPTYPE>()[ir].real();
            }
        }
    }
    ModuleBase::timer::tick(this->classname, "recip2real");
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
}
