#include <complex>
#include "pw_basis_k.h"
#include <cassert>
#include "../module_base/timer.h"
#include "pw_gatherscatter.h"
#include "module_pw/include/pw_multi_device.h"

namespace ModulePW
{

///
/// transform real space to reciprocal space
/// in: (nplane, ny, nx), complex<double> data
/// out: (nz, ns),  complex<double> data
///
void PW_Basis_K:: real2recip(const std::complex<double> * in, std::complex<double> * out, const int ik, const bool add, const double factor)
{
    ModuleBase::timer::tick(this->classname, "real2recip");

    assert(this->gamma_only == false);
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.auxr[ir] = in[ir];
    }
    this->ft.fftxyfor(ft.auxr,ft.auxr);

    this->gatherp_scatters(this->ft.auxr, this->ft.auxg);
    
    this->ft.fftzfor(ft.auxg,ft.auxg);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    if(add)
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] += factor / double(this->nxyz) * this->ft.auxg[this->igl2isz_k[igl+startig]];
    }
    else
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] = this->ft.auxg[this->igl2isz_k[igl+startig]] / double(this->nxyz);
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
    return;
}

///
/// transform real space to reciprocal space
/// in: (nplane, ny, nx), double data
/// out: (nz, ns), complex<double> data
///
void PW_Basis_K:: real2recip(const double * in, std::complex<double> * out, const int ik, const bool add, const double factor)
{
    ModuleBase::timer::tick(this->classname, "real2recip");
    assert(this->gamma_only == true);
    // for(int ir = 0 ; ir < this->nrxx ; ++ir)
    // {
    //     this->ft.r_rspace[ir] = in[ir];
    // }
    // r2c in place
    const int npy = this->ny * this->nplane;
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            this->ft.r_rspace[ixpy + ipy] = in[ixpy + ipy];
        }
    }

    this->ft.fftxyr2c(ft.r_rspace,ft.auxr);

    this->gatherp_scatters(this->ft.auxr, this->ft.auxg);
    
    this->ft.fftzfor(ft.auxg,ft.auxg);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    if(add)
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] += factor / double(this->nxyz) * this->ft.auxg[this->igl2isz_k[igl+startig]];
    }
    else
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] = this->ft.auxg[this->igl2isz_k[igl+startig]] / double(this->nxyz);
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<double>
/// out: (nplane, ny, nx), complex<double>
///
void PW_Basis_K:: recip2real(const std::complex<double> * in, std::complex<double> * out, const int ik, const bool add, const double factor)
{
    ModuleBase::timer::tick(this->classname, "recip2real");
    assert(this->gamma_only == false);
    ModuleBase::GlobalFunc::ZEROS(ft.auxg, this->nst * this->nz);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        this->ft.auxg[this->igl2isz_k[igl+startig]] = in[igl];
    }
    this->ft.fftzbac(ft.auxg, ft.auxg);

    this->gathers_scatterp(this->ft.auxg,this->ft.auxr);

    this->ft.fftxybac(ft.auxr,ft.auxr);
    
    if(add)
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] += factor * this->ft.auxr[ir];
    }
    else
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] = this->ft.auxr[ir];
    }
    ModuleBase::timer::tick(this->classname, "recip2real");

    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<double>
/// out: (nplane, ny, nx), double
///
void PW_Basis_K:: recip2real(const std::complex<double> * in, double * out, const int ik, const bool add, const double factor)
{
    ModuleBase::timer::tick(this->classname, "recip2real");
    assert(this->gamma_only == true);
    ModuleBase::GlobalFunc::ZEROS(ft.auxg, this->nst * this->nz);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        this->ft.auxg[this->igl2isz_k[igl+startig]] = in[igl];
    }
   this->ft.fftzbac(ft.auxg, ft.auxg);
    
    this->gathers_scatterp(this->ft.auxg, this->ft.auxr);

    this->ft.fftxyc2r(ft.auxr,ft.r_rspace);

    // for(int ir = 0 ; ir < this->nrxx ; ++ir)
    // {
    //     out[ir] = this->ft.r_rspace[ir] / this->nxyz;
    // }

    // r2c in place
    const int npy = this->ny * this->nplane;
    if(add)
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            out[ixpy + ipy] += factor * this->ft.r_rspace[ixpy + ipy];
        }
    }
    else
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            out[ixpy + ipy] = this->ft.r_rspace[ixpy + ipy];
        }
    }
    ModuleBase::timer::tick(this->classname, "recip2real");
    return;
}

void PW_Basis_K::real_to_recip(const psi::DEVICE_CPU * /*dev*/, const std::complex<double> * in, std::complex<double> * out, const int ik, const bool add, const double factor)
{
    this->real2recip(in, out, ik, add, factor);
}

void PW_Basis_K::recip_to_real(const psi::DEVICE_CPU * /*dev*/, const std::complex<double> * in, std::complex<double> * out, const int ik, const bool add, const double factor)
{
    this->recip2real(in, out, ik, add, factor);
}

#if defined(__CUDA) || defined(__UT_USE_CUDA)
void PW_Basis_K::real_to_recip(const psi::DEVICE_GPU * ctx, const std::complex<double> * in, std::complex<double> * out, const int ik, const bool add, const double factor)
{
    ModuleBase::timer::tick(this->classname, "real_to_recip gpu");
    assert(this->gamma_only == false);
    assert(this->poolnproc == 1);

    // for(int ir = 0; ir < this->nrxx ; ++ir) {
    //     this->ft.auxr[ir] = in[ir];
    // }
    // psi::DEVICE_CPU *cpu_ctx = {};
     psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_GPU>()(
         ctx, ctx,
         this->ft.auxr_3d, in,
         this->nrxx);
    // psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_GPU>()(
    //     ctx, ctx,
    //     this->ft.auxr_3d, in,
    //     this->nrxx);
    // cudaMemcpy(this->ft.auxr_3d, this->ft.auxr, sizeof(std::complex<double>) * this->nxyz, cudaMemcpyHostToDevice);
    this->ft.fft3D_forward(this->ft.auxr_3d, this->ft.auxr_3d);
    // cudaMemcpy(this->ft.auxr, this->ft.auxr_3d, sizeof(std::complex<double>) * this->nxyz, cudaMemcpyDeviceToHost);

    const int startig = ik*this->npwk_max;
    const int npw_k = this->npwk[ik];
    // if(add) {
    //     for(int ig = 0; ig < npw_k; ++ig) {
    //         out[ig] += factor / static_cast<double>(this->nxyz) * this->ft.auxr[this->ig2ixyz_k[ig + startig]];
    //     }
    // }
    // else {
    //     for(int ig = 0; ig < npw_k; ++ig) {
    //         out[ig] = this->ft.auxr[this->ig2ixyz_k[ig + startig]] / static_cast<double>(this->nxyz);
    //     }
    // }
    // std::complex<double> * d_out = nullptr;
    // cudaMalloc(reinterpret_cast<void**>(&d_out), sizeof(std::complex<double>) * this->nxyz);
    // cudaMemcpy(d_out, out, sizeof(std::complex<double>) * this->nxyz, cudaMemcpyHostToDevice);
    set_real_to_recip_output_op<double, psi::DEVICE_GPU>()(
        ctx, npw_k, this->nxyz, add, factor,  this->ig2ixyz_k + startig, this->ft.auxr_3d, out);
    // cudaMemcpy(out, d_out, sizeof(std::complex<double>) * this->nxyz, cudaMemcpyDeviceToHost);
    // cudaFree(d_out);
    ModuleBase::timer::tick(this->classname, "real_to_recip gpu");
}

void PW_Basis_K::recip_to_real(const psi::DEVICE_GPU * ctx, const std::complex<double> * in, std::complex<double> * out, const int ik, const bool add, const double factor)
{
    ModuleBase::timer::tick(this->classname, "recip_to_real gpu");
    assert(this->gamma_only == false);
    assert(this->poolnproc == 1);
    // ModuleBase::GlobalFunc::ZEROS(ft.auxr_3d, this->nxyz);
    psi::memory::set_memory_op<std::complex<double>, psi::DEVICE_GPU>()(
        ctx, this->ft.auxr_3d, 0, this->nxyz);

    const int startig = ik*this->npwk_max;
    const int npw_k = this->npwk[ik];

    // for(int ig = 0; ig < npw_k; ++ig)
    // {
    //     this->ft.auxr[this->ig2ixyz_k[ig + startig]] = in[ig];
    // }
    // std::complex<double> * d_in = nullptr;
    // cudaMalloc(reinterpret_cast<void**>(&d_in), sizeof(std::complex<double>) * this->nxyz);
    // cudaMemcpy(d_in, in, sizeof(std::complex<double>) * this->nxyz, cudaMemcpyHostToDevice);
    // set_3d_fft_box_op<double, psi::DEVICE_GPU>()(
    //     ctx, npw_k, d_in, this->ft.auxr_3d, this->ig2ixyz_k + startig);
    set_3d_fft_box_op<double, psi::DEVICE_GPU>()(
        ctx, npw_k, this->ig2ixyz_k + startig, in, this->ft.auxr_3d);
    //auxg should be "auxg = new complex<double>[nxyz]“
    // cudaMemcpy(this->ft.auxr_3d, this->ft.auxr, sizeof(std::complex<double>) * this->nxyz, cudaMemcpyHostToDevice);
    this->ft.fft3D_backward(this->ft.auxr_3d, this->ft.auxr_3d);
    // cudaMemcpy(this->ft.auxr, this->ft.auxr_3d, sizeof(std::complex<double>) * this->nxyz, cudaMemcpyDeviceToHost);

    // if(add) {
    //     for(int ir = 0 ; ir < this->nrxx ; ++ir) {
    //         out[ir] += factor * this->ft.auxr[ir];
    //     }
    // }
    // else {
    //     for(int ir = 0; ir < this->nrxx ; ++ir) {
    //         out[ir] = this->ft.auxr[ir];
    //     }
    // }
    // std::complex<double> * d_out = nullptr;
    // cudaMalloc(reinterpret_cast<void**>(&d_out), sizeof(std::complex<double>) * this->nxyz);
    // cudaMemcpy(d_out, out, sizeof(std::complex<double>) * this->nxyz, cudaMemcpyHostToDevice);
    set_recip_to_real_output_op<double, psi::DEVICE_GPU>()(
        ctx, this->nrxx, add, factor, this->ft.auxr_3d, out);
    // cudaMemcpy(out, d_out, sizeof(std::complex<double>) * this->nxyz, cudaMemcpyDeviceToHost);
    // cudaFree(d_in);
    // cudaFree(d_out);

    ModuleBase::timer::tick(this->classname, "recip_to_real gpu");
}
#endif //defined(__CUDA) || defined(__UT_USE_CUDA)

#ifdef __MIX_PRECISION
///
/// transform real space to reciprocal space
/// in: (nplane, ny, nx), complex<float> data
/// out: (nz, ns),  complex<float> data
///
void PW_Basis_K:: real2recip(const std::complex<float> * in, std::complex<float> * out, const int ik, const bool add, const float factor)
{
    ModuleBase::timer::tick(this->classname, "real2recip");
    assert(this->gamma_only == false);
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.auxfr[ir] = in[ir];
    }
    this->ft.fftfxyfor(ft.auxfr,ft.auxfr);

    this->gatherp_scatters(this->ft.auxfr, this->ft.auxfg);
    
    this->ft.fftfzfor(ft.auxfg,ft.auxfg);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    if(add)
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] += factor / float(this->nxyz) * this->ft.auxfg[this->igl2isz_k[igl+startig]];
    }
    else
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] = this->ft.auxfg[this->igl2isz_k[igl+startig]] / float(this->nxyz);
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
    return;
}

///
/// transform real space to reciprocal space
/// in: (nplane, ny, nx), float data
/// out: (nz, ns), complex<float> data
///
void PW_Basis_K:: real2recip(const float * in, std::complex<float> * out, const int ik, const bool add, const float factor)
{
    ModuleBase::timer::tick(this->classname, "real2recip");
    assert(this->gamma_only == true);
    const int npy = this->ny * this->nplane;
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            this->ft.rf_rspace[ixpy + ipy] = in[ixpy + ipy];
        }
    }

    this->ft.fftfxyr2c(ft.rf_rspace,ft.auxfr);

    this->gatherp_scatters(this->ft.auxfr, this->ft.auxfg);
    
    this->ft.fftfzfor(ft.auxfg,ft.auxfg);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    if(add)
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] += factor / float(this->nxyz) * this->ft.auxfg[this->igl2isz_k[igl+startig]];
    }
    else
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        out[igl] = this->ft.auxfg[this->igl2isz_k[igl+startig]] / float(this->nxyz);
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<float>
/// out: (nplane, ny, nx), complex<float>
///
void PW_Basis_K:: recip2real(const std::complex<float> * in, std::complex<float> * out, const int ik, const bool add, const float factor)
{
    ModuleBase::timer::tick(this->classname, "recip2real");
    assert(this->gamma_only == false);
    ModuleBase::GlobalFunc::ZEROS(ft.auxfg, this->nst * this->nz);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        this->ft.auxfg[this->igl2isz_k[igl+startig]] = in[igl];
    }
    this->ft.fftfzbac(ft.auxfg, ft.auxfg);

    this->gathers_scatterp(this->ft.auxfg,this->ft.auxfr);

    this->ft.fftfxybac(ft.auxfr,ft.auxfr);
    
    if(add)
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] += factor * this->ft.auxfr[ir];
    }
    else
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        out[ir] = this->ft.auxfr[ir];
    }
    ModuleBase::timer::tick(this->classname, "recip2real");
    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<float>
/// out: (nplane, ny, nx), float
///
void PW_Basis_K:: recip2real(const std::complex<float> * in, float * out, const int ik, const bool add, const float factor)
{
    ModuleBase::timer::tick(this->classname, "recip2real");
    assert(this->gamma_only == true);
    ModuleBase::GlobalFunc::ZEROS(ft.auxfg, this->nst * this->nz);

    const int startig = ik*this->npwk_max;
    const int npwk = this->npwk[ik];
    for(int igl = 0 ; igl < npwk ; ++igl)
    {
        this->ft.auxfg[this->igl2isz_k[igl+startig]] = in[igl];
    }
   this->ft.fftfzbac(ft.auxfg, ft.auxfg);
    
    this->gathers_scatterp(this->ft.auxfg, this->ft.auxfr);

    this->ft.fftfxyc2r(ft.auxfr,ft.rf_rspace);

    const int npy = this->ny * this->nplane;
    if(add)
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            out[ixpy + ipy] += factor * this->ft.rf_rspace[ixpy + ipy];
        }
    }
    else
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        const int ixpy = ix*npy;
        for(int ipy = 0 ; ipy < npy ; ++ipy)
        {
            out[ixpy + ipy] = this->ft.rf_rspace[ixpy + ipy];
        }
    }
    ModuleBase::timer::tick(this->classname, "recip2real");
    return;
}

#endif
}