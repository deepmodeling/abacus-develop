#include "fft.h"
#include <complex>
#include "pw_basis.h"
#include <cassert>
#include "../module_base/global_function.h"
#include "../module_base/timer.h"
#include "pw_gatherscatter.h"

namespace ModulePW
{

/// 
/// transform real space to reciprocal space
/// c(k,g)=\int dr*f(r)*exp(-ig*r)
/// c(k,g)=c_k(g)*exp(ik*r)
/// c_k(g)=\int dr*f(r)*exp(-i(g+k)*r)
/// Here we calculate c(k,g)
/// in: (nplane,ny,nx), complex<double> data
/// out: (nz, ns),  complex<double> data
/// 
void PW_Basis:: real2recip(const std::complex<double> * in, std::complex<double> * out, const bool add, const double factor)
{
    ModuleBase::timer::tick(this->classname, "real2recip");

    assert(this->gamma_only == false);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.auxr[ir] = in[ir];
    }
    this->ft.fftxyfor(ft.auxr,ft.auxr);

    this->gatherp_scatters(this->ft.auxr, this->ft.auxg);
    
    this->ft.fftzfor(ft.auxg,ft.auxg);

    if(add)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] += factor / double(this->nxyz) * this->ft.auxg[this->ig2isz[ig]];
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] = this->ft.auxg[this->ig2isz[ig]] / double(this->nxyz);
        }
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
    return;
}

///
/// transform real space to reciprocal space
/// in: (nplane,ny,nx), double data
/// out: (nz, ns), complex<double> data
///
void PW_Basis:: real2recip(const double * in, std::complex<double> * out, const bool add, const double factor)
{
    ModuleBase::timer::tick(this->classname, "real2recip");
    if(this->gamma_only)
    {
        const int npy = this->ny * this->nplane;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
        for(int ix = 0 ; ix < this->nx ; ++ix)
        {
            for(int ipy = 0 ; ipy < npy ; ++ipy)
            {
                this->ft.r_rspace[ix*npy + ipy] = in[ix*npy + ipy];
            }
        }

        this->ft.fftxyr2c(ft.r_rspace,ft.auxr);
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir = 0 ; ir < this->nrxx ; ++ir)
        {
            this->ft.auxr[ir] = std::complex<double>(in[ir],0);
        }
        this->ft.fftxyfor(ft.auxr,ft.auxr);
    }
    this->gatherp_scatters(this->ft.auxr, this->ft.auxg);
    
    this->ft.fftzfor(ft.auxg,ft.auxg);

    if(add)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] += factor / double(this->nxyz) * this->ft.auxg[this->ig2isz[ig]];
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] = this->ft.auxg[this->ig2isz[ig]] / double(this->nxyz);
        }
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
    return;
}

/// 
/// transform reciprocal space to real space
/// f(r)=1/V * \sum_{g} c(k,g)*exp(ig*r)
/// c(k,g)=c_k(g)*exp(ik*r)
/// f(r)=1/V * \sum_{g} c_k(g)*exp(i(g+k)*r)
/// Here we use c(k,g)
/// in: (nz,ns), complex<double>
/// out: (nplane, ny, nx), complex<double>
/// 
void PW_Basis:: recip2real(const std::complex<double> * in, std::complex<double> * out, const bool add, const double factor)
{
    ModuleBase::timer::tick(this->classname, "recip2real");
    assert(this->gamma_only == false);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int i = 0 ; i < this->nst * this->nz ; ++i)
    {
        ft.auxg[i] = std::complex<double>(0, 0);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        this->ft.auxg[this->ig2isz[ig]] = in[ig];
    }
    this->ft.fftzbac(ft.auxg, ft.auxg);

    this->gathers_scatterp(this->ft.auxg,this->ft.auxr);

    this->ft.fftxybac(ft.auxr,ft.auxr);
    
    if(add)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir = 0 ; ir < this->nrxx ; ++ir)
        {
            out[ir] += factor * this->ft.auxr[ir];
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir = 0 ; ir < this->nrxx ; ++ir)
        {
            out[ir] = this->ft.auxr[ir];
        }
    }
    ModuleBase::timer::tick(this->classname, "recip2real");

    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<double>
/// out: (nplane, ny, nx), double
///
void PW_Basis:: recip2real(const std::complex<double> * in, double * out, const bool add, const double factor)
{
    ModuleBase::timer::tick(this->classname, "recip2real");
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int i = 0 ; i < this->nst * this->nz ; ++i)
    {
        ft.auxg[i] = std::complex<double>(0, 0);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        this->ft.auxg[this->ig2isz[ig]] = in[ig];
    }
    this->ft.fftzbac(ft.auxg, ft.auxg);

    this->gathers_scatterp(this->ft.auxg, this->ft.auxr);

    if(this->gamma_only)
    {
        this->ft.fftxyc2r(ft.auxr,ft.r_rspace);

        // r2c in place
        const int npy = this->ny * this->nplane;

        if(add)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
            for(int ix = 0 ; ix < this->nx ; ++ix)
            {
                for(int ipy = 0 ; ipy < npy ; ++ipy)
                {
                    out[ix*npy + ipy] += factor * this->ft.r_rspace[ix*npy + ipy];
                }
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
            for(int ix = 0 ; ix < this->nx ; ++ix)
            {
                for(int ipy = 0 ; ipy < npy ; ++ipy)
                {
                    out[ix*npy + ipy] = this->ft.r_rspace[ix*npy + ipy];
                }
            }
        }
    }
    else
    {
        this->ft.fftxybac(ft.auxr,ft.auxr);
        if(add)
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for(int ir = 0 ; ir < this->nrxx ; ++ir)
            {
                out[ir] += factor * this->ft.auxr[ir].real();
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for(int ir = 0 ; ir < this->nrxx ; ++ir)
            {
                out[ir] = this->ft.auxr[ir].real();
            }
        }
    }
    ModuleBase::timer::tick(this->classname, "recip2real");
    return;
}

#ifdef __MIX_PRECISION
///
/// transform real space to reciprocal space
/// in: (nplane,ny,nx), complex<float> data
/// out: (nz, ns),  complex<float> data
///
void PW_Basis:: real2recip(const std::complex<float> * in, std::complex<float> * out, const bool add, const float factor)
{
    ModuleBase::timer::tick(this->classname, "real2recip");
    assert(this->gamma_only == false);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int ir = 0 ; ir < this->nrxx ; ++ir)
    {
        this->ft.auxfr[ir] = in[ir];
    }
    this->ft.fftfxyfor(ft.auxfr,ft.auxfr);

    this->gatherp_scatters(this->ft.auxfr, this->ft.auxfg);
    
    this->ft.fftfzfor(ft.auxfg,ft.auxfg);

    if(add)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] += factor / float(this->nxyz) * this->ft.auxfg[this->ig2isz[ig]];
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] = this->ft.auxfg[this->ig2isz[ig]] / float(this->nxyz);
        }
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
    return;
}

///
/// transform real space to reciprocal space
/// in: (nplane,ny,nx), float data
/// out: (nz, ns), complex<float> data
///
void PW_Basis:: real2recip(const float * in, std::complex<float> * out, const bool add, const float factor)
{
    ModuleBase::timer::tick(this->classname, "real2recip");
    if(this->gamma_only)
    {
        const int npy = this->ny * this->nplane;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
        for(int ix = 0 ; ix < this->nx ; ++ix)
        {
            for(int ipy = 0 ; ipy < npy ; ++ipy)
            {
                this->ft.rf_rspace[ix*npy + ipy] = in[ix*npy + ipy];
            }
        }

        this->ft.fftfxyr2c(ft.rf_rspace,ft.auxfr);
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir = 0 ; ir < this->nrxx ; ++ir)
        {
            this->ft.auxfr[ir] = std::complex<float>(in[ir], 0);
        }
        this->ft.fftfxyfor(ft.auxfr,ft.auxfr);
    }

    this->gatherp_scatters(this->ft.auxfr, this->ft.auxfg);
    
    this->ft.fftfzfor(ft.auxfg,ft.auxfg);

    if(add)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] += factor / float(this->nxyz) * this->ft.auxfg[this->ig2isz[ig]];
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig = 0 ; ig < this->npw ; ++ig)
        {
            out[ig] = this->ft.auxfg[this->ig2isz[ig]] / float(this->nxyz);
        }
    }
    ModuleBase::timer::tick(this->classname, "real2recip");
    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<float>
/// out: (nplane, ny,nx), complex<float>
///
void PW_Basis:: recip2real(const std::complex<float> * in, std::complex<float> * out, const bool add, const float factor)
{
    ModuleBase::timer::tick(this->classname, "recip2real");
    assert(this->gamma_only == false);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int i = 0 ; i < this->nst * this->nz ; ++i)
    {
        ft.auxfg[i] = std::complex<float>(0, 0);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        this->ft.auxfg[this->ig2isz[ig]] = in[ig];
    }
    this->ft.fftfzbac(ft.auxfg, ft.auxfg);

    this->gathers_scatterp(this->ft.auxfg,this->ft.auxfr);

    this->ft.fftfxybac(ft.auxfr,ft.auxfr);
    
    if(add)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir = 0 ; ir < this->nrxx ; ++ir)
        {
            out[ir] += factor * this->ft.auxfr[ir];
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir = 0 ; ir < this->nrxx ; ++ir)
        {
            out[ir] = this->ft.auxfr[ir];
        }
    }
    ModuleBase::timer::tick(this->classname, "recip2real");
    return;
}

///
/// transform reciprocal space to real space
/// in: (nz,ns), complex<float>
/// out: (nplane, ny,nx), float
///
void PW_Basis:: recip2real(const std::complex<float> * in, float * out, const bool add, const float factor)
{
    ModuleBase::timer::tick(this->classname, "recip2real");
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int i = 0 ; i < this->nst * this->nz ; ++i)
    {
        ft.afuxg[i] = std::complex<float>(0, 0);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int ig = 0 ; ig < this->npw ; ++ig)
    {
        this->ft.auxfg[this->ig2isz[ig]] = in[ig];
    }
    this->ft.fftfzbac(ft.auxfg, ft.auxfg);
    
    this->gathers_scatterp(this->ft.auxfg, this->ft.auxfr);

    if(this->gamma_only)
    {
        this->ft.fftfxyc2r(ft.auxfr,ft.rf_rspace);

        const int npy = this->ny * this->nplane;
        if(add)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
            for(int ix = 0 ; ix < this->nx ; ++ix)
            {
                for(int ipy = 0 ; ipy < npy ; ++ipy)
                {
                    out[ix*npy + ipy] += factor * this->ft.rf_rspace[ix*npy + ipy];
                }
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
            for(int ix = 0 ; ix < this->nx ; ++ix)
            {
                for(int ipy = 0 ; ipy < npy ; ++ipy)
                {
                    out[ix*npy + ipy] = this->ft.rf_rspace[ix*npy + ipy];
                }
            }
        }
    }
    else
    {
        this->ft.fftfxybac(ft.auxfr,ft.auxfr);

        if(add)
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for(int ir = 0 ; ir < this->nrxx ; ++ir)
            {
                out[ir] += factor * this->ft.auxfr[ir].real();
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for(int ir = 0 ; ir < this->nrxx ; ++ir)
            {
                out[ir] = this->ft.auxfr[ir].real();
            }
        }
    }
    ModuleBase::timer::tick(this->classname, "recip2real");
    return;
}

#endif
}