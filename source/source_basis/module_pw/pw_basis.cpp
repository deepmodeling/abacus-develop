#include "pw_basis.h"

#include <utility>
#include "source_base/mymath.h"
#include "source_base/timer.h"
#include "source_base/global_function.h"

#include <vector>

namespace ModulePW
{
PW_Basis::PW_Basis()
{
    classname="PW_Basis";
}

PW_Basis::PW_Basis(std::string device_, std::string precision_) : device(std::move(device_)), precision(std::move(precision_)) {
    classname="PW_Basis";
    this->fft_bundle.setfft("cpu",this->precision);
    this->double_data_ = (this->precision == "double") || (this->precision == "mixing");
    this->float_data_ = (this->precision == "single")  || (this->precision == "mixing");
}

PW_Basis:: ~PW_Basis()
{
    delete[] ig2isz;
    delete[] istot2ixy;
    delete[] is2fftixy;
    delete[] fftixy2ip;
    delete[] nst_per;
    delete[] npw_per;
    delete[] startz;
    delete[] numz;
    delete[] numg;
    delete[] numr;
    delete[] startg;
    delete[] startr;
    this->clear_owned_cache();
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == "gpu")
    {
        delmem_int_op()(this->d_is2fftixy);
        delmem_int_op()(this->ig2ixyz_gpu);
    }
#endif
}

void PW_Basis::clear_owned_cache()
{
    std::lock_guard<std::mutex> guard(this->cache_mutex);
    this->invalidate_cache_unlocked();
}

PW_Basis::CacheStats PW_Basis::get_cache_stats() const
{
    std::lock_guard<std::mutex> guard(this->cache_mutex);
    return this->get_cache_stats_unlocked();
}

PW_Basis::CacheStats PW_Basis::get_cache_stats_unlocked() const
{
    CacheStats stats;
    stats.local_pw_hits = this->local_pw_cache_hits.load();
    stats.local_pw_misses = this->local_pw_cache_misses.load();
    stats.uniqgg_hits = this->uniqgg_cache_hits.load();
    stats.uniqgg_misses = this->uniqgg_cache_misses.load();
    const bool has_local_pw_cache = this->local_pw_cache_valid.load()
                                    && this->npw > 0
                                    && this->gg != nullptr
                                    && this->gdirect != nullptr
                                    && this->gcar != nullptr;
    const bool has_uniqgg_cache = this->uniqgg_cache_valid.load()
                                  && this->ngg > 0
                                  && this->ig2igg != nullptr
                                  && this->gg_uniq != nullptr;
    if (has_local_pw_cache)
    {
        stats.cache_bytes += sizeof(double) * this->npw;
        stats.cache_bytes += sizeof(ModuleBase::Vector3<double>) * this->npw * 2;
    }
    if (has_uniqgg_cache)
    {
        stats.cache_bytes += sizeof(int) * this->npw;
        stats.cache_bytes += sizeof(double) * this->ngg;
    }
    return stats;
}

void PW_Basis::reset_cache_stats()
{
    this->local_pw_cache_hits.store(0);
    this->local_pw_cache_misses.store(0);
    this->uniqgg_cache_hits.store(0);
    this->uniqgg_cache_misses.store(0);
}

/// 
/// distribute plane wave basis and real-space grids to different processors
/// set up maps for fft and create arrays for MPI_Alltoall
/// set up ffts
///
void PW_Basis::setuptransform()
{
    ModuleBase::timer::start(this->classname, "setuptransform");
    this->distribute_r();
    this->distribute_g();
    this->getstartgr();
    this->fft_bundle.clear();
    
    if(this->xprime)    
    {
        this->fft_bundle.initfft(this->nx,this->ny,this->nz,this->lix,this->rix,this->nst,this->nplane,this->poolnproc,this->gamma_only, this->xprime);
    }
    else                
    {
        this->fft_bundle.initfft(this->nx,this->ny,this->nz,this->liy,this->riy,this->nst,this->nplane,this->poolnproc,this->gamma_only, this->xprime);
    }
    this->fft_bundle.setupFFT();
    ModuleBase::timer::end(this->classname, "setuptransform");
}

void PW_Basis::getstartgr()
{
    if(this->gamma_only)    
    {
        this->nmaxgr = ( this->npw > (this->nrxx+1)/2 ) ? this->npw : (this->nrxx+1)/2;
    }
    else
    {
        this->nmaxgr = ( this->npw > this->nrxx ) ? this->npw : this->nrxx;
    }
    
    //---------------------------------------------
	// sum : starting plane of FFT box.
	//---------------------------------------------
    delete[] this->numg; this->numg = new int[poolnproc];
	delete[] this->startg; this->startg = new int[poolnproc];
	delete[] this->startr; this->startr = new int[poolnproc];
	delete[] this->numr; this->numr = new int[poolnproc];

	// Each processor has a set of full sticks,
	// 'rank_use' processor send a piece(npps[ip]) of these sticks(nst_per[rank_use])
	// to all the other processors in this pool
	for (int ip = 0;ip < poolnproc; ++ip)
    {
        this->numg[ip] = this->nst_per[poolrank] * this->numz[ip];
    }


	// Each processor in a pool send a piece of each stick(nst_per[ip]) to
	// other processors in this pool
	// rank_use processor receive datas in npps[rank_p] planes.
	for (int ip = 0;ip < poolnproc; ++ip)
    {
        this->numr[ip] = this->nst_per[ip] * this->numz[poolrank];
    }


	// startg record the starting 'numg' position in each processor.
	this->startg[0] = 0;
	for (int ip = 1;ip < poolnproc; ++ip)
    {
        this->startg[ip] = this->startg[ip-1] + this->numg[ip-1];
    }


	// startr record the starting 'numr' position
	this->startr[0] = 0;
	for (int ip = 1;ip < poolnproc; ++ip)
    {
        this->startr[ip] = this->startr[ip-1] + this->numr[ip-1];
    }
    return;
}

///
/// Collect planewaves on current core, and construct gg, gdirect, gcar according to ig2isz and is2fftixy.
/// known: ig2isz, is2fftixy
/// output: gg, gdirect, gcar
/// 
void PW_Basis::collect_local_pw()
{
    if(this->npw <= 0)
    {
        return;
    }
    std::lock_guard<std::mutex> guard(this->cache_mutex);
    if (this->local_pw_cache_valid.load())
    {
        this->local_pw_cache_hits.fetch_add(1);
        return;
    }
    this->local_pw_cache_misses.fetch_add(1);
    this->ig_gge0 = -1;
    this->gg_cache_storage.reset(new double[this->npw]);
    this->gdirect_cache_storage.reset(new ModuleBase::Vector3<double>[this->npw]);
    this->gcar_cache_storage.reset(new ModuleBase::Vector3<double>[this->npw]);
    this->gg = this->gg_cache_storage.get();
    this->gdirect = this->gdirect_cache_storage.get();
    this->gcar = this->gcar_cache_storage.get();
    // Unique-G data depends on gg, so rebuilding local G data invalidates it.
    this->uniqgg_cache_valid.store(false);
    this->ig2igg_cache_storage.reset();
    this->gg_uniq_cache_storage.reset();
    this->ig2igg = nullptr;
    this->gg_uniq = nullptr;
    this->ngg = 0;

    ModuleBase::Vector3<double> f;
    int gamma_num = 0;
    for(int ig = 0 ; ig < this-> npw ; ++ig)
    {
        int isz = this->ig2isz[ig];
        int iz = isz % this->nz;
        int is = isz / this->nz;
        int ixy = this->is2fftixy[is];
        int ix = ixy / this->fftny;
        int iy = ixy % this->fftny;
        if (ix >= int(this->nx/2) + 1)
        {
            ix -= this->nx;
        }
        if (iy >= int(this->ny/2) + 1)
        {
            iy -= this->ny;
        }
        if (iz >= int(this->nz/2) + 1)
        {
            iz -= this->nz;
        }
        f.x = ix;
        f.y = iy;
        f.z = iz;
        this->gg[ig] = f * (this->GGT * f);
        this->gdirect[ig] = f;
        this->gcar[ig] = f * this->G;
        if(this->gg[ig] < 1e-8)
        {
            this->ig_gge0 = ig;
            ++gamma_num;
            if (gamma_num > 1)
            {
                ModuleBase::WARNING_QUIT("PW_Basis::collect_local_pw", 
                                        "More than one gamma point found in the plane wave basis set.\n");
            }
        }
    }
    this->local_pw_cache_valid.store(true);
    return;
}

///
/// Collect modulus of planewaves on current cores
/// known: ig2isz, is2fftixy
/// output: ig2igg, gg_uniq, ngg
/// 
void PW_Basis::collect_uniqgg()
{
    if(this->npw <= 0)
    {
        return;
    }
    std::lock_guard<std::mutex> guard(this->cache_mutex);
    if (this->uniqgg_cache_valid.load())
    {
        this->uniqgg_cache_hits.fetch_add(1);
        return;
    }
    this->uniqgg_cache_misses.fetch_add(1);
    this->ig_gge0 = -1;
    this->ig2igg_cache_storage.reset(new int[this->npw]);
    this->ig2igg = this->ig2igg_cache_storage.get();
    
    std::vector<int> sortindex(this->npw); // Reconstruct the plane-wave index mapping after sorting by energy.
    std::vector<double> tmpgg(this->npw);
    std::vector<double> tmpgg2(this->npw);
    // Reuse gg when collect_local_pw has already built the same G^2 values.
    if (this->local_pw_cache_valid.load() && this->gg != nullptr)
    {
        for(int ig = 0 ; ig < this-> npw ; ++ig)
        {
            tmpgg[ig] = this->gg[ig];
            if(tmpgg[ig] < 1e-8)
            {
                this->ig_gge0 = ig;
            }
        }
    }
    else
    {
        ModuleBase::Vector3<double> f;
        for(int ig = 0 ; ig < this-> npw ; ++ig)
        {
            int isz = this->ig2isz[ig];
            int iz = isz % this->nz;
            int is = isz / this->nz;
            int ixy = this->is2fftixy[is];
            int ix = ixy / this->fftny;
            int iy = ixy % this->fftny;
            if (ix >= int(this->nx/2) + 1)
            {
                ix -= this->nx;
            }
            if (iy >= int(this->ny/2) + 1)
            {
                iy -= this->ny;
            }
            if (iz >= int(this->nz/2) + 1)
            {
                iz -= this->nz;
            }
            f.x = ix;
            f.y = iy;
            f.z = iz;
            tmpgg[ig] = f * (this->GGT * f);
            if(tmpgg[ig] < 1e-8)
            {
                this->ig_gge0 = ig;
            }
        }
    }

    ModuleBase::GlobalFunc::ZEROS(sortindex.data(), this->npw);
    ModuleBase::heapsort(this->npw, tmpgg.data(), sortindex.data());
   

    int igg = 0;
    this->ig2igg[sortindex[0]] = 0;
    tmpgg2[0] = tmpgg[0];
    double avg_gg = tmpgg2[igg];//For waves with similar energy,take the average
    int avg_n = 1;//The number of waves required to take the average
    for (int ig = 1; ig < this->npw; ++ig)
    {
        if (std::abs(tmpgg[ig] - tmpgg2[igg]) > 1.0e-8)
        {
            tmpgg2[igg] = avg_gg / double(avg_n);
            ++igg;
            tmpgg2[igg] = tmpgg[ig];
            avg_gg = tmpgg2[igg];
            avg_n = 1;   
        }
        else
        {
            avg_n++;
            avg_gg += tmpgg[ig];
        }
        this->ig2igg[sortindex[ig]] = igg;
    }
    tmpgg2[igg] = avg_gg / double(avg_n);
    this->ngg = igg + 1;
    this->gg_uniq_cache_storage.reset(new double[this->ngg]);
    this->gg_uniq = this->gg_uniq_cache_storage.get();
    for(int igg = 0 ; igg < this->ngg ; ++igg)
    {
            gg_uniq[igg] = tmpgg2[igg];
    }
    this->uniqgg_cache_valid.store(true);
}

void PW_Basis::getfftixy2is(int * fftixy2is) const
{
//Note: please assert when is1 >= is2, fftixy2is[is1] >= fftixy2is[is2]!
    for(int ixy = 0 ; ixy < this->fftnxy ; ++ixy)
    {
        fftixy2is[ixy] = -1;
    }
    int ixy = 0;
    for(int is = 0; is < this->nst; ++is)
    {
        for(; ixy < this->fftnxy ; ++ixy)
        {
            if(this->is2fftixy[is] == ixy)
            {
                fftixy2is[ixy] = is;
                ++ixy;
                break;
            }
        }
    }
}

void PW_Basis::set_device(std::string device_) {
    this->device = std::move(device_);
    this->invalidate_cache();
}

void PW_Basis::set_precision(std::string precision_) {
    this->precision = std::move(precision_);
    this->invalidate_cache();
}

}
