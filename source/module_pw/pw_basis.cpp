#include "pw_basis.h"
#include "../module_base/mymath.h"
#include "../module_base/timer.h"
#include "../module_base/global_function.h"


namespace ModulePW
{
#ifdef __MPI
MPI_Datatype PW_Basis::mpi_dcomplex;
int PW_Basis::member = 0;
#endif
PW_Basis::PW_Basis()
{
    classname="PW_Basis";
#ifdef __MPI
    if(member == 0)
    {
        MPI_Datatype block[2];
	    block[0]=MPI_DOUBLE;
	    block[1]=MPI_DOUBLE;
	    int ac[2]={1,1};
	    MPI_Aint dipc[2]={0,sizeof(double)};
	    MPI_Type_create_struct(2,ac,dipc,block,&mpi_dcomplex);
	    MPI_Type_commit(&mpi_dcomplex);
        ++member;
    }
#endif
}

PW_Basis:: ~PW_Basis()
{
    delete[] ig2isz;
    delete[] istot2ixy;
    delete[] is2fftixy;
    delete[] fftixy2ip;
    delete[] nst_per;
    delete[] npw_per;
    delete[] gdirect;
    delete[] gcar;
    delete[] gg;
    delete[] startz;
    delete[] numz;
    delete[] numg;
    delete[] numr;
    delete[] startg;
    delete[] startr;
    delete[] ig2igg;
    delete[] gg_uniq;
#ifdef __MPI
    --member;
    if(member==0)
    {
        MPI_Type_free(&mpi_dcomplex);
    }
#endif
}

/// 
/// distribute plane wave basis and real-space grids to different processors
/// set up maps for fft and create arrays for MPI_Alltoall
/// set up ffts
///
void PW_Basis::setuptransform()
{
    ModuleBase::timer::tick(this->classname, "setuptransform");
    this->distribute_r();
    this->distribute_g();
    this->getstartgr();
    this->ft.clear();
    if(this->xprime)    this->ft.initfft(this->nx,this->ny,this->nz,this->lix,this->rix,this->nst,this->nplane,this->poolnproc,this->gamma_only, this->xprime);
    else                this->ft.initfft(this->nx,this->ny,this->nz,this->liy,this->riy,this->nst,this->nplane,this->poolnproc,this->gamma_only, this->xprime);
    this->ft.setupFFT();
    ModuleBase::timer::tick(this->classname, "setuptransform");
}

void PW_Basis::getstartgr()
{
    if(this->gamma_only)    this->nmaxgr = ( this->npw > (this->nrxx+1)/2 ) ? this->npw : (this->nrxx+1)/2;
    else                    this->nmaxgr = ( this->npw > this->nrxx ) ? this->npw : this->nrxx;
    this->nmaxgr = (this->nz * this->nst > this->nxy * nplane) ? this->nz * this->nst : this->nxy * nplane;
    
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
	for (int ip = 0;ip < poolnproc; ++ip) this->numg[ip] = this->nst_per[poolrank] * this->numz[ip];


	// Each processor in a pool send a piece of each stick(nst_per[ip]) to
	// other processors in this pool
	// rank_use processor receive datas in npps[rank_p] planes.
	for (int ip = 0;ip < poolnproc; ++ip) this->numr[ip] = this->nst_per[ip] * this->numz[poolrank];


	// startg record the starting 'numg' position in each processor.
	this->startg[0] = 0;
	for (int ip = 1;ip < poolnproc; ++ip) this->startg[ip] = this->startg[ip-1] + this->numg[ip-1];


	// startr record the starting 'numr' position
	this->startr[0] = 0;
	for (int ip = 1;ip < poolnproc; ++ip) this->startr[ip] = this->startr[ip-1] + this->numr[ip-1];
    return;
}

///
/// Collect planewaves on current core, and construct gg, gdirect, gcar according to ig2isz and is2fftixy.
/// known: ig2isz, is2fftixy
/// output: gg, gdirect, gcar
/// 
void PW_Basis::collect_local_pw()
{
   delete[] this->gg; this->gg = new double[this->npw];
   delete[] this->gdirect; this->gdirect = new ModuleBase::Vector3<double>[this->npw];
   delete[] this->gcar; this->gcar = new ModuleBase::Vector3<double>[this->npw];

    ModuleBase::Vector3<double> f;
    for(int ig = 0 ; ig < this-> npw ; ++ig)
    {
        int isz = this->ig2isz[ig];
        int iz = isz % this->nz;
        int is = isz / this->nz;
        int ixy = this->is2fftixy[is];
        int ix = ixy / this->fftny;
        int iy = ixy % this->fftny;
        if (ix >= int(this->nx/2) + 1) ix -= this->nx;
        if (iy >= int(this->ny/2) + 1) iy -= this->ny;
        if (iz >= int(this->nz/2) + 1) iz -= this->nz;
        f.x = ix;
        f.y = iy;
        f.z = iz;
        this->gg[ig] = f * (this->GGT * f);
        this->gdirect[ig] = f;
        this->gcar[ig] = f * this->G;
        if(this->gg[ig] < 1e-8) this->ig_gge0 = ig;
    }
    return;
}

///
/// Collect modulus of planewaves on current cores
/// known: ig2isz, is2fftixy
/// output: ig2igg, gg_uniq, ngg
/// 
void PW_Basis::collect_uniqgg()
{
    delete[] this->ig2igg; this->ig2igg = new int [this->npw];
    int *sortindex = new int [this->npw];
    double *tmpgg = new double [this->npw];
    double *tmpgg2 = new double [this->npw];
    ModuleBase::Vector3<double> f;
    for(int ig = 0 ; ig < this-> npw ; ++ig)
    {
        int isz = this->ig2isz[ig];
        int iz = isz % this->nz;
        int is = isz / this->nz;
        int ixy = this->is2fftixy[is];
        int ix = ixy / this->fftny;
        int iy = ixy % this->fftny;
        if (ix >= int(this->nx/2) + 1) ix -= this->nx;
        if (iy >= int(this->ny/2) + 1) iy -= this->ny;
        if (iz >= int(this->nz/2) + 1) iz -= this->nz;
        f.x = ix;
        f.y = iy;
        f.z = iz;
        tmpgg[ig] = f * (this->GGT * f);
        if(tmpgg[ig] < 1e-8) this->ig_gge0 = ig;
    }

    ModuleBase::GlobalFunc::ZEROS(sortindex, this->npw);
    ModuleBase::heapsort(this->npw, tmpgg, sortindex);
   

    int igg = 0;
    this->ig2igg[sortindex[0]] = 0;
    tmpgg2[0] = tmpgg[0];
    double avg_gg = tmpgg2[igg];
    int avg_n = 1;
    for (int ig = 1; ig < this->npw; ++ig)
    {
        if (std::abs(tmpgg[ig] - tmpgg2[igg]) > 1.0e-8)
        {
            tmpgg2[igg] = avg_gg / double(avg_n) ;
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
        if(ig == this->npw)
        {
            tmpgg2[igg] = avg_gg / double(avg_n) ;
        }
    }
    this->ngg = igg + 1;
    delete[] this->gg_uniq; this->gg_uniq = new double [this->ngg];
    for(int igg = 0 ; igg < this->ngg ; ++igg)
    {
            gg_uniq[igg] = tmpgg2[igg];
    }
    delete[] sortindex;
    delete[] tmpgg;
    delete[] tmpgg2;
}

void PW_Basis::getfftixy2is(int * fftixy2is)
{
//Note: please assert when is1 >= is2, fftixy2is[is1] >= fftixy2is[is2]!
    for(int ixy = 0 ; ixy < this->fftnxy ; ++ixy)   fftixy2is[ixy] = -1;
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


}