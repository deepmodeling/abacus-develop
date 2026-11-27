#include "pw_basis.h"
#include "source_base/global_function.h"
#include "source_base/timer.h"
#include <algorithm>
#include <typeinfo>

namespace ModulePW
{
namespace detail
{
template <typename T>
inline void copy_complex_buffer(const std::complex<T>* in, std::complex<T>* out, const int count)
{
    if (count <= 0)
    {
        return;
    }

    std::copy_n(in, count, out);
}

// Top-level transform copies own the OpenMP parallel region; gather/scatter
// loops call the non-parallel helper inside their existing parallel regions.
template <typename T>
inline void copy_complex_buffer_parallel(const std::complex<T>* in, std::complex<T>* out, const int count)
{
    constexpr int chunk_size = 1024;
    if (count <= chunk_size)
    {
        copy_complex_buffer(in, out, count);
        return;
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
    for (int offset = 0; offset < count; offset += chunk_size)
    {
        const int chunk_count = std::min(chunk_size, count - offset);
        std::copy_n(in + offset, chunk_count, out + offset);
    }
#else
    copy_complex_buffer(in, out, count);
#endif
}
} // namespace detail

/**
 * @brief gather planes and scatter sticks
 * @param in: (nplane,fftny,fftnx)
 * @param out: (nz,nst)
 * @note in and out should be in different places
 * @note in[] will be changed
 */
template <typename T>
void PW_Basis::gatherp_scatters(std::complex<T>* in, std::complex<T>* out) const
{
    
    if(this->poolnproc == 1) //In this case nst=nstot, nz = nplane, 
    {
        const int nst_ = this->nst;
        const int nz_ = this->nz;
        const int* istot2ixy_ = this->istot2ixy;
        ModuleBase::timer::start(this->classname, "gatherp_copy_serial");
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int is = 0 ; is < nst_ ; ++is)
        {
            int ixy = istot2ixy_[is];
            std::complex<T>* outp = &out[is*nz_];
            const std::complex<T>* inp = &in[ixy*nz_];
            detail::copy_complex_buffer(inp, outp, nz_);
        }
        ModuleBase::timer::end(this->classname, "gatherp_copy_serial");
        return;
    }


#ifdef __MPI
    //change (nplane fftnxy) to (nplane,nstot)
    // Hence, we can send them at one time.
    const int nstot_gps = this->nstot;
    const int nplane_gps = this->nplane;
    const int* istot2ixy_gps = this->istot2ixy;
    ModuleBase::timer::start(this->classname, "gatherp_copy_pack");
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int istot = 0; istot < nstot_gps; ++istot)
    {
        int ixy = istot2ixy_gps[istot];
        std::complex<T>* outp = &out[istot * nplane_gps];
        const std::complex<T>* inp = &in[ixy * nplane_gps];
        detail::copy_complex_buffer(inp, outp, nplane_gps);
    }
    ModuleBase::timer::end(this->classname, "gatherp_copy_pack");

    //exchange data
    //(nplane,nstot) to (numz[ip],ns, poolnproc)
    if(typeid(T) == typeid(double))
    {
        MPI_Alltoallv(out, numr, startr, MPI_DOUBLE_COMPLEX, in, numg, startg, MPI_DOUBLE_COMPLEX, this->pool_world);
    }
    else if(typeid(T) == typeid(float))
    {
        MPI_Alltoallv(out, numr, startr, MPI_COMPLEX, in, numg, startg, MPI_COMPLEX, this->pool_world);
    }
    else
    {
        ModuleBase::WARNING_QUIT("PW_Basis::gatherp_scatters", "Unsupported data type for MPI_Alltoallv");
    }

    // change (nz,ns) to (numz[ip],ns, poolnproc)
    const int poolnproc_gps = this->poolnproc;
    const int nst_gps = this->nst;
    const int nz_gps = this->nz;
    const int* numz_gps = this->numz;
    const int* startg_gps = this->startg;
    const int* startz_gps = this->startz;
    ModuleBase::timer::start(this->classname, "gatherp_copy_unpack");
#ifdef _OPENMP
    #pragma omp parallel for collapse(2)
#endif
    for (int ip = 0; ip < poolnproc_gps ;++ip)
    {
        for (int is = 0; is < nst_gps; ++is)
        {
            int nzip = numz_gps[ip];
            std::complex<T> *outp0 = &out[startz_gps[ip]];
            std::complex<T> *inp0 = &in[startg_gps[ip]];
            std::complex<T>* outp = &outp0[is * nz_gps];
            const std::complex<T>* inp = &inp0[is * nzip ];
            detail::copy_complex_buffer(inp, outp, nzip);
        }
    }
    ModuleBase::timer::end(this->classname, "gatherp_copy_unpack");
#endif
    return;
}

/**
 * @brief gather sticks and scatter planes
 * @param in: (nz,nst)
 * @param out: (nplane,fftny,fftnx)
 * @note in and out should be in different places
 * @note in[] will be changed
 */
template <typename T>
void PW_Basis::gathers_scatterp(std::complex<T>* in, std::complex<T>* out) const
{
    if(this->poolnproc == 1) //In this case nrxx=fftnx*fftny*nz, nst = nstot, 
    {
        const int nrxx_ = this->nrxx;
        const int nst_ = this->nst;
        const int nz_ = this->nz;
        const int* istot2ixy_ = this->istot2ixy;
        ModuleBase::timer::start(this->classname, "gathers_zero_serial");
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(int i = 0; i < nrxx_; ++i)
        {
            out[i] = std::complex<T>(0, 0);
        }
        ModuleBase::timer::end(this->classname, "gathers_zero_serial");

        ModuleBase::timer::start(this->classname, "gathers_copy_serial");
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int is = 0 ; is < nst_ ; ++is)
        {
            int ixy = istot2ixy_[is];
            std::complex<T>* outp = &out[ixy*nz_];
            const std::complex<T>* inp = &in[is*nz_];
            detail::copy_complex_buffer(inp, outp, nz_);
        }
        ModuleBase::timer::end(this->classname, "gathers_copy_serial");
        return;
    }


#ifdef __MPI
    // change (nz,ns) to (numz[ip],ns, poolnproc)
    // Hence, we can send them at one time. 
    const int poolnproc_ = this->poolnproc;
    const int nst_ = this->nst;
    const int nz_ = this->nz;
    const int* numz_ = this->numz;
    const int* startg_ = this->startg;
    const int* startz_ = this->startz;
    ModuleBase::timer::start(this->classname, "gathers_copy_pack");
#ifdef _OPENMP
    #pragma omp parallel for collapse(2)
#endif
    for (int ip = 0; ip < poolnproc_ ;++ip)
    {
        for (int is = 0; is < nst_; ++is)
        {
            int nzip = numz_[ip];
            std::complex<T> *outp0 = &out[startg_[ip]];
            std::complex<T> *inp0 = &in[startz_[ip]];
            std::complex<T>* outp = &outp0[is * nzip];
            const std::complex<T>* inp = &inp0[is * nz_ ];
            detail::copy_complex_buffer(inp, outp, nzip);
        }
    }
    ModuleBase::timer::end(this->classname, "gathers_copy_pack");

    //exchange data
    //(numz[ip],ns, poolnproc) to (nplane,nstot)
    if(typeid(T) == typeid(double))
    {
        MPI_Alltoallv(out, numg, startg, MPI_DOUBLE_COMPLEX, in, numr, startr, MPI_DOUBLE_COMPLEX, this->pool_world);
    }
    else if(typeid(T) == typeid(float))
    {
        MPI_Alltoallv(out, numg, startg, MPI_COMPLEX, in, numr, startr, MPI_COMPLEX, this->pool_world);
    }
    else
    {
        ModuleBase::WARNING_QUIT("PW_Basis::gathers_scatterp", "Unsupported data type for MPI_Alltoallv");
    }

    const int nrxx_gsp = this->nrxx;
    ModuleBase::timer::start(this->classname, "gathers_zero_mpi");
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for(int i = 0; i < nrxx_gsp; ++i)
    {
        out[i] = std::complex<T>(0, 0);
    }
    ModuleBase::timer::end(this->classname, "gathers_zero_mpi");
    //change (nplane,nstot) to (nplane fftnxy)
    const int nstot = this->nstot;
    const int nplane = this->nplane;
    const int* istot2ixy = this->istot2ixy;
    ModuleBase::timer::start(this->classname, "gathers_copy_unpack");
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int istot = 0;istot < nstot; ++istot)
    {
        int ixy = istot2ixy[istot];
        //int ixy = (ixy / fftny)*ny + ixy % fftny;
        std::complex<T>* outp = &out[ixy * nplane];
        const std::complex<T>* inp = &in[istot * nplane];
        detail::copy_complex_buffer(inp, outp, nplane);
    }
    ModuleBase::timer::end(this->classname, "gathers_copy_unpack");
#endif
    return;
}



}
