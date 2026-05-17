#include "pw_basis.h"
#include "source_base/global_function.h"
#include "source_base/timer.h"
#include <typeinfo>

namespace ModulePW
{
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
    //ModuleBase::timer::start(this->classname, "gatherp_scatters");
    
    if(this->poolnproc == 1) //In this case nst=nstot, nz = nplane, 
    {
        const int nst = this->nst;
        const int nz = this->nz;
        const int* istot2ixy = this->istot2ixy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int is = 0 ; is < nst ; ++is)
        {
            int ixy = istot2ixy[is];
            std::complex<T> *outp = &out[is*nz];
            std::complex<T> *inp = &in[ixy*nz];
            for(int iz = 0 ; iz < nz ; ++iz)
            {
                outp[iz] = inp[iz];
            }
        }
        //ModuleBase::timer::end(this->classname, "gatherp_scatters");
        return;
    }
#ifdef __MPI
    //change (nplane fftnxy) to (nplane,nstot)
    // Hence, we can send them at one time.
    const int nstot_gps = this->nstot;
    const int nplane_gps = this->nplane;
    const int* istot2ixy_gps = this->istot2ixy;
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp for
#endif
        for (int istot = 0; istot < nstot_gps; ++istot)
        {
            int ixy = istot2ixy_gps[istot];
            std::complex<T> *outp = &out[istot * nplane_gps];
            std::complex<T> *inp = &in[ixy * nplane_gps];
            for (int iz = 0; iz < nplane_gps; ++iz)
            {
                outp[iz] = inp[iz];
            }
        }
#ifdef _OPENMP
        #pragma omp barrier
    }
#endif

    //exchange data
    //(nplane,nstot) to (numz[ip],ns, poolnproc)
    // OMP barrier: Ensure all OMP threads have finished writing to buffers
    // before MPI communication reads from them. This prevents data race
    // between OMP parallel write (e.g., FFT calculation) and MPI read.
    if(typeid(T) == typeid(double))
    {
        MPI_Alltoallv(out, numr, startr, MPI_DOUBLE_COMPLEX, in, numg, startg, MPI_DOUBLE_COMPLEX, this->pool_world);
    }
    else if(typeid(T) == typeid(float))
    {
        MPI_Alltoallv(out, numr, startr, MPI_COMPLEX, in, numg, startg, MPI_COMPLEX, this->pool_world);
    }

    // change (nz,ns) to (numz[ip],ns, poolnproc)
    const int poolnproc_gps = this->poolnproc;
    const int nst_gps = this->nst;
    const int nz_gps = this->nz;
    const int* numz_gps = this->numz;
    const int* startg_gps = this->startg;
    const int* startz_gps = this->startz;
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
            std::complex<T> *outp = &outp0[is * nz_gps];
            std::complex<T> *inp = &inp0[is * nzip ];
            for (int izip = 0; izip < nzip; ++izip)
            {
                outp[izip] = inp[izip];
            }
        }
    }
#endif
    //ModuleBase::timer::start(this->classname, "gatherp_scatters");
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
    // ModuleBase::timer::start(this->classname, "gathers_scatterp");
    if(this->poolnproc == 1) //In this case nrxx=fftnx*fftny*nz, nst = nstot, 
    {
        const int nrxx = this->nrxx;
        const int nst = this->nst;
        const int nz = this->nz;
        const int* istot2ixy = this->istot2ixy;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(int i = 0; i < nrxx; ++i)
        {
            out[i] = std::complex<T>(0, 0);
        }

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int is = 0 ; is < nst ; ++is)
        {
            int ixy = istot2ixy[is];
            std::complex<T> *outp = &out[ixy*nz];
            std::complex<T> *inp = &in[is*nz];
            for(int iz = 0 ; iz < nz ; ++iz)
            {
                outp[iz] = inp[iz];
            }
        }
        // ModuleBase::timer::end(this->classname, "gathers_scatterp");
        return;
    }
#ifdef __MPI
    // change (nz,ns) to (numz[ip],ns, poolnproc)
    // Hence, we can send them at one time. 
    const int poolnproc = this->poolnproc;
    const int nst = this->nst;
    const int nz = this->nz;
    const int* numz = this->numz;
    const int* startg = this->startg;
    const int* startz = this->startz;
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp for collapse(2)
#endif
        for (int ip = 0; ip < poolnproc ;++ip)
        {
            for (int is = 0; is < nst; ++is)
            {
                int nzip = numz[ip];
                std::complex<T> *outp0 = &out[startg[ip]];
                std::complex<T> *inp0 = &in[startz[ip]];
                std::complex<T> *outp = &outp0[is * nzip];
                std::complex<T> *inp = &inp0[is * nz ];
                for (int izip = 0; izip < nzip; ++izip)
                {
                    outp[izip] = inp[izip];
                }
            }
        }
#ifdef _OPENMP
        #pragma omp barrier
    }
#endif

    //exchange data
    //(numz[ip],ns, poolnproc) to (nplane,nstot)
    // OMP barrier: Ensure all OMP threads have finished writing to buffers
    // before MPI communication reads from them. This prevents data race
    // between OMP parallel write (e.g., FFT calculation) and MPI read.
    if(typeid(T) == typeid(double))
    {
        MPI_Alltoallv(out, numg, startg, MPI_DOUBLE_COMPLEX, in, numr, startr, MPI_DOUBLE_COMPLEX, this->pool_world);
    }
    else if(typeid(T) == typeid(float))
    {
        MPI_Alltoallv(out, numg, startg, MPI_COMPLEX, in, numr, startr, MPI_COMPLEX, this->pool_world);
    }

    const int nrxx_gsp = this->nrxx;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for(int i = 0; i < nrxx_gsp; ++i)
    {
        out[i] = std::complex<T>(0, 0);
    }
    //change (nplane,nstot) to (nplane fftnxy)
    const int nstot = this->nstot;
    const int nplane = this->nplane;
    const int* istot2ixy = this->istot2ixy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int istot = 0;istot < nstot; ++istot)
    {
        int ixy = istot2ixy[istot];
        //int ixy = (ixy / fftny)*ny + ixy % fftny;
        std::complex<T> *outp = &out[ixy * nplane];
        std::complex<T> *inp = &in[istot * nplane];
        for (int iz = 0; iz < nplane; ++iz)
        {
            outp[iz] = inp[iz];
        }
    }
#endif
    // ModuleBase::timer::start(this->classname, "gathers_scatterp");
    return;
}



}
