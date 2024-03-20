#include "pw_basis.h"
#include "module_base/global_function.h"
#include "module_base/timer.h"
#include "typeinfo"
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
    ModuleBase::timer::tick(this->classname, "gatherp_scatters");
    
    if(this->poolnproc == 1) //In this case nst=nstot, nz = nplane, 
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int is = 0 ; is < this->nst ; ++is)
        {
            int ixy = this->istot2ixy[is];
            //int ixy = (ixy / fftny)*ny + ixy % fftny;
            std::complex<T> *outp = &out[is*nz];
            std::complex<T> *inp = &in[ixy*nz];
            for(int iz = 0 ; iz < this->nz ; ++iz)
            {
                outp[iz] = inp[iz];
            }
        }
        ModuleBase::timer::tick(this->classname, "gatherp_scatters");
        return;
    }
#ifdef __MPI
    //change (nplane fftnxy) to (nplane,nstot)
    // Hence, we can send them at one time.
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int istot = 0;istot < nstot; ++istot)
	{
		int ixy = this->istot2ixy[istot];
        //int ixy = (ixy / fftny)*ny + ixy % fftny;
        std::complex<T> *outp = &out[istot*nplane];
        std::complex<T> *inp = &in[ixy*nplane];
		for (int iz = 0; iz < nplane; ++iz)
		{
			outp[iz] = inp[iz];
		}
	}

    //exchange data
    //(nplane,nstot) to (numz[ip],ns, poolnproc)
    if(typeid(T) == typeid(double))
	    MPI_Alltoallv(out, numr, startr, MPI_DOUBLE_COMPLEX, in, numg, startg, MPI_DOUBLE_COMPLEX, this->pool_world);
    else if(typeid(T) == typeid(float))
        MPI_Alltoallv(out, numr, startr, MPI_COMPLEX, in, numg, startg, MPI_COMPLEX, this->pool_world);
    // change (nz,ns) to (numz[ip],ns, poolnproc)
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int ip = 0; ip < this->poolnproc ;++ip)
	{
		for (int is = 0; is < this->nst; ++is)
		{
            int nzip = this->numz[ip];
            std::complex<T> *outp0 = &out[startz[ip]];
            std::complex<T> *inp0 = &in[startg[ip]];
            std::complex<T> *outp = &outp0[is * nz];
            std::complex<T> *inp = &inp0[is * nzip ];
			for (int izip = 0; izip < nzip; ++izip)
			{
				outp[izip] = inp[izip];
			}
		}
	}
   
#endif
    ModuleBase::timer::tick(this->classname, "gatherp_scatters");
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
    ModuleBase::timer::tick(this->classname, "gathers_scatterp");
    
    if(this->poolnproc == 1) //In this case nrxx=fftnx*fftny*nz, nst = nstot, 
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(T))
#endif
        for(int i = 0; i < this->nrxx; ++i)
        {
            out[i] = std::complex<T>(0, 0);
        }

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int is = 0 ; is < this->nst ; ++is)
        {
            int ixy = istot2ixy[is];
            //int ixy = (ixy / fftny)*ny + ixy % fftny;
            std::complex<T> *outp = &out[ixy*nz];
            std::complex<T> *inp = &in[is*nz];
            for(int iz = 0 ; iz < this->nz ; ++iz)
            {
                outp[iz] = inp[iz];
            }
        }
        ModuleBase::timer::tick(this->classname, "gathers_scatterp");
        return;
    }
#ifdef __MPI
    // change (nz,ns) to (numz[ip],ns, poolnproc)
    // Hence, we can send them at one time. 
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int ip = 0; ip < this->poolnproc ;++ip)
	{
		for (int is = 0; is < this->nst; ++is)
		{
            int nzip = this->numz[ip];
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

	//exchange data
    //(numz[ip],ns, poolnproc) to (nplane,nstot)
    if(typeid(T) == typeid(double))
	    MPI_Alltoallv(out, numg, startg, MPI_DOUBLE_COMPLEX, in, numr, startr, MPI_DOUBLE_COMPLEX, this->pool_world);
    else if(typeid(T) == typeid(float))
        MPI_Alltoallv(out, numg, startg, MPI_COMPLEX, in, numr, startr, MPI_COMPLEX, this->pool_world);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(T))
#endif
    for(int i = 0; i < this->nrxx; ++i)
    {
        out[i] = std::complex<T>(0, 0);
    }
    //change (nplane,nstot) to (nplane fftnxy)
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int istot = 0;istot < nstot; ++istot)
	{
		int ixy = this->istot2ixy[istot];
        //int ixy = (ixy / fftny)*ny + ixy % fftny;
        std::complex<T> *outp = &out[ixy * nplane];
        std::complex<T> *inp = &in[istot * nplane];
		for (int iz = 0; iz < nplane; ++iz)
		{
			outp[iz] = inp[iz];
		}
	}

#endif
    ModuleBase::timer::tick(this->classname, "gathers_scatterp");
    return;
}



}