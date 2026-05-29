#include "pw_basis.h"
#include "source_base/global_function.h"
#include "source_base/timer.h"
#include <typeinfo>
#include <vector>

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
    ModuleBase::timer::start(this->classname, "gatherp_scatters");

    if(this->poolnproc == 1) //In this case nst=nstot, nz = nplane,
    {
        const int nst_ = this->nst;
        const int nz_ = this->nz;
        const int* istot2ixy_ = this->istot2ixy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int is = 0 ; is < nst_ ; ++is)
        {
            int ixy = istot2ixy_[is];
            std::complex<T> *outp = &out[is*nz_];
            std::complex<T> *inp = &in[ixy*nz_];
            for(int iz = 0 ; iz < nz_ ; ++iz)
            {
                outp[iz] = inp[iz];
            }
        }
        ModuleBase::timer::end(this->classname, "gatherp_scatters");
        return;
    }


#ifdef __MPI
    //change (nplane fftnxy) to (nplane,nstot)
    // Hence, we can send them at one time.
    ModuleBase::timer::start(this->classname, "gatherp_pack");
    const int nstot_gps = this->nstot;
    const int nplane_gps = this->nplane;
    const int* istot2ixy_gps = this->istot2ixy;
    const int* numg_gps = this->numg;
    const int* numr_gps = this->numr;
    const int* startg_gps = this->startg;
    const int* startr_gps = this->startr;
    const int poolrank_gps = this->poolrank;
    const int poolnproc_gps = this->poolnproc;
    const int send_count_gps = startr_gps[poolnproc_gps - 1] + numr_gps[poolnproc_gps - 1];
    const int recv_count_gps = startg_gps[poolnproc_gps - 1] + numg_gps[poolnproc_gps - 1];
    std::complex<T>* commbuf = this->acquire_comm_workbuf<T>(send_count_gps + recv_count_gps);
    std::complex<T>* sendbuf = commbuf;
    // Keep a dedicated receive slice so ranks with zero local planes do not
    // need their logical input array to also satisfy the receive-buffer bound.
    std::complex<T>* recvbuf = commbuf + send_count_gps;
    if (nplane_gps > 0)
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int istot = 0; istot < nstot_gps; ++istot)
        {
            int ixy = istot2ixy_gps[istot];
            std::complex<T> *outp = &sendbuf[istot * nplane_gps];
            std::complex<T> *inp = &in[ixy * nplane_gps];
            for (int iz = 0; iz < nplane_gps; ++iz)
            {
                outp[iz] = inp[iz];
            }
        }
    }
    ModuleBase::timer::end(this->classname, "gatherp_pack");

    //exchange data
    //(nplane,nstot) to (numz[ip],ns, poolnproc)
    MPI_Datatype mpi_type = MPI_DATATYPE_NULL;
    if(typeid(T) == typeid(double))
    {
        mpi_type = MPI_DOUBLE_COMPLEX;
    }
    else if(typeid(T) == typeid(float))
    {
        mpi_type = MPI_COMPLEX;
    }
    else
    {
        ModuleBase::WARNING_QUIT("PW_Basis::gatherp_scatters", "Unsupported data type for MPI gather/scatter");
    }
    std::vector<MPI_Request> recv_requests(poolnproc_gps, MPI_REQUEST_NULL);
    std::vector<MPI_Request> send_requests(poolnproc_gps, MPI_REQUEST_NULL);
    std::vector<MPI_Status> recv_status(poolnproc_gps);
    std::vector<int> recv_indices(poolnproc_gps, MPI_UNDEFINED);
    int active_recvs = 0;
    int active_sends = 0;

    ModuleBase::timer::start(this->classname, "gatherp_alltoallv");
    for (int ip = 0; ip < poolnproc_gps; ++ip)
    {
        if (ip == poolrank_gps || numg_gps[ip] == 0)
        {
            continue;
        }
        MPI_Irecv(&recvbuf[startg_gps[ip]], numg_gps[ip], mpi_type, ip, 0, this->pool_world, &recv_requests[ip]);
        ++active_recvs;
    }
    for (int ip = 0; ip < poolnproc_gps; ++ip)
    {
        if (ip == poolrank_gps || numr_gps[ip] == 0)
        {
            continue;
        }
        MPI_Isend(&sendbuf[startr_gps[ip]], numr_gps[ip], mpi_type, ip, 0, this->pool_world, &send_requests[ip]);
        ++active_sends;
    }
    ModuleBase::timer::end(this->classname, "gatherp_alltoallv");

    // change (nz,ns) to (numz[ip],ns, poolnproc)
    const int nst_gps = this->nst;
    const int nz_gps = this->nz;
    const int* numz_gps = this->numz;
    const int* startz_gps = this->startz;
    auto unpack_peer = [&](const int ip)
    {
        const int nzip = numz_gps[ip];
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int is = 0; is < nst_gps; ++is)
        {
            std::complex<T> *outp = &out[is * nz_gps + startz_gps[ip]];
            std::complex<T> *inp = &recvbuf[startg_gps[ip] + is * nzip];
            for (int izip = 0; izip < nzip; ++izip)
            {
                outp[izip] = inp[izip];
            }
        }
    };

    ModuleBase::timer::start(this->classname, "gatherp_unpack");
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < numg_gps[poolrank_gps]; ++i)
    {
        recvbuf[startg_gps[poolrank_gps] + i] = sendbuf[startr_gps[poolrank_gps] + i];
    }
    unpack_peer(poolrank_gps);
    ModuleBase::timer::end(this->classname, "gatherp_unpack");

    while (active_recvs > 0)
    {
        int outcount = 0;
        ModuleBase::timer::start(this->classname, "gatherp_alltoallv");
        MPI_Waitsome(poolnproc_gps,
                     recv_requests.data(),
                     &outcount,
                     recv_indices.data(),
                     recv_status.data());
        ModuleBase::timer::end(this->classname, "gatherp_alltoallv");
        if (outcount == MPI_UNDEFINED)
        {
            break;
        }
        for (int idx = 0; idx < outcount; ++idx)
        {
            ModuleBase::timer::start(this->classname, "gatherp_unpack");
            unpack_peer(recv_indices[idx]);
            ModuleBase::timer::end(this->classname, "gatherp_unpack");
        }
        active_recvs -= outcount;
    }

    if (active_sends > 0)
    {
        ModuleBase::timer::start(this->classname, "gatherp_alltoallv");
        MPI_Waitall(poolnproc_gps, send_requests.data(), MPI_STATUSES_IGNORE);
        ModuleBase::timer::end(this->classname, "gatherp_alltoallv");
    }
#endif
    ModuleBase::timer::end(this->classname, "gatherp_scatters");
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
    ModuleBase::timer::start(this->classname, "gathers_scatterp");
    if(this->poolnproc == 1) //In this case nrxx=fftnx*fftny*nz, nst = nstot,
    {
        const int nrxx_ = this->nrxx;
        const int nst_ = this->nst;
        const int nz_ = this->nz;
        const int* istot2ixy_ = this->istot2ixy;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(int i = 0; i < nrxx_; ++i)
        {
            out[i] = std::complex<T>(0, 0);
        }

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int is = 0 ; is < nst_ ; ++is)
        {
            int ixy = istot2ixy_[is];
            std::complex<T> *outp = &out[ixy*nz_];
            std::complex<T> *inp = &in[is*nz_];
            for(int iz = 0 ; iz < nz_ ; ++iz)
            {
                outp[iz] = inp[iz];
            }
        }
        ModuleBase::timer::end(this->classname, "gathers_scatterp");
        return;
    }


#ifdef __MPI
    // change (nz,ns) to (numz[ip],ns, poolnproc)
    // Hence, we can send them at one time.
    ModuleBase::timer::start(this->classname, "gathers_pack");
    const int poolnproc_ = this->poolnproc;
    const int nst_ = this->nst;
    const int nz_ = this->nz;
    const int* numz_ = this->numz;
    const int* startg_ = this->startg;
    const int* startz_ = this->startz;
    const int* nst_per_ = this->nst_per;
    const int* startr_ = this->startr;
    const int poolrank_ = this->poolrank;
    const int send_count_ = startg_[poolnproc_ - 1] + this->numg[poolnproc_ - 1];
    const int recv_count_ = startr_[poolnproc_ - 1] + this->numr[poolnproc_ - 1];
    std::complex<T>* commbuf = this->acquire_comm_workbuf<T>(send_count_ + recv_count_);
    std::complex<T>* sendbuf = commbuf;
    std::complex<T>* recvbuf = commbuf + send_count_;
#ifdef _OPENMP
    #pragma omp parallel for collapse(2)
#endif
    for (int ip = 0; ip < poolnproc_ ;++ip)
    {
        for (int is = 0; is < nst_; ++is)
        {
            int nzip = numz_[ip];
            std::complex<T> *outp0 = &sendbuf[startg_[ip]];
            std::complex<T> *inp0 = &in[startz_[ip]];
            std::complex<T> *outp = &outp0[is * nzip];
            std::complex<T> *inp = &inp0[is * nz_ ];
            for (int izip = 0; izip < nzip; ++izip)
            {
                outp[izip] = inp[izip];
            }
        }
    }
    ModuleBase::timer::end(this->classname, "gathers_pack");

    //exchange data
    //(numz[ip],ns, poolnproc) to (nplane,nstot)
    MPI_Datatype mpi_type = MPI_DATATYPE_NULL;
    if(typeid(T) == typeid(double))
    {
        mpi_type = MPI_DOUBLE_COMPLEX;
    }
    else if(typeid(T) == typeid(float))
    {
        mpi_type = MPI_COMPLEX;
    }
    else
    {
        ModuleBase::WARNING_QUIT("PW_Basis::gathers_scatterp", "Unsupported data type for MPI gather/scatter");
    }
    std::vector<MPI_Request> recv_requests(poolnproc_, MPI_REQUEST_NULL);
    std::vector<MPI_Request> send_requests(poolnproc_, MPI_REQUEST_NULL);
    std::vector<MPI_Status> recv_status(poolnproc_);
    std::vector<int> recv_indices(poolnproc_, MPI_UNDEFINED);
    int active_recvs = 0;
    int active_sends = 0;

    ModuleBase::timer::start(this->classname, "gathers_alltoallv");
    for (int ip = 0; ip < poolnproc_; ++ip)
    {
        if (ip == poolrank_ || this->numr[ip] == 0)
        {
            continue;
        }
        MPI_Irecv(&recvbuf[startr_[ip]], this->numr[ip], mpi_type, ip, 0, this->pool_world, &recv_requests[ip]);
        ++active_recvs;
    }
    for (int ip = 0; ip < poolnproc_; ++ip)
    {
        if (ip == poolrank_ || this->numg[ip] == 0)
        {
            continue;
        }
        MPI_Isend(&sendbuf[startg_[ip]], this->numg[ip], mpi_type, ip, 0, this->pool_world, &send_requests[ip]);
        ++active_sends;
    }
    ModuleBase::timer::end(this->classname, "gathers_alltoallv");

    ModuleBase::timer::start(this->classname, "gathers_clear");
    const int nrxx_gsp = this->nrxx;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for(int i = 0; i < nrxx_gsp; ++i)
    {
        out[i] = std::complex<T>(0, 0);
    }
    ModuleBase::timer::end(this->classname, "gathers_clear");

    //change (nplane,nstot) to (nplane fftnxy)
    const int nplane = this->nplane;
    const int* istot2ixy = this->istot2ixy;
    std::vector<int> istot_offsets(poolnproc_, 0);
    for (int ip = 1; ip < poolnproc_; ++ip)
    {
        istot_offsets[ip] = istot_offsets[ip - 1] + nst_per_[ip - 1];
    }
    auto unpack_peer = [&](const int ip)
    {
        const int peer_nst = nst_per_[ip];
        if (peer_nst == 0 || nplane == 0)
        {
            return;
        }
        const int istot0 = istot_offsets[ip];
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int is = 0; is < peer_nst; ++is)
        {
            const int istot = istot0 + is;
            const int ixy = istot2ixy[istot];
            std::complex<T> *outp = &out[ixy * nplane];
            std::complex<T> *inp = &recvbuf[startr_[ip] + is * nplane];
            for (int iz = 0; iz < nplane; ++iz)
            {
                outp[iz] = inp[iz];
            }
        }
    };

    ModuleBase::timer::start(this->classname, "gathers_unpack");
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < this->numr[poolrank_]; ++i)
    {
        recvbuf[startr_[poolrank_] + i] = sendbuf[startg_[poolrank_] + i];
    }
    unpack_peer(poolrank_);
    ModuleBase::timer::end(this->classname, "gathers_unpack");

    while (active_recvs > 0)
    {
        int outcount = 0;
        ModuleBase::timer::start(this->classname, "gathers_alltoallv");
        MPI_Waitsome(poolnproc_,
                     recv_requests.data(),
                     &outcount,
                     recv_indices.data(),
                     recv_status.data());
        ModuleBase::timer::end(this->classname, "gathers_alltoallv");
        if (outcount == MPI_UNDEFINED)
        {
            break;
        }
        for (int idx = 0; idx < outcount; ++idx)
        {
            ModuleBase::timer::start(this->classname, "gathers_unpack");
            unpack_peer(recv_indices[idx]);
            ModuleBase::timer::end(this->classname, "gathers_unpack");
        }
        active_recvs -= outcount;
    }

    if (active_sends > 0)
    {
        ModuleBase::timer::start(this->classname, "gathers_alltoallv");
        MPI_Waitall(poolnproc_, send_requests.data(), MPI_STATUSES_IGNORE);
        ModuleBase::timer::end(this->classname, "gathers_alltoallv");
    }
#endif
    ModuleBase::timer::end(this->classname, "gathers_scatterp");
    return;
}



}
