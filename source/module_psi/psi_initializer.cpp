#include "psi_initializer.h"

#include "module_base/memory.h"
// basic functions support
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
// three global variables definition
#include "module_base/global_variable.h"
#include "module_parameter/parameter.h"

#ifdef __MPI
#include "module_base/parallel_reduce.h"
#endif

template <typename T>
void psi_initializer<T>::initialize(const Structure_Factor* sf,
                                    const ModulePW::PW_Basis_K* pw_wfc,
                                    const UnitCell* p_ucell,
                                    const K_Vectors* p_kv_in,
                                    const int& random_seed,
                                    const pseudopot_cell_vnl* p_pspot_nl,
                                    const int& rank)
{
    this->sf_ = sf;
    this->pw_wfc_ = pw_wfc;
    this->p_ucell_ = p_ucell;
    this->p_kv = p_kv_in;
    this->random_seed_ = random_seed;
    this->p_pspot_nl_ = p_pspot_nl;
}

template <typename T>
void psi_initializer<T>::random_t(T* psi, const int iw_start, const int iw_end, const int ik, const int mode)
{
    ModuleBase::timer::tick("psi_initializer", "random_t");
    assert(mode <= 1);
    assert(iw_start >= 0);
    const int ng = this->pw_wfc_->npwk[ik];
    const int npwk_max = this->pw_wfc_->npwk_max;
    const int npol = PARAM.globalv.npol;

    // If random seed is specified, then generate random wavefunction satisfying that
    // it can generate the same results using different number of processors.
    if (this->random_seed_ > 0) // qianrui add 2021-8-13
    {
#ifdef __MPI
        srand(unsigned(this->random_seed_ + this->p_kv->ik2iktot[ik]));
#else
        srand(unsigned(this->random_seed_ + ik));
#endif
        const int nxy = this->pw_wfc_->fftnxy;
        const int nz = this->pw_wfc_->nz;
        const int nstnz = this->pw_wfc_->nst * nz;

        std::vector<Real> stickrr(nz);
        std::vector<Real> stickarg(nz);
        std::vector<Real> tmprr(nstnz);
        std::vector<Real> tmparg(nstnz);

        for (int iw = iw_start; iw < iw_end; iw++)
        {
            // get the starting memory address of iw band
            T* psi_slice = &(psi[iw * npwk_max * npol]);
            for (int ipol = 0; ipol < npol; ++ipol)
            {
                // loop over all fft (x,y), but actually loop over all sticks
                for (int ir = 0; ir < nxy; ir++)
                {
                    // if the stick is not on present processor, then skip
                    if (this->pw_wfc_->fftixy2ip[ir] < 0)
                    {
                        continue;
                    }
                    // otherwise
                    // the following code is very time-consuming, but it can be skipped with pw_seed = 0
                    if (GlobalV::RANK_IN_POOL == 0)
                    {
                        // generate random number for (x,y) and all z, the stick will must
                        // be filled, because length of stick can be no longer than nz
                        // with: rr*exp(i*arg) = rr*cos(arg) + i*rr*sin(arg)
                        for (int iz = 0; iz < nz; iz++)
                        {
                            stickrr[iz] = std::rand() / Real(RAND_MAX);  // amplitude
                            stickarg[iz] = std::rand() / Real(RAND_MAX); // phase
                        }
                    }
#ifdef __MPI // the stick-distribution is not used for non-MPI version
             // then distribute the data to all processors in the pool
                    stick_to_pool(stickrr.data(), ir, tmprr.data());
                    stick_to_pool(stickarg.data(), ir, tmparg.data());
#endif
                }
                // then for each g-component, initialize the wavefunction value
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(T))
#endif
                for (int ig = 0; ig < ng; ig++)
                {
                    // get the correct value of "rr" and "arg" by indexing map "getigl2isz"
                    const int isz = this->pw_wfc_->getigl2isz(ik, ig);
                    const double rr = tmprr[isz];
                    const double arg = ModuleBase::TWO_PI * tmparg[isz];
                    // initialize the wavefunction value with rr * exp(i*arg)
                    psi_slice[ig] = this->template cast_to_T<T>(std::complex<double>(rr * cos(arg), rr * sin(arg)));
                }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(T))
#endif
                for (int ig = ng; ig < npwk_max; ++ig)
                {
                    psi_slice[ig] = static_cast<T>(0.0);
                }
                psi_slice += npwk_max; // move to the next polarization
            }
        }
    }
    // If random seed is not specified, then generate random wavefunction directly
    // It does not guarantee the same results using different number of processors.
    else
    {
        for (int iw = iw_start; iw < iw_end; iw++)
        {
            T* psi_slice = &(psi[iw * npwk_max * npol]); // get the memory to write directly. For nspin 4, nbasis*2
            // donot use openmp here, because the random number generator is not thread-safe
            for (int ig = 0; ig < ng; ig++)
            {
                const double rr = std::rand() / double(RAND_MAX);
                const double arg = ModuleBase::TWO_PI * std::rand() / double(RAND_MAX);
                const double gk2 = this->pw_wfc_->getgk2(ik, ig);
                psi_slice[ig] = this->template cast_to_T<T>(
                    std::complex<double>(rr * cos(arg) / (gk2 + 1.0), rr * sin(arg) / (gk2 + 1.0)));
            }
            if (npol == 2)
            {
                for (int ig = npwk_max; ig < npwk_max + ng; ig++)
                {
                    const double rr = std::rand() / double(RAND_MAX);
                    const double arg = ModuleBase::TWO_PI * std::rand() / double(RAND_MAX);
                    const double gk2 = this->pw_wfc_->getgk2(ik, ig - npwk_max);
                    psi_slice[ig] = this->template cast_to_T<T>(
                        std::complex<double>(rr * cos(arg) / (gk2 + 1.0), rr * sin(arg) / (gk2 + 1.0)));
                }
            }
        }
    }
    if (mode == 1)
    {
        for (int iw = iw_start; iw < iw_end; iw++)
        {
            T* psi_slice = &(psi[iw * npwk_max * npol]);
            for (int ipol = 0; ipol < npol; ipol++)
            {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(T))
#endif
                for (int ig = 0; ig < ng; ig++)
                {
                    const double gk2 = this->pw_wfc_->getgk2(ik, ig);
                    const Real inv_gk2 = 1.0 / (gk2 + 1.0);
                    psi_slice[ig] *= inv_gk2;
                }
                psi_slice += npwk_max;
            }
        }
    }
    ModuleBase::timer::tick("psi_initializer_random", "random_t");
}

#ifdef __MPI
template <typename T>
void psi_initializer<T>::stick_to_pool(Real* stick, const int& ir, Real* out) const
{
    ModuleBase::timer::tick("psi_initializer", "stick_to_pool");
    MPI_Status ierror;
    const int is = this->ixy2is_[ir];
    const int ip = this->pw_wfc_->fftixy2ip[ir];
    const int nz = this->pw_wfc_->nz;

    if (ip == 0 && GlobalV::RANK_IN_POOL == 0)
    {
        for (int iz = 0; iz < nz; iz++)
        {
            out[is * nz + iz] = stick[iz];
        }
    }
    else if (ip == GlobalV::RANK_IN_POOL)
    {
        if (std::is_same<Real, double>::value)
        {
            MPI_Recv(stick, nz, MPI_DOUBLE, 0, ir, POOL_WORLD, &ierror);
        }
        else if (std::is_same<Real, float>::value)
        {
            MPI_Recv(stick, nz, MPI_FLOAT, 0, ir, POOL_WORLD, &ierror);
        }
        else
        {
            ModuleBase::WARNING_QUIT("psi_initializer", "stick_to_pool: Real type not supported");
        }
        for (int iz = 0; iz < nz; iz++)
        {
            out[is * nz + iz] = stick[iz];
        }
    }
    else if (GlobalV::RANK_IN_POOL == 0)
    {
        if (std::is_same<Real, double>::value)
        {
            MPI_Send(stick, nz, MPI_DOUBLE, ip, ir, POOL_WORLD);
        }
        else if (std::is_same<Real, float>::value)
        {
            MPI_Send(stick, nz, MPI_FLOAT, ip, ir, POOL_WORLD);
        }
        else
        {
            ModuleBase::WARNING_QUIT("psi_initializer", "stick_to_pool: Real type not supported");
        }
    }

    return;
    ModuleBase::timer::tick("psi_initializer", "stick_to_pool");
}
#endif

// explicit instantiation
template class psi_initializer<std::complex<double>>;
template class psi_initializer<std::complex<float>>;
// gamma point calculation
template class psi_initializer<double>;
template class psi_initializer<float>;
