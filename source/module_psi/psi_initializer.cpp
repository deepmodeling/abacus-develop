#include "psi_initializer.h"
#include "module_base/memory.h"
// basic functions support
#include "module_base/tool_quit.h"
#include "module_base/timer.h"
// three global variables definition
#include "module_base/global_variable.h"

template<typename T, typename Device>
psi::Psi<std::complex<double>>* psi_initializer<T, Device>::allocate(bool only_psig)
{
    ModuleBase::timer::tick("psi_initializer", "allocate");
    /*
        WARNING: when basis_type = "pw", the variable GlobalV::NLOCAL will also be set, in this case, it is set to
        9 = 1 + 3 + 5, which is the maximal number of orbitals spd, I don't think it is reasonable
        The way of calculating this->p_ucell_->natomwfc is, for each atom, read pswfc and for s, it is 1, for p, it is 3
        , then multiplied by the number of atoms, and then add them together.
    */
	int prefactor = 1;
    int nbands_actual = 0;
    if(this->method_ == "random") 
    {
        nbands_actual = GlobalV::NBANDS;
        this->nbands_complem_ = 0;
    }
    else
    {
        if(this->method_.substr(0, 6) == "atomic")
        {
            if(this->p_ucell_->natomwfc >= GlobalV::NBANDS)
            {
                nbands_actual = this->p_ucell_->natomwfc;
                this->nbands_complem_ = 0;
            }
            else
            {
                nbands_actual = GlobalV::NBANDS;
                this->nbands_complem_ = GlobalV::NBANDS - this->p_ucell_->natomwfc;
            }
        }
        else if(this->method_.substr(0, 3) == "nao")
        {
            /*
                previously GlobalV::NLOCAL is used here, however it is wrong. GlobalV::NLOCAL is fixed to 9*nat.
            */
            int nbands_local = 0;
            for(int it = 0; it < this->p_ucell_->ntype; it++)
            {
                for(int ia = 0; ia < this->p_ucell_->atoms[it].na; ia++)
                {
            /* FOR EVERY ATOM */
                    for(int l = 0; l < this->p_ucell_->atoms[it].nwl + 1; l++)
                    {
            /* EVERY ZETA FOR (2l+1) ORBS */
                        /*
                            non-rotate basis, nbands_local*=2 for GlobalV::NPOL = 2 is enough
                        */
                        //nbands_local += this->p_ucell_->atoms[it].l_nchi[l]*(2*l+1) * GlobalV::NPOL;
                        /*
                            rotate basis, nbands_local*=4 for p, d, f,... orbitals, and nbands_local*=2 for s orbitals
                            risky when NSPIN = 4, problematic psi value, needed to be checked
                        */
                        if(l == 0) nbands_local += this->p_ucell_->atoms[it].l_nchi[l] * GlobalV::NPOL;
                        else nbands_local += this->p_ucell_->atoms[it].l_nchi[l]*(2*l+1) * GlobalV::NPOL;
                    }
                }
            }
            if(nbands_local >= GlobalV::NBANDS)
            {
                nbands_actual = nbands_local;
                this->nbands_complem_ = 0;
            }
            else
            {
                nbands_actual = GlobalV::NBANDS;
                this->nbands_complem_ = GlobalV::NBANDS - nbands_local;
            }
        }
    }
	int nkpts_actual = (GlobalV::CALCULATION == "nscf" && this->mem_saver_ == 1)? 1 : this->pw_wfc_->nks;
    int nbasis_actual = this->pw_wfc_->npwk_max * GlobalV::NPOL;
    psi::Psi<std::complex<double>>* psi_out = nullptr;
    if(!only_psig)
    {
        psi_out = new psi::Psi<std::complex<double>>(nkpts_actual, 
                                                     GlobalV::NBANDS, // because no matter what, the wavefunction finally needed has GlobalV::NBANDS bands
                                                     nbasis_actual, 
                                                     this->pw_wfc_->npwk);
        /*
            WARNING: this will cause DIRECT MEMORY LEAK, psi is not properly deallocated
        */
        const size_t memory_cost_psi = 
                nkpts_actual*
                    GlobalV::NBANDS * this->pw_wfc_->npwk_max * GlobalV::NPOL*
                        sizeof(std::complex<double>);
        std::cout << " MEMORY FOR PSI PER PROCESSOR (MB)  : " << double(memory_cost_psi)/1024.0/1024.0 << std::endl;
        ModuleBase::Memory::record("Psi_PW", memory_cost_psi);
    }
    this->psig_ = std::make_shared<psi::Psi<T, Device>>(nkpts_actual, 
                                                        nbands_actual, 
                                                        nbasis_actual, 
                                                        this->pw_wfc_->npwk);
    const size_t memory_cost_psig = 
            nkpts_actual*
                nbands_actual * this->pw_wfc_->npwk_max * GlobalV::NPOL*
                    sizeof(T);
    std::cout << " MEMORY FOR AUXILLARY PSI PER PROCESSOR (MB)  : " << double(memory_cost_psig)/1024.0/1024.0 << std::endl;

    GlobalV::ofs_running << "Allocate memory for psi and psig done.\n"
                         << "Print detailed information of dimension of psi and psig:\n"
                         << "psi: (" << nkpts_actual << ", " << GlobalV::NBANDS << ", " << nbasis_actual << ")\n"
                         << "psig: (" << nkpts_actual << ", " << nbands_actual << ", " << nbasis_actual << ")\n"
                         << "nkpts_actual = " << nkpts_actual << "\n"
                         << "GlobalV::NBANDS = " << GlobalV::NBANDS << "\n"
                         << "nbands_actual = " << nbands_actual << "\n"
                         << "nbands_complem = " << this->nbands_complem_ << "\n"
                         << "nbasis_actual = " << nbasis_actual << "\n"
                         << "npwk_max = " << this->pw_wfc_->npwk_max << "\n"
                         << "npol = " << GlobalV::NPOL << "\n";
    ModuleBase::Memory::record("psigPW", memory_cost_psig);
    ModuleBase::timer::tick("psi_initializer", "allocate");
    return psi_out;
}

template<typename T, typename Device>
void psi_initializer<T, Device>::random_t(T* psi, const int iw_start, const int iw_end, const int ik)
{
    ModuleBase::timer::tick("psi_initializer", "random_t");
    assert(iw_start >= 0);
    const int ng = this->pw_wfc_->npwk[ik];
#ifdef __MPI
    if (this->random_seed_ > 0) // qianrui add 2021-8-13
    {
        srand(unsigned(this->random_seed_ + this->p_parakpts_->startk_pool[GlobalV::MY_POOL] + ik));
        const int nxy = this->pw_wfc_->fftnxy;
        const int nz = this->pw_wfc_->nz;
        const int nstnz = this->pw_wfc_->nst*nz;

        std::vector<Real> stickrr(nz);
        std::vector<Real> stickarg(nz);
        std::vector<Real> tmprr(nstnz);
        std::vector<Real> tmparg(nstnz);
        for (int iw = iw_start; iw < iw_end; iw++)
        {   
            // get the starting memory address of iw band
            T* psi_slice = &(psi[iw * this->pw_wfc_->npwk_max * GlobalV::NPOL]);
            int startig = 0;
            for(int ipol = 0 ; ipol < GlobalV::NPOL ; ++ipol)
            {
                    
                for(int ir=0; ir < nxy; ir++)
                {
                    if(this->pw_wfc_->fftixy2ip[ir] < 0) continue;
                    if(GlobalV::RANK_IN_POOL==0)
                    {
                        for(int iz=0; iz<nz; iz++)
                        {
                            stickrr[iz] = std::rand()/Real(RAND_MAX);
                            stickarg[iz] = std::rand()/Real(RAND_MAX);
                        }
                    }
                    stick_to_pool(stickrr.data(), ir, tmprr.data());
                    stick_to_pool(stickarg.data(), ir, tmparg.data());
                }

                for (int ig = 0;ig < ng;ig++)
                {
                    const double rr = tmprr[this->pw_wfc_->getigl2isz(ik,ig)];
                    const double arg= ModuleBase::TWO_PI * tmparg[this->pw_wfc_->getigl2isz(ik,ig)];
                    const double gk2 = this->pw_wfc_->getgk2(ik,ig);
                    psi_slice[ig+startig] = this->template cast_to_T<T>(std::complex<double>(rr*cos(arg)/(gk2 + 1.0), rr*sin(arg)/(gk2 + 1.0)));
                }
                startig += this->pw_wfc_->npwk_max;
            }
        }
    }
    else
    {
#else  // !__MPI
        if (this->random_seed_ > 0) // qianrui add 2021-8-13
        {
            srand(unsigned(this->random_seed_ + ik));
        }
#endif
        for (int iw = iw_start ;iw < iw_end; iw++)
        {
            T* psi_slice = &(psi[iw * this->pw_wfc_->npwk_max * GlobalV::NPOL]);
            for (int ig = 0; ig < ng; ig++)
            {
                const double rr = std::rand()/double(RAND_MAX); //qianrui add RAND_MAX
                const double arg= ModuleBase::TWO_PI * std::rand()/double(RAND_MAX);
                const double gk2 = this->pw_wfc_->getgk2(ik,ig);
                psi_slice[ig] = this->template cast_to_T<T>(std::complex<double>(rr*cos(arg)/(gk2 + 1.0), rr*sin(arg)/(gk2 + 1.0)));
            }
            if(GlobalV::NPOL==2)
            {
                for (int ig = this->pw_wfc_->npwk_max; ig < this->pw_wfc_->npwk_max + ng; ig++)
                {
                    const double rr = std::rand()/double(RAND_MAX);
                    const double arg= ModuleBase::TWO_PI * std::rand()/double(RAND_MAX);
                    const double gk2 = this->pw_wfc_->getgk2(ik,ig-this->pw_wfc_->npwk_max);
                    psi_slice[ig] = this->template cast_to_T<T>(std::complex<double>(rr*cos(arg)/(gk2 + 1.0), rr*sin(arg)/(gk2 + 1.0)));
                }
            }
        }
#ifdef __MPI
    }
#endif
    ModuleBase::timer::tick("psi_initializer_random", "random_t");
}

#ifdef __MPI
template<typename T, typename Device>
void psi_initializer<T, Device>::stick_to_pool(Real* stick, const int& ir, Real* out) const
{	
    ModuleBase::timer::tick("psi_initializer", "stick_to_pool");
	MPI_Status ierror;
    const int is = this->ixy2is_[ir];
	const int ip = this->pw_wfc_->fftixy2ip[ir];
    const int nz = this->pw_wfc_->nz;

	if(ip == 0 && GlobalV::RANK_IN_POOL ==0)
	{
		for(int iz=0; iz<nz; iz++)
		{
			out[is*nz+iz] = stick[iz];
		}
	}
	else if(ip == GlobalV::RANK_IN_POOL )
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
		for(int iz=0; iz<nz; iz++)
		{
			out[is*nz+iz] = stick[iz];
		}
	}
	else if(GlobalV::RANK_IN_POOL==0)
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
template class psi_initializer<std::complex<double>, psi::DEVICE_CPU>;
template class psi_initializer<std::complex<float>, psi::DEVICE_CPU>;
// gamma point calculation
template class psi_initializer<double, psi::DEVICE_CPU>;
template class psi_initializer<float, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class psi_initializer<std::complex<double>, psi::DEVICE_GPU>;
template class psi_initializer<std::complex<float>, psi::DEVICE_GPU>;
// gamma point calculation
template class psi_initializer<double, psi::DEVICE_GPU>;
template class psi_initializer<float, psi::DEVICE_GPU>;
#endif