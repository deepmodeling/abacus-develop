#include "module_psi/psi_initializer/psi_initializer.h"

template<typename T, typename Device>
psi::Psi<std::complex<double>>* psi_initializer<T, Device>::allocate()
{
    ModuleBase::timer::tick("psi_initializer", "allocate");
    /*
        WARNING: when basis_type = "pw", the variable GlobalV::NLOCAL will also be set, in this case, it is set to
        9 = 1 + 3 + 5, which is the maximal number of orbitals spd, I don't think it is reasonable
        The way of calculating this->p_ucell->natomwfc is, for each atom, read pswfc and for s, it is 1, for p, it is 3
        , then multiplied by the number of atoms, and then add them together.
    */

    if(this->p_rep->psig != nullptr) delete this->p_rep->psig;
	int prefactor = 1;
    int nbands_actual = 0;
    if(GlobalV::init_wfc == "random") 
    {
        nbands_actual = GlobalV::NBANDS;
        this->nbands_complem = 0;
    }
    else
    {
        if(GlobalV::init_wfc.substr(0, 6) == "atomic")
        {
            if(this->p_ucell->natomwfc >= GlobalV::NBANDS)
            {
                nbands_actual = this->p_ucell->natomwfc;
                this->nbands_complem = 0;
            }
            else
            {
                nbands_actual = GlobalV::NBANDS;
                this->nbands_complem = GlobalV::NBANDS - this->p_ucell->natomwfc;
            }
        }
        else if(GlobalV::init_wfc.substr(0, 3) == "nao")
        {
            /*
                previously GlobalV::NLOCAL is used here, however it is wrong. GlobalV::NLOCAL is fixed to 9*nat.
            */
            int nbands_local = 0;
            for(int it = 0; it < this->p_ucell->ntype; it++)
            {
                for(int ia = 0; ia < this->p_ucell->atoms[it].na; ia++)
                {
            /* FOR EVERY ATOM */
                    for(int l = 0; l < this->p_ucell->atoms[it].nwl + 1; l++)
                    {
            /* EVERY ZETA FOR (2l+1) ORBS, for NSPIN = 4, 4 times (because change to rotate base) */
                        nbands_local += this->p_ucell->atoms[it].l_nchi[l]*(2*l+1) * GlobalV::NPOL;
                        /* the following is for rotate base. However it cannot yield correct result because many zeros yielded and cause diag failure */
                        /*
                        if(l == 0)
                        {
                            nbands_local += this->p_ucell->atoms[it].l_nchi[l]*(2*l+1) * GlobalV::NPOL;
                        }
                        else
                        {
                            nbands_local += this->p_ucell->atoms[it].l_nchi[l]*(2*l+1) * GlobalV::NPOL * GlobalV::NPOL;
                        }
                        */
                    }
                }
            }
            if(nbands_local >= GlobalV::NBANDS)
            {
                nbands_actual = nbands_local;
                this->nbands_complem = 0;
            }
            else
            {
                nbands_actual = GlobalV::NBANDS;
                this->nbands_complem = GlobalV::NBANDS - nbands_local;
            }
        }
    }
	int nkpts_actual = (GlobalV::CALCULATION == "nscf" && this->mem_saver == 1)? 
                            1 : this->pw_wfc->nks;
    int nbasis_actual = this->pw_wfc->npwk_max * GlobalV::NPOL;
    psi::Psi<std::complex<double>>* psi_out = nullptr;
    psi_out = new psi::Psi<std::complex<double>>(
        nkpts_actual, 
            GlobalV::NBANDS, // because no matter what, the wavefunction finally needed has GlobalV::NBANDS bands
                nbasis_actual, 
                    this->pw_wfc->npwk);
    this->p_rep->allocate_psig(nkpts_actual, nbands_actual, nbasis_actual, this->pw_wfc->npwk);
    // memory cost estimation per processor
    const size_t memory_cost_psi = 
            nkpts_actual*
                GlobalV::NBANDS * this->pw_wfc->npwk_max * GlobalV::NPOL*
                    sizeof(std::complex<double>);
	std::cout << " MEMORY FOR PSI PER PROCESSOR (MB)  : " << double(memory_cost_psi)/1024.0/1024.0 << std::endl;
    const size_t memory_cost_psig = 
            nkpts_actual*
                nbands_actual * this->pw_wfc->npwk_max * GlobalV::NPOL*
                    sizeof(T);
    std::cout << " MEMORY FOR AUXILLARY PSI PER PROCESSOR (MB)  : " << double(memory_cost_psig)/1024.0/1024.0 << std::endl;
	ModuleBase::Memory::record("Psi_PW", memory_cost_psi);
    ModuleBase::Memory::record("PsiG_PW", memory_cost_psig);
    ModuleBase::timer::tick("psi_initializer", "allocate");
    return psi_out;
}

template<typename T, typename Device>
void psi_initializer<T, Device>::random_t(T* psi, const int iw_start, const int iw_end, const int ik)
{
    ModuleBase::timer::tick("psi_initializer", "random_t");
    assert(iw_start >= 0);
    const int ng = this->pw_wfc->npwk[ik];
#ifdef __MPI
    if (this->random_seed > 0) // qianrui add 2021-8-13
    {
        srand(unsigned(this->random_seed + this->p_parakpts->startk_pool[GlobalV::MY_POOL] + ik));
        const int nxy = this->pw_wfc->fftnxy;
        const int nz = this->pw_wfc->nz;
        const int nstnz = this->pw_wfc->nst*nz;

        Real *stickrr = new Real[nz];
        Real *stickarg = new Real[nz];
        Real *tmprr = new Real[nstnz];
        Real *tmparg = new Real[nstnz];
        for (int iw = iw_start; iw < iw_end; iw++)
        {   
            // get the starting memory address of iw band
            T* psi_slice = &(psi[iw * this->pw_wfc->npwk_max * GlobalV::NPOL]);
            int startig = 0;
            for(int ipol = 0 ; ipol < GlobalV::NPOL ; ++ipol)
            {
                    
                for(int ir=0; ir < nxy; ir++)
                {
                    if(this->pw_wfc->fftixy2ip[ir] < 0) continue;
                    if(GlobalV::RANK_IN_POOL==0)
                    {
                        for(int iz=0; iz<nz; iz++)
                        {
                            stickrr[ iz ] = std::rand()/Real(RAND_MAX);
                            stickarg[ iz ] = std::rand()/Real(RAND_MAX);
                        }
                    }
                    stick_to_pool(stickrr, ir, tmprr);
                    stick_to_pool(stickarg, ir, tmparg);
                }

                for (int ig = 0;ig < ng;ig++)
                {
                    const double rr = tmprr[this->pw_wfc->getigl2isz(ik,ig)];
                    const double arg= ModuleBase::TWO_PI * tmparg[this->pw_wfc->getigl2isz(ik,ig)];
                    const double gk2 = this->pw_wfc->getgk2(ik,ig);
                    psi_slice[ig+startig] = this->p_rep->repin->template cast_to_T<T>(
                        std::complex<double>(rr*cos(arg)/(gk2 + 1.0), rr*sin(arg)/(gk2 + 1.0)));
                }
                startig += this->pw_wfc->npwk_max;
            }
        }
        delete[] stickrr;
        delete[] stickarg;
        delete[] tmprr;
        delete[] tmparg;
    }
    else
    {
#else  // !__MPI
        if (this->random_seed > 0) // qianrui add 2021-8-13
        {
            srand(unsigned(this->random_seed + ik));
        }
#endif
        for (int iw = iw_start ;iw < iw_end; iw++)
        {
            T* psi_slice = &(psi[iw * this->pw_wfc->npwk_max * GlobalV::NPOL]);
            for (int ig = 0; ig < ng; ig++)
            {
                const double rr = std::rand()/double(RAND_MAX); //qianrui add RAND_MAX
                const double arg= ModuleBase::TWO_PI * std::rand()/double(RAND_MAX);
                const double gk2 = this->pw_wfc->getgk2(ik,ig);
                psi_slice[ig] = this->p_rep->repin->template cast_to_T<T>(
                    std::complex<double>(rr*cos(arg)/(gk2 + 1.0), rr*sin(arg)/(gk2 + 1.0)));
            }
            if(GlobalV::NPOL==2)
            {
                for (int ig = this->pw_wfc->npwk_max; ig < this->pw_wfc->npwk_max + ng; ig++)
                {
                    const double rr = std::rand()/double(RAND_MAX);
                    const double arg= ModuleBase::TWO_PI * std::rand()/double(RAND_MAX);
                    const double gk2 = this->pw_wfc->getgk2(ik,ig-this->pw_wfc->npwk_max);
                    psi_slice[ig] = this->p_rep->repin->template cast_to_T<T>(
                        std::complex<double>(rr*cos(arg)/(gk2 + 1.0), rr*sin(arg)/(gk2 + 1.0)));
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
    const int is = this->ixy2is[ir];
	const int ip = this->pw_wfc->fftixy2ip[ir];
    const int nz = this->pw_wfc->nz;

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
