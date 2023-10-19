#include "module_psi/psi_initializer/psi_initializer.h"
// basic functions support
#include "module_base/tool_quit.h"
#include "module_base/timer.h"
#include "module_base/memory.h"
// basic data support
#include "module_base/global_variable.h"

template<typename T, typename Device>
#ifdef __MPI
psi_initializer<T, Device>::psi_initializer(
    UnitCell* p_ucell_in,
    ModulePW::PW_Basis_K* pw_wfc_in,
    Parallel_Kpoints* p_parakpts_in,
    Representation<T, Device>* p_rep_in
    ): p_ucell(p_ucell_in),
       pw_wfc(pw_wfc_in),
       p_parakpts(p_parakpts_in),
       p_rep(p_rep_in)
#else
psi_initializer<T, Device>::psi_initializer(
    UnitCell* p_ucell_in,
    ModulePW::PW_Basis_K* pw_wfc_in,
    Representation<T, Device>* p_rep_in
    ): p_ucell(p_ucell_in),
       pw_wfc(pw_wfc_in),
       p_rep(p_rep_in)
#endif
{
    this->ixy2is = new int[this->pw_wfc->fftnxy];
    this->pw_wfc->getfftixy2is(this->ixy2is);
    this->method = GlobalV::init_wfc;
}

template<typename T, typename Device>
psi_initializer<T, Device>::~psi_initializer()
{
    if (this->ixy2is != nullptr) delete[] this->ixy2is;
}

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

    this->p_rep->deallocate_psig();

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
                    psi_slice[ig+startig] = RepIn<T, Device>::template cast_to_T<T>(
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
                psi_slice[ig] = RepIn<T, Device>::template cast_to_T<T>(
                    std::complex<double>(rr*cos(arg)/(gk2 + 1.0), rr*sin(arg)/(gk2 + 1.0)));
            }
            if(GlobalV::NPOL==2)
            {
                for (int ig = this->pw_wfc->npwk_max; ig < this->pw_wfc->npwk_max + ng; ig++)
                {
                    const double rr = std::rand()/double(RAND_MAX);
                    const double arg= ModuleBase::TWO_PI * std::rand()/double(RAND_MAX);
                    const double gk2 = this->pw_wfc->getgk2(ik,ig-this->pw_wfc->npwk_max);
                    psi_slice[ig] = RepIn<T, Device>::template cast_to_T<T>(
                        std::complex<double>(rr*cos(arg)/(gk2 + 1.0), rr*sin(arg)/(gk2 + 1.0)));
                }
            }
        }
#ifdef __MPI
    }
#endif
    ModuleBase::timer::tick("psi_initializer", "random_t");
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

template<typename T, typename Device>
void psi_initializer<T, Device>::initialize()
{
    ModuleBase::timer::tick("psi_initializer", "initialize");
    // get dimension of psig for once
    int psig_nbands = this->p_rep->get_psig()->get_nbands();
    int psig_nbasis = this->p_rep->get_psig()->get_nbasis();
    // register transform, for random, initialize directly
    if(GlobalV::init_wfc.substr(0, 3) == "nao")
    {
        this->p_rep->add_transform_pair("nao", "pw");
    }
    else if(GlobalV::init_wfc.substr(0, 6) == "atomic")
    {
        this->p_rep->add_transform_pair("pao", "pw");
    }
    else
    {
        for(int ik = 0; ik < this->pw_wfc->nks; ik++)
        {
            this->p_rep->get_psig()->fix_k(ik);
            this->random_t(this->p_rep->get_psig_pointer(), 
                           0,                      // range, starting band index (included)
                           psig_nbands,                 // range, ending band index (excluded)
                           ik);                    // kpoint index
        }
        return;
    }
    // do transform
    for(int ik = 0; ik < this->pw_wfc->nks; ik++)
    {
        this->p_rep->align_kpoint(ik);
        // it is a special use of transform function, in parameter list the input wavefunction is set to nullptr, so that
        // transform function will recognize it is just a half transform, from basis function to its pw representation
        // this transform result is stored in representation::psig.
        // For other case in which input wavefunction is not nullptr, transform function will recognize it is a full transform
        // then psig will not be the final product but a intermediate one, for example:
        // a transform from lcao to qo should be:
        // ```
        //     add_transform_pair("nao", "qo");
        //     this->p_rep->transform(psi_in, psi_out);
        // ```
        // , where the psi_in is wavefunction in lcao representation, and psi_out is wavefunction in qo representation, can be
        // created temporarily: 
        // ```
        //     psi::Psi<T, Device> psi_out = psi::Psi<T, Device>(nkpts, nbands, nbasis, npwk);
        // ```
        // and then delete.
        this->p_rep->transform(nullptr, nullptr); // HERE psig is calculated
        // if not enough number of bands initialized, fill with random number
        int nbands_complem = this->get_nbands_complem();
        // random number mixing
        if(GlobalV::init_wfc.find_first_of('+') != std::string::npos)
        {
            if(
                (GlobalV::init_wfc == "atomic+random")
              ||(GlobalV::init_wfc == "nao+random")
                )
            {
                int iband_start = 0;
                int iband_end = GlobalV::NBANDS - nbands_complem;
                psi::Psi<T, Device> psi_random = psi::Psi<T, Device>(1, 
                                                                     iband_end - iband_start, 
                                                                     psig_nbasis, 
                                                                     this->pw_wfc->npwk);
                psi_random.fix_k(0);
                this->random_t(psi_random.get_pointer(), 
                               iband_start, // range, starting band index (included)
                               iband_end,   // range, ending band index (excluded)
                               ik);         // kpoint index
                for(int iband = iband_start; iband < iband_end; iband++)
                {
                    for(int ibasis = 0; ibasis < psig_nbasis; ibasis++)
                    {
                        this->p_rep->get_psig_pointer()[iband * psig_nbasis + ibasis] 
                        = 
                        this->random_mix*
                            psi_random.get_pointer()[iband * psig_nbasis + ibasis]
                        + Real(1.0 - this->random_mix)*
                            this->p_rep->get_psig_pointer()[iband * psig_nbasis + ibasis];
                    }
                }
            }
            else
            {
                ModuleBase::WARNING_QUIT("psi_initializer", "initialize: unknown mixing method");
            }
        }
        // bands complement with random number
        if(nbands_complem != 0)
        {
            int iband_start = GlobalV::NBANDS - nbands_complem;
            int iband_end = GlobalV::NBANDS;
            this->random_t(this->p_rep->get_psig_pointer(), 
                           iband_start, // range, starting band index (included)
                           iband_end,   // range, ending band index (excluded)
                           ik);         // kpoint index
        }
    }
    ModuleBase::timer::tick("psi_initializer", "initialize");
}

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