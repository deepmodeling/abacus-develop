#include "psi_initializer.h"
#include "module_base/memory.h"
// basic functions support
#include "module_base/tool_quit.h"
#include "module_base/timer.h"
// three global variables definition
#include "module_base/global_variable.h"
#include "module_parameter/parameter.h"

#ifdef __MPI
#include "module_base/parallel_reduce.h"
#endif

template<typename T, typename Device>
psi::Psi<std::complex<double>>* psi_initializer<T, Device>::allocate(const bool only_psig)
{
    ModuleBase::timer::tick("psi_initializer", "allocate");
    /*
        WARNING: when basis_type = "pw", the variable PARAM.globalv.nlocal will also be set, in this case, it is set to
        9 = 1 + 3 + 5, which is the maximal number of orbitals spd, I don't think it is reasonable
        The way of calculating this->p_ucell_->natomwfc is, for each atom, read pswfc and for s, it is 1, for p, it is 3
        , then multiplied by the number of atoms, and then add them together.
    */
    int nbands_actual = 0;
    if(this->method_ == "random") 
    {
        nbands_actual = PARAM.inp.nbands;
        this->nbands_complem_ = 0;
    }
    else
    {
        if(this->method_.substr(0, 6) == "atomic")
        {
            nbands_actual = std::max(this->p_ucell_->natomwfc, PARAM.inp.nbands);
            this->nbands_complem_ = nbands_actual - this->p_ucell_->natomwfc;
        }
        else if(this->method_.substr(0, 3) == "nao")
        {
            /*
                previously PARAM.globalv.nlocal is used here, however it is wrong. PARAM.globalv.nlocal is fixed to 9*nat.
            */
            int nbands_local = 0;
            for(int it = 0; it < this->p_ucell_->ntype; it++)
            {
                for(int l = 0; l < this->p_ucell_->atoms[it].nwl + 1; l++)
                {
                    /* EVERY ZETA FOR (2l+1) ORBS */
                    const int nchi = this->p_ucell_->atoms[it].l_nchi[l];
                    const int degen_l = (l == 0)? 1 : 2*l+1;
                    nbands_local += nchi * degen_l * PARAM.globalv.npol * this->p_ucell_->atoms[it].na;
                    /*
                        non-rotate basis, nbands_local*=2 for PARAM.globalv.npol = 2 is enough
                    */
                    //nbands_local += this->p_ucell_->atoms[it].l_nchi[l]*(2*l+1) * PARAM.globalv.npol;
                    /*
                        rotate basis, nbands_local*=4 for p, d, f,... orbitals, and nbands_local*=2 for s orbitals
                        risky when NSPIN = 4, problematic psi value, needed to be checked
                    */
                }
            }
            nbands_actual = std::max(nbands_local, PARAM.inp.nbands);
            this->nbands_complem_ = nbands_actual - nbands_local;
        }
    }
    assert(this->nbands_complem_ >= 0);

	const int nks_psi = (PARAM.inp.calculation == "nscf" && this->mem_saver_ == 1)? 1 : this->pw_wfc_->nks;
    const int nbasis_actual = this->pw_wfc_->npwk_max * PARAM.globalv.npol;
    psi::Psi<std::complex<double>>* psi_out = nullptr;
    if(!only_psig)
    {
        psi_out = new psi::Psi<std::complex<double>>(nks_psi, 
                                                     PARAM.inp.nbands, // because no matter what, the wavefunction finally needed has PARAM.inp.nbands bands
                                                     nbasis_actual, 
                                                     this->pw_wfc_->npwk);
        double memory_cost_psi = sizeof(std::complex<double>) * nks_psi * PARAM.inp.nbands 
                                                              * this->pw_wfc_->npwk_max * PARAM.globalv.npol;
#ifdef __MPI
        // get the correct memory cost for psi by all-reduce sum
        Parallel_Reduce::reduce_all(memory_cost_psi);
#endif
        // std::cout << " MEMORY FOR PSI PER PROCESSOR (MB)  : " << double(memory_cost_psi)/1024.0/1024.0 << std::endl;
        ModuleBase::Memory::record("Psi_PW", memory_cost_psi);
    }
    // psi_initializer also works for basis transformation tasks. In this case, psig needs to allocate memory for 
    // each kpoint, otherwise, for initializing pw wavefunction, only one kpoint's space is enough.
    const int nks_psig = (PARAM.inp.basis_type == "pw")? 1 : nks_psi;
    this->psig_ = std::make_shared<psi::Psi<T, Device>>(nks_psig, 
                                                        nbands_actual, 
                                                        nbasis_actual, 
                                                        this->pw_wfc_->npwk);

    double memory_cost_psig = 
            nks_psig * nbands_actual * this->pw_wfc_->npwk_max * PARAM.globalv.npol * sizeof(T);
#ifdef __MPI
    // get the correct memory cost for psig by all-reduce sum
    Parallel_Reduce::reduce_all(memory_cost_psig);
#endif
    // std::cout << " MEMORY FOR AUXILLARY PSI PER PROCESSOR (MB)  : " << double(memory_cost_psig)/1024.0/1024.0 << std::endl;

    GlobalV::ofs_running << "Allocate memory for psi and psig done.\n"
                         << "Print detailed information of dimension of psi and psig:\n"
                         << "psi: (" << nks_psi << ", " << PARAM.inp.nbands << ", " << nbasis_actual << ")\n"
                         << "psig: (" << nks_psig << ", " << nbands_actual << ", " << nbasis_actual << ")\n"
                         << "nks (psi) = " << nks_psi << "\n"
                         << "nks (psig) = " << nks_psig << "\n"
                         << "PARAM.inp.nbands = " << PARAM.inp.nbands << "\n"
                         << "nbands_actual = " << nbands_actual << "\n"
                         << "nbands_complem = " << this->nbands_complem_ << "\n"
                         << "nbasis_actual = " << nbasis_actual << "\n"
                         << "npwk_max = " << this->pw_wfc_->npwk_max << "\n"
                         << "npol = " << PARAM.globalv.npol << "\n";
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

    // if random seed is specified, then based on this seed to generate random wavefunction
    if (this->random_seed_ > 0) // qianrui add 2021-8-13
    {
#ifdef __MPI // if compile with MPI, then let the seed include the kpoint information
        srand(unsigned(this->random_seed_ + this->p_parakpts_->startk_pool[GlobalV::MY_POOL] + ik));
#else // otherwise, it is the run in serial, without the Parallel_Kpoints class
        srand(unsigned(this->random_seed_ + ik));
#endif
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
            T* psi_slice = &(psi[iw * this->pw_wfc_->npwk_max * PARAM.globalv.npol]);
            int startig = 0;
            for(int ipol = 0 ; ipol < PARAM.globalv.npol ; ++ipol)
            {
                // loop over all fft (x,y), but actually loop over all sticks
                for(int ir = 0; ir < nxy; ir++)
                {
                    // if the stick is not on present processor, then skip
                    if(this->pw_wfc_->fftixy2ip[ir] < 0) { continue; }
                    // otherwise
                    // the following code is very time-consuming, but it can be skipped with pw_seed = 0
                    if(GlobalV::RANK_IN_POOL == 0)
                    {
                        // generate random number for (x,y) and all z, the stick will must
                        // be filled, because length of stick can be no longer than nz
                        // with: rr*exp(i*arg) = rr*cos(arg) + i*rr*sin(arg)
                        for(int iz = 0; iz < nz; iz++)
                        {
                            stickrr[iz] = std::rand()/Real(RAND_MAX);  // amplitude
                            stickarg[iz] = std::rand()/Real(RAND_MAX); // phase
                        }
                    }
#ifdef __MPI // the stick-distribution is not used for non-MPI version
                    // then distribute the data to all processors in the pool 
                    stick_to_pool(stickrr.data(), ir, tmprr.data());
                    stick_to_pool(stickarg.data(), ir, tmparg.data());
#endif
                }
                // then for each g-component, initialize the wavefunction value
                #pragma omp parallel for schedule(static, 4096/sizeof(T))
                for (int ig = 0; ig < ng; ig++)
                {
                    // get the correct value of "rr" and "arg" by indexing map "getigl2isz"
                    const double rr = tmprr[this->pw_wfc_->getigl2isz(ik, ig)];
                    const double arg = ModuleBase::TWO_PI * tmparg[this->pw_wfc_->getigl2isz(ik,ig)];
                    // get the |G+k|^2
                    const double gk2 = this->pw_wfc_->getgk2(ik, ig);
                    // initialize the wavefunction value with rr/(gk2 + 1.0) * exp(i*arg)
                    psi_slice[ig+startig] = this->template cast_to_T<T>(std::complex<double>(rr*cos(arg)/(gk2 + 1.0), rr*sin(arg)/(gk2 + 1.0)));
                }
                startig += this->pw_wfc_->npwk_max; // move to the next polarization
            }
        }
    }
    else // if random seed is not specified, then generate random wavefunction directly
    {
        for (int iw = iw_start ;iw < iw_end; iw++)
        {
            T* psi_slice = &(psi[iw * this->pw_wfc_->npwk_max * PARAM.globalv.npol]); // get the memory to write directly. For nspin 4, nbasis*2

            #pragma omp parallel for schedule(static, 4096/sizeof(T))
            for (int ig = 0; ig < ng; ig++)
            {
                const double rr = std::rand()/double(RAND_MAX); //qianrui add RAND_MAX
                const double arg = ModuleBase::TWO_PI * std::rand()/double(RAND_MAX);
                const double gk2 = this->pw_wfc_->getgk2(ik, ig);
                psi_slice[ig] = this->template cast_to_T<T>(std::complex<double>(rr*cos(arg)/(gk2 + 1.0), rr*sin(arg)/(gk2 + 1.0)));
            }
            if(PARAM.globalv.npol == 2) // additionally for nspin 4...
            {
                #pragma omp parallel for schedule(static, 4096/sizeof(T))
                for (int ig = this->pw_wfc_->npwk_max; ig < this->pw_wfc_->npwk_max + ng; ig++)
                {
                    const double rr = std::rand()/double(RAND_MAX);
                    const double arg = ModuleBase::TWO_PI * std::rand()/double(RAND_MAX);
                    const double gk2 = this->pw_wfc_->getgk2(ik, ig - this->pw_wfc_->npwk_max);
                    psi_slice[ig] = this->template cast_to_T<T>(std::complex<double>(rr*cos(arg)/(gk2 + 1.0), rr*sin(arg)/(gk2 + 1.0)));
                }
            }
        }
    }
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
template class psi_initializer<std::complex<double>, base_device::DEVICE_CPU>;
template class psi_initializer<std::complex<float>, base_device::DEVICE_CPU>;
// gamma point calculation
template class psi_initializer<double, base_device::DEVICE_CPU>;
template class psi_initializer<float, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class psi_initializer<std::complex<double>, base_device::DEVICE_GPU>;
template class psi_initializer<std::complex<float>, base_device::DEVICE_GPU>;
// gamma point calculation
template class psi_initializer<double, base_device::DEVICE_GPU>;
template class psi_initializer<float, base_device::DEVICE_GPU>;
#endif
