#include "psi_initializer.h"
#include "module_base/memory.h"
// basic functions support
#include "module_base/tool_quit.h"
#include "module_base/timer.h"
// three global variables definition
#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

psi_initializer::psi_initializer(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in): sf(sf_in), pw_wfc(pw_wfc_in)
{
    this->ixy2is = new int[this->pw_wfc->fftnxy];
    this->pw_wfc->getfftixy2is(this->ixy2is);
}


psi_initializer::~psi_initializer()
{
    delete[] this->ixy2is;
    if (this->psig != nullptr) delete this->psig;
}

psi::Psi<std::complex<double>>* psi_initializer::allocate()
{
    ModuleBase::timer::tick("psi_initializer", "allocate");
    /*
        WARNING: when basis_type = "pw", the variable GlobalV::NLOCAL will also be set, in this case, it is set to
        9 = 1 + 3 + 5, which is the maximal number of orbitals spd, I don't think it is reasonable
        The way of calculating GlobalC::ucell.natomwfc is, for each atom, read pswfc and for s, it is 1, for p, it is 3
        , then multiplied by the number of atoms, and then add them together.
    */
    if(this->psig != nullptr) delete this->psig;
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
            if(GlobalC::ucell.natomwfc >= GlobalV::NBANDS)
            {
                nbands_actual = GlobalC::ucell.natomwfc;
                this->nbands_complem = 0;
            }
            else
            {
                nbands_actual = GlobalV::NBANDS;
                this->nbands_complem = GlobalV::NBANDS - GlobalC::ucell.natomwfc;
            }
        }
        else if(GlobalV::init_wfc.substr(0, 3) == "nao")
        {
            /*
                previously GlobalV::NLOCAL is used here, however it is wrong. GlobalV::NLOCAL is fixed to 9*nat.
            */
            int nbands_local = 0;
            for(int it = 0; it < GlobalC::ucell.ntype; it++)
            {
                for(int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
                {
            /* FOR EVERY ATOM */
                    for(int l = 0; l < GlobalC::ucell.atoms[it].nwl + 1; l++)
                    {
            /* EVERY ZETA FOR (2l+1) ORBS, for NSPIN = 4, DOUBLE */
                        nbands_local += GlobalC::ucell.atoms[it].l_nchi[l]*(2*l+1) * GlobalV::NPOL;
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
    this->psig = new psi::Psi<std::complex<double>>(
        nkpts_actual, 
            nbands_actual, 
                nbasis_actual, 
                    this->pw_wfc->npwk);
    const size_t memory_cost = 
        nkpts_actual*
            nkpts_actual*
                this->pw_wfc->npwk_max * GlobalV::NPOL*
                    sizeof(std::complex<double>);
	std::cout << " MEMORY FOR PSI (MB)  : " << double(memory_cost)/1024.0/1024.0 << std::endl;
	ModuleBase::Memory::record("Psi_PW", memory_cost);
    ModuleBase::timer::tick("psi_initializer", "allocate");
    return psi_out;
}

void psi_initializer::write_psig() const
{
    for(int ik = 0; ik < this->psig->get_nk(); ik++)
    {
        std::string filename = "psig_"+std::to_string(ik);
        std::ofstream ofs_psig;
        ofs_psig.open(filename+"_kpt.out");
        ofs_psig << "N.B.: output data is complex, therefore every data will be enclosed by parenthesis." << std::endl;
        ofs_psig << "psig information" << std::endl;
        ofs_psig << "number of kpoints: " << this->psig->get_nk() << std::endl;
        ofs_psig << "number of bands: " << this->psig->get_nbands() << std::endl;
        ofs_psig << "number of planewaves: " << this->psig->get_nbasis() << std::endl;
        ofs_psig << "Calculation information" << std::endl;
        ofs_psig << "method of psi initialization: " << GlobalV::init_wfc << std::endl;
        ofs_psig << "method of diagonalization: " << GlobalV::KS_SOLVER << std::endl;
        this->psig->fix_k(ik);
        ofs_psig << "k point No. " << ik << std::endl;
        for(int iband = 0; iband < this->psig->get_nbands(); iband++)
        {
            ofs_psig << "energy band No. " << iband << std::endl;
            for(int ibasis = 0; ibasis < this->psig->get_nbasis(); ibasis++)
            {
                ofs_psig << std::setprecision(10) << std::fixed << (*(this->psig))(iband, ibasis) << " ";
            }
            ofs_psig << std::endl;
        }
        ofs_psig << std::endl;
        ofs_psig.close();
    }
}

void psi_initializer::write_psig(int ik) const
{
    std::string filename = "psig_"+std::to_string(ik);
    std::ofstream ofs_psig;
    ofs_psig.open(filename+"_kpt.out");
    ofs_psig << "N.B.: output data is complex, therefore every data will be enclosed by parenthesis." << std::endl;
    ofs_psig << "psig information" << std::endl;
    ofs_psig << "number of kpoints: " << this->psig->get_nk() << std::endl;
    ofs_psig << "number of bands: " << this->psig->get_nbands() << std::endl;
    ofs_psig << "number of planewaves: " << this->psig->get_nbasis() << std::endl;
    ofs_psig << "Calculation information" << std::endl;
    ofs_psig << "method of psi initialization: " << GlobalV::init_wfc << std::endl;
    ofs_psig << "method of diagonalization: " << GlobalV::KS_SOLVER << std::endl;
    this->psig->fix_k(ik);
    ofs_psig << "k point No. " << ik << std::endl;
    for(int iband = 0; iband < this->psig->get_nbands(); iband++)
    {
        ofs_psig << "energy band No. " << iband << std::endl;
        for(int ibasis = 0; ibasis < this->psig->get_nbasis(); ibasis++)
        {
            ofs_psig << std::setprecision(10) << std::fixed << (*(this->psig))(iband, ibasis) << " ";
        }
        ofs_psig << std::endl;
    }
    ofs_psig << std::endl;
    ofs_psig.close();
}

void psi_initializer::print_status(psi::Psi<std::complex<double>>& psi) const
{
    std::cout << "Current method: " << this->method << std::endl;
    std::cout << "Psi status:" << std::endl;
    std::cout << "  number of kpoints: " << psi.get_nk() << std::endl;
    std::cout << "  number of bands: " << psi.get_nbands() << std::endl;
    std::cout << "  number of planewaves: " << psi.get_nbasis() << std::endl;
}

void psi_initializer::random_t(std::complex<double>* psi, const int iw_start, const int iw_end, const int ik, const ModulePW::PW_Basis_K* wfc_basis)
{
    ModuleBase::timer::tick("psi_initializer", "random_t");
    assert(iw_start >= 0);
    const int ng = wfc_basis->npwk[ik];
#ifdef __MPI
    if (INPUT.pw_seed > 0) // qianrui add 2021-8-13
    {
        srand(unsigned(INPUT.pw_seed + GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL] + ik));
        const int nxy = wfc_basis->fftnxy;
        const int nz = wfc_basis->nz;
        const int nstnz = wfc_basis->nst*nz;

        double *stickrr = new double[nz];
        double *stickarg = new double[nz];
        double *tmprr = new double[nstnz];
        double *tmparg = new double[nstnz];
        for (int iw = iw_start; iw < iw_end; iw++)
        {   
            // get the starting memory address of iw band
            std::complex<double>* psi_slice = &(psi[iw * this->pw_wfc->npwk_max * GlobalV::NPOL]);
            int startig = 0;
            for(int ipol = 0 ; ipol < GlobalV::NPOL ; ++ipol)
            {
                    
                for(int ir=0; ir < nxy; ir++)
                {
                    if(wfc_basis->fftixy2ip[ir] < 0) continue;
                    if(GlobalV::RANK_IN_POOL==0)
                    {
                        for(int iz=0; iz<nz; iz++)
                        {
                            stickrr[ iz ] = std::rand()/double(RAND_MAX);
                            stickarg[ iz ] = std::rand()/double(RAND_MAX);
                        }
                    }
                    stick_to_pool(stickrr, ir, tmprr, wfc_basis);
                    stick_to_pool(stickarg, ir, tmparg, wfc_basis);
                }

                for (int ig = 0;ig < ng;ig++)
                {
                    const double rr = tmprr[wfc_basis->getigl2isz(ik,ig)];
                    const double arg= ModuleBase::TWO_PI * tmparg[wfc_basis->getigl2isz(ik,ig)];
                    const double gk2 = wfc_basis->getgk2(ik,ig);
                    psi_slice[ig+startig] = std::complex<double>(rr * cos(arg), rr * sin(arg)) / double(gk2 + 1.0);
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
        if (INPUT.pw_seed > 0) // qianrui add 2021-8-13
        {
            srand(unsigned(INPUT.pw_seed + ik));
        }
#endif
        for (int iw = iw_start ;iw < iw_end; iw++)
        {
            std::complex<double>* psi_slice = &(psi[iw * this->pw_wfc->npwk_max * GlobalV::NPOL]);
            for (int ig = 0; ig < ng; ig++)
            {
                const double rr = std::rand()/double(RAND_MAX); //qianrui add RAND_MAX
                const double arg= ModuleBase::TWO_PI * std::rand()/double(RAND_MAX);
                const double gk2 = wfc_basis->getgk2(ik,ig);
                psi_slice[ig] = std::complex<double>(rr * cos(arg), rr * sin(arg)) / double(gk2 + 1.0);
            }
            if(GlobalV::NPOL==2)
            {
                for (int ig = this->pw_wfc->npwk_max; ig < this->pw_wfc->npwk_max + ng; ig++)
                {
                    const double rr = std::rand()/double(RAND_MAX);
                    const double arg= ModuleBase::TWO_PI * std::rand()/double(RAND_MAX);
                    const double gk2 = wfc_basis->getgk2(ik,ig-this->pw_wfc->npwk_max);
                    psi_slice[ig] = std::complex<double>(rr * cos(arg), rr * sin(arg)) / double(gk2 + 1.0);
                }
            }
        }
#ifdef __MPI
    }
#endif
    ModuleBase::timer::tick("psi_initializer_random", "random_t");
}

#ifdef __MPI
void psi_initializer::stick_to_pool(double* stick, const int& ir, double* out, const ModulePW::PW_Basis_K* wfc_basis) const
{	
    ModuleBase::timer::tick("psi_initializer", "stick_to_pool");
	MPI_Status ierror;
    const int is = this->ixy2is[ir];
	const int ip = wfc_basis->fftixy2ip[ir];
    const int nz = wfc_basis->nz;

	if(ip == 0 && GlobalV::RANK_IN_POOL ==0)
	{
		for(int iz=0; iz<nz; iz++)
		{
			out[is*nz+iz] = stick[iz];
		}
	}
	else if(ip == GlobalV::RANK_IN_POOL )
	{
		MPI_Recv(stick, nz, MPI_DOUBLE, 0, ir, POOL_WORLD,&ierror);
		for(int iz=0; iz<nz; iz++)
		{
			out[is*nz+iz] = stick[iz];
		}
	}
	else if(GlobalV::RANK_IN_POOL==0)
	{
		MPI_Send(stick, nz, MPI_DOUBLE, ip, ir, POOL_WORLD);
	}

	return;	
    ModuleBase::timer::tick("psi_initializer", "stick_to_pool");
}
#endif
