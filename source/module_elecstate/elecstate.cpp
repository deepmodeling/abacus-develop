#include "elecstate.h"

#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_base/memory.h"
#include "src_parallel/parallel_reduce.h"
#include "src_pw/global.h"
#include "src_pw/occupy.h"

namespace elecstate
{

const double* ElecState::getRho(int spin) const
{
    // hamilt::MatrixBlock<double> temp{&(this->charge->rho[spin][0]), 1, this->charge->nrxx}; //
    // this->chr->get_nspin(), this->chr->get_nrxx()};
    return &(this->charge->rho[spin][0]);
}

void ElecState::fixed_weights(const double * const ocp_kb)
{
    for (int ik = 0; ik < this->wg.nr; ik++)
    {
        for (int ib = 0; ib < this->wg.nc; ib++)
        {
            this->wg(ik, ib) = ocp_kb[ik * this->wg.nc + ib];
        }
    }
    this->skip_weights = true;
}

void ElecState::init_nelec_spin()
{
    this->nelec_spin.resize(GlobalV::NSPIN);
    if(GlobalV::NSPIN==2)
    {
        //in fact, when TWO_EFERMI(nupdown in INPUT is not 0.0), nelec_spin will be fixed.
        this->nelec_spin[0] = (GlobalV::nelec + GlobalV::nupdown) / 2.0 ;
        this->nelec_spin[1] = (GlobalV::nelec - GlobalV::nupdown) / 2.0 ;
    }
}

void ElecState::calculate_weights()
{
    ModuleBase::TITLE("ElecState", "calculate_weights");
    if(this->skip_weights)
    {
        return;
    }

    double** ekb_tmp = new double*[this->ekb.nr];
    for (int i = 0; i < this->ekb.nr; ++i)
    {
        ekb_tmp[i] = &(this->ekb(i, 0));
    }
    int nbands = this->ekb.nc;
    int nks = this->ekb.nr;

    if (!Occupy::use_gaussian_broadening && !Occupy::fixed_occupations)
    {
        if (GlobalV::TWO_EFERMI)
        {
            Occupy::iweights(nks,
                             this->klist->wk,
                             nbands,
                             this->nelec_spin[0],
                             ekb_tmp,
                             GlobalC::en.ef_up,
                             this->wg,
                             0,
                             this->klist->isk);
            Occupy::iweights(nks,
                             this->klist->wk,
                             nbands,
                             this->nelec_spin[1],
                             ekb_tmp,
                             GlobalC::en.ef_dw,
                             this->wg,
                             1,
                             this->klist->isk);
            // ef = ( ef_up + ef_dw ) / 2.0_dp need??? mohan add 2012-04-16
        }
        else
        {
            // -1 means don't need to consider spin.
            Occupy::iweights(nks,
                             this->klist->wk,
                             nbands,
                             GlobalV::nelec,
                             ekb_tmp,
                             this->ef,
                             this->wg,
                             -1,
                             this->klist->isk);
        }
    }
    else if (Occupy::use_gaussian_broadening)
    {
        if (GlobalV::TWO_EFERMI)
        {
            double demet_up = 0.0;
            double demet_dw = 0.0;
            Occupy::gweights(nks,
                             this->klist->wk,
                             nbands,
                             this->nelec_spin[0],
                             Occupy::gaussian_parameter,
                             Occupy::gaussian_type,
                             ekb_tmp,
                             GlobalC::en.ef_up,
                             demet_up,
                             this->wg,
                             0,
                             this->klist->isk);
            Occupy::gweights(nks,
                             this->klist->wk,
                             nbands,
                             this->nelec_spin[1],
                             Occupy::gaussian_parameter,
                             Occupy::gaussian_type,
                             ekb_tmp,
                             GlobalC::en.ef_dw,
                             demet_dw,
                             this->wg,
                             1,
                             this->klist->isk);
            this->demet = demet_up + demet_dw;
        }
        else
        {
            // -1 means is no related to spin.
            Occupy::gweights(nks,
                             this->klist->wk,
                             nbands,
                             GlobalV::nelec,
                             Occupy::gaussian_parameter,
                             Occupy::gaussian_type,
                             ekb_tmp,
                             this->ef,
                             this->demet,
                             this->wg,
                             -1,
                             this->klist->isk);
        }

        // qianrui fix a bug on 2021-7-21
        Parallel_Reduce::reduce_double_allpool(this->demet);
    }
    else if (Occupy::fixed_occupations)
    {
        // fix occupations need nelup and neldw.
        // mohan add 2011-04-03
        this->ef = -1.0e+20;
        for (int ik = 0; ik < nks; ik++)
        {
            for (int ibnd = 0; ibnd < nbands; ibnd++)
            {
                if (this->wg(ik, ibnd) > 0.0)
                {
                    this->ef = std::max(this->ef, ekb_tmp[ik][ibnd]);
                }
            }
        }
    }

    if (GlobalV::TWO_EFERMI)
    {
        Parallel_Reduce::gather_max_double_all(GlobalC::en.ef_up);
        Parallel_Reduce::gather_max_double_all(GlobalC::en.ef_dw);
    }
    else
    {
        double ebotom = ekb_tmp[0][0];
        double etop = ekb_tmp[0][0];
        for (int ik = 0; ik < nks; ik++)
        {
            for (int ib = 0; ib < nbands; ib++)
            {
                ebotom = min(ebotom, ekb_tmp[ik][ib]);
                etop = max(etop, ekb_tmp[ik][ib]);
            }
        }

        // parallel
        Parallel_Reduce::gather_max_double_all(this->ef);
        Parallel_Reduce::gather_max_double_all(etop);
        Parallel_Reduce::gather_min_double_all(ebotom);

        // not parallel yet!
        //		OUT(GlobalV::ofs_running,"Top    Energy (eV)", etop * ModuleBase::Ry_to_eV);
        //      OUT(GlobalV::ofs_running,"Fermi  Energy (eV)", this->ef * ModuleBase::Ry_to_eV);
        //		OUT(GlobalV::ofs_running,"Bottom Energy (eV)", ebotom * ModuleBase::Ry_to_eV);
        //		OUT(GlobalV::ofs_running,"Range  Energy (eV)", etop-ebotom * ModuleBase::Ry_to_eV);
    }

    delete[] ekb_tmp;

    return;
}

void ElecState::calEBand()
{
    ModuleBase::TITLE("ElecState", "calEBand");
    // calculate ebands using wg and ekb
    this->eband = 0.0;
    for (int ik = 0; ik < this->ekb.nr; ++ik)
    {
        for (int ibnd = 0; ibnd < this->ekb.nc; ibnd++)
        {
            this->eband += this->ekb(ik, ibnd) * this->wg(ik, ibnd);
        }
    }
    if (GlobalV::KPAR != 1 && GlobalV::ESOLVER_TYPE != "sdft")
    {
        //==================================
        // Reduce all the Energy in each cpu
        //==================================
        this->eband /= GlobalV::NPROC_IN_POOL;
        Parallel_Reduce::reduce_double_all(this->eband);
    }
    return;
}

void ElecState::print_band(const int& ik, const int& printe, const int& iter)
{
    //check the band energy.
    bool wrong = false;
	for(int ib=0; ib<GlobalV::NBANDS; ++ib)
	{
		if( abs( this->ekb(ik, ib) ) > 1.0e10)
		{
			GlobalV::ofs_warning << " ik=" << ik+1 << " ib=" << ib+1 << " " << this->ekb(ik, ib) << " Ry" << std::endl;
			wrong = true;
		}
	}
	if(wrong)
    {
        ModuleBase::WARNING_QUIT("Threshold_Elec::print_eigenvalue","Eigenvalues are too large!");
    }



	if(GlobalV::MY_RANK==0)
	{
        if( printe>0 && ((iter+1) % printe == 0))
        {            
            GlobalV::ofs_running << std::setprecision(6);
            GlobalV::ofs_running << " Energy (eV) & Occupations  for spin=" << GlobalV::CURRENT_SPIN+1 << " K-point=" << ik+1 << std::endl;
            GlobalV::ofs_running << std::setiosflags(ios::showpoint);
            for(int ib=0;ib<GlobalV::NBANDS;ib++)
            {
                GlobalV::ofs_running << " "<< std::setw(6) << ib+1  
                    << std::setw(15) << this->ekb(ik, ib) * ModuleBase::Ry_to_eV;
                // for the first electron iteration, we don't have the energy
                // spectrum, so we can't get the occupations. 
                GlobalV::ofs_running << std::setw(15) << this->wg(ik,ib);
                GlobalV::ofs_running << std::endl;
            }
        }
	}
	return;
}

void ElecState::init_scf(const int istep, const ModuleBase::ComplexMatrix& strucfac)
{
    //---------Charge part-----------------
    // core correction potential.
    this->charge->set_rho_core(strucfac);

    //--------------------------------------------------------------------
    // (2) other effective potentials need charge density,
    // choose charge density from ionic step 0.
    //--------------------------------------------------------------------
    if (istep == 0)
    {
        this->charge->init_rho();
    }

    // renormalize the charge density
    this->charge->renormalize_rho();

    //---------Potential part--------------
    this->pot->init_pot(istep, this->charge);
}

void ElecState::init_ks(
    Charge *chg_in, // pointer for class Charge
    const K_Vectors *klist_in,
    int nk_in
)
{
    this->charge = chg_in;
    this->klist = klist_in;
    //init nelec_spin with nelec and nupdown
    this->init_nelec_spin();
    //autoset and check GlobalV::NBANDS, nelec_spin is used when NSPIN==2
    this->cal_nbands();
    //initialize ekb and wg
    this->ekb.create(nk_in, GlobalV::NBANDS);
    this->wg.create(nk_in, GlobalV::NBANDS);
}

void ElecState::cal_nbands()
{
    if ( GlobalV::ESOLVER_TYPE == "sdft" ) //qianrui 2021-2-20
	{
        return;
    }
    //=======================================
	// calculate number of bands (setup.f90)
	//=======================================
	double occupied_bands = static_cast<double>(GlobalV::nelec/ModuleBase::DEGSPIN);	
	if(GlobalV::LSPINORB==1) occupied_bands = static_cast<double>(GlobalV::nelec);

	if( (occupied_bands - std::floor(occupied_bands)) > 0.0 )
	{
		occupied_bands = std::floor(occupied_bands) + 1.0; //mohan fix 2012-04-16
	}

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"occupied bands",occupied_bands);

	// mohan add 2010-09-04
    //std::cout << "nbands(GlobalC::ucell) = " <<GlobalV::NBANDS <<std::endl;
	if(GlobalV::NBANDS==occupied_bands)
	{
		if( Occupy::gauss() )
		{
			ModuleBase::WARNING_QUIT("ElecState::cal_nbands","for smearing, num. of bands > num. of occupied bands");
		}
	}
	
	if(GlobalV::NBANDS == 0)
	{
		if(GlobalV::NSPIN == 1)
		{
			const int nbands1 = static_cast<int>(occupied_bands) + 10;
			const int nbands2 = static_cast<int>(1.2 * occupied_bands) + 1;
			GlobalV::NBANDS = std::max(nbands1, nbands2);
			if(GlobalV::BASIS_TYPE!="pw") GlobalV::NBANDS = std::min(GlobalV::NBANDS, GlobalV::NLOCAL);
		}
		else if (GlobalV::NSPIN == 4)
		{
			const int nbands3 = GlobalV::nelec + 20;
			const int nbands4 = static_cast<int>(1.2 * GlobalV::nelec) + 1;
			GlobalV::NBANDS = std::max(nbands3, nbands4);
			if(GlobalV::BASIS_TYPE!="pw") GlobalV::NBANDS = std::min(GlobalV::NBANDS, GlobalV::NLOCAL);
		}
        else if (GlobalV::NSPIN == 2)
        {
            const double max_occ = std::max(this->nelec_spin[0], this->nelec_spin[1]);
            const int nbands3 = static_cast<int>(max_occ) + 11;
			const int nbands4 = static_cast<int>(1.2 * max_occ) + 1;
            GlobalV::NBANDS = std::max(nbands3, nbands4);
            if(GlobalV::BASIS_TYPE!="pw") GlobalV::NBANDS = std::min(GlobalV::NBANDS, GlobalV::NLOCAL);
        }
		ModuleBase::GlobalFunc::AUTO_SET("NBANDS",GlobalV::NBANDS);
	}
	//else if ( GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="md" || GlobalV::CALCULATION=="relax") //pengfei 2014-10-13
	else
	{
		if(GlobalV::NBANDS < occupied_bands) ModuleBase::WARNING_QUIT("unitcell","Too few bands!");
        if(GlobalV::NSPIN==2)
        {
            if(GlobalV::NBANDS < this->nelec_spin[0] ) 
            {
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nelec_up", this->nelec_spin[0]);
                ModuleBase::WARNING_QUIT("ElecState::cal_nbands","Too few spin up bands!");
            }
            if(GlobalV::NBANDS < this->nelec_spin[1] )
            {
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nelec_down", this->nelec_spin[1]);
                ModuleBase::WARNING_QUIT("ElecState::cal_nbands","Too few spin down bands!");
            }
        }
	}

	// mohan update 2021-02-19
    // mohan add 2011-01-5
    if(GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw")
    {
        if( GlobalV::NBANDS > GlobalV::NLOCAL )
        {
            ModuleBase::WARNING_QUIT("ElecState::cal_nbandsc","NLOCAL < NBANDS");
        }
        else
        {
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"NLOCAL",GlobalV::NLOCAL);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"NBANDS",GlobalV::NBANDS);
        }
    }

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"NBANDS",GlobalV::NBANDS);
}

} // namespace elecstate
