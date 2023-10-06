#include "charge_mixing.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/inverse_matrix.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/element_elec_config.h"

Charge_Mixing::Charge_Mixing()
{
	rstep = 0;
	dstep = rstep - 1;//alway like this.
    initp = false;
	initb = false;
}

Charge_Mixing::~Charge_Mixing()
{
	if (initp)
    {
		deallocate_Pulay();
    }
	if (initb)
	{
		deallocate_Broyden();
	}
}

void Charge_Mixing::set_mixing
(
    const std::string &mixing_mode_in,
    const double &mixing_beta_in,
    const int &mixing_ndim_in,
	const double &mixing_gg0_in,
	const bool &mixing_tau_in
)
{
    this->mixing_mode = mixing_mode_in;
    this->mixing_beta = mixing_beta_in;
    this->mixing_ndim = mixing_ndim_in;
	this->mixing_gg0 = mixing_gg0_in; //mohan add 2014-09-27
	this->mixing_tau = mixing_tau_in;

	if(mixing_tau && mixing_mode == "broyden")
	{
		GlobalV::ofs_running << "Note : mixing_tau has only been implemented for plain and pulay mixing" << std::endl;
	}

    return;
}

void Charge_Mixing::need_auto_set()
{
	this->autoset = true;
}

void Charge_Mixing::auto_set(const double& bandgap_in, const UnitCell& ucell_)
{
	//auto set parameters once
	if(!this->autoset)
	{
		return;
	}
	else {
		this->autoset = false;
	}
	GlobalV::ofs_running<<"--------------AUTO-SET---------------"<<std::endl;
	//0.2 for metal and 0.7 for others
	if(bandgap_in * ModuleBase::Ry_to_eV < 1.0)
	{
		this->mixing_beta = 0.2;
	}
	else
	{
		this->mixing_beta = 0.7;
	}
	GlobalV::ofs_running<<"      Autoset mixing_beta to "<<this->mixing_beta<<std::endl;

	bool has_trans_metal = false;
	// find elements of cell
	for(int it=0;it<ucell_.ntype;it++)
	{
		if(ModuleBase::IsTransMetal.find(ucell_.atoms[it].ncpp.psd) != ModuleBase::IsTransMetal.end())
		{
			if(ModuleBase::IsTransMetal.at(ucell_.atoms[it].ncpp.psd)) 
			{
				has_trans_metal = true;
			}
		}
	}
	// auto set kerker mixing for trans metal system
	if(has_trans_metal)
	{
		this->mixing_gg0 = 1.5;
	}
	else
	{
		this->mixing_gg0 = 0.0;
	}
	GlobalV::ofs_running<<"      Autoset mixing_gg0 to "<<this->mixing_gg0<<std::endl;
	GlobalV::ofs_running<<"-------------------------------------"<<std::endl;
	//auto set for inhomogeneous system

}

double Charge_Mixing::get_drho(Charge* chr, const double nelec)
{
	for (int is=0; is<GlobalV::NSPIN; is++)
    {
		ModuleBase::GlobalFunc::NOTE("Perform FFT on rho(r) to obtain rho(G).");
        this->rhopw->real2recip(chr->rho[is], chr->rhog[is]);

		ModuleBase::GlobalFunc::NOTE("Perform FFT on rho_save(r) to obtain rho_save(G).");
        this->rhopw->real2recip(chr->rho_save[is], chr->rhog_save[is]);


		ModuleBase::GlobalFunc::NOTE("Calculate the charge difference between rho(G) and rho_save(G)");
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for (int ig=0; ig<this->rhopw->npw; ig++)
        {
           chr->rhog[is][ig] -= chr->rhog_save[is][ig];
        }
    }

	ModuleBase::GlobalFunc::NOTE("Calculate the norm of the Residual std::vector: < R[rho] | R[rho_save] >");
    double scf_thr = this->rhog_dot_product( chr->rhog, chr->rhog);
	
	if(GlobalV::test_charge)GlobalV::ofs_running << " scf_thr from rhog_dot_product is " << scf_thr << std::endl;
	
	// This is a temporary method to choose the type of scf_thr, in order to test different scf_thr_type.
	// Charge mixing will be refactored later.
	if(GlobalV::SCF_THR_TYPE == 1 && !GlobalV::test_charge)
	{
		return scf_thr;
	}

	// scf_thr calculated from real space.
	double scf_thr2 = 0.0;
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		if(is!=0 && is!=3 && GlobalV::DOMAG_Z)
		{
			continue;
		}
#ifdef _OPENMP
#pragma omp parallel for reduction(+:scf_thr2)
#endif
		for(int ir=0; ir<this->rhopw->nrxx; ir++)
		{
			scf_thr2 += std::abs( chr->rho[is][ir] - chr->rho_save[is][ir] );
		}
	}

	Parallel_Reduce::reduce_pool( scf_thr2 );
	assert( nelec != 0);
	assert( GlobalC::ucell.omega > 0);
	assert( this->rhopw->nxyz > 0);
	scf_thr2 *= GlobalC::ucell.omega / static_cast<double>( this->rhopw->nxyz );
	scf_thr2 /= nelec;
	if(GlobalV::test_charge)GlobalV::ofs_running << " scf_thr from real space grid is " << scf_thr2 << std::endl;

	// mohan add 2011-01-22
	//if(LINEAR_SCALING && LOCAL_BASIS) xiaohui modify 2013-09-01
	if(GlobalV::SCF_THR_TYPE == 2)
	{
		scf_thr = scf_thr2;	
	}
	return scf_thr;
}

void Charge_Mixing::mix_rho(const int &iter, Charge* chr)
{
    ModuleBase::TITLE("Charge_Mixing","mix_rho");
	ModuleBase::timer::tick("Charge", "mix_rho");

	// the charge before mixing.
	double **rho123 = new double*[GlobalV::NSPIN];
	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		rho123[is] = new double[this->rhopw->nrxx];
		if(is!=0 && is!=3 && GlobalV::DOMAG_Z)
		{
			ModuleBase::GlobalFunc::ZEROS(rho123[is], this->rhopw->nrxx);
			continue;
		} 
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
		for(int ir=0; ir<this->rhopw->nrxx; ++ir)
		{
			rho123[is][ir] = chr->rho[is][ir];
		}
	}

	double **kin_r123;
	if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
    {
		kin_r123 = new double*[GlobalV::NSPIN];
		for(int is=0; is<GlobalV::NSPIN; ++is)
		{
			kin_r123[is] = new double[this->rhopw->nrxx];
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
			for(int ir=0; ir<this->rhopw->nrxx; ++ir)
			{
				kin_r123[is][ir] = chr->kin_r[is][ir];
			}
		}
	}

	if ( this->mixing_mode == "plain")
    {
		this->plain_mixing(chr);
    }
    else if ( this->mixing_mode == "pulay")
    {
        this->Pulay_mixing(chr);
    }
    else if ( this->mixing_mode == "broyden")
    {
		this->Simplified_Broyden_mixing(iter, chr);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing","Not implemended yet,coming soon.");
    }

	// mohan add 2012-06-05
	// rho_save is the charge before mixing
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		if(is!=0 && is!=3 && GlobalV::DOMAG_Z) continue;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
		for(int ir=0; ir<this->rhopw->nrxx; ++ir)
		{
			chr->rho_save[is][ir] = rho123[is][ir];
		}
    }

	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		delete[] rho123[is];
	}
	delete[] rho123;

	if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			for(int ir=0; ir<this->rhopw->nrxx; ++ir)
			{
				chr->kin_r_save[is][ir] = kin_r123[is][ir];
			}
		}	

		for(int is=0; is<GlobalV::NSPIN; ++is)
		{
			delete[] kin_r123[is];
		}
		delete[] kin_r123;
	}


	if(new_e_iteration) new_e_iteration = false;

    ModuleBase::timer::tick("Charge","mix_rho");
    return;
}

void Charge_Mixing::plain_mixing(Charge* chr) const
{
	ModuleBase::timer::tick("Charge", "plain_mixing");
    // if mixing_beta == 1, each electron iteration,
    // use all new charge density,
    // on the contrary, if mixing_beta == 0,
    // no new charge will be generated!
    const double mix_old = 1 - mixing_beta;

//xiaohui add 2014-12-09
	for (int is=0; is<GlobalV::NSPIN; is++)
	{
		if(this->mixing_gg0 > 0.0)
		{
			double* Rrho = new double[this->rhopw->nrxx];
			std::complex<double> *kerpulay = new std::complex<double>[this->rhopw->npw];
			double* kerpulayR = new double[this->rhopw->nrxx];
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
			for(int ir=0; ir<this->rhopw->nrxx; ir++)
			{
				Rrho[ir] = chr->rho[is][ir] - chr->rho_save[is][ir];
			}
			this->rhopw->real2recip(Rrho, kerpulay);

			const double fac = this->mixing_gg0;
			const double gg0 = std::pow(fac * 0.529177 / GlobalC::ucell.tpiba, 2);
			double* filter_g = new double[this->rhopw->npw];
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
			for(int ig=0; ig<this->rhopw->npw; ig++)
			{
				double gg = this->rhopw->gg[ig];
				filter_g[ig] = std::max(gg / (gg + gg0), 0.1);

				kerpulay[ig] = (1 - filter_g[ig]) * kerpulay[ig];
			}
			this->rhopw->recip2real(kerpulay, kerpulayR);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 128)
#endif
			for(int ir=0; ir<this->rhopw->nrxx; ir++)
			{
				Rrho[ir] = Rrho[ir] - kerpulayR[ir];
				chr->rho[is][ir] = Rrho[ir] * mixing_beta + chr->rho_save[is][ir];
				chr->rho_save[is][ir] = chr->rho[is][ir];
			}

			delete[] Rrho;
			delete[] kerpulay;
			delete[] kerpulayR;
			delete[] filter_g;
		}
		else
		{
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
			for (int ir=0; ir<this->rhopw->nrxx; ir++)
			{
				chr->rho[is][ir] = chr->rho[is][ir]*mixing_beta + mix_old*chr->rho_save[is][ir];
				chr->rho_save[is][ir] = chr->rho[is][ir];
			}
		}

		if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
		{
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
			for (int ir=0; ir<this->rhopw->nrxx; ir++)
			{
				chr->kin_r[is][ir] = chr->kin_r[is][ir]*mixing_beta + mix_old*chr->kin_r_save[is][ir];
				chr->kin_r_save[is][ir] = chr->kin_r[is][ir];
			}
		}

	}
	ModuleBase::timer::tick("Charge", "plain_mixing");
    return;
}

double Charge_Mixing::rhog_dot_product(
	const std::complex<double>*const*const rhog1,
	const std::complex<double>*const*const rhog2
) const
{
    ModuleBase::TITLE("Charge_Mixing","rhog_dot_product");
	ModuleBase::timer::tick("Charge_Mixing","rhog_dot_product");
    static const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / GlobalC::ucell.tpiba2;
    static const double fac2 = ModuleBase::e2 * ModuleBase::FOUR_PI / (ModuleBase::TWO_PI * ModuleBase::TWO_PI);

    double sum = 0.0;
	
	auto part_of_noncolin = [&]()			// Peize Lin change goto to function at 2020.01.31
	{
		double sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum)
#endif
		for (int ig=0; ig<this->rhopw->npw; ++ig)
		{
			if(this->rhopw->gg[ig]<1e-8) continue;
			sum += ( conj( rhog1[0][ig] )* rhog2[0][ig] ).real() / this->rhopw->gg[ig];
		}
		sum *= fac;
		return sum;
	};

    switch ( GlobalV::NSPIN )
    {
	case 1:
		sum += part_of_noncolin();
		break;

	case 2:
		{
			// (1) First part of density error.
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum)
#endif
			for (int ig=0; ig<this->rhopw->npw; ++ig)
			{
				if(this->rhopw->gg[ig]<1e-8) continue;
				sum += ( conj( rhog1[0][ig]+rhog1[1][ig] ) * (rhog2[0][ig]+rhog2[1][ig]) ).real() / this->rhopw->gg[ig];
			}
			sum *= fac;

			if(GlobalV::GAMMA_ONLY_PW)
			{
				sum *= 2.0;
			}

			// (2) Second part of density error.
			// including |G|=0 term.
			double sum2 = 0.0;

			sum2 += fac2 * ( conj( rhog1[0][0]-rhog1[1][0] ) * ( rhog2[0][0]-rhog2[1][0] ) ).real();

			double mag = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:mag)
#endif
			for (int ig=0; ig<this->rhopw->npw; ig++)
			{
				mag += ( conj( rhog1[0][ig]-rhog1[1][ig] ) * ( rhog2[0][ig]-rhog2[1][ig] ) ).real();
			}
			mag *= fac2;

			//if(GlobalV::GAMMA_ONLY_PW);
			if(GlobalV::GAMMA_ONLY_PW)			// Peize Lin delete ; 2020.01.31
			{
				mag *= 2.0;
			}

			//std::cout << " sum=" << sum << " mag=" << mag << std::endl;
			sum2 += mag;
			sum += sum2;
			break;
		}
	case 4:
		// non-collinear spin, added by zhengdy
		if(!GlobalV::DOMAG&&!GlobalV::DOMAG_Z)
			sum += part_of_noncolin();
		else
		{
			//another part with magnetization
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum)
#endif
			for (int ig=0; ig<this->rhopw->npw; ig++)
			{
				if(ig==this->rhopw->ig_gge0) continue;
				sum += ( conj( rhog1[0][ig] )* rhog2[0][ig] ).real() / this->rhopw->gg[ig];
			}
			sum *= fac;
			const int ig0 = this->rhopw->ig_gge0;
			if(ig0 > 0)
			{
				sum += fac2 * ((conj( rhog1[1][ig0])*rhog2[1][ig0]).real() +
					(conj( rhog1[2][ig0])*rhog2[2][ig0]).real() +
					(conj( rhog1[3][ig0])*rhog2[3][ig0]).real());
			}
			double fac3 = fac2;
			if(GlobalV::GAMMA_ONLY_PW)
			{
				fac3 *= 2.0;
			}
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum)
#endif
			for (int ig=0; ig<this->rhopw->npw; ig++)
			{
				if(ig == ig0) continue;
				sum += fac3 * ((conj( rhog1[1][ig])*rhog2[1][ig]).real() +
					(conj( rhog1[2][ig])*rhog2[2][ig]).real() +
					(conj( rhog1[3][ig])*rhog2[3][ig]).real());
			}
		}
		break;
    }

    Parallel_Reduce::reduce_pool( sum );

	ModuleBase::timer::tick("Charge_Mixing","rhog_dot_product");

	sum *= GlobalC::ucell.omega * 0.5;

    return sum;
}

