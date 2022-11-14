#include "charge_mixing.h"
#include "global.h"
#include "../module_base/inverse_matrix.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"

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
	const double &mixing_gg0_in
)
{
    this->mixing_mode = mixing_mode_in;
    this->mixing_beta = mixing_beta_in;
    this->mixing_ndim = mixing_ndim_in;
	this->mixing_gg0 = mixing_gg0_in; //mohan add 2014-09-27

    return;
}

double Charge_Mixing::get_drho(double** rho, double** rho_save,
	std::complex<double>** rhog, std::complex<double>** rhog_save, const double nelec)
{
	double scf_thr;
	for (int is=0; is<GlobalV::NSPIN; is++)
    {
		ModuleBase::GlobalFunc::NOTE("Perform FFT on rho(r) to obtain rho(G).");
        GlobalC::rhopw->real2recip(rho[is], rhog[is]);

		ModuleBase::GlobalFunc::NOTE("Perform FFT on rho_save(r) to obtain rho_save(G).");
        GlobalC::rhopw->real2recip(rho_save[is], rhog_save[is]);


		ModuleBase::GlobalFunc::NOTE("Calculate the charge difference between rho(G) and rho_save(G)");
        for (int ig=0; ig<GlobalC::rhopw->npw; ig++)
        {
            rhog[is][ig] -= rhog_save[is][ig];
        }

    }

	ModuleBase::GlobalFunc::NOTE("Calculate the norm of the Residual std::vector: < R[rho] | R[rho_save] >");
    scf_thr = this->rhog_dot_product( rhog, rhog);
	
	if(GlobalV::test_charge)GlobalV::ofs_running << " scf_thr from rhog_dot_product is " << scf_thr << std::endl;

	// scf_thr calculated from real space.
	double scf_thr2 = 0.0;
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
		{
			scf_thr2 += abs( rho[is][ir] - rho_save[is][ir] );
		}
	}

	Parallel_Reduce::reduce_double_pool( scf_thr2 );
	assert( nelec != 0);
	assert( GlobalC::ucell.omega > 0);
	assert( GlobalC::rhopw->nxyz > 0);
	scf_thr2 *= GlobalC::ucell.omega / static_cast<double>( GlobalC::rhopw->nxyz );
	scf_thr2 /= nelec;
	if(GlobalV::test_charge)GlobalV::ofs_running << " scf_thr from real space grid is " << scf_thr2 << std::endl;

	// mohan add 2011-01-22
	//if(LINEAR_SCALING && LOCAL_BASIS) xiaohui modify 2013-09-01
	if(GlobalV::BASIS_TYPE=="lcao" )
	{
		scf_thr = scf_thr2;	
	}
	return scf_thr;
}

void Charge_Mixing::mix_rho
(
    const int &iter,
	double** rho,
	double** rho_save,
	std::complex<double>** rhog,
	std::complex<double>** rhog_save
)
{
    ModuleBase::TITLE("Charge_Mixing","mix_rho");
	ModuleBase::timer::tick("Charge", "mix_rho");

	// the charge before mixing.
	double **rho123 = new double*[GlobalV::NSPIN];
	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		rho123[is] = new double[GlobalC::rhopw->nrxx];
		ModuleBase::GlobalFunc::ZEROS(rho123[is], GlobalC::rhopw->nrxx);
		for(int ir=0; ir<GlobalC::rhopw->nrxx; ++ir)
		{
			rho123[is][ir] = rho[is][ir];
		}
	}
	
	if ( this->mixing_mode == "plain")
    {
        // calculate mixing change, and save it in rho1.
        for (int is=0; is<GlobalV::NSPIN; is++)
        {
            this->plain_mixing( rho[is], rho_save[is]);
        }
    }
    else if ( this->mixing_mode == "pulay")
    {
        this->Pulay_mixing(rho, rho_save);
    }
    else if ( this->mixing_mode == "broyden")
    {
		this->Simplified_Broyden_mixing(iter, rho, rho_save, rhog, rhog_save);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing","Not implemended yet,coming soon.");
    }

	// mohan add 2012-06-05
	// rho_save is the charge before mixing
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int ir=0; ir<GlobalC::rhopw->nrxx; ++ir)
		{
			rho_save[is][ir] = rho123[is][ir];
		}
    }

	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		delete[] rho123[is];
	}
	delete[] rho123;

	if(new_e_iteration) new_e_iteration = false;

    ModuleBase::timer::tick("Charge","mix_rho");
    return;
}

void Charge_Mixing::plain_mixing( double *rho, double *rho_save_in ) const
{
    // if mixing_beta == 1, each electron iteration,
    // use all new charge density,
    // on the contrary, if mixing_beta == 0,
    // no new charge will be generated!
    const double mix_old = 1 - mixing_beta;

//xiaohui add 2014-12-09
	if(this->mixing_gg0 > 0.0)
	{
		double* Rrho = new double[GlobalC::rhopw->nrxx];
		std::complex<double> *kerpulay = new std::complex<double>[GlobalC::rhopw->npw];
		double* kerpulayR = new double[GlobalC::rhopw->nrxx];

		for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
		{
			Rrho[ir] = rho[ir] - rho_save_in[ir];
		}
		GlobalC::rhopw->real2recip(Rrho, kerpulay);

		const double fac = this->mixing_gg0;
		const double gg0 = std::pow(fac * 0.529177 / GlobalC::ucell.tpiba, 2);
		double* filter_g = new double[GlobalC::rhopw->npw];
		for(int ig=0; ig<GlobalC::rhopw->npw; ig++)
		{
			double gg = GlobalC::rhopw->gg[ig];
			filter_g[ig] = max(gg / (gg + gg0), 0.1);

			kerpulay[ig] = (1 - filter_g[ig]) * kerpulay[ig];
		}
		GlobalC::rhopw->recip2real(kerpulay, kerpulayR);

		for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
		{
			Rrho[ir] = Rrho[ir] - kerpulayR[ir];
			rho[ir] = Rrho[ir] * mixing_beta + rho_save_in[ir];
		}

		delete[] Rrho;
		delete[] kerpulay;
		delete[] kerpulayR;
		delete[] filter_g;
	}
	else
	{
		for (int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
		{
			rho[ir] = rho[ir]*mixing_beta + mix_old*rho_save_in[ir];
		}
	}

	ModuleBase::GlobalFunc::DCOPY( rho, rho_save_in, GlobalC::rhopw->nrxx);

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
		for (int ig=0; ig<GlobalC::rhopw->npw; ++ig)
		{
			if(GlobalC::rhopw->gg[ig]<1e-8) continue;
			sum += ( conj( rhog1[0][ig] )* rhog2[0][ig] ).real() / GlobalC::rhopw->gg[ig];
		}
		sum *= fac;
	};

    switch ( GlobalV::NSPIN )
    {
	case 1:
		part_of_noncolin();
		break;

	case 2:
		{
			// (1) First part of density error.
			for (int ig=0; ig<GlobalC::rhopw->npw; ++ig)
			{
				if(GlobalC::rhopw->gg[ig]<1e-8) continue;
				sum += ( conj( rhog1[0][ig]+rhog1[1][ig] ) * (rhog2[0][ig]+rhog2[1][ig]) ).real() / GlobalC::rhopw->gg[ig];
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
			for (int ig=0; ig<GlobalC::rhopw->npw; ig++)
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
			part_of_noncolin();
		else
		{
			//another part with magnetization
			for (int ig=0; ig<GlobalC::rhopw->npw; ig++)
			{
				if(ig==GlobalC::rhopw->ig_gge0) continue;
				sum += ( conj( rhog1[0][ig] )* rhog2[0][ig] ).real() / GlobalC::rhopw->gg[ig];
			}
			sum *= fac;
			const int ig0 = GlobalC::rhopw->ig_gge0;
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
			for (int ig=0; ig<GlobalC::rhopw->npw; ig++)
			{
				if(ig == ig0) continue;
				sum += fac3 * ((conj( rhog1[1][ig])*rhog2[1][ig]).real() +
					(conj( rhog1[2][ig])*rhog2[2][ig]).real() +
					(conj( rhog1[3][ig])*rhog2[3][ig]).real());
			}
		}
		break;
    }

    Parallel_Reduce::reduce_double_pool( sum );

	ModuleBase::timer::tick("Charge_Mixing","rhog_dot_product");

	sum *= GlobalC::ucell.omega * 0.5;

    return sum;
}

