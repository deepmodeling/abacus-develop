#include "charge_mixing.h"
#include "global.h"
#include "../module_base/inverse_matrix.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"

Charge_Mixing::Charge_Mixing(){}
Charge_Mixing::~Charge_Mixing(){}

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

//    if (mixing_mode != "plain") // "TF","local-TF","potential"
//    {
//        ModuleBase::WARNING_QUIT("set_mixing","only plain mixing availabel.");
//    }

    return;
}

    
// Fourier transform of rho(g) to rho(r) in real space.
void Charge_Mixing::set_rhor(std::complex<double> *rhog, double *rho)const
{
    if (GlobalV::test_charge)ModuleBase::TITLE("Charge_Mixing","set_rhor");
    for (int is=0; is < GlobalV::NSPIN; is++)
    {
		GlobalC::UFFT.ToRealSpace(rhog, rho);
    }
    return;
}


void Charge_Mixing::set_rhog(const double *rho_in, std::complex<double> *rhog_in)const
{
    if (GlobalV::test_charge)ModuleBase::TITLE("Charge_Mixing","set_rhog");
	GlobalC::UFFT.ToReciSpace(rho_in, rhog_in);
    return;
}


void Charge_Mixing::plain_mixing( double *rho, double *rho_save_in ) const
{
    // if mixing_beta == 1, each electron iteration,
    // use all new charge density,
    // on the contrary, if mixing_beta == 0,
    // no new charge will be generated!
    const double mix_old = 1 - mixing_beta;
	
//	this->check_ne(rho);
//	this->check_ne(rho_save_in);

    // in real space
	// mohan modify 2010-02-05
	// after mixing, the charge density become 
	// the input charge density of next iteration.
    //double* rho_tmp = new double[GlobalC::pw.nrxx];
    //ModuleBase::GlobalFunc::DCOPY( rho, rho_tmp, GlobalC::pw.nrxx);

//xiaohui add 2014-12-09
	if(this->mixing_gg0 > 0.0)
	{
		double* Rrho = new double[GlobalC::pw.nrxx];
		std::complex<double> *kerpulay = new std::complex<double>[GlobalC::pw.ngmc];
		double* kerpulayR = new double[GlobalC::pw.nrxx];

		for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
		{
			Rrho[ir] = rho[ir] - rho_save_in[ir];
		}
		set_rhog(Rrho, kerpulay);

		const double fac = this->mixing_gg0;
		const double gg0 = std::pow(fac * 0.529177 / GlobalC::ucell.tpiba, 2);
		double* filter_g = new double[GlobalC::pw.ngmc];
		for(int ig=0; ig<GlobalC::pw.ngmc; ig++)
		{
			double gg = GlobalC::pw.get_NormG_cartesian(ig);
			filter_g[ig] = max(gg / (gg + gg0), 0.1);

			kerpulay[ig] = (1 - filter_g[ig]) * kerpulay[ig];
		}
		set_rhor(kerpulay, kerpulayR);

		for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
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
		for (int ir=0; ir<GlobalC::pw.nrxx; ir++)
		{
			rho[ir] = rho[ir]*mixing_beta + mix_old*rho_save_in[ir];
		}
	}

	ModuleBase::GlobalFunc::DCOPY( rho, rho_save_in, GlobalC::pw.nrxx);
//    delete[] rho_tmp;

    return;
}


void Charge_Mixing::Kerker_mixing( double *rho, const std::complex<double> *residual_g, double *rho_save)
{
//  ModuleBase::TITLE("Charge_Mixing","Kerker");
    ModuleBase::timer::tick("Charge_Mixing","Kerker");

//	std::cout << " here is Kerker" << std::endl;
//	this->check_ne(rho);
//	this->check_ne(rho_save);

    // (1) do kerker mixing in reciprocal space.
    std::complex<double> *rhog = new std::complex<double>[GlobalC::pw.ngmc];
    ModuleBase::GlobalFunc::ZEROS(rhog, GlobalC::pw.ngmc);

	// mohan modify 2010-02-03, rhog should store the old
	// charge density. " rhog = FFT^{-1}(rho_save) "
    this->set_rhog(rho_save, rhog);

    // (2) set up filter
    //const double a = 0.8; // suggested by VASP.

	// mohan fixed bug 2010/03/25
	// suggested by VASP, 1.5(angstrom^-1) is always satisfied.
    const double gg0 = std::pow(1.5 * 0.529177 / GlobalC::ucell.tpiba, 2);
    double *filter_g = new double[GlobalC::pw.ngmc];
    for (int ig=0; ig<GlobalC::pw.ngmc; ig++)
    {
        double gg = GlobalC::pw.get_NormG_cartesian(ig);
//      filter_g[ig] = a * gg / (gg+gg0);
		filter_g[ig] = mixing_beta * gg / (gg+gg0);//mohan modify 2010/03/25
    }

    // (3)
    for (int ig=0; ig<GlobalC::pw.ngmc; ig++)
    {
        rhog[ig] += filter_g[ig] * residual_g[ig];
    }

    // (4) transform rhog from G-space to real space.
	// get the new charge. " FFT(rhog) = rho "
    this->set_rhor( rhog, rho);

    // (5)
	// mohan change the order of (4) (5), 2010-02-05
    ModuleBase::GlobalFunc::DCOPY(rho, rho_save, GlobalC::pw.nrxx);

    //this->renormalize_rho();


    delete[] filter_g;
    delete[] rhog;
    ModuleBase::timer::tick("Charge_Mixing","Kerker");
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
		for (int ig=GlobalC::pw.gstart; ig<GlobalC::pw.ngmc; ig++)
		{
			sum += ( conj( rhog1[0][ig] )* rhog2[0][ig] ).real() / GlobalC::pw.gg[ig];
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
			for (int ig=GlobalC::pw.gstart; ig<GlobalC::pw.ngmc; ig++)
			{
				sum += ( conj( rhog1[0][ig]+rhog1[1][ig] ) * (rhog2[0][ig]+rhog2[1][ig]) ).real() / GlobalC::pw.gg[ig];
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
			for (int ig=0; ig<GlobalC::pw.ngmc; ig++)
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
			for (int ig=GlobalC::pw.gstart; ig<GlobalC::pw.ngmc; ig++)
			{
				sum += ( conj( rhog1[0][ig] )* rhog2[0][ig] ).real() / GlobalC::pw.gg[ig];
			}
			sum *= fac;
			if(GlobalC::pw.gstart == 2)
			{
				sum += fac2 * ((conj( rhog1[1][0])*rhog2[1][0]).real() +
					(conj( rhog1[2][0])*rhog2[2][0]).real() +
					(conj( rhog1[3][0])*rhog2[3][0]).real());
			}
			double fac3 = fac2;
			if(GlobalV::GAMMA_ONLY_PW)
			{
				fac3 *= 2.0;
			}
			for (int ig=GlobalC::pw.gstart; ig<GlobalC::pw.ngmc; ig++)
			{
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

