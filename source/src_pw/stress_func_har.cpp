#include "./stress_func.h"
#include "./myfunc.h"
#include "./H_Hartree_pw.h"
#include "../module_base/timer.h"
#include "global.h"

//calculate the Hartree part in PW or LCAO base
void Stress_Func::stress_har(ModuleBase::matrix& sigma, const bool is_pw)
{
	ModuleBase::timer::tick("Stress_Func","stress_har");
	double shart;

	std::complex<double> *Porter = GlobalC::UFFT.porter;

	//  Hartree potential VH(r) from n(r)
	ModuleBase::GlobalFunc::ZEROS( Porter, GlobalC::pw.nrxx );
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for (int ir=0; ir<GlobalC::pw.nrxx; ir++)
		{
			Porter[ir] += std::complex<double>( GlobalC::CHR.rho[is][ir], 0.0 );
		}
	}
	//=============================
	//  bring rho (aux) to G space
	//=============================
	GlobalC::pw.FFT_chg.FFT3D(Porter, -1);

	std::complex<double> *psic = new std::complex<double> [GlobalC::pw.nrxx];
	double *psic0 = new double[GlobalC::pw.nrxx];
	ModuleBase::GlobalFunc::ZEROS( psic0, GlobalC::pw.nrxx);
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		daxpy (GlobalC::pw.nrxx, 1.0, GlobalC::CHR.rho[is], 1, psic0, 2);
		for (int ir=0; ir<GlobalC::pw.nrxx; ir++)
		{
			psic[ir] = std::complex<double>(psic0[ir], 0.0);
		}
	}

	GlobalC::pw.FFT_chg.FFT3D(psic, -1) ;

	std::complex<double> *vh_g  = new std::complex<double>[GlobalC::pw.ngmc];
	ModuleBase::GlobalFunc::ZEROS(vh_g, GlobalC::pw.ngmc);

//	double ehart=0;
	for (int ig = GlobalC::pw.gstart; ig<GlobalC::pw.ngmc; ig++)
	{
		const int j = GlobalC::pw.ig2fftc[ig];
		//const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (GlobalC::ucell.tpiba2 * GlobalC::pw.gg [ig]);
		//ehart += ( conj( Porter[j] ) * Porter[j] ).real() * fac;
		//vh_g[ig] = fac * Porter[j];
		shart= ( conj( Porter[j] ) * Porter[j] ).real()/(GlobalC::ucell.tpiba2 * GlobalC::pw.gg [ig]);
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<l+1;m++)
			{
				sigma(l, m) += shart * 2 * GlobalC::pw.get_G_cartesian_projection(ig, l) * GlobalC::pw.get_G_cartesian_projection(ig, m) / GlobalC::pw.gg[ig];
			}
		}
	}
	//	Parallel_Reduce::reduce_double_pool( GlobalC::en.ehart );
	//	ehart *= 0.5 * GlobalC::ucell.omega;
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<l+1;m++)
		{
			Parallel_Reduce::reduce_double_pool( sigma(l,m) );
		}
	}

//        Parallel_Reduce::reduce_double_pool( ehart );
//        ehart *= 0.5 * GlobalC::ucell.omega;
        //psic(:)=(0.0,0.0)
	if(is_pw&&INPUT.gamma_only)
	{
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<3;m++)
			{
				sigma(l,m) *= ModuleBase::e2 * ModuleBase::FOUR_PI;
			}
		}
	}
	else
	{
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<3;m++)
			{
				sigma(l,m) *= 0.5 * ModuleBase::e2 * ModuleBase::FOUR_PI;
			}
		}
	}
	
	for(int l=0;l<3;l++)
	{
		if(is_pw) sigma(l,l) -= H_Hartree_pw::hartree_energy /GlobalC::ucell.omega;
		else sigma(l,l) += H_Hartree_pw::hartree_energy /GlobalC::ucell.omega;
		for(int m=0;m<l;m++)
		{
			sigma(m,l)=sigma(l,m);
		}
	}
	
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<3;m++)
		{
			sigma(l,m)*=-1;
		}
	}

	delete[] vh_g;
	delete[] psic;
	delete[] psic0;
	ModuleBase::timer::tick("Stress_Func","stress_har");
	return;
}
