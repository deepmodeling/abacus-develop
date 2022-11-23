#include "./stress_func.h"
#include "./myfunc.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "../module_base/timer.h"
#include "global.h"

//calculate the Hartree part in PW or LCAO base
void Stress_Func::stress_har(ModuleBase::matrix& sigma, ModulePW::PW_Basis* rho_basis, const bool is_pw, const Charge* const chr)
{
	ModuleBase::timer::tick("Stress_Func","stress_har");
	double shart;

	std::complex<double> *aux = new std::complex<double>[rho_basis->nmaxgr];

	//  Hartree potential VH(r) from n(r)
	ModuleBase::GlobalFunc::ZEROS( aux, rho_basis->nrxx );
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for (int ir=0; ir<rho_basis->nrxx; ir++)
		{
			aux[ir] += std::complex<double>( chr->rho[is][ir], 0.0 );
		}
	}
	//=============================
	//  bring rho (aux) to G space
	//=============================
	rho_basis->real2recip(aux, aux);


//	double ehart=0;
	for (int ig = 0 ; ig < rho_basis->npw ; ++ig)
	{
		const double g2 = rho_basis->gg[ig];
		if(g2 < 1e-8) continue;
		//const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (GlobalC::ucell.tpiba2 * GlobalC::sf.gg [ig]);
		//ehart += ( conj( Porter[j] ) * Porter[j] ).real() * fac;
		//vh_g[ig] = fac * Porter[j];
		shart= ( conj( aux[ig] ) * aux[ig] ).real()/(GlobalC::ucell.tpiba2 * g2);
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<l+1;m++)
			{
				sigma(l, m) += shart * 2 * rho_basis->gcar[ig][l] * rho_basis->gcar[ig][m] / g2;
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
		if(is_pw) sigma(l,l) -= elecstate::H_Hartree_pw::hartree_energy /GlobalC::ucell.omega;
		else sigma(l,l) += elecstate::H_Hartree_pw::hartree_energy /GlobalC::ucell.omega;
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

	delete[] aux;
	ModuleBase::timer::tick("Stress_Func","stress_har");
	return;
}
