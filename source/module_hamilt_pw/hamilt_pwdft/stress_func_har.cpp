#include "stress_func.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

//calculate the Hartree part in PW or LCAO base
template<typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_har(ModuleBase::matrix& sigma, ModulePW::PW_Basis* rho_basis, const bool is_pw, const Charge* const chr)
{
    ModuleBase::TITLE("Stress_Func","stress_har");
	ModuleBase::timer::tick("Stress_Func","stress_har");

	std::complex<FPTYPE> *aux = new std::complex<FPTYPE>[rho_basis->nmaxgr];

	//  Hartree potential VH(r) from n(r)
    /*
        blocking rho_basis->nrxx for data locality.

        By blocking aux with block size 1024,
        we can keep the blocked aux in L1 cache when iterating GlobalV::NSPIN loop
        performance will be better when number of atom is quite huge
    */
    const int block_ir = 1024;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int irb = 0; irb < rho_basis->nrxx; irb += block_ir)
	{
		// calculate the actual task length of this block
 		int ir_end = std::min(irb + block_ir, rho_basis->nrxx);

		{ // is = 0
			for (int ir = irb; ir < ir_end; ++ir)
			{ // initialize aux
				aux[ir] = std::complex<FPTYPE>( chr->rho[0][ir], 0.0 );
			}
		}
		for (int is = 1; is < GlobalV::NSPIN; is++)
		{
			for (int ir = irb; ir < ir_end; ++ir)
			{ // accumulate aux
				aux[ir] += std::complex<FPTYPE>( chr->rho[is][ir], 0.0 );
			}
		}
	}
	//=============================
	//  bring rho (aux) to G space
	//=============================
	rho_basis->real2recip(aux, aux);

//	FPTYPE ehart=0;
#ifndef _OPENMP
	ModuleBase::matrix& local_sigma = sigma;
#else
#pragma omp parallel
{
	ModuleBase::matrix local_sigma(3, 3);
	#pragma omp for
#endif
	for (int ig = 0 ; ig < rho_basis->npw ; ++ig)
	{
		const FPTYPE g2 = rho_basis->gg[ig];
		if(g2 < 1e-8) continue;
		//const FPTYPE fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (GlobalC::ucell.tpiba2 * GlobalC::sf.gg [ig]);
		//ehart += ( conj( Porter[j] ) * Porter[j] ).real() * fac;
		//vh_g[ig] = fac * Porter[j];
		FPTYPE shart= ( conj( aux[ig] ) * aux[ig] ).real()/(GlobalC::ucell.tpiba2 * g2);
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<l+1;m++)
			{
				local_sigma(l, m) += shart * 2 * rho_basis->gcar[ig][l] * rho_basis->gcar[ig][m] / g2;
			}
		}
	}
#ifdef _OPENMP
	#pragma omp critical(stress_har_reduce)
	{
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<l+1;m++)
			{
				sigma(l,m) += local_sigma(l,m);
			}
		}
	}
}
#endif

    for(int l=0;l<3;l++)
	{
		for(int m=0;m<l+1;m++)
		{
            Parallel_Reduce::reduce_pool(sigma(l, m));
		}
	}

    //        Parallel_Reduce::reduce_pool( ehart );
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

template class Stress_Func<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, psi::DEVICE_GPU>;
#endif