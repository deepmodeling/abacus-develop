#include "stress_func.h"
#include "module_base/math_integral.h"
#include "module_base/tool_threading.h"
#include "module_base/timer.h"
#include "module_base/libm/libm.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

//calculate local pseudopotential stress in PW or VL_dVL stress in LCAO
template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_loc(ModuleBase::matrix& sigma,
                                             ModulePW::PW_Basis* rho_basis,
                                             const Structure_Factor* p_sf,
                                             const bool is_pw,
                                             const Charge* const chr)
{
    ModuleBase::TITLE("Stress_Func","stress_loc");
    ModuleBase::timer::tick("Stress_Func","stress_loc");

    FPTYPE *dvloc = new FPTYPE[rho_basis->npw];
    FPTYPE evloc=0.0;
	FPTYPE fact=1.0;

	if (INPUT.gamma_only && is_pw) fact=2.0;

    

	std::complex<FPTYPE> *aux = new std::complex<FPTYPE> [rho_basis->nmaxgr];

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
				aux[ir] = std::complex<FPTYPE>(chr->rho[0][ir], 0.0 );
			}
		}
		for (int is = 1; is < GlobalV::NSPIN; is++)
		{
			for (int ir = irb; ir < ir_end; ++ir)
			{ // accumulate aux
				aux[ir] += std::complex<FPTYPE>(chr->rho[is][ir], 0.0 );
			}
		}
 	}
	rho_basis->real2recip(aux,aux);

//    if(INPUT.gamma_only==1) fact=2.0;
//    else fact=1.0;

	if(is_pw)
	{
#pragma omp parallel for collapse(2) reduction(+:evloc)
		for (int it=0; it<GlobalC::ucell.ntype; it++)
		{
			for (int ig=0; ig<rho_basis->npw; ig++)
			{
                if (rho_basis->ig_gge0 == ig)
                    evloc += GlobalC::ppcell.vloc(it, rho_basis->ig2igg[ig])
                             * (p_sf->strucFac(it, ig) * conj(aux[ig])).real();
                else
                    evloc += GlobalC::ppcell.vloc(it, rho_basis->ig2igg[ig])
                             * (p_sf->strucFac(it, ig) * conj(aux[ig]) * fact).real();
            }
		}
	}
	for(int nt = 0;nt< GlobalC::ucell.ntype; nt++)
	{
		const Atom* atom = &GlobalC::ucell.atoms[nt];
		//mark by zhengdy for check
		// if ( GlobalC::ppcell.vloc == NULL ){
		if(0)
		{
		//
		// special case: pseudopotential is coulomb 1/r potential
		//
			this->dvloc_coul (atom->ncpp.zv, dvloc, rho_basis);
		//
		}
		else
		{
		//
		// normal case: dvloc contains dV_loc(G)/dG
		//
			this->dvloc_of_g ( atom->ncpp.msh, atom->ncpp.rab, atom->ncpp.r,
					atom->ncpp.vloc_at, atom->ncpp.zv, dvloc, rho_basis);
		//
		}
#ifndef _OPENMP
		ModuleBase::matrix &local_sigma = sigma;
#else
#pragma omp parallel
{
		ModuleBase::matrix local_sigma(3, 3);
		#pragma omp for
#endif
		for(int ig = 0;ig< rho_basis->npw;ig++)
		{
			for (int l = 0;l< 3;l++)
			{
				for (int m = 0; m<l+1;m++)
				{
                    local_sigma(l, m) = local_sigma(l, m)
                                        + (conj(aux[ig]) * p_sf->strucFac(nt, ig)).real() * 2.0
                                              * dvloc[rho_basis->ig2igg[ig]] * GlobalC::ucell.tpiba2
                                              * rho_basis->gcar[ig][l] * rho_basis->gcar[ig][m] * fact;
                }
			}
		}
#ifdef _OPENMP
		#pragma omp critical(stress_loc_reduce)
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
	}

	if(is_pw)
	{
		for(int l = 0;l< 3;l++)
		{
			sigma(l,l) += evloc;
		}
	}
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<l+1;m++)
		{
            Parallel_Reduce::reduce_pool(sigma(l, m));
		}
	}
	for(int l = 0;l< 3;l++)
	{
		for (int m = 0; m<l+1; m++)
		{
			sigma(m,l) = sigma(l,m);
		}
	}
	delete[] dvloc;
	delete[] aux;


	ModuleBase::timer::tick("Stress_Func","stress_loc");
	return;
}

template<typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::dvloc_of_g
(
const int& msh,
const FPTYPE* rab,
const FPTYPE* r,
const FPTYPE* vloc_at,
const FPTYPE& zp,
FPTYPE*  dvloc,
ModulePW::PW_Basis* rho_basis
)
{
  //----------------------------------------------------------------------
  //
  // dvloc = D Vloc (g^2) / D g^2 = (1/2g) * D Vloc(g) / D g
  //

  //
  //FPTYPE  dvloc[ngl];
  // the fourier transform dVloc/dG
  //
	FPTYPE  *aux1;

	int igl0;
	// counter on erf functions or gaussians
	// counter on g shells vectors
	// first shell with g != 0

	aux1 = new FPTYPE[msh];

	// the  G=0 component is not computed
	if (rho_basis->gg_uniq[0] < 1.0e-8)
	{
		dvloc[0] = 0.0;
		igl0 = 1;
	}
	else
	{
		igl0 = 0;
	}

	// Pseudopotentials in numerical form (Vloc contains the local part)
	// In order to perform the Fourier transform, a term erf(r)/r is
	// subtracted in real space and added again in G space

#ifdef _OPENMP
#pragma omp parallel
{
	#pragma omp for
#endif
	//
	//   This is the part of the integrand function
	//   indipendent of |G| in real space
	//
	for(int i = 0;i< msh; i++)
	{
		aux1[i] = r [i] * vloc_at [i] + zp * ModuleBase::e2 * erf(r[i]);
	}

	FPTYPE  *aux;
	aux = new FPTYPE[msh];
	aux[0] = 0.0;
#ifdef _OPENMP
	#pragma omp for
#endif
	for(int igl = igl0;igl< rho_basis->ngg;igl++)
	{
		const FPTYPE g2 = rho_basis->gg_uniq[igl];
		const FPTYPE gx = sqrt (g2 * GlobalC::ucell.tpiba2);
		const FPTYPE gx2 = g2 * GlobalC::ucell.tpiba2;
		//
		//    and here we perform the integral, after multiplying for the |G|
		//    dependent  part
		//
		// DV(g)/Dg = ModuleBase::Integral of r (Dj_0(gr)/Dg) V(r) dr
		for(int i = 1;i< msh;i++)
		{
			FPTYPE sinp, cosp;
            ModuleBase::libm::sincos(gx * r [i], &sinp, &cosp);
			aux [i] = aux1 [i] * (r [i] * cosp / gx - sinp / pow(gx,2));
		}
		FPTYPE vlcp=0;
		// simpson (msh, aux, rab, vlcp);
		ModuleBase::Integral::Simpson_Integral(msh, aux, rab, vlcp );
		// DV(g^2)/Dg^2 = (DV(g)/Dg)/2g
		vlcp *= ModuleBase::FOUR_PI / GlobalC::ucell.omega / 2.0 / gx;
		// subtract the long-range term
		FPTYPE g2a = gx2 / 4.0;
		vlcp += ModuleBase::FOUR_PI / GlobalC::ucell.omega * zp * ModuleBase::e2 * ModuleBase::libm::exp ( - g2a) * (g2a + 1) / pow(gx2 , 2);
		dvloc [igl] = vlcp;
	}
	delete[] aux;
#ifdef _OPENMP
}
#endif
	delete[] aux1;

	return;
}

template<typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::dvloc_coul
(
const FPTYPE& zp,
FPTYPE* dvloc,
ModulePW::PW_Basis* rho_basis
)
{
	//----------------------------------------------------------------------
	//
	//    Fourier transform of the Coulomb potential - For all-electron
	//    calculations, in specific cases only, for testing purposes
	//


	// fourier transform: dvloc = D Vloc (g^2) / D g^2 = 4pi e^2/omegai /G^4

	int  igl0;
	// first shell with g != 0

	// the  G=0 component is 0
	if (rho_basis->gg_uniq[0] < 1.0e-8)
	{
		dvloc[0] = 0.0;
		igl0 = 1;
	}
	else
	{
		igl0 = 0;
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(int i=igl0;i<rho_basis->ngg;i++)
	{
		dvloc[i] = ModuleBase::FOUR_PI * zp * ModuleBase::e2 / GlobalC::ucell.omega / pow(( GlobalC::ucell.tpiba2 * rho_basis->gg_uniq[i] ),2);
	}

	return;
}

template class Stress_Func<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, psi::DEVICE_GPU>;
#endif