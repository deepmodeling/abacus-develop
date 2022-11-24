#include "stress_func.h"
#include "../module_base/math_integral.h"
#include "../module_base/timer.h"
#include "global.h"

//calculate local pseudopotential stress in PW or VL_dVL stress in LCAO
void Stress_Func::stress_loc(ModuleBase::matrix& sigma, ModulePW::PW_Basis* rho_basis, const bool is_pw, const Charge* const chr)
{
    ModuleBase::timer::tick("Stress_Func","stress_loc");

    double *dvloc = new double[rho_basis->npw];
    double evloc=0.0;
	double fact=1.0;

	if (INPUT.gamma_only && is_pw) fact=2.0;

    

	std::complex<double> *aux = new std::complex<double> [rho_basis->nmaxgr];

	ModuleBase::GlobalFunc::ZEROS( aux, rho_basis->nrxx );
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for (int ir=0; ir<rho_basis->nrxx; ir++)
		{
			aux[ir] += std::complex<double>(chr->rho[is][ir], 0.0 );
		}
	}
	rho_basis->real2recip(aux,aux);

//    if(INPUT.gamma_only==1) fact=2.0;
//    else fact=1.0;

	evloc=0.0;

	for (int it=0; it<GlobalC::ucell.ntype; it++)
	{
		for (int ig=0; ig<rho_basis->npw; ig++)
		{
			if(rho_basis->ig_gge0==ig)
				evloc += GlobalC::ppcell.vloc(it, rho_basis->ig2igg[ig]) * (GlobalC::sf.strucFac(it,ig) * conj(aux[ig])).real();
			else
				evloc += GlobalC::ppcell.vloc(it, rho_basis->ig2igg[ig]) * (GlobalC::sf.strucFac(it,ig) * conj(aux[ig]) * fact).real();
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

		for(int ig = 0;ig< rho_basis->npw;ig++)
		{
			for (int l = 0;l< 3;l++)
			{
				for (int m = 0; m<l+1;m++)
				{
					sigma(l, m) = sigma(l, m) + (conj(aux[ig]) * GlobalC::sf.strucFac(nt, ig)).real() 
						* 2.0 * dvloc[rho_basis->ig2igg[ig]] * GlobalC::ucell.tpiba2 * 
						rho_basis->gcar[ig][l] * rho_basis->gcar[ig][m] * fact;
				}
			}
		}
	}

    for(int l = 0;l< 3;l++)
	{
		if(is_pw) sigma(l,l) += evloc;
		for (int m = 0; m<l+1; m++)
		{
			sigma(m,l) = sigma(l,m);
		}
	}
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<3;m++)
		{
			Parallel_Reduce::reduce_double_pool( sigma(l,m) );
		}
	}
	delete[] dvloc;
	delete[] aux;


	ModuleBase::timer::tick("Stress_Func","stress_loc");
	return;
}

void Stress_Func::dvloc_of_g 
(
const int& msh,
const double* rab,
const double* r,
const double* vloc_at,
const double& zp,
double*  dvloc,
ModulePW::PW_Basis* rho_basis
)
{
  //----------------------------------------------------------------------
  //
  // dvloc = D Vloc (g^2) / D g^2 = (1/2g) * D Vloc(g) / D g
  //

  //
  //double  dvloc[ngl];
  // the fourier transform dVloc/dG
  //
	double vlcp=0;
	double  *aux, *aux1;

	int igl0;
	// counter on erf functions or gaussians
	// counter on g shells vectors
	// first shell with g != 0

	aux = new double[msh];
	aux1 = new double[msh];
	ModuleBase::GlobalFunc::ZEROS(aux, msh);
	ModuleBase::GlobalFunc::ZEROS(aux1, msh);

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

	//
	//   This is the part of the integrand function
	//   indipendent of |G| in real space
	//
	for(int i = 0;i< msh; i++)
	{
		aux1[i] = r [i] * vloc_at [i] + zp * ModuleBase::e2 * erf(r[i]);
	}
	for(int igl = igl0;igl< rho_basis->ngg;igl++)
	{
		const double g2 = rho_basis->gg_uniq[igl];
		const double gx = sqrt (g2 * GlobalC::ucell.tpiba2);
		const double gx2 = g2 * GlobalC::ucell.tpiba2;
		//
		//    and here we perform the integral, after multiplying for the |G|
		//    dependent  part
		//
		// DV(g)/Dg = ModuleBase::Integral of r (Dj_0(gr)/Dg) V(r) dr
		for(int i = 1;i< msh;i++)
		{
			aux [i] = aux1 [i] * (r [i] * cos (gx * r [i] ) / gx - sin (gx * r [i] ) / pow(gx,2));
		}
		// simpson (msh, aux, rab, vlcp);
		ModuleBase::Integral::Simpson_Integral(msh, aux, rab, vlcp );
		// DV(g^2)/Dg^2 = (DV(g)/Dg)/2g
		vlcp *= ModuleBase::FOUR_PI / GlobalC::ucell.omega / 2.0 / gx;
		// subtract the long-range term
		double g2a = gx2 / 4.0;
		vlcp += ModuleBase::FOUR_PI / GlobalC::ucell.omega * zp * ModuleBase::e2 * exp ( - g2a) * (g2a + 1) / pow(gx2 , 2);
		dvloc [igl] = vlcp;
	}
	delete[] aux1;
	delete[] aux;

	return;
}

void Stress_Func::dvloc_coul 
(
const double& zp,
double* dvloc,
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
	for(int i=igl0;i<rho_basis->ngg;i++)
	{
		dvloc[i] = ModuleBase::FOUR_PI * zp * ModuleBase::e2 / GlobalC::ucell.omega / pow(( GlobalC::ucell.tpiba2 * rho_basis->gg_uniq[i] ),2);
	}
	
	return;
}
