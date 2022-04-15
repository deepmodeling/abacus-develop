#include "stress_func.h"
#include "../module_base/math_integral.h"
#include "../module_base/timer.h"
#include "global.h"

//calculate local pseudopotential stress in PW or VL_dVL stress in LCAO
void Stress_Func::stress_loc(ModuleBase::matrix& sigma, const bool is_pw)
{
    ModuleBase::timer::tick("Stress_Func","stress_loc");

    double *dvloc;
    double evloc=0.0;
	double fact=1.0;

	if (INPUT.gamma_only && is_pw) fact=2.0;

    dvloc = new double[GlobalC::pw.ngmc];

	std::complex<double> *Porter = GlobalC::UFFT.porter;

	ModuleBase::GlobalFunc::ZEROS( Porter, GlobalC::pw.nrxx );
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for (int ir=0; ir<GlobalC::pw.nrxx; ir++)
		{
			Porter[ir] += std::complex<double>(GlobalC::CHR.rho[is][ir], 0.0 );
		}
	}
	GlobalC::pw.FFT_chg.FFT3D(Porter, -1);

//    if(INPUT.gamma_only==1) fact=2.0;
//    else fact=1.0;

	evloc=0.0;

	std::complex<double> *vg = new std::complex<double>[GlobalC::pw.ngmc];
	ModuleBase::GlobalFunc::ZEROS( vg, GlobalC::pw.ngmc );
	for (int it=0; it<GlobalC::ucell.ntype; it++)
	{
		if (GlobalC::pw.gstart==1) evloc += GlobalC::ppcell.vloc(it, GlobalC::pw.ig2ngg[0]) * (GlobalC::pw.strucFac(it,0) * conj(Porter[GlobalC::pw.ig2fftc[0]])).real();
		for (int ig=GlobalC::pw.gstart; ig<GlobalC::pw.ngmc; ig++)
		{
			const int j = GlobalC::pw.ig2fftc[ig];
			evloc += GlobalC::ppcell.vloc(it, GlobalC::pw.ig2ngg[ig]) * (GlobalC::pw.strucFac(it,ig) * conj(Porter[j]) * fact).real();
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
			this->dvloc_coul (atom->zv, dvloc);
		//
		}
		else
		{
		//
		// normal case: dvloc contains dV_loc(G)/dG
		//
			this->dvloc_of_g ( atom->msh, atom->rab, atom->r,
					atom->vloc_at, atom->zv, dvloc);
		//
		}

		for(int ng = 0;ng< GlobalC::pw.ngmc;ng++)
		{
			const int j = GlobalC::pw.ig2fftc[ng];
			for (int l = 0;l< 3;l++)
			{
				for (int m = 0; m<l+1;m++)
				{
					sigma(l, m) = sigma(l, m) + (conj(Porter[j]) * GlobalC::pw.strucFac(nt, ng)).real() 
						* 2.0 * dvloc[GlobalC::pw.ig2ngg[ng]] * GlobalC::ucell.tpiba2 * 
						GlobalC::pw.get_G_cartesian_projection(ng, l) * GlobalC::pw.get_G_cartesian_projection(ng, m) * fact;
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
	delete[] vg;


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
double*  dvloc
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
	if (GlobalC::pw.ggs[0] < 1.0e-8)
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
	for(int igl = igl0;igl< GlobalC::pw.nggm;igl++)
	{
		double gx = sqrt (GlobalC::pw.ggs [igl] * GlobalC::ucell.tpiba2);
		double gx2 = GlobalC::pw.ggs [igl] * GlobalC::ucell.tpiba2;
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
double* dvloc
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
	if (GlobalC::pw.ggs[0] < 1.0e-8)
	{
		dvloc[0] = 0.0;
		igl0 = 1;
	}
	else
	{
		igl0 = 0;
	}
	for(int i=igl0;i<GlobalC::pw.nggm;i++)
	{
		dvloc[i] = ModuleBase::FOUR_PI * zp * ModuleBase::e2 / GlobalC::ucell.omega / pow(( GlobalC::ucell.tpiba2 * GlobalC::pw.ggs[i] ),2);
	}
	
	return;
}
