#include "./stress_func.h"
#include "../module_xc/xc_functional.h"
#include "../module_base/math_integral.h"
#include "../module_base/timer.h"
#include "global.h"

//NLCC term, need to be tested
void Stress_Func::stress_cc(ModuleBase::matrix& sigma, ModulePW::PW_Basis* rho_basis, const bool is_pw, const Charge* const chr)
{
	ModuleBase::timer::tick("Stress_Func","stress_cc");
        
	double fact=1.0;

	if(is_pw&&INPUT.gamma_only) 
	{
		fact = 2.0; //is_pw:PW basis, gamma_only need to double.
	}

	std::complex<double> sigmadiag;
	double* rhocg;

	int judge=0;
	for(int nt=0;nt<GlobalC::ucell.ntype;nt++)
	{
		if(GlobalC::ucell.atoms[nt].ncpp.nlcc) 
		{
			judge++;
		}
	}

	if(judge==0) 
	{
		ModuleBase::timer::tick("Stress_Func","stress_cc");
		return;
	}

	//recalculate the exchange-correlation potential
	ModuleBase::matrix vxc;
	if(XC_Functional::get_func_type() == 3)
	{
#ifdef USE_LIBXC
    	const auto etxc_vtxc_v = XC_Functional::v_xc_meta(
            rho_basis->nrxx, rho_basis->nxyz, GlobalC::ucell.omega,
            chr);
        
        GlobalC::en.etxc = std::get<0>(etxc_vtxc_v);
        GlobalC::en.vtxc = std::get<1>(etxc_vtxc_v);
        vxc = std::get<2>(etxc_vtxc_v);
#else
        ModuleBase::WARNING_QUIT("cal_force_cc","to use mGGA, compile with LIBXC");
#endif
	}
	else
	{
    	const auto etxc_vtxc_v = XC_Functional::v_xc(rho_basis->nrxx, rho_basis->nxyz, GlobalC::ucell.omega, chr);
		GlobalC::en.etxc    = std::get<0>(etxc_vtxc_v);			// may delete?
		GlobalC::en.vtxc    = std::get<1>(etxc_vtxc_v);			// may delete?
		vxc = std::get<2>(etxc_vtxc_v);
	}

	std::complex<double> * psic = new std::complex<double> [rho_basis->nmaxgr];

	ModuleBase::GlobalFunc::ZEROS(psic, rho_basis->nrxx);

	if(GlobalV::NSPIN==1||GlobalV::NSPIN==4)
	{
		for(int ir=0;ir<rho_basis->nrxx;ir++)
		{
			// psic[ir] = vxc(0,ir);
			psic[ir] = std::complex<double>(vxc(0, ir),  0.0);
		}
	}
	else
	{
		for(int ir=0;ir<rho_basis->nrxx;ir++)
		{
			psic[ir] = 0.5 * (vxc(0, ir) + vxc(1, ir));
		}
	}

	// to G space
	rho_basis->real2recip(psic, psic); 

	//psic cantains now Vxc(G)
	rhocg= new double [rho_basis->ngg];
	ModuleBase::GlobalFunc::ZEROS(rhocg, rho_basis->ngg);

	sigmadiag=0.0;
	for(int nt=0;nt<GlobalC::ucell.ntype;nt++)
	{
		if(GlobalC::ucell.atoms[nt].ncpp.nlcc)
		{
			//drhoc();
			chr->non_linear_core_correction(
				GlobalC::ppcell.numeric,
				GlobalC::ucell.atoms[nt].ncpp.msh,
				GlobalC::ucell.atoms[nt].ncpp.r,
				GlobalC::ucell.atoms[nt].ncpp.rab,
				GlobalC::ucell.atoms[nt].ncpp.rho_atc,
				rhocg,
				rho_basis);


			//diagonal term 
			for(int ig = 0;ig< rho_basis->npw;ig++)
			{
				if(rho_basis->ig_gge0==ig)
					sigmadiag += conj(psic[ig] ) * GlobalC::sf.strucFac (nt, ig) * rhocg[rho_basis->ig2igg[ig]];
				else
					sigmadiag += conj(psic[ig] ) * GlobalC::sf.strucFac (nt, ig) * rhocg[rho_basis->ig2igg[ig]] * fact;
			}
			this->deriv_drhoc (
				GlobalC::ppcell.numeric,
				GlobalC::ucell.atoms[nt].ncpp.msh,
				GlobalC::ucell.atoms[nt].ncpp.r,
				GlobalC::ucell.atoms[nt].ncpp.rab,
				GlobalC::ucell.atoms[nt].ncpp.rho_atc,
				rhocg,
				rho_basis);
			// non diagonal term (g=0 contribution missing)
			for(int ig = 0;ig< rho_basis->npw;ig++)
			{
				const double norm_g = sqrt(rho_basis->gg[ig]);
				if(norm_g < 1e-4) 	continue;
				for (int l = 0; l < 3; l++)
				{
					for (int m = 0;m< 3;m++)
					{
						const std::complex<double> t = conj(psic[ig]) * GlobalC::sf.strucFac(nt, ig) * rhocg[rho_basis->ig2igg[ig]] * GlobalC::ucell.tpiba *
												  rho_basis->gcar[ig][l] * rho_basis->gcar[ig][m] / norm_g * fact;
						//						sigmacc [l][ m] += t.real();
						sigma(l,m) += t.real();
					}//end m
				}//end l
			}//end ng
		}//end if
	}//end nt

	for(int l = 0;l< 3;l++)
	{
		sigma(l,l) += sigmadiag.real();
//		sigmacc [l][ l] += sigmadiag.real();
	}
	for(int l = 0;l< 3;l++)
	{
		for (int m = 0;m< 3;m++)
		{
			Parallel_Reduce::reduce_double_pool( sigma(l,m) );
		}
	}

	delete[] rhocg;
	delete[] psic;

	ModuleBase::timer::tick("Stress_Func","stress_cc");
	return;
}


void Stress_Func::deriv_drhoc 
(
	const bool &numeric,
	const int mesh,
	const double *r,
	const double *rab,
	const double *rhoc,
	double *drhocg,
	ModulePW::PW_Basis* rho_basis
)
{

	double gx = 0, rhocg1 = 0;
	// the modulus of g for a given shell
	// the fourier transform
	double *aux = new double[ mesh];
	// auxiliary memory for integration

	int  igl0;
	// counter on radial mesh points
	// counter on g shells
	// lower limit for loop on ngl

	//
	// G=0 term
	//
	if (rho_basis->gg_uniq[0] < 1.0e-8)
	{
		drhocg [0] = 0.0;
		igl0 = 1;
	}
	else
	{
		igl0 = 0;
	}
	//
	// G <> 0 term
	//
	
	for(int igl = igl0;igl< rho_basis->ngg;igl++)
	{
		gx = sqrt(rho_basis->gg_uniq[igl] * GlobalC::ucell.tpiba2);
		for( int ir = 0;ir< mesh; ir++)
		{
			aux [ir] = r [ir] * rhoc [ir] * (r [ir] * cos (gx * r [ir] ) / gx - sin (gx * r [ir] ) / pow(gx,2));
		}//ir
		ModuleBase::Integral::Simpson_Integral(mesh, aux, rab, rhocg1);
		drhocg [igl] = ModuleBase::FOUR_PI / GlobalC::ucell.omega * rhocg1;
	}//igl
	
	delete [] aux;

	return;
}
