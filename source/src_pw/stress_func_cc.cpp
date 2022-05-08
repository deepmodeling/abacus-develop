#include "./stress_func.h"
#include "../module_xc/xc_functional.h"
#include "../module_base/math_integral.h"
#include "../module_base/timer.h"
#include "global.h"

//NLCC term, need to be tested
void Stress_Func::stress_cc(ModuleBase::matrix& sigma, const bool is_pw)
{
	ModuleBase::timer::tick("Stress_Func","stress_cc");
        
	int nt,ng,l,m,ir;
	double fact=1.0;

	if(is_pw&&INPUT.gamma_only) 
	{
		fact = 2.0; //is_pw:PW basis, gamma_only need to double.
	}

	std::complex<double> sigmadiag;
	double* rhocg;

	int judge=0;
	for(nt=0;nt<GlobalC::ucell.ntype;nt++)
	{
		if(GlobalC::ucell.atoms[nt].nlcc) 
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
    const auto etxc_vtxc_v = XC_Functional::v_xc(GlobalC::pw.nrxx, GlobalC::pw.ncxyz, GlobalC::ucell.omega, GlobalC::CHR.rho, GlobalC::CHR.rho_core);
	GlobalC::en.etxc    = std::get<0>(etxc_vtxc_v);			// may delete?
	GlobalC::en.vtxc    = std::get<1>(etxc_vtxc_v);			// may delete?
	const ModuleBase::matrix vxc = std::get<2>(etxc_vtxc_v);

	std::complex<double> * psic = new std::complex<double> [GlobalC::pw.nrxx];

	ModuleBase::GlobalFunc::ZEROS(psic, GlobalC::pw.nrxx);

	if(GlobalV::NSPIN==1||GlobalV::NSPIN==4)
	{
		for(ir=0;ir<GlobalC::pw.nrxx;ir++)
		{
			// psic[ir] = vxc(0,ir);
			psic[ir] = std::complex<double>(vxc(0, ir),  0.0);
		}
	}
	else
	{
		for(ir=0;ir<GlobalC::pw.nrxx;ir++)
		{
			psic[ir] = 0.5 * (vxc(0, ir) + vxc(1, ir));
		}
	}

	// to G space
	GlobalC::pw.FFT_chg.FFT3D(psic, -1);

	//psic cantains now Vxc(G)
	rhocg= new double [GlobalC::pw.nggm];
	ModuleBase::GlobalFunc::ZEROS(rhocg, GlobalC::pw.nggm);

	sigmadiag=0.0;
	for(nt=0;nt<GlobalC::ucell.ntype;nt++)
	{
		if(GlobalC::ucell.atoms[nt].nlcc)
		{
			//drhoc();
			GlobalC::CHR.non_linear_core_correction(
				GlobalC::ppcell.numeric,
				GlobalC::ucell.atoms[nt].msh,
				GlobalC::ucell.atoms[nt].r,
				GlobalC::ucell.atoms[nt].rab,
				GlobalC::ucell.atoms[nt].rho_atc,
				rhocg);


			//diagonal term 
			if (GlobalC::pw.gstart==1) 
			{
				sigmadiag += conj(psic [GlobalC::pw.ig2fftc[0]] ) * GlobalC::pw.strucFac (nt, 0) * rhocg [GlobalC::pw.ig2ngg[0] ];
			}
			for( ng = GlobalC::pw.gstart;ng< GlobalC::pw.ngmc;ng++)
			{
				sigmadiag +=  conj(psic[GlobalC::pw.ig2fftc[ng]] ) *
					GlobalC::pw.strucFac (nt, ng) * rhocg [GlobalC::pw.ig2ngg[ng] ] * fact;
			}
			this->deriv_drhoc (
				GlobalC::ppcell.numeric,
				GlobalC::ucell.atoms[nt].msh,
				GlobalC::ucell.atoms[nt].r,
				GlobalC::ucell.atoms[nt].rab,
				GlobalC::ucell.atoms[nt].rho_atc,
				rhocg);
			// non diagonal term (g=0 contribution missing)
			for( ng = GlobalC::pw.gstart;ng< GlobalC::pw.ngmc;ng++)
			{
				double norm_g = sqrt(GlobalC::pw.get_NormG_cartesian(ng));
				assert(norm_g != 0.000);
				for (l = 0; l < 3; l++)
				{
					for (m = 0;m< 3;m++)
					{
						const std::complex<double> t = conj(psic[GlobalC::pw.ig2fftc[ng]]) * GlobalC::pw.strucFac(nt, ng) * rhocg[GlobalC::pw.ig2ngg[ng]] * GlobalC::ucell.tpiba *
												  GlobalC::pw.get_G_cartesian_projection(ng, l) * GlobalC::pw.get_G_cartesian_projection(ng, m) 
												  / norm_g * fact;
						//						sigmacc [l][ m] += t.real();
						sigma(l,m) += t.real();
					}//end m
				}//end l
			}//end ng
		}//end if
	}//end nt

	for( l = 0;l< 3;l++)
	{
		sigma(l,l) += sigmadiag.real();
//		sigmacc [l][ l] += sigmadiag.real();
	}
	for( l = 0;l< 3;l++)
	{
		for (m = 0;m< 3;m++)
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
	double *drhocg
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
	if (GlobalC::pw.ggs[0] < 1.0e-8)
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
	
	for(int igl = igl0;igl< GlobalC::pw.nggm;igl++)
	{
		gx = sqrt(GlobalC::pw.ggs [igl] * GlobalC::ucell.tpiba2);
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
