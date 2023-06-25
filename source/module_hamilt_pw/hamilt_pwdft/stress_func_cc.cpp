#include "stress_func.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_base/math_integral.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

//NLCC term, need to be tested
template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_cc(ModuleBase::matrix& sigma,
                                            ModulePW::PW_Basis* rho_basis,
                                            const Structure_Factor* p_sf,
                                            const bool is_pw,
                                            const Charge* const chr)
{
    ModuleBase::TITLE("Stress_Func","stress_cc");
	ModuleBase::timer::tick("Stress_Func","stress_cc");
        
	FPTYPE fact=1.0;

	if(is_pw&&INPUT.gamma_only) 
	{
		fact = 2.0; //is_pw:PW basis, gamma_only need to FPTYPE.
	}

	FPTYPE sigmadiag;
	FPTYPE* rhocg;

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
	if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
	{
#ifdef USE_LIBXC
        const auto etxc_vtxc_v
            = XC_Functional::v_xc_meta(rho_basis->nrxx, GlobalC::ucell.omega, GlobalC::ucell.tpiba, chr);

        // etxc = std::get<0>(etxc_vtxc_v);
        // vtxc = std::get<1>(etxc_vtxc_v);
        vxc = std::get<2>(etxc_vtxc_v);
#else
        ModuleBase::WARNING_QUIT("cal_force_cc","to use mGGA, compile with LIBXC");
#endif
	}
	else
	{
		if(GlobalV::NSPIN==4) GlobalC::ucell.cal_ux();
        const auto etxc_vtxc_v = XC_Functional::v_xc(rho_basis->nrxx, chr, &GlobalC::ucell);
        // etxc = std::get<0>(etxc_vtxc_v); // may delete?
        // vtxc = std::get<1>(etxc_vtxc_v); // may delete?
        vxc = std::get<2>(etxc_vtxc_v);
    }

    std::complex<FPTYPE>* psic = new std::complex<FPTYPE>[rho_basis->nmaxgr];

    if(GlobalV::NSPIN==1||GlobalV::NSPIN==4)
	{
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
		for(int ir=0;ir<rho_basis->nrxx;ir++)
		{
			// psic[ir] = vxc(0,ir);
			psic[ir] = std::complex<FPTYPE>(vxc(0, ir),  0.0);
		}
	}
	else
	{
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
		for(int ir=0;ir<rho_basis->nrxx;ir++)
		{
			psic[ir] = 0.5 * (vxc(0, ir) + vxc(1, ir));
		}
	}

	// to G space
	rho_basis->real2recip(psic, psic); 

	//psic cantains now Vxc(G)
	rhocg= new FPTYPE [rho_basis->ngg];

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
				rhocg);


			//diagonal term 
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sigmadiag) schedule(static, 256)
#endif
			for(int ig = 0;ig< rho_basis->npw;ig++)
			{
                std::complex<double> local_sigmadiag;
                if (rho_basis->ig_gge0 == ig)
                    local_sigmadiag = conj(psic[ig]) * p_sf->strucFac(nt, ig) * rhocg[rho_basis->ig2igg[ig]];
                else
                    local_sigmadiag = conj(psic[ig]) * p_sf->strucFac(nt, ig) * rhocg[rho_basis->ig2igg[ig]] * fact;
                sigmadiag += local_sigmadiag.real();
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
#ifdef _OPENMP
#pragma omp parallel
{
			ModuleBase::matrix local_sigma(3, 3);
			#pragma omp for
#else
			ModuleBase::matrix& local_sigma = sigma;
#endif
			for(int ig = 0;ig< rho_basis->npw;ig++)
			{
				const FPTYPE norm_g = sqrt(rho_basis->gg[ig]);
				if(norm_g < 1e-4) 	continue;
				for (int l = 0; l < 3; l++)
				{
					for (int m = 0;m< 3;m++)
					{
                        const std::complex<FPTYPE> t
                            = conj(psic[ig]) * p_sf->strucFac(nt, ig) * rhocg[rho_basis->ig2igg[ig]]
                              * GlobalC::ucell.tpiba * rho_basis->gcar[ig][l] * rho_basis->gcar[ig][m] / norm_g * fact;
                        //						sigmacc [l][ m] += t.real();
                        local_sigma(l,m) += t.real();
					}//end m
				}//end l
			}//end ng
#ifdef _OPENMP
			#pragma omp critical(stress_cc_reduce)
			{
				for(int l=0;l<3;l++)
				{
					for(int m=0;m<3;m++)
					{
						sigma(l,m) += local_sigma(l,m);
					}
				}
			}
}
#endif
		}//end if
	}//end nt

	for(int l = 0;l< 3;l++)
	{
		sigma(l,l) += sigmadiag;
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


template<typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::deriv_drhoc
(
	const bool &numeric,
	const int mesh,
	const FPTYPE *r,
	const FPTYPE *rab,
	const FPTYPE *rhoc,
	FPTYPE *drhocg,
	ModulePW::PW_Basis* rho_basis
)
{
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
#ifdef _OPENMP
#pragma omp parallel
{
#endif
	double gx = 0, rhocg1 = 0;
	// the modulus of g for a given shell
	// the fourier transform
	double *aux = new double[ mesh];
	// auxiliary memory for integration

	//
	// G <> 0 term
	//
#ifdef _OPENMP
#pragma omp for
#endif
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
#ifdef _OPENMP
}
#endif
	return;
}

template class Stress_Func<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, psi::DEVICE_GPU>;
#endif