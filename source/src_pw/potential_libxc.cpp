//==========================================================
// AUTHOR : Peize Lin
// DATE :   2017-09-14
// UPDATE : 2021-02-28
//==========================================================

#ifdef USE_LIBXC

#include "potential_libxc.h"
#include "global.h"
#include "tools.h"
#include "xc_gga_pw.h"
#include "../module_base/global_function.h"
#ifdef __LCAO
#include "../src_lcao/global_fp.h"
#endif

//the interfacd to libxc xc_lda_exc_vxc and xc_gga_exc_vxc
//XC_POLARIZED, XC_UNPOLARIZED: internal flags used in LIBXC, denote the polarized(nspin=1) or unpolarized(nspin=2) calculations, definition can be found in xc.h from LIBXC
//XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_CORRELATION: internal flags used in LIBXC, denote the types of functional associated with a certain functional ID, definition can be found in xc.h from LIBXC

// [etxc, vtxc, v] = Potential_Libxc::v_xc(...)
std::tuple<double,double,ModuleBase::matrix> Potential_Libxc::v_xc(
	const double * const * const rho_in,
	const double * const rho_core_in)
{
    ModuleBase::TITLE("Potential_Libxc","v_xc");
    ModuleBase::timer::tick("Potential_Libxc","v_xc");

    double etxc = 0.0;
    double vtxc = 0.0;
	ModuleBase::matrix v(GlobalV::NSPIN,GlobalC::pw.nrxx);

	if(GlobalV::VXC_IN_H == 0 )
	{
    	ModuleBase::timer::tick("Potential_Libxc","v_xc");
		return std::make_tuple( etxc, vtxc, std::move(v) );
	}
	//----------------------------------------------------------
	// xc_func_type is defined in Libxc package
	// to understand the usage of xc_func_type,
	// use can check on website, for example:
	// https://www.tddft.org/programs/libxc/manual/libxc-5.1.x/
	//----------------------------------------------------------
	std::vector<xc_func_type> funcs = init_func();

	// the type of rho_sigma_gdr is automatically set to 'tuple'
	// [rho, sigma, gdr] = cal_input( funcs, rho_in );
	const auto rho_sigma_gdr = cal_input( funcs, rho_in, rho_core_in );
	const std::vector<double> &rho = std::get<0>(rho_sigma_gdr);
	const std::vector<double> &sigma = std::get<1>(rho_sigma_gdr);

	for( xc_func_type &func : funcs )
	{
		// jiyy add for threshold
		constexpr double rho_threshold = 1E-6;
		constexpr double grho_threshold = 1E-10;
		xc_func_set_dens_threshold(&func, rho_threshold);
		// sgn for threshold mask
		const std::vector<double> sgn = [&]() -> std::vector<double>
		{
			std::vector<double> sgn( GlobalC::pw.nrxx * nspin0(), 1.0);
			if(nspin0()==2 && func.info->family != XC_FAMILY_LDA && func.info->kind==XC_CORRELATION)
			{
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
				{
					if ( rho[ir*2]<rho_threshold || sqrt(abs(sigma[ir*3]))<grho_threshold ) 
						sgn[ir*2] = 0.0;
					if ( rho[ir*2+1]<rho_threshold || sqrt(abs(sigma[ir*3+2]))<grho_threshold ) 
						sgn[ir*2+1] = 0.0;
				}
			}
			return sgn;
		}();

		std::vector<double> exc   ( GlobalC::pw.nrxx                       );
		std::vector<double> vrho  ( GlobalC::pw.nrxx * nspin0()            );
		std::vector<double> vsigma( GlobalC::pw.nrxx * ((1==nspin0())?1:3) );
		
		// cal etxc from rho, exc
		auto process_exc = [&]()
		{
			for( size_t is=0; is!=nspin0(); ++is )
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
					etxc += ModuleBase::e2 * exc[ir] * rho[ir*nspin0()+is] * sgn[ir*nspin0()+is];
		};

		// cal vtx, v from rho_in, vrho
		auto process_vrho = [&]()
		{
			if(nspin0()==1 || GlobalV::NSPIN==2)
			{
				for( size_t is=0; is!=nspin0(); ++is )
				{
					for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
					{
						const double v_tmp = ModuleBase::e2 * vrho[ir*nspin0()+is] * sgn[ir*nspin0()+is];
						v(is,ir) += v_tmp;
						vtxc += v_tmp * rho_in[is][ir];
					}
				}
			}
			else // may need updates for SOC
			{
				constexpr double vanishing_charge = 1.0e-12;
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
				{
					std::vector<double> v_tmp(4);
					v_tmp[0] = ModuleBase::e2 * (0.5 * (vrho[ir*2] + vrho[ir*2+1]));
					const double vs = 0.5 * (vrho[ir*2] - vrho[ir*2+1]);
					const double amag = sqrt( pow(rho_in[1][ir],2) 
						+ pow(rho_in[2][ir],2) 
						+ pow(rho_in[3][ir],2) );

					if(amag>vanishing_charge)
					{
						for(int ipol=1; ipol<4; ++ipol)
						{
							v_tmp[ipol] = ModuleBase::e2 * vs * rho_in[ipol][ir] / amag;
						}
					}
					for(int ipol=0; ipol<4; ++ipol)
					{
						v(ipol, ir) += v_tmp[ipol];
						vtxc += v_tmp[ipol] * rho_in[ipol][ir];
					}
				}
			}
			
		};

		// cal vtxc, v from rho_in, rho, gdr, vsigma
		auto process_vsigma = [&]()
		{
			const std::vector<std::vector<ModuleBase::Vector3<double>>> &gdr = std::get<2>(rho_sigma_gdr);
			
			std::vector<std::vector<ModuleBase::Vector3<double>>> h( nspin0(), std::vector<ModuleBase::Vector3<double>>(GlobalC::pw.nrxx) );
			if( 1==nspin0() )
			{
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
				{
					h[0][ir] = ModuleBase::e2 * gdr[0][ir] * vsigma[ir] * 2.0 * sgn[ir];
				}
			}
			else
			{
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
				{
					h[0][ir] = ModuleBase::e2 * (gdr[0][ir] * vsigma[ir*3  ] * 2.0 * sgn[ir*2  ]
						           + gdr[1][ir] * vsigma[ir*3+1]       * sgn[ir*2]   * sgn[ir*2+1]);
					h[1][ir] = ModuleBase::e2 * (gdr[1][ir] * vsigma[ir*3+2] * 2.0 * sgn[ir*2+1]
					               + gdr[0][ir] * vsigma[ir*3+1]       * sgn[ir*2]   * sgn[ir*2+1]);
				}
			}

			// define two dimensional array dh [ nspin, GlobalC::pw.nrxx ]
			std::vector<std::vector<double>> dh(nspin0(), std::vector<double>(GlobalC::pw.nrxx));
			for( size_t is=0; is!=nspin0(); ++is )
				GGA_PW::grad_dot( ModuleBase::GlobalFunc::VECTOR_TO_PTR(h[is]), ModuleBase::GlobalFunc::VECTOR_TO_PTR(dh[is]) );

			for( size_t is=0; is!=nspin0(); ++is )
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
					vtxc -= dh[is][ir] * rho[ir*nspin0()+is];

			if(nspin0()==1 || GlobalV::NSPIN==2)
			{
				for( size_t is=0; is!=nspin0(); ++is )
					for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
						v(is,ir) -= dh[is][ir];
			}
			else // may need updates for SOC
			{
				constexpr double vanishing_charge = 1.0e-12;
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
				{
					v(0,ir) -= 0.5 * (dh[0][ir] + dh[1][ir]);
					const double amag = sqrt( pow(rho_in[1][ir],2) + pow(rho_in[2][ir],2) + pow(rho_in[3][ir],2) );
					const double neg = (GlobalC::ucell.magnet.lsign_
					                 && rho_in[1][ir]*GlobalC::ucell.magnet.ux_[0]+rho_in[2][ir]*GlobalC::ucell.magnet.ux_[1]+rho_in[3][ir]*GlobalC::ucell.magnet.ux_[2]<=0)
									 ? -1 : 1;
					if(amag > vanishing_charge)
					{
						for(int i=1;i<4;i++)
							v(i,ir) -= neg * 0.5 * (dh[0][ir]-dh[1][ir]) * rho_in[i][ir] / amag;
					}
				}
			}
		};

		//------------------------------------------------------------------
		// "func.info->family" includes LDA, GGA, Hybrid functional (HYB).
		//------------------------------------------------------------------
		switch( func.info->family )
		{
			case XC_FAMILY_LDA:
				// call Libxc function: xc_lda_exc_vxc
				xc_lda_exc_vxc( &func, GlobalC::pw.nrxx, rho.data(),
					exc.data(), vrho.data() );
				process_exc();
				process_vrho();
				break;
			case XC_FAMILY_GGA:
			case XC_FAMILY_HYB_GGA:
				// call Libxc function: xc_gga_exc_vxc
				xc_gga_exc_vxc( &func, GlobalC::pw.nrxx, rho.data(), sigma.data(),
					exc.data(), vrho.data(), vsigma.data() );
				process_exc();
				process_vrho();
				process_vsigma();
				break;
			default:
				throw std::domain_error("func.info->family ="+ModuleBase::GlobalFunc::TO_STRING(func.info->family)
					+" unfinished in "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
				break;
		}
		// Libxc function: Deallocate memory
		xc_func_end(&func);
	}

	//-------------------------------------------------
	// for MPI, reduce the exchange-correlation energy
	//-------------------------------------------------
	Parallel_Reduce::reduce_double_pool( etxc );
	Parallel_Reduce::reduce_double_pool( vtxc );

    etxc *= GlobalC::ucell.omega / GlobalC::pw.ncxyz;
    vtxc *= GlobalC::ucell.omega / GlobalC::pw.ncxyz;

    ModuleBase::timer::tick("Potential_Libxc","v_xc");
	return std::make_tuple( etxc, vtxc, std::move(v) );
}


//----------------------------------------------------------------------------
// for adding new xc functionals, only this function needs to be updated
// now we use iexch, igcx, icorr, igcc in xcf to characterize XC functionals
//----------------------------------------------------------------------------
std::vector<xc_func_type> Potential_Libxc::init_func()
{
	// 'funcs' is the return value
	std::vector<xc_func_type> funcs;

	const int xc_polarized = (1==nspin0()) ? XC_UNPOLARIZED : XC_POLARIZED;

	//-------------------------------------------
	// define a function named 'add_func', which 
	// will be called in the following codes
	//-------------------------------------------
	auto add_func = [&]( const int function )
	{
		funcs.push_back({});
		// 'xc_func_init' is defined in Libxc
		xc_func_init( &funcs.back(), function, xc_polarized );
	};

	//--------------------------------------
	// for the exchange energy 
	//--------------------------------------
	if( 6==GlobalC::xcf.iexch_now && 
		8==GlobalC::xcf.igcx_now && 
		4==GlobalC::xcf.icorr_now && 
		4==GlobalC::xcf.igcc_now ) // GGA functional
	{
		add_func( XC_HYB_GGA_XC_PBEH );		
		double parameter_hse[3] = { GlobalC::exx_global.info.hybrid_alpha, 
			GlobalC::exx_global.info.hse_omega, 
			GlobalC::exx_global.info.hse_omega };
		xc_func_set_ext_params(&funcs.back(), parameter_hse);
		return funcs;	
	}
	else if( 9==GlobalC::xcf.iexch_now && 
			12==GlobalC::xcf.igcx_now && 
			4==GlobalC::xcf.icorr_now && 
			4==GlobalC::xcf.igcc_now ) // HSE06 hybrid functional
	{
		add_func( XC_HYB_GGA_XC_HSE06 );	
		double parameter_hse[3] = { GlobalC::exx_global.info.hybrid_alpha, 
			GlobalC::exx_global.info.hse_omega, 
			GlobalC::exx_global.info.hse_omega };
		xc_func_set_ext_params(&funcs.back(), parameter_hse);
		return funcs;	
	}
	
	if( 1==GlobalC::xcf.iexch_now && 
		0==GlobalC::xcf.igcx_now ) // LDA functional
	{
		add_func( XC_LDA_X );
	}
	else if( 1==GlobalC::xcf.iexch_now && 
		3==GlobalC::xcf.igcx_now ) // GGA functional
	{
		add_func( XC_GGA_X_PBE );
	}
	else if(GlobalC::xcf.igcx_now == 13 ) //SCAN_X
	{
		add_func(263);
	}
	else
	{
		throw std::domain_error("iexch="+ModuleBase::GlobalFunc::TO_STRING(GlobalC::xcf.iexch_now)+", igcx="+ModuleBase::GlobalFunc::TO_STRING(GlobalC::xcf.igcx_now)
			+" unfinished in "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
	}

	//--------------------------------------
	// for the correlation energy part
	//--------------------------------------
	if( 1==GlobalC::xcf.icorr_now && 0==GlobalC::xcf.igcc_now )
	{
		add_func( XC_LDA_C_PZ );
	}
	else if( 4==GlobalC::xcf.icorr_now && 4==GlobalC::xcf.igcc_now )
	{
		add_func( XC_GGA_C_PBE );
	}
	else if (GlobalC::xcf.igcc_now == 9)
	{
		add_func( 267 );
	}
	else
	{
		throw std::domain_error("icorr="+ModuleBase::GlobalFunc::TO_STRING(GlobalC::xcf.icorr_now)+", igcc="+ModuleBase::GlobalFunc::TO_STRING(GlobalC::xcf.igcc_now)
			+" unfinished in "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
	}

	return funcs;
}



// [rho, sigma, gdr] = Potential_Libxc::cal_input(...)
std::tuple< std::vector<double>, 
			std::vector<double>, 
			std::vector<std::vector<ModuleBase::Vector3<double>>> > 
Potential_Libxc::cal_input( 
	const std::vector<xc_func_type> &funcs, 
	const double * const * const rho_in,
	const double * const rho_core_in )
{
	// ..., ↑_{i},↓_{i}, ↑_{i+1},↓_{i+1}, ...
	std::vector<double> rho;
	bool finished_rho   = false;

	// here we assume rho_core equally exists in different spins, may need double check
	auto cal_rho = [&]()
	{
		if(!finished_rho)
		{
			rho.resize(GlobalC::pw.nrxx*nspin0());
			if(nspin0()==1 || GlobalV::NSPIN==2)
			{
				for( size_t is=0; is!=nspin0(); ++is )
					for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
						rho[ir*nspin0()+is] = rho_in[is][ir] + 1.0/nspin0()*rho_core_in[ir];
			}
			else // may need updates for SOC
			{
				if(GlobalC::xcf.igcx||GlobalC::xcf.igcc)
				{
					GlobalC::ucell.cal_ux();
				}
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
				{
					const double amag = sqrt( pow(rho_in[1][ir],2) + pow(rho_in[2][ir],2) + pow(rho_in[3][ir],2) );
					const double neg = (GlobalC::ucell.magnet.lsign_ && rho_in[1][ir]*GlobalC::ucell.magnet.ux_[0]
					                    +rho_in[2][ir]*GlobalC::ucell.magnet.ux_[1]
					                    +rho_in[3][ir]*GlobalC::ucell.magnet.ux_[2]<=0)
					                   ? -1 : 1;
					rho[ir*2]   = 0.5 * (rho_in[0][ir] + neg * amag) + 0.5 * rho_core_in[ir];
					rho[ir*2+1] = 0.5 * (rho_in[0][ir] - neg * amag) + 0.5 * rho_core_in[ir];
				}
			}
		}
		finished_rho = true;
	};		
		
	// [...,↑_{i},↑_{i+1},...], [...,↓_{i},↓_{i+1},...]
	std::vector<std::vector<ModuleBase::Vector3<double>>> gdr;
	bool finished_gdr = false;
	auto cal_gdr = [&]()
	{
		if(!finished_gdr)
		{
			gdr.resize( nspin0() );
			for( size_t is=0; is!=nspin0(); ++is )
			{
				std::vector<double> rhor(GlobalC::pw.nrxx);
				for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
					rhor[ir] = rho[ir*nspin0()+is];

				//------------------------------------------
				// initialize the charge density arry in reciprocal space
				// bring electron charge density from real space to reciprocal space
				//------------------------------------------
				std::vector<std::complex<double>> rhog(GlobalC::pw.ngmc);
				GlobalC::CHR.set_rhog(rhor.data(), rhog.data());

				//-------------------------------------------
				// compute the gradient of charge density and
				// store the gradient in gdr[is]
				//-------------------------------------------
				gdr[is].resize(GlobalC::pw.nrxx);
				GGA_PW::grad_rho(rhog.data(), gdr[is].data());
			}
		}
		finished_gdr = true;
	};
	
	// ..., ↑↑_{i},↑↓_{i},↓↓_{i}, ↑↑_{i+1},↑↓_{i+1},↓↓_{i+1}, ...
	std::vector<double> sigma;
	bool finished_sigma = false;
	auto cal_sigma = [&]()
	{
		if(!finished_sigma)
		{
			sigma.resize( GlobalC::pw.nrxx * ((1==nspin0())?1:3) );

			if( 1==nspin0() )
			{
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
				{
					sigma[ir] = gdr[0][ir]*gdr[0][ir];
				}
			}
			else
			{
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
				{
					sigma[ir*3]   = gdr[0][ir]*gdr[0][ir];
					sigma[ir*3+1] = gdr[0][ir]*gdr[1][ir];
					sigma[ir*3+2] = gdr[1][ir]*gdr[1][ir];
				}
			}
		}
		finished_sigma = true;
	};

	for( const xc_func_type &func : funcs )
	{
		switch( func.info->family )
		{
			case XC_FAMILY_LDA:
				cal_rho();
				break;
			case XC_FAMILY_GGA:
			case XC_FAMILY_HYB_GGA:
			case XC_FAMILY_MGGA:
				cal_rho();
				cal_gdr();
				cal_sigma();
				break;
			default:
				throw std::domain_error("func.info->family ="+ModuleBase::GlobalFunc::TO_STRING(func.info->family)
					+" unfinished in "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
				break;
		}
	}
	
	return std::make_tuple( std::move(rho), std::move(sigma), std::move(gdr) );
}

#endif	//ifdef USE_LIBXC
