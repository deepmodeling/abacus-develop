//==========================================================
// AUTHOR : Peize Lin
// DATE :   2017-09-14
// UPDATE : 2021-02-28
//==========================================================

#ifdef USE_LIBXC

#include "potential_libxc.h"
#include "src_pw/global.h"
#include "src_lcao/global_fp.h"
#include "src_pw/tools.h"
#include "src_global/global_function.h"
#include "src_pw/xc_gga_pw.h"

void Potential_Libxc::v_xc(
	const double * const * const rho_in, // CHR.rho_core may be needed in future
    double &etxc,
    double &vtxc,
    matrix &v)
{
    TITLE("Potential_Libxc","v_xc");
    timer::tick("Potential_Libxc","v_xc");

    etxc = 0.0;
    vtxc = 0.0;
	v.zero_out();

	//----------------------------------------------------------
	// xc_func_type is defined in Libxc package
	// to understand the usage of xc_func_type,
	// use can check on website, for example:
	// https://www.tddft.org/programs/libxc/manual/libxc-5.1.x/
	//----------------------------------------------------------
	vector<xc_func_type> funcs = init_func();

	// the type of input_tmp is automatically set to 'tuple'
	// [rho, sigma, gdr] = cal_input( funcs, rho_in );
	const auto input_tmp = cal_input( funcs, rho_in );

	// obtain the electron charge density rho from 'get' function (first variable)
	const vector<double> &rho = std::get<0>(input_tmp);

	// obtain the sigma from 'get' function (second variable)
	const vector<double> &sigma = std::get<1>(input_tmp);

	for( xc_func_type &func : funcs )
	{
		vector<double> exc   ( pw.nrxx                       );
		vector<double> vrho  ( pw.nrxx * nspin0()            );
		vector<double> vsigma( pw.nrxx * ((1==nspin0())?1:3) );
		
		// cal etxc from rho, exc
		auto process_exc = [&](vector<double> &sgn)
		{
			for( size_t is=0; is!=nspin0(); ++is )
			{
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
				{
					etxc += e2 * exc[ir] * rho[ir*nspin0()+is] * sgn[ir*nspin0()+is];
				}
			}
		};

		// cal vtx, v from rho_in, vrho
		auto process_vrho = [&](vector<double> &sgn)
		{
			if(nspin0()==1 || NSPIN==2)
			{
				for( size_t is=0; is!=nspin0(); ++is )
				{
					for( size_t ir=0; ir!=pw.nrxx; ++ir )
					{
						const double v_tmp = e2 * vrho[ir*nspin0()+is] * sgn[ir*nspin0()+is];
						v(is,ir) += v_tmp;
						vtxc += v_tmp * rho_in[is][ir];
					}
				}
			}
			else // may need updates for SOC
			{
				constexpr double vanishing_charge = 1.0e-12;
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
				{
					vector<double> v_tmp(4);
					v_tmp[0] = e2 * (0.5 * (vrho[ir*2] + vrho[ir*2+1]));
					const double vs = 0.5 * (vrho[ir*2] - vrho[ir*2+1]);
					const double amag = sqrt( pow(rho_in[1][ir],2) 
						+ pow(rho_in[2][ir],2) 
						+ pow(rho_in[3][ir],2) );

					if(amag>vanishing_charge)
					{
						for(int ipol=1; ipol<4; ++ipol)
						{
							v_tmp[ipol] = e2 * vs * rho_in[ipol][ir] / amag;
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
		auto process_vsigma = [&](vector<double> &sgn)
		{
			const std::vector<std::vector<Vector3<double>>> &gdr = std::get<2>(input_tmp);
			
			vector<vector<Vector3<double>>> h( nspin0(), vector<Vector3<double>>(pw.nrxx) );
			if( 1==nspin0() )
			{
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
				{
					h[0][ir] = e2 * gdr[0][ir] * vsigma[ir] * 2.0 * sgn[ir];
				}
			}
			else
			{
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
				{
					h[0][ir] = e2 * (gdr[0][ir] * vsigma[ir*3  ] * 2.0 
						* sgn[ir*2  ] + gdr[1][ir] * vsigma[ir*3+1] * sgn[ir*2] * sgn[ir*2+1]);
					h[1][ir] = e2 * (gdr[1][ir] * vsigma[ir*3+2] * 2.0 
						* sgn[ir*2+1] + gdr[0][ir] * vsigma[ir*3+1] * sgn[ir*2] * sgn[ir*2+1]);
				}
			}

			// define two dimensional array dh [ nspin, pw.nrxx ]
			vector<vector<double>> dh(nspin0(), vector<double>(pw.nrxx));
			for( size_t is=0; is!=nspin0(); ++is )
			{
				GGA_PW::grad_dot( VECTOR_TO_PTR(h[is]), VECTOR_TO_PTR(dh[is]) );
			}


			for( size_t is=0; is!=nspin0(); ++is )
			{
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
				{
					vtxc -= dh[is][ir] * rho[ir*nspin0()+is];
				}
			}

			if(nspin0()==1 || NSPIN==2)
			{
				for( size_t is=0; is!=nspin0(); ++is )
				{
					for( size_t ir=0; ir!=pw.nrxx; ++ir )
					{
						v(is,ir) -= dh[is][ir];
					}
				}
			}
			else // may need updates for SOC
			{
				constexpr double vanishing_charge = 1.0e-12;
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
				{
					v(0,ir) -= 0.5 * (dh[0][ir] + dh[1][ir]);
					const double amag = sqrt( pow(rho_in[1][ir],2) + pow(rho_in[2][ir],2) + pow(rho_in[3][ir],2) );
					const double neg = (soc.lsign && rho_in[1][ir]*soc.ux[0]
						+rho_in[2][ir]*soc.ux[1]+rho_in[3][ir]*soc.ux[2]<=0) ? -1 : 1;
					if(amag > vanishing_charge)
					{
						for(int i=1;i<4;i++)
						{
							v(i,ir) -= neg * 0.5 * (dh[0][ir]-dh[1][ir]) * rho_in[i][ir] / amag;
						}
					}
				}
			}
		};
		
		// jiyy add for threshold
		constexpr double rho_threshold = 1E-6;
		constexpr double grho_threshold = 1E-10;
		xc_func_set_dens_threshold(&func, rho_threshold);

		// LPZ: Explain what is sgn
		vector<double> sgn( pw.nrxx * nspin0(), 1.0);
		if(nspin0()==2 && func.info->family != XC_FAMILY_LDA && func.info->kind==XC_CORRELATION)
		{
			for( size_t ir=0; ir!=pw.nrxx; ++ir )
			{
				if ( rho[ir*2]<rho_threshold || sqrt(abs(sigma[ir*3]))<grho_threshold ) 
				{
					sgn[ir*2] = 0.0;
				}
				if ( rho[ir*2+1]<rho_threshold || sqrt(abs(sigma[ir*3+2]))<grho_threshold ) 
				{
					sgn[ir*2+1] = 0.0;
				}
			}
		}

		//------------------------------------------------------------------
		// "func.info->family" includes LDA, GGA, Hybrid functional (HYB).
		//------------------------------------------------------------------
		switch( func.info->family )
		{
			case XC_FAMILY_LDA:
				// call Libxc function: xc_lda_exc_vxc
				xc_lda_exc_vxc( &func, pw.nrxx, 
					VECTOR_TO_PTR(rho), VECTOR_TO_PTR(exc), VECTOR_TO_PTR(vrho) );
				process_exc(sgn);
				process_vrho(sgn);
				break;
			case XC_FAMILY_GGA:
			case XC_FAMILY_HYB_GGA:
				// call Libxc function: xc_gga_exc_vxc
				xc_gga_exc_vxc( &func, pw.nrxx, 
					VECTOR_TO_PTR(rho), VECTOR_TO_PTR(sigma), VECTOR_TO_PTR(exc), 
					VECTOR_TO_PTR(vrho), VECTOR_TO_PTR(vsigma) );
				process_exc(sgn);
				process_vrho(sgn);
				process_vsigma(sgn);
				break;
			default:
				throw domain_error("func.info->family ="+TO_STRING(func.info->family)
					+" unfinished in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
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

    etxc *= ucell.omega / pw.ncxyz;
    vtxc *= ucell.omega / pw.ncxyz;

    timer::tick("Potential_Libxc","v_xc");
    return;
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
	if( 6==xcf.iexch_now && 
		8==xcf.igcx_now && 
		4==xcf.icorr_now && 
		4==xcf.igcc_now ) // GGA functional
	{
		add_func( XC_HYB_GGA_XC_PBEH );		
		double parameter_hse[3] = { exx_global.info.hybrid_alpha, 
			exx_global.info.hse_omega, 
			exx_global.info.hse_omega };
		xc_func_set_ext_params(&funcs.back(), parameter_hse);
		return funcs;	
	}
	else if( 9==xcf.iexch_now && 
			12==xcf.igcx_now && 
			4==xcf.icorr_now && 
			4==xcf.igcc_now ) // HSE06 hybrid functional
	{
		add_func( XC_HYB_GGA_XC_HSE06 );	
		double parameter_hse[3] = { exx_global.info.hybrid_alpha, 
			exx_global.info.hse_omega, 
			exx_global.info.hse_omega };
		xc_func_set_ext_params(&funcs.back(), parameter_hse);
		return funcs;	
	}
	
	if( 1==xcf.iexch_now && 
		0==xcf.igcx_now ) // LDA functional
	{
		add_func( XC_LDA_X );
	}
	else if( 1==xcf.iexch_now && 
		3==xcf.igcx_now ) // GGA functional
	{
		add_func( XC_GGA_X_PBE );
	}
	else
	{
		throw domain_error("iexch="+TO_STRING(xcf.iexch_now)+", igcx="
			+TO_STRING(xcf.igcx_now)+" unfinished in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
	}

	//--------------------------------------
	// for the correlation energy part
	//--------------------------------------
	if( 1==xcf.icorr_now && 0==xcf.igcc_now )
	{
		add_func( XC_LDA_C_PZ );
	}
	else if( 4==xcf.icorr_now && 4==xcf.igcc_now )
	{
		add_func( XC_GGA_C_PBE );
	}
	else
	{
		throw domain_error("icorr="+TO_STRING(xcf.icorr_now)+", igcc="
			+TO_STRING(xcf.igcc_now)+" unfinished in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
	}

	return funcs;
}



// [rho, sigma, gdr] = cal_input( funcs, rho_in )
std::tuple< std::vector<double>, 
			std::vector<double>, 
			std::vector<std::vector<Vector3<double>>> > 
Potential_Libxc::cal_input( 
	const std::vector<xc_func_type> &funcs, 
	const double * const * const rho_in )
{
	// ..., ↑_{i},↓_{i}, ↑_{i+1},↓_{i+1}, ...
	std::vector<double> rho;
	bool finished_rho   = false;

	// here we assume CHR.rho_core equally exists in different spins, may need double check
	auto cal_rho = [&]()
	{
		if(!finished_rho)
		{
			rho.resize(pw.nrxx*nspin0());
			if(nspin0()==1 || NSPIN==2)
			{
				for( size_t is=0; is!=nspin0(); ++is )
				{
					for( size_t ir=0; ir!=pw.nrxx; ++ir )
					{
						rho[ir*nspin0()+is] = rho_in[is][ir] + 1.0/nspin0()*CHR.rho_core[ir];
					}
				}
			}
			else // may need updates for SOC
			{
				if(xcf.igcx||xcf.igcc)
				{
					soc.cal_ux(ucell.ntype);
				}
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
				{
					const double amag = sqrt( pow(rho_in[1][ir],2) + pow(rho_in[2][ir],2) + pow(rho_in[3][ir],2) );
					const double neg = (soc.lsign && rho_in[1][ir]*soc.ux[0]
						+rho_in[2][ir]*soc.ux[1]
						+rho_in[3][ir]*soc.ux[2]<=0) ? -1 : 1;
					rho[ir*2]   = 0.5 * (rho_in[0][ir] + neg * amag) + 0.5 * CHR.rho_core[ir];
					rho[ir*2+1] = 0.5 * (rho_in[0][ir] - neg * amag) + 0.5 * CHR.rho_core[ir];
				}
			}
		}
		finished_rho = true;
	};		
		
	// [...,↑_{i},↑_{i+1},...], [...,↓_{i},↓_{i+1},...]
	std::vector<std::vector<Vector3<double>>> gdr;
	bool finished_gdr = false;
	auto cal_gdr = [&]()
	{
		if(!finished_gdr)
		{
			gdr.resize( nspin0() );
			for( size_t is=0; is!=nspin0(); ++is )
			{
				vector<double> rhor(pw.nrxx);
				for(int ir=0; ir<pw.nrxx; ++ir)
				{
					rhor[ir] = rho[ir*nspin0()+is];
				}

				//------------------------------------------
				// initialize the charge density arry in 
				// reciprocal space
				//------------------------------------------
				vector<complex<double>> rhog(pw.ngmc);

				//-------------------------------------------
				// bring electron charge density from real
				// space to reciprocal space
				//-------------------------------------------
				CHR.set_rhog(rhor.data(), rhog.data());

				//-------------------------------------------
				// compute the gradient of charge density and
				// store the gradient in gdr[is]
				//-------------------------------------------
				gdr[is].resize(pw.nrxx);
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
			sigma.resize( pw.nrxx * ((1==nspin0())?1:3) );

			if( 1==nspin0() )
			{
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
				{
					sigma[ir] = gdr[0][ir]*gdr[0][ir];
				}
			}
			else
			{
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
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
				cal_rho();
				cal_gdr();
				cal_sigma();
				break;
			default:
				throw domain_error("func.info->family ="
				+TO_STRING(func.info->family)+" unfinished in "
				+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
				break;
		}
	}
	
	return std::make_tuple( rho, sigma, gdr );
}

#endif	//ifdef USE_LIBXC
