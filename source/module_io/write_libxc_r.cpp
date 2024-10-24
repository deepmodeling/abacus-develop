//======================
// AUTHOR : Peize Lin
// DATE :   2024-09-12
//======================

#ifdef USE_LIBXC

#include "write_libxc_r.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_general/module_xc/xc_functional_libxc.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_esolver/esolver_fp.h"
#include "module_io/cube_io.h"
#include "module_base/global_variable.h"
#include "module_parameter/parameter.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"

#include <xc.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <unistd.h>

void ModuleIO::write_libxc_r(
	const int order,
	const std::vector<int> &func_id,
	const int &nrxx, // number of real-space grid
	const double &omega, // volume of cell
	const double tpiba,
	const Charge*const chr,
	const ModuleESolver::ESolver_FP*const esolver)
{
	ModuleBase::TITLE("ModuleIO","write_libxc_r");
	ModuleBase::timer::tick("ModuleIO","write_libxc_r");

	const int nspin =
		(PARAM.inp.nspin == 1 || ( PARAM.inp.nspin ==4 && !PARAM.globalv.domag && !PARAM.globalv.domag_z))
		? 1 : 2;

	//----------------------------------------------------------
	// xc_func_type is defined in Libxc package
	// to understand the usage of xc_func_type,
	// use can check on website, for example:
	// https://www.tddft.org/programs/libxc/manual/libxc-5.1.x/
	//----------------------------------------------------------

	std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func( func_id, (1==nspin) ? XC_UNPOLARIZED:XC_POLARIZED );

	const bool is_gga = [&funcs]()
	{
		for( xc_func_type &func : funcs )
		{
			switch( func.info->family )
			{
				case XC_FAMILY_GGA:
				case XC_FAMILY_HYB_GGA:
					return true;
			}
		}
		return false;
	}();

	// converting rho
	std::vector<double> rho;
	std::vector<double> amag;
	if(1==nspin || 2==PARAM.inp.nspin)
	{
		rho = XC_Functional_Libxc::convert_rho(nspin, nrxx, chr);
	}
	else
	{
		std::tuple<std::vector<double>,std::vector<double>> rho_amag = XC_Functional_Libxc::convert_rho_amag_nspin4(nspin, nrxx, chr);
		rho = std::get<0>(std::move(rho_amag));
		amag = std::get<1>(std::move(rho_amag));
	}

	std::vector<double> sigma;
	if(is_gga)
	{
		const std::vector<std::vector<ModuleBase::Vector3<double>>> gdr = XC_Functional_Libxc::cal_gdr(nspin, nrxx, rho, tpiba, chr);
		sigma = XC_Functional_Libxc::convert_sigma(gdr);
	}

	std::vector<double> exc;
	std::vector<double> vrho;
	std::vector<double> vsigma;
	std::vector<double> v2rho2;
	std::vector<double> v2rhosigma;
	std::vector<double> v2sigma2;
	std::vector<double> v3rho3;
	std::vector<double> v3rho2sigma;
	std::vector<double> v3rhosigma2;
	std::vector<double> v3sigma3;
	std::vector<double> v4rho4;
	std::vector<double> v4rho3sigma;
	std::vector<double> v4rho2sigma2;
	std::vector<double> v4rhosigma3;
	std::vector<double> v4sigma4;
	// attention: order 4321 don't break
	switch( order )
	{
		case 4:		v4rho4.resize( nrxx * ((1==nspin)?1:5) );
		case 3:		v3rho3.resize( nrxx * ((1==nspin)?1:4) );
		case 2:		v2rho2.resize( nrxx * ((1==nspin)?1:3) );
		case 1:		vrho  .resize( nrxx * nspin            );
		case 0:		exc   .resize( nrxx                    );
					break;
		default:	throw std::domain_error("order ="+std::to_string(order)
						+" unfinished in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
					break;
	}
	if(is_gga)
	{
		switch( order )
		{
			case 4:		v4rho3sigma .resize( nrxx * ((1==nspin)?1:12) );
						v4rho2sigma2.resize( nrxx * ((1==nspin)?1:15) );
						v4rhosigma3 .resize( nrxx * ((1==nspin)?1:20) );
						v4sigma4    .resize( nrxx * ((1==nspin)?1:15) );
			case 3:		v3rho2sigma .resize( nrxx * ((1==nspin)?1:9)  );
						v3rhosigma2 .resize( nrxx * ((1==nspin)?1:12) );
						v3sigma3    .resize( nrxx * ((1==nspin)?1:10) );
			case 2:		v2rhosigma  .resize( nrxx * ((1==nspin)?1:6)  );
						v2sigma2    .resize( nrxx * ((1==nspin)?1:6)  );
			case 1:		vsigma      .resize( nrxx * ((1==nspin)?1:3)  );
			case 0:		break;
			default:	throw std::domain_error("order ="+std::to_string(order)
							+" unfinished in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
						break;
		}
	}

	for( xc_func_type &func : funcs )
	{
		// jiyy add for threshold
		constexpr double rho_threshold = 1E-6;
		constexpr double grho_threshold = 1E-10;

		xc_func_set_dens_threshold(&func, rho_threshold);

		// sgn for threshold mask
		const std::vector<double> sgn = XC_Functional_Libxc::cal_sgn(rho_threshold, grho_threshold, func, nspin, nrxx, rho, sigma);

		// call libxc function
		// attention: order 432 don't break
		switch( func.info->family )
		{
			case XC_FAMILY_LDA:
			{
				switch( order )
				{
					case 4:		xc_lda_lxc    ( &func, nrxx, rho.data(), v4rho4.data() );
					case 3:		xc_lda_kxc    ( &func, nrxx, rho.data(), v3rho3.data() );
					case 2:		xc_lda_fxc    ( &func, nrxx, rho.data(), v2rho2.data() );
					case 1:		xc_lda_exc_vxc( &func, nrxx, rho.data(), exc.data(), vrho.data() );
								break;
					case 0:		xc_lda_exc    ( &func, nrxx, rho.data(), exc.data() );
								break;
					default:	throw std::domain_error("order ="+std::to_string(order)
									+" unfinished in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
								break;
				}
				break;
			}
			case XC_FAMILY_GGA:
			case XC_FAMILY_HYB_GGA:
			{
				switch( order )
				{
					case 4:		xc_gga_lxc    ( &func, nrxx, rho.data(), sigma.data(), v4rho4.data(), v4rho3sigma.data(), v4rho2sigma2.data(), v4rhosigma3.data(), v4sigma4.data() );
					case 3:		xc_gga_kxc    ( &func, nrxx, rho.data(), sigma.data(), v3rho3.data(), v3rho2sigma.data(), v3rhosigma2.data(), v3sigma3.data() );
					case 2:		xc_gga_fxc    ( &func, nrxx, rho.data(), sigma.data(), v2rho2.data(), v2rhosigma.data(), v2sigma2.data() );
					case 1:		xc_gga_exc_vxc( &func, nrxx, rho.data(), sigma.data(), exc.data(), vrho.data(), vsigma.data() );
								break;
					case 0:		xc_gga_exc    ( &func, nrxx, rho.data(), sigma.data(), exc.data() );
								break;
					default:	throw std::domain_error("order ="+std::to_string(order)
									+" unfinished in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
								break;
				}
				break;
			}
			default:
			{
				throw std::domain_error("func.info->family ="+std::to_string(func.info->family)
					+" unfinished in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
				break;
			}
		} // end switch( func.info->family )
	} // end for( xc_func_type &func : funcs )

	auto write_data = [&esolver](
		const std::string data_name,
		const std::vector<double> &data,
		const int number_spin)
	{
		for(int is=0; is<number_spin; ++is)
		{
			std::ofstream ofs;
			if(GlobalV::MY_RANK==0)
			{
				const std::string folder_name = PARAM.globalv.global_out_dir + "/xc_r/";
				const std::string command0 =  "test -d " + folder_name + " || mkdir " + folder_name;
				system( command0.c_str() );
				const std::string file_name = folder_name + data_name+"_"+std::to_string(is);
				ofs.open(file_name);

				ofs.unsetf(std::ostream::fixed);
				ofs << std::setprecision(PARAM.inp.out_xc_r[1]);
				ofs << std::scientific;
			}

		  #ifdef __MPI
			ModuleIO::write_cube_core(
				ofs,
				esolver->pw_big->bz,
				esolver->pw_big->nbz,
				esolver->pw_rhod->nplane * number_spin,
				esolver->pw_rhod->startz_current,
				data.data()+is,
				esolver->pw_rhod->nx * esolver->pw_rhod->ny,
				esolver->pw_rhod->nz,
				number_spin,
				esolver->pw_rhod->nz);
		  #else
		  	if(nspin!=1)
				{ throw std::invalid_argument("nspin="+std::to_string(nspin)+" is invalid for ModuleIO::write_cube_core without MPI. see "+std:string(__FILE__)+" line "+std::to_string(__LINE__)); }
			ModuleIO::write_cube_core(
				ofs,
				data.data()+is,
				esolver->pw_rhod->nx * esolver->pw_rhod->ny,
				esolver->pw_rhod->nz,
				esolver->pw_rhod->nz);
		  #endif
		}
	};

	write_data( "rho", rho, nspin );

	if(1!=nspin && 2!=PARAM.inp.nspin)
		write_data( "amag", amag, 1 );

	if(is_gga)
		write_data( "sigma", sigma, (1==nspin)?1:3 );

	switch( order )
	{
		case 4:		write_data( "v4rho4", v4rho4, (1==nspin)?1:5 );
		case 3:		write_data( "v3rho3", v3rho3, (1==nspin)?1:4 );
		case 2:		write_data( "v2rho2", v2rho2, (1==nspin)?1:3 );
		case 1:		write_data( "vrho"  , vrho  , nspin );
		case 0:		write_data( "exc"   , exc   , 1 );
					break;
		default:	throw std::domain_error("order ="+std::to_string(order)
						+" unfinished in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
					break;
	}
	if(is_gga)
	{
		switch( order )
		{
			case 4:		write_data( "v4rho3sigma" , v4rho3sigma , (1==nspin)?1:12 );
						write_data( "v4rho2sigma2", v4rho2sigma2, (1==nspin)?1:15 );
						write_data( "v4rhosigma3" , v4rhosigma3 , (1==nspin)?1:20 );
						write_data( "v4sigma4"    , v4sigma4    , (1==nspin)?1:15 );
			case 3:		write_data( "v3rho2sigma" , v3rho2sigma , (1==nspin)?1:9  );
						write_data( "v3rhosigma2" , v3rhosigma2 , (1==nspin)?1:12 );
						write_data( "v3sigma3"    , v3sigma3    , (1==nspin)?1:10 );
			case 2:		write_data( "v2rhosigma"  , v2rhosigma  , (1==nspin)?1:6  );
						write_data( "v2sigma2"    , v2sigma2    , (1==nspin)?1:6  );
			case 1:		write_data( "vsigma"      , vsigma      , (1==nspin)?1:3  );
			case 0:		break;
			default:	throw std::domain_error("order ="+std::to_string(order)
							+" unfinished in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
						break;
		}
	}

	XC_Functional_Libxc::finish_func(funcs);

	ModuleBase::timer::tick("ModuleIO","write_libxc_r");
}

#endif // USE_LIBXC