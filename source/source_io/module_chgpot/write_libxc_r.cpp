//======================
// AUTHOR : Peize Lin
// DATE :   2024-09-12
//======================

#ifdef USE_LIBXC

#include "write_libxc_r.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_hamilt/module_xc/xc_functional_libxc.h"
#include "source_estate/module_charge/charge.h"
#include "source_basis/module_pw/pw_basis_big.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_io/cube_io.h"
#include "source_base/global_variable.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"

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
	const Charge &chr,
	const ModulePW::PW_Basis_Big &pw_big,
	const ModulePW::PW_Basis &pw_rhod)
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
		rho = XC_Functional_Libxc::convert_rho(nspin, nrxx, &chr);
	}
	else
	{
		std::tuple<std::vector<double>,std::vector<double>> rho_amag = XC_Functional_Libxc::convert_rho_amag_nspin4(nspin, nrxx, &chr);
		rho = std::get<0>(std::move(rho_amag));
		amag = std::get<1>(std::move(rho_amag));
	}

	std::vector<double> sigma;
	if(is_gga)
	{
		const std::vector<std::vector<ModuleBase::Vector3<double>>> gdr = XC_Functional_Libxc::cal_gdr(nspin, nrxx, rho, tpiba, &chr);
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
					case 4:		xc_gga_lxc ( &func, nrxx, rho.data(), sigma.data(), v4rho4.data(), v4rho3sigma.data(), v4rho2sigma2.data(), v4rhosigma3.data(), v4sigma4.data() );
					case 3:		xc_gga_kxc ( &func, nrxx, rho.data(), sigma.data(), v3rho3.data(), v3rho2sigma.data(), v3rhosigma2.data(), v3sigma3.data() );
					case 2:		xc_gga_fxc ( &func, nrxx, rho.data(), sigma.data(), v2rho2.data(), v2rhosigma.data(), v2sigma2.data() );
					case 1:		xc_gga_exc_vxc( &func, nrxx, rho.data(), sigma.data(), exc.data(), vrho.data(), vsigma.data() );
								break;
					case 0:		xc_gga_exc ( &func, nrxx, rho.data(), sigma.data(), exc.data() );
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

	auto write_data = [&pw_big, &pw_rhod](
		const std::string data_name,
		const std::vector<double> &data,
		const int number_spin)
	{
		for(int is=0; is<number_spin; ++is)
		{
			std::ofstream ofs;
			if(GlobalV::MY_RANK==0)
			{
				const std::string file_name = PARAM.globalv.global_out_dir + "xc_"+data_name+"_s"+std::to_string(is+1)+".cube";
				ofs.open(file_name);

				ofs.unsetf(std::ostream::fixed);
				ofs << std::setprecision(PARAM.inp.out_xc_r[1]);
				ofs << std::scientific;
			}

		  #ifdef __MPI
			ModuleIO::write_cube_core(
				ofs,
				pw_big.bz,
				pw_big.nbz,
				pw_rhod.nplane * number_spin,
				pw_rhod.startz_current,
				data.data()+is,
				pw_rhod.nx * pw_rhod.ny,
				pw_rhod.nz,
				number_spin,
				pw_rhod.nz);
		  #else
			if(nspin!=1)
				{ throw std::invalid_argument("nspin="+std::to_string(nspin)+" is invalid for ModuleIO::write_cube_core without MPI. see "+std:string(__FILE__)+" line "+std::to_string(__LINE__)); }
			ModuleIO::write_cube_core(
				ofs,
				data.data()+is,
				pw_rhod.nx * pw_rhod.ny,
				pw_rhod.nz,
				pw_rhod.nz);
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


#ifdef __MPI

void ModuleIO::write_cube_core(
	std::ofstream &ofs_cube,
	const int bz,
	const int nbz,
	const int nplane,
	const int startz_current,
	const double*const data,
	const int nxy,
	const int nz,
	const int nld,
	const int n_data_newline)
{
	ModuleBase::TITLE("ModuleIO", "write_cube_core");

	const int my_rank = GlobalV::MY_RANK;
	const int my_pool = GlobalV::MY_POOL;
	const int rank_in_pool = GlobalV::RANK_IN_POOL;
	const int nproc_in_pool = GlobalV::NPROC_IN_POOL;

	// only do in the first pool.
	if (my_pool == 0)
	{
		/// for cube file
		const int nxyz = nxy * nz;
		std::vector<double> data_cube(nxyz, 0.0);

		// num_z: how many planes on processor 'ip'
		std::vector<int> num_z(nproc_in_pool, 0);
		for (int iz = 0; iz < nbz; iz++)
		{
			const int ip = iz % nproc_in_pool;
			num_z[ip] += bz;
		}

		// start_z: start position of z in
		// processor ip.
		std::vector<int> start_z(nproc_in_pool, 0);
		for (int ip = 1; ip < nproc_in_pool; ip++)
		{
			start_z[ip] = start_z[ip - 1] + num_z[ip - 1];
		}

		// which_ip: found iz belongs to which ip.
		std::vector<int> which_ip(nz, 0);
		for (int iz = 0; iz < nz; iz++)
		{
			for (int ip = 0; ip < nproc_in_pool; ip++)
			{
				if (iz >= start_z[nproc_in_pool - 1])
				{
					which_ip[iz] = nproc_in_pool - 1;
					break;
				}
				else if (iz >= start_z[ip] && iz < start_z[ip + 1])
				{
					which_ip[iz] = ip;
					break;
				}
			}
		}

		int count = 0;
		std::vector<double> zpiece(nxy, 0.0);

		// save the rho one z by one z.
		for (int iz = 0; iz < nz; iz++)
		{
			zpiece.assign(nxy, 0.0);

			// tag must be different for different iz.
			const int tag = iz;
			MPI_Status ierror;

			// case 1: the first part of rho in processor 0.
			if (which_ip[iz] == 0 && rank_in_pool == 0)
			{
				for (int ixy = 0; ixy < nxy; ixy++)
				{
					// mohan change to rho_save on 2012-02-10
					// because this can make our next restart calculation lead
					// to the same scf_thr as the one saved.
					zpiece[ixy] = data[ixy * nplane + (iz - startz_current) * nld];
				}
			}
			// case 2: > first part rho: send the rho to
			// processor 0.
			else if (which_ip[iz] == rank_in_pool)
			{
				for (int ixy = 0; ixy < nxy; ixy++)
				{
					zpiece[ixy] = data[ixy * nplane + (iz - startz_current) * nld];
				}
				MPI_Send(zpiece.data(), nxy, MPI_DOUBLE, 0, tag, POOL_WORLD);
			}

			// case 2: > first part rho: processor 0 receive the rho
			// from other processors
			else if (rank_in_pool == 0)
			{
				MPI_Recv(zpiece.data(), nxy, MPI_DOUBLE, which_ip[iz], tag, POOL_WORLD, &ierror);
			}

			if (my_rank == 0)
			{
				/// for cube file
				for (int ixy = 0; ixy < nxy; ixy++)
				{
					data_cube[ixy * nz + iz] = zpiece[ixy];
				}
				/// for cube file
			}
		} // end iz

		// for cube file
		if (my_rank == 0)
		{
			for (int ixy = 0; ixy < nxy; ixy++)
			{
				for (int iz = 0; iz < nz; iz++)
				{
					ofs_cube << " " << data_cube[ixy * nz + iz];
					if ((iz % n_data_newline == n_data_newline-1) && (iz != nz - 1))
					{
						ofs_cube << "\n";
					}
				}
				ofs_cube << "\n";
			}
		}
		/// for cube file
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

#else   // #ifdef __MPI

void ModuleIO::write_cube_core(
	std::ofstream &ofs_cube,
	const double*const data,
	const int nxy,
	const int nz,
	const int n_data_newline)
{
	ModuleBase::TITLE("ModuleIO", "write_cube_core");
	for (int ixy = 0; ixy < nxy; ixy++)
	{
		for (int iz = 0; iz < nz; iz++)
		{
			ofs_cube << " " << data[iz * nxy + ixy];
			// ++count_cube;
			if ((iz % n_data_newline == n_data_newline-1) && (iz != nz - 1))
			{
				ofs_cube << "\n";
			}
		}
		ofs_cube << "\n";
	}
}

#endif  // #ifdef __MPI

#endif // USE_LIBXC