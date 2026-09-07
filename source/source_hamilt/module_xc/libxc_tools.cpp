#ifdef __LIBXC

#include "libxc_abacus.h"
#include "xc_functional.h"
#include "source_estate/module_charge/charge.h"
#include "source_io/module_parameter/parameter.h"

// converting rho (abacus=>libxc)
std::vector<double> XC_Functional_Libxc::convert_rho(
	const int nspin,
	const std::size_t nrxx,
	const Charge* const chr)
{
	std::vector<double> rho(nrxx*nspin);
	#ifdef _OPENMP
	#pragma omp parallel for collapse(2) schedule(static, 1024)
	#endif
	for( int is=0; is<nspin; ++is )
	{
		for( int ir=0; ir<nrxx; ++ir )
		{
			rho[ir*nspin+is] = chr->rho[is][ir] + 1.0/nspin*chr->rho_core[ir];
		}
	}
	return rho;
}

// converting rho (abacus=>libxc)
std::tuple<std::vector<double>, std::vector<double>>
XC_Functional_Libxc::convert_rho_amag_nspin4(
	const int nspin,
	const std::size_t nrxx,
	const Charge* const chr)
{
	assert(PARAM.inp.nspin==4);
	std::vector<double> rho(nrxx*nspin);
	std::vector<double> amag(nrxx);
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for( int ir=0; ir<nrxx; ++ir )
	{
		const double arhox = std::abs( chr->rho[0][ir] + chr->rho_core[ir] );
		amag[ir] = std::sqrt( std::pow(chr->rho[1][ir],2)
							+ std::pow(chr->rho[2][ir],2)
							+ std::pow(chr->rho[3][ir],2) );
		const double amag_clip = (amag[ir]<arhox) ? amag[ir] : arhox;
		rho[ir*nspin+0] = (arhox + amag_clip) / 2.0;
		rho[ir*nspin+1] = (arhox - amag_clip) / 2.0;
	}
	return std::make_tuple(std::move(rho), std::move(amag));
}

// calculating grho
std::vector<std::vector<ModuleBase::Vector3<double>>>
XC_Functional_Libxc::cal_gdr(
	const int nspin,
	const std::size_t nrxx,
	const std::vector<double> &rho,
	const double tpiba,
	const Charge* const chr)
{
	std::vector<std::vector<ModuleBase::Vector3<double>>> gdr(nspin);
	for( int is=0; is!=nspin; ++is )
	{
		std::vector<double> rhor(nrxx);
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 1024)
		#endif
		for(std::size_t ir=0; ir<nrxx; ++ir)
		{
			rhor[ir] = rho[ir*nspin+is];
		}
		//------------------------------------------
		// initialize the charge density array in reciprocal space
		// bring electron charge density from real space to reciprocal space
		//------------------------------------------
		std::vector<std::complex<double>> rhog(chr->rhopw->npw);
		chr->rhopw->real2recip(rhor.data(), rhog.data());

		//-------------------------------------------
		// compute the gradient of charge density and
		// store the gradient in gdr[is]
		//-------------------------------------------
		gdr[is].resize(nrxx);
		XC_Functional::grad_rho(rhog.data(), gdr[is].data(), chr->rhopw, tpiba);
	} // end for(is)
	return gdr;
}

void XC_Functional_Libxc::cal_gdr_and_lapl(
	const int nspin,
	const std::size_t nrxx,
	const std::vector<double> &rho,
	const double tpiba,
	const Charge* const chr,
	std::vector<std::vector<ModuleBase::Vector3<double>>> &gdr,
	std::vector<double> &lapl,
	const bool need_laplacian)
{
	gdr.resize(nspin);
	lapl.assign(nrxx * nspin, 0.0);
	for( int is=0; is!=nspin; ++is )
	{
		std::vector<double> rhor(nrxx);
		for(std::size_t ir=0; ir<nrxx; ++ir)
			rhor[ir] = rho[ir*nspin+is];
		std::vector<std::complex<double>> rhog(chr->rhopw->npw);
		chr->rhopw->real2recip(rhor.data(), rhog.data());
		gdr[is].resize(nrxx);
		XC_Functional::grad_rho(rhog.data(), gdr[is].data(), chr->rhopw, tpiba);
		if (need_laplacian)
		{
			std::vector<double> lapl_spin(nrxx);
			XC_Functional::laplacian_rho(rhog.data(), lapl_spin.data(), chr->rhopw, tpiba);
			for(std::size_t ir=0; ir<nrxx; ++ir)
				lapl[ir*nspin+is] = lapl_spin[ir];
		}
	}
}

// converting grho (abacus=>libxc)
std::vector<double> XC_Functional_Libxc::convert_sigma(
	const std::vector<std::vector<ModuleBase::Vector3<double>>> &gdr)
{
	const std::size_t nspin = gdr.size();
	assert(nspin>0);
	const std::size_t nrxx = gdr[0].size();
	for(std::size_t is=1; is<nspin; ++is)
	{
		assert(nrxx==gdr[is].size());
	}

	std::vector<double> sigma( nrxx * ((1==nspin)?1:3) );
	if( 1==nspin )
	{
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 1024)
		#endif
		for( std::size_t ir=0; ir<nrxx; ++ir )
			sigma[ir] = gdr[0][ir]*gdr[0][ir];
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 256)
		#endif
		for( std::size_t ir=0; ir<nrxx; ++ir )
		{
			sigma[ir*3]   = gdr[0][ir]*gdr[0][ir];
			sigma[ir*3+1] = gdr[0][ir]*gdr[1][ir];
			sigma[ir*3+2] = gdr[1][ir]*gdr[1][ir];
		}
	}
	return sigma;
}

// sgn for threshold mask
std::vector<double> XC_Functional_Libxc::cal_sgn(
	const double rho_threshold,
	const double grho_threshold,
	const xc_func_type &func,
	const int nspin,
	const std::size_t nrxx,
	const std::vector<double> &rho,
	const std::vector<double> &sigma)
{
    //assert(nrxx>0); // adding this once will cause error in examples
	std::vector<double> sgn(nrxx*nspin, 1.0);
	// in the case of GGA correlation for polarized case,
	// a cutoff for grho is required to ensure that libxc gives reasonable results
	if(nspin==2 && func.info->family != XC_FAMILY_LDA && func.info->kind==XC_CORRELATION)
	{
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 512)
		#endif
		for( int ir=0; ir<nrxx; ++ir )
		{
			if ( rho[ir*2]<rho_threshold || std::sqrt(std::abs(sigma[ir*3]))<grho_threshold )
			{
				sgn[ir*2] = 0.0;
			}
			if ( rho[ir*2+1]<rho_threshold || std::sqrt(std::abs(sigma[ir*3+2]))<grho_threshold )
			{
				sgn[ir*2+1] = 0.0;
			}
		}
	}
	return sgn;
}

// converting etxc from exc (libxc=>abacus)
double XC_Functional_Libxc::convert_etxc(
	const int nspin,
	const std::size_t nrxx,
	const std::vector<double> &sgn,
	const std::vector<double> &rho,
	std::vector<double> exc)
{
	double etxc = 0.0;
	#ifdef _OPENMP
	#pragma omp parallel for collapse(2) reduction(+:etxc) schedule(static, 256)
	#endif
	for( int is=0; is<nspin; ++is )
	{
		for( int ir=0; ir<nrxx; ++ir )
		{
			etxc += ModuleBase::e2 * exc[ir] * rho[ir*nspin+is] * sgn[ir*nspin+is];
		}
	}
	return etxc;
}

XC_Functional_Libxc::LibxcWeightedDerivatives
XC_Functional_Libxc::make_libxc_weighted_derivatives(
	const xc_func_type &func,
	const int nspin,
	const std::size_t nrxx,
	const std::vector<double> &sgn,
	const std::vector<double> &rho,
	const std::vector<double> &sigma,
	const std::vector<double> &exc,
	const std::vector<double> &vrho,
	const std::vector<double> &vsigma)
{
	assert(nspin == 1 || nspin == 2);
	assert(sgn.size() == nrxx * nspin);
	assert(rho.size() == nrxx * nspin);
	assert(exc.size() == nrxx);
	assert(vrho.size() == nrxx * nspin);
	assert(func.nspin == nspin);

	const bool is_gga
		= func.info->family == XC_FAMILY_GGA || func.info->family == XC_FAMILY_HYB_GGA;
	const std::size_t nsigma = nspin == 1 ? 1 : 3;
	if (is_gga)
	{
		assert(sigma.size() == nrxx * nsigma);
		assert(vsigma.size() == nrxx * nsigma);
	}

	LibxcWeightedDerivatives weighted;
	weighted.energy_sum = 0.0;
	weighted.drho.assign(nrxx * nspin, 0.0);
	if (is_gga)
	{
		weighted.dsigma.assign(nrxx * nsigma, 0.0);
	}

	const double density_floor = func.dens_threshold;
	const double sigma_floor = func.sigma_threshold * func.sigma_threshold;
	double energy_sum = 0.0;
	#ifdef _OPENMP
	#pragma omp parallel for reduction(+:energy_sum) schedule(static, 512)
	#endif
	for (std::size_t ir = 0; ir < nrxx; ++ir)
	{
		double raw_density_sum = 0.0;
		double sanitized_density_sum = 0.0;
		double energy_weight = 0.0;
		for (int is = 0; is < nspin; ++is)
		{
			const std::size_t index = ir * nspin + is;
			raw_density_sum += rho[index];
			sanitized_density_sum += std::max(density_floor, rho[index]);
			energy_weight += sgn[index] * rho[index];
		}

		// Libxc leaves all outputs zero below the total-density threshold.
		// Inside that branch, the ABACUS weighted energy is locally constant.
		if (raw_density_sum < density_floor)
		{
			continue;
		}

		// ABACUS accumulates M*eps while Libxc differentiates Y*eps after
		// y_s=max(T,rho_s). Hence d eps/d y_s=(vrho_s-eps)/Y.
		energy_sum += energy_weight * exc[ir];
		const double libxc_weight = energy_weight / sanitized_density_sum;
		for (int is = 0; is < nspin; ++is)
		{
			const std::size_t index = ir * nspin + is;
			const double floor_jacobian = rho[index] > density_floor ? 1.0 : 0.0;
			weighted.drho[index]
				= sgn[index] * exc[ir]
				  + libxc_weight * floor_jacobian * (vrho[index] - exc[ir]);
		}

		if (!is_gga)
		{
			continue;
		}

		if (nspin == 1)
		{
			const double floor_jacobian = sigma[ir] > sigma_floor ? 1.0 : 0.0;
			weighted.dsigma[ir] = libxc_weight * floor_jacobian * vsigma[ir];
			continue;
		}

		const std::size_t sigma_index = 3 * ir;
		const double sigma_uu = sigma[sigma_index];
		const double sigma_ud = sigma[sigma_index + 1];
		const double sigma_dd = sigma[sigma_index + 2];
		const double jacobian_uu = sigma_uu > sigma_floor ? 1.0 : 0.0;
		const double jacobian_dd = sigma_dd > sigma_floor ? 1.0 : 0.0;
		const double sanitized_uu = std::max(sigma_floor, sigma_uu);
		const double sanitized_dd = std::max(sigma_floor, sigma_dd);
		const double cross_limit = 0.5 * (sanitized_uu + sanitized_dd);

		double cross_to_diagonal = 0.0;
		double cross_jacobian = 1.0;
		if (sigma_ud < -cross_limit)
		{
			cross_to_diagonal = -0.5;
			cross_jacobian = 0.0;
		}
		else if (sigma_ud > cross_limit)
		{
			cross_to_diagonal = 0.5;
			cross_jacobian = 0.0;
		}

		weighted.dsigma[sigma_index]
			= libxc_weight * jacobian_uu
			  * (vsigma[sigma_index] + cross_to_diagonal * vsigma[sigma_index + 1]);
		weighted.dsigma[sigma_index + 1]
			= libxc_weight * cross_jacobian * vsigma[sigma_index + 1];
		weighted.dsigma[sigma_index + 2]
			= libxc_weight * jacobian_dd
			  * (vsigma[sigma_index + 2] + cross_to_diagonal * vsigma[sigma_index + 1]);
	}
	weighted.energy_sum = energy_sum;
	return weighted;
}

// converting vtxc and v from vrho and vsigma (libxc=>abacus)
std::pair<double,ModuleBase::matrix> XC_Functional_Libxc::convert_vtxc_v(
	const xc_func_type &func,
	const int nspin,
	const std::size_t nrxx,
	const std::vector<double> &sgn,
	const std::vector<double> &rho,
	const std::vector<std::vector<ModuleBase::Vector3<double>>> &gdr,
	const std::vector<double> &vrho,
	const std::vector<double> &vsigma,
	const double tpiba,
	const Charge* const chr)
{
    // assert(nrxx>0); // will cause error
	double vtxc = 0.0;
	ModuleBase::matrix v(nspin, nrxx);

	#ifdef _OPENMP
	#pragma omp parallel for collapse(2) reduction(+:vtxc) schedule(static, 256)
	#endif
	for( int is=0; is<nspin; ++is )
	{
		for( std::size_t ir=0; ir<nrxx; ++ir )
		{
            const std::size_t index = ir*nspin+is;
			const double v_tmp = ModuleBase::e2 * vrho[index] * sgn[index];
			v(is,ir) += v_tmp;
			vtxc += v_tmp * rho[index];
		}
	}

	if(func.info->family == XC_FAMILY_GGA || func.info->family == XC_FAMILY_HYB_GGA)
	{
		const std::vector<std::vector<double>> dh = XC_Functional_Libxc::cal_dh(nspin, nrxx, sgn, gdr, vsigma, tpiba, chr);

		double rvtxc = 0.0;
		#ifdef _OPENMP
		#pragma omp parallel for collapse(2) reduction(+:rvtxc) schedule(static, 256)
		#endif
		for( int is=0; is<nspin; ++is )
		{
			for( std::size_t ir=0; ir<nrxx; ++ir )
			{
				rvtxc += dh[is][ir] * rho[ir*nspin+is];
				v(is,ir) -= dh[is][ir];
			}
		}

		vtxc -= rvtxc;
	} // end if(func.info->family == XC_FAMILY_GGA || func.info->family == XC_FAMILY_HYB_GGA))

	return std::make_pair(vtxc, std::move(v));
}


// dh for gga v
std::vector<std::vector<double>> XC_Functional_Libxc::cal_dh(
	const int nspin,
	const std::size_t nrxx,
	const std::vector<double> &sgn,
	const std::vector<std::vector<ModuleBase::Vector3<double>>> &gdr,
	const std::vector<double> &vsigma,
	const double tpiba,
	const Charge* const chr)
{
    //assert(nrxx>0); // this line will cause bug
	std::vector<std::vector<ModuleBase::Vector3<double>>> h(
		nspin,
		std::vector<ModuleBase::Vector3<double>>(nrxx) );

	if( nspin==1 )
	{
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 1024)
		#endif
		for( std::size_t ir=0; ir<nrxx; ++ir )
		{
			h[0][ir] = 2.0 * gdr[0][ir] * vsigma[ir] * 2.0 * sgn[ir];
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 1024)
		#endif
		for( std::size_t ir=0; ir< nrxx; ++ir )
		{
			h[0][ir] = 2.0 * (gdr[0][ir] * vsigma[ir*3  ] * sgn[ir*2  ] * 2.0
							+ gdr[1][ir] * vsigma[ir*3+1] * sgn[ir*2]   * sgn[ir*2+1]);
			h[1][ir] = 2.0 * (gdr[1][ir] * vsigma[ir*3+2] * sgn[ir*2+1] * 2.0
							+ gdr[0][ir] * vsigma[ir*3+1] * sgn[ir*2]   * sgn[ir*2+1]);
		}
	}

	// define two dimensional array dh [ nspin, nrxx ]
	std::vector<std::vector<double>> dh(nspin, std::vector<double>(nrxx));
	for( int is=0; is!=nspin; ++is )
	{
		XC_Functional::grad_dot( h[is].data(), dh[is].data(), chr->rhopw, tpiba);
	}

	return dh;
}


// convert v for NSPIN=4
ModuleBase::matrix XC_Functional_Libxc::convert_v_nspin4(
	const std::size_t nrxx,
	const Charge* const chr,
	const std::vector<double> &amag,
	const ModuleBase::matrix &v)
{
    //assert(nrxx>0);
	assert(PARAM.inp.nspin==4);
	constexpr double vanishing_charge = 1.0e-10;
	ModuleBase::matrix v_nspin4(PARAM.inp.nspin, nrxx);
	for( int ir=0; ir<nrxx; ++ir )
	{
		v_nspin4(0,ir) = 0.5 * (v(0,ir)+v(1,ir));
	}
	if(PARAM.globalv.domag || PARAM.globalv.domag_z)
	{
		for( int ir=0; ir<nrxx; ++ir )
		{
			if ( amag[ir] > vanishing_charge )
			{
				const double vs = 0.5 * (v(0,ir)-v(1,ir));
				for(int ipol=1; ipol<PARAM.inp.nspin; ++ipol)
				{
					v_nspin4(ipol,ir) = vs * chr->rho[ipol][ir] / amag[ir];
				}
			}
		}
	}
	return v_nspin4;
}

#endif
