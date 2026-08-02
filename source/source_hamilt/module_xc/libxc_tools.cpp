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

std::vector<double> XC_Functional_Libxc::compute_mag_part_nspin4(
	const std::size_t nrxx,
	const Charge* const chr)
{
	std::vector<double> mag_part(3 * nrxx, 0.0);
	#ifdef _OPENMP
	#pragma omp parallel for schedule(static, 1024)
	#endif
	for (std::size_t ir = 0; ir < nrxx; ++ir)
	{
		double mx = chr->rho[1][ir], my = chr->rho[2][ir], mz = chr->rho[3][ir];
		double amag = std::sqrt(mx * mx + my * my + mz * mz);
		if (amag > 1e-12)
		{
			mag_part[ir] = mx / amag;
			mag_part[ir + nrxx] = my / amag;
			mag_part[ir + 2 * nrxx] = mz / amag;
		}
	}
	return mag_part;
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
		if(PARAM.inp.nspin==4 && (PARAM.inp.gga_grad == 2 || PARAM.inp.gga_grad == 3) && (PARAM.globalv.domag || PARAM.globalv.domag_z))
		{
			std::vector<double> mag_part_tmp = XC_Functional_Libxc::compute_mag_part_nspin4(nrxx, chr);
			const std::vector<std::vector<double>> dh = XC_Functional_Libxc::cal_dh_sf(nspin, nrxx, sgn, gdr, vsigma, mag_part_tmp, tpiba, chr);

			constexpr int nspin4 = 4;
			double rvtxc = 0.0;
			#ifdef _OPENMP
			#pragma omp parallel for collapse(2) reduction(+:rvtxc) schedule(static, 256)
			#endif
			for (int is = 0; is < nspin4; ++is)
			{
				for (std::size_t ir = 0; ir < nrxx; ++ir)
				{
					double rho_ir = 0.0;
					if (is == 0) { rho_ir = chr->rho[0][ir]; }
					else { rho_ir = chr->rho[is][ir]; }
					rvtxc += dh[is][ir] * rho_ir;
					v(is, ir) -= dh[is][ir];
				}
			}
			vtxc -= rvtxc;
		}
		else
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
		}
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

std::vector<std::vector<ModuleBase::Vector3<double>>> XC_Functional_Libxc::cal_gdr_sf(
	const int nspin,
	const std::size_t nrxx,
	const std::vector<double> &rho,
	const std::vector<double> &mag_part,
	const double tpiba,
	const Charge* const chr)
{
	std::vector<std::vector<ModuleBase::Vector3<double>>> gdr(nspin);
	std::vector<double> rhor(nrxx);
	std::vector<std::complex<double>> rhog(chr->rhopw->npw);
	std::vector<ModuleBase::Vector3<double>> gdr_tmp(nrxx);

	#ifdef _OPENMP
	#pragma omp parallel for schedule(static, 1024)
	#endif
	for (std::size_t ir = 0; ir < nrxx; ++ir)
	{
		rhor[ir] = rho[ir * nspin + 0];
	}
	chr->rhopw->real2recip(rhor.data(), rhog.data());
	gdr[0].resize(nrxx);
	XC_Functional::grad_rho(rhog.data(), gdr[0].data(), chr->rhopw, tpiba);

	gdr[1].resize(nrxx);
	#ifdef _OPENMP
	#pragma omp parallel for schedule(static, 1024)
	#endif
	for (std::size_t ir = 0; ir < nrxx; ++ir)
	{
		gdr_tmp[ir] = gdr[0][ir];
		gdr[0][ir] = 0.5 * gdr_tmp[ir];
		gdr[1][ir] = 0.5 * gdr_tmp[ir];
	}

	for (int is = 1; is <= 3; ++is)
	{
		chr->rhopw->real2recip(chr->rho[is], rhog.data());
		XC_Functional::grad_rho(rhog.data(), gdr_tmp.data(), chr->rhopw, tpiba);
		const double* mp = mag_part.data() + (is - 1) * nrxx;
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 1024)
		#endif
		for (std::size_t ir = 0; ir < nrxx; ++ir)
		{
			const ModuleBase::Vector3<double> g = 0.5 * gdr_tmp[ir] * mp[ir];
			gdr[0][ir] += g;
			gdr[1][ir] -= g;
		}
	}

	return gdr;
}

std::vector<std::vector<double>> XC_Functional_Libxc::cal_dh_sf(
	const int nspin,
	const std::size_t nrxx,
	const std::vector<double> &sgn,
	const std::vector<std::vector<ModuleBase::Vector3<double>>> &gdr,
	const std::vector<double> &vsigma,
	const std::vector<double> &mag_part,
	const double tpiba,
	const Charge* const chr)
{
	std::vector<ModuleBase::Vector3<double>> h1(nrxx), h2(nrxx);
	#ifdef _OPENMP
	#pragma omp parallel for schedule(static, 1024)
	#endif
	for (std::size_t ir = 0; ir < nrxx; ++ir)
	{
		h1[ir] = 2.0 * (gdr[0][ir] * vsigma[ir*3  ] * sgn[ir*2  ] * 2.0
					  + gdr[1][ir] * vsigma[ir*3+1] * sgn[ir*2]   * sgn[ir*2+1]);
		h2[ir] = 2.0 * (gdr[1][ir] * vsigma[ir*3+2] * sgn[ir*2+1] * 2.0
					  + gdr[0][ir] * vsigma[ir*3+1] * sgn[ir*2]   * sgn[ir*2+1]);
	}

	std::vector<ModuleBase::Vector3<double>> tmp_h(nrxx);
	std::vector<double> dh0(nrxx), dh_mu(nrxx);

	#ifdef _OPENMP
	#pragma omp parallel for schedule(static, 1024)
	#endif
	for (std::size_t ir = 0; ir < nrxx; ++ir)
	{
		tmp_h[ir] = 0.5 * (h1[ir] + h2[ir]);
	}
	XC_Functional::grad_dot(tmp_h.data(), dh0.data(), chr->rhopw, tpiba);

	std::vector<std::vector<double>> dh_total(4, std::vector<double>(nrxx, 0.0));
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
	for (std::size_t ir = 0; ir < nrxx; ++ir)
	{
		dh_total[0][ir] = dh0[ir];
	}

	for (int mu = 1; mu < 4; ++mu)
	{
		const double* mp_mu = mag_part.data() + (mu - 1) * nrxx;
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 1024)
		#endif
		for (std::size_t ir = 0; ir < nrxx; ++ir)
		{
			tmp_h[ir] = 0.5 * (h1[ir] - h2[ir]) * mp_mu[ir];
		}
		XC_Functional::grad_dot(tmp_h.data(), dh_mu.data(), chr->rhopw, tpiba);
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 1024)
		#endif
		for (std::size_t ir = 0; ir < nrxx; ++ir)
		{
			dh_total[mu][ir] = dh_mu[ir];
		}
	}

	return dh_total;
}

ModuleBase::matrix XC_Functional_Libxc::convert_v_nspin4_sf(
	const std::size_t nrxx,
	const Charge* const chr,
	const std::vector<double> &mag_part,
	const ModuleBase::matrix &v)
{
	assert(PARAM.inp.nspin==4);
	constexpr double vanishing_charge = 1.0e-10;
	ModuleBase::matrix v_nspin4(PARAM.inp.nspin, nrxx);

	for (std::size_t ir = 0; ir < nrxx; ++ir)
	{
		v_nspin4(0, ir) = 0.5 * (v(0, ir) + v(1, ir));
	}

	if (PARAM.globalv.domag || PARAM.globalv.domag_z)
	{
		for (std::size_t ir = 0; ir < nrxx; ++ir)
		{
			double amag = std::sqrt(std::pow(chr->rho[1][ir], 2)
								  + std::pow(chr->rho[2][ir], 2)
								  + std::pow(chr->rho[3][ir], 2));
			if (amag > vanishing_charge)
			{
				const double vs = 0.5 * (v(0, ir) - v(1, ir));
				for (int ipol = 1; ipol < PARAM.inp.nspin; ++ipol)
				{
					v_nspin4(ipol, ir) = vs * mag_part[ir + (ipol - 1) * nrxx];
				}
			}
		}
	}
	return v_nspin4;
}

#endif
