#ifdef USE_LIBXC

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

// calculating laplacian of rho
std::vector<double> XC_Functional_Libxc::cal_lapl(
	const int nspin,
	const std::size_t nrxx,
	const std::vector<double> &rho,
	const double tpiba,
	const Charge* const chr)
{
	std::vector<double> lapl(nrxx * nspin);
	for (int is = 0; is != nspin; ++is)
	{
		std::vector<double> rhor(nrxx);
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 1024)
		#endif
		for (std::size_t ir = 0; ir < nrxx; ++ir)
		{
			rhor[ir] = rho[ir * nspin + is];
		}
		std::vector<std::complex<double>> rhog(chr->rhopw->npw);
		chr->rhopw->real2recip(rhor.data(), rhog.data());

		std::vector<double> lapl_spin(nrxx);
		XC_Functional::laplacian_rho(rhog.data(), lapl_spin.data(), chr->rhopw, tpiba);

		#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 1024)
		#endif
		for (std::size_t ir = 0; ir < nrxx; ++ir)
		{
			lapl[ir * nspin + is] = lapl_spin[ir];
		}
	}
	return lapl;
}

// calculating laplacian of rho using FD kernel (consistent with vlapl potential)
std::vector<double> XC_Functional_Libxc::cal_lapl_fd(
	const int nspin,
	const std::size_t nrxx,
	const std::vector<double> &rho,
	const Charge* const chr)
{
	const int nx = chr->rhopw->nx;
	const int ny = chr->rhopw->ny;
	const int nz = chr->rhopw->nz;
	const double tpiba = chr->rhopw->tpiba;
	const double tpiba2 = tpiba * tpiba;
	const double twopi = 2.0 * M_PI;

	std::vector<double> gg_fd(chr->rhopw->npw);
	#ifdef _OPENMP
	#pragma omp parallel for schedule(static, 1024)
	#endif
	for (int ig = 0; ig < chr->rhopw->npw; ++ig)
	{
		int m0 = chr->rhopw->gdirect[ig].x;
		int m1 = chr->rhopw->gdirect[ig].y;
		int m2 = chr->rhopw->gdirect[ig].z;

		double fd00 = 2.0 * nx * nx * (1.0 - std::cos(twopi * m0 / nx));
		double fd11 = 2.0 * ny * ny * (1.0 - std::cos(twopi * m1 / ny));
		double fd22 = 2.0 * nz * nz * (1.0 - std::cos(twopi * m2 / nz));
		double fd01 = nx * ny * std::sin(twopi * m0 / nx) * std::sin(twopi * m1 / ny);
		double fd02 = nx * nz * std::sin(twopi * m0 / nx) * std::sin(twopi * m2 / nz);
		double fd12 = ny * nz * std::sin(twopi * m1 / ny) * std::sin(twopi * m2 / nz);

		double ggt00 = chr->rhopw->GGT.e11;
		double ggt01 = chr->rhopw->GGT.e12;
		double ggt02 = chr->rhopw->GGT.e13;
		double ggt11 = chr->rhopw->GGT.e22;
		double ggt12 = chr->rhopw->GGT.e23;
		double ggt22 = chr->rhopw->GGT.e33;

		gg_fd[ig] = (ggt00 * fd00 + ggt11 * fd11 + ggt22 * fd22
		           + 2.0 * ggt01 * fd01 + 2.0 * ggt02 * fd02 + 2.0 * ggt12 * fd12)
		          / (twopi * twopi);
	}

	std::vector<double> lapl(nrxx * nspin);
	for (int is = 0; is != nspin; ++is)
	{
		std::vector<double> rhor(nrxx);
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 1024)
		#endif
		for (std::size_t ir = 0; ir < nrxx; ++ir)
		{
			rhor[ir] = rho[ir * nspin + is];
		}
		std::vector<std::complex<double>> rhog(chr->rhopw->npw);
		chr->rhopw->real2recip(rhor.data(), rhog.data());

		std::vector<std::complex<double>> lapl_g(chr->rhopw->nmaxgr);
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 1024)
		#endif
		for (int ig = 0; ig < chr->rhopw->npw; ++ig)
		{
			lapl_g[ig] = -rhog[ig] * gg_fd[ig] * tpiba2;
		}
		for (int ig = chr->rhopw->npw; ig < chr->rhopw->nmaxgr; ++ig)
		{
			lapl_g[ig] = 0.0;
		}

		std::vector<std::complex<double>> lapl_r(chr->rhopw->nmaxgr);
		chr->rhopw->recip2real(lapl_g.data(), lapl_r.data());

		#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 1024)
		#endif
		for (std::size_t ir = 0; ir < nrxx; ++ir)
		{
			lapl[ir * nspin + is] = lapl_r[ir].real();
		}
	}
	return lapl;
}

// calculating Hessian of rho: d²rho / dr_alpha dr_beta
// returns 6 independent components per spin: xx, yy, zz, xy, yz, zx
std::vector<std::vector<double>> XC_Functional_Libxc::cal_rho_hessian(
	const int nspin,
	const std::size_t nrxx,
	const std::vector<double> &rho,
	const Charge* const chr)
{
	// 6 independent components: (0)xx (1)yy (2)zz (3)xy (4)yz (5)zx
	std::vector<std::vector<double>> hessian(nspin, std::vector<double>(6 * nrxx, 0.0));

	for (int is = 0; is != nspin; ++is)
	{
		std::vector<double> rhor(nrxx);
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 1024)
		#endif
		for (std::size_t ir = 0; ir < nrxx; ++ir)
		{
			rhor[ir] = rho[ir * nspin + is];
		}
		std::vector<std::complex<double>> rhog(chr->rhopw->npw);
		chr->rhopw->real2recip(rhor.data(), rhog.data());

		// Hessian in G-space: H_ab(G) = -G_a * G_b * rho(G)
		// We compute each of the 6 independent components
		int ab_map[6][2] = {{0,0}, {1,1}, {2,2}, {0,1}, {1,2}, {0,2}};

		for (int ab = 0; ab < 6; ++ab)
		{
			int a = ab_map[ab][0];
			int b = ab_map[ab][1];
			std::vector<std::complex<double>> hess_g(chr->rhopw->npw);

			#ifdef _OPENMP
			#pragma omp parallel for schedule(static, 1024)
			#endif
			for (int ig = 0; ig < chr->rhopw->npw; ++ig)
			{
				hess_g[ig] = -rhog[ig] * chr->rhopw->gcar[ig][a] * chr->rhopw->gcar[ig][b];
			}

			std::vector<std::complex<double>> hess_r(chr->rhopw->nmaxgr);
			chr->rhopw->recip2real(hess_g.data(), hess_r.data());

			const double tpiba2 = chr->rhopw->tpiba * chr->rhopw->tpiba;
			#ifdef _OPENMP
			#pragma omp parallel for schedule(static, 1024)
			#endif
			for (std::size_t ir = 0; ir < nrxx; ++ir)
			{
				hessian[is][ab * nrxx + ir] = hess_r[ir].real() * tpiba2;
			}
		}
	}
	return hessian;
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
