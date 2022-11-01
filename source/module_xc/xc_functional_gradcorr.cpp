// This file contains subroutines realted to gradient calculations
// it contains 5 subroutines:
// 1. gradcorr, which calculates gradient correction
// 2. grad_wfc, which calculates gradient of wavefunction
//		it is used in stress_func_mgga.cpp
// 3. grad_rho, which calculates gradient of density
// 4. grad_dot, which calculates divergence of something
// 5. noncolin_rho, which diagonalizes the spin density matrix
//  and gives the spin up and spin down components of the charge.

#include "../src_pw/global.h"
#include "xc_functional.h"

// from gradcorr.f90
void XC_Functional::gradcorr(double &etxc, double &vtxc, ModuleBase::matrix &v, std::vector<double> &stress_gga, const bool is_stress)
{
	ModuleBase::TITLE("XC_Functional","gradcorr");
	
	if(func_type == 0 || func_type == 1) return; // none or LDA functional

	bool igcc_is_lyp = false;
	if( func_id[1] == XC_GGA_C_LYP) igcc_is_lyp = true;

	int nspin0 = GlobalV::NSPIN;
	if(GlobalV::NSPIN==4) nspin0 =1;
	if(GlobalV::NSPIN==4&&(GlobalV::DOMAG||GlobalV::DOMAG_Z)) nspin0 = 2;
	if(GlobalV::NSPIN==4) GlobalC::ucell.cal_ux();

	assert(nspin0>0);
	const double fac = 1.0/ nspin0;

	if(is_stress)
	{
		stress_gga.resize(9);
		for(int i=0;i<9;i++)
		{
			stress_gga[i] = 0.0;
		}
	}

	// doing FFT to get rho in G space: rhog1 
    GlobalC::CHR.set_rhog(GlobalC::CHR.rho[0], GlobalC::CHR.rhog[0]);
	if(GlobalV::NSPIN==2)//mohan fix bug 2012-05-28
	{
		GlobalC::CHR.set_rhog(GlobalC::CHR.rho[1], GlobalC::CHR.rhog[1]);
	}
    GlobalC::CHR.set_rhog(GlobalC::CHR.rho_core, GlobalC::CHR.rhog_core);
		
	// sum up (rho_core+rho) for each spin in real space
	// and reciprocal space.
	double* rhotmp1 = nullptr;
	double* rhotmp2 = nullptr;
	std::complex<double>* rhogsum1 = nullptr;
	std::complex<double>* rhogsum2 = nullptr;
	ModuleBase::Vector3<double>* gdr1 = nullptr;
	ModuleBase::Vector3<double>* gdr2 = nullptr;
	ModuleBase::Vector3<double>* h1 = nullptr;
	ModuleBase::Vector3<double>* h2 = nullptr;
	double* neg = nullptr;
	double** vsave = nullptr;
	double** vgg = nullptr;
	
	// for spin unpolarized case, 
	// calculate the gradient of (rho_core+rho) in reciprocal space.
	rhotmp1 = new double[GlobalC::rhopw->nrxx];
	rhogsum1 = new std::complex<double>[GlobalC::rhopw->npw];
	ModuleBase::GlobalFunc::ZEROS(rhotmp1, GlobalC::rhopw->nrxx);
	ModuleBase::GlobalFunc::ZEROS(rhogsum1, GlobalC::rhopw->npw);
	for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++) rhotmp1[ir] = GlobalC::CHR.rho[0][ir] + fac * GlobalC::CHR.rho_core[ir];
	for(int ig=0; ig<GlobalC::rhopw->npw; ig++) rhogsum1[ig] = GlobalC::CHR.rhog[0][ig] + fac * GlobalC::CHR.rhog_core[ig];

	gdr1 = new ModuleBase::Vector3<double>[GlobalC::rhopw->nrxx];
	if(!is_stress)	h1 = new ModuleBase::Vector3<double>[GlobalC::rhopw->nrxx];
	
	XC_Functional::grad_rho( rhogsum1 , gdr1, GlobalC::rhopw);

	// for spin polarized case;
	// calculate the gradient of (rho_core+rho) in reciprocal space.
	if(GlobalV::NSPIN==2)
	{
		rhotmp2 = new double[GlobalC::rhopw->nrxx];
		rhogsum2 = new std::complex<double>[GlobalC::rhopw->npw];
		ModuleBase::GlobalFunc::ZEROS(rhotmp2, GlobalC::rhopw->nrxx);
		ModuleBase::GlobalFunc::ZEROS(rhogsum2, GlobalC::rhopw->npw);
		for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++) rhotmp2[ir] = GlobalC::CHR.rho[1][ir] + fac * GlobalC::CHR.rho_core[ir];
		for(int ig=0; ig<GlobalC::rhopw->npw; ig++) rhogsum2[ig] = GlobalC::CHR.rhog[1][ig] + fac * GlobalC::CHR.rhog_core[ig];

		gdr2 = new ModuleBase::Vector3<double>[GlobalC::rhopw->nrxx];
		if(!is_stress) h2 = new ModuleBase::Vector3<double>[GlobalC::rhopw->nrxx];
		
		XC_Functional::grad_rho( rhogsum2 , gdr2, GlobalC::rhopw);
	}

	if(GlobalV::NSPIN == 4&&(GlobalV::DOMAG||GlobalV::DOMAG_Z))
	{
		ModuleBase::GlobalFunc::ZEROS(rhotmp1, GlobalC::rhopw->nrxx);
		ModuleBase::GlobalFunc::ZEROS(rhogsum1, GlobalC::rhopw->npw);
		rhotmp2 = new double[GlobalC::rhopw->nrxx];
		ModuleBase::GlobalFunc::ZEROS(rhotmp2, GlobalC::rhopw->nrxx);
 		neg = new double [GlobalC::rhopw->nrxx];
		ModuleBase::GlobalFunc::ZEROS(neg, GlobalC::rhopw->nrxx);
		if(!is_stress)
		{
			vsave = new double* [GlobalV::NSPIN];
			for(int is = 0;is<GlobalV::NSPIN;is++) {
				vsave[is]= new double [GlobalC::rhopw->nrxx];
				for(int ir =0;ir<GlobalC::rhopw->nrxx;ir++){
					vsave[is][ir] = v(is,ir);
					v(is,ir) = 0;
				}
			}
			vgg = new double* [nspin0];
			for(int is = 0;is<nspin0;is++)vgg[is] = new double[GlobalC::rhopw->nrxx];
		}
		noncolin_rho(rhotmp1,rhotmp2,neg);

		rhogsum2 = new std::complex<double>[GlobalC::rhopw->npw];
		ModuleBase::GlobalFunc::ZEROS(rhogsum2, GlobalC::rhopw->npw);
		GlobalC::CHR.set_rhog(rhotmp1, rhogsum1);
		GlobalC::CHR.set_rhog(rhotmp2, rhogsum2);
		for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
		{
			rhotmp2[ir] += fac * GlobalC::CHR.rho_core[ir];
			rhotmp1[ir] += fac * GlobalC::CHR.rho_core[ir];
		}
		for(int ig=0; ig<GlobalC::rhopw->npw; ig++)
		{
			rhogsum2[ig] += fac * GlobalC::CHR.rhog_core[ig];
			rhogsum1[ig] += fac * GlobalC::CHR.rhog_core[ig];
		}

		gdr2 = new ModuleBase::Vector3<double>[GlobalC::rhopw->nrxx];
		h2 = new ModuleBase::Vector3<double>[GlobalC::rhopw->nrxx];

		XC_Functional::grad_rho( rhogsum1 , gdr1, GlobalC::rhopw );
		XC_Functional::grad_rho( rhogsum2 , gdr2, GlobalC::rhopw );

	}
	
	const double epsr = 1.0e-6;
	const double epsg = 1.0e-10;

	double grho2a = 0.0;
	double grho2b = 0.0;
	double sxc = 0.0;
	double v1xc = 0.0;
	double v2xc = 0.0;
	double vtxcgc = 0.0;
	double etxcgc = 0.0;

	if(nspin0==1)
	{
		double segno;
		for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
		{
			const double arho = std::abs( rhotmp1[ir] );
			if(!is_stress) h1[ir].x = h1[ir].y = h1[ir].z = 0.0;

			if(arho > epsr)
			{
				grho2a = gdr1[ir].norm2();

				if( rhotmp1[ir] >= 0.0 ) segno = 1.0;
				if( rhotmp1[ir] < 0.0 ) segno = -1.0;
				if(func_type == 3 && is_stress) //the gradcorr part to stress of mGGA
				{
#ifdef USE_LIBXC
					double v3xc;
					double atau = GlobalC::CHR.kin_r[0][ir]/2.0;
					XC_Functional::tau_xc( arho, grho2a, atau, sxc, v1xc, v2xc, v3xc);
#endif
				}
				else
				{
					XC_Functional::gcxc( arho, grho2a, sxc, v1xc, v2xc);
				}
				if(is_stress)
				{
					double tt[3];
					tt[0] = gdr1[ir].x;
					tt[1] = gdr1[ir].y;
					tt[2] = gdr1[ir].z;
					for(int l = 0;l< 3;l++)
					{
						for(int m = 0;m< l+1;m++)
						{
							int ind = l*3 + m;
							stress_gga[ind] += tt[l] * tt[m] * ModuleBase::e2 * v2xc;
						}
					}
				}
				else
				{
					// first term of the gradient correction:
					// D(rho*Exc)/D(rho)
					v(0, ir) += ModuleBase::e2 * v1xc;
					// cout << "v    " << v(0, ir) << endl;
					
					// h contains
					// D(rho*Exc) / D(|grad rho|) * (grad rho) / |grad rho|
					h1[ir] = ModuleBase::e2 * v2xc * gdr1[ir];
					
					vtxcgc += ModuleBase::e2* v1xc * ( rhotmp1[ir] - GlobalC::CHR.rho_core[ir] );
					etxcgc += ModuleBase::e2* sxc  * segno;
				}
			} // end arho > epsr
		}
	}// end nspin0 == 1
	else // spin polarized case
	{
		for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
		{
			if(use_libxc)
			{
				double sxc, v1xcup, v1xcdw, v2xcup, v2xcdw, v2xcud;
				XC_Functional::gcxc_spin_libxc(rhotmp1[ir], rhotmp2[ir], gdr1[ir], gdr2[ir], 
					sxc, v1xcup, v1xcdw, v2xcup, v2xcdw, v2xcud);
				if(is_stress)
				{
					double tt1[3],tt2[3];
					{
						tt1[0] = gdr1[ir].x;
						tt1[1] = gdr1[ir].y;
						tt1[2] = gdr1[ir].z;
						tt2[0] = gdr2[ir].x;
						tt2[1] = gdr2[ir].y;
						tt2[2] = gdr2[ir].z;
					}
					for(int l = 0;l< 3;l++)
					{
						for(int m = 0;m< l+1;m++)
						{
							int ind = l*3 + m;
							stress_gga [ind] += ( tt1[l] * tt1[m] * v2xcup + 
									tt2[l] * tt2[m] * v2xcdw + 
									(tt1[l] * tt2[m] +
									tt2[l] * tt1[m] ) * v2xcud ) * ModuleBase::e2;
						}
					}
				}
				else
				{
					// first term of the gradient correction : D(rho*Exc)/D(rho)
					v(0,ir) += ModuleBase::e2 * v1xcup;
					v(1,ir) += ModuleBase::e2 * v1xcdw;
				
					// h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
					h1[ir] += ModuleBase::e2 * ( v2xcup * gdr1[ir] + v2xcud * gdr2[ir] );
					h2[ir] += ModuleBase::e2 * ( v2xcdw * gdr2[ir] + v2xcud * gdr1[ir] );

					vtxcgc = vtxcgc + ModuleBase::e2 * v1xcup * ( rhotmp1[ir] - GlobalC::CHR.rho_core[ir] * fac );
					vtxcgc = vtxcgc + ModuleBase::e2 * v1xcdw * ( rhotmp2[ir] - GlobalC::CHR.rho_core[ir] * fac );
					etxcgc = etxcgc + ModuleBase::e2 * sxc;
				}
			}
			else
			{
				double v1cup = 0.0;
				double v1cdw = 0.0;
				double v2cup = 0.0;
				double v2cdw = 0.0;
				double v1xup = 0.0;
				double v1xdw = 0.0;
				double v2xup = 0.0;
				double v2xdw = 0.0;
				double v2cud = 0.0;
				double v2c = 0.0;
				double sx = 0.0;
				double sc = 0.0;
				double rh = rhotmp1[ir] + rhotmp2[ir];
				grho2a = gdr1[ir].norm2();
				grho2b = gdr2[ir].norm2();
				XC_Functional::gcx_spin(rhotmp1[ir], rhotmp2[ir], grho2a, grho2b,
					sx, v1xup, v1xdw, v2xup, v2xdw);
				
				if(rh > epsr)
				{
					if(igcc_is_lyp)
					{
						ModuleBase::WARNING_QUIT("XC_Functional","igcc_is_lyp is not available now.");
					}
					else
					{
						double zeta = ( rhotmp1[ir] - rhotmp2[ir] ) / rh;
						if(GlobalV::NSPIN==4&&(GlobalV::DOMAG||GlobalV::DOMAG_Z)) zeta = fabs(zeta) * neg[ir];
						const double grh2 = (gdr1[ir]+gdr2[ir]).norm2();
						XC_Functional::gcc_spin(rh, zeta, grh2, sc, v1cup, v1cdw, v2c);
						v2cup = v2c;
						v2cdw = v2c;
						v2cud = v2c;
					}
				}
				else
				{
					sc = 0.0;
					v1cup = 0.0;
					v1cdw = 0.0;
					v2c = 0.0;
					v2cup = 0.0;
					v2cdw = 0.0;
					v2cud = 0.0;
				}

				if(is_stress)
				{
					double tt1[3],tt2[3];
					{
						tt1[0] = gdr1[ir].x;
						tt1[1] = gdr1[ir].y;
						tt1[2] = gdr1[ir].z;
						tt2[0] = gdr2[ir].x;
						tt2[1] = gdr2[ir].y;
						tt2[2] = gdr2[ir].z;
					}
					for(int l = 0;l< 3;l++)
					{
						for(int m = 0;m< l+1;m++)
						{
							int ind = l*3 + m;
							//    exchange
							stress_gga [ind] += tt1[l] * tt1[m] * ModuleBase::e2 * v2xup + 
									tt2[l] * tt2[m] * ModuleBase::e2 * v2xdw;
							//    correlation
							stress_gga [ind] += ( tt1[l] * tt1[m] * v2cup + 
									tt2[l] * tt2[m] * v2cdw + 
									(tt1[l] * tt2[m] +
									tt2[l] * tt1[m] ) * v2cud ) * ModuleBase::e2;
						}
					}
				}
				else
				{
					// first term of the gradient correction : D(rho*Exc)/D(rho)
					v(0,ir) = v(0,ir) + ModuleBase::e2 * ( v1xup + v1cup );
					v(1,ir) = v(1,ir) + ModuleBase::e2 * ( v1xdw + v1cdw );
				
					// h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
					h1[ir] = ModuleBase::e2 * ( ( v2xup + v2cup ) * gdr1[ir] + v2cud * gdr2[ir] );
					h2[ir] = ModuleBase::e2 * ( ( v2xdw + v2cdw ) * gdr2[ir] + v2cud * gdr1[ir] );

					vtxcgc = vtxcgc + ModuleBase::e2 * ( v1xup + v1cup ) * ( rhotmp1[ir] - GlobalC::CHR.rho_core[ir] * fac );
					vtxcgc = vtxcgc + ModuleBase::e2 * ( v1xdw + v1cdw ) * ( rhotmp2[ir] - GlobalC::CHR.rho_core[ir] * fac );
					etxcgc = etxcgc + ModuleBase::e2 * ( sx + sc );
				}
			}
		}// end ir

	}

	//std::cout << "\n vtxcgc=" << vtxcgc;
	//std::cout << "\n etxcgc=" << etxcgc << std::endl;

	if(!is_stress)
	{
		for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++) rhotmp1[ir] -= fac * GlobalC::CHR.rho_core[ir];
		if(nspin0==2) for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++) rhotmp2[ir] -= fac * GlobalC::CHR.rho_core[ir];
		
		// second term of the gradient correction :
		// \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )

		// dh is in real sapce.
		double* dh = new double[GlobalC::rhopw->nrxx];

		for(int is=0; is<nspin0; is++)
		{
			ModuleBase::GlobalFunc::ZEROS(dh, GlobalC::rhopw->nrxx);
			if(is==0)XC_Functional::grad_dot(h1,dh,GlobalC::rhopw);
			if(is==1)XC_Functional::grad_dot(h2,dh,GlobalC::rhopw);

			for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
				v(is, ir) -= dh[ir];
		
			double sum = 0.0;
			if(is==0)
				for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
					sum += dh[ir] * rhotmp1[ir];
			else if(is==1)
				for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
					sum += dh[ir] * rhotmp2[ir];
			
			vtxcgc -= sum;
		}
		
		delete[] dh;

		vtxc += vtxcgc;
		etxc += etxcgc;

		if(GlobalV::NSPIN == 4 && (GlobalV::DOMAG||GlobalV::DOMAG_Z))
		{
			for(int is=0;is<GlobalV::NSPIN;is++)
			{
				for(int ir=0;ir<GlobalC::rhopw->nrxx;ir++)
				{
					if(is<nspin0) vgg[is][ir] = v(is,ir);
					v(is,ir) = vsave[is][ir];
				}
			}
			for(int ir=0;ir<GlobalC::rhopw->nrxx;ir++)
			{
				v(0,ir) += 0.5 * (vgg[0][ir] + vgg[1][ir]);
				double amag = sqrt(pow(GlobalC::CHR.rho[1][ir],2)+pow(GlobalC::CHR.rho[2][ir],2)+pow(GlobalC::CHR.rho[3][ir],2));
				if(amag>1e-12)
				{
					for(int i=1;i<4;i++)
						v(i,ir)+= neg[ir] * 0.5 *(vgg[0][ir]-vgg[1][ir])*GlobalC::CHR.rho[i][ir]/amag;
				}
			}
		}
	}
	// deacllocate
	delete[] rhotmp1;
	delete[] rhogsum1;
	delete[] gdr1;
	if(!is_stress) delete[] h1;

	if(GlobalV::NSPIN==2)
	{
		delete[] rhotmp2;
		delete[] rhogsum2;
		delete[] gdr2;
		if(!is_stress) delete[] h2;
	}
	if(GlobalV::NSPIN == 4 && (GlobalV::DOMAG||GlobalV::DOMAG_Z))
	{
		delete[] neg;
		if(!is_stress) 
		{
			for(int i=0; i<nspin0; i++) delete[] vgg[i];
			delete[] vgg;
			for(int i=0; i<GlobalV::NSPIN; i++) delete[] vsave[i];
			delete[] vsave;
			delete[] h2;
		}
		delete[] rhotmp2;
		delete[] rhogsum2;
		delete[] gdr2;
	}

	return;
}

void XC_Functional::grad_wfc( const std::complex<double> *rhog, const int ik, std::complex<double> **grad, ModulePW::PW_Basis_K *wfc_basis)
{
	double *kplusg;
	const int npw_k = wfc_basis->npwk[ik];
	kplusg = new double[npw_k];
	ModuleBase::GlobalFunc::ZEROS(kplusg, npw_k);

	std::complex<double> *Porter = new std::complex<double> [wfc_basis->nmaxgr];

	for(int ipol=0; ipol<3; ipol++)
	{
		// the formula is : rho(r)^prime = \int iG * rho(G)e^{iGr} dG
		for(int ig=0; ig<npw_k; ++ig)
			kplusg[ig] = wfc_basis->getgpluskcar(ik,ig)[ipol] * GlobalC::ucell.tpiba;

		// calculate the charge density gradient in reciprocal space.
		for(int ig=0; ig<npw_k; ig++)
			Porter[ig] = complex<double>(0.0,kplusg[ig]) * rhog[ig];

		// bring the gdr from G --> R
		wfc_basis->recip2real(Porter, Porter, ik);

		for(int ir=0; ir<wfc_basis->nrxx; ++ir)
			grad[ir][ipol]= Porter[ir];
	}//end loop ipol
	delete[] kplusg;
	delete[] Porter;
	return;
}

void XC_Functional::grad_rho( const std::complex<double> *rhog, ModuleBase::Vector3<double> *gdr, ModulePW::PW_Basis *rho_basis)
{
	std::complex<double> *gdrtmp = new std::complex<double>[rho_basis->nmaxgr];
	ModuleBase::GlobalFunc::ZEROS(gdrtmp, rho_basis->nmaxgr);

	

	// the formula is : rho(r)^prime = \int iG * rho(G)e^{iGr} dG
	for(int i = 0 ; i < 3 ; ++i)
	{
		// calculate the charge density gradient in reciprocal space.
		for(int ig=0; ig<rho_basis->npw; ig++)
			gdrtmp[ig] = ModuleBase::IMAG_UNIT * rhog[ig] * rho_basis->gcar[ig][i];

		// bring the gdr from G --> R
		rho_basis->recip2real(gdrtmp, gdrtmp);

		// remember to multily 2pi/a0, which belongs to G vectors.
		for(int ir=0; ir<rho_basis->nrxx; ir++)
			gdr[ir][i] = gdrtmp[ir].real() * GlobalC::ucell.tpiba;
	}

	delete[] gdrtmp;
	return;
}


void XC_Functional::grad_dot(const ModuleBase::Vector3<double> *h, double *dh, ModulePW::PW_Basis *rho_basis)
{
	std::complex<double> *aux = new std::complex<double>[rho_basis->nmaxgr];
	std::complex<double> *gaux = new std::complex<double>[rho_basis->npw];
	ModuleBase::GlobalFunc::ZEROS(gaux, rho_basis->npw);

	for(int i = 0 ; i < 3 ; ++i)
	{
		for(int ir = 0; ir < rho_basis->nrxx; ++ir)
			aux[ir] = std::complex<double>( h[ir][i], 0.0);

		// bring to G space.
		rho_basis->real2recip(aux,aux);
		for(int ig = 0; ig < rho_basis->npw; ++ig)
			gaux[ig] +=  ModuleBase::IMAG_UNIT * aux[ig] * rho_basis->gcar[ig][i];
	}

	// bring back to R space
	rho_basis->recip2real(gaux,aux);

	for(int ir=0; ir<rho_basis->nrxx; ir++)
		dh[ir] = aux[ir].real() * GlobalC::ucell.tpiba;
	
	delete[] aux;	
	delete[] gaux;
	return;
}

void XC_Functional::noncolin_rho(double *rhoout1,double *rhoout2, double *neg)
{
	//this function diagonalizes the spin density matrix and gives as output the
	//spin up and spin down components of the charge.
	//If lsign is true up and dw are with respect to the fixed quantization axis 
	//ux, otherwise rho + |m| is always rhoup and rho-|m| is always rhodw.
	double amag=0;
	for(int ir = 0;ir<GlobalC::rhopw->nrxx;ir++) neg[ir] = 1.0;
	if(GlobalC::ucell.magnet.lsign_)
	{
		for(int ir = 0;ir<GlobalC::rhopw->nrxx;ir++)
		{
			if(GlobalC::CHR.rho[1][ir]*GlobalC::ucell.magnet.ux_[0] + GlobalC::CHR.rho[2][ir]*GlobalC::ucell.magnet.ux_[1] + GlobalC::CHR.rho[3][ir]*GlobalC::ucell.magnet.ux_[2]>0) neg[ir] = 1.0;
			else neg[ir] = -1.0;
		}
	}
	for(int ir = 0;ir<GlobalC::rhopw->nrxx;ir++)
	{
		amag = sqrt(pow(GlobalC::CHR.rho[1][ir],2)+pow(GlobalC::CHR.rho[2][ir],2)+pow(GlobalC::CHR.rho[3][ir],2));
		rhoout1[ir] = 0.5 * (GlobalC::CHR.rho[0][ir] + neg[ir] * amag);
		rhoout2[ir] = 0.5 * (GlobalC::CHR.rho[0][ir] - neg[ir] * amag);
	}
	return;
}

