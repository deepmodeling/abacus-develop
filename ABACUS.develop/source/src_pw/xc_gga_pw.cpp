#include "xc_gga_pw.h"
#include "global.h"
#include "xc_functional.h"
#include "../src_pw/myfunc.h"

// from gradcorr.f90
void GGA_PW::gradcorr(double &etxc, double &vtxc, matrix &v)
{
	TITLE("GGA_PW","gradcorr");
	
	if (xcf.igcx_now == 0  &&  xcf.igcc_now == 0)
	{
		return;
	}

	bool igcc_is_lyp = false;
	if( xcf.igcc_now == 3 || xcf.igcc_now == 7)
	{
		igcc_is_lyp = true;
	}

	int nspin0 = NSPIN;
	if(NSPIN==4) nspin0 =1;
	if(NSPIN==4&&(DOMAG||DOMAG_Z)) nspin0 = 2;
	if(NSPIN==4)
	{
		if(xcf.igcx != 0  ||  xcf.igcc != 0) mag.cal_ux(ucell.ntype);
	}

	assert(nspin0>0);
	const double fac = 1.0/ nspin0;

	// doing FFT to get rho in G space: rhog1 
    CHR.set_rhog(CHR.rho[0], CHR.rhog[0]);
	if(NSPIN==2)//mohan fix bug 2012-05-28
	{
		CHR.set_rhog(CHR.rho[1], CHR.rhog[1]);
	}
    CHR.set_rhog(CHR.rho_core, CHR.rhog_core);
		
	// sum up (rho_core+rho) for each spin in real space
	// and reciprocal space.
	double* rhotmp1 = nullptr;
	double* rhotmp2 = nullptr;
	complex<double>* rhogsum1 = nullptr;
	complex<double>* rhogsum2 = nullptr;
	Vector3<double>* gdr1 = nullptr;
	Vector3<double>* gdr2 = nullptr;
	Vector3<double>* h1 = nullptr;
	Vector3<double>* h2 = nullptr;
	double* neg = nullptr;
	double** vsave = nullptr;
	double** vgg = nullptr;
	
	// for spin unpolarized case, 
	// calculate the gradient of (rho_core+rho) in reciprocal space.
	rhotmp1 = new double[pw.nrxx];
	rhogsum1 = new complex<double>[pw.ngmc];
	ZEROS(rhotmp1, pw.nrxx);
	ZEROS(rhogsum1, pw.ngmc);
	for(int ir=0; ir<pw.nrxx; ir++) rhotmp1[ir] = CHR.rho[0][ir] + fac * CHR.rho_core[ir];
	for(int ig=0; ig<pw.ngmc; ig++) rhogsum1[ig] = CHR.rhog[0][ig] + fac * CHR.rhog_core[ig];

	gdr1 = new Vector3<double>[pw.nrxx];
	h1 = new Vector3<double>[pw.nrxx];
	
	GGA_PW::grad_rho( rhogsum1 , gdr1 );

	// for spin polarized case;
	// calculate the gradient of (rho_core+rho) in reciprocal space.
	if(NSPIN==2)
	{
		rhotmp2 = new double[pw.nrxx];
		rhogsum2 = new complex<double>[pw.ngmc];
		ZEROS(rhotmp2, pw.nrxx);
		ZEROS(rhogsum2, pw.ngmc);
		for(int ir=0; ir<pw.nrxx; ir++) rhotmp2[ir] = CHR.rho[1][ir] + fac * CHR.rho_core[ir];
		for(int ig=0; ig<pw.ngmc; ig++) rhogsum2[ig] = CHR.rhog[1][ig] + fac * CHR.rhog_core[ig];

		gdr2 = new Vector3<double>[pw.nrxx];
		h2 = new Vector3<double>[pw.nrxx];
		
		GGA_PW::grad_rho( rhogsum2 , gdr2 );
	}

	if(NSPIN == 4&&(DOMAG||DOMAG_Z))
	{
		ZEROS(rhotmp1, pw.nrxx);
		ZEROS(rhogsum1, pw.ngmc);
		rhotmp2 = new double[pw.nrxx];
		ZEROS(rhotmp2, pw.nrxx);
 		neg = new double [pw.nrxx];
		ZEROS(neg, pw.nrxx);
		vsave = new double* [NSPIN];
		for(int is = 0;is<NSPIN;is++) {
			vsave[is]= new double [pw.nrxx];
			for(int ir =0;ir<pw.nrxx;ir++){
				vsave[is][ir] = v(is,ir);
				v(is,ir) = 0;
			}
		}
		vgg = new double* [nspin0];
		for(int is = 0;is<nspin0;is++)vgg[is] = new double[pw.nrxx];

		noncolin_rho(rhotmp1,rhotmp2,neg);

		rhogsum2 = new complex<double>[pw.ngmc];
		ZEROS(rhogsum2, pw.ngmc);
		CHR.set_rhog(rhotmp1, rhogsum1);
		CHR.set_rhog(rhotmp2, rhogsum2);
		for(int ir=0; ir<pw.nrxx; ir++)
		{
			rhotmp2[ir] += fac * CHR.rho_core[ir];
			rhotmp1[ir] += fac * CHR.rho_core[ir];
		}
		for(int ig=0; ig<pw.ngmc; ig++)
		{
			rhogsum2[ig] += fac * CHR.rhog_core[ig];
			rhogsum1[ig] += fac * CHR.rhog_core[ig];
		}

		gdr2 = new Vector3<double>[pw.nrxx];
		h2 = new Vector3<double>[pw.nrxx];

		GGA_PW::grad_rho( rhogsum1 , gdr1 );
		GGA_PW::grad_rho( rhogsum2 , gdr2 );

	}

	// for test
	/*
	double sum[6]={0,0,0,0,0,0};
	for(int ir=0; ir<pw.nrxx; ir++)
	{
		sum[0] += abs(gdr1[0][ir]);
		sum[1] += abs(gdr1[1][ir]);
		sum[2] += abs(gdr1[2][ir]);	
		sum[3] += abs(rhotmp1[ir]);	
		sum[4] += rhotmp1[ir]*rhotmp1[ir];	
	}
	*/
	
	/*
	cout << "\n sum grad 1= " << sum[0] << " "  << sum[1] << " " << sum[2] << endl;
	cout << " sum rho = " << sum[3] << " "  << sum[4] << endl;
	ZEROS(sum,6);
	for(int ir=0; ir<pw.nrxx; ir++)
	{
		sum[0] += abs(gdr2[0][ir]);
		sum[1] += abs(gdr2[1][ir]);
		sum[2] += abs(gdr2[2][ir]);	
		sum[3] += abs(rhotmp2[ir]);	
		sum[4] += rhotmp2[ir]*rhotmp2[ir];	
	}
	cout << "\n sum grad 2= " << sum[0] << " "  << sum[1] << " " << sum[2] << endl;
	cout << " sum rho = " << sum[3] << " "  << sum[4] << endl;
	*/
	
	const double epsr = 1.0e-6;
	const double epsg = 1.0e-10;

	double grho2a = 0.0;
	double grho2b = 0.0;
	double sx = 0.0;
	double sc = 0.0;
	double v1x = 0.0;
	double v2x = 0.0;
	double v1c = 0.0;
	double v2c = 0.0;
	double vtxcgc = 0.0;
	double etxcgc = 0.0;

	if(nspin0==1)
	{
		double segno;
		for(int ir=0; ir<pw.nrxx; ir++)
		{
			const double arho = std::abs( rhotmp1[ir] );
			h1[ir].x = h1[ir].y = h1[ir].z = 0.0;
			if(arho > epsr)
			{
				grho2a = gdr1[ir].norm2();
				if( grho2a > epsg )
				{
					if( rhotmp1[ir] >= 0.0 ) segno = 1.0;
					if( rhotmp1[ir] < 0.0 ) segno = -1.0;
					
					XC_Functional::gcxc( arho, grho2a, sx, sc, v1x, v2x, v1c, v2c);

					// first term of the gradient correction:
					// D(rho*Exc)/D(rho)
					v(0, ir) += e2 * ( v1x + v1c );
					
					// h contains
					// D(rho*Exc) / D(|grad rho|) * (grad rho) / |grad rho|
					h1[ir] = e2 * ( v2x + v2c ) * gdr1[ir];
					
					vtxcgc += e2*( v1x + v1c ) * ( rhotmp1[ir] - CHR.rho_core[ir] );
					etxcgc += e2*( sx + sc ) * segno;
				}
			} // end arho > epsr
		}
	}// end nspin0 == 1
	else // spin polarized case
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
		for(int ir=0; ir<pw.nrxx; ir++)
		{
			double rh = rhotmp1[ir] + rhotmp2[ir]; 
			grho2a = gdr1[ir].norm2();
			grho2b = gdr2[ir].norm2();
			//XC_Functional::gcx_spin();
			gcx_spin(rhotmp1[ir], rhotmp2[ir], grho2a, grho2b,
				sx, v1xup, v1xdw, v2xup, v2xdw);
			
			if(rh > epsr)
			{
				if(igcc_is_lyp)
				{
					WARNING_QUIT("gga_pw","igcc_is_lyp is not available now.");
				}
				else
				{
					double zeta = ( rhotmp1[ir] - rhotmp2[ir] ) / rh;
					if(NSPIN==4&&(DOMAG||DOMAG_Z)) zeta = fabs(zeta) * neg[ir];
					const double grh2 = (gdr1[ir]+gdr2[ir]).norm2();
					//XC_Functional::gcc_spin(rh, zeta, grh2, sc, v1cup, v1cdw, v2c);
					gcc_spin(rh, zeta, grh2, sc, v1cup, v1cdw, v2c);
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

			// first term of the gradient correction : D(rho*Exc)/D(rho)
			v(0,ir) = v(0,ir) + e2 * ( v1xup + v1cup );
			v(1,ir) = v(1,ir) + e2 * ( v1xdw + v1cdw );

//			continue; //mohan tmp
			
			// h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
			h1[ir] = e2 * ( ( v2xup + v2cup ) * gdr1[ir] + v2cud * gdr2[ir] );
			h2[ir] = e2 * ( ( v2xdw + v2cdw ) * gdr2[ir] + v2cud * gdr1[ir] );

			vtxcgc = vtxcgc + e2 * ( v1xup + v1cup ) * ( rhotmp1[ir] - CHR.rho_core[ir] * fac );
			vtxcgc = vtxcgc + e2 * ( v1xdw + v1cdw ) * ( rhotmp2[ir] - CHR.rho_core[ir] * fac );
			etxcgc = etxcgc + e2 * ( sx + sc );
			

		}// end ir

	}

	//cout << "\n vtxcgc=" << vtxcgc;
	//cout << "\n etxcgc=" << etxcgc << endl;

	for(int ir=0; ir<pw.nrxx; ir++) rhotmp1[ir] -= fac * CHR.rho_core[ir];
	if(nspin0==2) for(int ir=0; ir<pw.nrxx; ir++) rhotmp2[ir] -= fac * CHR.rho_core[ir];
	
	// second term of the gradient correction :
	// \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )

	// dh is in real sapce.
	double* dh = new double[pw.nrxx];

	for(int is=0; is<nspin0; is++)
	{
		ZEROS(dh, pw.nrxx);
		if(is==0)GGA_PW::grad_dot(h1,dh);
		if(is==1)GGA_PW::grad_dot(h2,dh);

		for(int ir=0; ir<pw.nrxx; ir++)
			v(is, ir) -= dh[ir];
		
		double sum = 0.0;
		if(is==0)
			for(int ir=0; ir<pw.nrxx; ir++)
				sum += dh[ir] * rhotmp1[ir];
		else if(is==1)
			for(int ir=0; ir<pw.nrxx; ir++)
				sum += dh[ir] * rhotmp2[ir];
		
		vtxcgc -= sum;
	}

	delete[] dh;

	vtxc += vtxcgc;
	etxc += etxcgc;

	if(NSPIN == 4 && (DOMAG||DOMAG_Z))
	{
		for(int is=0;is<NSPIN;is++)
		{
			for(int ir=0;ir<pw.nrxx;ir++)
			{
				if(is<nspin0) vgg[is][ir] = v(is,ir);
				v(is,ir) = vsave[is][ir];
			}
		}
		for(int ir=0;ir<pw.nrxx;ir++)
		{
			v(0,ir) += 0.5 * (vgg[0][ir] + vgg[1][ir]);
			double amag = sqrt(pow(CHR.rho[1][ir],2)+pow(CHR.rho[2][ir],2)+pow(CHR.rho[3][ir],2));
			if(amag>1e-12)
			{
				for(int i=1;i<4;i++)
					v(i,ir)+= neg[ir] * 0.5 *(vgg[0][ir]-vgg[1][ir])*CHR.rho[i][ir]/amag;
			}
		}
	}
	
	// deacllocate
	delete[] rhotmp1;
	delete[] rhogsum1;
	delete[] gdr1;
	delete[] h1;

	if(NSPIN==2)
	{
		delete[] rhotmp2;
		delete[] rhogsum2;
		delete[] gdr2;
		delete[] h2;
	}
	if(NSPIN == 4 && (DOMAG||DOMAG_Z))
	{
		delete[] neg;
		for(int i=0; i<nspin0; i++) delete[] vgg[i];
		delete[] vgg;
		for(int i=0; i<NSPIN; i++) delete[] vsave[i];
		delete[] vsave;
		delete[] rhotmp2;
		delete[] rhogsum2;
		delete[] gdr2;
		delete[] h2;
	}

	return;
}

void GGA_PW::grad_rho( const complex<double> *rhog, Vector3<double> *gdr )
{
	complex<double> *gdrtmpg = new complex<double>[pw.ngmc];
	ZEROS(gdrtmpg, pw.ngmc);

	complex<double> *Porter = UFFT.porter;

	// the formula is : rho(r)^prime = \int iG * rho(G)e^{iGr} dG
	for(int ig=0; ig<pw.ngmc; ig++)
		gdrtmpg[ig] = IMAG_UNIT * rhog[ig];

	// calculate the charge density gradient in reciprocal space.
	ZEROS(Porter, pw.nrxx);
	for(int ig=0; ig<pw.ngmc; ig++)
		Porter[ pw.ig2fftc[ig] ] = gdrtmpg[ig]* complex<double>(pw.get_G_cartesian_projection(ig, 0), 0.0);
	// bring the gdr from G --> R
	pw.FFT_chg.FFT3D(Porter, 1);
	// remember to multily 2pi/a0, which belongs to G vectors.
	for(int ir=0; ir<pw.nrxx; ir++)
		gdr[ir].x = Porter[ir].real() * ucell.tpiba;

	// calculate the charge density gradient in reciprocal space.
	ZEROS(Porter, pw.nrxx);
	for(int ig=0; ig<pw.ngmc; ig++)
		Porter[pw.ig2fftc[ig]] = gdrtmpg[ig] * complex<double>(pw.get_G_cartesian_projection(ig, 1), 0.0);
	// bring the gdr from G --> R
	pw.FFT_chg.FFT3D(Porter, 1);
	// remember to multily 2pi/a0, which belongs to G vectors.
	for(int ir=0; ir<pw.nrxx; ir++)
		gdr[ir].y = Porter[ir].real() * ucell.tpiba;

	// calculate the charge density gradient in reciprocal space.
	ZEROS(Porter, pw.nrxx);
	for(int ig=0; ig<pw.ngmc; ig++)
		Porter[pw.ig2fftc[ig]] = gdrtmpg[ig] * complex<double>(pw.get_G_cartesian_projection(ig, 2), 0.0);
	// bring the gdr from G --> R
	pw.FFT_chg.FFT3D(Porter, 1);
	// remember to multily 2pi/a0, which belongs to G vectors.
	for(int ir=0; ir<pw.nrxx; ir++)
		gdr[ir].z = Porter[ir].real() * ucell.tpiba;

	delete[] gdrtmpg;
	return;
}


void GGA_PW::grad_dot(const Vector3<double> *h, double *dh)
{
	complex<double> *aux = new complex<double>[pw.nrxx];
	complex<double> *gaux = new complex<double>[pw.ngmc];
	ZEROS(gaux, pw.ngmc);
	
	ZEROS(aux, pw.nrxx);
	for(int ir=0; ir<pw.nrxx; ir++)
		aux[ir] = complex<double>( h[ir].x, 0.0);
	// bring to G space.
	pw.FFT_chg.FFT3D(aux, -1);
	for(int ig=0; ig<pw.ngmc; ig++)
		gaux[ig] += pw.get_G_cartesian_projection(ig, 0) * IMAG_UNIT * aux[pw.ig2fftc[ig]];

	ZEROS(aux, pw.nrxx);
	for(int ir=0; ir<pw.nrxx; ir++)
		aux[ir] = complex<double>( h[ir].y, 0.0);
	// bring to G space.
	pw.FFT_chg.FFT3D(aux, -1);
	for(int ig=0; ig<pw.ngmc; ig++)
		gaux[ig] += pw.get_G_cartesian_projection(ig, 1) * IMAG_UNIT * aux[pw.ig2fftc[ig]];

	ZEROS(aux, pw.nrxx);
	for(int ir=0; ir<pw.nrxx; ir++)
		aux[ir] = complex<double>( h[ir].z, 0.0);
	// bring to G space.
	pw.FFT_chg.FFT3D(aux, -1);
	for(int ig=0; ig<pw.ngmc; ig++)
		gaux[ig] += pw.get_G_cartesian_projection(ig, 2) * IMAG_UNIT * aux[pw.ig2fftc[ig]];

	ZEROS(aux, pw.nrxx);
	for(int ig=0; ig<pw.ngmc; ig++)
		aux[ pw.ig2fftc[ig] ] = gaux[ig];
	// bring back to R space
	pw.FFT_chg.FFT3D(aux, 1);
	for(int ir=0; ir<pw.nrxx; ir++)
		dh[ir] = aux[ir].real() * ucell.tpiba;
	
	delete[] aux;	
	delete[] gaux; //mohan fix 2012-04-02
	return;
}

void GGA_PW::noncolin_rho(double *rhoout1,double *rhoout2, double *neg)
{
	//this function diagonalizes the spin density matrix and gives as output the
	//spin up and spin down components of the charge.
	//If lsign is true up and dw are with respect to the fixed quantization axis 
	//ux, otherwise rho + |m| is always rhoup and rho-|m| is always rhodw.
	double amag=0;
	for(int ir = 0;ir<pw.nrxx;ir++) neg[ir] = 1.0;
	if(mag.lsign_)
	{
		for(int ir = 0;ir<pw.nrxx;ir++)
		{
			if(CHR.rho[1][ir]*mag.ux_[0] + CHR.rho[2][ir]*mag.ux_[1] + CHR.rho[3][ir]*mag.ux_[2]>0) neg[ir] = 1.0;
			else neg[ir] = -1.0;
		}
	}
	for(int ir = 0;ir<pw.nrxx;ir++)
	{
		amag = sqrt(pow(CHR.rho[1][ir],2)+pow(CHR.rho[2][ir],2)+pow(CHR.rho[3][ir],2));
		rhoout1[ir] = 0.5 * (CHR.rho[0][ir] + neg[ir] * amag);
		rhoout2[ir] = 0.5 * (CHR.rho[0][ir] - neg[ir] * amag);
	}
	return;
}

