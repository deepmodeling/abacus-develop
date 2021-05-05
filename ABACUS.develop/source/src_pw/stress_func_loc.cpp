#include"stress_func.h"

//calculate local pseudopotential stress in PW or VL_dVL stress in LCAO
void Stress_Func::stress_loc(matrix& sigma, const bool is_pw)
{
    timer::tick("Stress_Func","stress_loc",'F');

    double *dvloc;
    double evloc,fact=1.0;
    int ng,nt,l,m,is;

	if (INPUT.gamma_only && is_pw) fact=2.0;

    dvloc = new double[pw.ngmc];

	complex<double> *Porter = UFFT.porter;

	ZEROS( Porter, pw.nrxx );
	for(int is=0; is<NSPIN; is++)
	{
		for (int ir=0; ir<pw.nrxx; ir++)
		{
			Porter[ir] += complex<double>(CHR.rho[is][ir], 0.0 );
		}
	}
	pw.FFT_chg.FFT3D(Porter, -1);

//    if(INPUT.gamma_only==1) fact=2.0;
//    else fact=1.0;

	evloc=0.0;
	double g[3]={0,0,0};


	complex<double> *vg = new complex<double>[pw.ngmc];
	ZEROS( vg, pw.ngmc );
	for (int it=0; it<ucell.ntype; it++)
	{
		if (pw.gstart==1) evloc += ppcell.vloc(it, pw.ig2ngg[0]) * (pw.strucFac(it,0) * conj(Porter[pw.ig2fftc[0]])).real();
		for (int ig=pw.gstart; ig<pw.ngmc; ig++)
		{
			const int j = pw.ig2fftc[ig];
			evloc += ppcell.vloc(it, pw.ig2ngg[ig]) * (pw.strucFac(it,ig) * conj(Porter[j]) * fact).real();
		}
	}
	for( nt = 0;nt< ucell.ntype; nt++)
	{
		const Atom* atom = &ucell.atoms[nt];
		//mark by zhengdy for check
		// if ( ppcell.vloc == NULL ){
		if(0)
		{
		//
		// special case: pseudopotential is coulomb 1/r potential
		//
			this->dvloc_coul (atom->zv, dvloc);
		//
		}
		else
		{
		//
		// normal case: dvloc contains dV_loc(G)/dG
		//
			this->dvloc_of_g ( atom->msh, atom->rab, atom->r,
					atom->vloc_at, atom->zv, dvloc);
		//
		}

		for( ng = 0;ng< pw.ngmc;ng++)
		{
			const int j = pw.ig2fftc[ng];
			g[0]=pw.gcar[ng].x;
			g[1]=pw.gcar[ng].y;
			g[2]=pw.gcar[ng].z;
			for (l = 0;l< 3;l++)
			{
				for (m = 0; m<l+1;m++)
				{
					sigma(l,m) = sigma(l,m) +  ( conj( Porter[j] )
							* pw.strucFac (nt, ng) ).real() * 2.0 * dvloc [ pw.ig2ngg[ng] ]
							* ucell.tpiba2 * g[l] * g[m] * fact;
				}
			}
		}
	}

    for( l = 0;l< 3;l++)
	{
		if(is_pw) sigma(l,l) += evloc;
		for (m = 0; m<l+1; m++)
		{
			sigma(m,l) = sigma(l,m);
		}
	}
	for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			Parallel_Reduce::reduce_double_pool( sigma(l,m) );
		}
	}
	delete[] dvloc;
	delete[] vg;


	timer::tick("Stress_Func","stress_loc");
	return;
}

void Stress_Func::dvloc_of_g 
(
const int& msh,
const double* rab,
const double* r,
const double* vloc_at,
const double& zp,
double*  dvloc
)
{
  //----------------------------------------------------------------------
  //
  // dvloc = D Vloc (g^2) / D g^2 = (1/2g) * D Vloc(g) / D g
  //

  //
  //double  dvloc[ngl];
  // the fourier transform dVloc/dG
  //
	double vlcp=0;
	double  *aux, *aux1;

	int i, igl, igl0;
	// counter on erf functions or gaussians
	// counter on g shells vectors
	// first shell with g != 0

	aux = new double[msh];
	aux1 = new double[msh];
	ZEROS(aux, msh);
	ZEROS(aux1, msh);

	// the  G=0 component is not computed
	if (pw.ggs[0] < 1.0e-8)
	{
		dvloc[0] = 0.0;
		igl0 = 1;
	}
	else
	{
		igl0 = 0;
	}

	// Pseudopotentials in numerical form (Vloc contains the local part)
	// In order to perform the Fourier transform, a term erf(r)/r is
	// subtracted in real space and added again in G space

	//
	//   This is the part of the integrand function
	//   indipendent of |G| in real space
	//
	for( i = 0;i< msh; i++)
	{
		aux1[i] = r [i] * vloc_at [i] + zp * e2 * erf(r[i]);
	}
	for( igl = igl0;igl< pw.nggm;igl++)
	{
		double gx = sqrt (pw.ggs [igl] * ucell.tpiba2);
		double gx2 = pw.ggs [igl] * ucell.tpiba2;
		//
		//    and here we perform the integral, after multiplying for the |G|
		//    dependent  part
		//
		// DV(g)/Dg = Integral of r (Dj_0(gr)/Dg) V(r) dr
		for( i = 1;i< msh;i++)
		{
			aux [i] = aux1 [i] * (r [i] * cos (gx * r [i] ) / gx - sin (gx * r [i] ) / pow(gx,2));
		}
		// simpson (msh, aux, rab, vlcp);
		Mathzone::Simpson_Integral(msh, aux, rab, vlcp );
		// DV(g^2)/Dg^2 = (DV(g)/Dg)/2g
		vlcp *= FOUR_PI / ucell.omega / 2.0 / gx;
		// subtract the long-range term
		double g2a = gx2 / 4.0;
		vlcp += FOUR_PI / ucell.omega * zp * e2 * exp ( - g2a) * (g2a + 1) / pow(gx2 , 2);
		dvloc [igl] = vlcp;
	}
	delete[] aux1;
	delete[] aux;

	return;
}

void Stress_Func::dvloc_coul 
(
const double& zp,
double* dvloc
)
{
	//----------------------------------------------------------------------
	//
	//    Fourier transform of the Coulomb potential - For all-electron
	//    calculations, in specific cases only, for testing purposes
	//


	// fourier transform: dvloc = D Vloc (g^2) / D g^2 = 4pi e^2/omegai /G^4

	int  igl0;
	// first shell with g != 0

	// the  G=0 component is 0
	if (pw.ggs[0] < 1.0e-8)
	{
		dvloc[0] = 0.0;
		igl0 = 1;
	}
	else
	{
		igl0 = 0;
	}
	for(int i=igl0;i<pw.nggm;i++)
	{
		dvloc[i] = FOUR_PI * zp * e2 / ucell.omega / pow(( ucell.tpiba2 * pw.ggs[i] ),2);
	}
	
	return;
}
