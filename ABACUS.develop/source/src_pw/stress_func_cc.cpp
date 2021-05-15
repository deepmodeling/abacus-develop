#include "./stress_func.h"
#include "./H_XC_pw.h"
#include "../src_global/math_integral.h"

//NLCC term, need to be tested
void Stress_Func::stress_cc(matrix& sigma, const bool is_pw)
{
	timer::tick("Stress_Func","stress_cc",'F');
        
	int nt,ng,l,m,ir;
	double fact=1.0;

	if(is_pw&&INPUT.gamma_only) 
	{
		fact = 2.0; //is_pw:PW basis, gamma_only need to double.
	}

	complex<double> sigmadiag;
	double* rhocg;
	double g[3];

	int judge=0;
	for(nt=0;nt<ucell.ntype;nt++)
	{
		if(ucell.atoms[nt].nlcc) 
		{
			judge++;
		}
	}

	if(judge==0) 
	{
		return;
	}

	//recalculate the exchange-correlation potential
    const auto etxc_vtxc_v = H_XC_pw::v_xc(pw.nrxx, pw.ncxyz, ucell.omega, CHR.rho, CHR.rho_core);
	H_XC_pw::etxc    = std::get<0>(etxc_vtxc_v);			// may delete?
	H_XC_pw::vtxc    = std::get<1>(etxc_vtxc_v);			// may delete?
	const matrix vxc = std::get<2>(etxc_vtxc_v);

	complex<double> * psic = new complex<double> [pw.nrxx];

	ZEROS(psic, pw.nrxx);

	if(NSPIN==1||NSPIN==4)
	{
		for(ir=0;ir<pw.nrxx;ir++)
		{
			// psic[ir] = vxc(0,ir);
			psic[ir] = complex<double>(vxc(0, ir),  0.0);
		}
	}
	else
	{
		for(ir=0;ir<pw.nrxx;ir++)
		{
			psic[ir] = 0.5 * (vxc(0, ir) + vxc(1, ir));
		}
	}

	// to G space
	pw.FFT_chg.FFT3D(psic, -1);

	//psic cantains now Vxc(G)
	rhocg= new double [pw.nggm];
	ZEROS(rhocg, pw.nggm);

	sigmadiag=0.0;
	for(nt=0;nt<ucell.ntype;nt++)
	{
		if(ucell.atoms[nt].nlcc)
		{
			//drhoc();
			CHR.non_linear_core_correction(
				ppcell.numeric,
				ucell.atoms[nt].msh,
				ucell.atoms[nt].r,
				ucell.atoms[nt].rab,
				ucell.atoms[nt].rho_atc,
				rhocg);


			//diagonal term 
			if (pw.gstart==1) 
			{
				sigmadiag += conj(psic [pw.ig2fftc[0]] ) * pw.strucFac (nt, 0) * rhocg [pw.ig2ngg[0] ];
			}
			for( ng = pw.gstart;ng< pw.ngmc;ng++)
			{
				sigmadiag +=  conj(psic[pw.ig2fftc[ng]] ) *
					pw.strucFac (nt, ng) * rhocg [pw.ig2ngg[ng] ] * fact;
			}
			this->deriv_drhoc (
				ppcell.numeric,
				ucell.atoms[nt].msh,
				ucell.atoms[nt].r,
				ucell.atoms[nt].rab,
				ucell.atoms[nt].rho_atc,
				rhocg);
			// non diagonal term (g=0 contribution missing)
			for( ng = pw.gstart;ng< pw.ngmc;ng++)
			{
				g[0] = pw.gcar[ng].x;
				g[1] = pw.gcar[ng].y;
				g[2] = pw.gcar[ng].z;
				for( l = 0;l< 3;l++)
				{
					for (m = 0;m< 3;m++)
					{
						const complex<double> t = conj(psic[pw.ig2fftc[ng]] )
							* pw.strucFac (nt, ng) * rhocg [pw.ig2ngg[ng] ] * ucell.tpiba *
							g [l] * g [m] / pw.gcar[ng].norm() * fact;
//						sigmacc [l][ m] += t.real();
						sigma(l,m) += t.real();
					}//end m
				}//end l
			}//end ng
		}//end if
	}//end nt

	for( l = 0;l< 3;l++)
	{
		sigma(l,l) += sigmadiag.real();
//		sigmacc [l][ l] += sigmadiag.real();
	}
	for( l = 0;l< 3;l++)
	{
		for (m = 0;m< 3;m++)
		{
			Parallel_Reduce::reduce_double_pool( sigma(l,m) );
		}
	}

	delete[] rhocg;
	delete[] psic;

	timer::tick("Stress_Func","stress_cc");
	return;
}


void Stress_Func::deriv_drhoc 
(
	const bool &numeric,
	const int mesh,
	const double *r,
	const double *rab,
	const double *rhoc,
	double *drhocg
)
{

	double gx = 0, rhocg1 = 0;
	// the modulus of g for a given shell
	// the fourier transform
	double *aux = new double[ mesh];
	// auxiliary memory for integration

	int  igl0;
	// counter on radial mesh points
	// counter on g shells
	// lower limit for loop on ngl

	//
	// G=0 term
	//
	if (pw.ggs[0] < 1.0e-8)
	{
		drhocg [0] = 0.0;
		igl0 = 1;
	}
	else
	{
		igl0 = 0;
	}
	//
	// G <> 0 term
	//
	
	for(int igl = igl0;igl< pw.nggm;igl++)
	{
		gx = sqrt(pw.ggs [igl] * ucell.tpiba2);
		for( int ir = 0;ir< mesh; ir++)
		{
			aux [ir] = r [ir] * rhoc [ir] * (r [ir] * cos (gx * r [ir] ) / gx - sin (gx * r [ir] ) / pow(gx,2));
		}//ir
		Integral::Simpson_Integral(mesh, aux, rab, rhocg1);
		drhocg [igl] = FOUR_PI / ucell.omega * rhocg1;
	}//igl
	
	delete [] aux;

	return;
}
