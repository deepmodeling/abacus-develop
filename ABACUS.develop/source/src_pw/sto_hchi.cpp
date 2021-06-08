#include "tools.h"
#include "global.h"
#include "sto_hchi.h" 

int Stochastic_hchi::nrxx;
int Stochastic_hchi::nx;
int Stochastic_hchi::ny;
int Stochastic_hchi::nz;

fftw_plan Stochastic_hchi::pb;
fftw_plan Stochastic_hchi::pf;

double Stochastic_hchi::Emin;
double Stochastic_hchi::Emax;

bool Stochastic_hchi::initplan;
bool Stochastic_hchi::ortho;

complex<double>* Stochastic_hchi::rp_chi;
complex<double>* Stochastic_hchi::rl_chi;

int * Stochastic_hchi:: GRA_index;

Stochastic_hchi::Stochastic_hchi()
{
	initplan = false;
	ortho = false;
	nrxx = 0;
	rp_chi = new complex<double> [1];
	rl_chi = new complex<double> [1];
	GRA_index = new int [1];
}

Stochastic_hchi::~Stochastic_hchi()
{
	if(initplan)
	{
		fftw_destroy_plan(pb);
		fftw_destroy_plan(pf);
		delete[] rp_chi;
		delete[] rl_chi;
		delete[] GRA_index;
	}
}

void Stochastic_hchi:: init()
{
	//wait for init--------------------------------------
		//nrxx
	//---------------------------------------------------
	nrxx = pw.nrxx;
	nx = pw.nx;
	ny = pw.ny;
	nz = pw.nz;
    if(nrxx != 0)
    {
        delete[] rp_chi;
        delete[] rl_chi;
		delete[] GRA_index;
		GRA_index = new int [wf.npw];
		if(STO_WF.stotype != "pw")
		{
			rp_chi = new complex<double> [nrxx];
			rl_chi = new complex<double> [nrxx];

			pf=fftw_plan_dft_3d(nx,ny,nz,
				(fftw_complex *)rp_chi,
				(fftw_complex *)rp_chi, 
				FFTW_FORWARD, 
				FFTW_MEASURE);

			pb=fftw_plan_dft_3d(nx,ny,nz,
				(fftw_complex *)rl_chi,
				(fftw_complex *)rl_chi, 
				FFTW_BACKWARD, 
				FFTW_MEASURE);

			initplan = true;
		}
    }
    else
    {
        WARNING_QUIT("Stochastic_hchi", "Number of grids should be at least one!");
    }

}

void Stochastic_hchi::get_GRA_index()
{
	
	if(STO_WF.stotype == "pw")
	{
		for(int ig = 0 ; ig < wf.npw; ++ig)
		{
			GRA_index[ig] = pw.ig2fftw[wf.igk(0,ig)]; //GAMMA POINT temporarily
		}
	}
	else
	{
		int ix,iy,iz;
		int ir;
		ZEROS(GRA_index,wf.npw);
		for(int ig = 0 ; ig < wf.npw; ++ig)
		{
			ix = floor(pw.gcar[wf.igk(0, ig)].x+0.1);
			iy = floor(pw.gcar[wf.igk(0, ig)].y+0.1);
			iz = floor(pw.gcar[wf.igk(0, ig)].z+0.1);
			if(ix < 0) ix += nx;
			if(iy < 0) iy += ny;
			if(iz < 0) iz += nz;
			ir = ix * ny * nz + iy * nz + iz;
			GRA_index[ig] = ir;
		}
	}
}

void Stochastic_hchi::orthogonal_to_psi_real(complex<double> *wfin, complex<double> *wfout, int &ikk)
{

	TITLE("Stochastic_hchi","orthogonal_to_psi0");
	if(!initplan) WARNING_QUIT("Stochastic_hchi", "Please init hchi first!");

	DCOPY(wfin,rp_chi,nrxx);
	//LapackConnector::copy(nrxx,wfin,1,rp_chi,1);
	fftw_execute(pf);
	
	complex<double> * chig = new complex<double> [wf.npw];
	for(int ig = 0; ig < wf.npw; ++ig)
	{
		chig[ig] = rp_chi[GRA_index[ig]]; 
	}

	//orthogonal part
	complex<double> sum = 0;
	int inc=1;
	for(int iksb = 0; iksb < NBANDS; ++iksb)
	{
		complex<double> *kswf = &wf.evc[ikk](iksb,0); 
		zdotc_(&sum,&wf.npw,kswf,&inc,chig,&inc);

//
		//LapackConnector::axpy(wf.npw,-sum,kswf,1,chig,1);
		for(int ig = 0; ig < wf.npw; ++ig)
		{
			chig[ig] -= sum*kswf[ig];
		}
	}

	//test orthogonal in reciprocal space
	//complex<double> overlap;
	//for(int iksb = 0; iksb < NBANDS; ++iksb)
	//{
	//	complex<double> *kswf = &wf.evc[ikk](iksb,0); 
	//	overlap=0;
	//	for(int ig = 0; ig < wf.npw; ++ig)
	//	{
	//		overlap += conj(kswf[ig]) * chig[ig];
	//	}
	//	cout<<"OVERLAP "<<overlap<<endl;
	//}

	for(int ig = 0; ig < wf.npw; ++ig)
	{
		rp_chi[GRA_index[ig]] = chig[ig];
	}
	DCOPY(rp_chi,rl_chi,nrxx);
	//LapackConnector::copy(nrxx,rp_chi,1,rl_chi,1);

	fftw_execute(pb);

	//LapackConnector::copy(nrxx,rl_chi,1,wfout,1);
	//LapackConnector::scal(nrxx,1/double(nrxx),wfout,1);
	for(int ir = 0; ir < nrxx ; ++ir)
	{
		wfout[ir] = rl_chi[ir] / nrxx;
	}
	

	/*//test orthogonal in real space
	complex<double> overlap;
	complex<double> * kswf = new complex<double> [nrxx];
	fftw_plan pp=fftw_plan_dft_3d(pw.nx,pw.ny,pw.nz,(fftw_complex *)kswf,(fftw_complex *)kswf, FFTW_BACKWARD, FFTW_ESTIMATE);
	for(int iksb = 0; iksb < NBANDS; ++iksb)
	{
		ZEROS(kswf,nrxx);
		for(int ig = 0 ; ig < wf.npw; ++ig)
		{
			kswf[GRA_index[ig]] = wf.evc[ikk](iksb,ig);
		}
		fftw_execute(pp);
		overlap=0;
		for(int ir = 0; ir < nrxx; ++ir)
		{
			overlap += conj(kswf[ir]) * wfout[ir];
		}
		cout<<"OVERLAP "<<overlap<<endl;
	}
	delete [] kswf;*/

	ortho = true;
	delete[] chig;
	return;
}

void Stochastic_hchi::orthogonal_to_psi_reciprocal(complex<double> *wfgin, complex<double> *wfgout, int &ikk)
{

	TITLE("Stochastic_hchi","orthogonal_to_psi0");
	int nchip=STO_WF.nchip;
	int npw = wf.npw;
	for(int ig = 0 ; ig < npw * nchip; ++ig)
	{
		wfgout[ig] = wfgin[ig];
	}

	//orthogonal part
	
	complex<double> *sum = new complex<double> [NBANDS * nchip];
	char transC='C';
	char transN='N';
	
	//sum(b<NBANDS, a<nchi) = < psi_b | chi_a >
	zgemm_(&transC, &transN, &NBANDS, &nchip, &npw, &ONE, wf.evc[ikk].c, &wf.npwx, wfgout, &npw, &ZERO, sum, &NBANDS);
	Parallel_Reduce::reduce_complex_double_pool(sum, NBANDS * nchip);
	
	//psi -= psi * sum
	zgemm_(&transN, &transN, &npw, &nchip, &NBANDS, &NEG_ONE, wf.evc[ikk].c, &wf.npwx, sum, &NBANDS, &ONE, wfgout, &npw);
	
	ortho = true;
	delete[] sum;
	return;
}



void Stochastic_hchi::hchi_real(complex<double>*chi_in, complex<double> *hchi, const int m)
{
	
	double*vr = pot.vr_eff1;  //vr= pot.vrs1 temporarily use cutoff vr.

	//wait for init--------------------------------------	
	double dk1,dk2,dk3;
	dk1 = ucell.tpiba;
	dk2 = ucell.tpiba;
	dk3 = ucell.tpiba;
		
	//---------------------------------------------------
	if(!initplan) WARNING_QUIT("Stochastic_hchi", "Please init hchi first!");

	ZEROS(hchi,nrxx);
	DCOPY(chi_in, rp_chi, nrxx);
	//LapackConnector::copy(nrxx,chi_in,1,rp_chi,1);
	fftw_execute(pf);

	complex<double> * chig = new complex<double> [wf.npw];
	ZEROS(chig,wf.npw);
	for(int ig = 0; ig < wf.npw; ++ig)
	{
		chig[ig] = rp_chi[GRA_index[ig]]; 
	}

	//------------------------------------
	//(1) the local potential.
	//------------------------------------

	if(VL_IN_H)
	{
		for(int ir = 0; ir < nrxx; ++ir)
		{
			hchi[ir] += chi_in[ir] * vr[ir] ;
		}
		
	}
	/*cout<<"HCHI-------------------------"<<endl;

		for(int ir  = 63800 ; ir < 64000; ++ir)
		{
			cout<<chi_in[ir] * vr[ir]<<" ";
		}
		cout<<endl;*/
	
	
	
	//------------------------------------
	//(2) the kinetical energy.
	//------------------------------------
	if(T_IN_H)
	{
		Vector3<double> gg;
		int gx,gy,gz;
		for(int ig1 = 0, i = 0; ig1 < nx; ++ig1)
		{
			for(int ig2 = 0; ig2 < ny; ++ig2)
			{
				for(int ig3 = 0; ig3 < nz; ++ig3)
				{
					gx = ig1;
					gy = ig2;
					gz = ig3;
					if(ig1 > nx/2) gx -= nx;
					if(ig2 > ny/2) gy -= ny;
					if(ig3 > nz/2) gz -= nz;
					gg.set(gx*dk1, gy*dk2, gz*dk3);
					rl_chi[i] = gg.norm2()*rp_chi[i];
					++i;
				}
			}
		
		}
	}

	
	

	//------------------------------------
	// (3) the nonlocal pseudopotential.
	//------------------------------------
	int inc = 1;
	if(VNL_IN_H)
	{
		if ( ppcell.nkb > 0)
		{
			complex<double> *becp = new complex<double>[ ppcell.nkb * NPOL ];
			ZEROS(becp,ppcell.nkb * NPOL);

			for (int i=0;i< ppcell.nkb;++i)
			{
				const complex<double>* p = &ppcell.vkb(i,0);
				zdotc_(&becp[i],&wf.npw,p,&inc,chig,&inc);
				//for (int ig=0; ig< wf.npw; ++ig)
				//{
				//	if(NSPIN!=4) becp[i] += chig[ig] * conj( p[ig] );
				//} 
			}

			//Parallel_Reduce::reduce_complex_double_pool( becp, ppcell.nkb * NPOL);
			complex<double> * Ps = new complex<double> [ppcell.nkb * NPOL];
			ZEROS( Ps, ppcell.nkb * NPOL );
			int sum = 0;
    		int iat = 0;
    		// this function sum up each non-local pseudopotential located in each atom,
    		// all we need to do is put the right Dij coefficient to each becp, which
    		// is calculated before.
    		for (int it=0; it<ucell.ntype; ++it)
    		{
    		    const int Nprojs = ucell.atoms[it].nh;
    		    for (int ia=0; ia<ucell.atoms[it].na; ++ia)
    		    {
    		        // each atom has Nprojs, means this is with structure factor;
    		        // each projector (each atom) must multiply coefficient
    		        // with all the other projectors.
    		        for (int ip=0; ip<Nprojs; ++ip)
    		        {
    		            for (int ip2=0; ip2<Nprojs; ++ip2)
    		            {
							if(NSPIN!=4)
								Ps[sum+ip2] += ppcell.deeq(CURRENT_SPIN, iat, ip, ip2) * becp[sum+ip];
    		            }// end ih
    		        }//end jh
					if(NSPIN!=4) sum += Nprojs;
					++iat;
    		    } //end na
    		} //end nt

			// use simple method.
			ZEROS(chig, wf.npw);
			if(NSPIN!=4)
				for(int i=0; i<ppcell.nkb; ++i)
				{
					complex<double>* p = &ppcell.vkb(i,0);
					//LapackConnector::axpy(wf.npw,Ps[i],p,1,chig,1);
					for(int ig=0; ig< wf.npw; ++ig)
					{
						chig[ig] += Ps[i] * p[ig];
					}
				}
			delete[] becp;
			delete[] Ps;
		}
		for(int ig = 0; ig < wf.npw; ++ig)
		{
			rl_chi[GRA_index[ig]] += chig[ig];
		}
	}

	
	//------------------------------------
	// (4) Conver (2) & (3) in Reciprocal space to Real one
	//------------------------------------
	fftw_execute(pb);

	double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;

	//LapackConnector::axpy(nrxx,1/double(nrxx),rl_chi,1,hchi,1);
	//LapackConnector::axpy(nrxx,-Ebar,chi_in,1,hchi,1);
	//LapackConnector::scal(nrxx,1/DeltaE,hchi,1);
	for(int i = 0; i < nrxx; ++i)
	{
		hchi[i] += rl_chi[i] / nrxx;
	}
	for(int i = 0; i < nrxx; ++i)
	{
		hchi[i] = (hchi[i] - Ebar * chi_in[i]) / DeltaE;
	}

	delete [] chig;

	//test Emax & Emin
	//------------------------------------------------------------
	//double sum1 = 0;
	//double sum2 = 0;
	//for(int i = 0 ; i < nrxx; ++i)
	//{
	//	sum1 += norm(chi_in[i]);
	//	sum2 += real(conj(chi_in[i]) * hchi[i]);
	//}
	//cout<<setw(15)<<sum2 <<setw(15)<<sum1<<setw(15)<<sum2/sum1<<endl;
	//------------------------------------------------------------

	//test hermit property
	//------------------------------------------------------------
	//complex<double> sum=0;
	//for(int i = 0 ; i < nrxx; ++i)
	//{
	//	sum+=conj(chi_in[i]) * hchi[i];
	//}
	//cout<<sum<<" must be real numebr."<<endl; //sum must be a real number
	//------------------------------------------------------------
	return;
}

void Stochastic_hchi:: hchi_reciprocal(complex<double> *chig, complex<double> *hchig, const int m)
{
	timer::tick("Stochastic_hchi","hchi_reciprocal",'H');
	
	//---------------------------------------------------

	int npw = wf.npw;
	int npm = NPOL * m;
	int inc = 1;
	//------------------------------------
	//(1) the kinetical energy.
	//------------------------------------
	complex<double> *chibg = chig;
	complex<double> *hchibg = hchig;
	if(T_IN_H)
	{
		for (int ib = 0; ib < m ; ++ib)
		{
			for (int ig = 0; ig < npw; ++ig)
			{
				hchibg[ig] = wf.g2kin[ig] * chibg[ig];
			}
			chibg += npw;
			hchibg += npw;
		}
	}
	
	//------------------------------------
	//(2) the local potential.
	//------------------------------------
	timer::tick("Stochastic_hchi","vloc",'H');
	if(VL_IN_H)
	{
		chibg = chig;
		hchibg = hchig;
		for(int ib = 0 ; ib < m ; ++ib)
		{
			ZEROS( UFFT.porter, pw.nrxx);
			UFFT.RoundTrip( chibg, pot.vr_eff1, GRA_index, UFFT.porter );
			for (int ig = 0; ig < npw; ++ig)
			{
				hchibg[ig] += UFFT.porter[ GRA_index[ig] ];
			}
			chibg += npw;
			hchibg += npw;
		}
			
	}
	timer::tick("Stochastic_hchi","vloc",'H');


	//------------------------------------
	// (3) the nonlocal pseudopotential.
	//------------------------------------
	timer::tick("Stochastic_hchi","vnl",'H');
	if(VNL_IN_H)
	{
		if ( ppcell.nkb > 0)
		{
			int nkb = ppcell.nkb;
			complex<double> *becp = new complex<double>[ nkb * NPOL * m ];
			char transc = 'C';
			char transn = 'N';
			char transt = 'T';
			if(m==1 && NPOL ==1)
			{
				zgemv_(&transc, &npw, &nkb, &ONE, ppcell.vkb.c, &wf.npwx, chig, &inc, &ZERO, becp, &inc);
			}
			else
			{
				zgemm_(&transc,&transn,&nkb,&npm,&npw,&ONE,ppcell.vkb.c,&wf.npwx,chig,&npw,&ZERO,becp,&nkb);
			}
			Parallel_Reduce::reduce_complex_double_pool( becp, nkb * NPOL * m);

			complex<double> *Ps  = new complex<double> [nkb * NPOL * m];
   			ZEROS( Ps, NPOL * m * nkb);
			
			int sum = 0;
    		int iat = 0;
    		for (int it=0; it<ucell.ntype; it++)
    		{
    		    const int Nprojs = ucell.atoms[it].nh;
    		    for (int ia=0; ia<ucell.atoms[it].na; ia++)
    		    {
    		        // each atom has Nprojs, means this is with structure factor;
    		        // each projector (each atom) must multiply coefficient
    		        // with all the other projectors.
    		        for (int ip=0; ip<Nprojs; ip++)
    		        {
    		            for (int ip2=0; ip2<Nprojs; ip2++)
    		            {
							for(int ib = 0; ib < m ; ++ib)
							{
								Ps[(sum + ip2) * m + ib] += 
								ppcell.deeq(CURRENT_SPIN, iat, ip, ip2) * becp[ib * nkb + sum + ip];
							}//end ib
    		            }// end ih
    		        }//end jh 
					sum += Nprojs;
					++iat;
    		    } //end na
    		} //end nt

			if(NPOL==1 && m==1)
			{
				zgemv_(&transn, &npw, &nkb, &ONE, ppcell.vkb.c, &wf.npwx, Ps, &inc, &ONE, hchig, &inc);
			}
			else
			{
				zgemm_(&transn,&transt,&npw,&npm,&nkb,&ONE,ppcell.vkb.c,&wf.npwx,Ps,&npm,&ONE,hchig,&npw);
			}

			delete[] becp;
			delete[] Ps;
		}
	}
	timer::tick("Stochastic_hchi","vnl",'H');



	double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;

	for(int ig = 0; ig < npw * m; ++ig)
	{
		hchig[ig] = (hchig[ig] - Ebar * chig[ig]) / DeltaE;
	}
	
	timer::tick("Stochastic_hchi","hchi_reciprocal",'H');


	return;
}
