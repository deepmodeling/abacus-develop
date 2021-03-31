#include "tools.h"
#include "global.h"
#include "sto_hchi.h" 

int Stochastic_hchi:: nrxx;
int Stochastic_hchi:: nx,Stochastic_hchi::ny,Stochastic_hchi::nz;
fftw_plan Stochastic_hchi:: pb, Stochastic_hchi::pf;
double Stochastic_hchi:: Emin, Stochastic_hchi:: Emax;
bool Stochastic_hchi:: initplan, Stochastic_hchi::ortho;
complex<double>* Stochastic_hchi:: rp_chi, * Stochastic_hchi::rl_chi;
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
	}
	delete[] rp_chi;
	delete[] rl_chi;
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
        rp_chi = new complex<double> [nrxx];
		rl_chi = new complex<double> [nrxx];
		GRA_index = new int [wf.npw];
		pf=fftw_plan_dft_3d(nx,ny,nz,(fftw_complex *)rp_chi,(fftw_complex *)rp_chi, FFTW_FORWARD, FFTW_MEASURE);
		pb=fftw_plan_dft_3d(nx,ny,nz,(fftw_complex *)rl_chi,(fftw_complex *)rl_chi, FFTW_BACKWARD, FFTW_MEASURE);
		initplan = true;
    }
    else
    {
        WARNING_QUIT("Stochastic_hchi", "Number of grids should be at least one!");
    }

}

void Stochastic_hchi::get_GRA_index()
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

void Stochastic_hchi::orthogonal_to_psi_real(complex<double> *wfin, complex<double> *wfout, int &ikk)
{

	TITLE("Stochastic_hchi","orthogonal_to_psi0");
	if(!initplan) WARNING_QUIT("Stochastic_hchi", "Please init hchi first!");

	DCOPY(wfin,rp_chi,nrxx);
	fftw_execute(pf);
	
	complex<double> * chig = new complex<double> [wf.npw];
	for(int ig = 0; ig < wf.npw; ++ig)
	{
		chig[ig] = rp_chi[GRA_index[ig]]; 
	}

	//orthogonal part
	complex<double> sum;
	for(int iksb = 0; iksb < NBANDS; ++iksb)
	{
		complex<double> *kswf = &wf.evc[ikk](iksb,0); 
		sum=0;
		for(int ig = 0; ig < wf.npw; ++ig)
		{
			sum += conj(kswf[ig]) * chig[ig];
		}
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

	fftw_execute(pb);
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




void Stochastic_hchi::hchi_real(complex<double>*chi_in, complex<double> *hchi)
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
	if(VNL_IN_H)
	{
		if ( ppcell.nkb > 0)
		{
			complex<double> *becp = new complex<double>[ ppcell.nkb * NPOL ];
			ZEROS(becp,ppcell.nkb * NPOL);

			for (int i=0;i< ppcell.nkb;++i)
			{
				const complex<double>* p = &ppcell.vkb(i,0);
				for (int ig=0; ig< wf.npw; ++ig)
				{
					if(NSPIN!=4) becp[i] += chig[ig] * conj( p[ig] );
					else
					{
						//We didnot consider it temporarily.
					}	
				} 
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
	for(int i = 0; i < nrxx; ++i)
	{
		hchi[i] += rl_chi[i] / nrxx;
	}

	
	double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
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
	//cout<<sum2 << " "<<sum1<<" "<<sum2/sum1<<"  dddd"<<endl;
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
