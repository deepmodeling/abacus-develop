#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "global.h"
#include "sto_hchi.h" 
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"

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

std::complex<double>* Stochastic_hchi::rp_chi;
std::complex<double>* Stochastic_hchi::rl_chi;

int * Stochastic_hchi:: GRA_index;

Stochastic_hchi::Stochastic_hchi()
{
	initplan = false;
	ortho = false;
	nrxx = 0;
	rp_chi = new std::complex<double> [1];
	rl_chi = new std::complex<double> [1];
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
	nrxx = GlobalC::pw.nrxx;
	nx = GlobalC::pw.nx;
	ny = GlobalC::pw.ny;
	nz = GlobalC::pw.nz;
    if(nrxx != 0)
    {
        delete[] rp_chi;
        delete[] rl_chi;
		delete[] GRA_index;
		GRA_index = new int [GlobalC::wf.npw];
		if(GlobalC::sto_wf.stotype != "pw")
		{
			rp_chi = new std::complex<double> [nrxx];
			rl_chi = new std::complex<double> [nrxx];

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
        ModuleBase::WARNING_QUIT("Stochastic_hchi", "Number of grids should be at least one!");
    }

}

void Stochastic_hchi::get_GRA_index()
{
	
	if(GlobalC::sto_wf.stotype == "pw")
	{
		for(int ig = 0 ; ig < GlobalC::wf.npw; ++ig)
		{
			GRA_index[ig] = GlobalC::pw.ig2fftw[GlobalC::wf.igk(0,ig)]; //GAMMA POINT temporarily
		}
	}
	else
	{
		int ix,iy,iz;
		int ir;
		ModuleBase::GlobalFunc::ZEROS(GRA_index,GlobalC::wf.npw);
		for(int ig = 0 ; ig < GlobalC::wf.npw; ++ig)
		{
			ix = floor(GlobalC::pw.get_G_cartesian_projection(GlobalC::wf.igk(0, ig), 0) + 0.1);
			iy = floor(GlobalC::pw.get_G_cartesian_projection(GlobalC::wf.igk(0, ig), 1) + 0.1);
			iz = floor(GlobalC::pw.get_G_cartesian_projection(GlobalC::wf.igk(0, ig), 2) + 0.1);
			if(ix < 0) ix += nx;
			if(iy < 0) iy += ny;
			if(iz < 0) iz += nz;
			ir = ix * ny * nz + iy * nz + iz;
			GRA_index[ig] = ir;
		}
	}
}

void Stochastic_hchi::orthogonal_to_psi_real(std::complex<double> *wfin, std::complex<double> *wfout, int &ikk)
{

	ModuleBase::TITLE("Stochastic_hchi","orthogonal_to_psi0");
	if(!initplan) ModuleBase::WARNING_QUIT("Stochastic_hchi", "Please init hchi first!");

	ModuleBase::GlobalFunc::DCOPY(wfin,rp_chi,nrxx);
	//BlasConnector::copy(nrxx,wfin,1,rp_chi,1);
	fftw_execute(pf);
	
	std::complex<double> * chig = new std::complex<double> [GlobalC::wf.npw];
	for(int ig = 0; ig < GlobalC::wf.npw; ++ig)
	{
		chig[ig] = rp_chi[GRA_index[ig]]; 
	}

	//orthogonal part
	std::complex<double> sum = 0;
	int inc=1;
	for(int iksb = 0; iksb < GlobalV::NBANDS; ++iksb)
	{
		std::complex<double> *kswf = &GlobalC::wf.evc[ikk](iksb,0); 
		zdotc_(&sum,&GlobalC::wf.npw,kswf,&inc,chig,&inc);

		for(int ig = 0; ig < GlobalC::wf.npw; ++ig)
		{
			chig[ig] -= sum*kswf[ig];
		}
	}

	//test orthogonal in reciprocal space
	//std::complex<double> overlap;
	//for(int iksb = 0; iksb < GlobalV::NBANDS; ++iksb)
	//{
	//	std::complex<double> *kswf = &GlobalC::wf.evc[ikk](iksb,0); 
	//	overlap=0;
	//	for(int ig = 0; ig < GlobalC::wf.npw; ++ig)
	//	{
	//		overlap += conj(kswf[ig]) * chig[ig];
	//	}
	//	std::cout<<"OVERLAP "<<overlap<<std::endl;
	//}

	for(int ig = 0; ig < GlobalC::wf.npw; ++ig)
	{
		rp_chi[GRA_index[ig]] = chig[ig];
	}
	ModuleBase::GlobalFunc::DCOPY(rp_chi,rl_chi,nrxx);

	fftw_execute(pb);

	for(int ir = 0; ir < nrxx ; ++ir)
	{
		wfout[ir] = rl_chi[ir] / static_cast<double>(nrxx);
	}
	

	/*//test orthogonal in real space
	std::complex<double> overlap;
	std::complex<double> * kswf = new std::complex<double> [nrxx];
	fftw_plan pp=fftw_plan_dft_3d(GlobalC::pw.nx,GlobalC::pw.ny,GlobalC::pw.nz,(fftw_complex *)kswf,(fftw_complex *)kswf, FFTW_BACKWARD, FFTW_ESTIMATE);
	for(int iksb = 0; iksb < GlobalV::NBANDS; ++iksb)
	{
		ModuleBase::GlobalFunc::ZEROS(kswf,nrxx);
		for(int ig = 0 ; ig < GlobalC::wf.npw; ++ig)
		{
			kswf[GRA_index[ig]] = GlobalC::wf.evc[ikk](iksb,ig);
		}
		fftw_execute(pp);
		overlap=0;
		for(int ir = 0; ir < nrxx; ++ir)
		{
			overlap += conj(kswf[ir]) * wfout[ir];
		}
		std::cout<<"OVERLAP "<<overlap<<std::endl;
	}
	delete [] kswf;*/

	ortho = true;
	delete[] chig;
	return;
}

void Stochastic_hchi::orthogonal_to_psi_reciprocal(std::complex<double> *wfgin, std::complex<double> *wfgout, int &ikk)
{

	ModuleBase::TITLE("Stochastic_hchi","orthogonal_to_psi0");
	int nchip=GlobalC::sto_wf.nchip;
	int npw = GlobalC::wf.npw;
	for(int ig = 0 ; ig < npw * nchip; ++ig)
	{
		wfgout[ig] = wfgin[ig];
	}

	//orthogonal part
	
	std::complex<double> *sum = new std::complex<double> [GlobalV::NBANDS * nchip];
	char transC='C';
	char transN='N';
	
	//sum(b<GlobalV::NBANDS, a<nchi) = < psi_b | chi_a >
	zgemm_(&transC, &transN, &GlobalV::NBANDS, &nchip, &npw, &ModuleBase::ONE, GlobalC::wf.evc[ikk].c, &GlobalC::wf.npwx, wfgout, &npw, &ModuleBase::ZERO, sum, &GlobalV::NBANDS);
	Parallel_Reduce::reduce_complex_double_pool(sum, GlobalV::NBANDS * nchip);
	
	//psi -= psi * sum
	zgemm_(&transN, &transN, &npw, &nchip, &GlobalV::NBANDS, &ModuleBase::NEG_ONE, GlobalC::wf.evc[ikk].c, &GlobalC::wf.npwx, sum, &GlobalV::NBANDS, &ModuleBase::ONE, wfgout, &npw);
	
	ortho = true;
	delete[] sum;
	return;
}



void Stochastic_hchi::hchi_real(std::complex<double>*chi_in, std::complex<double> *hchi, const int m)
{
	
	double*vr = GlobalC::pot.vr_eff1;  //vr= GlobalC::pot.vrs1 temporarily use cutoff vr.

	//wait for init--------------------------------------	
	double dk1,dk2,dk3;
	dk1 = GlobalC::ucell.tpiba;
	dk2 = GlobalC::ucell.tpiba;
	dk3 = GlobalC::ucell.tpiba;
		
	//---------------------------------------------------
	if(!initplan) ModuleBase::WARNING_QUIT("Stochastic_hchi", "Please init hchi first!");

	ModuleBase::GlobalFunc::ZEROS(hchi,nrxx);
	ModuleBase::GlobalFunc::DCOPY(chi_in, rp_chi, nrxx);
	//BlasConnector::copy(nrxx,chi_in,1,rp_chi,1);
	fftw_execute(pf);

	std::complex<double> * chig = new std::complex<double> [GlobalC::wf.npw];
	ModuleBase::GlobalFunc::ZEROS(chig,GlobalC::wf.npw);
	for(int ig = 0; ig < GlobalC::wf.npw; ++ig)
	{
		chig[ig] = rp_chi[GRA_index[ig]]; 
	}

	//------------------------------------
	//(1) the local potential.
	//------------------------------------

	if(GlobalV::VL_IN_H)
	{
		for(int ir = 0; ir < nrxx; ++ir)
		{
			hchi[ir] += chi_in[ir] * vr[ir] ;
		}
		
	}
	/*std::cout<<"HCHI-------------------------"<<std::endl;

		for(int ir  = 63800 ; ir < 64000; ++ir)
		{
			std::cout<<chi_in[ir] * vr[ir]<<" ";
		}
		std::cout<<std::endl;*/
	
	
	
	//------------------------------------
	//(2) the kinetical energy.
	//------------------------------------
	if(GlobalV::T_IN_H)
	{
		ModuleBase::Vector3<double> gg;
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
	if(GlobalV::VNL_IN_H)
	{
		if ( GlobalC::ppcell.nkb > 0)
		{
			std::complex<double> *becp = new std::complex<double>[ GlobalC::ppcell.nkb * GlobalV::NPOL ];
			ModuleBase::GlobalFunc::ZEROS(becp,GlobalC::ppcell.nkb * GlobalV::NPOL);

			for (int i=0;i< GlobalC::ppcell.nkb;++i)
			{
				const std::complex<double>* p = &GlobalC::ppcell.vkb(i,0);
				zdotc_(&becp[i],&GlobalC::wf.npw,p,&inc,chig,&inc);
				//for (int ig=0; ig< GlobalC::wf.npw; ++ig)
				//{
				//	if(GlobalV::NSPIN!=4) becp[i] += chig[ig] * conj( p[ig] );
				//} 
			}

			//Parallel_Reduce::reduce_complex_double_pool( becp, GlobalC::ppcell.nkb * GlobalV::NPOL);
			std::complex<double> * Ps = new std::complex<double> [GlobalC::ppcell.nkb * GlobalV::NPOL];
			ModuleBase::GlobalFunc::ZEROS( Ps, GlobalC::ppcell.nkb * GlobalV::NPOL );
			int sum = 0;
    		int iat = 0;
    		// this function sum up each non-local pseudopotential located in each atom,
    		// all we need to do is put the right Dij coefficient to each becp, which
    		// is calculated before.
    		for (int it=0; it<GlobalC::ucell.ntype; ++it)
    		{
    		    const int Nprojs = GlobalC::ucell.atoms[it].nh;
    		    for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ++ia)
    		    {
    		        // each atom has Nprojs, means this is with structure factor;
    		        // each projector (each atom) must multiply coefficient
    		        // with all the other projectors.
    		        for (int ip=0; ip<Nprojs; ++ip)
    		        {
    		            for (int ip2=0; ip2<Nprojs; ++ip2)
    		            {
							if(GlobalV::NSPIN!=4)
								Ps[sum+ip2] += GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, ip, ip2) * becp[sum+ip];
    		            }// end ih
    		        }//end jh
					if(GlobalV::NSPIN!=4) sum += Nprojs;
					++iat;
    		    } //end na
    		} //end nt

			// use simple method.
			ModuleBase::GlobalFunc::ZEROS(chig, GlobalC::wf.npw);
			if(GlobalV::NSPIN!=4)
				for(int i=0; i<GlobalC::ppcell.nkb; ++i)
				{
					std::complex<double>* p = &GlobalC::ppcell.vkb(i,0);
					for(int ig=0; ig< GlobalC::wf.npw; ++ig)
					{
						chig[ig] += Ps[i] * p[ig];
					}
				}
			delete[] becp;
			delete[] Ps;
		}
		for(int ig = 0; ig < GlobalC::wf.npw; ++ig)
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

	for(int i = 0; i < nrxx; ++i)
	{
		hchi[i] += rl_chi[i] / static_cast<double>(nrxx);
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
	//std::cout<<std::setw(15)<<sum2 <<std::setw(15)<<sum1<<std::setw(15)<<sum2/sum1<<std::endl;
	//------------------------------------------------------------

	//test hermit property
	//------------------------------------------------------------
	//std::complex<double> sum=0;
	//for(int i = 0 ; i < nrxx; ++i)
	//{
	//	sum+=conj(chi_in[i]) * hchi[i];
	//}
	//std::cout<<sum<<" must be real numebr."<<std::endl; //sum must be a real number
	//------------------------------------------------------------
	return;
}

void Stochastic_hchi:: hchi_reciprocal(std::complex<double> *chig, std::complex<double> *hchig, const int m)
{
	ModuleBase::timer::tick("Stochastic_hchi","hchi_reciprocal");
	
	//---------------------------------------------------

	int npw = GlobalC::wf.npw;
	int npm = GlobalV::NPOL * m;
	int inc = 1;
	//------------------------------------
	//(1) the kinetical energy.
	//------------------------------------
	std::complex<double> *chibg = chig;
	std::complex<double> *hchibg = hchig;
	if(GlobalV::T_IN_H)
	{
		for (int ib = 0; ib < m ; ++ib)
		{
			for (int ig = 0; ig < npw; ++ig)
			{
				hchibg[ig] = GlobalC::wf.g2kin[ig] * chibg[ig];
			}
			chibg += npw;
			hchibg += npw;
		}
	}
	
	//------------------------------------
	//(2) the local potential.
	//------------------------------------
	ModuleBase::timer::tick("Stochastic_hchi","vloc");
	if(GlobalV::VL_IN_H)
	{
		chibg = chig;
		hchibg = hchig;
		for(int ib = 0 ; ib < m ; ++ib)
		{
			ModuleBase::GlobalFunc::ZEROS( GlobalC::UFFT.porter, GlobalC::pw.nrxx);
			GlobalC::UFFT.RoundTrip( chibg, GlobalC::pot.vr_eff1, GRA_index, GlobalC::UFFT.porter );
			for (int ig = 0; ig < npw; ++ig)
			{
				hchibg[ig] += GlobalC::UFFT.porter[ GRA_index[ig] ];
			}
			chibg += npw;
			hchibg += npw;
		}
			
	}
	ModuleBase::timer::tick("Stochastic_hchi","vloc");


	//------------------------------------
	// (3) the nonlocal pseudopotential.
	//------------------------------------
	ModuleBase::timer::tick("Stochastic_hchi","vnl");
	if(GlobalV::VNL_IN_H)
	{
		if ( GlobalC::ppcell.nkb > 0)
		{
			int nkb = GlobalC::ppcell.nkb;
			std::complex<double> *becp = new std::complex<double>[ nkb * GlobalV::NPOL * m ];
			char transc = 'C';
			char transn = 'N';
			char transt = 'T';
			if(m==1 && GlobalV::NPOL ==1)
			{
				zgemv_(&transc, &npw, &nkb, &ModuleBase::ONE, GlobalC::ppcell.vkb.c, &GlobalC::wf.npwx, chig, &inc, &ModuleBase::ZERO, becp, &inc);
			}
			else
			{
				zgemm_(&transc,&transn,&nkb,&npm,&npw,&ModuleBase::ONE,GlobalC::ppcell.vkb.c,&GlobalC::wf.npwx,chig,&npw,&ModuleBase::ZERO,becp,&nkb);
			}
			Parallel_Reduce::reduce_complex_double_pool( becp, nkb * GlobalV::NPOL * m);

			std::complex<double> *Ps  = new std::complex<double> [nkb * GlobalV::NPOL * m];
   			ModuleBase::GlobalFunc::ZEROS( Ps, GlobalV::NPOL * m * nkb);
			
			int sum = 0;
    		int iat = 0;
    		for (int it=0; it<GlobalC::ucell.ntype; it++)
    		{
    		    const int Nprojs = GlobalC::ucell.atoms[it].nh;
    		    for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
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
								GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, ip, ip2) * becp[ib * nkb + sum + ip];
							}//end ib
    		            }// end ih
    		        }//end jh 
					sum += Nprojs;
					++iat;
    		    } //end na
    		} //end nt

			if(GlobalV::NPOL==1 && m==1)
			{
				zgemv_(&transn, &npw, &nkb, &ModuleBase::ONE, GlobalC::ppcell.vkb.c, &GlobalC::wf.npwx, Ps, &inc, &ModuleBase::ONE, hchig, &inc);
			}
			else
			{
				zgemm_(&transn,&transt,&npw,&npm,&nkb,&ModuleBase::ONE,GlobalC::ppcell.vkb.c,&GlobalC::wf.npwx,Ps,&npm,&ModuleBase::ONE,hchig,&npw);
			}

			delete[] becp;
			delete[] Ps;
		}
	}
	ModuleBase::timer::tick("Stochastic_hchi","vnl");



	double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;

	for(int ig = 0; ig < npw * m; ++ig)
	{
		hchig[ig] = (hchig[ig] - Ebar * chig[ig]) / DeltaE;
	}
	
	ModuleBase::timer::tick("Stochastic_hchi","hchi_reciprocal");


	return;
}
