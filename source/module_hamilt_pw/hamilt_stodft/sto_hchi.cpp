#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "sto_hchi.h" 
#include "module_base/tool_title.h"
#include "module_base/timer.h"
#include "module_base/parallel_reduce.h"
#include "module_esolver/esolver_sdft_pw.h"


Stochastic_hchi::Stochastic_hchi()
{
	Emin = INPUT.emin_sto;	
	Emax = INPUT.emax_sto;
}

Stochastic_hchi::~Stochastic_hchi()
{
}

void Stochastic_hchi:: init(ModulePW::PW_Basis_K* wfc_basis, K_Vectors* pkv_in)
{
	wfcpw = wfc_basis;
	pkv = pkv_in;
}


void Stochastic_hchi:: hchi(complex<double> *chig, complex<double> *hchig, const int m)
{
	
	
	//---------------------------------------------------

	const int ik = this->current_ik;
	const int current_spin = pkv->isk[ik];
	const int npwx = this->wfcpw->npwk_max;
	const int npw = this->wfcpw->npwk[ik];
	const int npm = GlobalV::NPOL * m;
	const int inc = 1;
	const double tpiba2 = GlobalC::ucell.tpiba2;
	const int nrxx = this->wfcpw->nrxx;
	//------------------------------------
	//(1) the kinetical energy.
	//------------------------------------
	complex<double> *chibg = chig;
	complex<double> *hchibg = hchig;
	if(GlobalV::T_IN_H)
	{
		for (int ib = 0; ib < m ; ++ib)
		{
			for (int ig = 0; ig < npw; ++ig)
			{
				hchibg[ig] = this->wfcpw->getgk2(ik,ig) * tpiba2 * chibg[ig];
			}
			chibg += npwx;
			hchibg += npwx;
		}
	}
	
	//------------------------------------
	//(2) the local potential.
	//------------------------------------
	ModuleBase::timer::tick("Stochastic_hchi","vloc");
	std::complex<double>* porter = new std::complex<double>[nrxx];
	if(GlobalV::VL_IN_H)
	{
		chibg = chig;
		hchibg = hchig;
		const double* pveff = &((*GlobalTemp::veff)(current_spin, 0));
		for(int ib = 0 ; ib < m ; ++ib)
		{
			this->wfcpw->recip2real(chibg, porter, ik);
			for (int ir=0; ir< nrxx; ir++)
			{
				porter[ir] *=  pveff[ir];
			}
			this->wfcpw->real2recip(porter, hchibg, ik, true);
			
			chibg += npwx;
			hchibg += npwx;
		}
			
	}
	delete[] porter;
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
			complex<double> *becp = new complex<double>[ nkb * GlobalV::NPOL * m ];
			char transc = 'C';
			char transn = 'N';
			char transt = 'T';
			if(m==1 && GlobalV::NPOL ==1)
			{
				zgemv_(&transc, &npw, &nkb, &ModuleBase::ONE, GlobalC::ppcell.vkb.c, &npwx, chig, &inc, &ModuleBase::ZERO, becp, &inc);
			}
			else
			{
				zgemm_(&transc,&transn,&nkb,&npm,&npw,&ModuleBase::ONE,GlobalC::ppcell.vkb.c,&npwx,chig,&npwx,&ModuleBase::ZERO,becp,&nkb);
			}
            Parallel_Reduce::reduce_pool(becp, nkb * GlobalV::NPOL * m);

			complex<double> *Ps  = new complex<double> [nkb * GlobalV::NPOL * m];
   			ModuleBase::GlobalFunc::ZEROS( Ps, GlobalV::NPOL * m * nkb);
			
			int sum = 0;
    		int iat = 0;
    		for (int it=0; it<GlobalC::ucell.ntype; it++)
    		{
    		    const int Nprojs = GlobalC::ucell.atoms[it].ncpp.nh;
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
								GlobalC::ppcell.deeq(current_spin, iat, ip, ip2) * becp[ib * nkb + sum + ip];
							}//end ib
    		            }// end ih
    		        }//end jh 
					sum += Nprojs;
					++iat;
    		    } //end na
    		} //end nt

			if(GlobalV::NPOL==1 && m==1)
			{
				zgemv_(&transn, &npw, &nkb, &ModuleBase::ONE, 
						GlobalC::ppcell.vkb.c, &npwx, Ps, &inc, &ModuleBase::ONE, hchig, &inc);
			}
			else
			{
				zgemm_(&transn,&transt,&npw,&npm,&nkb,&ModuleBase::ONE,
						GlobalC::ppcell.vkb.c,&npwx,Ps,&npm,&ModuleBase::ONE,hchig,&npwx);
			}

			delete[] becp;
			delete[] Ps;
		}
	}
	ModuleBase::timer::tick("Stochastic_hchi","vnl");

	return;
}
void Stochastic_hchi:: hchi_norm(complex<double> *chig, complex<double> *hchig, const int m)
{
	ModuleBase::timer::tick("Stochastic_hchi","hchi_norm");

	this->hchi(chig,hchig,m);

	const int ik = this->current_ik;
	const int npwx = this->wfcpw->npwk_max;
	const int npw = this->wfcpw->npwk[ik];
	const double Ebar = (Emin + Emax)/2;
	const double DeltaE = (Emax - Emin)/2;
	for(int ib = 0 ; ib < m ; ++ib)
	{
		for(int ig = 0; ig < npw; ++ig)
		{
			hchig[ib*npwx+ig] = (hchig[ib*npwx+ig] - Ebar * chig[ib*npwx+ig]) / DeltaE;
		}
	}
	ModuleBase::timer::tick("Stochastic_hchi","hchi_norm");
}