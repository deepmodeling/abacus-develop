#include "build_st_pw.h"
#include "../src_pw/global.h"
#include "../module_base/timer.h"

#include "global_fp.h" // mohan add 2021-01-30

Build_ST_pw::Build_ST_pw(LCAO_Matrix* lm):
LM(lm)
{

}

Build_ST_pw::~Build_ST_pw()
{

}

// be called in LCAO_Hamilt::calculate_STNR_k()
// FUNCTION: calculate the overlap and kinetic matrix
// in localized basis (expanded in plane wave basis).
void Build_ST_pw::set_ST(const int &ik, const char& dtype)
{
	switch (dtype)
	{
		case 'S':
		{
			for(int i=0; i<GlobalV::NLOCAL; i++)
			{
				const int mu = this->LM->ParaV->trace_loc_row[i];
				if(mu < 0)continue;
				for(int j=0; j<GlobalV::NLOCAL; j++)
				{
					const int nu = this->LM->ParaV->trace_loc_col[j];
					if(nu < 0)continue;
					
					if(GlobalV::NSPIN!=4)
					{
						std::complex<double> v = ModuleBase::ZERO;
						for (int ig = 0; ig < GlobalC::kv.ngk[ik]; ig++) 
						{
							v += conj(GlobalC::wf.wanf2[ik](mu, ig)) * GlobalC::wf.wanf2[ik](nu, ig);
						}
						
		//				std::cout << "i=" << i << " j=" << j << " v=" << v << std::endl;
						//-----------------------------------
						// The results are saved in Sloc2.
						// 2 stands for k points.
						//-----------------------------------
						this->LM->Sloc2[ mu * this->LM->ParaV->ncol + nu ] = v;
					}
					else//added by zhengdy-soc
					{
/*						std::complex<double> v0 = ModuleBase::ZERO, v1 = ModuleBase::ZERO, v2 = ModuleBase::ZERO, v3 = ModuleBase::ZERO;
						for (int ig = 0; ig < GlobalC::kv.ngk[ik]; ig++)
						{
							v0 += conj(GlobalC::wf.wanf2[ik](mu, ig)) * GlobalC::wf.wanf2[ik](nu, ig);
							v1 += conj(GlobalC::wf.wanf2[ik](mu, ig)) * GlobalC::wf.wanf2[ik](nu, ig + GlobalC::wf.npwx);
							v2 += conj(GlobalC::wf.wanf2[ik](mu, ig + GlobalC::wf.npwx)) * GlobalC::wf.wanf2[ik](nu, ig);
							v3 += conj(GlobalC::wf.wanf2[ik](mu, ig + GlobalC::wf.npwx)) * GlobalC::wf.wanf2[ik](nu, ig + GlobalC::wf.npwx);
						}
						this->LM->Sloc2_soc(0, mu * this->LM->ParaV->ncol + nu) = v0;
						this->LM->Sloc2_soc(1, mu * this->LM->ParaV->ncol + nu) = v1;
						this->LM->Sloc2_soc(2, mu * this->LM->ParaV->ncol + nu) = v2;
						this->LM->Sloc2_soc(3, mu * this->LM->ParaV->ncol + nu) = v3;*/
						std::complex<double> v0 = ModuleBase::ZERO;
						for (int ig = 0; ig < GlobalC::wf.npwx*GlobalV::NPOL; ig++)
							v0 += conj(GlobalC::wf.wanf2[ik](mu, ig)) * GlobalC::wf.wanf2[ik](nu, ig);
                        this->LM->Sloc2[mu * this->LM->ParaV->ncol + nu] = v0;

					}
				}
			}
			break;
		}
		case 'T':
		{
			//------------------------------------
			//calculate the kinetic energy of ik.
			//------------------------------------

			for(int i=0; i<GlobalV::NLOCAL; i++)
			{
				const int mu = this->LM->ParaV->trace_loc_row[i];
				if(mu < 0)continue;
				for(int j=0; j<GlobalV::NLOCAL; j++)
				{
					const int nu = this->LM->ParaV->trace_loc_col[j];
					if(nu < 0)continue;
					
					std::complex<double> v = ModuleBase::ZERO;
					for (int ig = 0; ig < GlobalC::kv.ngk[ik]; ig++) 
					{
						v += conj(GlobalC::wf.wanf2[ik](mu, ig)) * GlobalC::wf.wanf2[ik](nu, ig) * GlobalC::wfcpw->getgk2(ik,ig) * GlobalC::ucell.tpiba2;
					}
					if(GlobalV::NSPIN==4)
					for (int ig = 0; ig < GlobalC::kv.ngk[ik]; ig++)
					{
						v += conj(GlobalC::wf.wanf2[ik](mu, ig + GlobalC::wf.npwx)) * GlobalC::wf.wanf2[ik](nu, ig + GlobalC::wf.npwx) * GlobalC::wfcpw->getgk2(ik,ig) * GlobalC::ucell.tpiba2;
					}
					
	//				std::cout << "i=" << i << " j=" << j << " v=" << v << std::endl;
					//-----------------------------------------
					// The results are saved in Hloc_fixed2.
					//-----------------------------------------
                    this->LM->Hloc_fixed2[mu * this->LM->ParaV->ncol + nu] = v;
                }
			}
			break;
		}
	}

	return;
}

/*void Build_ST_pw::set_local(const int &ik)
{
	ModuleBase::TITLE("Build_ST_pw","set_local");
	ModuleBase::timer::tick("Build_ST_pw","set_local");
	assert(GlobalV::NLOCAL>0);
	assert(!GlobalV::GAMMA_ONLY_LOCAL);

    const int npw = GlobalC::kv.ngk[ik];
    std::complex<double> *hpsi = new std::complex<double>[npw];
	std::complex<double> *psic = new std::complex<double>[GlobalC::wfcpw->nrxx];

//	ModuleBase::ComplexMatrix vij(GlobalV::NLOCAL, GlobalV::NLOCAL);

	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		if(GlobalV::NSPIN!=4)
		{
			// fft to real space and doing things.
			GlobalC::wfcpw->recip2real(&GlobalC::wf.wanf2[ik](i, 0), psic, ik);
			for (int ir=0; ir< GlobalC::wfcpw->nrxx; ir++)
			{
				psic[ir] *= GlobalC::pot.vr_eff1[ir];
			}

			// (3) fft back to G space.
			GlobalC::wfcpw->real2recip(psic, hpsi, ik);

			for(int j=i; j<GlobalV::NLOCAL; j++)
			{
				std::complex<double> v = ModuleBase::ZERO;
				for(int ig=0; ig<npw; ig++)
				{
					v += conj( GlobalC::wf.wanf2[ik](j,ig) ) * hpsi[ig];
				}
	//			vij(j, i) = v;
				this->LM->set_HSk(j,i,v,'L');
				if(i!=j)
				{
					this->LM->set_HSk(i,j,conj(v),'L');
				}
			}
		}
		else//noncolinear case
		{
			std::complex<double> *psic1 = new std::complex<double>[GlobalC::wfcpw->nrxx];
			delete[] hpsi;
			hpsi = new std::complex<double> [GlobalC::wf.npwx*GlobalV::NPOL];

			// fft to real space and doing things.
			GlobalC::wfcpw->recip2real(&GlobalC::wf.wanf2[ik](i, 0), psic, ik);
			GlobalC::wfcpw->recip2real(&GlobalC::wf.wanf2[ik](i, GlobalC::wf.npwx), psic1, ik);
			std::complex<double> sup,sdown;
			for (int ir=0; ir< GlobalC::wfcpw->nrxx; ir++)
			{
				sup = psic[ir] * (GlobalC::pot.vr_eff(0,ir) + GlobalC::pot.vr_eff(3,ir)) +
					psic1[ir] * (GlobalC::pot.vr_eff(1,ir) - std::complex<double>(0.0,1.0) * GlobalC::pot.vr_eff(2,ir));
				sdown = psic1[ir] * (GlobalC::pot.vr_eff(0,ir) - GlobalC::pot.vr_eff(3,ir)) +
					psic[ir] * (GlobalC::pot.vr_eff(1,ir) + std::complex<double>(0.0,1.0) * GlobalC::pot.vr_eff(2,ir));
				
				psic[ir] = sup;
				psic1[ir] = sdown;
			}
	
			// (3) fft back to G space.
			GlobalC::wfcpw->real2recip(psic, hpsi, ik);
			GlobalC::wfcpw->real2recip(psic, hpsi+GlobalC::wf.npwx, ik);

			for(int j=i; j<GlobalV::NLOCAL; j++)
			{
				std::complex<double> v = ModuleBase::ZERO;
				for(int ig=0; ig<npw; ig++)
				{
					v += conj( GlobalC::wf.wanf2[ik](j,ig) ) * hpsi[ig];
					v += conj( GlobalC::wf.wanf2[ik](j,ig + GlobalC::wf.npwx) ) * hpsi[ig + GlobalC::wf.npwx];
				}
//			vij(j, i) = v;
				this->LM->set_HSk(j,i,v,'L');
				if(i!=j)
				{
					this->LM->set_HSk(i,j,conj(v),'L');
				}
			}
			delete[] psic1;
		}
	}

//	out.printcm_norm("vij",vij,1.0e-5);

    delete[] hpsi;
	delete[] psic;
	ModuleBase::timer::tick("Build_ST_pw","set_local");
	return;
}*/
