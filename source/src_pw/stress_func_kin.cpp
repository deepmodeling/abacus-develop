#include"stress_func.h"
#include "global.h"
#include "module_base/timer.h"

//calculate the kinetic stress in PW base
template<typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_kin(ModuleBase::matrix& sigma, const ModuleBase::matrix& wg, const psi::Psi<complex<FPTYPE>>* psi_in)
{
    ModuleBase::TITLE("Stress_Func","stress_kin");
	ModuleBase::timer::tick("Stress_Func","stress_kin");
	
	FPTYPE **gk;
	gk=new FPTYPE* [3];
	int npw;
	FPTYPE s_kin[3][3];
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<3;m++)
		{
			s_kin[l][m]=0.0;
		}
	}
		
	int npwx=0;
	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{
		if(npwx<GlobalC::kv.ngk[ik])npwx=GlobalC::kv.ngk[ik];
	}
		
	gk[0]= new FPTYPE[npwx];
	gk[1]= new FPTYPE[npwx];
	gk[2]= new FPTYPE[npwx];
	FPTYPE factor=ModuleBase::TWO_PI/GlobalC::ucell.lat0;

	for(int ik=0;ik<GlobalC::kv.nks;ik++)
	{
		npw = GlobalC::kv.ngk[ik];
		for(int i=0;i<npw;i++)
		{
			gk[0][i] = GlobalC::wfcpw->getgpluskcar(ik,i)[0] * factor;
			gk[1][i] = GlobalC::wfcpw->getgpluskcar(ik,i)[1] * factor;
			gk[2][i] = GlobalC::wfcpw->getgpluskcar(ik,i)[2] * factor;
		}

		//kinetic contribution

		for(int l=0;l<3;l++)
		{
			for(int m=0;m<l+1;m++)
			{
				for(int ibnd=0;ibnd<GlobalV::NBANDS;ibnd++)
				{
					const std::complex<FPTYPE>* ppsi=nullptr;
					if(psi_in!=nullptr)
					{
						ppsi = &(psi_in[0](ik, ibnd, 0));
					}
					else
					{
						ppsi = &(GlobalC::wf.evc[ik](ibnd, 0));
					}
					for(int i=0;i<npw;i++)
					{
						s_kin[l][m] +=
							wg(ik, ibnd)*gk[l][i]*gk[m][i]
							*(FPTYPE((conj(ppsi[i]) * ppsi[i]).real()));
					}
				}
			}
		}
		   
		//contribution from the nonlocal part
		   
		//stres_us(ik, gk, npw);
	}
		
	//add the US term from augmentation charge derivatives
		
	// addussstres(sigmanlc);
	
	//mp_cast
		
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<l;m++)
		{
			s_kin[m][l]=s_kin[l][m];
		}
	}

	if(INPUT.gamma_only)
	{
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<3;m++)
			{
				s_kin[l][m] *= 2.0*ModuleBase::e2/GlobalC::ucell.omega;
			}
		}
	}
	else 
	{
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<3;m++)
			{
				s_kin[l][m] *= ModuleBase::e2/GlobalC::ucell.omega;
			}
		}
	}

	for(int l=0;l<3;l++)
	{
		for(int m=0;m<3;m++)
		{
			Parallel_Reduce::reduce_double_all( s_kin[l][m] ); //qianrui fix a bug for kpar > 1
		}
	}


	for(int l=0;l<3;l++)
	{
		for(int m=0;m<3;m++)
		{
			sigma(l,m) = s_kin[l][m];
		}
	}
	//do symmetry
	if(ModuleSymmetry::Symmetry::symm_flag == 1)
	{
		GlobalC::symm.stress_symmetry(sigma, GlobalC::ucell);
	}//end symmetry
	
	delete[] gk[0];
	delete[] gk[1];
	delete[] gk[2];
	delete[] gk;
		
	ModuleBase::timer::tick("Stress_Func","stress_kin");
	return;
}

template class Stress_Func<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, psi::DEVICE_GPU>;
#endif