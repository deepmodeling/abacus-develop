#include "stress_func.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/timer.h"

//calculate the kinetic stress in PW base
template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_kin(ModuleBase::matrix& sigma,
                                             const ModuleBase::matrix& wg,
                                             ModuleSymmetry::Symmetry* p_symm,
                                             K_Vectors* p_kv,
                                             ModulePW::PW_Basis_K* wfc_basis,
                                             const psi::Psi<complex<FPTYPE>>* psi_in)
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
    for (int ik = 0; ik < p_kv->nks; ik++)
    {
        if (npwx < p_kv->ngk[ik])
            npwx = p_kv->ngk[ik];
    }

    gk[0]= new FPTYPE[npwx];
	gk[1]= new FPTYPE[npwx];
	gk[2]= new FPTYPE[npwx];
    FPTYPE tpiba = ModuleBase::TWO_PI / GlobalC::ucell.lat0;
    FPTYPE twobysqrtpi = 2.0 / std::sqrt(ModuleBase::PI);
    FPTYPE* kfac = new FPTYPE[npwx];

    for (int ik = 0; ik < p_kv->nks; ik++)
    {
        npw = p_kv->ngk[ik];
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < npw; i++)
        {
            gk[0][i] = wfc_basis->getgpluskcar(ik, i)[0] * tpiba;
            gk[1][i] = wfc_basis->getgpluskcar(ik, i)[1] * tpiba;
            gk[2][i] = wfc_basis->getgpluskcar(ik, i)[2] * tpiba;
            if (wfc_basis->erf_height > 0)
            {
                FPTYPE gk2 = gk[0][i] * gk[0][i] + gk[1][i] * gk[1][i] + gk[2][i] * gk[2][i];
                FPTYPE arg = (gk2 - wfc_basis->erf_ecut) / wfc_basis->erf_sigma;
                kfac[i] = 1.0 + wfc_basis->erf_height / wfc_basis->erf_sigma * twobysqrtpi * std::exp(-arg * arg);
            }
            else
            {
                kfac[i] = 1.0;
            }
        }

        // kinetic contribution

        for (int l = 0; l < 3; l++)
        {
            for (int m = 0; m < l + 1; m++)
            {
                for (int ibnd = 0; ibnd < GlobalV::NBANDS; ibnd++)
                {
                    if (std::fabs(wg(ik, ibnd)) < ModuleBase::threshold_wg * wg(ik, 0))
                        continue;
                    const std::complex<FPTYPE>* ppsi = nullptr;
                    ppsi = &(psi_in[0](ik, ibnd, 0));

                    FPTYPE sum = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
                    for (int i = 0; i < npw; i++)
                    {
                        sum += wg(ik, ibnd) * gk[l][i] * gk[m][i] * kfac[i]
                               * (FPTYPE((conj(ppsi[i]) * ppsi[i]).real()));
                    }
                    s_kin[l][m] += sum;
                }
            }
        }

        //contribution from the nonlocal part
		   
		//stres_us(ik, gk, npw);
    }

    // add the US term from augmentation charge derivatives

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
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        p_symm->stress_symmetry(sigma, GlobalC::ucell);
    } // end symmetry

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