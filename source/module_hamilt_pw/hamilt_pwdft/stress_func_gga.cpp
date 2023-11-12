#include "stress_func.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

//calculate the GGA stress correction in PW and LCAO
template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_gga(ModuleBase::matrix& sigma,
                                             ModulePW::PW_Basis* rho_basis,
                                             const Charge* const chr)
{
    ModuleBase::TITLE("Stress_Func","stress_gga");
	ModuleBase::timer::tick("Stress_Func","stress_gga");
     
	int func_type = XC_Functional::get_func_type();
	if (func_type == 0 || func_type == 1)
	{
		ModuleBase::timer::tick("Stress_Func","stress_gga");
		return;
	}

	FPTYPE sigma_gradcorr[3][3];
	std::vector<FPTYPE> stress_gga;
	FPTYPE dum1, dum2;
	ModuleBase::matrix dum3;
	// call gradcorr to evaluate gradient correction to stress
	// the first three terms are etxc, vtxc and v, which
	// is not used here, so dummy variables are used.
    XC_Functional::gradcorr(dum1, dum2, dum3, chr, rho_basis, &GlobalC::ucell, stress_gga, 1);

    for(int l = 0;l< 3;l++)
	{
		for(int m = 0;m< l+1;m++)
		{
			int ind = l*3 + m;
			sigma_gradcorr[l][m] = stress_gga [ind];
			sigma_gradcorr[m][l] = sigma_gradcorr[l][m];
		}
	}

	for(int l = 0;l<3;l++)
	{
		for(int m = 0;m<3;m++)
		{
            Parallel_Reduce::reduce_pool(sigma_gradcorr[l][m]);
		}
	}
		
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
            sigma(i, j) += sigma_gradcorr[i][j] / rho_basis->nxyz;
        }
	}

	ModuleBase::timer::tick("Stress_Func","stress_gga");
	return;
}

template class Stress_Func<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, psi::DEVICE_GPU>;
#endif