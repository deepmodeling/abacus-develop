#include "./stress_func.h"
#include "../module_xc/xc_functional.h"
#include "../module_base/timer.h"
#include "global.h"

//calculate the GGA stress correction in PW and LCAO
void Stress_Func::stress_gga(ModuleBase::matrix& sigma) 
{
	ModuleBase::timer::tick("Stress_Func","stress_gga");
     
	int func_type = XC_Functional::get_func_type();
	if (func_type == 0 || func_type == 1)
	{
		ModuleBase::timer::tick("Stress_Func","stress_gga");
		return;
	}

	double sigma_gradcorr[3][3];
	std::vector<double> stress_gga;
	double dum1, dum2;
	ModuleBase::matrix dum3;
	// call gradcorr to evaluate gradient correction to stress
	// the first three terms are etxc, vtxc and v, which
	// is not used here, so dummy variables are used.
	XC_Functional::gradcorr(dum1, dum2, dum3, stress_gga, 1);

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
			Parallel_Reduce::reduce_double_pool( sigma_gradcorr[l][m] );
		}
	}
		
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			sigma(i,j) += sigma_gradcorr[i][j] / GlobalC::pw.ncxyz;
		}
	}

	ModuleBase::timer::tick("Stress_Func","stress_gga");
	return;
}
