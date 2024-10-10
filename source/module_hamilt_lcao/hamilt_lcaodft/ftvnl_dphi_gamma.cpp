#include "FORCE.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_parameter/parameter.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include <unordered_map>
#include "module_base/timer.h"

template<>
void Force_LCAO<double>::cal_ftvnl_dphi(
    const elecstate::DensityMatrix<double, double>* dm,
    const Parallel_Orbitals& pv,
    const UnitCell& ucell,
    ForceStressArrays& fsr,
    const bool isforce, 
	const bool isstress, 
	ModuleBase::matrix& ftvnl_dphi, 
    ModuleBase::matrix& stvnl_dphi,
    Record_adj* ra)
{
    ModuleBase::TITLE("Force_LCAO","cal_ftvnl_dphi");
    ModuleBase::timer::tick("Force_LCAO","cal_ftvnl_dphi");

    const int nlocal = PARAM.globalv.nlocal;
    const int nspin = PARAM.inp.nspin;

    for(int i=0; i<nlocal; i++)
    {
        const int iat = ucell.iwt2iat[i];
        for(int j=0; j<nlocal; j++)
        {
            const int mu = pv.global2local_row(j);
            const int nu = pv.global2local_col(i);

            if (mu >= 0 && nu >= 0 )
            {
                const int index = mu * pv.ncol + nu;
                //contribution from deriv of AO's in T+VNL term

                double sum = 0.0;
                for(int is=0; is<nspin; ++is)
                {
                    //sum += dm2d[is](nu, mu);
                    sum += dm->get_DMK(is+1, 0, nu, mu);
                }
                sum *= 2.0;

                if(isforce)
				{
					ftvnl_dphi(iat,0) += sum * fsr.DHloc_fixed_x[index];
					ftvnl_dphi(iat,1) += sum * fsr.DHloc_fixed_y[index];
					ftvnl_dphi(iat,2) += sum * fsr.DHloc_fixed_z[index];
				}
                if(isstress)
                {
                    stvnl_dphi(0,0) += sum/2.0 * fsr.DHloc_fixed_11[index];
                    stvnl_dphi(0,1) += sum/2.0 * fsr.DHloc_fixed_12[index];
                    stvnl_dphi(0,2) += sum/2.0 * fsr.DHloc_fixed_13[index];
                    stvnl_dphi(1,1) += sum/2.0 * fsr.DHloc_fixed_22[index];
                    stvnl_dphi(1,2) += sum/2.0 * fsr.DHloc_fixed_23[index];
                    stvnl_dphi(2,2) += sum/2.0 * fsr.DHloc_fixed_33[index];   
                }
            }
        }
    }

	if(isstress)
	{
        StressTools::stress_fill(ucell.lat0, ucell.omega, stvnl_dphi);
	}

    ModuleBase::timer::tick("Force_LCAO","cal_ftvnl_dphi");
    return;
}
