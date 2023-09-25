#include "FORCE_gamma.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/timer.h"

void Force_LCAO_gamma::cal_fvl_dphi(
	double*** DM_in,
	const bool isforce, 
    const bool isstress,
    const elecstate::Potential* pot_in,
    ModuleBase::matrix& fvl_dphi,
	ModuleBase::matrix& svl_dphi)
{   
    ModuleBase::TITLE("Force_LCAO_gamma","cal_fvl_dphi");
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_fvl_dphi");
    int istep = 1;
    fvl_dphi.zero_out();
    svl_dphi.zero_out();
    for(int is=0; is<GlobalV::NSPIN; ++is)
    {
        GlobalV::CURRENT_SPIN = is;
        const double* vr_eff1 = pot_in->get_effective_v(GlobalV::CURRENT_SPIN);
        const double* vofk_eff1 = nullptr;
        if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
        {
            vofk_eff1 = pot_in->get_effective_vofk(GlobalV::CURRENT_SPIN);
        }

        if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
        {
            Gint_inout inout(DM_in, is, vr_eff1, vofk_eff1, isforce, isstress, &fvl_dphi, &svl_dphi, Gint_Tools::job_type::force_meta);
            this->UHM->GG.cal_gint(&inout);
        }
        else
        {
            Gint_inout inout(DM_in, is, vr_eff1, isforce, isstress, &fvl_dphi, &svl_dphi, Gint_Tools::job_type::force);
            this->UHM->GG.cal_gint(&inout);
        }
        
    }

    if(isstress)
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(i<j) svl_dphi(j,i) = svl_dphi(i,j);
				svl_dphi(i,j) /= -GlobalC::ucell.omega;
            }
        }
    }
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_fvl_dphi");
}