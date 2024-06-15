#include "FORCE.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/timer.h"

template<>
void Force_LCAO<double>::cal_fvl_dphi(
	const bool isforce, 
    const bool isstress,
    const elecstate::Potential* pot_in,
    TGint<double>::type& gint,
    ModuleBase::matrix& fvl_dphi,
	ModuleBase::matrix& svl_dphi)
{   
    ModuleBase::TITLE("Force_LCAO","cal_fvl_dphi");
    ModuleBase::timer::tick("Force_LCAO","cal_fvl_dphi");

    const int nspin = GlobalV::NSPIN;

    fvl_dphi.zero_out();
    svl_dphi.zero_out();

    for(int is=0; is<nspin; ++is)
    {
        const double* vr_eff1 = pot_in->get_effective_v(is);
        const double* vofk_eff1 = nullptr;

        if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
        {
            vofk_eff1 = pot_in->get_effective_vofk(is);

			Gint_inout inout(
					is,
					vr_eff1, 
					vofk_eff1, 
					isforce, 
					isstress, 
					&fvl_dphi, 
					&svl_dphi, 
					Gint_Tools::job_type::force_meta);

			gint.cal_gint(&inout); 
        }
        else
        {
			Gint_inout inout(
					is, 
					vr_eff1, 
					isforce, 
					isstress, 
					&fvl_dphi, 
					&svl_dphi, 
					Gint_Tools::job_type::force);

			gint.cal_gint(&inout); 
        }
    }

    if(isstress)
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(i<j) 
				{
					svl_dphi(j,i) = svl_dphi(i,j);
				}
				svl_dphi(i,j) /= -GlobalC::ucell.omega;
            }
        }
    }

    ModuleBase::timer::tick("Force_LCAO","cal_fvl_dphi");
}
