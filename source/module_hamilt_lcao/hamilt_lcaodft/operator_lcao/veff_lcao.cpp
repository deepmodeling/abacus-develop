#include "veff_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace hamilt
{

template class Veff<OperatorLCAO<double>>;

template class Veff<OperatorLCAO<std::complex<double>>>;

template<>
Veff<OperatorLCAO<double>>::~Veff()
{
}

template<>
Veff<OperatorLCAO<std::complex<double>>>::~Veff()
{
    if(this->allocated_pvpR)
    {
        GK->destroy_pvpR();
    }
}

template<>
void Veff<OperatorLCAO<double>>::contributeHR()
{

}

template<>
void Veff<OperatorLCAO<std::complex<double>>>::contributeHR()
{
    ModuleBase::TITLE("Veff", "contributeHR");
    return;
}

template<>
void Veff<OperatorLCAO<double>>::contributeHk(int ik)
{
    ModuleBase::TITLE("Veff", "contributeHk");
    ModuleBase::timer::tick("Veff", "contributeHk");
    // this->hk_fixed_mock(ik);
    // this->hk_update_mock(ik);

    //-----------------------------------------
    //(1) prepare data for this k point.
    // copy the local potential from array.
    //-----------------------------------------
    const double* vr_eff1 = this->pot->get_effective_v(GlobalV::CURRENT_SPIN);
    const double* vofk_eff1 = this->pot->get_effective_vofk(GlobalV::CURRENT_SPIN);

    //--------------------------------------------
    // (3) folding matrix,
    // and diagonalize the H matrix (T+Vl+Vnl).
    //--------------------------------------------

    if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
    {
        Gint_inout inout(vr_eff1, vofk_eff1, Gint_Tools::job_type::vlocal_meta);
        this->GG->cal_vlocal(&inout, this->LM, new_e_iteration);
    }
    else
    {
        Gint_inout inout(vr_eff1, Gint_Tools::job_type::vlocal);
        this->GG->cal_vlocal(&inout, this->LM, new_e_iteration);
    }

    this->new_e_iteration = false;

    ModuleBase::timer::tick("Veff", "contributeHk");
}

template<>
void Veff<OperatorLCAO<std::complex<double>>>::contributeHk(int ik)
{
    ModuleBase::TITLE("Veff", "contributeHk");
    ModuleBase::timer::tick("Veff", "contributeHk");
    //-----------------------------------------
    //(1) prepare data for this k point.
    // copy the local potential from array.
    //-----------------------------------------
    double* vr_eff1 = this->pot->get_effective_v(GlobalV::CURRENT_SPIN);
    double* vofk_eff1 = this->pot->get_effective_vofk(GlobalV::CURRENT_SPIN);

    //--------------------------------------------
    //(2) check if we need to calculate
    // pvpR = < phi0 | v(spin) | phiR> for a new spin.
    //--------------------------------------------
    if (GlobalV::CURRENT_SPIN == this->GK->get_spin())
    {
        // GlobalV::ofs_running << " Same spin, same vlocal integration." << std::endl;
    }
    else
    {
        // GlobalV::ofs_running << " (spin change)" << std::endl;
        this->GK->reset_spin(GlobalV::CURRENT_SPIN);

        // if you change the place of the following code,
        // rememeber to delete the #include
        if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
        {
            Gint_inout inout(vr_eff1, vofk_eff1, 0, Gint_Tools::job_type::vlocal_meta);
            this->GK->cal_gint(&inout);
        }
        else
        {
            // vlocal = Vh[rho] + Vxc[rho] + Vl(pseudo)
            Gint_inout inout(vr_eff1, 0, Gint_Tools::job_type::vlocal);
            this->GK->cal_gint(&inout);
        }

        // added by zhengdy-soc, for non-collinear case
        // integral 4 times, is there any method to simplify?
        if (GlobalV::NSPIN == 4)
        {
            for (int is = 1; is < 4; is++)
            {
                vr_eff1 = this->pot->get_effective_v(is);
                if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
                {
                    vofk_eff1 = this->pot->get_effective_vofk(is);
                }
                
                if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
                {
                    Gint_inout inout(vr_eff1, vofk_eff1, is, Gint_Tools::job_type::vlocal_meta);
                    this->GK->cal_gint(&inout);
                }
                else
                {
                    Gint_inout inout(vr_eff1, is, Gint_Tools::job_type::vlocal);
                    this->GK->cal_gint(&inout);
                }
            }
        }
    }

    this->GK->folding_vl_k(ik, this->LM, kvec_d);
    ModuleBase::timer::tick("Veff", "contributeHk");
}

}