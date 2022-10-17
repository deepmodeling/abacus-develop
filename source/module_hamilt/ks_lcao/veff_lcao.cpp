#include "veff_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "src_pw/global.h"

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

    for (int ir = 0; ir < GlobalC::rhopw->nrxx; ir++)
    {
        GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir);
        if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
        {
            GlobalC::pot.vofk_eff1[ir] = GlobalC::pot.vofk(GlobalV::CURRENT_SPIN, ir);
        }
    }

    //--------------------------------------------
    // (3) folding matrix,
    // and diagonalize the H matrix (T+Vl+Vnl).
    //--------------------------------------------

    if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
    {
        Gint_inout inout(GlobalC::pot.vr_eff1, GlobalC::pot.vofk_eff1, this->LM, Gint_Tools::job_type::vlocal_meta);
        this->GG->cal_vlocal(&inout);
    }
    else
    {
        Gint_inout inout(GlobalC::pot.vr_eff1, this->LM, Gint_Tools::job_type::vlocal);
        this->GG->cal_vlocal(&inout);
    }

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
    for (int ir = 0; ir < GlobalC::rhopw->nrxx; ir++)
    {
        GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir);
        if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
        {
            GlobalC::pot.vofk_eff1[ir] = GlobalC::pot.vofk(GlobalV::CURRENT_SPIN, ir);
        }
    }

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
            Gint_inout inout(GlobalC::pot.vr_eff1, GlobalC::pot.vofk_eff1, 0, Gint_Tools::job_type::vlocal_meta);
            this->GK->cal_gint(&inout);
        }
        else
        {
            // vlocal = Vh[rho] + Vxc[rho] + Vl(pseudo)
            Gint_inout inout(GlobalC::pot.vr_eff1, 0, Gint_Tools::job_type::vlocal);
            this->GK->cal_gint(&inout);
        }

        // added by zhengdy-soc, for non-collinear case
        // integral 4 times, is there any method to simplify?
        if (GlobalV::NSPIN == 4)
        {
            for (int is = 1; is < 4; is++)
            {
                for (int ir = 0; ir < GlobalC::rhopw->nrxx; ir++)
                {
                    GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(is, ir);
                    if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
                    {
                        GlobalC::pot.vofk_eff1[ir] = GlobalC::pot.vofk(is, ir);
                    }
                }
                
                if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
                {
                    Gint_inout inout(GlobalC::pot.vr_eff1, GlobalC::pot.vofk_eff1, is, Gint_Tools::job_type::vlocal_meta);
                    this->GK->cal_gint(&inout);
                }
                else
                {
                    Gint_inout inout(GlobalC::pot.vr_eff1, is, Gint_Tools::job_type::vlocal);
                    this->GK->cal_gint(&inout);
                }
            }
        }
    }

    this->GK->folding_vl_k(ik, this->LM);
    ModuleBase::timer::tick("Veff", "contributeHk");
}

}