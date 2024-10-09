#include "esolver_ks_lcao.h"

#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_parameter/parameter.h"

//--------------temporary----------------------------
#include <memory>

#include "module_base/global_function.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"


//-----HSolver ElecState Hamilt--------
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hsolver/hsolver_lcao.h"
//---------------------------------------------------


namespace ModuleESolver
{

//------------------------------------------------------------------------------
//! the 4th function of ESolver_KS_LCAO: init_after_vc
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::init_after_vc(const Input_para& inp, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "init_after_vc");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "init_after_vc");

    ESolver_KS<TK>::init_after_vc(inp, ucell);
    if (inp.mdp.md_prec_level == 2)
    {
        delete this->pelec;
        this->pelec = new elecstate::ElecStateLCAO<TK>(
            &(this->chr),
            &(this->kv),
            this->kv.get_nks(),
            &(this->GG), // mohan add 2024-04-01
            &(this->GK), // mohan add 2024-04-01
            this->pw_rho,
            this->pw_big);

        dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)
            ->init_DM(&this->kv, &this->pv, PARAM.inp.nspin);

        GlobalC::ppcell.init_vloc(GlobalC::ppcell.vloc, this->pw_rho);

        this->pelec->charge->allocate(PARAM.inp.nspin);
        this->pelec->omega = GlobalC::ucell.omega;

        // Initialize the potential.
        if (this->pelec->pot == nullptr)
        {
            this->pelec->pot = new elecstate::Potential(this->pw_rhod,
                                                        this->pw_rho,
                                                        &GlobalC::ucell,
                                                        &(GlobalC::ppcell.vloc),
                                                        &(this->sf),
                                                        &(this->pelec->f_en.etxc),
                                                        &(this->pelec->f_en.vtxc));
        }
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "init_after_vc");
    return;
}


template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
}
