#include "source_lcao/module_lr/potentials/lr_gs_pot.h"

#include "source_base/parallel_grid.h"
#include "source_cell/stru_fac.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_hamilt/module_surchem/surchem.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_pw/module_pwdft/vl_pw.h"

#include <vector>

namespace LR
{

GsPot::GsPot() = default;
GsPot::~GsPot() = default;

void GsPot::refresh(const Input_para& inp, const std::string& xc_kernel, const UnitCell& ucell,
                    ModulePW::PW_Basis* pw_rhod, ModulePW::PW_Basis* pw_rho, const Parallel_Grid& pgrid,
                    const Charge& chg, Structure_Factor& sf, pseudopot_cell_vl& vloc, surchem& solvent,
                    bool refresh_sfac, bool has_float_data)
{
    std::vector<std::string> pot_register;
    if (inp.vl_in_h)
    {
        if (refresh_sfac)
        {   // on the `ks-lr` path `sf`/`vloc` alias the ground-state solver's, which
            // `ESolver_FP::before_scf` already refreshed for the current geometry
            sf.setup(&ucell, pgrid, pw_rhod, has_float_data);
            vloc.init_vloc(ucell, pw_rho);
        }
        pot_register.push_back("local");
    }
    if (inp.vh_in_h)
    {
        pot_register.push_back("hartree");
    }
    pot_register.push_back("xc");

    // initialize the ground state potential
    pot_ = LR_Util::make_unique<elecstate::Potential>(
        pw_rhod, pw_rho, &ucell, &vloc.vloc, &sf, &solvent, &etxc_, &vtxc_);
    pot_->pot_register(pot_register);
    XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func); // set XC type of the ground state
    pot_->init_pot(&chg);                                    // call update_from_charge inside
    if (LR_Util::has_local_xc(xc_kernel))
    {
        XC_Functional::set_xc_type(xc_kernel); // recover the excited state xc kernel type
    }
    if (inp.test_force)
    {
        pot_hartree_ = LR_Util::make_unique<elecstate::Potential>(
            pw_rhod, pw_rho, &ucell, &vloc.vloc, &sf, &solvent, &etxc_, &vtxc_);
        pot_hartree_->pot_register({"hartree"});
        pot_hartree_->init_pot(&chg); // call update_from_charge inside
    }
}

elecstate::Potential* GsPot::pot_hartree() const
{
    return pot_hartree_.get();
}

} // namespace LR
