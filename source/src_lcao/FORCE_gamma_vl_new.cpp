#include "FORCE_gamma.h"
#include "../src_pw/global.h"
#include "../module_base/timer.h"

void Force_LCAO_gamma::cal_fvl_dphi_new(
	double*** DM_in,
	const bool isforce, 
    const bool isstress,
    ModuleBase::matrix& fvl_dphi,
	ModuleBase::matrix& svl_dphi)
{   
    ModuleBase::TITLE("Force_LCAO_gamma","cal_fvl_dphi_new");
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_fvl_dphi_new");
    int istep = 1;
    GlobalC::pot.init_pot(istep, GlobalC::pw.strucFac);
    fvl_dphi.zero_out();
    for(int is=0; is<GlobalV::NSPIN; ++is)
    {
        GlobalV::CURRENT_SPIN = is;
        for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
        {
            GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir);
        }

        this->UHM->GG.cal_force_new(DM_in[is], GlobalC::pot.vr_eff1, fvl_dphi);
    }
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_fvl_dphi_new");
}