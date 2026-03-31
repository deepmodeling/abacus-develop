//=======================
// AUTHOR : Peize Lin
// DATE :   2025-10-01
//=======================

#include "pot_xc_fdm.h"
#include "source_hamilt/module_xc/xc_functional.h"

namespace elecstate
{

PotXC_FDM::PotXC_FDM(
	const ModulePW::PW_Basis* rho_basis_in,
	const Charge*const chg_0_in,
	const UnitCell*const ucell)
	: chg_0(chg_0_in)
{
	this->rho_basis_ = rho_basis_in;
	this->dynamic_mode = true;
	this->fixed_mode = false;

	const std::tuple<double, double, ModuleBase::matrix> etxc_vtxc_v_0
		= XC_Functional::v_xc(this->chg_0->nrxx, this->chg_0, ucell);
	this->v_xc_0 = std::get<2>(etxc_vtxc_v_0);
}

void PotXC_FDM::cal_v_eff(
	const Charge*const chg_1,
	const UnitCell*const ucell,
	ModuleBase::matrix& v_eff)
{
	ModuleBase::TITLE("PotXC_FDM", "cal_veff");
	ModuleBase::timer::start("PotXC_FDM", "cal_veff");

	assert(this->chg_0->nrxx == chg_1->nrxx);
	assert(this->chg_0->nspin == chg_1->nspin);

	Charge chg_01;
	chg_01.set_rhopw(chg_1->rhopw);
	chg_01.allocate(chg_1->nspin, chg_01.kin_density());

	for(int ir=0; ir<chg_01.nrxx; ++ir)
	{
		for(int is=0; is<chg_01.nspin; ++is)
			{ chg_01.rho[is][ir] = chg_0->rho[is][ir] + chg_1->rho[is][ir]; }
		chg_01.rho_core[ir] = chg_0->rho_core[ir] + chg_1->rho_core[ir];
	}

	const std::tuple<double, double, ModuleBase::matrix> etxc_vtxc_v_01
		= XC_Functional::v_xc(chg_01.nrxx, &chg_01, ucell);
	const ModuleBase::matrix &v_xc_01 = std::get<2>(etxc_vtxc_v_01);

	v_eff += v_xc_01 - this->v_xc_0;

	ModuleBase::timer::end("PotXC_FDM", "cal_veff");
}

} // namespace elecstate

