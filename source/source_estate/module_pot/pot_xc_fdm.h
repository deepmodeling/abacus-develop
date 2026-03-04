//=======================
// AUTHOR : Peize Lin
// DATE :   2025-10-01
//=======================

#ifndef POTXC_FDM_H
#define POTXC_FDM_H

#include "pot_base.h"

namespace elecstate
{

class PotXC_FDM : public PotBase
{
public:

	PotXC_FDM(
		const ModulePW::PW_Basis* rho_basis_in,
		const Charge*const chg_0_in,
		const UnitCell*const ucell);

	void cal_v_eff(
		const Charge*const chg_1,
		const UnitCell*const ucell,
		ModuleBase::matrix& v_eff) override;

	const Charge*const chg_0 = nullptr;
	ModuleBase::matrix v_xc_0;
};

} // namespace elecstate

#endif