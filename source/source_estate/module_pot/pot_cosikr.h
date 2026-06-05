//=======================
// AUTHOR : Peize Lin
// DATE :   2025-10-05
//=======================

#ifndef POT_COSIKR_H
#define POT_COSIKR_H

#include "pot_base.h"
#include "source_base/vector3.h"


namespace elecstate
{

// ampitude * cos( 2pi*( k * r + phase ) )
class Pot_Cosikr : public PotBase
{
  public:
	Pot_Cosikr(
		const ModulePW::PW_Basis* rho_basis_in,
		const ModuleBase::Vector3<double> &kvec_d_in,
		const std::vector<double> &phase_in,
		const std::vector<double> &amplitude_in);
	
	void cal_v_eff(const Charge*const chg, const UnitCell*const ucell, ModuleBase::matrix &v_eff) override;

  private:
	const ModuleBase::Vector3<double> kvec_d;
	const std::vector<double> phase;
	const std::vector<double> amplitude;
};

}

#endif