//=======================
// AUTHOR : Peize Lin
// DATE :   2025-10-05
//=======================

#include "pot_cosikr.h"

#include <cmath>

namespace elecstate
{

Pot_Cosikr::Pot_Cosikr(
		const ModulePW::PW_Basis* rho_basis_in,
		const ModuleBase::Vector3<double> &kvec_d_in,
		const std::vector<double> &phase_in,
		const std::vector<double> &amplitude_in)
	:kvec_d(kvec_d_in),
	 phase(phase_in),
	 amplitude(amplitude_in)
{
	this->rho_basis_ = rho_basis_in;
	this->dynamic_mode = true;
	this->fixed_mode = false;
}


void Pot_Cosikr::cal_v_eff(const Charge*const chg, const UnitCell*const ucell, ModuleBase::matrix &v_eff)
{
	ModuleBase::TITLE("Pot_Cosikr", "cal_v_eff");
	ModuleBase::timer::start("Pot_Cosikr", "cal_veff");
	assert(v_eff.nr == this->phase.size());
	assert(v_eff.nr == this->amplitude.size());
	int ir = 0;
	for (int ix = 0; ix < this->rho_basis_->nx; ++ix)
	{
		const double phase_x = this->kvec_d.x * ix / this->rho_basis_->nx;
		for (int iy = 0; iy < this->rho_basis_->ny; ++iy)
		{
			const double phase_xy = phase_x + this->kvec_d.y * iy / this->rho_basis_->ny;
			for (int iz = this->rho_basis_->startz_current; iz < this->rho_basis_->startz_current + this->rho_basis_->nplane; ++iz)
			{
				const double phase_xyz = phase_xy + this->kvec_d.z * iz / this->rho_basis_->nz;
				for(int is=0; is<v_eff.nr; ++is)
					v_eff(is,ir) += this->amplitude[is] * std::cos((phase_xyz + this->phase[is]) * ModuleBase::TWO_PI);
				++ir;
			}
		}
	}
	ModuleBase::timer::end("Pot_Cosikr", "cal_veff");
}

}