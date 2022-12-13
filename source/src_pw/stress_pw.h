#ifndef STRESS_PW_H
#define STRESS_PW_H

#include "stress_func.h"
#include "module_elecstate/elecstate.h"

template <typename FPTYPE, typename Device = psi::DEVICE_CPU>
class Stress_PW: public Stress_Func<FPTYPE, Device>
{
	public :
	
	Stress_PW (const elecstate::ElecState* pelec_in):pelec(pelec_in){};

	//calculate the stress in PW basis
	void cal_stress(ModuleBase::matrix& smearing_sigmatot, const psi::Psi<complex<FPTYPE>>* psi_in=nullptr, const psi::Psi<complex<FPTYPE>, Device>* d_psi_in=nullptr);

	protected :

	//call the vdw stress
	void stress_vdw(ModuleBase::matrix& smearing_sigma);   //force and stress calculated in vdw together.

	const elecstate::ElecState* pelec = nullptr;

};
#endif
