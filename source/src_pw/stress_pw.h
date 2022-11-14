#ifndef STRESS_PW_H
#define STRESS_PW_H

#include "stress_func.h"
#include "module_elecstate/elecstate.h"

class Stress_PW: public Stress_Func
{
	public :
	
	Stress_PW (const elecstate::ElecState* pelec_in):pelec(pelec_in){};

	//calculate the stress in PW basis
	void cal_stress(ModuleBase::matrix& smearing_sigmatot, const psi::Psi<complex<double>>* psi_in=nullptr);

	protected :

	//call the vdw stress
	void stress_vdw(ModuleBase::matrix& smearing_sigma);   //force and stress calculated in vdw together.

	const elecstate::ElecState* pelec = nullptr;

};
#endif
