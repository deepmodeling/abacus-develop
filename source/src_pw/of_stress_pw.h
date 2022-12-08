#ifndef OF_STRESS_PW_H
#define OF_STRESS_PW_H

#include "stress_func.h"
#include "module_elecstate/elecstate.h"

class OF_Stress_PW: public Stress_Func<double>
{
	public :
	
	OF_Stress_PW (const elecstate::ElecState* pelec_in):pelec(pelec_in){};

	//calculate the stress in OFDFT
	void cal_stress(ModuleBase::matrix& sigmatot, ModuleBase::matrix& kinetic_stress, const psi::Psi<complex<double>>* psi_in=nullptr);

	protected :

	//call the vdw stress
	void stress_vdw(ModuleBase::matrix& smearing_sigma);   //force and stress calculated in vdw together.

	const elecstate::ElecState* pelec = nullptr;

};
#endif
