#ifndef OF_STRESS_PW_H
#define OF_STRESS_PW_H

#include "stress_func.h"

class OF_Stress_PW: public Stress_Func
{
	public :
	
	OF_Stress_PW (){};

	//calculate the stress in OFDFT
	void cal_stress(ModuleBase::matrix& sigmatot, const ModuleBase::matrix& wg, ModuleBase::matrix& kinetic_stress, const psi::Psi<complex<double>>* psi_in=nullptr);

	protected :

	//call the vdw stress
	void stress_vdw(ModuleBase::matrix& smearing_sigma);   //force and stress calculated in vdw together.

};
#endif
