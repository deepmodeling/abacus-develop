#ifndef STRESS_PW_H
#define STRESS_PW_H

#include "stress_func.h"

class Stress_PW:public Stress_Func
{
	public :
	
	Stress_PW (){};
	~Stress_PW (){};

	//calculate the stress in PW basis
	void cal_stress(matrix& sigmatot);

	private :
	//call the vdw stress
	void stress_vdw(matrix& sigma);   //force and stress calculated in vdw together.
	


};
#endif
