//==========================================================
// AUTHOR : Peize Lin
// DATE : 2014-04-25
// UPDATE : 2019-04-26
//==========================================================
#ifndef VDWD2_H
#define VDWD2_H

#include "vdwd2_parameters.h"
#include "module_cell/unitcell_pseudo.h"
#include "src_global/vector3.h"
#include <vector>

class Vdwd2
{
public:

	Vdwd2(const UnitCell_pseudo &unit_in, Vdwd2_Parameters &para_in);
	
	void cal_energy();	
	void cal_force();
	void cal_stress();

	const double                       &get_energy()const{ return energy; }
	const std::vector<Vector3<double>> &get_force() const{ return force;  }
	const Matrix3                      &get_stress()const{ return stress; }

private:

	const UnitCell_pseudo &ucell;
	Vdwd2_Parameters &para;

	double energy = 0;
	std::vector<Vector3<double>> force;
	Matrix3 stress;
};

#endif
