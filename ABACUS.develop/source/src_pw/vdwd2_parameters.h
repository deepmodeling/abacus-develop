//==========================================================
// AUTHOR:	Peize Lin
// DATE:	2021-03-09
//==========================================================

#ifndef VDWD2_PARAMETERS_H
#define VDWD2_PARAMETERS_H

#include <map>
#include <string>
#include "src_global/vector3.h"
#include "src_pw/unitcell_pseudo.h"

class Vdwd2_Parameters
{
public:
	Vdwd2_Parameters();

	void init_C6();
	void init_R0();
	void C6_input(const std::string &file, const std::string &unit);
	void R0_input(const std::string &file, const std::string &unit);	

	bool flag_vdwd2;	

	map<std::string,double> C6;
	map<std::string,double> R0;	
	
	double scaling;
	double damping;
	
	std::string model;
	double radius;
	Vector3<int> period;

	void initset(const UnitCell_pseudo &ucell);
};

#endif // define VDWD2_PARAMETERS_H