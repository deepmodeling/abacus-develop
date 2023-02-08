#ifndef SYMM_OTHER_H
#define SYMM_OTHER_H

#include "module_base/vector3.h"
#include "module_base/global_function.h"
namespace ModuleSymmetry
{
namespace Symm_Other
{
	void print1(const int &ibrav, const double *cel_const, std::ofstream &ofs_running);
	bool right_hand_sense(ModuleBase::Vector3<double> &v1,ModuleBase::Vector3<double> &v2,ModuleBase::Vector3<double> &v3);
	double celvol(const ModuleBase::Vector3<double> &a, 
		const ModuleBase::Vector3<double> &b, const ModuleBase::Vector3<double> &c);


}
}

#endif
