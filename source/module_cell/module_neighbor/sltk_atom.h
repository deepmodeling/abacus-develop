#ifndef INCLUDE_FATOM
#define INCLUDE_FATOM

#include <memory>
#include "sltk_util.h"
#include "module_base/timer.h"
#include <vector>

// a class contains the atom position, 
// the type and the index,
class FAtom
{
public:
	double x;
	double y;
	double z;

	int type;
	int natom;

	int cell_x;
	int cell_y;
	int cell_z;

	FAtom();
	FAtom(const double& x_in, const double& y_in, const double& z_in, 
			const int& type_in, const int& natom_in, 
			const int& cell_x_in, const int& cell_y_in, const int& cell_z_in)
	{
		x = x_in;
		y = y_in;
		z = z_in;
		type = type_in;
		natom = natom_in;
		cell_x = cell_x_in;
		cell_y = cell_y_in;
		cell_z = cell_z_in;
	}
	~FAtom()
	{
	}
};

#endif
