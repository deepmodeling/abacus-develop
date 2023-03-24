#ifndef DM_IO_H
#define DM_IO_H

#include <string>
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_base/abfs-vector3_order.h"

namespace ModuleIO
{
	void read_dm(const int &is, const std::string &fn, double*** DM, double** DM_R);
	void write_dm(const int &is, const int &iter, const std::string &fn, const int &precision, const int &out_dm, double*** DM);
}

#endif

