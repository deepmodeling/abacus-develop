#ifndef DMK_IO_H
#define DMK_IO_H

#include <string>
#include "module_cell/klist.h"
#include "module_base/complexmatrix.h"
#include "module_base/matrix.h"

namespace ModuleIO
{
	void write_dmk(const K_Vectors& kv,const int& ik,const std::string &fn, const int &precision,std::vector<ModuleBase::ComplexMatrix> &dm_k);
	void read_dmk(const K_Vectors& kv,const int& ik,const std::string &fn,std::vector<ModuleBase::ComplexMatrix> &dm_k);
};

#endif

