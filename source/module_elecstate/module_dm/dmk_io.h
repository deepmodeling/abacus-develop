#ifndef DMK_IO_H
#define DMK_IO_H

#include <string>
#include "module_cell/klist.h"
#include "module_base/complexmatrix.h"
#include "module_base/matrix.h"

namespace elecstate
{
	// write density matrix dm(k) into *.dmk
	void write_dmk(const K_Vectors& kv,const int& ik,const int& nlocal,const std::string &fn,std::vector<ModuleBase::ComplexMatrix> &dm_k);
	// read *.dmk into density matrix dm(k)
	void read_dmk(const K_Vectors& kv,const int& ik,const int& nlocal,const std::string &fn,std::vector<ModuleBase::ComplexMatrix> &dm_k);
};

#endif

