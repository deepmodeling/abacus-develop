#ifndef WRITE_ISTATE_INFO_H
#define WRITE_ISTATE_INFO_H
#include "module_base/matrix.h"
#include "module_cell/klist.h"
#include "module_cell/parallel_kpoints.h"

namespace ModuleIO
{
	void write_istate_info(const ModuleBase::matrix &ekb,
		const ModuleBase::matrix &wg,
		const K_Vectors& kv,
		const Parallel_Kpoints* Pkpoints);
}

#endif
