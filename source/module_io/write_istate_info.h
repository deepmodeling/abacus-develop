#ifndef WRITE_ISTATE_INFO_H
#define WRITE_ISTATE_INFO_H
#include "module_elecstate/elecstate.h"
#include "module_cell/klist.h"
#include "src_parallel/parallel_kpoints.h"

namespace ModuleIO
{
	void write_istate_info(const elecstate::ElecState* pelec,
		const K_Vectors* kv,
		const Parallel_Kpoints* Pkpoints);
}

#endif
