#ifndef WRITE_EIG_OCC_H
#define WRITE_EIG_OCC_H
#include "source_base/matrix.h"
#include "source_cell/klist.h"
#include "source_cell/parallel_kpoints.h"

namespace ModuleIO
{
	void write_eig_iter(const ModuleBase::matrix &ekb,
		const ModuleBase::matrix &wg,
		const K_Vectors& kv);

	void write_eig_file(const ModuleBase::matrix &ekb,
		const ModuleBase::matrix &wg,
		const K_Vectors& kv);
}

#endif
