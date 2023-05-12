#ifndef NSCF_BAND_H
#define NSCF_BAND_H
#include "module_base/matrix.h"
#include "module_cell/klist.h"
#include "module_cell/parallel_kpoints.h"

namespace ModuleIO
{
	void nscf_band(
		const int &is,
		const std::string &out_band_dir, 
		const int &nks, 
		const int &nband, 
		const double &fermie,
		const ModuleBase::matrix &ekb,
		const K_Vectors& kv,
		const Parallel_Kpoints* Pkpoints);
}

#endif
