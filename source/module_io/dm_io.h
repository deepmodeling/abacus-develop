#ifndef DM_IO_H
#define DM_IO_H

#include <string>
#include "module_cell/unitcell.h"

namespace ModuleIO
{
	
void read_dm(
#ifdef __MPI
	const int nnrg,
	const int* trace_lo,
#endif
	const bool gamma_only_local,
	const int nlocal,
	const int nspin,
	const int &is,
	const std::string &fn,
	double*** DM,
	double** DM_R,
	double& ef,
	const UnitCell* ucell);

void write_dm(
#ifdef __MPI
	const int* trace_lo,
#endif
	const int &is,
	const int &iter,
	const std::string &fn,
	int precision,
	int out_dm,
	double*** DM,
	const double& ef,
	const UnitCell* ucell);

}

#endif

