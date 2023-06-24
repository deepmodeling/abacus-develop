#ifndef CUBE_IO_H
#define CUBE_IO_H
#include <string>
#include "module_cell/unitcell.h"
#ifdef __MPI
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#endif

namespace ModuleIO
{
bool read_cube(
#ifdef __MPI
    Parallel_Grid* Pgrid,
#endif
    const int& is,
    const int& nspin,
    const std::string& fn,
    double* data,
    const int& nx,
    const int& ny,
    const int& nz,
    double& ef,
    const UnitCell* ucell,
    int& prenspin);

	void write_cube(
#ifdef __MPI
		const int& bz,
		const int& nbz,
		const int& nplane,
		const int& startz_current,
#endif
		const double* data,
		const int& is,
		const int& nspin,
		const int& iter,
		const std::string& fn,
		const int& nx,
		const int& ny,
		const int& nz,
		const double& ef,
		const UnitCell* ucell,
		const int &precision = 11,
		const int &out_fermi = 1);//mohan add 2007-10-17
}

#endif
