#include "module_io/rho_io.h"
#include "module_io/cube_io.h"


bool ModuleIO::read_rho(
#ifdef __MPI
		Parallel_Grid* Pgrid,
#endif
		const int &is,
		const int &nspin,
		const std::string &fn,
		double* rho,
		int& nx,
		int& ny,
		int& nz,
		double& ef,
		const UnitCell* ucell,
		int &prenspin)
{
    return ModuleIO::read_cube(
#ifdef __MPI
		Pgrid,
#endif
		is,
		nspin,
		fn,
		rho,
		nx,
		ny,
		nz,
		ef,
		ucell,
        prenspin);
}
