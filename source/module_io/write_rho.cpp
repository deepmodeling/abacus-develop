#include "module_io/cube_io.h"
#include "module_io/rho_io.h"

void ModuleIO::write_rho(
#ifdef __MPI
	const int& bz,
	const int& nbz,
	const int& nplane,
	const int& startz_current,
#endif
	const double* rho_save,
	const int& is,
	const int& nspin,
	const int& iter,
	const std::string& fn,
	const int& nx,
	const int& ny,
	const int& nz,
	const double& ef,
	const UnitCell* ucell,
	const int &precision)
{
	ModuleIO::write_cube(
#ifdef __MPI
		bz,
		nbz,
		nplane,
		startz_current,
#endif
		rho_save,
		is,
		nspin,
		iter,
		fn,
		nx,
		ny,
		nz,
		ef,
		ucell,
		precision,
		1);

    return;
}
