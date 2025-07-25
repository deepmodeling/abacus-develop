#ifndef EXX_ABFS_JLE_H
#define EXX_ABFS_JLE_H

#include "exx_abfs.h"
#include "../../source_hamilt/module_xc/exx_info.h"
#include "../../source_basis/module_ao/ORB_atomic_lm.h"

#include <vector>

	class LCAO_Orbitals;
	class UnitCell;

class Exx_Abfs::Jle
{
public:
	static std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
	init_jle(
		const Exx_Info::Exx_Info_Opt_ABFs &info,
		const double kmesh_times, 
		const UnitCell& ucell,
		const LCAO_Orbitals& orb);
};

#endif	// EXX_ABFS_JLE_H
