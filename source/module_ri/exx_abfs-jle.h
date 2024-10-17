#ifndef EXX_ABFS_JLE_H
#define EXX_ABFS_JLE_H

#include "exx_abfs.h"
#include "../module_hamilt_general/module_xc/exx_info.h"
#include "../module_basis/module_ao/ORB_atomic_lm.h"
#include <vector>

	class LCAO_Orbitals;

namespace Exx_Abfs
{
namespace Jle
{
	// jle[T][L][E]
	extern std::vector< std::vector< std::vector< Numerical_Orbital_Lm>>>
	init_jle(
		const Exx_Info::Exx_Info_Opt_ABFs &info,
		const double kmesh_times,
		const LCAO_Orbitals &orb );
}
}

#endif	// EXX_ABFS_JLE_H
